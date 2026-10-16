#!/usr/bin/env python3
"""
Parametric OpenFOAM study: inlet velocity sweep on the pitzDaily
backward-facing step (steady RANS, k-epsilon).

What it does
  1. Copies the pitzDaily tutorial once per inlet velocity.
  2. Edits the inlet boundary conditions:
       U       -> the new velocity
       k       -> scaled with U^2 (same turbulence intensity)
       epsilon -> scaled with U^3 (same turbulence length scale)
  3. Runs blockMesh and the steady solver for each case.
  4. Reads the mesh and result files in Python and works out, at the
     inlet and outlet:
       - static pressure (area-averaged)
       - total pressure p + 0.5|U|^2 (mass-flow-weighted)
  5. Reports the total pressure loss, the loss coefficient
     K = dp0 / (0.5 rho U^2), the pumping power dp0 * Q and how
     each scales with inlet velocity.
  6. Writes results.csv and a plot (if matplotlib is installed).

Usage
  Source your OpenFOAM environment first, then:
    python3 pitzdaily_sweep.py
    python3 pitzdaily_sweep.py --velocities 5 10 15 20 25
    python3 pitzdaily_sweep.py --case /path/to/pitzDaily
"""

import argparse
import csv
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

RHO = 1.2  # kg/m^3, air: converts OpenFOAM's kinematic pressure to Pa
DEPTH = 0.001  # m, thickness of the 2D pitzDaily mesh


# ---------------------------------------------------------------- setup ----
def find_tutorial(user_case):
    """Return the path of a steady (SIMPLE) pitzDaily case."""
    if user_case:
        case = Path(user_case).expanduser()
        if not (case / "system" / "fvSolution").exists():
            sys.exit(f"Not an OpenFOAM case: {case}")
        return case

    # A copy of the case sitting next to this script wins
    local = Path(__file__).resolve().parent / "pitzDaily"
    if (local / "system" / "fvSolution").exists():
        return local

    roots = [os.environ.get("FOAM_TUTORIALS"),
             os.path.join(os.environ.get("WM_PROJECT_DIR", ""), "tutorials")]
    roots = [Path(r) for r in roots if r and Path(r).is_dir()]

    preferred = [
        "incompressible/simpleFoam/pitzDaily",   # openfoam.com, .org <= v10
        "incompressibleFluid/pitzDailySteady",   # openfoam.org v11 and later
    ]
    candidates = [r / p for r in roots for p in preferred]
    for r in roots:                               # any other layout
        candidates += sorted(r.rglob("pitzDaily*"))

    for case in candidates:
        fvsol = case / "system" / "fvSolution"
        # Only accept steady cases, which have a SIMPLE dictionary
        if (fvsol.exists() and "SIMPLE" in fvsol.read_text()
                and (case / "0" / "U").exists()):
            return case
    sys.exit("Could not find a steady pitzDaily tutorial in this OpenFOAM "
             "install. Unzip pitzDaily.zip next to this script and run again.")


def solver_command():
    """Pick the steady solver for this OpenFOAM version."""
    if shutil.which("foamRun"):        # openfoam.org v11 and later
        return ["foamRun"]
    if shutil.which("simpleFoam"):     # openfoam.com, .org v10 and earlier
        return ["simpleFoam"]
    sys.exit("No OpenFOAM solver found on PATH. Source OpenFOAM first.")


def run(cmd, case, log_name):
    """Run an OpenFOAM utility in a case folder, logging its output."""
    log = case / log_name
    env = dict(os.environ, PWD=str(case))
    with open(log, "w") as f:
        result = subprocess.run(cmd, cwd=case, env=env,
                                stdout=f, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        sys.exit(f"'{' '.join(cmd)}' failed in {case}. See {log}")


# ------------------------------------------------------- edit the case -----
def patch_block(text, patch):
    """Return (start, end) of a patch block in a field file."""
    m = re.search(rf"\b{patch}\s*\{{", text)
    if not m:
        raise ValueError(f"no '{patch}' patch found")
    depth, i = 1, m.end()
    while depth:
        depth += {"{": 1, "}": -1}.get(text[i], 0)
        i += 1
    return m.start(), i


def read_inlet_speed(path):
    text = path.read_text()
    s, e = patch_block(text, "inlet")
    m = re.search(r"value\s+uniform\s*\(\s*([-+0-9.eE]+)", text[s:e])
    if not m:
        raise ValueError(f"no fixed inlet velocity in {path}")
    return float(m.group(1))


def set_inlet_velocity(path, v):
    text = path.read_text()
    s, e = patch_block(text, "inlet")
    block = re.sub(r"(value\s+uniform\s*)\([^)]*\)",
                   lambda m: f"{m.group(1)}({v:g} 0 0)", text[s:e], count=1)
    path.write_text(text[:s] + block + text[e:])


def scale_inlet_scalar(path, factor):
    """Multiply the fixed inlet value of a scalar field by factor."""
    if not path.exists():
        return None
    text = path.read_text()
    s, e = patch_block(text, "inlet")
    m = re.search(r"(value\s+uniform\s+)([-+0-9.eE]+)", text[s:e])
    if not m:
        return None
    new = float(m.group(2)) * factor
    block = text[s:e].replace(m.group(0), f"{m.group(1)}{new:.6g}", 1)
    path.write_text(text[:s] + block + text[e:])
    return new


# --------------------------------------------- read OpenFOAM ASCII files ---
def strip_comments(text):
    return re.sub(r"/\*.*?\*/|//[^\n]*", "", text, flags=re.S)


def first_list(text):
    """Body of the first top-level 'N ( ... )' list in a file."""
    m = re.search(r"\n\s*(\d+)\s*\n\s*\(", text)
    if not m:
        raise ValueError("no ASCII list found (is writeFormat binary?)")
    return m.end()


def read_mesh(case):
    mesh = case / "constant" / "polyMesh"

    btext = strip_comments((mesh / "boundary").read_text())
    patches = {}
    for name, body in re.findall(r"(\w+)\s*\{([^}]*)\}", btext):
        n = re.search(r"nFaces\s+(\d+);", body)
        s = re.search(r"startFace\s+(\d+);", body)
        if n and s:
            patches[name] = (int(s.group(1)), int(n.group(1)))

    text = strip_comments((mesh / "owner").read_text())
    i = first_list(text)
    owner = [int(x) for x in text[i:text.index(")", i)].split()]

    text = strip_comments((mesh / "points").read_text())
    pts = [tuple(map(float, t.split()))
           for t in re.findall(r"\(([-+0-9.eE\s]+)\)", text[first_list(text):])]

    text = strip_comments((mesh / "faces").read_text())
    faces = [list(map(int, f.split()))
             for f in re.findall(r"\d+\(([\d\s]+)\)", text[first_list(text):])]

    return patches, owner, pts, faces


def read_internal_field(path):
    """Cell values of a volScalarField or volVectorField."""
    text = strip_comments(path.read_text())
    m = re.search(r"internalField\s+nonuniform\s+List<(\w+)>\s*(\d+)\s*\(", text)
    if not m:
        raise ValueError(f"internalField in {path} is not a nonuniform list")
    body = text[m.end():]
    n = int(m.group(2))
    if m.group(1) == "scalar":
        return [float(x) for x in body.split(None, n)[:n]]
    vecs = re.findall(r"\(([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\)", body)
    return [tuple(map(float, v)) for v in vecs[:n]]


def face_area_vector(face, pts):
    """Area vector of a polygon face (magnitude = area)."""
    x0 = pts[face[0]]
    ax = ay = az = 0.0
    for a, b in zip(face[1:-1], face[2:]):
        u = [pts[a][k] - x0[k] for k in range(3)]
        w = [pts[b][k] - x0[k] for k in range(3)]
        ax += u[1] * w[2] - u[2] * w[1]
        ay += u[2] * w[0] - u[0] * w[2]
        az += u[0] * w[1] - u[1] * w[0]
    return 0.5 * ax, 0.5 * ay, 0.5 * az


def patch_values(case, time, v_inlet, mesh):
    """Static and total pressure at inlet and outlet (kinematic, m^2/s^2).

    Inlet:  U fixed, p zeroGradient  -> face p = next cell's p
    Outlet: p fixed at 0, U zeroGradient -> face U = next cell's U
    """
    patches, owner, pts, faces = mesh
    tdir = case / f"{time:g}"
    p = read_internal_field(tdir / "p")
    U = read_internal_field(tdir / "U")

    out = {}
    for name in ("inlet", "outlet"):
        start, n = patches[name]
        area = p_area = flow = p0_flow = 0.0
        for f in range(start, start + n):
            sf = face_area_vector(faces[f], pts)
            a = math.sqrt(sum(c * c for c in sf))
            c = owner[f]
            if name == "inlet":
                pf, uf = p[c], (v_inlet, 0.0, 0.0)
            else:
                pf, uf = 0.0, U[c]
            flux = abs(sum(uf[k] * sf[k] for k in range(3)))   # volume flow
            p0 = pf + 0.5 * sum(x * x for x in uf)
            area += a
            p_area += pf * a
            flow += flux
            p0_flow += p0 * flux
        out[name] = {"p": p_area / area, "p0": p0_flow / flow, "Q": flow}
    return out


def latest_time(case):
    times = []
    for d in case.iterdir():
        try:
            t = float(d.name)
        except ValueError:
            continue
        if d.is_dir() and t > 0:
            times.append(t)
    return max(times) if times else 0


# ----------------------------------------------------------- analysis ------
def power_law_exponent(x, y):
    """Least-squares slope of ln(y) against ln(x)."""
    lx = [math.log(a) for a in x]
    ly = [math.log(b) for b in y]
    xm, ym = sum(lx) / len(lx), sum(ly) / len(ly)
    return (sum((a - xm) * (b - ym) for a, b in zip(lx, ly)) /
            sum((a - xm) ** 2 for a in lx))


def make_plot(rows, out):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed, so no plot "
              "(install it with: pip install matplotlib)")
        return

    v = [r["U_inlet_m_s"] for r in rows]
    dp0 = [r["total_p_loss_Pa"] for r in rows]
    k = [r["loss_coeff_K"] for r in rows]
    a = sum(x * x * y for x, y in zip(v, dp0)) / sum(x ** 4 for x in v)
    vf = [min(v) + i * (max(v) - min(v)) / 50 for i in range(51)]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    ax1.plot(v, dp0, "o", label="CFD")
    ax1.plot(vf, [a * x * x for x in vf], "-",
             label=f"fit: $\\Delta p_0$ = {a:.3f} U$^2$")
    ax1.set_xlabel("Inlet velocity U (m/s)")
    ax1.set_ylabel("Total pressure loss (Pa)")
    ax1.set_title("Total pressure loss vs inlet velocity")
    ax1.legend()
    ax1.grid(alpha=0.3)

    ax2.plot(v, k, "s-")
    ax2.set_xlabel("Inlet velocity U (m/s)")
    ax2.set_ylabel("Loss coefficient K")
    ax2.set_title("K = $\\Delta p_0$ / (0.5 $\\rho$ U$^2$)")
    ax2.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print(f"Plot saved: {out}")


# --------------------------------------------------------------- main ------
def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--velocities", type=float, nargs="+",
                    default=[5, 10, 15, 20], help="inlet velocities in m/s")
    ap.add_argument("--case", help="path to a steady pitzDaily case")
    ap.add_argument("--out", default="runs", help="output folder")
    args = ap.parse_args()

    template = find_tutorial(args.case)
    solver = solver_command()
    out = Path(args.out).resolve()
    out.mkdir(exist_ok=True)
    print(f"Template case : {template}")
    print(f"Solver        : {' '.join(solver)}\n")

    u0 = read_inlet_speed(template / "0" / "U")
    rows, mesh = [], None

    for v in args.velocities:
        case = out / f"U_{v:g}"
        if case.exists():
            shutil.rmtree(case)
        shutil.copytree(template, case)

        # Streamline output is not needed here and slows the run down
        cdict = case / "system" / "controlDict"
        cdict.write_text(re.sub(r"^\s*#includeFunc\s+streamlines.*$", "",
                                cdict.read_text(), flags=re.M))

        set_inlet_velocity(case / "0" / "U", v)
        k_in = scale_inlet_scalar(case / "0" / "k", (v / u0) ** 2)
        e_in = scale_inlet_scalar(case / "0" / "epsilon", (v / u0) ** 3)

        print(f"U = {v:>5g} m/s  k_in = {k_in:.4g}  eps_in = {e_in:.4g} ... ",
              end="", flush=True)
        run(["blockMesh"], case, "log.blockMesh")
        run(solver, case, "log.solver")

        if mesh is None:              # same mesh for every case
            mesh = read_mesh(case)
        t = latest_time(case)
        res = patch_values(case, t, v, mesh)

        dp0 = RHO * (res["inlet"]["p0"] - res["outlet"]["p0"])
        rows.append({
            "U_inlet_m_s": v,
            "iterations": int(t),
            "static_p_rise_Pa": round(RHO * (res["outlet"]["p"] - res["inlet"]["p"]), 3),
            "total_p_loss_Pa": round(dp0, 3),
            "loss_coeff_K": round(dp0 / (0.5 * RHO * v * v), 4),
            "power_W_per_m": round(dp0 * res["inlet"]["Q"] / DEPTH, 3),
            "mass_balance_error_pct": round(
                100 * abs(res["inlet"]["Q"] - res["outlet"]["Q"]) / res["inlet"]["Q"], 3),
        })
        print(f"{int(t)} iterations, total pressure loss {dp0:.2f} Pa, "
              f"K = {rows[-1]['loss_coeff_K']}")

    csv_path = out / "results.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    print(f"\nResults saved: {csv_path}")

    if len(rows) >= 2:
        n = power_law_exponent([r["U_inlet_m_s"] for r in rows],
                               [r["total_p_loss_Pa"] for r in rows])
        print(f"Total pressure loss scales as U^{n:.2f} "
              "(2.0 = fully turbulent, Reynolds-independent losses)")
        m = power_law_exponent([r["U_inlet_m_s"] for r in rows],
                               [r["power_W_per_m"] for r in rows])
        print(f"Pumping power scales as U^{m:.2f} (about 3 expected)")
        make_plot(rows, out / "loss_vs_velocity.png")


if __name__ == "__main__":
    main()
