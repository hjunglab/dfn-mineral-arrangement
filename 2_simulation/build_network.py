#!/usr/bin/env python3
"""Generate and mesh one stochastic discrete fracture network with dfnWorks.

    python build_network.py <seed> [min_fractures=12]

Output: $WORK/networks/s<seed>/  (conforming mesh, PFLOTRAN grid files, boundary files, pinned.ex)
Prints 'SUCCESS <dir> <n_fractures>' or 'FAIL <reason>' (exit code 2), so that a driver loop can
simply try the next seed when a realization does not percolate or cannot be meshed.

Requires the dfnWorks Python package (pydfnworks, v2.9); run inside the dfnWorks Docker container.
"""
import os, sys, json

HERE = os.path.dirname(os.path.abspath(__file__))
INPUTS = os.path.join(HERE, "..", "1_inputs")
WORK = os.environ.get("WORK", os.path.join(HERE, "..", "work"))
FAM = json.load(open(os.path.join(INPUTS, "fracture_families.json")))


def fail(msg):
    print("FAIL %s" % msg); sys.exit(2)


def build(seed):
    from pydfnworks import DFNWORKS
    job = os.path.join(WORK, "networks", "s%d" % seed)
    os.makedirs(os.path.dirname(job), exist_ok=True)
    dfn = DFNWORKS(job, dfnFlow_file=os.path.join(INPUTS, "flow.in"), ncpu=int(os.environ.get("NCPU", "4")))
    p = dfn.params
    p["domainSize"]["value"] = FAM["domain_m"]
    p["h"]["value"] = FAM["mesh_edge_length_m"]
    p["seed"]["value"] = seed
    p["stopCondition"]["value"] = 1                      # generate until the target P32 is reached
    p["domainSizeIncrease"]["value"] = [1.0, 1.0, 0.5]
    p["numOfLayers"]["value"] = 1
    p["layers"]["value"] = [FAM["layer_z_m"]]
    p["keepOnlyLargestCluster"]["value"] = 1             # keep the percolating cluster only
    p["ignoreBoundaryFaces"]["value"] = False
    p["boundaryFaces"]["value"] = [1, 1, 0, 0, 0, 0]     # flow along x: inlet (west) and outlet (east)
    p["orientationOption"]["value"] = 2                  # dip / strike

    ap, size = FAM["aperture"], FAM["size"]
    mu_h = ap["measured_geometric_mean_m"] * ap["hydraulic_factor"]   # hydraulic aperture, geometric mean ~0.10 mm
    for fam in FAM["families"].values():
        dfn.add_fracture_family(
            shape="ell", distribution="tpl", alpha=size["alpha"],
            min_radius=size["min_radius_m"], max_radius=size["max_radius_m"],
            kappa=fam["fisher_kappa"], dip=fam["mean_dip_deg"], strike=fam["mean_strike_deg"],
            p32=fam["p32_per_m"], layer=1, aspect=size["aspect_ratio"], number_of_points=8,
            hy_variable="aperture", hy_function="log-normal",
            hy_params={"mu": mu_h, "sigma": ap["ln_variance"]})
    dfn.make_working_directory(delete=True)
    dfn.check_input()
    dfn.create_network()
    dfn.mesh_network(uniform_mesh=True)
    dfn.lagrit2pflotran(boundary_cell_area=1)            # unit boundary cell area keeps the rate-BC matrix well conditioned
    return job


if __name__ == "__main__":
    seed = int(sys.argv[1]); min_frac = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    try:
        job = build(seed)
    except SystemExit as e:                              # dfnGen exits when the network does not connect the boundaries
        fail("not_connected_%s" % e)
    except Exception as e:
        fail("generation_%s" % repr(e)[:60])

    # pressure datum: first cell of the bottom boundary
    lines = open(os.path.join(job, "boundary_bottom.ex")).read().splitlines()
    first = next((ln for ln in lines[1:] if ln.strip()), None)
    if first is None:
        fail("empty_bottom_boundary")
    open(os.path.join(job, "pinned.ex"), "w").write("CONNECTIONS\t1\n" + first + "\n")

    mat = open(os.path.join(job, "materialid.dat")).read().splitlines()[3:]
    nfrac = len(set(int(l.split()[0]) for l in mat))
    if nfrac < min_frac:
        fail("too_few_fractures_%d" % nfrac)
    print("SUCCESS %s %d" % (job, nfrac))
