"""
This script creates a crustal thickness map of Mars.
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pyshtools as pysh

pysh.utils.figstyle()


# ==== MAIN FUNCTION ====

def main():
    TestCrustalThickness()

# ==== TEST FUNCTIONS ====


def TestCrustalThickness():
    """
    Example routine that calculates the crustal thickness of Mars
    """
    delta_max = 5.0
    nmax = 6
    degmax = 50
    lmax = 200
    rho_c = 2900.0
    rho_m = 3500.0
    filter_type = 0
    half = 0

    gravfile = 'gmm3_120_sha.tab'
    if len(sys.argv) > 1:
        gravfile = os.path.join(sys.argv[1], gravfile)
    else:
        gravfile = os.path.join('../../ExampleDataFiles', gravfile)

    pot, lmaxp, header = pysh.shio.shread(gravfile, lmax=degmax, header=True)
    gm = float(header[1]) * 1.e9
    mass = gm / pysh.constants.G.value
    r_grav = float(header[0]) * 1.e3
    print(r_grav, gm, mass, lmaxp)

    topofile = 'MarsTopo719.shape'
    if len(sys.argv) > 1:
        topofile = os.path.join(sys.argv[1], topofile)
    else:
        topofile = os.path.join('../../ExampleDataFiles', topofile)

    hlm, lmaxt = pysh.shio.shread(topofile)
    r0 = hlm[0, 0, 0]
    d = r0 - 45.217409924028445e3
    print(r0, lmaxt)

    for l in range(2, lmaxp + 1):
        pot[:, l, :l + 1] = pot[:, l, :l + 1] * (r_grav / r0)**l

    topo_grid = pysh.expand.MakeGridDH(hlm, lmax=lmax, sampling=2,
                                       lmax_calc=degmax)
    print("Maximum radius (km) = ", topo_grid.max() / 1.e3)
    print("Minimum radius (km) = ", topo_grid.min() / 1.e3)

    bc, r0 = pysh.gravmag.CilmPlusDH(topo_grid, nmax, mass, rho_c, lmax=degmax)

    ba = pot - bc

    moho_c = np.zeros([2, degmax + 1, degmax + 1], dtype=np.float64)
    moho_c[0, 0, 0] = d

    for l in range(1, degmax + 1):
        if filter_type == 0:
            moho_c[:, l, :l + 1] = ba[:, l, :l + 1] * mass * (2 * l + 1) * \
                ((r0 / d)**l) / (4.0 * np.pi * (rho_m - rho_c) * d**2)
        elif filter_type == 1:
            moho_c[:, l, :l + 1] = pysh.gravmag.DownContFilterMA(
                l, half, r0, d) * ba[:, l, :l + 1] * mass * (2 * l + 1) * \
                ((r0 / d)**l) / (4.0 * np.pi * (rho_m - rho_c) * d**2)
        else:
            moho_c[:, l, :l + 1] = pysh.gravmag.DownContFilterMC(
                l, half, r0, d) * ba[:, l, :l + 1] * mass * (2 * l + 1) *\
                ((r0 / d)**l) / (4.0 * np.pi * (rho_m - rho_c) * d**2)

    moho_grid3 = pysh.expand.MakeGridDH(moho_c, lmax=lmax, sampling=2,
                                        lmax_calc=degmax)
    print('Maximum Crustal thickness (km) = ',
          (topo_grid - moho_grid3).max() / 1.e3)
    print('Minimum Crustal thickness (km) = ',
          (topo_grid - moho_grid3).min() / 1.e3)

    moho_c = pysh.gravmag.BAtoHilmDH(ba, moho_grid3, nmax, mass, r0,
                                     (rho_m - rho_c), lmax=lmax,
                                     filter_type=filter_type, filter_deg=half,
                                     lmax_calc=degmax)

    moho_grid2 = pysh.expand.MakeGridDH(moho_c, lmax=lmax, sampling=2,
                                        lmax_calc=degmax)
    print('Delta (km) = ', abs(moho_grid3 - moho_grid2).max() / 1.e3)

    temp_grid = topo_grid - moho_grid2
    print('Maximum Crustal thickness (km) = ', temp_grid.max() / 1.e3)
    print('Minimum Crustal thickness (km) = ', temp_grid.min() / 1.e3)

    iter = 0
    delta = 1.0e9

    while delta > delta_max:
        iter += 1
        print('Iteration ', iter)

        moho_grid = (moho_grid2 + moho_grid3) / 2.0
        print("Delta (km) = ", abs(moho_grid - moho_grid2).max() / 1.e3)

        temp_grid = topo_grid - moho_grid
        print('Maximum Crustal thickness (km) = ', temp_grid.max() / 1.e3)
        print('Minimum Crustal thickness (km) = ', temp_grid.min() / 1.e3)

        moho_grid3 = moho_grid2
        moho_grid2 = moho_grid

        iter += 1
        print('Iteration ', iter)

        moho_c = pysh.gravmag.BAtoHilmDH(ba, moho_grid2, nmax, mass, r0,
                                         rho_m - rho_c, lmax=lmax,
                                         filter_type=filter_type,
                                         filter_deg=half,
                                         lmax_calc=degmax)

        moho_grid = pysh.expand.MakeGridDH(moho_c, lmax=lmax, sampling=2,
                                           lmax_calc=degmax)
        delta = abs(moho_grid - moho_grid2).max()
        print('Delta (km) = ', delta / 1.e3)

        temp_grid = topo_grid - moho_grid
        print('Maximum Crustal thickness (km) = ', temp_grid.max() / 1.e3)
        print('Minimum Crustal thickness (km) = ', temp_grid.min() / 1.e3)

        moho_grid3 = moho_grid2
        moho_grid2 = moho_grid

        if temp_grid.max() > 100.e3:
            print('Not converging')
            exit(1)

    fig_map = plt.figure()
    plt.imshow(temp_grid)
    fig_map.savefig('Mars_CrustalThicknes.png')


# ==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
