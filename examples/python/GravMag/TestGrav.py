"""
This script tests the gravity and magnetics routines.
"""
import numpy as np
import matplotlib.pyplot as plt
import pyshtools as pysh

pysh.utils.figstyle()


# ==== MAIN FUNCTION ====

def main():
    TestMakeGravGrid()
    TestNormalGravity()
    TestGravGrad()
    TestFilter()
    TestMakeMagGrid()


# ==== TEST FUNCTIONS ====

def TestMakeGravGrid():
    infile = '../../ExampleDataFiles/gmm3_120_sha.tab'
    clm, lmax, header = pysh.shio.shread(infile, header=True)
    r0 = float(header[0]) * 1.e3
    gm = float(header[1]) * 1.e9
    clm[0, 0, 0] = 1.0
    print(gm, r0)

    geoid = pysh.gravmag.MakeGeoidGridDH(clm, r0, gm,
                                         pysh.constants.Mars.u0.value,
                                         a=pysh.constants.Mars.a.value,
                                         f=pysh.constants.Mars.f.value,
                                         omega=pysh.constants.Mars.omega.value)
    geoid = geoid / 1.e3  # convert to meters
    fig_map = plt.figure()
    plt.imshow(geoid)
    fig_map.savefig('MarsGeoid.png')

    rad, theta, phi, total, pot = pysh.gravmag.MakeGravGridDH(
        clm, gm, r0, lmax=719, a=pysh.constants.Mars.a.value,
        f=pysh.constants.Mars.f.value, lmax_calc=85,
        omega=pysh.constants.Mars.omega.value, normal_gravity=1)
    fig, axes = plt.subplots(2, 2)

    for num, vv, s in ((0, rad, "$g_{r}$"), (1, theta, "$g_{\\theta}$"),
                       (2, phi, "$g_{\phi}$"),
                       (3, total, "Gravity disturbance")):
        if (num == 3):
            axes.flat[num].imshow(vv * 1.e5, vmin=-400, vmax=550)
            # Convert to mGals
        else:
            axes.flat[num].imshow(vv)
        axes.flat[num].set_title(s)
        axes.flat[num].set_xticks(())
        axes.flat[num].set_yticks(())

    fig.savefig('Mars_Grav.png')


def TestNormalGravity():
    gm = pysh.constants.Mars.gm.value
    omega = pysh.constants.Mars.omega.value
    a = pysh.constants.Mars.a.value
    b = pysh.constants.Mars.b.value
    lat = np.arange(-90., 90., 1.)
    ng = np.array(
        [pysh.gravmag.NormalGravity(x, gm, omega, a, b) for x in lat])
    fig = plt.figure()
    plt.plot(lat, ng, '-')
    plt.xlim(-90, 90)
    plt.xlabel('latitude')
    plt.ylabel('$g, m s^{-2}$')
    fig.savefig('Mars_normalgravity.png')


def TestGravGrad():
    # ---- input parameters ----
    lmax = 100
    clm = np.zeros((2, lmax + 1, lmax + 1), dtype=np.float64)
    clm[0, 2, 2] = 1.0
    gm = 1.0
    r0 = 1.0
    a = 1.0
    f = 0.0

    vxx, vyy, vzz, vxy, vxz, vyz = pysh.gravmag.MakeGravGradGridDH(clm, gm, r0,
                                                                   a=a, f=f)

    print("Maximum Trace(Vxx+Vyy+Vzz) = ", np.max(vxx + vyy + vzz))
    print("Minimum Trace(Vxx+Vyy+Vzz) = ", np.min(vxx + vyy + vzz))

    fig, axes = plt.subplots(2, 3)
    fig.suptitle("Gravity gradient tensor", fontsize=10)

    for num, vv, s in ((0, vxx, "$V_{xx}$"), (1, vyy, "$V_{yy}$"),
                       (2, vzz, "$V_{zz}$"), (3, vxy, "$V_{xy}$"),
                       (4, vxz, "$V_{xz}$"), (5, vyz, "$V_{yz}$")):
        axes.flat[num].imshow(vv, vmin=-5, vmax=5)
        axes.flat[num].set_title(s)
        axes.flat[num].set_xticks(())
        axes.flat[num].set_yticks(())

    fig.savefig('GravGrad_C22.png')


def TestFilter():
    half = 80
    r = pysh.constants.Moon.mean_radius.value
    d = r - 40.e3
    deglist = np.arange(1, 200, 1)
    wl = np.zeros(len(deglist) + 1)
    wlcurv = np.zeros(len(deglist) + 1)

    for l in deglist:
        wl[l] = pysh.gravmag.DownContFilterMA(l, half, r, d)
        wlcurv[l] = pysh.gravmag.DownContFilterMC(l, half, r, d)

    fig = plt.figure()
    plt.plot(deglist, wl[1:], 'b-', label='Minimum amplitude')
    plt.plot(deglist, wlcurv[1:], 'r-', label='Minimum curvature')
    plt.xlabel('degree, l')
    plt.ylabel('W')
    plt.legend()

    fig.savefig('Filter.png')


def TestMakeMagGrid():
    infile = '../../ExampleDataFiles/FSU_mars90.sh'
    clm, lmax, header = pysh.shio.shread(infile, header=True, skip=1)
    r0 = float(header[0]) * 1.e3
    a = pysh.constants.Mars.mean_radius.value + 145.0e3
    # radius to evaluate the field

    rad, theta, phi, total, pot = pysh.gravmag.MakeMagGridDH(
        clm, r0, lmax=719, a=a, f=pysh.constants.Mars.f.value, lmax_calc=90)
    fig, axes = plt.subplots(2, 2)

    for num, vv, s in ((0, rad, "$B_{r}$"), (1, theta, "$B_{\\theta}$"),
                       (2, phi, "$B_{\phi}$"), (3, total, "$|B|$")):
        if (num == 3):
            axes.flat[num].imshow(vv, vmin=0, vmax=700)
        else:
            axes.flat[num].imshow(vv)
        axes.flat[num].set_title(s)
        axes.flat[num].set_xticks(())
        axes.flat[num].set_yticks(())

    fig.savefig('Mars_Mag.png')

    ls = np.arange(lmax + 1)
    pspectrum = pysh.gravmag.mag_spectrum(clm, r0, r0)

    fig_spectrum, ax = plt.subplots(1, 1)
    ax.set_xscale('linear')
    ax.set_yscale('log')
    ax.set_xlabel('degree, l')
    ax.set_ylabel('Power')
    ax.grid(True, which='both')

    ax.plot(ls[1:], pspectrum[1:], label='Magnetic power spectrum')
    ax.legend()

    fig_spectrum.savefig('Mars_MagPowerSpectrum.png')


# ==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
