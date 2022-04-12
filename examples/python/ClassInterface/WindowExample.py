"""
This script tests the python class interface
"""
import pyshtools as pysh

pysh.utils.figstyle()


# ==== MAIN FUNCTION ====

def main():
    example1()
    example2()


# ==== EXAMPLES ====
def example1():
    # generate cap window
    lmax = 20
    nwin = 20
    theta = 25.
    cap = pysh.SHWindow.from_cap(theta, lmax, nwin=nwin)
    cap.info()
    cap.plot_windows(20, show=False, fname='cap_tapers.png')
    cap.plot_coupling_matrix(30, k=5, show=False, fname='cap_coupling.png')


# ==== EXAMPLES ====
def example2():
    # generate cap window
    lmax = 15
    nwins = 15

    coeffs = pysh.SHCoeffs.from_file(
        '../../ExampleDataFiles/srtmp300.msl')
    topo = coeffs.expand(grid='DH2')
    dh_mask = topo.data > 0.
    print(dh_mask.shape)
    region = pysh.SHWindow.from_mask(dh_mask, lmax, nwins)
    region.info()
    region.plot_windows(nwins, show=False, fname='continent_tapers.png')
    region.plot_coupling_matrix(30, k=5, show=False,
                                fname='continent_coupling.png')


# ==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
