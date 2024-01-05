"""
This script tests the gravity and magnetics routines.
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pyshtools as pysh

pysh.utils.figstyle()


# ==== MAIN FUNCTION ====

def main():
    TestCircle()
    TestWigner()
    TestSHExpandWLSQ()


# ==== TEST FUNCTIONS ====

def TestCircle():
    coord = pysh.utils.MakeCircleCoord(30, 10, 30)
    lat = coord[:, 0]
    lon = coord[:, 1]
    fig = plt.figure()
    plt.axis([-180, 180, -90, 90])
    plt.plot(lon, lat, 'r-', 10, 30, 'r.')
    coord = pysh.utils.MakeCircleCoord(-75, -45, 10)
    plt.plot(coord[:, 1], coord[:, 0], 'b-', -45, -75, 'b.')
    coord = pysh.utils.MakeEllipseCoord(0, 45, 20, 30, 10)
    plt.plot(coord[:, 1], coord[:, 0], 'g-', 45, 0, 'g.')

    fig.savefig('Circles.png')


def TestWigner():
    w3j, jmin, jmax = pysh.utils.Wigner3j(4, 2, 0, 0, 0)
    print("< J, 4, 2 / 0, 0, 0 >")
    print("jmin = ", jmin)
    print("jmax = ", jmax)
    print(w3j)
    w3j, jmin, jmax = pysh.utils.Wigner3j(10, 14, -1, -4, 5)
    print("< J, 10, 14 / -1, -4, 5 >")
    print("jmin = ", jmin)
    print("jmax = ", jmax)
    print(w3j)


def TestSHExpandWLSQ():
    topofile = 'MarsTopo719.shape'
    if len(sys.argv) > 1:
        topofile = os.path.join(sys.argv[1], topofile)
    else:
        topofile = os.path.join('../../ExampleDataFiles', topofile)

    clm = pysh.SHCoeffs.from_file(topofile, lmax=9)
    nmax = 100
    np.random.seed(seed=123456)
    x = 2*np.random.random_sample(nmax)-1
    y = 2*np.random.random_sample(nmax)-1
    z = 2*np.random.random_sample(nmax)-1

    lat = np.rad2deg(np.arctan2(z, np.sqrt(x**2 + y**2)))
    lon = np.rad2deg(np.arctan2(y, x))
    d = clm.expand(lat=lat, lon=lon)
    print('Minimum and maximum d (km) ', d.max()/1.e3, d.min()/1.e3)
    w = np.ones(nmax)

    print('Least squares inversion misit: chi2, chi2' +
          '(weighted, uniform weights),')
    for l in range(10):
        hilm, chi2 = pysh.expand.SHExpandLSQ(d, lat, lon, l)
        # hlm = pyshtools.SHCoeffs.from_array(hilm)

        hilmw, chi2w = pysh.expand.SHExpandWLSQ(d, w, lat, lon, l)
        # hlmw = pyshtools.SHCoeffs.from_array(hilmw)

        print('L = {:d}, chi2 = {:e}, chi2w = {:e}'
              .format(l, chi2, chi2w))


# ==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
