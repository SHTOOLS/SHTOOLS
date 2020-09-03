#include <iomanip>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "cWrapper.h"

/*
 * This is a very simple example, explaining how to use the c/cpp interface
 *
 * By default arguments in Fortran are passed by reference and in c/cpp by
 * value. To pass a reference to Fortran get the pointer with the address-of
 * operator: &
 *
 * A lot of the fortran subroutines have optional arguments. In c/c++ optional
 * arguments in arbitrary order are not possible. Thus, always all arguments
 * need to be specified. However, you can use nullptr, which is equivalent to
 * not specifying the argument in fortran.
 *
 * Currently, the function calls are quiet ugly by reason of the large amount of
 * arguments. I am also working on a interface, internally using the same calls
 * but with improved signature. However, the development may take some time.
 * 
 *
 */

int
main(int argc, char** argv)
{

  std::string infile = "../ExampleDataFiles/MarsTopo719.shape";

  int lmax = 15;
  int cilm_dim = lmax + 1;

  // In Fortran Cilm has the dimension 2 x cilm_dim x cilm_dim
  // In the C interface we always use 1 D arrays with the same number of elements
  std::vector<double> cilm(2 * cilm_dim * cilm_dim);

  shtools::cSHRead(infile.c_str(),
                   infile.size(),
                   &cilm[0],
                   cilm_dim,
                   &lmax,
                   nullptr,
                   nullptr,
                   0,
                   nullptr,
                   0,
                   0,
                   0,
                   nullptr);

  double lat = 10.0;
  double lon = 30.0;

  double val = shtools::cMakeGridPoint(
    &cilm[0], cilm_dim, lmax, lat, lon, nullptr, nullptr, nullptr);

  std::cout << "diff to python " << val - 3395259.548270001 << std::endl;

  return 0;
}
