
#include <iostream>
#include <iomanip>
#include <math.h>
#include <memory>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

#include "shtoolswrapper.h"


int
main(int argc, char** argv)
{

  std::string infile = "../examples/ExampleDataFiles/MarsTopo719.shape";


  int lmax = 15;
  int n = lmax+1; 
  std::vector<double> mars = cpp_sh_read(infile, lmax );
    
  Eigen::TensorMap<Cilm> mars_tensor(&mars[0], 2, n, n);
  
  for(int i = 0; i < 2; ++i){
      for(int j = 0; j < n; ++j){
          for(int k = 0; k < n; ++k){
              std::cout << std::setw(12) << std::setprecision(3) << mars_tensor(i,j,k) << " ";
          }
          std::cout << std::endl;
      }
      std::cout << std::string(13*n,'-') << std::endl;
  }
  
  double val = cpp_make_grid_point( mars, 10.0, 30.0);
  std::cout <<  std::setprecision(16) << val << std::endl;
  std::cout <<  "diff to python " << val-3395259.548270001 << std::endl;
  
  

  
//   std::vector<double> cindex(90,1);
//   std::vector<double> cilm(2*8*8);
// 
//   
//   cpp_sh_cindex_to_cilm(cindex, cilm, degmax);
  
  return 0;
  
  
 
}
