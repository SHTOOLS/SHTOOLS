#pragma once

constexpr int n2 = 2;

extern "C"
{
     void cbind_sh_cindex_to_cilm( const double* cindex, const int* cindex_d0, const int* cindex_d1,
                                   double* cilm, const int* cilm_d0,
                                   const int* degmax,
                                   int* exitstatus);
     
     void cbind_sh_read(    const char* filename, 
                            int* flen, 
                            double* cilm, const int* cilm_d0,
                            int* lmax,
                            int* skip,
                            double* header, const int* header_d1,
                            double* error, const int* error_d0, const int* error_d1, const int* error_d2,
                            int* exitstatus);
     
     double cbind_make_grid_point( const double* cilm, const int* cilm_d0,
                                   const int* lmax,
                                   const double* lat, const double* lon,
                                   const int* norm,
                                   const int* csphase,
                                   const int* dealloc);

}

inline int deg_2_n(int deg){
    return std::sqrt( 1./4. + 2.*deg ) - 3./2.;
}


inline double cpp_make_grid_point( const std::vector<double>& cilm, double lat, double lon, 
                                     int norm=1, int csphase=1, int dealloc=0)
{

   int n = cilm.size();
   int cilmd = std::sqrt(n/2);
   int lmax = cilmd-1;
      
   return cbind_make_grid_point (&cilm[0], &cilmd, &lmax, &lat, &lon, &norm, &csphase, &dealloc);

 }


 inline void cpp_sh_cindex_to_cilm( std::vector<double> cindex, std::vector<double> cilm, int degmax )
 {
  

    int exitstatus;
    
    int n = cindex.size()/2;
    int m = deg_2_n(n);
    
    cbind_sh_cindex_to_cilm( &cindex[0], &n2, &n, &cilm[0], &m, &degmax, &exitstatus );
    


 }
 
  inline std::vector<double> cpp_sh_read( const std::string& filename, int degree  )
  {
  

    int exitstatus;       
    int cilm_dim = degree+1;

    std::vector<double> cilm(2*cilm_dim*cilm_dim);
    std::vector<double> error(2*cilm_dim*cilm_dim);
    
    double* header = nullptr;
    int* header_d = nullptr;
    int* skip = nullptr;

    
    int s = filename.size();
    
    cbind_sh_read(
                      filename.c_str(),
                      &s,
                      &cilm[0], &cilm_dim,
                      &degree,
                      skip,
                      header, header_d,
                      &error[0], &n2, &cilm_dim, &cilm_dim,
                      &exitstatus
                     );
    
    return cilm;

 }
 
 
 typedef Eigen::Matrix<double, 2, Eigen::Dynamic, Eigen::ColMajor> Cindex;
 typedef Eigen::Tensor<double, 3, Eigen::ColMajor> Cilm;
 
//  inline void cpp_shc_index_to_cilm1( Cindex cindex, Cilm cilm, int degmax )
//  {
//      
//  }


