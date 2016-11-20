#ifndef MEXPLUSPLUS_EIGEN_HPP
#define MEXPLUSPLUS_EIGEN_HPP

#include <mexplus.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

#include "mexplusplus/sparse_matlab_eigen.hpp"

template <typename Scalar>
using Vector3 = Eigen::Matrix<Scalar,3,1>; 

namespace mexplus  {
  // Define two template specializations.

  template <> 
  mxArray* MxArray::from(const std::vector<Vector3<double> >& u) 
  {
    //    typedef typename Derived::Scalar Scalar;
    typedef double Scalar;
    size_t n = u.size();
    if (n > 0) {
      size_t m = u[0].rows();
      MxArray numeric(MxArray::Numeric<Scalar>(m,n));
      size_t j = 0;
      for (auto u_i : u) {
	for  (size_t i=0;i<m;i++)  // rows
	  numeric.set(i,j,u_i.derived().coeff(i,0));
	j++;
      }
      return numeric.release();
    }
  }

  template <>
  void MxArray::to(const mxArray* array, std::vector<Vector3<double> >* u) 
  {
    //typedef typename Derived::Scalar Scalar;
    // Write your conversion code. For example,
    typedef double Scalar;

    MxArray numeric(array);
    size_t m = numeric.rows();
    size_t n = numeric.cols();
    u->reserve(m);
    for (size_t j = 0; j < n; j++) {
      Vector3<double> u_i;
      for (size_t i = 0; i < m; i++) 
	u_i[i] = numeric.at<double>(i,j);    
      u->push_back(u_i);
    }
  }

  template <>  
  mxArray* 
  MxArray::from(const std::vector<Vector3<float> >& u) 
  {
    //    typedef typename Derived::Scalar Scalar;
    typedef float Scalar;
    size_t n = u.size();
    if (n > 0) {
      size_t m = u[0].rows();
      MxArray numeric(MxArray::Numeric<Scalar>(m,n));
      size_t j = 0;
      for (auto u_i : u) {
	for (size_t i=0; i<m; i++)  // rows
	  numeric.set(i,j,u_i.derived().coeff(i,0));
	j++;
      }
      return numeric.release();
    }
  }

  template <>
  void 
  MxArray::to(const mxArray* array, std::vector<Vector3<float> >* u) 
  {
    //typedef typename Derived::Scalar Scalar;
    // Write your conversion code. For example,
    typedef float Scalar;

    MxArray numeric(array);
    size_t m = numeric.rows();
    size_t n = numeric.cols();
    u->reserve(m);
    for (size_t j = 0; j < n; j++) {
      Vector3<float> u_i;
      for (size_t i = 0; i < m; i++) 
	u_i[i] = numeric.at<float>(i,j);    
      u->push_back(u_i);
    }
  }

  // Define template specializations for doubles.
  template <>
  mxArray* 
  MxArray::from(const Eigen::SparseMatrix<double, Eigen::ColMajor>& u) 
  {
    return eigen2matlab<double, Eigen::ColMajor>(u);
  }

  template <>
  void 
  MxArray::to(const mxArray* array, Eigen::SparseMatrix<double, Eigen::ColMajor>* u) 
  {
    matlab2eigen<double, Eigen::ColMajor>(array,u);
  }

  template <>
  mxArray* 
  MxArray::from(const Eigen::SparseMatrix<double, Eigen::RowMajor>& u) 
  {
    return eigen2matlab<double,Eigen::RowMajor>(u);
  }

  template <>
  void 
  MxArray::to(const mxArray* array, Eigen::SparseMatrix<double, Eigen::RowMajor>* u) 
  {
    matlab2eigen<double,Eigen::RowMajor>(array,u);
  }

  // Define template specializations for ints.
  template <>
  mxArray* 
  MxArray::from(const Eigen::SparseMatrix<int, Eigen::ColMajor>& u) 
  {
    return eigen2matlab<int, Eigen::ColMajor>(u);
  }

  template <>
  void 
  MxArray::to(const mxArray* array, Eigen::SparseMatrix<int, Eigen::ColMajor>* u) 
  {
    matlab2eigen<int, Eigen::ColMajor>(array,u);
  }

  template <>
  mxArray* 
  MxArray::from(const Eigen::SparseMatrix<int, Eigen::RowMajor>& u) 
  {
    return eigen2matlab<int,Eigen::RowMajor>(u);
  }

  template <>
  void 
  MxArray::to(const mxArray* array, Eigen::SparseMatrix<int, Eigen::RowMajor>* u) 
  {
    matlab2eigen<int,Eigen::RowMajor>(array,u);
  }

  // Define template specializations for long long.
  template <>
  mxArray* 
  MxArray::from(const Eigen::SparseMatrix<long long, Eigen::ColMajor>& u) 
  {
    return eigen2matlab<long long, Eigen::ColMajor>(u);
  }

  template <>
  void 
  MxArray::to(const mxArray* array, Eigen::SparseMatrix<long long, Eigen::ColMajor>* u) 
  {
    matlab2eigen<long long, Eigen::ColMajor>(array,u);
  }

  template <>
  mxArray* 
  MxArray::from(const Eigen::SparseMatrix<long long, Eigen::RowMajor>& u) 
  {
    return eigen2matlab<long long,Eigen::RowMajor>(u);
  }

  template <>
  void 
  MxArray::to(const mxArray* array, Eigen::SparseMatrix<long long, Eigen::RowMajor>* u) 
  {
    matlab2eigen<long long,Eigen::RowMajor>(array,u);
  }

} // namespace mexplus

#endif
