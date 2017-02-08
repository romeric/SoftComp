#ifndef COMMONS_H
#define COMMONS_H

#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <algorithm>
#include <functional>
#include <numeric>
//#include <tuple>
#include <type_traits>
#include <utility>
#include <cmath>
#include <cstdlib>


#ifdef __SSE2__
#include <emmintrin.h>
#endif
#ifdef __AVX__
#include <immintrin.h>
#endif



#ifdef __GNUC__
    #ifndef __clang__
        #ifndef __INTEL_COMPILER
            #define FASTOR_GCC
        #endif
    #endif
#endif

#ifdef __INTEL_COMPILER
    #define SC_INTEL
#endif

#ifdef __clang__
    #define SC_CLANG
#endif

#if defined(_MSC_VER)
    #define SC_MSVC
#endif

#if defined(_MSC_VER)
    #if _MSC_VER < 1800
       #error SOFTCOMP REQUIRES AN ISO C++11 COMPLIANT COMPILER
    #endif
#elif defined(__GNUC__) || defined(__GNUG__)
    #if __cplusplus <= 199711L
        #error SOFTCOMP REQUIRES AN ISO C++11 COMPLIANT COMPILER
    #endif
#endif

#if defined(__GNUC__) || defined(__GNUG__)
    #define SC_INLINE inline __attribute__((always_inline))
    #define SC_NOINLINE __attribute__((noinline))
#elif defined(_MSC_VER)
    #define FASTOR_INLINE __forceinline
    #define FASTOR_NOINLINE __declspec(noinline)
#endif

#if defined(__GNUC__) || defined(__GNUG__)
    #define SC_ALIGN __attribute__((aligned(0x20))) __restrict__
#elif defined(_MSC_VER)
    #define SC_ALIGN __declspec(align(32))
#endif


//#ifndef HAS_EIGEN
//#ifndef HAS_ARMA
//#ifndef HAS_BLAZE
//#define HAS_EIGEN
//#endif
//#endif
//#endif


//#define HAS_EIGEN
//#define HAS_ARMA
//#define HAS_BLAZE

#if defined(HAS_EIGEN)

#define EIGEN_VECTORIZE
#define EIGEN_DEFAULT_TO_ROW_MAJOR
#define EIGEN_HAVE_RVALUE_REFERENCES

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#elif defined(HAS_BLAZE)

#include <blaze/Blaze.h>
#include <blaze/Math.h>

#elif defined(HAS_ARMA)

//#define ARMA_DONT_USE_WRAPPER
#include <armadillo>

#endif

// Define these after Eigen/Blaze includes
using real    = double;
using integer = long long int;
using std::size_t;
using character = char;
using std::string;
using boolean = bool;

// macros
#define None 0
constexpr real EPS = std::numeric_limits<real>::epsilon();
//constexpr real NAN = std::numeric_limits<real>::quiet_NaN();
#define PI           3.14159265358979323846  /* pi */



#if defined(HAS_EIGEN)

template<typename T, size_t M, size_t N>
using smatrix = Eigen::Matrix<T,M,N>;

//template<typename T>
//using matrix = Eigen::Matrix<T,-1,-1,Eigen::ColMajor,-1,-1>;

template<typename T>
using matrix = Eigen::Matrix<T,-1,-1,Eigen::RowMajor,-1,-1>;

template<typename T, size_t M>
using svector = Eigen::Matrix<T,M,1>;

template<typename T>
using vector = Eigen::Matrix<T,-1,1>;

// Un-comment once sparse matrices are used
template<typename T>
using spmatrix = Eigen::SparseMatrix<T,Eigen::RowMajor>;

#elif defined(HAS_BLAZE)

template<typename T, size_t M, size_t N>
using smatrix = blaze::StaticMatrix<T,M,N,blaze::rowMajor>;

template<typename T>
using matrix = blaze::DynamicMatrix<T,blaze::rowMajor>;

template<typename T, size_t M>
using svector = blaze::StaticVector<T,M>;

template<typename T>
using vector = blaze::DynamicVector<T>;

#elif defined(HAS_ARMA)

//template<typename T, size_t M, size_t N>
//using smatrix = arma::Mat<T>::fixed<M,N>;

template<typename T>
using matrix = arma::Mat<T>;

template<typename T>
using vector = arma::Col<T>;

#endif


void SC_ASSERT(boolean cond, const string &x) {
    if (cond==true) {
        return;
    }
    else {
        std::cout << x << "\n";
        exit(EXIT_FAILURE);
    }
}

#define SC_STATIC_ASSERT(cond, msg) static_assert(cond,msg);




//#include "utils.h"



#endif // COMMONS_H
