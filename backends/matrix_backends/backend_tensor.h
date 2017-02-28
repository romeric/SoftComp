#ifndef BACKEND_TENSOR_H
#define BACKEND_TENSOR_H

#include "commons/commons.h"

#if defined(HAS_EIGEN)
#include <unsupported/Eigen/CXX11/Tensor>
template<typename T, integer rank>
using tensor = Eigen::Tensor<T, rank, Eigen::RowMajor>;
#endif

namespace SoftComp {

using tpair = Eigen::IndexPair<integer>;

// template<typename Derived, int AccessLevel>
// SC_INLINE integer size(const Eigen::TensorBase<Derived, AccessLevel> &tens, integer i) {
//    return tens.dimension(i);
// }
template<typename T, int rank>
SC_INLINE integer size(const Eigen::Tensor<T,rank,Eigen::RowMajor> &tens) {
   integer si = 1;
   const auto &dims = tens.dimensions(); 
   for (auto &k: dims) si *=k;
   return si;
}

template<typename T, int rank>
SC_INLINE integer size(const Eigen::Tensor<T,rank,Eigen::RowMajor> &tens, integer i) {
   return tens.dimension(i);
}

template<typename T>
SC_INLINE T tensor2scalar(const Eigen::Tensor<T,0,Eigen::RowMajor> &tens) {
   return tens.data()[0];
}

template<typename T>
SC_INLINE boolean isempty(const tensor<T,1> &a) {
   return a.dimension(0)==0 ? true : false;
}
template<typename T>
SC_INLINE boolean isempty(const tensor<T,2> &a) {
   return (a.dimension(0)==0 && a.dimension(1)==0) ? true : false;
}
template<typename T>
SC_INLINE boolean isempty(const tensor<T,3> &a) {
   return (a.dimension(0)==0 && a.dimension(1)==0 && a.dimension(2)==0) ? true : false;
}
template<typename T>
SC_INLINE boolean isempty(const tensor<T,4> &a) {
   return (a.dimension(0)==0 && a.dimension(1)==0 && 
    a.dimension(2)==0 && a.dimension(3)==0) ? true : false;
}
template<typename T>
SC_INLINE boolean isempty(const tensor<T,5> &a) {
   return (a.dimension(0)==0 && a.dimension(1)==0 && 
    a.dimension(2)==0 && a.dimension(3)==0 && a.dimension(4)==0) ? true : false;
}


// MATH
template<typename T, int rank>
SC_INLINE auto abs(const tensor<T,rank> a) -> decltype(a.abs()) {
    return a.abs();
}

template<typename T, int rank>
SC_INLINE auto log(const tensor<T,rank> a) -> decltype(a.log()) {
    return a.log();
}

template<typename T, int rank>
SC_INLINE auto log1p(const tensor<T,rank> a) -> decltype(a.log1p()) {
    return a.log1p();
}

template<typename T, int rank>
SC_INLINE auto exp(const tensor<T,rank> a) -> decltype(a.exp()) {
    return a.exp();
}

template<typename T, int rank>
SC_INLINE auto sum(const tensor<T,rank> a) -> decltype(a.sum()) {
    return a.sum();
}
template<typename T, int rank>
SC_INLINE tensor<T,rank-1> sum(const tensor<T,rank> a, integer axis) {
    std::array<integer,1> idx = {axis};
    return a.sum(idx);
}

template<typename T, int rank>
SC_INLINE auto mean(const tensor<T,rank> a) -> decltype(a.mean()) {
    return a.mean();
}
template<typename T, int rank>
SC_INLINE tensor<T,rank-1> mean(const tensor<T,rank> a, integer axis) {
    std::array<integer,1> idx = {axis};
    return a.mean(idx);
}

template<typename T, int rank>
SC_INLINE auto max(const tensor<T,rank> a) -> decltype(a.maximum()) {
    return a.maximum();
}
template<typename T, int rank>
SC_INLINE tensor<T,rank-1> max(const tensor<T,rank> a, integer axis) {
    std::array<integer,1> idx = {axis};
    return a.maximum(idx);
}

template<typename T, int rank>
SC_INLINE auto min(const tensor<T,rank> a) -> decltype(a.minimum()) {
    return a.minimum();
}
template<typename T, int rank>
SC_INLINE tensor<T,rank-1> min(const tensor<T,rank> a, integer axis) {
    std::array<integer,1> idx = {axis};
    return a.minimum(idx);
}

// template<typename T, int rank>
// SC_INLINE real l2norm(const tensor<T,rank> &a) {
    // integer ncomponents = size(a.data()); // this is wrong
    // real norm           =  0;
    // for (auto i=0; i<ncomponents; i++){
    //     norm            +=  a.data()[i]*a.data()[i];
    // }
    // norm                =  sqrt(norm);
    // return norm;
// }

template<typename T, int rank>
SC_INLINE T norm(const tensor<T,rank> &a) {
    return sqrt(std::inner_product(a.data(),a.data()+size(a),a.data(),0));
}

// template<typename T>
// SC_INLINE real l2norm(const vector<T> a) {
//     integer ncomponents = size(a);
//     real norm           =  0;
//     for (auto i=0; i<ncomponents; i++){
//         norm            +=  a.data()[i]*a.data()[i];
//     }
//     norm                =  sqrt(norm);
//     return norm;
// }

// END

//
template<typename T>
SC_INLINE matrix<T> tomatrix(const tensor<T,2> a) {
    matrix<T> out(size(a,0),size(a,1)); 
    std::copy(a.data(),a.data()+size(a),out.data());
    return out;
}
//

// Einsum functionality
//------------------------------------------------------------------------------
template<int ... rest>
struct CIndex {
    static constexpr integer no_indices = sizeof...(rest);
};

template<class Ind0, class Ind1>
struct einsum_extractor {};

template<int ... idx0, int ... idx1>
struct einsum_extractor<CIndex<idx0...>,CIndex<idx1...>> {

template<typename T, int rank0, int rank1,
        typename std::enable_if<sizeof...(idx0)==sizeof...(idx1),bool>::type=0>
static SC_INLINE tensor<T,rank0+rank1-(int)sizeof...(idx0)-(int)sizeof...(idx0)> 
einsum_impl(const tensor<T,rank0> &a, const tensor<T,rank1> &b) {

    std::array<Eigen::IndexPair<int>,sizeof...(idx0)> idx;
    constexpr int arr0[sizeof...(idx0)] = {idx0...};
    constexpr int arr1[sizeof...(idx1)] = {idx1...};

    for (auto i=0; i<sizeof...(idx0); ++i) {
        idx[i] = Eigen::IndexPair<int>(arr0[i],arr1[i]);
    }
    return a.contract(b,idx);
}
};
template<class Ind0, class Ind1, typename T, int rank0, int rank1>
SC_INLINE auto einsum(const tensor<T,rank0> &a, const tensor<T,rank1> &b) 
-> decltype(einsum_extractor<Ind0,Ind1>::einsum_impl(a,b)) {
    return einsum_extractor<Ind0,Ind1>::einsum_impl(a,b);
}
//------------------------------------------------------------------------------


template<typename T>
void print(const Eigen::Tensor<T,3,Eigen::RowMajor> &a) {
    std::string s;
    for (auto i=0; i<a.dimension(0); ++i) {
        Eigen::Tensor<T,2,Eigen::RowMajor> a_chip = a.chip(i,0);
        s = "["+std::to_string(i)+",:,:]";
        print(s,a_chip); print();
    }
}
template<typename T>
void print(const Eigen::Tensor<T,4,Eigen::RowMajor> &a) {
    std::string s;
    for (auto i=0; i<a.dimension(0); ++i) {
        Eigen::Tensor<T,3,Eigen::RowMajor> a_chip3 = a.chip(i,0);
        for (auto j=0; j<a.dimension(1); ++j) {
            Eigen::Tensor<T,2,Eigen::RowMajor> a_chip2 = a_chip3.chip(j,0);
            s = "["+std::to_string(i)+","+std::to_string(j)+",:,:]";
            print(s,a_chip2); print();
        }
    }
}
template<typename T>
void print(const Eigen::Tensor<T,5,Eigen::RowMajor> &a) {
    std::string s;
    for (auto i=0; i<a.dimension(0); ++i) {
        Eigen::Tensor<T,4,Eigen::RowMajor> a_chip4 = a.chip(i,0);
        for (auto j=0; j<a.dimension(1); ++j) {
            Eigen::Tensor<T,3,Eigen::RowMajor> a_chip3 = a_chip4.chip(j,0);
            for (auto k=0; k<a.dimension(2); ++k) {
                Eigen::Tensor<T,2,Eigen::RowMajor> a_chip2 = a_chip3.chip(k,0);
                s = "["+std::to_string(i)+","+std::to_string(j)+","+std::to_string(k)+",:,:]";
                print(s,a_chip2); print();
            }
        }
    }
}
template<typename T>
void print(const Eigen::Tensor<T,6,Eigen::RowMajor> &a) {
    std::string s;
    for (auto i=0; i<a.dimension(0); ++i) {
        Eigen::Tensor<T,5,Eigen::RowMajor> a_chip5 = a.chip(i,0);
        for (auto j=0; j<a.dimension(1); ++j) {
            Eigen::Tensor<T,4,Eigen::RowMajor> a_chip4 = a_chip5.chip(j,0);
            for (auto k=0; k<a.dimension(2); ++k) {
                Eigen::Tensor<T,3,Eigen::RowMajor> a_chip3 = a_chip4.chip(k,0);
                for (auto l=0; l<a.dimension(3); ++l) {
                    Eigen::Tensor<T,2,Eigen::RowMajor> a_chip2 = a_chip3.chip(k,0);
                    s = "["+std::to_string(i)+","+std::to_string(j)+","+std::to_string(k)+","+std::to_string(l)+",:,:]";
                    print(s,a_chip2); print();
                }
            }
        }
    }
}


// template<int ... all>
// struct TPair {};


// template<typename T, int rank0, int rank1, TPair<int> ... pairs>
// contract(tensor<T,rank0> &a, rank<T,rank1> &b) {
//     std::array<Eigen::IndexPair<int>,sizeof...(tpair)> idx;
//     return a.contract(b,idx);
// }
 // contract<TPair<0,0>,TPair<1,1>>(a,b);


SC_INLINE static tensor<real,2> identity2d(integer ndim) {
    //--------------------------------------------------------------
    // 2D and 3D identity matrix for tensor (only stored as 2nd
    // order tensors) (Why called static?Roman)
    //--------------------------------------------------------------
    tensor<real, 2> tens(ndim,ndim);
    tens.setZero();
    for (integer i=0; i<ndim; i++){
        tens(i,i)  =  1;
    }
    return tens;
}

SC_INLINE static tensor<real,3> levi_civita(integer ndim) {
    //--------------------------------------------------------------
    // Levi Civita third order tensor
    //--------------------------------------------------------------
    tensor<real, 3> epsilon(ndim,ndim,ndim);
    epsilon.setZero();
    epsilon(0,1,2)   =   1;
    epsilon(1,2,0)   =   1;
    epsilon(2,0,1)   =   1;
    epsilon(1,0,2)   =  -1;
    epsilon(2,1,0)   =  -1;
    epsilon(0,2,1)   =  -1;
    return epsilon;
}

SC_INLINE static tensor<real,4> identity4d(integer ndim) {
    //--------------------------------------------------------------
    // Fourth order identity tensor A(i,I,j,J)=delta(i,j)*delta(I,J)
    //--------------------------------------------------------------
    matrix<real> identity_2d = eye(ndim,ndim);
    tensor<real, 4> identity_4d(ndim,ndim,ndim,ndim);
    for (integer i=0; i<ndim; ++i) {
        for (integer I=0; I<ndim; ++I) {
            for (integer j=0; j<ndim; ++j) {
                for (integer J=0; J<ndim; ++J) {
                    identity_4d(i,I,j,J) = identity_2d(i,j)*identity_2d(I,J);
                }
            }
        }
    }
    return identity_4d;
}

template<typename T>
SC_INLINE void __cross__(const T *__restrict__ A_data, const T *__restrict__ B_data, T *__restrict__ out_data) {
    //--------------------------------------------------------------
    // Tensor contraction function (only for 3D)(*Roman)
    //--------------------------------------------------------------
       T A00        =  A_data[0];
       T A11        =  A_data[4];
       T A22        =  A_data[8];
       T A01        =  A_data[1];
       T A02        =  A_data[2];
       T A12        =  A_data[5];
       T A10        =  A_data[3];
       T A20        =  A_data[6];
       T A21        =  A_data[7];

       T B00        =  B_data[0];
       T B11        =  B_data[4];
       T B22        =  B_data[8];
       T B01        =  B_data[1];
       T B02        =  B_data[2];
       T B12        =  B_data[5];
       T B10        =  B_data[3];
       T B20        =  B_data[6];
       T B21        =  B_data[7];

       out_data[0]  =  A11*B22 - A12*B21 - A21*B12 + A22*B11;
       out_data[1]  =  A12*B20 - A10*B22 + A20*B12 - A22*B10;
       out_data[2]  =  A10*B21 - A11*B20 - A20*B11 + A21*B10;
       out_data[3]  =  A02*B21 - A01*B22 + A21*B02 - A22*B01;
       out_data[4]  =  A00*B22 - A02*B20 - A20*B02 + A22*B00;
       out_data[5]  =  A01*B20 - A00*B21 + A20*B01 - A21*B00;
       out_data[6]  =  A01*B12 - A02*B11 - A11*B02 + A12*B01;
       out_data[7]  =  A02*B10 - A00*B12 + A10*B02 - A12*B00;
       out_data[8]  =  A00*B11 - A01*B10 - A10*B01 + A11*B00;
}

SC_INLINE void __inv__(const real *__restrict__ A_data, real *__restrict__ out, int ndim) {
    //--------------------------------------------------------------
    // Inverse of a matrix (3D and 2D)
    //--------------------------------------------------------------
    if (ndim==3) {
        real A00      =  A_data[0];
        real A11      =  A_data[4];
        real A22      =  A_data[8];
        real A01      =  A_data[1];
        real A02      =  A_data[2];
        real A12      =  A_data[5];
        real A10      =  A_data[3];
        real A20      =  A_data[6];
        real A21      =  A_data[7];

        real detA     =  A_data[0]*A_data[4]*A_data[8] + A_data[7]*A_data[3]*A_data[2] + A_data[6]*A_data[1]*A_data[5] \
                        -A_data[6]*A_data[4]*A_data[2] - A_data[7]*A_data[5]*A_data[0] - A_data[8]*A_data[1]*A_data[3] ;
        real adj00   =  A11*A22 - A21*A12;
        real adj01   =  A02*A21 - A22*A01;
        real adj02   =  A01*A12 - A11*A02;
        real adj10   =  A12*A20 - A22*A10;
        real adj11   =  A00*A22 - A20*A02;
        real adj12   =  A02*A10 - A12*A00;
        real adj20   =  A10*A21 - A20*A11;
        real adj21   =  A01*A20 - A21*A00;
        real adj22   =  A00*A11 - A10*A01;

        out[0]       =  adj00/detA;
        out[1]       =  adj01/detA;
        out[2]       =  adj02/detA;
        out[3]       =  adj10/detA;
        out[4]       =  adj11/detA;
        out[5]       =  adj12/detA;
        out[6]       =  adj20/detA;
        out[7]       =  adj21/detA;
        out[8]       =  adj22/detA;
    }
    else {
        real A00      =  A_data[0];
        real A01      =  A_data[1];
        real A10      =  A_data[2];
        real A11      =  A_data[3];

        real detA     =  A00*A11 - A01*A10;
        out[0]        =  A11/detA;
        out[1]        = -A01/detA;
        out[2]        = -A10/detA;
        out[3]        =  A00/detA;
    }
}


template<typename T, int rank0, int rank1>
tensor<T,rank0+rank1> outer(tensor<T,rank0> a, const tensor<T,rank1> &b) {
    //--------------------------------------------------------------
    // Concatenate out product
    //--------------------------------------------------------------
    const T *__restrict__ a_data = a.data();
    const T *__restrict__ b_data = b.data();
    std::array<int,rank0+rank1> dims;
    int i=0;
    for (; i<rank0; ++i) {
        dims[i] = a.dimension(i);
    }
    for (; i<rank0+rank1; ++i) {
        dims[i] = b.dimension(i);
    }

    // out.resize(dims);
    tensor<T,rank0+rank1> out;
    T *__restrict__ out_data = out.data();


    for (auto i=0; i<a.size(); ++i) {
        T a0 = a_data[i];
        for (auto j=0; j<b.size(); ++j) {
          out_data[j+i*b.size()] = a0*b_data[j];
        }
    }

    return out;
}

tensor<real,2> cross(const tensor<real,2> &a, const tensor<real,2> &b) {
    //--------------------------------------------------------------
    // tensor cross product of two rank 2 tensors
    //--------------------------------------------------------------
    tensor<real,2> out(size(b,0),size(b,1));
    __cross__(a.data(),b.data(),out.data());
    return out;
}


tensor<real,3> cross(const tensor<real,3> &a, const tensor<real,3> &b) {
    //--------------------------------------------------------------
    // Concatenate tensor contraction function for every Gauss point(*Roman)
    //--------------------------------------------------------------
    auto n3           = b.dimension(0);
    tensor<real,3> out(b.dimension(0),b.dimension(1),b.dimension(2));
    for (auto i=0; i<n3; ++i) {
        tensor<real,2> a_2d = a.chip(i,0);
        tensor<real,2> b_2d = b.chip(i,0);
        __cross__(a_2d.data(), b_2d.data(),out.data()+a_2d.dimension(0)*a_2d.dimension(1)*i);
    }
    return out;
}

real det(const tensor<real,2> &a) {
    //--------------------------------------------------------------
    // determinant of a second order tensor
    //--------------------------------------------------------------
    const real *__restrict__ a_data = a.data();
    auto ndim    =  size(a,0);
    real out = 0;
    if (ndim==3){
       out =  a_data[0]*a_data[4]*a_data[8] + a_data[7]*a_data[3]*a_data[2] + a_data[6]*a_data[1]*a_data[5]
              -a_data[6]*a_data[4]*a_data[2] - a_data[7]*a_data[5]*a_data[0] - a_data[8]*a_data[1]*a_data[3] ;
    }
    else{
        out = a_data[0]*a_data[3] - a_data[1]*a_data[2];
    }
    return out;
}

tensor<real,2> inv(const tensor<real,2> &a) {
    //--------------------------------------------------------------
    // inverse of a second order tensor
    //--------------------------------------------------------------
    tensor<real,2> out(a.dimension(0),a.dimension(1));
    __inv__(a.data(),out.data(),size(a,0)); 
    return out;
}
tensor<real,3> inv(const tensor<real,3> &a) {
    //--------------------------------------------------------------
    // overload for third order tensors (computes the inverse in the xy plane)
    //--------------------------------------------------------------
    auto n3   =  a.dimension(0);
    tensor<real,3> out(a.dimension(0),a.dimension(1),a.dimension(2));
    for (auto i=0; i<n3; ++i) {
        tensor<real,2> a_2d = a.chip(i,0);
        __inv__(a_2d.data(),out.data()+size(a_2d,0)*size(a_2d,1)*i,n3);
    }
    return out;
}

tensor<real,2> cofactor(const tensor<real,2> &a) {
    //--------------------------------------------------------------
    // Concatenate Cofactor for tensors
    //--------------------------------------------------------------
    tensor<real,2>  out;
    auto ndim    =  a.dimension(1);
    if (ndim==3){
        out =  cross(a,a)*0.5;
    }
    else {
        out = tensor<real,2>(2,2);
        real *__restrict__ out_data = out.data();
        const real *__restrict__ a_data = a.data();
        out_data[0]  =   a_data[3];
        out_data[1]  =  -a_data[2];
        out_data[2]  =  -a_data[1];
        out_data[3]  =   a_data[0];
    }
    return out;
}


tensor<real,3> cofactor(const tensor<real,3> &a) {
    //--------------------------------------------------------------
    // Concatenate Cofactor for tensors
    //--------------------------------------------------------------
    tensor<real,3>  out;
    auto ndim    =  a.dimension(1);
    if (ndim==3){
        out =  cross(a,a)*0.5;
    }
    else {
        out = tensor<real,3>(2,2,size(a,0));
        real *__restrict__ out_data = out.data();
        const real *__restrict__ a_data = a.data();
        for (auto i=0; i<size(a,0); ++i) {
            out_data[i*4]   =   a_data[i*4+3];
            out_data[i*4+1] =  -a_data[i*4+2];
            out_data[i*4+2] =  -a_data[i*4+1];
            out_data[i*4+3] =   a_data[i*4];
        }
    }
    return out;
}


tensor<real,5> chiped_outer(const tensor<real,3> &a, const tensor<real,3> &b) {
    //--------------------------------------------------------------
    // Chip outer
    //--------------------------------------------------------------
    auto n3 = a.dimension(0);
    auto ndim = a.dimension(1);
    tensor<real,5> out(ndim,ndim,ndim,ndim,n3);
    for (integer i=0; i<n3; ++i) {
        for (auto j=0; j<ndim; ++j) {
            for (auto k=0; k<ndim; ++k) {
                for (auto l=0; l<ndim; ++l) {
                    for (auto m=0; m<ndim; ++m) {
                        out(j,l,k,m,i) = a(j,k,i)*b(l,m,i);
                    }
                }
            }
        }
    }
    return out;
}





}

#endif // BACKEND_TENSOR_H
