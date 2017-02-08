#ifndef BACKEND_ARMA_H
#define BACKEND_ARMA_H

#if defined(HAS_ARMA)

#include "commons/commons.h"
#include <cassert>

namespace SoftComp {



//template<typename MatType>
//SC_INLINE integer size(const MatType &mat) {
//    return mat.n_elem;
//}
//template<typename MatType>
//SC_INLINE integer size(MatType &&mat) {
//    return mat.n_nelem;
//}

//template<typename MatType>
//SC_INLINE integer size(const MatType &mat, integer idx) {
//    return idx==0 ? mat.n_rows : mat.n_cols;
//}
//template<typename MatType>
//SC_INLINE integer size(MatType &&mat, integer idx) {
//    return idx==0 ? mat.n_rows : mat.n_cols;
//}

template<typename MatType>
SC_INLINE integer rows(const MatType &mat) {
    return mat.n_rows;
}
template<typename MatType>
SC_INLINE integer rows(MatType &&mat) {
    return mat.n_rows;
}
template<typename MatType>
SC_INLINE integer cols(const MatType &mat) {
    return mat.n_cols;
}
template<typename MatType>
SC_INLINE integer cols(MatType &&mat) {
    return mat.n_cols;
}

// auto returns a arma::subview_col<T> which does not bind to arma::Col<T>

//template<typename T>
//SC_INLINE auto row(const arma::Mat<T> &mat, integer idx) -> decltype(mat.row(idx)) {
//    return mat.row(idx);
//}
//template<typename T>
//SC_INLINE auto row(arma::Mat<T> &&mat, integer idx) -> decltype(mat.row(idx)) {
//    return mat.row(idx);
//}

//template<typename T>
//SC_INLINE auto col(const arma::Mat<T> &mat, integer idx) -> decltype(mat.col(idx)) {
//    return mat.col(idx);
//}
//template<typename T>
//SC_INLINE auto col(arma::Mat<T> &&mat, integer idx) -> decltype(mat.col(idx)) {
//    return mat.col(idx);
//}
template<typename T>
SC_INLINE arma::Row<T> row(const arma::Mat<T> &mat, integer idx) {
    return mat.row(idx);
}
template<typename T>
SC_INLINE arma::Row<T> row(arma::Mat<T> &&mat, integer idx) {
    return mat.row(idx);
}

template<typename T>
SC_INLINE arma::Col<T> col(const arma::Mat<T> &mat, integer idx) {
    return mat.col(idx);
}
template<typename T>
SC_INLINE arma::Col<T> col(arma::Mat<T> &&mat, integer idx) {
    return mat.col(idx);
}

SC_INLINE matrix<real> rand(integer m, integer n) {
    matrix<real> out; out.randu(m,n);
    return out;
}

template<typename T>
SC_INLINE vector<T> zeros(integer m) {
    vector<T> out(m); out.zeros();
}
template<typename T>
SC_INLINE matrix<T> zeros(integer m, integer n) {
    matrix<T> out(m,n); out.zeros();
    return out;
}
SC_INLINE vector<real> zeros(integer m) {
    vector<real> out(m); out.zeros();
    return out;
}
SC_INLINE matrix<real> zeros(integer m, integer n) {
    matrix<real> out(m,n); out.zeros();
    return out;
}

template<typename T>
SC_INLINE vector<T> ones(integer m) {
    vector<T> out(m); out.ones();
    return out;
}
template<typename T>
SC_INLINE matrix<T> ones(integer m, integer n) {
    matrix<T> out(m,n); out.ones();
    return out;
}
SC_INLINE vector<real> ones(integer m) {
    vector<real> out(m); out.ones();
    return out;
}
SC_INLINE matrix<real> ones(integer m, integer n) {
    matrix<real> out(m,n); out.ones();
    return out;
}

SC_INLINE matrix<real> eye(integer m, integer n) {
    matrix<real> out(m,n); out.eye();
    return out;
}

template<typename T>
SC_INLINE void fill(matrix<T> &mat, integer num) {
    mat.fill(num);
}
template<typename T>
SC_INLINE vector<T> fill(integer cap, integer num) {
    vector<T> mat(cap);
    std::fill(mat.memptr(),mat.memptr()+mat.n_elem,num);
    return mat;
}

template<typename T>
SC_INLINE matrix<T> arange(integer m, integer n) {
    matrix<T> out(n-m,1);
    std::iota(out.memptr(),out.memptr()+n-m,m);
    return out;
}
template<typename T>
SC_INLINE matrix<T> arange(integer n) {
    matrix<T> out(n,1);
    std::iota(out.memptr(),out.memptr()+n,0);
    return out;
}

SC_INLINE vector<real> linspace(real low, real high, integer num=50)
{
    //! Linearly spaced points. Dont parametrise!
    return arma::linspace<vector<real>>(low,high,num);
}

SC_INLINE vector<real> logspace(real low, real high, integer num=50)
{
    //! Linearly spaced points. Dont parametrise!
    return arma::logspace<vector<real>>(low,high,num);
}


template<typename T>
SC_INLINE T max(const matrix<T> &mat) {
    return mat.max();
}
template<typename T>
SC_INLINE T min(const matrix<T> &mat) {
    return mat.min();
}
template<typename T>
SC_INLINE T argmax(const matrix<T> &mat) {
    return mat.index_max();
}
template<typename T>
SC_INLINE T argmin(const matrix<T> &mat) {
    return mat.index_min();
}

template<typename T>
SC_INLINE boolean isempty(const matrix<T> &mat) {
    return mat.n_elem==0;
}

template<typename T>
SC_INLINE boolean isempty(const vector<T> &mat) {
    return mat.n_rows==0;
}

template<typename T, typename U = T>
SC_INLINE vector<T>
append(const vector<T> &arr, U num) {
    //! Append to the end of Eigen vector (similar to push_back)

    vector<T> new_arr(rows(arr)+1); new_arr.zeros();
    new_arr(0,arma::span(0,rows(arr))) = arr;
    new_arr(rows(arr)) = static_cast<T>(num);
    return new_arr;
}

template<typename T>
SC_INLINE matrix<T> reshape(const matrix<T> &mat, integer m, integer n) {
    return mat.reshape(m,n);
}
template<typename T>
SC_INLINE matrix<T> reshape(const vector<T> &mat, integer m, integer n) {
    SC_ASSERT(m*n==size(mat),"CANNOT CHANGE SIZE ON RESHAPE");
    matrix<T> out(m,n);
    std::copy(out.memptr(),out.memptr()+size(out),mat.memptr());
    return out;
}

template<typename T>
SC_INLINE matrix<T>
ravel(const matrix<T> &arr) {
    //! Ravel/flatten a matrix. Makes a copy
    return reshape(arr,size(arr),1);
}


//template<typename T, typename U = T>
//std::tuple<matrix<integer>,matrix<integer> >
//SC_INLINE find(const Eigen::PlainObjectBase<T> &arr,
//         U num, real tolerance=1e-14)
//{
//    //! find the occurence of a value in a matrix
//    std::vector<integer> idx_rows;
//    std::vector<integer> idx_cols;
//    idx_rows.clear(); idx_cols.clear();
//    for (auto i=0; i<arr.rows();++i)
//    {
//        for (auto j=0; j<arr.cols();++j)
//        {
//            if (static_cast<real>(abs(arr(i,j)-num))<tolerance)
//            {
//                idx_rows.push_back(i);
//                idx_cols.push_back(j);
//            }
//        }
//    }

//    return std::make_tuple(
//                Eigen::Map<matrix<integer>>
//                (idx_rows.data(),idx_rows.size(),1),
//                Eigen::Map<matrix<integer>>
//                (idx_cols.data(),idx_cols.size(),1));
//}


template<typename T>
SC_INLINE matrix<T> fliplr(const matrix<T> &mat) {
    //! Reverse a matrix left to right
    SC_ASSERT(mat.n_cols!=1,"CANNOT FLIP A COLUMN VECTOR ROW-WISE");
    return arma::fliplr(mat);
}

template<typename T>
SC_INLINE matrix<T> flipud(const matrix<T> &mat) {
    //! Reverse a matrix left to right
    SC_ASSERT(mat.n_rows!=1,"CANNOT FLIP A COLUMN VECTOR ROW-WISE");
    return arma::flipud(mat);
}

template<typename T>
SC_INLINE matrix<T>
copy(const matrix<T> &arr) {
    return out(arr);
}

template<typename Derived>
SC_INLINE auto sort(const Derived &mat, integer axis=0) -> decltype(arma::sort(mat,"ascend",axis)) {
    //! Sort a matrix
    return arma::sort(mat,"ascend",axis);
}

template<typename Derived>
SC_INLINE auto argsort(const Derived &mat, integer axis=0) -> decltype(arma::sort(mat,"ascend",axis)) {
    //! Sort a matrix
    return arma::sort_index(mat,"ascend");
}

template<template<typename> class Mat, typename T>
SC_INLINE matrix<T> unique(const Mat<T> &arr, integer axis=-1,
                           boolean keep_order=true, boolean consider_sort=false) {
    //! Find unique values of a matrix as a whole, row-wise or column-wise
    //!
    //! input:
    //!     arr:                array whose unique values need to be found
    //!     axis:               axis along which to take unique values:
    //!                             0 - row-wise
    //!                             1 - column-wise
    //!                            -1 - whole array (default)
    //!
    //!     keep_order:         keep the original order of the array (default false). If
    //!                         true, returns a sorted array
    //!     consider_sort:      only holds for taking uniques along an axis, in that if permutation of
    //!                         array within an axis matters that is two rows/columns can be equal to
    //!                         each other but with different arrangement of elements. If consider_sort
    //!                         is true, those rows/columns would be treated equal

    if (axis==-1) return arma::unique(arr);
    else {
        matrix<T> arr_ = arr; if (axis==1) arr_ = arma::trans(arr);

        std::vector<std::vector<T>> input(rows(arr_));
        for(auto i=0; i<rows(arr_); ++i) {
            input[i] = std::vector<T>(cols(arr_));
            for (auto j=0; j<cols(arr_); ++j) {
                input[i][j] = arr_(i,j);
            }
        }

        if (consider_sort) {
            for (auto &k: input)
                std::sort(k.begin(),k.end());
        }

        if (keep_order) unsort_unique(input);
        else {
            std::sort(input.begin(), input.end());
            input.erase(std::unique(input.begin(), input.end()), input.end());
        }


        matrix<T> out(input.size(),cols(arr_));
        for(size_t i=0; i<input.size(); ++i) {
            for (size_t j=0; j<input[i].size(); ++j) {
                out(i,j) = input[i][j];
            }
        }

        if (axis==1) arma::inplace_trans(out);

        return out;
    }
}


template<typename T>
SC_INLINE matrix<T>
hstack(const matrix<T> &mat0, const matrix<T> &mat1) {
    SC_ASSERT(rows(mat0)==rows(mat1),"CONCATENATING MATRICES SHOULD HAVE THE SAME SIZE NUMBER OF ROWS");
    return arma::join_horiz(mat0,mat1);
}

//template<template<typename> class Mat0, template<typename> class Mat1, typename T>
//SC_INLINE matrix<T>
//vstack(const Mat0<T> &mat0, const Mat1<T> &mat1) {
//    SC_ASSERT(cols(mat0)==cols(mat1),"CONCATENATING MATRICES SHOULD HAVE THE SAME SIZE NUMBER OF COLUMNS");
//    return arma::join_vert(mat0,mat1);
//}

template<typename Derived0, typename Derived1>
SC_INLINE auto
vstack(const Derived0 &mat0, const Derived1 &mat1) -> decltype(arma::join_vert(mat0,mat1)) {
    return arma::join_vert(mat0,mat1);
}

template<template<typename> class Mat, typename T>
SC_INLINE Mat<T>
repmat(const Mat<T> &mat, integer a, integer b=None) {
    //! Equivalent to MATLAB's repmat
    b = b==None ? a : b;
    return arma::repmat(mat,a,b);
}




// Math
template<template<typename> class Mat, typename T>
SC_INLINE matrix<real> abs(const Mat<T> &mat) {
    return arma::abs(mat);
}
template<template<typename> class Mat, typename T, typename U>
SC_INLINE matrix<real> pow(const Mat<T> &mat, U num) {
    return arma::pow(mat,num);
}
//template<typename Derived, typename T>
//SC_INLINE matrix<T> pow(T num, const Eigen::MatrixBase<Derived> &mat) {
//    matrix<T> out = Eigen::pow(num,mat.array());
//    return out;
//}
//template<typename Derived0, typename Derived1>
//SC_INLINE matrix<real> pow(const Eigen::MatrixBase<Derived0> &mat, const Eigen::MatrixBase<Derived1> &nums) {
//    return mat.array().pow(nums.array());
//}
//template<typename Derived>
//SC_INLINE matrix<real> log(const Eigen::MatrixBase<Derived> &mat) {
//    return mat.array().log();
//}
//template<typename Derived>
//SC_INLINE matrix<real> log10(const Eigen::MatrixBase<Derived> &mat) {
//    return mat.array().log10();
//}
//template<typename Derived>
//SC_INLINE matrix<real> log1p(const Eigen::MatrixBase<Derived> &mat) {
//    return mat.array().log1p();
//}
//template<typename Derived>
//SC_INLINE matrix<real> exp(const Eigen::MatrixBase<Derived> &mat) {
//    return mat.array().exp();
//}

//// Trig
//template<typename Derived>
//SC_INLINE matrix<real> sin(const Eigen::MatrixBase<Derived> &mat) {
//    return mat.array().sin();
//}
//template<typename Derived>
//SC_INLINE matrix<real> cos(const Eigen::MatrixBase<Derived> &mat) {
//    return mat.array().cos();
//}
//template<typename Derived>
//SC_INLINE matrix<real> tan(const Eigen::MatrixBase<Derived> &mat) {
//    return mat.array().tan();
//}
//template<typename Derived>
//SC_INLINE matrix<real> sinh(const Eigen::MatrixBase<Derived> &mat) {
//    return mat.array().sinh();
//}
//template<typename Derived>
//SC_INLINE matrix<real> cosh(const Eigen::MatrixBase<Derived> &mat) {
//    return mat.array().cosh();
//}
//template<typename Derived>
//SC_INLINE matrix<real> tanh(const Eigen::MatrixBase<Derived> &mat) {
//    return mat.array().tanh();
//}




// Linear algebra subroutines
//template<template<typename> class Mat0, template<typename> class Mat1, typename T>
//SC_INLINE auto
//multiply(const Mat0<T> &mat0, const Mat1<T> &mat1) -> decltype(mat0 %mat1) {
//    return mat0 % mat1;
//}

template<typename Derived0, typename Derived1>
SC_INLINE auto
multiply(const Derived0 &mat0, const Derived1 &mat1) -> decltype(mat0 % mat1) {
    return mat0 % mat1;
}

//template<typename Derived>
//SC_INLINE matrix<real> transpose(const Eigen::PlainObjectBase<Derived> &mat) {
//    return mat.transpose();
//}

//template<typename Derived>
//SC_INLINE matrix<real> inv(const Eigen::PlainObjectBase<Derived> &mat) {
//    return mat.inverse();
//}

//template<typename Derived>
//SC_INLINE real det(const Eigen::PlainObjectBase<Derived> &mat) {
//    return mat.determinant();
//}

//template<typename Derived>
//SC_INLINE real norm(const Eigen::MatrixBase<Derived> &mat) {
//    return mat.squaredNorm();
//}

template<typename T, typename Derived1>
SC_INLINE matrix<T> solve(const matrix<T> &mat, const vector<T> &vec, integer info=linear_solver::non_symmetric) {
    info++;
    return solve(mat,vec);
}

}

#endif


#endif // BACKEND_ARMA_H
