#ifndef BACKEND_EIGEN_H
#define BACKEND_EIGEN_H

#if defined(HAS_EIGEN)

#include "commons/commons.h"


namespace SoftComp {


template<typename Derived>
SC_INLINE integer size(const Eigen::PlainObjectBase<Derived> &mat) {
    return mat.rows()*mat.cols();
}
template<typename Derived>
SC_INLINE integer size(Eigen::PlainObjectBase<Derived> &&mat) {
    return mat.rows()*mat.cols();
}

template<typename Derived>
SC_INLINE integer size(const Eigen::PlainObjectBase<Derived> &mat, integer idx) {
    return idx==0 ? mat.rows() : mat.cols();
}
template<typename Derived>
SC_INLINE integer size(Eigen::PlainObjectBase<Derived> &&mat, integer idx) {
    return idx==0 ? mat.rows() : mat.cols();
}

template<typename Derived>
SC_INLINE integer rows(const Eigen::PlainObjectBase<Derived> &mat) {
    return mat.rows();
}
template<typename Derived>
SC_INLINE integer rows(Eigen::PlainObjectBase<Derived> &&mat) {
    return mat.rows();
}

template<typename Derived>
SC_INLINE integer cols(const Eigen::PlainObjectBase<Derived> &mat) {
    return mat.cols();
}
template<typename Derived>
SC_INLINE integer cols(Eigen::PlainObjectBase<Derived> &&mat) {
    return mat.cols();
}

#ifdef USE_EIGEN_EXPR
template<typename Derived>
SC_INLINE auto row(const Derived &mat, integer idx) -> decltype(mat.row(idx)) {
    return mat.row(idx);
}
template<typename Derived>
SC_INLINE auto row(Derived &&mat, integer idx) -> decltype(mat.row(idx)) {
    return mat.row(idx);
}

template<typename Derived>
SC_INLINE auto col(const Derived &mat, integer idx) -> decltype(mat.col(idx)) {
    return mat.col(idx);
}
template<typename Derived>
SC_INLINE auto col(Derived &&mat, integer idx) -> decltype(mat.col(idx)) {
    return mat.col(idx);
}
#else
template<typename Derived>
SC_INLINE vector<typename Derived::Scalar> row(const Derived &mat, integer idx) {
    return mat.row(idx);
}
template<typename Derived>
SC_INLINE vector<typename Derived::Scalar> row(Derived &&mat, integer idx) {
    return mat.row(idx);
}

template<typename Derived>
SC_INLINE vector<typename Derived::Scalar> col(const Derived &mat, integer idx) {
    return mat.col(idx);
}
template<typename Derived>
SC_INLINE vector<typename Derived::Scalar> col(Derived &&mat, integer idx) {
    return mat.col(idx);
}

// Spans
//------------------------------------------------------------------------------------------------------------------

template<typename Derived, typename T>
SC_INLINE matrix<typename Derived::Scalar> seq(const Derived &mat,
                                               const std::array<T,3> &rows_,
                                               const std::array<T,3> &cols_) {

    integer starting_row = rows_[0];
    integer starting_col = cols_[0];
    integer span_rows = (rows_[2]-starting_row)/rows_[1];
    integer span_cols = (cols_[2]-starting_col)/cols_[1];
    matrix<typename Derived::Scalar> out(span_rows,span_cols);
    for (integer i=0; i<span_rows; ++i) {
        for (integer j=0; j<span_cols; ++j) {
            out(i,j) = mat(starting_row+i,starting_col+j);
        }
    }
    return out;
}


template<typename Derived, typename T>
SC_INLINE matrix<typename Derived::Scalar> seq(const Derived &mat,
                                               const vector<T> &rows_,
                                               const vector<T> &cols_) {

    matrix<typename Derived::Scalar> out(size(rows_),size(cols_));
    for (integer i=0; i<size(rows_); ++i) {
        for (integer j=0; j<size(cols_); ++j) {
            out(i,j) = mat(rows_(i),cols_(j));
        }
    }
    return out;
}


//template<typename Derived0, typename Derived1>
//SC_INLINE matrix<typename Derived::Scalar> seq(const Derived0 &mat,
//                                               const Derived1 &mat_idx) {

//    matrix<typename Derived::Scalar> out(rows(mat_idx0),cols(mat_idx1));
//    for (integer i=0; i<rows(mat_idx0); ++i) {
//        for (integer j=0; j<cols(mat_idx1); ++j) {
//            out(i,j) = mat(mat_idx(i,j));
//        }
//    }
//    return out;
//}


//template<typename Derived, typename T, typename U>
//SC_INLINE auto seq(Derived &mat,
//                   const vector<U> & rows_,
//                   const vector<U> & cols_)
// -> typename Eigen::Map<matrix<typename Derived::Scalar>,0,
//         Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>> {

//    SC_ASSERT(rows_.size()<4 && cols_.size()<4, "INVALID INDEX FOR MATRIX");

//    integer istride=1,ostride=1;
//    if (rows_.size()==3)
//        istride = *(rows_.begin()+1);
//    if (cols_.size()==3)
//        ostride = *(cols_.begin()+1);

//    integer starting_row = *rows_.begin();
//    integer starting_col = *cols_.begin();
//    integer span_rows, span_cols;
//    if (rows_.size()==3)
//        span_rows = (*(rows_.begin()+2)-starting_row)/istride;
//    else
//        span_rows = (*(rows_.begin()+1)-starting_row)/istride;
//    if (cols_.size()==3)
//        span_cols = (*(cols_.begin()+2)-starting_col)/ostride;
//    else
//        span_cols = (*(cols_.begin()+1)-starting_col)/ostride;

//    Eigen::Map<matrix<typename Derived::Scalar>,0,
//            Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic> > out(mat.data()+starting_row*mat.cols()+starting_col,
//                                                  span_rows, span_cols,
//                                                  Eigen::Stride<Eigen::Dynamic,
//                                                  Eigen::Dynamic>(ostride*mat.outerStride(),
//                                                    istride*mat.innerStride()));

//    return out;
//}

//template<typename Derived, typename T>
//SC_INLINE matrix<typename Derived::Scalar> seq(const Derived &mat,
//                                               const std::initializer_list<T> &rows_,
//                                               const std::initializer_list<T> &cols_) {

//    SC_ASSERT(rows_.size()<4 && cols_.size()<4, "INVALID INDEX FOR MATRIX");

//    integer starting_row = *rows_.begin();
//    integer starting_col = *cols_.begin();
//    integer span_rows = (*(rows_.begin()+2)-starting_row);
//    integer span_cols = (*(cols_.begin()+2)-starting_col);
//    if (rows_.size()==3)
//        span_rows /= *(rows_.begin()+1);
//    if (cols_.size()==3)
//        span_cols /= *(cols_.begin()+1);

//    matrix<typename Derived::Scalar> out(span_rows,span_cols);
//    for (integer i=0; i<span_rows; ++i) {
//        for (integer j=0; j<span_cols; ++j) {
//            out(i,j) = mat(starting_row+i,starting_col+j);
//        }
//    }

//    return out;
//}


template<typename Derived, typename T, typename U>
SC_INLINE auto seq(Derived &mat,
                   std::initializer_list<U> &&  rows_,
                   std::initializer_list<U> &&cols_)
 -> typename Eigen::Map<matrix<typename Derived::Scalar>,0,
         Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>> {

    SC_ASSERT(rows_.size()<4 && cols_.size()<4, "INVALID INDEX FOR MATRIX");

    integer istride=1,ostride=1;
    if (rows_.size()==3)
        istride = *(rows_.begin()+1);
    if (cols_.size()==3)
        ostride = *(cols_.begin()+1);

    integer starting_row = *rows_.begin();
    integer starting_col = *cols_.begin();
    integer span_rows, span_cols;
    if (rows_.size()==3)
        span_rows = (*(rows_.begin()+2)-starting_row)/istride;
    else
        span_rows = (*(rows_.begin()+1)-starting_row)/istride;
    if (cols_.size()==3)
        span_cols = (*(cols_.begin()+2)-starting_col)/ostride;
    else
        span_cols = (*(cols_.begin()+1)-starting_col)/ostride;

    Eigen::Map<matrix<typename Derived::Scalar>,0,
            Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic> > out(mat.data()+starting_row*mat.cols()+starting_col,
                                                  span_rows, span_cols,
                                                  Eigen::Stride<Eigen::Dynamic,
                                                  Eigen::Dynamic>(ostride*mat.outerStride(),
                                                    istride*mat.innerStride()));

    return out;
}


template<typename Derived, typename T, typename U>
SC_INLINE auto seq(Derived &mat,
                  T idx,
                  std::initializer_list<U> &&cols_)
-> typename Eigen::Map<matrix<typename Derived::Scalar>,0,
        Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>> {

   SC_ASSERT(cols_.size()<4, "INVALID INDEX FOR MATRIX");

   integer istride=1,ostride=1;
   if (cols_.size()==3)
       ostride = *(cols_.begin()+1);

   integer starting_row = idx;
   integer starting_col = *cols_.begin();
   integer span_rows = 1, span_cols;
   if (cols_.size()==3)
      span_cols = (*(cols_.begin()+2)-starting_col)/ostride;
   else
      span_cols = (*(cols_.begin()+1)-starting_col)/ostride;

   Eigen::Map<matrix<typename Derived::Scalar>,0,
           Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>> out(mat.data()+starting_row*mat.cols()+starting_col,
                                                 span_rows, span_cols,
                                                 Eigen::Stride<Eigen::Dynamic,
                                                 Eigen::Dynamic>(ostride*mat.outerStride(),
                                                   istride*mat.innerStride()));

   return out;
}


template<typename Derived, typename T, typename U>
SC_INLINE auto seq(Derived &mat,
                 std::initializer_list<U> &&rows_,
                 T idx)
-> typename Eigen::Map<matrix<typename Derived::Scalar>,0,
       Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>> {

  SC_ASSERT(rows_.size()<4, "INVALID INDEX FOR MATRIX");

  integer istride=1,ostride=1;
  if (rows_.size()==3)
      istride = *(rows_.begin()+1);

  integer starting_row = *rows_.begin();
  integer starting_col = idx;
  integer span_rows, span_cols=1;
  if (rows_.size()==3)
      span_rows = (*(rows_.begin()+2)-starting_row)/istride;
  else
      span_rows = (*(rows_.begin()+1)-starting_row)/istride;

  Eigen::Map<matrix<typename Derived::Scalar>,0,
          Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic> > out(mat.data()+starting_row*mat.cols()+starting_col,
                                                span_rows, span_cols,
                                                Eigen::Stride<Eigen::Dynamic,
                                                Eigen::Dynamic>(ostride*mat.outerStride(),
                                                  istride*mat.innerStride()));

  return out;
}



//------------------------------------------------------------------------------------------------------------------

#endif

SC_INLINE matrix<real> rand(integer m, integer n) {
    return matrix<real>::Random(m,n);
}

template<typename T>
SC_INLINE vector<T> zeros(integer m) {
    return vector<T>::Zero(m);
}
template<typename T>
SC_INLINE matrix<T> zeros(integer m, integer n) {
    return matrix<T>::Zero(m,n);
}
SC_INLINE vector<real> zeros(integer m) {
    return vector<real>::Zero(m);
}
SC_INLINE matrix<real> zeros(integer m, integer n) {
    return matrix<real>::Zero(m,n);
}

template<typename T>
SC_INLINE vector<T> ones(integer m) {
    return vector<T>::Ones(m);
}
template<typename T>
SC_INLINE matrix<T> ones(integer m, integer n) {
    return matrix<T>::Ones(m,n);
}
SC_INLINE vector<real> ones(integer m) {
    return vector<real>::Ones(m);
}
SC_INLINE matrix<real> ones(integer m, integer n) {
    return matrix<real>::Ones(m,n);
}

SC_INLINE matrix<real> eye(integer m, integer n) {
    return matrix<real>::Identity(m,n);
}

template<typename Derived, typename std::enable_if<!std::is_arithmetic<Derived>::value,bool>::type=0>
SC_INLINE void fill(Eigen::PlainObjectBase<Derived> &mat, integer num) {
    //! In-place
    std::fill(mat.data(),mat.data()+size(mat),num);
}
template<typename T>
SC_INLINE vector<T> fill(integer cap, integer num) {
    vector<T> mat(cap);
    std::fill(mat.data(),mat.data()+size(mat),num);
    return mat;
}

template<typename T>
SC_INLINE matrix<T> arange(integer m, integer n) {
    matrix<T> out(n-m,1);
    std::iota(out.data(),out.data()+n-m,m);
    return out;
}
template<typename T>
SC_INLINE vector<T> arange(integer n) {
    vector<T> out(n);
    std::iota(out.data(),out.data()+n,0);
    return out;
}

SC_INLINE vector<real> linspace(real low, real high, integer num=50)
{
    //! Linearly spaced points. Dont parametrise!
    return vector<real>::LinSpaced(num,low,high);
}


template<typename Derived, typename T>
SC_INLINE matrix<T> pow(T num, const Eigen::MatrixBase<Derived> &mat); // Forward declare

SC_INLINE vector<real> logspace(real low, real high, integer num=50)
{
    //! Linearly spaced points. Dont parametrise!
    vector<real> out = vector<real>::LinSpaced(Eigen::Sequential,num,low,high);
    return pow(10.,out);
}

template<typename Derived>
SC_INLINE real max(const Eigen::DenseBase<Derived> &mat) {
    return mat.maxCoeff();
}
template<typename Derived>
SC_INLINE real min(const Eigen::DenseBase<Derived> &mat) {
    return mat.minCoeff();
}
template<typename Derived>
SC_INLINE integer argmax(const Eigen::DenseBase<Derived> &mat) {
    integer i; mat.minCoeff(&i);
    return i;
}
template<typename Derived>
SC_INLINE integer argmin(const Eigen::DenseBase<Derived> &mat) {
    integer i; mat.minCoeff(&i);
    return i;
}

template<typename Derived>
SC_INLINE boolean isempty(const Eigen::PlainObjectBase<Derived> &mat) {
    return mat.rows()==0 && mat.cols()==0;
}

template<typename T, typename U>
SC_INLINE vector<T>
append(const vector<T> &arr, U num) {
    //! Append to the end of Eigen vector (similar to push_back)

    vector<T> new_arr(arr.rows()+1);
    new_arr.head(arr.rows()) = arr;
    new_arr[arr.rows()] = static_cast<T>(num);
    return new_arr;
}

template<typename Derived>
SC_INLINE matrix<typename Derived::Scalar> reshape(const Eigen::PlainObjectBase<Derived> &mat, integer m, integer n) {
    //! Reshape an existing matrix
    SC_ASSERT(m*n==size(mat),"CANNOT RESHAPE MATRICES TO A DIFFERENT SIZE");
    matrix<typename Derived::Scalar> out(m,n);
    std::copy(mat.data(),mat.data()+m*n,out.data());
    return out;
}

template<typename Derived>
SC_INLINE matrix<typename Derived::Scalar>
ravel(const Eigen::PlainObjectBase<Derived> &arr) {
    //! Ravel/flatten a matrix. Makes a copy
    return reshape(arr,size(arr),1);
}

template<typename Derived>
SC_INLINE matrix<typename Derived::Scalar>
flatten(const Eigen::PlainObjectBase<Derived> &arr) {
    //! Ravel/flatten a matrix. Makes a copy
    return reshape(arr,size(arr),1);
}

template<typename T, typename U = T>
std::tuple<matrix<integer>,matrix<integer> >
SC_INLINE find(const Eigen::PlainObjectBase<T> &arr,
         U num, real tolerance=1e-14) {
    //! find the occurence of a value in a matrix
    std::vector<integer> idx_rows;
    std::vector<integer> idx_cols;
    idx_rows.clear(); idx_cols.clear();
    for (auto i=0; i<arr.rows();++i)
    {
        for (auto j=0; j<arr.cols();++j)
        {
            if (static_cast<real>(abs(arr(i,j)-num))<tolerance)
            {
                idx_rows.push_back(i);
                idx_cols.push_back(j);
            }
        }
    }

    return std::make_tuple(
                Eigen::Map<matrix<integer>>
                (idx_rows.data(),idx_rows.size(),1),
                Eigen::Map<matrix<integer>>
                (idx_cols.data(),idx_cols.size(),1));
}


template<typename Derived>
SC_INLINE matrix<typename Derived::Scalar> fliplr(const Eigen::DenseBase<Derived> &mat) {
    //! Reverse a matrix left to right
    SC_ASSERT(mat.cols()!=1,"CANNOT FLIP A COLUMN VECTOR COLUMN-WISE");
    matrix<typename Derived::Scalar> out(mat.rows(),mat.cols());
    if (mat.rows()==1) {
        return mat.reverse();
    }
    else {
        for (auto i=0; i<mat.rows(); ++i) {
            out.row(i) = mat.row(i).reverse();
        }
    }
    return out;
}
template<typename T>
SC_INLINE vector<T> fliplr(const vector<T> &vec) {
    //! Reverse a vector
    return vec.reverse();
}

template<typename Derived>
SC_INLINE matrix<typename Derived::Scalar> flipud(const Eigen::DenseBase<Derived> &mat) {
    //! Reverse a matrix top to bottom
    SC_ASSERT(mat.rows()!=1,"CANNOT FLIP A ROW VECTOR ROW-WISE");
    matrix<typename Derived::Scalar> out(mat.rows(),mat.cols());
    if (mat.cols()==1) {
        return mat.reverse();
    }
    else {
        for (auto i=0; i<mat.cols(); ++i) {
            out.col(i) = mat.col(i).reverse();
        }
    }
    return out;
}
template<typename T>
SC_INLINE vector<T> flipud(const vector<T> &vec) {
    //! Reverse a vector
    return vec.reverse();
}

template<typename Derived>
SC_INLINE matrix<real> copy(const Eigen::DenseBase<Derived> &mat) {
    // Make a copy of an existing matrix
//    matrix<real> out(mat.rows(),mat.cols());
//    std::copy(mat.data(),mat.data()+mat.rows()*mat.cols(),out.data());
//    return out;
    return matrix<real>(mat);
}

template<typename Derived>
SC_INLINE matrix<real> copy(Eigen::DenseBase<Derived> &&mat) {
    // Move overload (named copy) for simplicity
    return matrix<real>(std::move(mat));
}

template<typename Derived>
SC_INLINE matrix<real> sort(const Eigen::PlainObjectBase<Derived> &mat) {
    //! Sort a matrix
    auto cmat = copy(mat);
    std::sort(cmat.data(),cmat.data()+size(cmat));
    return cmat;
}

template <typename T>
SC_INLINE matrix<integer> argsort(const matrix<T> &v) {
    //! Get indices of a sorted vector
    // INITIALIZE ORIGINAL INDEX LOCATIONS
    matrix<integer> idx(size(v));
    fill(idx,0);
    // SORT INDICES BY COMPARING VALUES IN V USING LAMBDA FUNCTION
    std::sort(idx.data(), idx.data()+size(idx),[&v](integer i1, integer i2) {return v[i1] < v[i2];});
    return idx;
}


//template<typename Derived>
//SC_INLINE vector<typename Derived::Scalar>
//unique(const Eigen::PlainObjectBase<Derived> &arr, boolean keep_order=true) {
//    std::vector<typename Derived::Scalar> vec(arr.data(),arr.data()+size(arr));
//    if (keep_order) unsort_unique(vec);
//    else {
//        std::sort( vec.begin(), vec.end() );
//        vec.erase(std::unique( vec.begin(), vec.end() ), vec.end() );
//    }
//    return Eigen::Map<vector<typename Derived::Scalar>>(vec.data(),vec.size());
//}


//enum do_unique {return_index, return_inverse};
//enum {return_index, return_inverse};

// THIS BECOMES AMBIGUOUS
//template<typename T>
//SC_INLINE std::tuple<std::vector<typename Eigen::PlainObjectBase<T>::Scalar>,std::vector<integer> >
//unique(const Eigen::PlainObjectBase<T> &arr, boolean keep_order=true, boolean retrun_index=true) {
//    //! RETURNS UNIQUE VALUES AND UNIQUE INDICES OF AN EIGEN MATRIX
//    std::vector<typename Eigen::PlainObjectBase<T>::Scalar> uniques;
//    std::vector<integer> idx;

//    for (auto i=0; i<arr.rows(); ++i) {
//        bool isunique = true;
//        for (auto j=0; j<=i; ++j) {
//            if (arr(i)==arr(j) && i!=j) {
//                isunique = false;
//                break;
//            }
//        }

//        if (isunique==true) {
//            uniques.push_back(arr(i));
//            idx.push_back(i);
//        }
//    }

//    // SORT UNIQUE VALUES
//    auto sorter = argsort(uniques);
//    std::sort(uniques.begin(),uniques.end());
//    std::vector<integer> idx_sorted(idx.size());
//    for (size_t i=0; i<uniques.size();++i) {
//        idx_sorted[i] = idx[sorter[i]];
//    }

//    std::tuple<std::vector<typename Eigen::PlainObjectBase<T>::Scalar>,std::vector<integer> >
//            uniques_idx = std::make_tuple(uniques,idx_sorted);

//    return uniques_idx;
//}


template <typename Derived>
matrix<typename Derived::Scalar> unique(const Eigen::PlainObjectBase<Derived> &arr, integer axis=-1,
                                          boolean keep_order=true, boolean consider_sort=false)
{
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

    if (axis==-1) {
        std::vector<typename Derived::Scalar> vec(arr.data(),arr.data()+size(arr));
        if (keep_order) unsort_unique(vec);
        else {
            std::sort( vec.begin(), vec.end() );
            vec.erase(std::unique( vec.begin(), vec.end() ), vec.end() );
        }
        return Eigen::Map<vector<typename Derived::Scalar>>(vec.data(),vec.size());
    }

    matrix<typename Derived::Scalar> &&arr_ = arr; if (axis==1) arr_ = arr.transpose();

    std::vector<std::vector<typename Derived::Scalar>> input(arr_.rows());
    for(auto i=0; i<arr_.rows(); ++i) {
        input[i] = std::vector<typename Derived::Scalar>(arr_.cols());
        for (auto j=0; j<arr_.cols(); ++j) {
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


    matrix<typename Derived::Scalar> out(input.size(),arr_.cols());
    for(size_t i=0; i<input.size(); ++i) {
        for (size_t j=0; j<input[i].size(); ++j) {
            out(i,j) = input[i][j];
        }
    }

    if (axis==1) out.transposeInPlace();

    return out;

}

template<typename Derived0, typename Derived1>
SC_INLINE matrix<typename Derived0::Scalar>
hstack(const Eigen::PlainObjectBase<Derived0> &mat0, const Eigen::PlainObjectBase<Derived1> &mat1) {
    SC_ASSERT(mat0.rows()==mat1.rows(),"CONCATENATING MATRICES SHOULD HAVE THE SAME NUMBER OF ROWS");
    matrix<typename Derived0::Scalar> out(mat0.rows(),mat0.cols()+mat1.cols());
    for (integer i=0; i<mat0.rows(); ++i) {
        integer j=0;
        for (; j<mat0.cols(); ++j) {
            out(i,j) = mat0(i,j);
        }
        integer counter = 0;
        for (; j<mat0.cols()+mat1.cols(); ++j) {
            out(i,j) = mat1(i,counter);
            counter++;
        }
    }
    return out;
}
template<typename Derived0, typename Derived1>
SC_INLINE matrix<typename Derived0::Scalar>
hstack(Eigen::PlainObjectBase<Derived0> &&mat0, Eigen::PlainObjectBase<Derived1> &&mat1) {
    SC_ASSERT(mat0.rows()==mat1.rows(),"CONCATENATING MATRICES SHOULD HAVE THE SAME NUMBER OF ROWS");
    matrix<typename Derived0::Scalar> out(mat0.rows(),mat0.cols()+mat1.cols());
    for (integer i=0; i<mat0.rows(); ++i) {
        integer j=0;
        for (; j<mat0.cols(); ++j) {
            out(i,j) = mat0(i,j);
        }
        integer counter = 0;
        for (; j<mat0.cols()+mat1.cols(); ++j) {
            out(i,j) = mat1(i,counter);
            counter++;
        }
    }
    return out;
}


template<typename Derived0, typename Derived1>
SC_INLINE matrix<typename Derived0::Scalar>
vstack(const Eigen::PlainObjectBase<Derived0> &mat0, const Eigen::PlainObjectBase<Derived1> &mat1) {
    SC_ASSERT(mat0.cols()==mat1.cols(),"CONCATENATING MATRICES SHOULD HAVE THE SAME NUMBER OF COLUMNS");
    matrix<typename Derived0::Scalar> out(mat0.rows()+mat1.rows(),mat0.cols());
    std::copy(mat0.data(),mat0.data()+size(mat0),out.data());
    std::copy(mat1.data(),mat1.data()+size(mat1),out.data()+size(mat0));
    return out;
}
template<typename Derived0, typename Derived1>
SC_INLINE matrix<typename Derived0::Scalar>
vstack(Eigen::PlainObjectBase<Derived0> &&mat0, Eigen::PlainObjectBase<Derived1> &&mat1) {
    SC_ASSERT(mat0.cols()==mat1.cols(),"CONCATENATING MATRICES SHOULD HAVE THE SAME NUMBER OF COLUMNS");
    matrix<typename Derived0::Scalar> out(mat0.rows()+mat1.rows(),mat0.cols());
    std::copy(mat0.data(),mat0.data()+size(mat0),out.data());
    std::copy(mat1.data(),mat1.data()+size(mat1),out.data()+size(mat0));
    return out;
}


template<typename Derived>
SC_INLINE matrix<typename Derived::Scalar>
repmat(const Eigen::PlainObjectBase<Derived> &mat, integer a, integer b=None) {
    //! Equivalent to MATLAB's repmat
    b = b==None ? a : b;
    return mat.replicate(a,b);
//    matrix<typename Derived::Scalar> out(mat.rows()*a,mat.cols()*b);
//    for (integer i=0; i<a; ++i) {
//        out.block(i*mat.rows(),0,mat.rows(),mat.cols()) = mat;
//    }
//    for (integer j=1; j<b; ++j) {
//        out.block(0,j*mat.cols(),out.rows(),mat.cols()) = out.block(0,0,out.rows(),mat.cols());
//    }
//    return out;
}
template<typename Derived>
SC_INLINE matrix<typename Derived::Scalar>
repmat(Eigen::PlainObjectBase<Derived> &&mat, integer a, integer b=None) {
    //! Equivalent to MATLAB's repmat
    b = b==None ? a : b;
    return mat.replicate(a,b);
//    matrix<typename Derived::Scalar> out(mat.rows()*a,mat.cols()*b);
//    for (integer i=0; i<a; ++i) {
//        out.block(i*mat.rows(),0,mat.rows(),mat.cols()) = mat;
//    }
//    for (integer j=1; j<b; ++j) {
//        out.block(0,j*mat.cols(),out.rows(),mat.cols()) = out.block(0,0,out.rows(),mat.cols());
//    }
//    return out;
}

template<typename Derived>
SC_INLINE matrix<typename Derived::Scalar>
tile(const Eigen::PlainObjectBase<Derived> &mat, integer reps, integer axis=0) {
    SC_ASSERT(axis<2,"AXIS ARGUMENT SHOULD BE EITHER 0 OR 1");
    matrix<typename Derived::Scalar> out;
    if (axis==0) {
        out = zeros<typename Derived::Scalar>(mat.rows()*reps,mat.cols());
        for (integer i=0; i<reps; ++i) {
            out.block(i*mat.rows(),0,mat.rows(),mat.cols()) = mat;
        }
    }
    else if (axis==1) {
        out = zeros<typename Derived::Scalar>(mat.rows(),mat.cols()*reps);
        for (integer i=0; i<reps; ++i) {
            out.block(0,i*mat.cols(),mat.rows(),mat.cols()) = mat;
        }
    }
    return out;
}
template<typename Derived>
SC_INLINE matrix<typename Derived::Scalar>
tile(Eigen::PlainObjectBase<Derived> &&mat, integer reps, integer axis=0) {
    SC_ASSERT(axis<2,"AXIS ARGUMENT SHOULD BE EITHER 0 OR 1");
    matrix<typename Derived::Scalar> out;
    if (axis==0) {
        out = zeros<typename Derived::Scalar>(mat.rows()*reps,mat.cols());
        for (integer i=0; i<reps; ++i) {
            out.block(i*mat.rows(),0,mat.rows(),mat.cols()) = mat;
        }
    }
    else if (axis==1) {
        out = zeros<typename Derived::Scalar>(mat.rows(),mat.cols()*reps);
        for (integer i=0; i<reps; ++i) {
            out.block(0,i*mat.cols(),mat.rows(),mat.cols()) = mat;
        }
    }
    return out;
}

template<typename Derived>
SC_INLINE vector<typename Derived::Scalar>
repeat(const Eigen::PlainObjectBase<Derived> &mat, integer reps) {
    SC_ASSERT(mat.cols()==1 || mat.rows()==1,"REPEAT ONLY WORKS FOR VECTORS AND/OR MATRICES WITH EITHER cols OR rows EQUAL TO 1");
    return mat.rows()==1 ? ravel(tile((matrix<typename Derived::Scalar>)mat.transpose(),reps,1)) :
                           ravel(tile(mat,reps,1));
}

template<typename Derived, typename T=real>
SC_INLINE matrix<typename Derived::Scalar>
zero_less_than(Eigen::PlainObjectBase<Derived> &mat, T tol=1e-14) {
    //! Make elements of a matrix which are less than tolerance, 0 (in-place)
    for (auto i=0; i<mat.rows(); ++i) {
        for (auto j=0; j<mat.cols(); ++j) {
            if (std::abs(mat(i,j))<tol) mat(i,j) = (typename Derived::Scalar)0.;
        }
    }
}

// Math
template<typename Derived>
SC_INLINE typename Derived::Scalar sum(const Eigen::DenseBase<Derived> &mat) {
    return mat.sum();
}
template<typename Derived>
SC_INLINE vector<typename Derived::Scalar> sum(const Eigen::PlainObjectBase<Derived> &mat, integer axis) {
    if (mat.cols()==1 || mat.rows()==1) {
      vector<typename Derived::Scalar> out(1); out[0] = mat.sum();
      return out;
    }
    else {
        if (axis==0) {
          vector<typename Derived::Scalar> out(mat.cols());
          for (auto i=0; i<mat.cols(); ++i)
              out[i] = mat.col(i).sum(); 
          return out;
        }
        else if (axis==1) {
          vector<typename Derived::Scalar> out(mat.rows());
          for (auto i=0; i<mat.rows(); ++i)
              out[i] = mat.row(i).sum(); 
          return out;
        }
        else {
            SC_ASSERT(false,"INVALID AXIS ARGUMENT");
            return vector<typename Derived::Scalar>{};
        }
    }
}

template<typename Derived>
SC_INLINE typename Derived::Scalar prod(const Eigen::DenseBase<Derived> &mat) {
    return mat.prod();
}
template<typename Derived>
SC_INLINE vector<typename Derived::Scalar> prod(const Eigen::PlainObjectBase<Derived> &mat, integer axis) {
    if (mat.cols()==1 || mat.rows()==1) {
      vector<typename Derived::Scalar> out(1); out[0] = mat.prod();
      return out;
    }
    else {
        if (axis==0) {
          vector<typename Derived::Scalar> out(mat.cols());
          for (auto i=0; i<mat.cols(); ++i)
              out[i] = mat.col(i).prod(); 
          return out;
        }
        else if (axis==1) {
          vector<typename Derived::Scalar> out(mat.rows());
          for (auto i=0; i<mat.rows(); ++i)
              out[i] = mat.row(i).prod(); 
          return out;
        }
        else {
            SC_ASSERT(false,"INVALID AXIS ARGUMENT");
            return vector<typename Derived::Scalar>{};
        }
    }
}

template<typename Derived>
SC_INLINE typename Derived::Scalar mean(const Eigen::DenseBase<Derived> &mat) {
    return mat.mean();
}
template<typename Derived>
SC_INLINE vector<typename Derived::Scalar> mean(const Eigen::PlainObjectBase<Derived> &mat, integer axis) {
    if (mat.cols()==1 || mat.rows()==1) {
      vector<typename Derived::Scalar> out(1); out[0] = mat.mean();
      return out;
    }
    else {
        if (axis==0) {
          vector<typename Derived::Scalar> out(mat.cols());
          for (auto i=0; i<mat.cols(); ++i)
              out[i] = mat.col(i).mean(); 
          return out;
        }
        else if (axis==1) {
          vector<typename Derived::Scalar> out(mat.rows());
          for (auto i=0; i<mat.rows(); ++i)
              out[i] = mat.row(i).mean(); 
          return out;
        }
        else {
            SC_ASSERT(false,"INVALID AXIS ARGUMENT");
            return vector<typename Derived::Scalar>{};
        }
    }
}

template<typename Derived>
SC_INLINE matrix<typename Derived::Scalar> abs(const Eigen::MatrixBase<Derived> &mat) {
    return mat.array().abs();
}
template<typename Derived>
SC_INLINE vector<typename Derived::Scalar> abs(const Eigen::PlainObjectBase<Derived> &mat, integer axis) {
    if (mat.cols()==1 || mat.rows()==1) {
      vector<typename Derived::Scalar> out(1); out[0] = mat.abs();
      return out;
    }
    else {
        if (axis==0) {
          vector<typename Derived::Scalar> out(mat.cols());
          for (auto i=0; i<mat.cols(); ++i)
              out[i] = mat.col(i).abs(); 
          return out;
        }
        else if (axis==1) {
          vector<typename Derived::Scalar> out(mat.rows());
          for (auto i=0; i<mat.rows(); ++i)
              out[i] = mat.row(i).abs(); 
          return out;
        }
        else {
            SC_ASSERT(false,"INVALID AXIS ARGUMENT");
            return vector<typename Derived::Scalar>{};
        }
    }
}

template<typename Derived, typename T,
         typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0>
SC_INLINE matrix<typename Derived::Scalar> pow(const Eigen::MatrixBase<Derived> &mat, T num) {
    return mat.array().pow(num);
}
template<typename Derived, typename T>
SC_INLINE matrix<T> pow(T num, const Eigen::MatrixBase<Derived> &mat) {
    matrix<T> out = Eigen::pow(num,mat.array());
    return out;
}
template<typename Derived0, typename Derived1>
SC_INLINE matrix<real> pow(const Eigen::MatrixBase<Derived0> &mat, const Eigen::MatrixBase<Derived1> &nums) {
    return mat.array().pow(nums.array());
}
template<typename Derived>
SC_INLINE matrix<real> log(const Eigen::MatrixBase<Derived> &mat) {
    return mat.array().log();
}
template<typename Derived>
SC_INLINE matrix<real> log10(const Eigen::MatrixBase<Derived> &mat) {
    return mat.array().log10();
}
template<typename Derived>
SC_INLINE matrix<real> log1p(const Eigen::MatrixBase<Derived> &mat) {
    return mat.array().log1p();
}
template<typename Derived>
SC_INLINE matrix<real> exp(const Eigen::MatrixBase<Derived> &mat) {
    return mat.array().exp();
}

// Trig
template<typename Derived>
SC_INLINE matrix<real> sin(const Eigen::MatrixBase<Derived> &mat) {
    return mat.array().sin();
}
template<typename Derived>
SC_INLINE vector<typename Derived::Scalar> sin(const Eigen::PlainObjectBase<Derived> &mat, integer axis) {
    if (mat.cols()==1 || mat.rows()==1) {
      vector<typename Derived::Scalar> out(1); out[0] = mat.sin();
      return out;
    }
    else {
        if (axis==0) {
          vector<typename Derived::Scalar> out(mat.cols());
          for (auto i=0; i<mat.cols(); ++i)
              out[i] = mat.col(i).sin(); 
          return out;
        }
        else if (axis==1) {
          vector<typename Derived::Scalar> out(mat.rows());
          for (auto i=0; i<mat.rows(); ++i)
              out[i] = mat.row(i).sin(); 
          return out;
        }
        else {
            SC_ASSERT(false,"INVALID AXIS ARGUMENT");
            return vector<typename Derived::Scalar>{};
        }
    }
}

template<typename Derived>
SC_INLINE matrix<real> cos(const Eigen::MatrixBase<Derived> &mat) {
    return mat.array().cos();
}
template<typename Derived>
SC_INLINE vector<typename Derived::Scalar> cos(const Eigen::PlainObjectBase<Derived> &mat, integer axis) {
    if (mat.cols()==1 || mat.rows()==1) {
      vector<typename Derived::Scalar> out(1); out[0] = mat.cos();
      return out;
    }
    else {
        if (axis==0) {
          vector<typename Derived::Scalar> out(mat.cols());
          for (auto i=0; i<mat.cols(); ++i)
              out[i] = mat.col(i).cos(); 
          return out;
        }
        else if (axis==1) {
          vector<typename Derived::Scalar> out(mat.rows());
          for (auto i=0; i<mat.rows(); ++i)
              out[i] = mat.row(i).cos(); 
          return out;
        }
        else {
            SC_ASSERT(false,"INVALID AXIS ARGUMENT");
            return vector<typename Derived::Scalar>{};
        }
    }
}

template<typename Derived>
SC_INLINE matrix<real> tan(const Eigen::MatrixBase<Derived> &mat) {
    return mat.array().tan();
}
template<typename Derived>
SC_INLINE vector<typename Derived::Scalar> tan(const Eigen::PlainObjectBase<Derived> &mat, integer axis) {
    if (mat.cols()==1 || mat.rows()==1) {
      vector<typename Derived::Scalar> out(1); out[0] = mat.tan();
      return out;
    }
    else {
        if (axis==0) {
          vector<typename Derived::Scalar> out(mat.cols());
          for (auto i=0; i<mat.cols(); ++i)
              out[i] = mat.col(i).tan(); 
          return out;
        }
        else if (axis==1) {
          vector<typename Derived::Scalar> out(mat.rows());
          for (auto i=0; i<mat.rows(); ++i)
              out[i] = mat.row(i).tan(); 
          return out;
        }
        else {
            SC_ASSERT(false,"INVALID AXIS ARGUMENT");
            return vector<typename Derived::Scalar>{};
        }
    }
}

template<typename Derived>
SC_INLINE matrix<real> sinh(const Eigen::MatrixBase<Derived> &mat) {
    return mat.array().sinh();
}
template<typename Derived>
SC_INLINE vector<typename Derived::Scalar> sinh(const Eigen::PlainObjectBase<Derived> &mat, integer axis) {
    if (mat.cols()==1 || mat.rows()==1) {
      vector<typename Derived::Scalar> out(1); out[0] = mat.sinh();
      return out;
    }
    else {
        if (axis==0) {
          vector<typename Derived::Scalar> out(mat.cols());
          for (auto i=0; i<mat.cols(); ++i)
              out[i] = mat.col(i).sinh(); 
          return out;
        }
        else if (axis==1) {
          vector<typename Derived::Scalar> out(mat.rows());
          for (auto i=0; i<mat.rows(); ++i)
              out[i] = mat.row(i).sinh(); 
          return out;
        }
        else {
            SC_ASSERT(false,"INVALID AXIS ARGUMENT");
            return vector<typename Derived::Scalar>{};
        }
    }
}

template<typename Derived>
SC_INLINE matrix<real> cosh(const Eigen::MatrixBase<Derived> &mat) {
    return mat.array().cosh();
}
template<typename Derived>
SC_INLINE vector<typename Derived::Scalar> cosh(const Eigen::PlainObjectBase<Derived> &mat, integer axis) {
    if (mat.cols()==1 || mat.rows()==1) {
      vector<typename Derived::Scalar> out(1); out[0] = mat.cosh();
      return out;
    }
    else {
        if (axis==0) {
          vector<typename Derived::Scalar> out(mat.cols());
          for (auto i=0; i<mat.cols(); ++i)
              out[i] = mat.col(i).cosh(); 
          return out;
        }
        else if (axis==1) {
          vector<typename Derived::Scalar> out(mat.rows());
          for (auto i=0; i<mat.rows(); ++i)
              out[i] = mat.row(i).cosh(); 
          return out;
        }
        else {
            SC_ASSERT(false,"INVALID AXIS ARGUMENT");
            return vector<typename Derived::Scalar>{};
        }
    }
}

template<typename Derived>
SC_INLINE matrix<real> tanh(const Eigen::MatrixBase<Derived> &mat) {
    return mat.array().tanh();
}
template<typename Derived>
SC_INLINE vector<typename Derived::Scalar> tanh(const Eigen::PlainObjectBase<Derived> &mat, integer axis) {
    if (mat.cols()==1 || mat.rows()==1) {
      vector<typename Derived::Scalar> out(1); out[0] = mat.tanh();
      return out;
    }
    else {
        if (axis==0) {
          vector<typename Derived::Scalar> out(mat.cols());
          for (auto i=0; i<mat.cols(); ++i)
              out[i] = mat.col(i).tanh(); 
          return out;
        }
        else if (axis==1) {
          vector<typename Derived::Scalar> out(mat.rows());
          for (auto i=0; i<mat.rows(); ++i)
              out[i] = mat.row(i).tanh(); 
          return out;
        }
        else {
            SC_ASSERT(false,"INVALID AXIS ARGUMENT");
            return vector<typename Derived::Scalar>{};
        }
    }
}


// Linear algebra subroutines
template<typename Derived0, typename T,
         typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0>
SC_INLINE matrix<typename Derived0::Scalar>
operator+(const Eigen::MatrixBase<Derived0> &mat0, T a) {
    return mat0.array()+a;
}
template<typename Derived0, typename T,
         typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0>
SC_INLINE matrix<typename Derived0::Scalar>
operator+(T a, const Eigen::MatrixBase<Derived0> &mat0) {
    return a+mat0.array();
}

template<typename Derived0, typename T,
         typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0>
SC_INLINE matrix<typename Derived0::Scalar>
operator-(const Eigen::MatrixBase<Derived0> &mat0, T a) {
    return mat0.array()-a;
}
template<typename Derived0, typename T,
         typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0>
SC_INLINE matrix<typename Derived0::Scalar>
operator-(T a, const Eigen::MatrixBase<Derived0> &mat0) {
    return a-mat0.array();
}

template<typename Derived0, typename Derived1>
SC_INLINE matrix<typename Derived0::Scalar>
multiply(const Eigen::MatrixBase<Derived0> &mat0, const Eigen::MatrixBase<Derived1> &mat1) {
//    SC_ASSERT(size(mat0)==size(mat1),"MATRICES SHOULD HAVE THE SAME SIZE");
    return mat0.cwiseProduct(mat1);
}
template<typename Derived0, typename Derived1>
SC_INLINE matrix<typename Derived0::Scalar>
operator/(const Eigen::MatrixBase<Derived0> &mat0, const Eigen::MatrixBase<Derived1> &mat1) {
//    SC_ASSERT(size(mat0)==size(mat1),"MATRICES SHOULD HAVE THE SAME SIZE");
    return mat0.cwiseQuotient(mat1);
}
template<typename Derived0, typename T,
    typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0>
SC_INLINE matrix<typename Derived0::Scalar>
operator/(const Eigen::PlainObjectBase<Derived0> &mat0, T a) {
    return mat0.array()/a;
}
template<typename Derived0, typename T,
    typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0>
SC_INLINE matrix<typename Derived0::Scalar>
operator/(T a, const Eigen::PlainObjectBase<Derived0> &mat0) {
    return a/mat0.array();
}

template<typename Derived>
SC_INLINE typename Derived::Scalar inner(const Eigen::PlainObjectBase<Derived> &mat) {
    return std::inner_product(mat.data(),mat.data()+size(mat),mat.data(),0);
}
template<typename Derived>
SC_INLINE typename Derived::Scalar inner(Eigen::PlainObjectBase<Derived> &&mat) {
    return std::inner_product(mat.data(),mat.data()+size(mat),mat.data(),0);
}

template<typename Derived>
SC_INLINE matrix<typename Derived::Scalar> transpose(const Eigen::DenseBase<Derived> &mat) {
    return mat.transpose();
}
template<typename Derived>
SC_INLINE matrix<typename Derived::Scalar> transpose(Eigen::DenseBase<Derived> &&mat) {
    mat.transposeInPlace();
    return mat;
}
template<typename T>
SC_INLINE matrix<T> transpose(const vector<T> &vec) {
    return vec.transpose();
}
template<typename T>
SC_INLINE matrix<T> transpose(vector<T> &&vec) {
    vec.transposeInPlace();
    return vec;
}

template<typename Derived>
SC_INLINE typename Derived::Scalar trace(const Eigen::MatrixBase<Derived> &mat) {
    return mat.trace();
}
template<typename Derived>
SC_INLINE typename Derived::Scalar trace(Eigen::MatrixBase<Derived> &&mat) {
    return mat.trace();
}

template<typename Derived>
SC_INLINE matrix<typename Derived::Scalar> inv(const Eigen::PlainObjectBase<Derived> &mat) {
    return mat.inverse();
}

template<typename Derived>
SC_INLINE typename Derived::Scalar det(Eigen::PlainObjectBase<Derived> &&mat) {
    return mat.determinant();
}
template<typename Derived>
SC_INLINE typename Derived::Scalar det(const Eigen::PlainObjectBase<Derived> &mat) {
    return mat.determinant();
}


template<typename Derived>
SC_INLINE auto diag(const Eigen::MatrixBase<Derived> &mat) -> decltype(mat.diagonal()) {
    return mat.diagonal();
}
template<typename Derived>
SC_INLINE auto diag(Eigen::MatrixBase<Derived> &&mat) -> decltype(mat.diagonal()) {
    return mat.diagonal();
}

template<typename Derived>
SC_INLINE typename Derived::Scalar norm(const Eigen::MatrixBase<Derived> &mat) {
    return mat.squaredNorm();
}
template<typename Derived>
SC_INLINE typename Derived::Scalar norm(Eigen::MatrixBase<Derived> &&mat) {
    return mat.squaredNorm();
}

template<typename Derived0, typename Derived1>
SC_INLINE auto matmul(const Eigen::MatrixBase<Derived0> &mat0,
                      const Eigen::MatrixBase<Derived1> &mat1) -> decltype(mat0*mat1) {
    return mat0*mat1;
}

template<typename T, typename Derived1>
SC_INLINE matrix<T> solve(const matrix<T> &mat, const Eigen::MatrixBase<Derived1> &vec, integer info=linear_solver::non_symmetric) {
    static_assert(std::is_same<T,typename Derived1::Scalar>::value,"DATA TYPES FOR LHS AND RHS SHOULD BE THE SAME");
    if (info==linear_solver::non_symmetric) {return mat.partialPivLu().solve(vec);}
    else if (info==linear_solver::semi_indefinite) {return mat.ldlt().solve(vec);}
    else {return mat.llt().solve(vec);}
}





// Sparse matrices



//template<typename Derived0, typename Derived1>
//SC_INLINE matrix<real> solve(const Eigen::EigenBase<Derived0> &mat, const Eigen::MatrixBase<Derived1> &vec) {
//    return Eigen::LDLT<Eigen::EigenBase<Derived0>>(mat);//.solve(vec);
//}

}


#endif

#endif // BACKEND_EIGEN_H
