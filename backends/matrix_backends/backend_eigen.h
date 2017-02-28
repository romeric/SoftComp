#ifndef BACKEND_EIGEN_H
#define BACKEND_EIGEN_H

#include "commons/commons.h"
#if defined(HAS_EIGEN)



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

// Spans/ Array views
//-----------------------------------------------------------------------------------------------

// #define all -101;
// #define fin -102;

template<typename Derived, typename U>
SC_INLINE auto seq(Derived &mat,
                   std::initializer_list<U> &&  rows_,
                   std::initializer_list<U> &&cols_)
 -> typename Eigen::Map<matrix<typename Derived::Scalar>,Eigen::Unaligned,
         Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>> {
    /* Array slicing
        A[0:5,0:5:2]  ---> seq(A,{0,5},{0,5,2})
    */

#ifndef NDEBUG
    SC_ASSERT(rows_.size()<4 && cols_.size()<4, "INVALID INDEX FOR MATRIX");
#endif


    integer istride=1,ostride=1;
    if (cols_.size()==3)
        istride = *(cols_.begin()+2);
    if (rows_.size()==3)
        ostride = *(rows_.begin()+2);

    integer starting_row = rows_.size()==0 ? 0 : *rows_.begin();
    integer starting_col = cols_.size()==0 ? 0 : *cols_.begin();
    integer end_row      = rows_.size()==0 ? mat.rows() : *(rows_.begin()+1);
    integer end_col      = cols_.size()==0 ? mat.cols() : *(cols_.begin()+1);

#ifndef NDEBUG
    SC_ASSERT(starting_row>=0 && starting_col>=0, "NEGATIVE INDICES NOT ALLOWED");
    SC_ASSERT(end_row<=size(mat,0) && end_col<=size(mat,1), "INDEX OUT OF RANGE");
#endif
    
    integer span_rows = 0; for (auto i=starting_row; i<end_row; i+=ostride) span_rows++;
    integer span_cols = 0; for (auto i=starting_col; i<end_col; i+=istride) span_cols++;

    Eigen::Map<matrix<typename Derived::Scalar>,Eigen::Unaligned,
            Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic> > out(mat.data()+starting_row*mat.cols()+starting_col,
                                                  span_rows, span_cols,
                                                  Eigen::Stride<Eigen::Dynamic,
                                                  Eigen::Dynamic>(ostride*mat.outerStride(),istride*mat.innerStride()));

    return out;
}


template<typename Derived, typename T, typename U>
SC_INLINE auto seq(Derived &mat,
                   T row,
                   std::initializer_list<U> &&cols_)
 -> typename Eigen::Map<matrix<typename Derived::Scalar>,Eigen::Unaligned,
         Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>> {
    /* Array slicing
        A[5,0:5:2]  ---> seq(A,5,{0,5,2})
    */
#ifndef NDEBUG
    SC_ASSERT(cols_.size()<4, "INVALID INDEX FOR MATRIX");
#endif

    integer istride=1,ostride=1;
    if (cols_.size()==3)
        istride = *(cols_.begin()+2);

    integer starting_row = row;
    integer starting_col = cols_.size()==0 ? 0 : *cols_.begin();
    integer end_row      = row+1;
    integer end_col      = cols_.size()==0 ? mat.cols() : *(cols_.begin()+1);

#ifndef NDEBUG
    SC_ASSERT(end_row<=size(mat,0) && end_col<=size(mat,1), "INDEX OUT OF RANGE");
#endif
    
    integer span_rows = 0; for (auto i=starting_row; i<end_row; i+=ostride) span_rows++;
    integer span_cols = 0; for (auto i=starting_col; i<end_col; i+=istride) span_cols++;

    Eigen::Map<matrix<typename Derived::Scalar>,Eigen::Unaligned,
            Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic> > out(mat.data()+starting_row*mat.cols()+starting_col,
                                                  span_rows, span_cols,
                                                  Eigen::Stride<Eigen::Dynamic,
                                                  Eigen::Dynamic>(ostride*mat.outerStride(),istride*mat.innerStride()));

    return out;
}


template<typename Derived, typename T, typename U>
SC_INLINE auto seq(Derived &mat,
                   std::initializer_list<U> &&  rows_,
                   T col)
 -> typename Eigen::Map<matrix<typename Derived::Scalar>,Eigen::Unaligned,
         Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>> {
    /* Array slicing
        A[0:5,3]  ---> seq(A,{0,5},3)
    */
#ifndef NDEBUG
    SC_ASSERT(rows_.size()<4, "INVALID INDEX FOR MATRIX");
#endif

    integer istride=1,ostride=1;
    if (rows_.size()==3)
        ostride = *(rows_.begin()+2);

    integer starting_row = rows_.size()==0 ? 0 : *rows_.begin();
    integer starting_col = col;
    integer end_row      = rows_.size()==0 ? mat.rows() : *(rows_.begin()+1);
    integer end_col      = col+1;

#ifndef NDEBUG
    SC_ASSERT(end_row<=size(mat,0) && end_col<=size(mat,1), "INDEX OUT OF RANGE");
#endif
    
    integer span_rows = 0; for (auto i=starting_row; i<end_row; i+=ostride) span_rows++;
    integer span_cols = 0; for (auto i=starting_col; i<end_col; i+=istride) span_cols++;

    Eigen::Map<matrix<typename Derived::Scalar>,Eigen::Unaligned,
            Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic> > out(mat.data()+starting_row*mat.cols()+starting_col,
                                                  span_rows, span_cols,
                                                  Eigen::Stride<Eigen::Dynamic,
                                                  Eigen::Dynamic>(ostride*mat.outerStride(),istride*mat.innerStride()));

    return out;
}



//------------------------------------------------------------------------------------------------------------------

#endif

SC_INLINE vector<real> rand(integer m) {
    return vector<real>::Random(m);
}
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
SC_INLINE matrix<typename Derived::Scalar> abs(const Eigen::MatrixBase<Derived> &mat) {
    return mat.array().abs();
}
template<typename Derived>
SC_INLINE typename Derived::Scalar sum(const Eigen::DenseBase<Derived> &mat) {
    return mat.sum();
}
template<typename Derived>
SC_INLINE vector<typename Derived::Scalar> sum(const Eigen::DenseBase<Derived> &mat, integer axis) {
    if (axis==0) return mat.rowwise().sum();
    else if (axis==1) return mat.colwise().sum();
    else {
      SC_ASSERT(false,"AXIS ARGUMENT SHOULD BE EITHER 0 OR 1");
      return vector<typename Derived::Scalar>{};
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
SC_INLINE matrix<real> cos(const Eigen::MatrixBase<Derived> &mat) {
    return mat.array().cos();
}
template<typename Derived>
SC_INLINE matrix<real> tan(const Eigen::MatrixBase<Derived> &mat) {
    return mat.array().tan();
}
template<typename Derived>
SC_INLINE matrix<real> sinh(const Eigen::MatrixBase<Derived> &mat) {
    return mat.array().sinh();
}
template<typename Derived>
SC_INLINE matrix<real> cosh(const Eigen::MatrixBase<Derived> &mat) {
    return mat.array().cosh();
}
template<typename Derived>
SC_INLINE matrix<real> tanh(const Eigen::MatrixBase<Derived> &mat) {
    return mat.array().tanh();
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

//template<typename T>
//matrix<T> outer(vector<T> a, vector<T> b) {
//    matrix<T> out = zeros<T>(size(a),size(b));
//    for (integer i=0; i<size(a); ++i) {
//        for (integer j=0; j<size(b); ++j) {
//            out(i,j) += a[i]*b[j];
//        }
//    }
//    return out;
//}
template<typename T>
SC_INLINE matrix<T> outer(vector<T> a, vector<T> b) {
    return a*b.transpose();
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
SC_INLINE Eigen::Matrix<T,1,-1> transpose(const vector<T> &vec) {
    return vec.transpose();
}
template<typename T>
SC_INLINE Eigen::Matrix<T,1,-1> transpose(vector<T> &&vec) {
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
    return mat.stableNorm();
}
template<typename Derived>
SC_INLINE typename Derived::Scalar norm(Eigen::MatrixBase<Derived> &&mat) {
    return mat.stableNorm();
}
// template<typename Derived>
// SC_INLINE typename Derived::Scalar norm(const Eigen::PlainObjectBase<Derived> &mat) {
//     return std::sqrt(std::inner_product(mat.data(),mat.data()+size(mat),mat.data(),0));
// }
// template<typename Derived>
// SC_INLINE typename Derived::Scalar norm(Eigen::PlainObjectBase<Derived> &&mat) {
//     return std::sqrt(std::inner_product(mat.data(),mat.data()+size(mat),mat.data(),0));
// }
// template<typename Derived>
// SC_INLINE typename Derived::Scalar norm(const Eigen::PlainObjectBase<Derived> &mat) {
//     real norm           =  0.;
//     for (auto i=0; i<size(mat); i++){
//         real  cur        =  mat.data()[i];
//         norm            +=  cur*cur;
//     }
//     return  sqrt(norm);
// }

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


static spmatrix<real> block_sparse_extractor(const integer *arr, integer free_size, integer total) {
    vector<integer> I(free_size), J(free_size);
//    for (auto ifree=0; ifree<free_size; ++ifree) {
//        I[ifree] = ifree;
//        J[ifree] = arr[ifree];
//        V[ifree] = 1;
//    }

//    vector<integer> I = arange<integer>(free_size);
//    vector<real>    V = ones(free_size);
    spmatrix<real> b_ones(free_size, total);

    using T = Eigen::Triplet<real>;
    std::vector<T> triplets;
    triplets.reserve(free_size);

    for (auto ifree=0; ifree<free_size; ++ifree) {
        triplets.push_back(T(ifree,arr[ifree],1.));
    }

    b_ones.setFromTriplets(triplets.begin(),triplets.end());
    return b_ones;
//    spmatrix<real>
}


//
template<class ForwardIt, class T>
inline ForwardIt binary_find(ForwardIt first, ForwardIt last, const T& value, std::less<int> comp={}) {
    first = std::lower_bound(first, last, value, comp);
    return first != last && !comp(value, *first) ? first : last;
}

template<typename T>
SC_INLINE int interpolation_search(const T *__restrict__ arr, int len, int key) {
    using int64 = long long int;

    int low = 0;
    int high = len - 1;
    int mid;

    int l = arr[low];
    int h = arr[high];

    while (l <= key && h >= key) {
        int64 high_low = (high - low);
        int64 key_l = (key - l);
        int64 product = high_low*key_l;
        int64 h_l = h-l;
        int64 step = product / h_l;
        mid = low + step;

        int m = arr[mid];

        if (m < key) {
            l = arr[low = mid + 1];
        } else if (m > key) {
            h = arr[high = mid - 1];
        } else {
            return mid;
        }
    }

    if (arr[low] == key)
        return low;
    else
        return len;
}


template<typename T>
spmatrix<T> 
seqNC(const spmatrix<T> &a, const std::vector<int> &krows, const std::vector<int> &kcols) {

#ifndef NDEBUG 
    assert(std::is_sorted(kcols.begin(),kcols.end()) && "ACCESSOR DATA MUST BE SORTED");
#endif

    const T   *__restrict__   valptr = a.valuePtr();
    const int *__restrict__ inptr    = a.innerIndexPtr();
    const int *__restrict__ outptr   = a.outerIndexPtr();
    int nnz                          = a.nonZeros();

    const int *keep_rows             = krows.data();
    const int *keep_cols             = kcols.data();
    int size_keep_rows               = krows.size();
    int size_keep_cols               = kcols.size();

    std::vector<int> o_inptr; std::vector<T> o_valptr;
    o_inptr.reserve(nnz); o_valptr.reserve(nnz);
    std::vector<int> o_outptr(size_keep_rows+1);
    o_outptr[0] = 0;

    const auto keep_cols_0 = keep_cols[0];
    const auto keep_cols_f = keep_cols[size_keep_cols-1];

    auto counter = 0;
    for (auto i=0; i<size_keep_rows; ++i) {
        const auto idx = keep_rows[i];
        auto out_start = outptr[idx];
        auto out_end = outptr[idx+1];  if (idx == a.rows() - 1) out_end = nnz;  
        for (auto j=out_start; j<out_end && inptr[j]<=keep_cols_f; ++j) {       
            if (inptr[j]>=keep_cols_0) {
                // const auto idx_into_keep_cols = binary_find(keep_cols,keep_cols+size_keep_cols,inptr[j]) - keep_cols;
                const auto idx_into_keep_cols = interpolation_search(keep_cols,size_keep_cols,inptr[j]); // 2x faster than simple binary search
                if (idx_into_keep_cols != size_keep_cols) {
                    o_inptr.push_back(keep_cols[idx_into_keep_cols] - keep_cols_0);
                    o_valptr.push_back(valptr[j]);
                    counter++;
                }
            }
        }
        o_outptr[i+1] = counter;
    }

    return Eigen::Map<spmatrix<T>>(size_keep_rows, size_keep_cols, 
        counter, o_outptr.data(), o_inptr.data(), o_valptr.data(),0);
}


template<typename T>
spmatrix<T> 
seqNC(const spmatrix<T> &a, const vector<integer> &krows, const vector<integer> &kcols) {

#ifndef NDEBUG 
    assert(std::is_sorted(kcols.data(),kcols.data()+size(kcols)) && "ACCESSOR DATA MUST BE SORTED");
#endif

    const T   *__restrict__   valptr = a.valuePtr();
    const int *__restrict__ inptr    = a.innerIndexPtr();
    const int *__restrict__ outptr   = a.outerIndexPtr();
    int nnz                          = a.nonZeros();

    const integer *keep_rows             = krows.data();
    const integer *keep_cols             = kcols.data();
    int size_keep_rows               = size(krows);
    int size_keep_cols               = size(kcols);

    std::vector<int> o_inptr; std::vector<T> o_valptr;
    o_inptr.reserve(nnz); o_valptr.reserve(nnz);
    std::vector<int> o_outptr(size_keep_rows+1);
    o_outptr[0] = 0;

    const auto keep_cols_0 = keep_cols[0];
    const auto keep_cols_f = keep_cols[size_keep_cols-1];

    auto counter = 0;
    for (auto i=0; i<size_keep_rows; ++i) {
        const auto idx = keep_rows[i];
        auto out_start = outptr[idx];
        auto out_end = outptr[idx+1];  if (idx == a.rows() - 1) out_end = nnz;  
        for (auto j=out_start; j<out_end && inptr[j]<=keep_cols_f; ++j) {       
            if (inptr[j]>=keep_cols_0) {
                // const auto idx_into_keep_cols = binary_find(keep_cols,keep_cols+size_keep_cols,inptr[j]) - keep_cols;
                const auto idx_into_keep_cols = interpolation_search(keep_cols,size_keep_cols,inptr[j]); // 2x faster than simple binary search
                if (idx_into_keep_cols != size_keep_cols) {
                    o_inptr.push_back(keep_cols[idx_into_keep_cols] - keep_cols_0);
                    o_valptr.push_back(valptr[j]);
                    counter++;
                }
            }
        }
        o_outptr[i+1] = counter;
    }

    return Eigen::Map<spmatrix<T>>(size_keep_rows, size_keep_cols, 
        counter, o_outptr.data(), o_inptr.data(), o_valptr.data(),0);
}


template<typename T>
spmatrix<T> 
seqNC(const spmatrix<T> &a, const std::vector<int> &krows) {

    std::vector<int> kcols(a.cols()); 
    std::iota(kcols.begin(),kcols.end(),0);

#ifndef NDEBUG 
    assert(std::is_sorted(kcols.begin(),kcols.end()) && "ACCESSOR DATA MUST BE SORTED");
#endif

    const T   *__restrict__   valptr = a.valuePtr();
    const int *__restrict__ inptr    = a.innerIndexPtr();
    const int *__restrict__ outptr   = a.outerIndexPtr();
    int nnz                          = a.nonZeros();

    const int *keep_rows             = krows.data();
    const int *keep_cols             = kcols.data();
    int size_keep_rows               = krows.size();
    int size_keep_cols               = kcols.size();

    std::vector<int> o_inptr; std::vector<T> o_valptr;
    o_inptr.reserve(nnz); o_valptr.reserve(nnz);
    std::vector<int> o_outptr(size_keep_rows+1);
    o_outptr[0] = 0;

    const auto keep_cols_0 = keep_cols[0];
    const auto keep_cols_f = keep_cols[size_keep_cols-1];

    auto counter = 0;
    for (auto i=0; i<size_keep_rows; ++i) {
        const auto idx = keep_rows[i];
        auto out_start = outptr[idx];
        auto out_end = outptr[idx+1];  if (idx == a.rows() - 1) out_end = nnz;  
        for (auto j=out_start; j<out_end && inptr[j]<=keep_cols_f; ++j) {       
            if (inptr[j]>=keep_cols_0) {
                // const auto idx_into_keep_cols = binary_find(keep_cols,keep_cols+size_keep_cols,inptr[j]) - keep_cols;
                const auto idx_into_keep_cols = interpolation_search(keep_cols,size_keep_cols,inptr[j]); // 2x faster than simple binary search
                if (idx_into_keep_cols != size_keep_cols) {
                    o_inptr.push_back(keep_cols[idx_into_keep_cols] - keep_cols_0);
                    o_valptr.push_back(valptr[j]);
                    counter++;
                }
            }
        }
        o_outptr[i+1] = counter;
    }

    return Eigen::Map<spmatrix<T>>(size_keep_rows, size_keep_cols, 
        counter, o_outptr.data(), o_inptr.data(), o_valptr.data(),0);
}


template<typename T>
spmatrix<T> 
seqNC(const spmatrix<T> &a, const integer *keep_rows, integer size_keep_rows) {

    std::vector<int> kcols(a.cols()); 
    std::iota(kcols.begin(),kcols.end(),0);

#ifndef NDEBUG 
    assert(std::is_sorted(kcols.begin(),kcols.end()) && "ACCESSOR DATA MUST BE SORTED");
#endif

    const T   *__restrict__   valptr = a.valuePtr();
    const int *__restrict__ inptr    = a.innerIndexPtr();
    const int *__restrict__ outptr   = a.outerIndexPtr();
    int nnz                          = a.nonZeros();

    const int *keep_cols             = kcols.data();
    int size_keep_cols               = kcols.size();

    std::vector<int> o_inptr; std::vector<T> o_valptr;
    o_inptr.reserve(nnz); o_valptr.reserve(nnz);
    std::vector<int> o_outptr(size_keep_rows+1);
    o_outptr[0] = 0;

    const auto keep_cols_0 = keep_cols[0];
    const auto keep_cols_f = keep_cols[size_keep_cols-1];

    auto counter = 0;
    for (auto i=0; i<size_keep_rows; ++i) {
        const auto idx = keep_rows[i];
        auto out_start = outptr[idx];
        auto out_end = outptr[idx+1];  if (idx == a.rows() - 1) out_end = nnz;  
        for (auto j=out_start; j<out_end && inptr[j]<=keep_cols_f; ++j) {       
            if (inptr[j]>=keep_cols_0) {
                // const auto idx_into_keep_cols = binary_find(keep_cols,keep_cols+size_keep_cols,inptr[j]) - keep_cols;
                const auto idx_into_keep_cols = interpolation_search(keep_cols,size_keep_cols,inptr[j]); // 2x faster than simple binary search
                if (idx_into_keep_cols != size_keep_cols) {
                    o_inptr.push_back(keep_cols[idx_into_keep_cols] - keep_cols_0);
                    o_valptr.push_back(valptr[j]);
                    counter++;
                }
            }
        }
        o_outptr[i+1] = counter;
    }

    return Eigen::Map<spmatrix<T>>(size_keep_rows, size_keep_cols, 
        counter, o_outptr.data(), o_inptr.data(), o_valptr.data(),0);
}





template<typename T, typename Derived>
SC_INLINE vector<T> spsolve(const spmatrix<T> &A, const Derived &b, T scale=1.0) {

    // // Scale global matrices
 //    Eigen::IterScaling<spmatrix<real>> scale;
 //    // Compute the left and right scaling vectors. The matrix is equilibrated at output
 //    scale.computeRef(A);
 //    // Scale the right hand side
 //    b = scale.LeftScaling().cwiseProduct(b);

    // Solve
    // Eigen::ConjugateGradient<spmatrix<real>> solver;
    Eigen::SparseLU<spmatrix<T>> solver;
    // Compute the ordering permutation vector from the structural pattern of A
    solver.analyzePattern(A);
    // Compute the numerical factorization
    solver.factorize(A);
    // Use the factors to solve the linear system
    auto sol = solver.solve(b);

    return sol;

    // Scale back the computed solution
    // sol = scale.RightScaling().cwiseProduct(sol);

// #include <unsupported/Eigen/src/IterativeSolvers/Scaling.h>
// #include <Eigen/IterativeLinearSolvers>
// #include <Eigen/SparseCholesky>
// #include <Eigen/UmfPackSupport>

                //     Eigen::SparseLU<spmatrix<real>> solver;
                // // Eigen::UmfPackLU<spmatrix<real>> solver(reducedmatrix);
                // // Eigen::BiCGSTAB<spmatrix<real>> solver;
                // // Eigen::SimplicialLDLT<spmatrix<real>> solver;
                // // Eigen::SimplicialLLT<spmatrix<real>> solver;
                // // Compute the ordering permutation vector from the structural pattern of A
                // solver.analyzePattern(reducedmatrix);
                // // Compute the numerical factorization
                // solver.factorize(reducedmatrix);
                // // Use the factors to solve the linear system
                // // auto sol = solver.solve(-reducedresidual);
                // vector<real> sol = solver.solve(-reducedresidual);
                // // std::cout << "#iterations:     " << solver.iterations() << std::endl;
                // // std::cout << "estimated error: " << solver.error()      << std::endl;


}



} // end of namespace SoftComp


#endif

#endif // BACKEND_EIGEN_H
