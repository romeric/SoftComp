#ifndef BACKEND_BLAZE_H
#define BACKEND_BLAZE_H

#if defined(HAS_BLAZE)

#include "commons/commons.h"


namespace SoftComp {

//template<typename Derived, boolean SO>
//SC_INLINE integer size(const blaze::Matrix<Derived,SO> &mat) {
//    return mat.rows()*mat.columns();
//}
//template<typename Derived, boolean SO>
//SC_INLINE integer size(blaze::Matrix<Derived,SO> &&mat) {
//    return mat.rows()*mat.columns();
//}

//template<typename Derived, boolean SO>
//SC_INLINE integer size(const blaze::Matrix<Derived,SO> &mat, integer idx) {
//    return idx==0 ? mat.rows() : mat.columns();
//}
//template<typename Derived, boolean SO>
//SC_INLINE integer size(blaze::Matrix<Derived,SO> &&mat, integer idx) {
//    return idx==0 ? mat.rows() : mat.columns();
//}

//template<typename Derived, boolean SO>
//SC_INLINE integer rows(const blaze::Matrix<Derived,SO> &mat) {
//    return mat.columns();
//}
//template<typename Derived, boolean SO>
//SC_INLINE integer rows(blaze::Matrix<Derived,SO> &&mat) {
//    return mat.columns();
//}

//template<typename Derived, boolean SO>
//SC_INLINE integer cols(const blaze::DenseMatrix<Derived,SO> &mat) {
//    return mat.columns();
//}
//template<typename Derived, boolean SO>
//SC_INLINE integer cols(blaze::Matrix<Derived,SO> &&mat) {
//    return mat.columns();
//}

//template<typename Derived>
//SC_INLINE integer size(const Derived &mat) {
//    return mat.rows()*mat.columns();
//}
template<typename T>
SC_INLINE integer size(const matrix<T> &mat) {
    return mat.rows()*mat.columns();
}
template<typename T>
SC_INLINE integer size(const vector<T> &mat) {
    return mat.size();
}

template<typename Derived>
SC_INLINE integer rows(const Derived &mat) {
    return mat.rows();
}
template<typename T>
SC_INLINE integer rows(const vector<T> &mat) {
    return mat.size();
}

template<typename Derived>
SC_INLINE integer cols(const Derived &mat) {
    return mat.columns();
}
template<typename T>
SC_INLINE integer cols(const vector<T> &) {
    return 1;
}

//SC_INLINE matrix<real> rand(integer m, integer n) {
//    return matrix<real>::Random(m,n);
//}

template<typename T>
SC_INLINE vector<T> zeros(integer m) {
    vector<T> out(m); std::fill(out.data(),out.data()+out.size(),0);
    return out;
}
template<typename T>
SC_INLINE matrix<T> zeros(integer m, integer n) {
    matrix<T> out(m,n); std::fill(out.data(),out.data()+size(out),0);
    return out;
}
SC_INLINE vector<real> zeros(integer m) {
    vector<real> out(m); std::fill(out.data(),out.data()+out.size(),0);
    return out;
}
SC_INLINE matrix<real> zeros(integer m, integer n) {
    matrix<real> out(m,n); std::fill(out.data(),out.data()+size(out),0);
    return out;
}

template<typename T>
SC_INLINE vector<T> ones(integer m) {
    vector<T> out(m); std::fill(out.data(),out.data()+out.size(),1);
    return out;
}
template<typename T>
SC_INLINE matrix<T> ones(integer m, integer n) {
    matrix<T> out(m,n); std::fill(out.data(),out.data()+size(out),1);
    return out;
}
SC_INLINE vector<real> ones(integer m) {
    vector<real> out(m); std::fill(out.data(),out.data()+out.size(),1);
    return out;
}
SC_INLINE matrix<real> ones(integer m, integer n) {
    matrix<real> out(m,n); std::fill(out.data(),out.data()+size(out),1);
    return out;
}


template<typename Derived, typename std::enable_if<!std::is_arithmetic<Derived>::value,bool>::type=0>
SC_INLINE void fill(Derived &mat, integer num) {
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
SC_INLINE vector<T> arange(integer m, integer n) {
    vector<T> out(n-m);
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
    vector<real> out(num);
    for (integer i=0; i<num; ++i) {
        out[i] = low + (high-low)/(real)num;
    }
    return out;
}


template<typename Derived>
SC_INLINE boolean isempty(const Derived &mat) {
    return mat.rows()==0 && mat.columns()==0;
}

template<typename T, typename U>
SC_INLINE vector<T>
append(const vector<T> &arr, U num) {
    //! Append to the end of Eigen vector (similar to push_back)

    vector<T> new_arr(arr.size()+1);
    new_arr = arr;
    new_arr[arr.size()] = static_cast<T>(num);
    return new_arr;
}


//template<typename Derived>
//SC_INLINE matrix<typename Derived::ElementType> flipud(const Derived &mat) {
//    //! Reverse a matrix top to bottom
//    SC_ASSERT(rows(mat)!=1,"CANNOT FLIP A ROW VECTOR ROW-WISE");
//    matrix<typename Derived::ElementType> out(rows(mat),cols(mat));
//    if (cols(mat)==1) {
//        std::reverse(mat.data(),mat.data()+size(mat));
//        return mat;
//    }
//    else {
//        for (auto i=0; i<cols(mat); ++i) {
//            out.col(i) = col(mat,i).reverse();
//        }
//    }
//    return out;
//}

template<typename T>
SC_INLINE vector<T> flipud(const vector<T> &mat) {
    //! Reverse a matrix top to bottom
    std::reverse(mat.data(),mat.data()+size(mat));
    return mat;
}

template<typename Derived, boolean SO>
SC_INLINE vector<typename Derived::ElementType> flipud(const blaze::DenseVector<Derived,SO> &mat) {
    //! Reverse a matrix top to bottom
    vector<typename Derived::ElementType> out(mat);
    std::reverse(out.data(),out.data()+size(mat));
    return mat;
}


template<typename T>
SC_INLINE matrix<T> hstack(const matrix<T> &mat0, const matrix<T> &mat1) {
    SC_ASSERT(mat0.rows()==mat1.rows(),"CONCATENATING MATRICES SHOULD HAVE THE SAME NUMBER OF ROWS");
    matrix<T> out(mat0.rows(),mat0.columns()+mat1.columns());
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
template<typename T>
SC_INLINE matrix<T> hstack(const vector<T> &mat0, const vector<T> &mat1) {
    SC_ASSERT(rows(mat0)==rows(mat1),"CONCATENATING MATRICES SHOULD HAVE THE SAME NUMBER OF ROWS");
    matrix<T> out(rows(mat0),cols(mat0)+cols(mat1));
    for (integer i=0; i<rows(mat0); ++i) {
        out(i,0) = mat0[i];
        out(i,1) = mat1[i];
    }
    return out;
}

template<typename T>
SC_INLINE matrix<T> vstack(const matrix<T> &mat0, const matrix<T> &mat1) {
    SC_ASSERT(cols(mat0)==cols(mat1),"CONCATENATING MATRICES SHOULD HAVE THE SAME NUMBER OF COLUMNS");
    matrix<T> out(rows(mat0)+rows(mat1),cols(mat0));
    std::copy(mat0.data(),mat0.data()+size(mat0),out.data());
    std::copy(mat1.data(),mat1.data()+size(mat1),out.data()+size(mat0));
    return out;
}
template<typename T>
SC_INLINE vector<T> vstack(const vector<T> &mat0, const vector<T> &mat1) {
    vector<T> out(rows(mat0)+rows(mat1));
    std::copy(mat0.data(),mat0.data()+size(mat0),out.data());
    std::copy(mat1.data(),mat1.data()+size(mat1),out.data()+size(mat0));
    return out;
}

//template<typename Derived0, typename Derived1>
//SC_INLINE matrix<typename Derived0::Scalar>
//vstack(const Eigen::PlainObjectBase<Derived0> &mat0, const Eigen::PlainObjectBase<Derived1> &mat1) {
//    SC_ASSERT(mat0.cols()==mat1.cols(),"CONCATENATING MATRICES SHOULD HAVE THE SAME NUMBER OF COLUMNS");
//    matrix<typename Derived0::Scalar> out(mat0.rows()+mat1.rows(),mat0.cols());
//    std::copy(mat0.data(),mat0.data()+size(mat0),out.data());
//    std::copy(mat1.data(),mat1.data()+size(mat1),out.data()+size(mat0));
//    return out;
//}



template<typename T>
SC_INLINE vector<T> col(const matrix<T> &mat, integer idx) {
    return blaze::column(mat,idx);
}

template<typename T>
SC_INLINE vector<T> row(const matrix<T> &mat, integer idx) {
    return blaze::row(mat,idx);
}




template<typename Derived>
SC_INLINE auto max(Derived &&mat) -> decltype(blaze::max(mat)) {
    return blaze::max(mat);
}





// Math
//template<typename Derived, boolean SO, typename std::enable_if<!std::is_arithmetic<Derived>::value,boolean>::type=0>
//SC_INLINE matrix<typename Derived::ElementType> abs(const blaze::DenseMatrix<Derived,SO> &mat) {
//    return blaze::abs(mat);
//}
//template<typename Derived, boolean SO, typename std::enable_if<!std::is_arithmetic<Derived>::value,boolean>::type=0>
//SC_INLINE vector<typename Derived::ElementType> abs(const blaze::DenseVector<Derived,SO> &mat) {
//    return blaze::abs(mat);
//}
template<typename Derived, typename std::enable_if<!std::is_arithmetic<Derived>::value,boolean>::type=0>
SC_INLINE auto abs(const Derived &mat) -> decltype(blaze::abs(mat)) {
    return blaze::abs(mat);
}

template<typename Derived, typename T, typename std::enable_if<!std::is_arithmetic<Derived>::value,boolean>::type=0>
SC_INLINE auto pow(const Derived &mat, T num) -> decltype(blaze::pow(mat,num)) {
    return blaze::pow(mat,num);
}

template<typename Derived, typename std::enable_if<!std::is_arithmetic<Derived>::value,boolean>::type=0>
SC_INLINE auto sin(const Derived &mat) -> decltype(blaze::sin(mat)) {
    return blaze::sin(mat);
}

template<typename Derived, typename std::enable_if<!std::is_arithmetic<Derived>::value,boolean>::type=0>
SC_INLINE auto cos(const Derived &mat) -> decltype(blaze::cos(mat)) {
    return blaze::cos(mat);
}


// Linear algebra subroutines
template<typename Derived, typename T, boolean SO,
         typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0>
SC_INLINE matrix<T> //blaze::DenseMatrix<Derived,SO>
operator+(const blaze::DenseMatrix<Derived,SO> &mat0, T a) {
    auto B = blaze::forEach(mat0, [&a](T d) { return a+d; } );
    return B;
}
template<typename Derived, typename T, boolean SO,
         typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0>
SC_INLINE vector<T> //blaze::DenseVector<Derived,SO>
operator+(const blaze::DenseVector<Derived,SO> &mat0, T a) {
    auto B = blaze::forEach(mat0, [&a](T d) { return a+d; } );
    return B;
}
template<typename Derived, typename T, boolean SO,
         typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0>
SC_INLINE matrix<T> //blaze::DenseMatrix<Derived,SO>
operator+(T a, const blaze::DenseMatrix<Derived,SO> &mat0) {
    auto B = blaze::forEach(mat0, [&a](T d) { return a+d; } );
    return B;
}
template<typename Derived, typename T, boolean SO,
         typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0>
SC_INLINE vector<T> //blaze::DenseVector<Derived,SO>
operator+(T a, const blaze::DenseVector<Derived,SO> &mat0) {
    auto B = blaze::forEach(mat0, [&a](T d) { return a+d; } );
    return B;
}

template<typename Derived, typename T, boolean SO,
         typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0>
SC_INLINE matrix<T> //blaze::DenseMatrix<Derived,SO>
operator-(const blaze::DenseMatrix<Derived,SO> &mat0, T a) {
    auto B = blaze::forEach(mat0, [&a](T d) { return d-a; } );
    return B;
}
template<typename Derived, typename T, boolean SO,
         typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0>
SC_INLINE vector<T> //blaze::DenseVector<Derived,SO>
operator-(const blaze::DenseVector<Derived,SO> &mat0, T a) {
    auto B = blaze::forEach(mat0, [&a](T d) { return d-a; } );
    return B;
}
template<typename Derived, typename T, boolean SO,
         typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0>
SC_INLINE matrix<T> //blaze::DenseMatrix<Derived,SO>
operator-(T a, const blaze::DenseMatrix<Derived,SO> &mat0) {
    auto B = blaze::forEach(mat0, [&a](T d) { return a-d; } );
    return B;
}
template<typename Derived, typename T, boolean SO,
         typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0>
SC_INLINE vector<T> //blaze::DenseVector<Derived,SO>
operator-(T a, const blaze::DenseVector<Derived,SO> &mat0) {
    auto B = blaze::forEach(mat0, [&a](T d) { return a-d; } );
    return B;
}

template<typename Derived, typename T, boolean SO,
         typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0>
SC_INLINE matrix<T> //blaze::DenseMatrix<Derived,SO>
multiply(const blaze::DenseMatrix<Derived,SO> &mat0, T a) {
    auto B = blaze::forEach(mat0, [&a](T d) { return a*d; } );
    return B;
}
//template<typename Derived, typename T, boolean SO,
//         typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0>
//SC_INLINE vector<T> //blaze::DenseVector<Derived,SO>
//multiply(const blaze::DenseVector<Derived,SO> &mat0, T a) {
//    auto B = blaze::forEach(mat0, [&a](T d) { return a*d; } );
//    return B;
//}
template<typename Derived, typename T, boolean SO,
         typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0>
SC_INLINE matrix<T> //blaze::DenseMatrix<Derived,SO>
multiply(T a, const blaze::DenseMatrix<Derived,SO> &mat0) {
    auto B = blaze::forEach(mat0, [&a](T d) { return a*d; } );
    return B;
}
//template<typename Derived, typename T, boolean SO,
//         typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0>
//SC_INLINE vector<T> //blaze::DenseVector<Derived,SO>
//multiply(T a, const blaze::DenseVector<Derived,SO> &mat0) {
//    auto B = blaze::forEach(mat0, [&a](T d) { return a*d; } );
//    return B;
//}

//template<typename Derived0, typename Derived1, typename T, boolean SO>
//SC_INLINE vector<typename Derived0::ElementType> //blaze::DenseMatrix<Derived,SO>
//multiply(const blaze::DenseVector<Derived0,SO> &mat0, const blaze::DenseVector<Derived1,SO> &mat1) {
//    vector<typename Derived0::ElementType> out;
//    for (integer i=0; i<size(mat0); ++i)
//        out[i] = mat0[i]*mat1[i];
//}
template<typename Derived0, typename Derived1, boolean SO>
SC_INLINE vector<typename Derived0::ElementType> //blaze::DenseMatrix<Derived,SO>
multiply(const blaze::DenseVector<Derived0,SO> &mat0, const blaze::DenseVector<Derived1,SO> &mat1) {
    vector<typename Derived0::ElementType> matt0(mat0);
    vector<typename Derived0::ElementType> matt1(mat1);
    vector<typename Derived0::ElementType> out;
    for (integer i=0; i<size(mat0); ++i)
        out[i] = matt0[i]*matt1[i];
    return out;
}



//template<typename Derived, typename T, boolean SO,
//         typename std::enable_if<std::is_arithmetic<T>::value && !std::is_arithmetic<Derived>::value,bool>::type=0>
//SC_INLINE matrix<T> //blaze::DenseMatrix<Derived,SO>
//operator/(const blaze::DenseMatrix<Derived,SO> &mat0, T a) {
//    auto B = blaze::forEach(mat0, [&a](T d) { return d/a; } );
//    return B;
//}
//template<typename Derived, typename T, boolean SO,
//         typename std::enable_if<std::is_arithmetic<T>::value && !std::is_arithmetic<Derived>::value,bool>::type=0>
//SC_INLINE vector<T> //blaze::DenseVector<Derived,SO>
//operator/(const blaze::DenseVector<Derived,SO> &mat0, T a) {
//    auto B = blaze::forEach(mat0, [&a](T d) { return d/a; } );
//    return B;
//}
//template<typename Derived, typename T, boolean SO,
//         typename std::enable_if<std::is_arithmetic<T>::value && !std::is_arithmetic<Derived>::value,bool>::type=0>
//SC_INLINE matrix<T> //blaze::DenseMatrix<Derived,SO>
//operator/(T a, const blaze::DenseMatrix<Derived,SO> &mat0) {
//    auto B = blaze::forEach(mat0, [&a](T d) { return a/d; } );
//    return B;
//}
//template<typename Derived, typename T, boolean SO,
//         typename std::enable_if<std::is_arithmetic<T>::value && !std::is_arithmetic<Derived>::value,bool>::type=0>
//SC_INLINE vector<T> //blaze::DenseVector<Derived,SO>
//operator/(T a, const blaze::DenseVector<Derived,SO> &mat0) {
//    auto B = blaze::forEach(mat0, [&a](T d) { return a/d; } );
//    return B;
//}

//template<typename T, typename U>
//SC_INLINE vector<U>
//operator/(T a, const vector<U> &mat0) {
//    auto B = blaze::forEach(mat0, [&a](U d) { return a/d; } );
//    return B;
//}

template<typename T, typename Derived,
         typename std::enable_if<std::is_arithmetic<T>::value && !std::is_arithmetic<Derived>::value,bool>::type=0>
SC_INLINE vector<T>
operator/(T a, const Derived &mat0) {
    auto B = blaze::forEach(mat0, [&a](T d) { return a/d; } );
    return B;
}


}
#endif

#endif // BACKEND_BLAZE_H
