#ifndef BACKEND_BUILTIN_H
#define BACKEND_BUILTIN_H

#include "commons/commons.h"
#if defined(HAS_NONE)


namespace SoftComp {

template<typename T, integer SO=1>
class matrix;

// vector class
template<typename T>
class vector {
public:
    using Scalar = T;
    SC_INLINE vector() = default;
    SC_INLINE vector(integer size): _data(size) {}
    SC_INLINE vector(const matrix<T> &a) {
        SC_ASSERT(a.rows()==1 || a.cols()==1, "CANNOT ASSIGN MATRIX TO VECTOR");
    }

    SC_INLINE integer size() const {return _size;}

    SC_INLINE T* data() {return _data.data();}
    SC_INLINE const T* data() const {return _data.data();}

    //! Get the underlying array
    SC_INLINE std::vector<T> array() {return _data;}

    SC_INLINE T& operator[](integer m) {
        return _data[m];
    }
    SC_INLINE T& operator()(integer m) {
        return _data[m];
    }

    SC_INLINE T operator+=(const vector<T>& other) {
        for (auto i=0; i<this->size(); ++i) {
            _data[i] += other.data()[i];
        }
    }
    SC_INLINE T operator-=(const vector<T>& other) {
        for (auto i=0; i<this->size(); ++i) {
            _data[i] *= other.data()[i];
        }
    }
    /* Use multiply instead
    SC_INLINE T operator*=(const vector<T>& other) {
        for (auto i=0; i<this->size(); ++i) {
            _data[i] *= other.data()[i];
        }
    }
    */
    SC_INLINE T operator/=(const vector<T>& other) {
        for (auto i=0; i<this->size(); ++i) {
            _data[i] /= other.data()[i];
        }
    }

private:
    std::vector<T> _data;
    integer _size;
};

// Overloads for vectors
template<typename T>
SC_INLINE integer size(const vector<T> &a) {
    return a.size();
}

template<typename T>
SC_INLINE integer rows(const vector<T> &a) {
    return a.size();
}

template<typename T>
SC_INLINE integer cols(const vector<T> &a) {
    return 1;
}

template<typename T>
SC_INLINE vector<T> zeros(integer m) {
    vector<T> out(m);
    std::fill(out.data(),out.data()+m,0);
    return 0;
}
SC_INLINE vector<real> zeros(integer m) {
    vector<real> out(m);
    std::fill(out.data(),out.data()+m,0);
    return 0;
}

template<typename T>
SC_INLINE vector<T> ones(integer m) {
    vector<T> out(m);
    std::fill(out.data(),out.data()+m,1);
    return 0;
}
SC_INLINE vector<real> ones(integer m) {
    vector<real> out(m);
    std::fill(out.data(),out.data()+m,1);
    return 0;
}

template<typename T>
SC_INLINE vector<T> arange(integer n) {
    vector<T> out(n);
    std::iota(out.data(),out.data()+n,0);
    return out;
}

template<typename T>
SC_INLINE vector<T> linspace(real low, real high, integer num=50) {
    vector<T> out(num);
    // CHCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCEEEECK
    std::for_each(out.data(),out.data()+num,
                  [&out,low,high,num](integer i){out[i]=(high-low)/(real)num;});
    return out;
}

//----------------------------------------------------------------------
//matrix class

/* Builtin dense matrix class
 * SO=1 for RowMajor
 * SO=0 for ColMajor
 */

template<typename T, integer SO>
class matrix {
public:
    using Scalar = T;
    SC_INLINE matrix() = default;
//    SC_INLINE matrix(integer size): _data(n) {}
    SC_INLINE matrix(integer rows, integer cols): _data(rows*cols), _rows(rows), _cols(cols) {}
    SC_INLINE matrix(vector<T> a) {
        _data = std::vector<T>(a.size());
        std::copy(a.data(),a.data()+a.size(),_data.data());
        _rows = a.size(), _cols = 1;
    }

    SC_INLINE integer rows() const {return _rows;}
    SC_INLINE integer cols() const {return _cols;}
    SC_INLINE integer size() const {return _rows*_cols;}

    SC_INLINE T* data() {return _data.data();}
    SC_INLINE const T* data() const {return _data.data();}

    //! Get the underlying array
    SC_INLINE std::vector<T> array() {return _data;}

    template<typename U>
    SC_INLINE T& operator()(U m) {
        SC_ASSERT(cols()==1, "CANNOT INDEX MATRIX WITH ONE INDEX");
        return _data[m];
    }
    template<typename U, typename V>
    SC_INLINE T& operator()(U m, V n) {
        return _data[m*_rows+n];
    }

    SC_INLINE T operator+=(const matrix<T>& other) {
        for (auto i=0; i<this->size(); ++i) {
            _data[i] += other.data()[i];
        }
    }
    SC_INLINE T operator-=(const matrix<T>& other) {
        for (auto i=0; i<this->size(); ++i) {
            _data[i] *= other.data()[i];
        }
    }
    /*
    SC_INLINE T operator*=(const matrix<T>& other) {
        for (auto i=0; i<this->size(); ++i) {
            _data[i] *= other.data()[i];
        }
    }
    */
    SC_INLINE T operator/=(const matrix<T>& other) {
        for (auto i=0; i<this->size(); ++i) {
            _data[i] /= other.data()[i];
        }
    }

    SC_INLINE vector<T>& row(integer m) {
        vector<T> out;
        for (integer i=0; i<this->cols(); ++i)
            out[i] = _data[m*_rows+i];
        return out;

    }
    SC_INLINE vector<T> col(integer m) {
        vector<T> out;
        for (integer i=0; i<this->rows(); ++i)
            out[i] = _data[i*_rows+m];
        return out;
    }

private:
    std::vector<T> _data;
    integer _rows, _cols;
};

template<typename T>
SC_INLINE integer size(const matrix<T> &a) {
    return a.rows()*a.cols();
}

template<typename T>
SC_INLINE integer rows(const matrix<T> &a) {
    return a.rows();
}

template<typename T>
SC_INLINE integer cols(const matrix<T> &a) {
    return a.cols();
}

template<typename T>
SC_INLINE matrix<T> zeros(integer m, integer n) {
    matrix<T> out(m,n);
    std::fill(out.data(),out.data()+m*n,0);
    return out;
}
SC_INLINE matrix<real> zeros(integer m, integer n) {
    matrix<real> out(m);
    std::fill(out.data(),out.data()+m*n,0);
    return out;
}

template<typename T>
SC_INLINE matrix<T> arange(integer m, integer n) {
    matrix<T> out(n-m,1);
    std::iota(out.data(),out.data()+n-m,m);
    return out;
}


// Generic operator overloaders
//--------------------------------------------------------------------------------------
template<typename MatType>
SC_INLINE MatType operator+(const MatType &a, const MatType &b) {
    MatType out(a.size());
    using T = typename MatType::Scalar;
    const T* a_data = a.data();
    const T* b_data = b.data();
    T* out_data = out.data();
    for (auto i=0; i<a.size(); ++i) {
        out_data[i] = a_data[i] + b_data[i];
    }
    return out;
}
template<typename MatType, typename U>
SC_INLINE MatType operator+(const MatType &a, U num) {
    MatType out(a.size());
    using T = typename MatType::Scalar;
    const T* a_data = a.data();
    T* out_data = out.data();
    for (auto i=0; i<a.size(); ++i) {
        out_data[i] = a_data[i] + num;
    }
    return out;
}
template<typename MatType, typename U>
SC_INLINE MatType operator+(U num, const MatType &a) {
    MatType out(a.size());
    using T = typename MatType::Scalar;
    const T* a_data = a.data();
    T* out_data = out.data();
    for (auto i=0; i<a.size(); ++i) {
        out_data[i] = num + a_data[i];
    }
    return out;
}

template<typename MatType>
SC_INLINE MatType operator-(const MatType &a, const MatType &b) {
    MatType out(a.size());
    using T = typename MatType::Scalar;
    const T* a_data = a.data();
    const T* b_data = b.data();
    T* out_data = out.data();
    for (auto i=0; i<a.size(); ++i) {
        out_data[i] = a_data[i] - b_data[i];
    }
    return out;
}
template<typename MatType, typename U>
SC_INLINE MatType operator-(const MatType &a, U num) {
    MatType out(a.size());
    using T = typename MatType::Scalar;
    const T* a_data = a.data();
    T* out_data = out.data();
    for (auto i=0; i<a.size(); ++i) {
        out_data[i] = a_data[i] - num;
    }
    return out;
}
template<typename MatType, typename U>
SC_INLINE MatType operator-(U num, const MatType &a) {
    MatType out(a.size());
    using T = typename MatType::Scalar;
    const T* a_data = a.data();
    T* out_data = out.data();
    for (auto i=0; i<a.size(); ++i) {
        out_data[i] = num - a_data[i];
    }
    return out;
}

template<typename MatType, typename U,
         typename std::enable_if<std::is_arithmetic<U>::value,boolean>::type=0>
SC_INLINE MatType operator*(const MatType &a, U num) {
    MatType out(a.size());
    using T = typename MatType::Scalar;
    const T* a_data = a.data();
    T* out_data = out.data();
    for (auto i=0; i<a.size(); ++i) {
        out_data[i] = a_data[i] * num;
    }
    return out;
}
template<typename MatType, typename U,
         typename std::enable_if<std::is_arithmetic<U>::value,boolean>::type=0>
SC_INLINE MatType operator*(U num, const MatType &a) {
    MatType out(a.size());
    using T = typename MatType::Scalar;
    const T* a_data = a.data();
    T* out_data = out.data();
    for (auto i=0; i<a.size(); ++i) {
        out_data[i] = num * a_data[i];
    }
    return out;
}

template<typename MatType, typename U>
SC_INLINE MatType multiply(const MatType &a, U num) {
    MatType out(a.size());
    using T = typename MatType::Scalar;
    const T* a_data = a.data();
    T* out_data = out.data();
    for (auto i=0; i<a.size(); ++i) {
        out_data[i] = a_data[i] * num;
    }
    return out;
}
template<typename MatType, typename U>
SC_INLINE MatType multiply(U num, const MatType &a) {
    MatType out(a.size());
    using T = typename MatType::Scalar;
    const T* a_data = a.data();
    T* out_data = out.data();
    for (auto i=0; i<a.size(); ++i) {
        out_data[i] = num * a_data[i];
    }
    return out;
}

template<typename MatType>
SC_INLINE MatType operator/(const MatType &a, const MatType &b) {
    MatType out(a.size());
    using T = typename MatType::Scalar;
    const T* a_data = a.data();
    const T* b_data = b.data();
    T* out_data = out.data();
    for (auto i=0; i<a.size(); ++i) {
        out_data[i] = a_data[i] / b_data[i];
    }
    return out;
}
template<typename MatType, typename U>
SC_INLINE MatType operator/(const MatType &a, U num) {
    MatType out(a.size());
    using T = typename MatType::Scalar;
    const T* a_data = a.data();
    T* out_data = out.data();
    for (auto i=0; i<a.size(); ++i) {
        out_data[i] = a_data[i] / num;
    }
    return out;
}
template<typename MatType, typename U>
SC_INLINE MatType operator/(U num, const MatType &a) {
    MatType out(a.size());
    using T = typename MatType::Scalar;
    const T* a_data = a.data();
    T* out_data = out.data();
    for (auto i=0; i<a.size(); ++i) {
        out_data[i] = num / a_data[i];
    }
    return out;
}


// Math
template<typename MatType>
SC_INLINE MatType abs(const MatType &a) {
    MatType out(a.size());
    using T = typename MatType::Scalar;
    const T* a_data = a.data();
    T* out_data = out.data();
    for (auto i=0; i<a.size(); ++i) {
        out_data[i] = std::abs(a_data[i]);
    }
}

template<typename MatType>
SC_INLINE typename MatType::Scalar max(const MatType &a) {
    MatType out(a.size());
    using T = typename MatType::Scalar;
    const T* a_data = a.data();
    T* out_data = out.data();
    return std::max(a.data().begin(),a.data().end());
}

template<typename MatType>
SC_INLINE typename MatType::Scalar min(const MatType &a) {
    MatType out(a.size());
    using T = typename MatType::Scalar;
    const T* a_data = a.data();
    T* out_data = out.data();
    return std::min(a.data().begin(),a.data().end());
}

template<typename MatType, typename U>
SC_INLINE MatType pow(const MatType &a, U num) {
    MatType out(a.size());
    using T = typename MatType::Scalar;
    const T* a_data = a.data();
    T* out_data = out.data();
    for (auto i=0; i<a.size(); ++i) {
        out_data[i] = std::pow(a_data[i],num);
    }
    return out;
}

template<typename MatType, typename U>
SC_INLINE MatType pow(U num, const MatType &a) {
    MatType out(a.size());
    using T = typename MatType::Scalar;
    const T* a_data = a.data();
    T* out_data = out.data();
    for (auto i=0; i<a.size(); ++i) {
        out_data[i] = std::pow(num, a_data[i]);
    }
    return out;
}


template<typename Derived, typename std::enable_if<!std::is_arithmetic<Derived>::value,bool>::type=0>
SC_INLINE void fill(Derived &mat, integer num) {
    //! In-place
    std::fill(mat.data(),mat.data()+size(mat),num);
}

template<typename T, typename U>
SC_INLINE vector<T> append(vector<T> &vec, U num) {
    vector<T> out(vec.size()+1);
    std::copy(vec.data(),vec.data()+vec.size(),out.data());
    out[vec.size()] = static_cast<T>(num);
    return out;
}

template<typename MatType>
SC_INLINE boolean isempty(const MatType &mat) {
    return mat.rows()==0 && mat.cols()==0;
}

template<typename T>
SC_INLINE matrix<T> fliplr(const matrix<T> &mat) {
    //! Reverse a matrix left to right
    SC_ASSERT(mat.cols()!=1,"CANNOT FLIP A COLUMN VECTOR COLUMN-WISE");
    matrix<T> out(mat);
    if (mat.rows()==1) {
        std::reverse(out.data(),out.data()+out.size());
        return out;
    }
    else {
        for (auto i=0; i<mat.rows(); ++i) {
            std::reverse(out.row(i).data(),out.row(i).data()+out.cols());
        }
    }
    return out;
}
template<typename T>
SC_INLINE vector<T> fliplr(const vector<T> &vec) {
    //! Reverse a vector
    vector<T> out(vec);
    std::reverse(out.data(),out.data()+out.size());
    return out;
}

template<typename T>
SC_INLINE matrix<T> flipud(const matrix<T> &mat) {
    //! Reverse a matrix left to right
    SC_ASSERT(mat.rows()!=1,"CANNOT FLIP A ROW VECTOR ROW-WISE");
    matrix<T> out(mat);
    if (mat.cols()==1) {
        std::reverse(out.data(),out.data()+out.size());
        return out;
    }
    else {
        for (auto i=0; i<mat.cols(); ++i) {
            std::reverse(out.col(i).data(),out.col(i).data()+out.rows());
        }
    }
    return out;
}
template<typename T>
SC_INLINE vector<T> flipud(const vector<T> &vec) {
    //! Reverse a vector
    vector<T> out(vec);
    std::reverse(out.data(),out.data()+out.size());
    return out;
}

template<typename MatType0, typename MatType1>
SC_INLINE matrix<typename MatType0::Scalar> vstack(const MatType0 &a, const MatType1 &b) {
    SC_ASSERT(rows(a)==rows(b),"NO VSTACK");
    matrix<typename MatType0::Scalar> out(rows(a)+rows(b),cols(a));
    std::copy(a.data(),a.data()+a.size(),out.data());
    std::copy(b.data(),b.data()+b.size(),out.data()+a.size());
    return out;
}

template<typename MatType0, typename MatType1>
SC_INLINE matrix<typename MatType0::Scalar> hstack(const MatType0 &a, const MatType1 &b) {
    SC_ASSERT(cols(a)==cols(b),"NO HSTACK");
    matrix<typename MatType0::Scalar> out(rows(a),cols(a)+cols(b));
    for (integer i=0; i<rows(a); ++i) {
        for (integer j=0; j<cols(a); ++j) {
            out(i,j) = a(i,j);
        }
    }
    std::copy(a.data(),a.data()+a.size(),out.data());
    std::copy(b.data(),b.data()+b.size(),out.data()+a.size());
    return out;
}

}

#endif

#endif // BACKEND_BUILTIN_H
