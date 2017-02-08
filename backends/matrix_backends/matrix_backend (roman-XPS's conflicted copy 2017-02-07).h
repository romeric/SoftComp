#ifndef MATRIX_BACKEND_H
#define MATRIX_BACKEND_H

#include <set>
#include <initializer_list>
#include <random>

namespace SoftComp {

enum linear_solver {non_symmetric, symmetric, symmetric_positive_definite, semi_indefinite};

// Aux
template<typename T>
void unsort_unique(std::vector<T>& vec)
{
    std::set<T> seen;
    auto newEnd = std::remove_if(vec.begin(), vec.end(), [&seen](const T& value)
    {
        if (seen.find(value) != std::end(seen))
            return true;
        seen.insert(value);
        return false;
    });
    vec.erase(newEnd, vec.end());
}
//
}

#if defined(HAS_EIGEN)
#include "backend_eigen.h"
#elif defined(HAS_ARMA)
#include "backend_arma.h"
#include "backend_blaze.h"
#else
#include "backend_builtin.h"
#endif

namespace SoftComp {

// For primitive types
template<typename T, typename U, typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0>
SC_INLINE T pow(T a, U p) {return std::pow(a,p);}

template<typename T, typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0> SC_INLINE T log(T a) {return std::log(a);}
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0> SC_INLINE T log10(T a) {return std::log10(a);}
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0> SC_INLINE T log1p(T a) {return std::log1p(a);}
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0> SC_INLINE T exp(T a) {return std::exp(a);}
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0> SC_INLINE T abs(T a) {return std::abs(a);}

#ifndef __SSE2__
template<typename T> SC_INLINE T sqrt(T a) {return std::sqrt(a);}
#else
template<typename T,
         typename std::enable_if<std::is_same<T,double>::value,bool>::type=0>
SC_INLINE T sqrt(T a) {
    return _mm_cvtsd_f64( _mm_sqrt_pd(_mm_set1_pd(a)));
}
template<typename T,
         typename std::enable_if<std::is_same<T,float>::value,bool>::type=0>
SC_INLINE T sqrt(T a) {
    return _mm_cvtss_f32( _mm_sqrt_ps(_mm_set1_ps(a)));
}
#endif

// Trig
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0> SC_INLINE T sin(T a) {return std::sin(a);}
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0> SC_INLINE T cos(T a) {return std::cos(a);}
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0> SC_INLINE T tan(T a) {return std::tan(a);}
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0> SC_INLINE T sinh(T a) {return std::sinh(a);}
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0> SC_INLINE T cosh(T a) {return std::cosh(a);}
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value,bool>::type=0> SC_INLINE T tanh(T a) {return std::tanh(a);}



/*---------------------------------------------------------------------------------------------------------------------------*/
/* Common impls */
SC_INLINE matrix<integer> randint(integer low, integer high, integer rows, integer cols) {
    std::random_device                      rd;
    std::mt19937                            generator(rd());
    std::uniform_int_distribution<integer>  distr(low, high);

    matrix<integer> out(rows,cols);
    for (auto i=0; i<rows; ++i) {
        for (auto j=0; j<cols; ++j) {
            out(i,j) = distr(generator);
        }
    }
    return out;
}

SC_INLINE vector<integer> randint(integer low, integer high, integer size) {
    std::random_device                      rd;
    std::mt19937                            generator(rd());
    std::uniform_int_distribution<integer>  distr(low, high);

    vector<integer> out(size);
    for (auto i=0; i<size; ++i) {
        out[i] = distr(generator);
    }
    return out;
}










template <typename T>
SC_INLINE
std::tuple<vector<T>,vector<integer>,vector<integer>>
unique_all(const vector<T> &a) {
    auto size_a = size(a);
    std::vector<T> uniques;
    std::vector<integer> idx;
    vector<integer> inv(size_a);

    for (auto i=0; i<size_a; ++i) {
        auto counter = 0;
        for (auto j=0; j<uniques.size(); ++j) {
            if (uniques[j]==a[i]) {
                counter +=1;
                break;
            }
        }
        if (counter==0) {
            uniques.push_back(a[i]);
            idx.push_back(i);
        }
    }

    for (auto i=0; i<size_a; ++i) {
        for (auto j=0; j<uniques.size(); ++j) {
            if (uniques[j]==a[i]) {
                inv[i] = j;
                break;
            }
        }
    }

    auto size_u = uniques.size();
    vector<T> uniques_e(size_u);
    vector<integer> idx_e(size_u);
#if defined(HAS_EIGEN)
    std::copy(idx.begin(),idx.end(),idx_e.data());
    std::copy(uniques.begin(),uniques.end(),uniques_e.data());
#elif defined(HAS_ARMA)
    std::copy(idx.begin(),idx.end(),idx_e.memptr());
    std::copy(uniques.begin(),uniques.end(),uniques_e.memptr());
#endif

    return std::make_tuple(uniques_e,idx_e,inv);
}


#include <map>

template <typename T>
inline std::tuple<vector<T>,
                  vector<integer>,
                  vector<integer>>
unique_all_2(const vector<T> &a)
 {
   integer               ind;
   std::map<T, integer>  m;
   std::vector<T>        uniques;
   std::vector<integer>  idx;
   std::vector<integer>  inv; inv.reserve(a.size());

   ind = 0;
   for(integer i=0; i <size(a); ++i) {
      auto e = m.insert(std::make_pair(a[i], ind));
      if (e.second) {
         uniques.push_back(a[i]);
         idx.push_back(i);
         ++ind;
      }
      inv.push_back(e.first->second);
    }

   auto size_u = uniques.size();
   vector<T> uniques_e(size_u);
   vector<integer> idx_e(size_u);
   vector<integer> inv_e(a.size());
#if defined(HAS_EIGEN)
   std::copy(idx.begin(),idx.end(),idx_e.data());
   std::copy(uniques.begin(),uniques.end(),uniques_e.data());
   std::copy(inv.begin(),inv.end(),inv_e.data());
#elif defined(HAS_ARMA)
   std::copy(idx.begin(),idx.end(),idx_e.memptr());
   std::copy(uniques.begin(),uniques.end(),uniques_e.memptr());
   std::copy(inv.begin(),inv.end(),inv_e.memptr());
#endif

   return std::make_tuple(uniques_e,idx_e,inv_e);
}
/*---------------------------------------------------------------------------------------------------------------------------*/

}


#endif // MATRIX_BACKEND_H
