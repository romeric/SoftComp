#ifndef EXAMPLES_H
#define EXAMPLES_H

#include "backends/matrix_backends/matrix_backend.h"
#include "commons/print.h"

namespace SoftComp {


//matrix<real> example_1(real a) {
//    matrix<real> b = ones(3)*a;
////    matrix<real> c = b*transpose(b);
////    print(ones(3));
//    matrix<real> c  =  b*transpose(b);
//    print(c);
////    return ones(3,3);
//    return c;
//}


template<typename T>
matrix<T> example_1(T a) {
    matrix<T> b = ones<T>(3)*a;
    matrix<T> c  =  b*transpose(b);
//    print(c);
//    return ones(3,3);
    return c;
}

}

#endif // EXAMPLES_H
