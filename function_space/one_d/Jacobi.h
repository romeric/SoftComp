#ifndef JACOBI_H
#define JACOBI_H

#include "backends/matrix_backends/matrix_backend.h"

namespace SoftComp {

class JacobiPolynomial {
public:
    integer n;
    matrix<real> bases;
    matrix<real> grad_bases;

    JacobiPolynomial() = default;
    JacobiPolynomial(integer degree) : n(degree) {}

    SC_INLINE void SetDegree(integer degree) {
        this->n = degree;
    }

    SC_INLINE void ComputeBases(real xi, real a=0., real b=0.)
    {
        bases = zeros(n+1,1);
        real a1n,a2n,a3n,a4n;

        bases(0,0)=1.0;
        if (n>0) {
            bases(1,0) = 0.5*((a-b)+(a+b+2)*xi);
        }
        if (n>1) {
            for (integer p=1; p<n; p++) {
                a1n = 2*(p+1)*(p+a+b+1)*(2*p+a+b);
                a2n = (2*p+a+b+1)*(a*a-b*b);
                a3n = (2*p+a+b)*(2*p+a+b+1)*(2*p+a+b+2);
                a4n = 2*(p+a)*(p+b)*(2*p+a+b+2);

                bases(p+1,0) = ((a2n+a3n*xi)*bases(p,0)-a4n*bases(p-1,0))/a1n;
            }
        }
    }


    SC_INLINE void ComputeGradients(real xi, real a=0., real b=0., integer opt=0)
    {
        SC_ASSERT(n!=0,"POLYNOMIAL DEGREE NOT SET");
        // Creating another instance is necessary to avoid Grad not to modify Bases
        JacobiPolynomial other(n);
        grad_bases = zeros(n+1,1);
        if (opt==1) {
            other.ComputeBases(xi,a+1,b+1);
        }
        else {
            other.ComputeBases(xi,a,b);
        }

        for (integer p=1; p<n+1; p++) {
            grad_bases(p,0) = 0.5*(a+b+p+1)*other.bases(p-1,0);
        }
    }
};






}

#endif // JACOBI_H
