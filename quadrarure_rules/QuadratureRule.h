#ifndef QUADRATURERULE_H
#define QUADRATURERULE_H

#include "backends/matrix_backends/matrix_backend.h"
#include "function_space/one_d/Jacobi.h"

namespace SoftComp {

class QuadratureRule {
public:
    real degree; // Corresponds to the interpolation degree
    real order;  // Corresponds to the quadrature order used for a given polynomial degree
    vector<real> points;
    vector<real> weights;

    QuadratureRule() = default;

#if !defined(HAS_BLAZE)

    void Gauss(integer N, real a=-1., real b=1.) {

        auto N0 = N-1;
        auto N1 = N0+1;
        auto N2 = N0+2;

        // Construct points on [-1,1] first
        auto xu         = linspace(-1.,1.,N1);
        // Legendre-Gauss-Vandermonde Matrix
        auto L          = zeros(N1,N2);
        // Derivative of Legendre-Gauss-Vandermonde Matrix
        auto Lp         = zeros(N1);
        auto dum        = linspace(0.0,N0,N1);
        vector<real> y  = cos((2.*dum+1.)*PI/(2.*N0+2.))+(0.27/N1)*sin(PI*xu*N0/N2);
        // Initial Guess
        auto y0         = fill<real>(N1,2);

        while (max(abs(y-y0)) > EPS) {

            L.col(0) = ones(N1);
            L.col(1) = y;

            for (auto k=1; k!=N1; k++) {
                L.col(k+1) = ((2*k+1)*multiply(L.col(k),y)-k*L.col(k-1))/(k+1);
            }

            Lp = N2*(L.col(N0)-multiply(L.col(N1),y))/(1.-pow(y,2));
            y0 = y;
            y  = y0 - L.col(N1)/Lp;
        }

        // Gauss Points
        this->points = flipud((a*(1-y)+b*(1+y))/2.);
        // Gauss Weights
        auto coeff = static_cast<real>(N2)/static_cast<real>(N1);
        this->weights = flipud((b-a)/(multiply((1-pow(y,2)),multiply(Lp,Lp)))*(coeff*coeff));
    }

    /*
    void Gauss(integer N, real a=-1., real b=1.) {

        auto N0 = N-1;
        auto N1 = N0+1;
        auto N2 = N0+2;

        // Construct points on [-1,1] first
        auto xu              = linspace(-1.,1.,N1);
        // Legendre-Gauss-Vandermonde Matrix
        auto L               = zeros(N1,N2);
        // Derivative of Legendre-Gauss-Vandermonde Matrix
        std::vector<real>    Lp(N1);
        auto dum             = linspace(0.0,N0,N1);
        vector<real> y       = cos((2.*dum+1.)*PI/(2.*N0+2.))+(0.27/N1)*sin(PI*xu*N0/N2);
        // Initial Guess
        auto y0              = fill<real>(N1,2);

        while (max(abs(y-y0)) > EPS) {
            for (auto i=0; i<N1; ++i) {
                L(i,0) = 1.;
                L(i,1) = y[i];
            }
            for (auto k=1; k!=N1; k++) {
                for (auto i=0; i<N1; ++i) {
                    L(i,k+1) = ((2.*k+1)*L(i,k)*y[i]-k*L(i,k-1))/(k+1);
                }
            }
            for (auto i=0; i<N1; ++i) {
                Lp[i] = N2*(L(i,N0)-L(i,N1)*y[i])/(1.-y[i]*y[i]);
                y0[i] = y[i];
                y[i]  = y0[i] - L(i,N1)/Lp[i];
            }
        }

        // Gauss Points
        points  = zeros(N1);
        weights = zeros(N1);
        for (auto i=N1-1; i>=0; --i) {
            points[N1-i-1] = (a*(1-y[i])+b*(1+y[i]))/2.;
        }
        // Gauss Weights
        auto coeff  = (real)N2/(real)N1;
        auto coeff0 = coeff*coeff;
        auto coeff1 = b-a;
        for (auto i=N1-1; i>=0; --i) {
            weights[N1-i-1] = coeff1/((1.-(y[i]*y[i]))*Lp[i]*Lp[i])*coeff0;
        }
    }
    */

#else
    void Gauss(integer N, real a=-1., real b=1.) {

        auto N0 = N-1;
        auto N1 = N0+1;
        auto N2 = N0+2;

        // Construct points on [-1,1] first
        auto xu = linspace(-1.,1.,N1);
        // Legendre-Gauss-Vandermonde Matrix
        auto L = zeros(N1,N2);
        // Derivative of Legendre-Gauss-Vandermonde Matrix
        auto Lp = zeros(N1);
        auto dum = linspace(0.0,N0,N1);
        vector<real> y = cos((2.*dum+1.)*PI/(2.*N0+2.))+(0.27/N1)*sin(PI*xu*N0/N2);
        // Initial Guess
        auto y0 = fill<real>(N1,2);

        while (max(abs(y-y0)) > EPS) {

            col(L,0) = ones(N1);
            Lp = zeros(N1);
            col(L,1) = y;

            for (auto k=1; k!=N1; k++) {
                col(L,k+1) = ((2*k+1)*multiply(col(L,k),y)-k*col(L,k-1))/(k+1);
            }

            Lp = N2*(col(L,N0)-multiply(col(L,N1),y))/(1.-pow(y,2));
            y0 = y;
            y  = y0 - col(L,N1)/Lp;
        }

        // Gauss Points
        this->points = flipud((a*(1-y)+b*(1+y))/2.);
        // Gauss Weights
        auto coeff = static_cast<real>(N2)/static_cast<real>(N1);
        this->weights = flipud((b-a)/(multiply((1-pow(y,2)),multiply(Lp,Lp)))*(coeff*coeff));
    }
#endif



    SC_INLINE void GaussLobatto(integer N) {

        N = N-1;
//        SC_ASSERT(N!=0,"THERE SHOULD BE ATLEAST TWO GAUSS LOBATTO POINTS");

//        if (N==0) {
//            points = zeros(2);
//            points(0,0) = -1.; points(1,0) = 1.;
//            weights = ones(2);
//            return;
//        }
//        else if (N==1) {
//            points = zeros(3);
//            weights = zeros(3);
//            points(0,0) = -1.; points(1,0) = 0.; points(2,0) = 1.;
//            weights(0,0) = 0.333333333333333315; weights(1,0) = 1.333333333333333259;
//            weights(2,0) = 0.333333333333333315;
//            return;
//        }
//        else if (N==2) {
//            points = zeros(4);
//            weights = zeros(4);
//            points(0,0) = -1.; points(1,0) = -0.447213595499957983;
//            points(2,0) = 0.447213595499957983; points(3,0) = 1.;
//            weights(0,0) = 0.166666666666666657; weights(1,0) = 0.83333333333333337 ;
//            weights(2,0) = 0.83333333333333337 ; weights(3,0) = 0.166666666666666657;
//            return;
//        }

        if (N==0) {
            points = zeros(2);
            points[0] = -1.; points[1] = 1.;
            weights = ones(2);
            return;
        }
        else if (N==1) {
            points = zeros(3);
            weights = zeros(3);
            points[0] = -1.; points[1] = 0.; points[2] = 1.;
            weights[0] = 0.333333333333333315; weights[1] = 1.333333333333333259;
            weights[2] = 0.333333333333333315;
            return;
        }
        else if (N==2) {
            points = zeros(4);
            weights = zeros(4);
            points[0] = -1.; points[1] = -0.447213595499957983;
            points[2] = 0.447213595499957983; points[3] = 1.;
            weights[0] = 0.166666666666666657; weights[1] = 0.83333333333333337 ;
            weights[2] = 0.83333333333333337 ; weights[3] = 0.166666666666666657;
            return;
        }


//        real a=1., b=1., x1;
//        auto jacobi = JacobiPolynomial();
//        // Initial Guess - Chebyshev-Gauss-Lobatto points
//        matrix<real> x = -cos(linspace(0,N+1,N+2)/(real)(N+1)*M_PI);
//        auto nrows = rows(x);
//        // Allocate space for points and weights
//        points  = zeros(nrows,1);
//        weights = zeros(nrows,1);
//        auto y  = zeros(nrows,1);
//        for (auto k=0; k<nrows; ++k) {
//            jacobi.SetDegree(N);
//            auto x0 = x(k,0);
//            real dell = 2.;
//            while (abs(dell) > EPS) {
//                // Polynomial Deflation: Exclude already determined roots
//                real s = 0.;
//                for (auto i=0; i<k; ++i) {
//                    s += 1.0/(x0-points(i,0));
//                }
//                // Compute Jacobi polynomial p(a,b)
//                jacobi.ComputeBases(x0,a,b);
//                // Compute derivative of Jacobi polynomial p(a,b)
//                jacobi.ComputeGradients(x0,a,b,1);
//                // Gauss-Lobatto points are roots of (1-x^2)*dp, hence
//                auto coeff = (1.0 - x0*x0);
//                auto nom = coeff*jacobi.bases(N,0);
//                auto dnom = -2.0*x0*jacobi.bases(N,0)+coeff*jacobi.grad_bases(N,0);
//                dell = - nom/(dnom-nom*s);
//                x1 = x0+dell;
//                x0 = x1;
//            }
//            points(k,0) = x1;
//            // Compute weights
//            jacobi.SetDegree(N+1);
//            jacobi.ComputeBases(x1,0.,0.);
//            weights(k,0) = 2.0/((N+1)*(N+2)*jacobi.bases(N+1,0)*jacobi.bases(N+1,0));
//        }


        real a=1., b=1., x1;
        auto jacobi = JacobiPolynomial();
        // Initial Guess - Chebyshev-Gauss-Lobatto points
        vector<real> x = -cos(linspace(0,N+1,N+2)/(real)(N+1)*PI);
        auto nrows = rows(x);
        // Allocate space for points and weights
        points  = zeros(nrows);
        weights = zeros(nrows);
        auto y  = zeros(nrows);
        for (auto k=0; k<nrows; ++k) {
            jacobi.SetDegree(N);
            auto x0 = x[k];
            real dell = 2.;
            while (abs(dell) > EPS) {
                // Polynomial Deflation: Exclude already determined roots
                real s = 0.;
                for (auto i=0; i<k; ++i) {
                    s += 1.0/(x0-points[i]);
                }
                // Compute Jacobi polynomial p(a,b)
                jacobi.ComputeBases(x0,a,b);
                // Compute derivative of Jacobi polynomial p(a,b)
                jacobi.ComputeGradients(x0,a,b,1);
                // Gauss-Lobatto points are roots of (1-x^2)*dp, hence
                auto coeff = (1.0 - x0*x0);
                auto nom = coeff*jacobi.bases(N,0);
                auto dnom = -2.0*x0*jacobi.bases(N,0)+coeff*jacobi.grad_bases(N,0);
                dell = - nom/(dnom-nom*s);
                x1 = x0+dell;
                x0 = x1;
            }
            points[k] = x1;
            // Compute weights
            jacobi.SetDegree(N+1);
            jacobi.ComputeBases(x1,0.,0.);
            weights[k] = 2.0/((N+1)*(N+2)*jacobi.bases(N+1,0)*jacobi.bases(N+1,0));
        }
    }
};

}

#endif // QUADRATURERULE_H
