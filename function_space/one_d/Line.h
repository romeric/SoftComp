#ifndef LINE_H
#define LINE_H

#include "quadrarure_rules/QuadratureRule.h"

namespace SoftComp {


class Line {
public:
    vector<real> bases;
    vector<real> grad_bases;

    Line() = default;

    SC_INLINE Line(integer C, real xi, boolean equally_spaced=false) {

        auto n = C+2;
        auto quad = QuadratureRule();
        vector<real> eps;
        if (!equally_spaced) {quad.GaussLobatto(C+1); eps = quad.points;}
        else {quad.Gauss(C+1,-1.,1.); eps = quad.points;}

        auto A = zeros(n,n);
        for (auto i=0; i<n; ++i) {
            A(i,0) = 1.;
        }

        for (auto j=0; j<n; ++j)
            for (auto i=1; i<n; ++i)
                A(j,i) = pow(eps[j],i);

        bases  = zeros(n);
        grad_bases = zeros(n);

        for (auto ishape=0; ishape<n; ishape++) {
            auto RHS = zeros(n);
            RHS[ishape] = 1.;

            // Solve linear system (dense LU)
            vector<real> coeff = solve(A,RHS);
            // Build shape functions
            for (auto incr=0; incr<n; incr++)
                bases[ishape] += coeff[incr]*pow(xi,incr);

            // Build derivate of shape functions
            for (auto incr=0; incr<n-1; incr++)
                grad_bases[ishape] += (incr+1)*coeff[incr+1]*pow(xi,incr);
        }


    }
};


#include "mesh/nodal_arrangement.h"


class Quad {
public:
    vector<real> bases;
    matrix<real> grad_bases;

    Quad() = default;

    SC_INLINE Quad(integer C, real zeta, real eta, boolean equally_spaced=false, boolean arrange=true) {

        // Allocate
        auto nsize = (C+2)*(C+2);
        bases = zeros(nsize);
        grad_bases = zeros(nsize,2);
        // Compute bases
        auto line_func_zeta = Line(C, zeta, equally_spaced);
        auto line_func_eta  = Line(C, eta, equally_spaced);
        auto &Nzeta = line_func_zeta.bases;
        auto &Neta = line_func_eta.bases;
        // Compute gradient of bases functions
        auto &gNzeta = line_func_zeta.grad_bases;
        auto &gNeta = line_func_eta.grad_bases;

        // Don't pass directly product expressions to functions receiving eigen generic expressions
//        matrix<real> dum0 = Neta*transpose(Nzeta);
//        vector<real> dum1 = flatten(dum0);
//        matrix<real> dum2 = gNeta*transpose(Nzeta);
//        matrix<real> dum3 = Neta*transpose(gNzeta);
//        vector<real> dum4 = flatten(dum2);
//        vector<real> dum5 = flatten(dum3);

        vector<real> dum1 = flatten(outer(Neta,Nzeta));
        vector<real> dum4 = flatten(outer(gNeta,Nzeta));
        vector<real> dum5 = flatten(outer(Neta,gNzeta));
        // Ternsorial product
        if (arrange==1) {
            auto node_arranger = NodalArrangement("quad",C+1);

            for (auto i=0; i< size(bases); ++i) {
                bases(i) = dum1(node_arranger.element_arrangement(i));
                grad_bases(i,0) = dum4(node_arranger.element_arrangement(i));
                grad_bases(i,1) = dum5(node_arranger.element_arrangement(i));
            }
        }
        else {
            bases.head(nsize) = dum1;
            grad_bases.col(0) = dum4;
            grad_bases.col(1) = dum5;
        }
    }

};

}

#endif // LINE_H
