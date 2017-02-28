#ifndef FUNCTIONSPACE_H
#define FUNCTIONSPACE_H

#include "function_space/one_d/Jacobi.h"
#include "function_space/one_d/Line.h"
#include "backends/matrix_backends/backend_tensor.h"

namespace SoftComp {

class FunctionSpace {
    /*!
    * FunctionSpace class, responsible for generating all kinds
    * of bases functions (i.e. for scalar & vector spaced variables)
    * for all kinds of elements (tris,quads,tets,hexes)
    */
public:
    integer degree;
    integer ndim;
    integer ngauss;
    string etype;
    tensor<real,2> bases;
    tensor<real,3> grad_bases;
    matrix<real> quadrature_points;
    vector<real> quadrature_weights;

    FunctionSpace() = default; /*Default constructor*/
    FunctionSpace(string etype, integer poly_degree, boolean equally_spaced=true) {
        /*! Custom Constructor
        *    input:
        *       etype:              [string] type of element (can be either "tri","quad","tet" or "hex")
        *       poly_degree:        [integer] polynomial degree (p) for the bases
        *       equally_spaced      [boolean] generate equally spaced/Gauss Lobatto/Fekete nodal Lagrangian bases
        *                           default is equally spaced
        */

        this->etype = etype;
        this->degree = poly_degree;
        this->_tol_ = 1e-10;

        auto quadrature_rule = QuadratureRule();
        quadrature_rule.Gauss(degree+1);
        if (etype == "quad") {
            this->ndim = 2;
            this->ngauss = size(quadrature_rule.points)*size(quadrature_rule.points);
            integer nsize = (degree+1)*(degree+1);
            quadrature_points = matrix<real>(ngauss,2);
            quadrature_weights = vector<real>(ngauss);
            auto counter = 0;
            for (auto i=0; i<size(quadrature_rule.points); ++i) {
                for (auto j=0; j<size(quadrature_rule.points); ++j) {
                    quadrature_points(counter,0) = quadrature_rule.points(i);
                    quadrature_points(counter,1) = quadrature_rule.points(j);
                    quadrature_weights(counter)  = quadrature_rule.weights(i)*quadrature_rule.weights(j);
                    counter++;
                }
            }
#ifndef NDEBUG
            // Sanity checks
            SC_ASSERT(abs(sum(quadrature_weights)-4)<_tol_,"QUADRATURE WEIGHTS DO NOT SUM UP TO THE AREA OF THE ELEMENT");
#endif
            auto etype_func = Quad();
            bases = tensor<real,2>(ngauss,nsize);
            grad_bases = tensor<real,3>(ngauss,nsize,2);
            // Evaluate bases at all integration points
            for (auto i=0; i<ngauss; ++i) {
                // Quad basis functions take continuity not degree as an input (i.e. C = degree -1)
                etype_func.Compute(degree-1,quadrature_points(i,0),quadrature_points(i,1));
#ifndef NDEBUG
                // Sanity checks
                SC_ASSERT(abs(sum(etype_func.bases)-1)<_tol_,"BASES DO NOT SATISFY PARTITION OF UNITY");
                SC_ASSERT(abs(sum(sum(etype_func.grad_bases,0)))<_tol_,"X DERIVATES OF BASES DO NOT SUM UP TO ZERO");
                SC_ASSERT(abs(sum(sum(etype_func.grad_bases,1)))<_tol_,"Y DERIVATES OF BASES DO NOT SUM UP TO ZERO");
#endif
                std::copy(etype_func.bases.data(),etype_func.bases.data()+nsize,bases.data()+i*nsize);
                std::copy(etype_func.grad_bases.data(),etype_func.grad_bases.data()+nsize*2,grad_bases.data()+i*nsize*2);
            }
        }

        else{
            SC_ASSERT(false, "BASIS FUNCTIONS FOR "+etype+" NOT YET IMPLEMENTED");
        }
    }

private:
    real _tol_;
};


}

#endif // FUNCTIONSPACE_H
