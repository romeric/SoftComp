#ifndef GRADIENTSPHASEFIELD_H
#define GRADIENTSPHASEFIELD_H


#include "backends/matrix_backends/backend_tensor.h"
#include "backends/matrix_backends/matrix_backend.h"


namespace SoftComp {

/*-------------------------------------------------------------------
 * ------------------------------------------------------------------
 * ------------------------------------------------------------------
This function computes the kinematcis
* -------------------------------------------------------------------
* -------------------------------------------------------------------
* -------------------------------------------------------------------
*/

class GradientsPhaseField {
public:
    tensor<real,1> s;
    tensor<real,2> material_gradient_s;
    tensor<real,3> material_gradient;

    GradientsPhaseField() = delete;

    GradientsPhaseField(const Kinematics &kin, const tensor<real,1> &s_elem, const QuadratureRule &quadrule, const ConstitutiveModel &model) {
        auto ngauss            =  size(grad_bases,0);
        auto ndim              =  size(grad_bases,2);
        auto nodeperelem       =  size(s_elem,0);
        //------------------------------------------
        //  Material gradient of shape functions
        //------------------------------------------
        material_gradient      =  tensor<real,3>(ngauss,nodeperelem,ndim);
        material_gradient_s    =  tensor<real,3>(ngauss,nodeperelem,ndim);
        s                      =  tensor<real,3>(ngauss,ndim,ndim);

        std::array<tpair, 1> idx0 = {tpair(0, 0)};
        std::array<tpair, 1> idx1 = {tpair(1, 0)};
        for (auto igauss=0;igauss<ngauss;igauss++){
            s(igauss)                          =  tensor2scalar(static_cast<tensor<real,0>>(s_elem.contract(quadrule.bases.chip(igauss,0),idx0)));
            material_gradient.chip(igauss,0)   =  quadrule.grad_bases.chip(igauss,0).contract(kin.inv_dX_chi.chip(igauss,0),idx1);
            material_gradient_s.chip(igauss,0) =  s_elem.contract(material_gradient.chip(igauss,0),idx0);
        }
    }
};

}

#endif // GRADIENTSPHASEFIELD_H
