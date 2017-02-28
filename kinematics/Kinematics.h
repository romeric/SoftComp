#ifndef KINEMATICS_H
#define KINEMATICS_H

#include "backends/matrix_backends/backend_tensor.h"
#include "backends/matrix_backends/matrix_backend.h"


namespace SoftComp {

class Kinematics {
    /*-------------------------------------------------------------------
     * ------------------------------------------------------------------
     * ------------------------------------------------------------------
    Kinematics class responsible for computing all kinematics quantity
    * -------------------------------------------------------------------
    * -------------------------------------------------------------------
    * -------------------------------------------------------------------
    */
public:
    tensor<real,3> F, H;
    tensor<real,1> J, material_iso_jacobian;
    tensor<real,3> material_gradient, inv_dX_chi;
    tensor<real,3> C, b;

    Kinematics() = default;

    Kinematics(const tensor<real,2> &X, const tensor<real,2> &x, const tensor<real,3> &grad_bases) {
        auto ngauss                =  size(grad_bases,0);
        auto ndim                  =  size(X,1);
        auto nodeperelem           =  size(X,0);
        //------------------------------------------
        //  Material gradient of shape functions
        //------------------------------------------
        material_gradient          =  tensor<real,3>(ngauss,nodeperelem,ndim);
        F                          =  tensor<real,3>(ngauss,ndim,ndim);
        H                          =  tensor<real,3>(ngauss,ndim,ndim);
        J                          =  tensor<real,1>(ngauss);
        material_iso_jacobian      =  tensor<real,1>(ngauss);
        inv_dX_chi                 =  tensor<real,3>(ngauss,ndim,ndim);

        std::array<tpair, 1> idx0  = {tpair(0, 0)};
        std::array<tpair, 1> idx1  = {tpair(1, 0)};
        std::array<tpair, 2> idx2  = {tpair(0, 0), tpair(1, 1)};
        for (auto igauss=0;igauss<ngauss;igauss++){
            tensor<real,2> dX_chi             =  X.contract(grad_bases.chip(igauss,0),idx0);
            inv_dX_chi.chip(igauss,0)         =  inv(dX_chi);
            material_gradient.chip(igauss,0)  =  grad_bases.chip(igauss,0).contract(inv_dX_chi.chip(igauss,0),idx1);
            F.chip(igauss,0)                  =  x.contract(material_gradient.chip(igauss,0),idx0);
            H.chip(igauss,0)                  =  cofactor(static_cast<tensor<real,2>>(F.chip(igauss,0)));
            J(igauss)                         =  tensor2scalar(static_cast<tensor<real,0>>
                                                 (F.chip(igauss,0).contract(H.chip(igauss,0),idx2)/(real)ndim));
            material_iso_jacobian(igauss)     =  abs(det(dX_chi));
        }
    }
};


}

#endif // KINEMATICS_H
