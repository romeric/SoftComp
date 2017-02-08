#ifndef DISPLACEMENTFORMULATION_H
#define DISPLACEMENTFORMULATION_H

#include "VariationalFormulation.h"

namespace SoftComp {

class DisplacementFormulation : public VariationalFormulation {
public:
    tensor<real,2> internal_traction;
    tensor<real,2> local_stiffness;

    void LocalResidual(const QuadratureRule &quadrule, const Kinematics &kin, const ConstitutiveModel &material) {

        auto ndim = kin.F.dimension(0);
        auto ngauss = kin.F.dimension(2);
        auto nodeperelem = kin.material_gradient.dimension(0);



        internal_traction = tensor<real,2>(ndim,nodeperelem);
        internal_traction.setZero();

        for (integer igauss=0; igauss<ngauss; ++igauss) {
            const tensor<real,2> &current_piola         = material.piola.chip(igauss,2);
            const tensor<real,2> &current_material_gradient = kin.material_gradient.chip(igauss,2);
            std::cout << current_material_gradient.dimension(0) << "\n";
//            const tensor<real,2> &current_material_gradient = kin.material_gradient.chip(igauss,2);

//            std::array<Eigen::IndexPair<integer>, 1> idx = { Eigen::IndexPair<integer>(0, 0) };
//            internal_traction += current_piola.contract(current_material_gradient)*quadrule.weights(igauss)*kin.material_iso_jacobian(igauss);
        }
    }

    void LocalStiffness(const QuadratureRule &quadrule, const Kinematics &kin, const ConstitutiveModel &material) {

        auto ndim = kin.F.dimension(0);
        auto ngauss = kin.F.dimension(2);
        auto nodeperelem = kin.material_gradient.dimension(0);

        local_stiffness = tensor<real,2>(ndim*nodeperelem,ndim*nodeperelem);
        local_stiffness.setZero();

        for (integer igauss=0; igauss<ngauss; ++igauss) {
            const tensor<real,2> &current_hessian = material.hessian.chip(igauss,4);
            const tensor<real,2> &current_material_gradient = kin.material_gradient.chip(igauss,2);

            std::array<Eigen::IndexPair<integer>, 1> idx = { Eigen::IndexPair<integer>(3, 1) };
            auto temp = current_hessian.contract(current_material_gradient,idx);
            idx = { Eigen::IndexPair<integer>(1, 1) };
//            local_stiffness += temp.contract(current_material_gradient)*quadrule.weights(igauss)*kin.material_iso_jacobian(igauss);
        }
    }

};

}

#endif // DISPLACEMENTFORMULATION_H
