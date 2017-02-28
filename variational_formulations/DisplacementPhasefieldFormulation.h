#ifndef DISPLACEMENTPHASEFIELDFORMULATION_H
#define DISPLACEMENTPHASEFIELDFORMULATION_H

#endif // DISPLACEMENTPHASEFIELDFORMULATION_H

#include "VariationalFormulation.h"

namespace SoftComp {

template<class MaterialType>
class DisplacementPhasefieldFormulation : public VariationalFormulation<MaterialType> {
public:
    tensor<real,2> local_internal_traction;
    tensor<real,2> local_stiffness;
    tensor<real,1> s;
    tensor<real,2> X, x;


    // Kinematics kinematics;

    DisplacementPhasefieldFormulation(integer ndim) {
        nvar = ndim;
    }

    tensor<real,2> LocalResidualDisplacement(const QuadratureRule &quadrule, 
        const Kinematics &kin, const ConstitutiveModel &material) {
        //----------------------------------------------------
        // Internal residual for linear momentum
        //----------------------------------------------------
        auto ndim                    =  kin.F.dimension(1);
        auto ngauss                  =  kin.F.dimension(0);
        auto nodeperelem             =  kin.material_gradient.dimension(1);

        tensor<real,2>  local_T_x(nodeperelem,ndim);
        local_Tx.setZero();

        std::array<tpair, 1> idx0    = {tpair(1,1)};
        for (integer igauss=0; igauss<ngauss; ++igauss) {
            const tensor<real,2> &current_piola              =  material.piola.chip(igauss,0);
            const tensor<real,2> &current_material_gradient  =  kin.material_gradient.chip(igauss,0);
            local_Tx                +=  current_material_gradient.contract(current_piola,idx0)*
                                         quadrule.displacement.weights(igauss)*kin.material_iso_jacobian(igauss);
         }
        return local_Tx;
    }
    tensor<real,4> LocalStiffnessKxx(const QuadratureRule &quadrule, 
        const Kinematics &kin, const ConstitutiveModel &material) {
        //----------------------------------------------------
        // Stiffness contribution Kxx
        //----------------------------------------------------
        auto ndim          =  kin.F.dimension(1);
        auto ngauss        =  kin.F.dimension(0);
        auto nodeperelem   =  kin.material_gradient.dimension(0);

        tensor<real,4> Kxx(nodeperelem,nodeperelem,ndim,ndim,);
        Kxx.setZero();

        std::array<tpair, 1> idx0 = {tpair(1,3)};

        for (integer igauss=0; igauss<ngauss; ++igauss) {
            const tensor<real,4> &current_hessian           = material.hessian.chip(igauss,0);
            const tensor<real,2> &current_material_gradient = kin.material_gradient.chip(igauss,0);

            auto temp      =  current_material_gradient.contract(current_hessian,idx0);
            Kxx           +=  current_material_gradient.contract(temp,idx0)*
                                quadrule.displacement.weights(igauss)*kin.material_iso_jacobian(igauss);
        }
        return Kxx;
    }
    tensor<real,3> LocalStiffnessKxs(const QuadratureRule &quadrule, const Kinematics &kin, 
        const GradientsPhasefield &grad, const tensor<real,1> &s, const ConstitutiveModel &material) {
        //----------------------------------------------------
        // Stiffness contribution Kxs
        //----------------------------------------------------
        auto ndim          =  kin.F.dimension(1);
        auto ngauss        =  kin.F.dimension(0);
        auto nodeperelemx  =  kin.material_gradient.dimension(0);
        auto nodeperelems  =  grad.material_gradient.dimension(0);

        Kxs                =  tensor<real,3>(nodeperelemx,ndim,nodeperelems);  // Not sure yet
        Kxs.setZero();
        // Unfinished
        return Kxs;
    }
    tensor<real,3> LocalStiffnessKsx(const QuadratureRule &quadrule, const Kinematics &kin, 
        const GradientsPhasefield &grad, const tensor<real,1> &s, const ConstitutiveModel &material) {
        //----------------------------------------------------
        // Stiffness contribution Ksx
        //----------------------------------------------------
        auto ndim          =  kin.F.dimension(1);
        auto ngauss        =  kin.F.dimension(0);
        auto nodeperelemx  =  kin.material_gradient.dimension(0);
        auto nodeperelems  =  grad.material_gradient.dimension(0);

        Ksx                =  tensor<real,3>(nodeperelems,nodeperelemx,ndim);  // Not sure yet
        Ksx.setZero();
        // Unfinished
        return Ksx;
    }
    tensor<real,2> LocalStiffnessKss(const QuadratureRule &quadrule, const GradientsPhasefield &grad, 
        const tensor<real,1> &s, const ConstitutiveModel &material) {
        //----------------------------------------------------
        // Stiffness contribution Ksx
        //----------------------------------------------------
        auto ngauss        =  kin.F.dimension(0);
        auto nodeperelems  =  grad.material_gradient.dimension(0);

        Kss                =  tensor<real,2>(nodeperelems,nodeperelems);
        Kss.setZero();
        // Unfinished
        return Kss;
    }


    tensor<real,1> LocalResidualPhasefield(const QuadratureRule &quadrule, 
        const GradientsPhasefield &grad, const ConstitutiveModel &material) {
        //----------------------------------------------------
        // Residual for the phase-field equation
        //----------------------------------------------------
        auto ngauss                  =  grad.material_gradient_s.dimension(0);
        auto nodeperelem             =  grad.material_gradient_s.dimension(0);

        tensor<real,2>  local_T_s(nodeperelem,ndim);
        local_Ts.setZero();

        //  Unfinished

        return local_Ts;

    }

public:
    void LocalResidual(const QuadratureRule &quadrule, const Kinematics &kin, 
        const ConstitutiveModel &material, const GradientsPhasefield &grad) {
        //----------------------------------------------------
        // Appropriate merge of both contributions in the residual,
        //namely linear momentum and the phase-field equations
        //----------------------------------------------------
        tensor<real,2>  local_internal_traction(nodeperelem,ndim);
        local_internal_traction.setZero();

        tensor<real,2> local_Tx  =  LocalResidualDisplacement(quadrule.displacement, kin, material);
        tensor<real,1> local_Ts  =  LocalResidualPhasefield(quadrule.phasefield, grad, material);

        // Unfinished

    }
    void LocalStiffness(const QuadratureRule &quadrule, const Kinematics &kin, const ConstitutiveModel &material) {
        //----------------------------------------------------
        // Internal stiffness in displacement-based formulation
        //----------------------------------------------------
        auto ndim         =  kin.F.dimension(1);
        auto ngauss       =  kin.F.dimension(0);
        auto nodeperelem  =  kin.material_gradient.dimension(0);

        local_stiffness   =  tensor<real,4>(nodeperelem,ndim,ndim,nodeperelem);
        local_stiffness.setZero();

        tensor<real,4> local_Kxx  =  LocalStiffnessKxx(quadrule, kin, material);
        tensor<real,3> local_Kxs  =  LocalStiffnessKxs(quadrule, kin, grad, s, material);
        tensor<real,3> local_Ksx  =  LocalStiffnessKsx(quadrule, kin, grad, s, material);
        tensor<real,2> local_Kss  =  LocalStiffnessKss(quadrule, grad, s, material);

        //  Unfinished
    }


    void Extractor () {}
    void LocalStiffness(const FunctionSpace &function_space, const Mesh &mesh, 
        const Material &material, const Unknowns &unk, integer &ielem=None) {

    }


    void InitialiseVariables(const Mesh &mesh) {
         nnode    =  mesh.phasefield.nnode;
         tensor<real,1> s(nnode);
         s.setZero();
    }

protected:
    MaterialType _material_;
    
    SC_INLINE void GetGaussInformation(const Mesh &mesh, const SolutionFields &solutionfields, const ConstitutiveModel &material, integer ielement){
       // Get x and X
       // Compute auto kin  =  Kinematics(X,x,grad_bases);
       // ComputeConjugate(kin);
       // ComputeHessianComponents(kin);
       // FirstPiola(kin);
       // Hessian(kin);
    }

};


}
