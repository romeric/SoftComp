#ifndef CONSTITUTIVEMODELS_H
#define CONSTITUTIVEMODELS_H

#include "backends/matrix_backends/matrix_backend.h"
#include "backends/matrix_backends/backend_tensor.h"
#include "backends/matrix_backends/backend_constitutivetensors.h"
#include "kinematics/Kinematics.h"


namespace SoftComp {


class ConstitutiveModel {
    /*-------------------------------------------------------------------
     * ------------------------------------------------------------------
     * ------------------------------------------------------------------
    Base class for computing the first Piola-Kirchhoff stress tensor
    and the equivalent Hessian coming from the derived constitutive model
    * -------------------------------------------------------------------
    * -------------------------------------------------------------------
    * -------------------------------------------------------------------
    */
public:
    tensor<real,3> sigma_F, sigma_H;
    tensor<real,1> sigma_J;

    tensor<real,5> W_FF, W_FH, W_HF, W_HH;
    tensor<real,3> W_FJ, W_HJ, W_JF, W_JH;
    tensor<real,1> W_JJ;
    tensor<real,5> geometrical_component;

    tensor<real,3> piola;
    tensor<real,5> hessian;


    ConstitutiveModel() {}

    //------------------------------------------------------------
    // Compute first Piola
    //------------------------------------------------------------
    SC_INLINE void FirstPiola(const Kinematics &kin) {
    #ifdef NDEBUG
    #endif
        piola      =  sigma_F + piolaH(sigma_H,kin.F) + piolaJ(sigma_J,kin.H);
    }
    //------------------------------------------------------------
    // Compute final Hessian (including all the components of the Hessian)
    //------------------------------------------------------------
    SC_INLINE void Hessian(const Kinematics &kin) {
        auto ndim  =  kin.F.dimension(1);
        hessian    =  W_FF + H_HH(kin.F,W_HH,ndim) + H_JJ(kin.H,W_JJ,ndim) + H_geom(kin.F,sigma_H,sigma_J,ndim);// + H_FH(kin.F,W_FH,ndim) + H_HH(kin.F,W_HH,ndim)+ H_FJ(kin.H,W_FJ,ndim) + H_HF(kin.F,W_HF,ndim) + H_HH(kin.F, W_HH,ndim);
    }

    //void ComputeConjugate(const Kinematics &kin) {}
    //void ComputeHessianComponents(const Kinematics &kin) {}
};

//--------------------------------------------------------------
//--------------------------------------------------------------
// Mooney-Rivlin constitutive tensor
//--------------------------------------------------------------
//--------------------------------------------------------------
class MooneyRivlin: public ConstitutiveModel {
public:
    real mu1, mu2, lambda;

    MooneyRivlin() = default;
    MooneyRivlin(real u1, real u2, real u3) : mu1(u1), mu2(u2), lambda(u3) {}

    MooneyRivlin(const MooneyRivlin &other) {
        mu1    = other.mu1;
        mu2    = other.mu2;
        lambda = other.lambda;
    }

    SC_INLINE void ComputeConjugate(const Kinematics &kin) {
        //------------------------------------------------------------
        // Compute all the work conjugates of the model (3D and 2D)
        //------------------------------------------------------------
        auto ndim    =  size(kin.F,1);
        if (ndim==3){
           sigma_F   =  mu1*kin.F;
           sigma_H   =  mu2*kin.H;
           sigma_J   =  (lambda*(kin.J - 1.)- (mu1+2.*mu2)/kin.J);
        }
        else{
            sigma_F  =  mu1*kin.F;
            sigma_H  =  mu2*kin.H;
            sigma_J  =  (lambda*(kin.J - 1.) - (mu1+2.*mu2)/kin.J) + mu2*kin.J;
        }
    }
    SC_INLINE void ComputeHessianComponents(const Kinematics &kin) {
        //------------------------------------------------------------
        // Compute all the components of the Hessian
        //------------------------------------------------------------
        auto ngauss  =  kin.F.dimension(0);
        auto ndim    =  kin.F.dimension(1);
        auto id4     =  identity4d(ndim);
        W_FF         =  tensor<real,5>(ngauss,ndim,ndim,ndim,ndim);
        W_HH         =  tensor<real,5>(ngauss,ndim,ndim,ndim,ndim);
        W_JJ         =  tensor<real,1>(ngauss);
        W_FH         =  tensor<real,5>(ngauss,ndim,ndim,ndim,ndim);
        W_HF         =  tensor<real,5>(ngauss,ndim,ndim,ndim,ndim);
        W_FJ         =  tensor<real,3>(ngauss,ndim,ndim);
        W_JF         =  tensor<real,3>(ngauss,ndim,ndim);
        W_HJ         =  tensor<real,3>(ngauss,ndim,ndim);
        W_JH         =  tensor<real,3>(ngauss,ndim,ndim);
        //***********************************
        // Non-zero components of the Hessian
        //***********************************
        if (ndim==3){
           for (integer igauss=0; igauss<ngauss; ++igauss) {
               W_FF.chip(igauss,0)  =  mu1*id4;
               W_HH.chip(igauss,0)  =  mu2*id4;
               W_JJ(igauss)         =  (mu1+2.*mu2)/(kin.J(igauss)*kin.J(igauss)) + lambda;
           }
        }
        else {
            for (integer igauss=0; igauss<ngauss; ++igauss) {
                W_FF.chip(igauss,0)  =  mu1*id4;
                W_HH.chip(igauss,0)  =  mu2*id4;
                W_JJ(igauss)         =  (mu1+2.*mu2)/(kin.J(igauss)*kin.J(igauss)) + mu2 + lambda;
            }
        }
    }

};

}

#endif // CONSTITUTIVEMODELS_H
