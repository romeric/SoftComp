#ifndef CONSTITUTIVEMODELS_H
#define CONSTITUTIVEMODELS_H

#include "backends/matrix_backends/matrix_backend.h"
#include "kinematics/Kinematics.h"


namespace SoftComp {



tensor<real,2> __cross__(const tensor<real,2> &A, const tensor<real,2> &B) {

    real A00=A(0,0);
    real A11=A(1,1);
    real A22=A(2,2);
    real A01=A(0,1);
    real A02=A(0,2);
    real A12=A(1,2);
    real A10=A(1,0);
    real A20=A(2,0);
    real A21=A(2,1);

    real B00=B(0,0);
    real B11=B(1,1);
    real B22=B(2,2);
    real B01=B(0,1);
    real B02=B(0,2);
    real B12=B(1,2);
    real B10=B(1,0);
    real B20=B(2,0);
    real B21=B(2,1);

    tensor<real,2> out;
    if (A.dimension(0)==3) {
        out = tensor<real,2>(3,3);
        out.setValues({{A11*B22 - A12*B21 - A21*B12 + A22*B11, A12*B20 - A10*B22 + A20*B12 - A22*B10, A10*B21 - A11*B20 - A20*B11 + A21*B10},
                   {A02*B21 - A01*B22 + A21*B02 - A22*B01, A00*B22 - A02*B20 - A20*B02 + A22*B00, A01*B20 - A00*B21 + A20*B01 - A21*B00},
                   {A01*B12 - A02*B11 - A11*B02 + A12*B01, A02*B10 - A00*B12 + A10*B02 - A12*B00, A00*B11 - A01*B10 - A10*B01 + A11*B00}});
    }
    else {
        out = tensor<real,2>(2,2);
        B22 = 1., A22 = 1., B02 = 0., B12=0., B20 = 0., B21 = 0., A02 = 0., A12=0., A20 = 0., A21 = 0.;
        out.setValues({{A11*B22 - A12*B21 - A21*B12 + A22*B11, A12*B20 - A10*B22 + A20*B12 - A22*B10},
                       {A02*B21 - A01*B22 + A21*B02 - A22*B01, A00*B22 - A02*B20 - A20*B02 + A22*B00}});
    }

    return out;
}

//template<typename T, integer rank>
//tensor<T,rank> cross(const tensor<T,rank> &a, const tensor<T,rank> &b) {
//    return tensor<T,rank>{};
//}


tensor<real,3> cross(const tensor<real,3> &a, const tensor<real,3> &b) {
    auto ngauss = b.dimension(2);
    tensor<real,3> out(b.dimension(0),b.dimension(1),b.dimension(2));
    for (auto igauss=0; igauss<ngauss; ++igauss) {
        out.chip(igauss,2) = __cross__(a.chip(igauss,2), b.chip(igauss,2));
    }
    return out;
}



static tensor<real,3> levi_civita(integer ndim) {
    tensor<real, 3> epsilon(ndim,ndim,ndim);
    epsilon.setZero();
    epsilon(0,1,2) = 1;
    epsilon(1,2,0) = 1;
    epsilon(2,0,1) = 1;
    epsilon(1,0,2) = -1;
    epsilon(2,1,0) = -1;
    epsilon(0,2,1) = -1;
    return epsilon;
}

static tensor<real,4> identity(integer ndim) {
    matrix<real> identity_2d = eye(ndim,ndim);
    tensor<real, 4> identity_4d(ndim,ndim,ndim,ndim);
    for (integer i=0; i<ndim; ++i) {
        for (integer I=0; I<ndim; ++I) {
            for (integer j=0; j<ndim; ++j) {
                for (integer J=0; J<ndim; ++J) {
                    identity_4d(i,I,j,J) = identity_2d(i,j)*identity_2d(I,J);
                }
            }
        }
    }
    return identity_4d;
}


//#######################

class ConstitutiveModel {
public:
    tensor<real,3> sigma_F, sigma_H;
    tensor<real,1> sigma_J;
    tensor<real,3> piola;

    tensor<real,5> W_FF, W_FH, W_HF, W_HH;
    tensor<real,3> W_FJ, W_HJ, W_JF, W_JH;
    tensor<real,1> W_JJ;
    tensor<real,5> hessian;

    tensor<real,5> geometrical_component;

    ConstitutiveModel() {}

    SC_INLINE void FirstPiola(const Kinematics &kin) {

//        std::array<Eigen::IndexPair<integer>, 1> idx = { Eigen::IndexPair<integer>(0, 2) };
//        tensor<real,3> temp = sigma_J.contract(kin.H, idx);

        auto ngauss = kin.F.dimension(2);
        auto ndim = kin.F.dimension(0);
        tensor<real,3> temp(ndim,ndim,ngauss);
        for (integer igauss=0; igauss<ngauss; ++igauss) {
            temp.chip(igauss,2)  = sigma_J(igauss)*kin.H.chip(igauss,2);
        }
        piola = sigma_F + cross(sigma_H,kin.F) + temp;
    }

    void Hessian(const Kinematics &kin) {

        auto ndim = kin.F.dimension(0);
        auto ngauss = kin.F.dimension(2);
        auto epsilon = levi_civita(ndim);
        tensor<real,5> hessian_HH;

//#pragma simd
        for (integer i=0; i<ndim; ++i) {
            for (integer p=0; p<ndim; ++p) {
                for (integer q=0; q<ndim; ++q) {
                    for (integer I=0; I<ndim; ++I) {
                        for (integer P=0; P<ndim; ++P) {
                            for (integer Q=0; Q<ndim; ++Q) {
                                for (integer j=0; j<ndim; ++j) {
                                    for (integer r=0; r<ndim; ++r) {
                                        for (integer s=0; s<ndim; ++s) {
                                            for (integer J=0; J<ndim; ++J) {
                                                for (integer R=0; R<ndim; ++R) {
                                                    for (integer S=0; S<ndim; ++S) {
                                                        for (integer igauss; igauss<ngauss; ++igauss) {
                                                            hessian_HH(i,I,j,J,igauss) += epsilon(i,p,q)*epsilon(I,P,Q)*\
                                                                    epsilon(j,r,s)*epsilon(J,R,S)*kin.F(p,P,igauss)*W_HH(q,Q,r,R,igauss)*kin.F(s,S,igauss);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        const auto &hessian_FF = W_FF;
        hessian = hessian_FF + hessian_HH + geometrical_component;
    }
};









class MooneyRivlin: public ConstitutiveModel {
public:
    real mu1, mu2, lambda;
//    tensor<real,3> piola;
//    tensor<real,3> sigma_F, sigma_H;
//    tensor<real,1> sigma_J;

//    tensor<real,5> W_FF, W_HH;
//    tensor<real,1> W_JJ;

    MooneyRivlin(real u1, real u2, real u3) : mu1(u1), mu2(u2), lambda(u3) {}

    SC_INLINE void ComputeConjugate(const Kinematics &kin) {

        sigma_F = mu1*kin.F;
        sigma_H = mu2*kin.H;
        sigma_J = (lambda*(kin.J- 1.)- 2.*(mu1+2.*mu2)/kin.J);
    }

    SC_INLINE void ComputeHessianComponents(const Kinematics &kin) {

        auto ngauss = kin.F.dimension(2);
        auto ndim = kin.F.dimension(0);
        auto id4 = identity(ndim);
        W_FF = tensor<real,5>(ndim,ndim,ndim,ndim,ngauss);
        W_HH = tensor<real,5>(ndim,ndim,ndim,ndim,ngauss);
        W_JJ = tensor<real,1>(ngauss);
        for (integer igauss=0; igauss<ngauss; ++igauss) {
            W_FF.chip(igauss,4) = mu1*id4;
            W_HH.chip(igauss,4) = mu2*id4;
            W_JJ(igauss) = 2.*(mu1+2.*mu2)/(kin.J(igauss)*kin.J(igauss)) + lambda;
        }


        W_FH = tensor<real,5>(ndim,ndim,ndim,ndim,ngauss);
        W_HF = tensor<real,5>(ndim,ndim,ndim,ndim,ngauss);
        W_FJ = tensor<real,3>(ndim,ndim,ngauss);
        W_JF = tensor<real,3>(ndim,ndim,ngauss);
        W_HJ = tensor<real,3>(ndim,ndim,ngauss);
        W_JH = tensor<real,3>(ndim,ndim,ngauss);

        // Compute Geometric term
        tensor<real,3> temp_b(ndim,ndim,ngauss);
        for (integer igauss=0; igauss<ngauss; ++igauss) {
            temp_b.chip(igauss,2)  = kin.H.chip(igauss,2) + sigma_J(igauss)*kin.F.chip(igauss,2);
        }
//        auto B = sigma_H;
        auto epsilon = levi_civita(ndim);
        std::array<Eigen::IndexPair<integer>, 1> idx = { Eigen::IndexPair<integer>(0, 0) };
        tensor<real,4> temp0 = epsilon.contract(temp_b, idx);

        idx = { Eigen::IndexPair<integer>(0, 2) };
        tensor<real,5> temp1 = epsilon.contract(temp0, idx);


    }



};

}

#endif // CONSTITUTIVEMODELS_H
