#ifndef KINEMATICS_H
#define KINEMATICS_H

// #include "function_space/one_d/Line.h"
#include "backends/matrix_backends/matrix_backend.h"


#if defined(HAS_EIGEN)
#include <unsupported/Eigen/CXX11/Tensor>
template<typename T, integer rank>
using tensor = Eigen::Tensor<T, rank, Eigen::RowMajor>;
#endif

namespace SoftComp {




class Kinematics {
//    using namespace Fastor;
public:
    tensor<real,3> F, H;
    tensor<real,1> J, material_iso_jacobian;
    tensor<real,3> material_gradient;
    tensor<real,3> C, b;

    Kinematics() = delete;

    Kinematics(const matrix<real> &X, const matrix<real> &x, const tensor<real,3> &grad_bases) {


        auto ngauss = grad_bases.dimension(2);
        auto ndim = size(X,1);
        auto nodeperelem = size(X,0);

        matrix<real> parent_gradient = zeros(ndim,ndim);
        material_iso_jacobian = tensor<real,1>(ngauss);
        for (integer igauss=0; igauss<ngauss; ++igauss) {

            for (auto i=0; i<nodeperelem; ++i) {
                for (auto j=0; j<ndim; ++j) {
                    for (auto k=0; k<ndim; ++k) {
                        parent_gradient(j,k) += X(i,j)*grad_bases(i,k,igauss);
                    }
                }
            }
            material_iso_jacobian(igauss) = abs(det(parent_gradient));

            matrix<real> inv_parent_gradient_t = inv(parent_gradient);

//            for (auto i=0; i<nodeperelem; ++i) {
//                for (auto j=0; j<ndim; ++j) {
//                    for (auto k=0; k<ngauss; ++k) {
//                        material_gradient(j,k,igauss) += grad_bases(i,k,igauss)*inv_parent_gradient_t(i,j); // CHECK
//                    }
//                }
//            }

            //     F = transpose(x)*material_gradient;
            //     J = det(F);
            //     H = J*transpose(inv(F));

            material_gradient = tensor<real,3>(4,2,8);
            F = tensor<real,3>(2,2,8);
            H = tensor<real,3>(2,2,8);


        }
    }

};



//class Kinematics {
//public:
//    matrix<real> F, H, J;
//    matrix<real> material_gradient;
//    matrix<real> C, b;

//    Kinematics() = delete;


////    Kinematics(const matrix<real> &X, const matrix<real> &x, const matrix<real> &grad_bases) {
        
////        auto ngauss = size(grad_bases,1);
////        auto ndim = size(X,1);
////        auto nodeperelem = size(X,0);

////        material_gradient = zeros(nodeperelem*ndim,ngauss);
////        F = zeros(ndim*ndim,ngauss);
////        H = zeros(ndim*ndim,ngauss);
////        J = zeros(ngauss);


////        for (auto igauss=0; igauss<ngauss; ++igauss) {
////            auto current_grad_bases = reshape(grad_bases.col(igauss),nodeperelem,ndim);
////            matrix<real> parent_gradient = transpose(X)*current_grad_bases;
////            matrix<real> inv_parent_gradient_t = inv(parent_gradient);
////            matrix<real> current_material_gradient = current_grad_bases*(inv_parent_gradient_t);
////            material_gradient.col(igauss) = flatten(current_material_gradient);
////            F.col(igauss) = flatten(transpose(x)*current_material_gradient);
////            J(igauss) = det(F.col(igauss));
////            H.col(igauss) = J(igauss)*transpose(inv(reshape(F.col(igauss),ndim,ndim)));
////        }

////    }

//};

}

#endif // KINEMATICS_H
