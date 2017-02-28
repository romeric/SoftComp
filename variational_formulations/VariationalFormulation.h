#ifndef VARIATIONALFORMULATION_H
#define VARIATIONALFORMULATION_H

#include "function_space/FunctionSpace.h"
#include "constitutive_models/ConstitutiveModels.h"
#include "quadrarure_rules/QuadratureRule.h"
#include "kinematics/Kinematics.h"
#include "boundary_condition/BoundaryCondition.h"

//#include "mesh/Mesh.h"

namespace SoftComp {

template<class Formulation>
class SolutionFields {};

template<class MaterialType>
class VariationalFormulation {
public:
    tensor<real,2> local_internal_traction;
    tensor<real,4> local_stiffness;
    tensor<real,4> local_mass;
    integer nvar;
    QuadratureRule quad_rule;
    vector<real> unknowns;


    VariationalFormulation() {
        quad_rule = QuadratureRule();
        // quad_rule.Gauss(2);
    }

    VariationalFormulation(const QuadratureRule &quadrule) {
        this->quad_rule = quadrule;
    }

    SC_INLINE void GetGaussInfo(const Mesh &mesh, const MaterialType &mat, integer ielement) {}
    SC_INLINE void LocalStiffness() {}
    SC_INLINE void LocalMass() {}
    SC_INLINE void LocalResidual() {}

protected:
    MaterialType _material_;
};



}

#endif // VARIATIONALFORMULATION_H
