#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H

#include "backends/matrix_backends/matrix_backend.h"
#include "mesh/Mesh.h"

namespace SoftComp {

class BoundaryCondition {
public:

    vector<real> dirichlet_values;
    vector<integer> dirichlet_fixed_dofs;
    vector<integer> dirichlet_free_dofs;

    vector<real> neumann_values;
    vector<integer> neumann_fixed_dofs;
    vector<integer> neumann_free_dofs;

    boolean is_solid, is_fluid, is_electro_mechanics, is_magneto_mechancics;


    BoundaryCondition() = default;

    template<typename T>
    void DirichletCriteria(vector<T> (*f)()) {
        /* Takes a user specified boundary condition
         * and modifies dirichlet_flags /fixed_dofs/free_dofs
         * accordingly
         *
         * input:
         *      f a function
        */

        dirichlet_values = f();
    }


    void DirichletCriteria(matrix<real> (*f)(Mesh), Mesh mesh) {
        /* Takes a user specified boundary condition
         * and modifies dirichlet_flags /fixed_dofs/free_dofs
         * accordingly
         *
         * input:
         *      f a function
        */

        dirichlet_values = ravel(f(mesh));
    }
};

}



#endif // BOUNDARYCONDITION_H
