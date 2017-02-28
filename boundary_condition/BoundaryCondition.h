#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H

#include "backends/matrix_backends/matrix_backend.h"
//#include "mesh/Mesh.h"
//#include "math.h"

namespace SoftComp {

class BoundaryCondition {
public:

    vector<real> dirichlet_values;
    vector<integer> dirichlet_fixed_dofs;
    vector<integer> dirichlet_free_dofs;
    integer n_free_dofs;
    integer n_fixed_dofs;

    vector<real> neumann_values;
    vector<integer> neumann_fixed_dofs;
    vector<integer> neumann_free_dofs;

    boolean is_solid, is_fluid, is_electro_mechanics, is_magneto_mechancics;


    BoundaryCondition() = default;

    template<typename T>
    void DirichletCriteria(matrix<T> (*f)(const Mesh&), const Mesh &mesh) {
        //------------------------------------------------
        // Dirichlet boundary conditions
        //------------------------------------------------
        dirichlet_values      =  ravel(f(mesh));
        auto nfixed           =  0;
        auto nfree            =  0;
        auto ndofs            =  size(dirichlet_values,0);
        for (auto i=0;  i<ndofs; ++i){
            if (dirichlet_values(i)==INFINITY){
               nfree ++;
            }
            else{
                nfixed ++;
            }
        }
        dirichlet_fixed_dofs  =  zeros<integer>(nfixed);
        dirichlet_free_dofs   =  zeros<integer>(nfree);
        auto ifree            =  0;
        auto ifixed           =  0;
        for (auto i=0; i<ndofs; ++i){
            if (dirichlet_values(i)==INFINITY){
                dirichlet_free_dofs(ifree)   =  i;
                ifree++;
            }
            else{
                dirichlet_fixed_dofs(ifixed) =  i;
                ifixed++;
            }
        }
        n_free_dofs           =  nfree;
        n_fixed_dofs          =  nfixed;
    }

    template<typename T>
    void NeumannCriteria(matrix<T> (*f)(const Mesh&), const Mesh &mesh) {
        //------------------------------------------------
        // Neumann boundary conditions
        //------------------------------------------------
        neumann_values        =  ravel(f(mesh));
        auto nfixed           =  0;
        auto nfree            =  0;
        auto ndofs            =  size(neumann_values,0);
        for (auto i=0;  i<ndofs; ++i){
            if (neumann_values(i)==0){
               nfree ++;
            }
            else{
                nfixed ++;
            }
        }
        neumann_fixed_dofs    =  zeros<integer>(nfixed);
        neumann_free_dofs     =  zeros<integer>(nfree);
        auto ifree            =  0;
        auto ifixed           =  0;
        for (auto i=0; i<ndofs; ++i){
            if (neumann_values(i)==0){
                neumann_free_dofs(ifree)   =  i;
                ifree++;
            }
            else{
                neumann_fixed_dofs(ifixed) =  i;
                ifixed++;
            }
        }
    }


};

}



#endif // BOUNDARYCONDITION_H
