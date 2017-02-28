#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include "mesh/Mesh.h"
#include "kinematics/Kinematics.h"
// #include "constitutive_models/ConstitutiveModels.h"
// #include "variational_formulations/VariationalFormulation.h"
// #include "variational_formulations/DisplacementFormulation.h"


namespace SoftComp {

class Assembly {
public:
    real tolerance;

    Assembly() {tolerance = 1e-14;};
    template<class MaterialType, template <class> class Formulation>
    std::tuple<spmatrix<real>,vector<real>>
    Assemble(const Mesh &mesh,
             Formulation<MaterialType> &formulation,
             MaterialType &material,
             SolutionFields<Formulation<MaterialType>> &fields) const {

        //---------------------------------------------------
        // Initialisation for stiffness matrix
        //---------------------------------------------------
        
        auto nodeperelem = mesh.nodeperelem;
        auto nvar = mesh.ndim;  // this needs to be changed!!
        auto local_capacity = nvar*nodeperelem*nvar*nodeperelem;
        auto capacity = local_capacity*mesh.nelem;

        using T = Eigen::Triplet<real>;
        std::vector<T> triplets;
        triplets.reserve(capacity);

        //---------------------------------------------------
        // Initialisation for global traction
        //---------------------------------------------------
        vector<real> gtraction = zeros(mesh.nnode*nvar);
        //---------------------------------------------------
        // Element loop
        //---------------------------------------------------
        //auto gcounter = 0;  //
        for (auto ielem=0; ielem<mesh.nelem; ++ielem) {
            //---------------------------------------------------
            // Get local Residual and Stiffness matrix
            //---------------------------------------------------
           formulation.LocalAssembly(mesh,material,fields,ielem);
            //---------------------------------------------------
            // Copy of local Residual and Stiffness matrix
            //---------------------------------------------------
           auto n_dofs               =  size(formulation.local_internal_traction);
           auto n_stiffness_entries  =  n_dofs*n_dofs;
           // matrix<real> local_stiffness(n_dofs,n_dofs);
           // std::copy(formulation.local_stiffness.data(), formulation.local_stiffness.data()+n_stiffness_entries,local_stiffness.data());
           // matrix<real> internal_traction(n_dofs,mesh.ndim);
           // std::copy(formulation.local_internal_traction.data(), formulation.local_internal_traction.data()+n_dofs,internal_traction.data());
            //---------------------------------------------------
            // Indices I and J and data including the non-zero
            // components of the stiffness
            //---------------------------------------------------
            for (auto i=0; i<n_stiffness_entries; ++i) {
               triplets.push_back(T(formulation.global_I_elem[i],
                                    formulation.global_J_elem[i],
                                    formulation.local_stiffness.data()[i]));
           }

           formulation.global_J_elem.clear();
           for (auto i=0; i<n_dofs; ++i) {
               gtraction[formulation.global_dofs[i]]  += formulation.local_internal_traction.data()[i];
           }
        }

        //---------------------------------------------------
        // Sparse matrix
        //---------------------------------------------------
        spmatrix<real> gstiffness(mesh.nnode*nvar,mesh.nnode*nvar);
        gstiffness.setFromTriplets(triplets.begin(),triplets.end());
        gstiffness.makeCompressed();

       return std::make_tuple(gstiffness,gtraction);
    }
};


}

#endif
