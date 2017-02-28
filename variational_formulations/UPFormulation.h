#ifndef UPFORMULATION_H
#define UPFORMULATION_H


#include "VariationalFormulation.h"
#include "backends/matrix_backends/backend_constitutivetensors.h"


namespace SoftComp {


template<class MaterialType>
class UPFormulation : public VariationalFormulation<MaterialType> {
public:
    tensor<real,2> local_internal_traction;
    tensor<real,4> local_stiffness;
    tensor<real,4> local_mass;
    integer nvar;
    QuadratureRule mass_quadrature, volume_quadrature, surface_quadrature;
    Kinematics kin;
    //boolean is_quadrature_computed;
    //boolean is_function_space_computed;
    std::vector<integer> global_dofs;
    FunctionSpace function_space, function_space_p;
    std::vector<integer> global_I_elem, global_J_elem;
    MaterialType material;

    UPFormulation(const string element_type, integer degree, integer degree_p) {
         //----------------------------------------------------
         // Get functional spaces
         //----------------------------------------------------
         function_space    =  FunctionSpace(element_type,degree);
         function_space_p  =  FunctionSpace(element_type,degree_p);
     }


    template<typename Formulation>
    SC_INLINE void GetGaussInfoInternalContribution(const Mesh &mesh, const MaterialType &mat,
        const SolutionFields<Formulation> & fields, integer ielement) {
         nvar = mesh.ndim;
        //----------------------------------------------------
        // Get x and X for a specific element (ielement)
        //----------------------------------------------------
        tensor<real,2>  Xelement(mesh.nodeperelem,mesh.ndim);
        tensor<real,2>  xelement(mesh.nodeperelem,mesh.ndim);
        tensor<real,1>  pelement();   // Get the pressure
        for (auto inode=0; inode<mesh.nodeperelem;  inode++){
            for (auto jdim=0; jdim<mesh.ndim; jdim++){
                Xelement(inode,jdim)      =  mesh.points(mesh.elements(ielement,inode),jdim);
                xelement(inode,jdim)      =  fields.x(mesh.elements(ielement,inode),jdim);
            }
        }
        //----------------------------------------------------
        // Get the kinematics for a specific element(ielement)
        //----------------------------------------------------
        kin  =  Kinematics(Xelement,xelement,function_space.grad_bases);
        //----------------------------------------------------
        // Get material properties (ielement)
        //----------------------------------------------------
        material =  mat;
        material.ComputeConjugate(kin);
        material.ComputeHessianComponents(kin);
        material.FirstPiola(kin);
        material.Hessian(kin);
        // Get geometrical term associated with the pressure!!!
    }

      void LocalResidualInternalContribution() {
      //----------------------------------------------------
      // Internal traction
      //----------------------------------------------------
#ifndef NDEBUG
        SC_ASSERT(!isempty(material.piola),"INVALID INITIALISATION. YOU NEED TO CALL GetGaussInfo FIRST");
        SC_ASSERT(!isempty(material.hessian),"INVALID INITIALISATION. YOU NEED TO CALL GetGaussInfo FIRST");
#endif

        auto ndim                    =  kin.F.dimension(1);
        auto ngauss                  =  kin.F.dimension(0);
        auto nodeperelem             =  kin.material_gradient.dimension(1);
        local_internal_traction      =  tensor<real,2>(nodeperelem,ndim);
        local_internal_traction.setZero();
        std::array<tpair, 1> idx0    = {tpair(1,1)};
        for (integer igauss=0; igauss<ngauss; ++igauss) {
            const tensor<real,2> &current_piola              =  material.piola.chip(igauss,0);
            const tensor<real,2> &current_material_gradient  =  kin.material_gradient.chip(igauss,0);
            local_internal_traction  +=  current_material_gradient.contract(current_piola,idx0)*
                                         function_space.quadrature_weights(igauss)*kin.material_iso_jacobian(igauss);
         }
    }

        void LocalStiffnessInternalContribution() {
        //----------------------------------------------------
        // Stiffness matrix for internal contribution
        //----------------------------------------------------
#ifndef NDEBUG
        SC_ASSERT(!isempty(material.piola),"INVALID INITIALISATION. YOU NEED TO CALL GetGaussInfo FIRST");
        SC_ASSERT(!isempty(material.hessian),"INVALID INITIALISATION. YOU NEED TO CALL GetGaussInfo FIRST");
#endif

        auto ndim         =  kin.F.dimension(1);
        auto ngauss       =  kin.F.dimension(0);
        auto nodeperelem  =  kin.material_gradient.dimension(1);

        local_stiffness   =  tensor<real,4>(nodeperelem,ndim,nodeperelem,ndim);
        tensor<real,4> temp_stiffness(nodeperelem,nodeperelem,ndim,ndim);
        temp_stiffness.setZero();
        std::array<tpair, 1> idx0 = {tpair(1,3)};
        std::array<tpair, 1> idx1 = {tpair(1,2)};
        for (integer igauss=0; igauss<ngauss; ++igauss) {
            const tensor<real,4> &current_hessian           = material.hessian.chip(igauss,0);
            const tensor<real,2> &current_material_gradient = kin.material_gradient.chip(igauss,0);
            auto temp        =  current_material_gradient.contract(current_hessian,idx0);
            temp_stiffness  +=  current_material_gradient.contract(temp,idx1)*
                                function_space.quadrature_weights(igauss)*kin.material_iso_jacobian(igauss);
        }
        //----------------------------------------------------
        // Permutation of 2nd and 3rd components of local stiffness
        // Very inefficient -- This needs to be permuated in-place instead of making a copy
        //----------------------------------------------------
        for (auto anode=0;anode<nodeperelem;anode++){
            for (auto idim=0;idim<ndim;idim++){
                for (auto bnode=0;bnode<nodeperelem;bnode++){
                    for (auto jdim=0;jdim<ndim;jdim++){
                        local_stiffness(anode,idim,bnode,jdim) = temp_stiffness(anode,bnode,idim,jdim);
                    }
                }
            }
        }
//        matrix<real> KK(8,8);
//        integer adof,bdof;
//        for (auto anode=0;anode<nodeperelem;anode++){
//            for (auto idim=0;idim<ndim;idim++){
//                for (auto bnode=0;bnode<nodeperelem;bnode++){
//                    for (auto jdim=0;jdim<ndim;jdim++){
//                        adof   =  idim + 2*anode;
//                        bdof   =  jdim + 2*bnode;
//                        KK(adof,bdof)  =  temp_stiffness(anode,bnode,idim,jdim);
//                    }
//                }
//            }
//        }
//        print(KK);
    }

    //void LocalMass() { // <-- Should use this instead
    void LocalMass(const Mesh &mesh, const Kinematics &kin, const ConstitutiveModel &mat) {
        //----------------------------------------------------
        // Mass matrix for an element
        //----------------------------------------------------
//        auto ndim         =  mesh.ndim;
//        auto ngauss       =  quadrule.mass_matrices.bases(0);  // Check this with Roman
//        auto nodeperelem  =  Mesh.nodeperelem;

//        local_mass         =  tensor<real,4>(nodeperelem,noderperelem,ndim,ndim); // Need to talk to Roman
//        local_mass.setZero();

//        tensor<real,2> bases(ngauss,nodeperelem);
//        bases.setRandom();

//        //tensor<real,ndim*2> outer(tensor<T,rank0> a, const tensor<T,rank1> &b)
//        for (auto anode=0;anode<nodeperelem;anode++){
//            for (auto bnode=0;bnode<nodeperelem;bnode++){
//                tensor<real,2> temp(ndim,ndim);
//                temp.SetZero();
//                for (auto igauss=0; igauss<ngauss; igauss++){
//                    //temp   =  quadrule.bases(igauss,anode)*quadrule.bases(igauss,bnode)*identity2d(ndim);
//                    //local_mass()  = local_mass() + material.density*temp*quadrule.mass.weights(igauss)*kin.material_iso_jacobian(igauss);
//                    temp   =  bases(igauss,anode)*bases(igauss,bnode)*identity2d(ndim);
//                    local_mass.chip(anode,0).chip(bnode,1)  = local_mass.chip(anode,0).chip(bnode,1) + material.density*temp*quadrule.mass.weights(igauss)*kin.material_iso_jacobian(igauss);
//                }
//            }
//        }


    }

//    void  LocalAssembly(const ){
    template<typename Formulation>
    void  LocalAssembly(const Mesh &mesh, const MaterialType &material, const SolutionFields<Formulation> & fields, integer ielem){
        //---------------------------------------------------
        // Get local Residual and Stiffness matrix
        //---------------------------------------------------
        GetGaussInfoInternalContribution(mesh,material,fields,ielem);
        LocalResidualInternalContribution();
        LocalStiffnessInternalContribution();
        //---------------------------------------------------
        // Initialise
        //---------------------------------------------------
        auto nodeperelem     =  kin.material_gradient.dimension(1);
        auto local_capacity  =  nvar*nodeperelem*nvar*nodeperelem;
        //---------------------------------------------------
        // Find the global degrees of freedom for the nodes of the current element
        //---------------------------------------------------
        global_dofs = std::vector<integer>(nvar*nodeperelem);
        for (auto counter=0; counter<nodeperelem; ++counter) {
            for (auto ncounter=0; ncounter<nvar; ++ncounter) {
                global_dofs[nvar*counter+ncounter] = nvar*mesh.elements(ielem,counter)+ncounter;
            }
        }
        //---------------------------------------------------
        // Indeces I and J
        //---------------------------------------------------
        global_I_elem  =  std::vector<integer>(local_capacity);
        global_J_elem  =  std::vector<integer>(local_capacity);
        // global_J_elem  =  std::vector<integer>;
        for (integer i=0; i<nvar*nodeperelem; ++i) {
            std::fill(global_I_elem.begin()+i*nvar*nodeperelem,global_I_elem.begin()+(i+1)*nvar*nodeperelem,global_dofs[i]);
            // global_J_elem.insert(global_J_elem.end(), global_dofs.begin(), global_dofs.end());
            for (integer j=0; j<nvar*nodeperelem; ++j) {
                global_J_elem[i*nvar*nodeperelem+j] = global_dofs[j];
            }
        }

//        //---------------------------------------------------
//        // Local counter accounting for the number of non-zero
//        // components of the stiffness
//        //---------------------------------------------------
//        auto lcounter = 0;
//        std::vector<integer> local_nnz; local_nnz.reserve(local_capacity);
//        for (auto i=0; i<nvar*nodeperelem; ++i) {
//            for (auto j=0; j<nvar*nodeperelem; ++j) {
//                if (abs(local_stiffness(i,j))>tolerance) {
//                    local_nnz.push_back(lcounter);
//                    lcounter++;
//                }
//            }
//        }
    }

//    void InitialisedSolutionFields(const Mesh &mesh){
//        //-------------------------------------------
//        //-------------------------------------------
//        x  =  tensor<real,2>(mesh.nnode,mesh.ndim);
//        X  =  tensor<real,2>(mesh.nnode,mesh.ndim);
//        for (auto inode=0;inode<mesh.nnode;inode++){
//            for (auto idim=0;idim<mesh.ndim;idim++){
//                x(inode,idim)  =  mesh.points(inode,idim);
//                X(inode,idim)  =  mesh.points(inode,idim);
//            }
//        }
//    }

// protected:
//     MaterialType _material_;

//    SC_INLINE void GetQuadratureRules() {
//        if (is_quadrature_computed) return;
//        is_quadrature_computed = true;
//        mass_quadrature.Gauss(2);
//        volume_quadrature.Gauss(1);
//        surface_quadrature.Gauss(1); // this has to be ndim-1
//    }

};

//--------------------------------------------------------------
// Specific unknown field of the displacement-based formulation:
// This include
//--------------------------------------------------------------
template<typename MaterialType>
class SolutionFields<DisplacementFormulation<MaterialType>> {
public:
    tensor<real,2> x;
    integer n_dofs_formulation;
    SolutionFields()           =  default;
    void InitialisedSolutionFields(const Mesh &mesh){
        //-------------------------------------------------------------------------------
       // Initialise the Solution Field x for the displacement formulation
        //-------------------------------------------------------------------------------
        x                      =  tensor<real,2>(mesh.nnode,mesh.ndim);
        for (auto inode=0;inode<mesh.nnode;inode++){
            for (auto idim=0;idim<mesh.ndim;idim++){
                x(inode,idim)  =  mesh.points(inode,idim);
            }
        }
        // Total number of dofs for this specific formulation
        n_dofs_formulation     =  size(x,0)*size(x,1);
    }
    void UpdateUnconstrainedSolutionFields(const Mesh &mesh,const BoundaryCondition &boundary_condition, const vector<real> &Dx){
        //-------------------------------------------------------------------------------
        // Update the Solution Field x for the displacement formulation (only for free dofs)
        //-------------------------------------------------------------------------------
        for (auto ifree=0;  ifree<boundary_condition.n_free_dofs;  ifree++){
            x.data()[boundary_condition.dirichlet_free_dofs(ifree)] += Dx(ifree);
        }
    }
    void UpdateConstrainedSolutionFields(const Mesh &mesh, const BoundaryCondition &boundary_condition, const real &accumulated_load_factor){
        //-------------------------------------------------------------------------------
        // Update the Solution Field x for degrees of freedom with imposed Dirichlet bc's
        //-------------------------------------------------------------------------------
        for (auto ifixed=0; ifixed<boundary_condition.n_fixed_dofs; ifixed++){
            x.data()[boundary_condition.dirichlet_fixed_dofs(ifixed)] = mesh.points.data()[ifixed] + accumulated_load_factor*boundary_condition.dirichlet_values(ifixed);
        }
    }
    tensor<real,1> operator-(const SolutionFields<DisplacementFormulation<MaterialType>> &oldfield){
        //-------------------------------------------------------------------------------
        // This function computes the difference bewteen the solution fields for two
        // different configurations
        //-------------------------------------------------------------------------------
        tensor<real,1> diff_fields(n_dofs_formulation);
        const real* nx_data = x.data();
        const real* ox_data = oldfield.x.data();
        for (auto idof=0; idof<n_dofs_formulation; idof++){
            diff_fields [idof]  =  nx_data[idof] - ox_data[idof];
        }
        return diff_fields;
    }

};


}
