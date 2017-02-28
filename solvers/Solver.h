#pragma once


#include "assembly/Assembly.h"
#include "commons/timeit.h"

namespace SoftComp {



class Solver {};


class NetwonRaphsonSolver: public Solver {    
public:
    real tolerance;
    integer number_of_increments;
    real load_factor;
    real accumulated_load_factor;
    integer max_iter;

    SC_INLINE NetwonRaphsonSolver(real tol, integer nincr=1) {
        tolerance = tol;
        number_of_increments = nincr;
        load_factor = 1./number_of_increments;
    }

    SC_INLINE void SetParameters(real tol, integer nincr=1) {
        tolerance = tol;
        number_of_increments = nincr;
        load_factor = 1./number_of_increments;
    }

    template<class MaterialType, template<class> class Formulation>
    void Solve(const Mesh &mesh, Formulation<MaterialType> &formulation,
            const BoundaryCondition &boundary_condition, MaterialType &material,
            SolutionFields<Formulation<MaterialType>> &fields, 
            const Assembly &assembly=Assembly()) {

            //---------------------------------------------
            // Sanity checks
            //---------------------------------------------
            print(FGRN(BOLD("Initialising the solver ...")));
            println(FGRN(BOLD("Number of elements:")), mesh.nelem, FGRN(BOLD(" Number of nodes")), mesh.nnode, '\n');

            //---------------------------------------------
            // Call the assembly
            //---------------------------------------------
            spmatrix<real> gstiffness; vector<real> gtraction;
            timer<real> t_assem, t_reduce, t_solve;
            t_assem.tic();
            std::tie(gstiffness,gtraction) = assembly.Assemble(mesh,formulation,material,fields);            
            t_assem.toc(FBLU(BOLD("Assembly time is: ")));
            //---------------------------------------------
            // Initialise the residual;
            //---------------------------------------------
            vector<real> residual    =  zeros(fields.n_dofs_formulation);
            accumulated_load_factor  =  0.0;
            //---------------------------------------------
            // Initialise the exterrnal force contribution
            //---------------------------------------------
            const vector<real> &total_nodal_forces   =   boundary_condition.neumann_values;
            //--------------------------------------------
            //--------------------------------------------
            // Loop for incremental load
            //--------------------------------------------
            //--------------------------------------------
            //while (accumulated_load_factor<1.0) {
            //--------------------------------------------
            // Update acumulated load factor
            //--------------------------------------------
            auto old_accumulated_load_factor =  accumulated_load_factor;
            accumulated_load_factor         +=  load_factor;
            auto delta_load_factor           =  accumulated_load_factor - old_accumulated_load_factor;
            auto oldfields                   =  fields;
            //--------------------------------------------
            // Update the value of the variables with associated
            // Dirichlet boundary conditions
            //--------------------------------------------
            fields.UpdateConstrainedSolutionFields(mesh,boundary_condition,accumulated_load_factor);
            tensor<real,1> diff_fields =  fields - oldfields;
            //--------------------------------------------
            // Add incremental load factor to the residual
            //--------------------------------------------
            // residual  += -delta_load_factor*total_nodal_forces;
            residual  -= delta_load_factor*total_nodal_forces;
            //--------------------------------------------
            // Get contribution of Dirichlet bc's in residual
            //--------------------------------------------
            t_reduce.tic();
            auto dirichlet_force     =  GetDirichletForceContribution(boundary_condition,fields.n_dofs_formulation,gstiffness,diff_fields);
            residual  +=  dirichlet_force;
            //--------------------------------------------
            // Reduced residual and stiffness matrix
            //--------------------------------------------
            auto reducedresidual      =  GetReducedVectors(residual,boundary_condition);
            auto reducedmatrix        =  GetReducedMatrix(gstiffness,boundary_condition);
            t_reduce.toc(FBLU(BOLD("Boundary condition time is: ")));
            //--------------------------------------------
            // Start Newton-Raphson algorithm for
            //--------------------------------------------
////******  while (norm(residual)<tolerance) {
            for (auto NRiteration=0; NRiteration<10; NRiteration++){
                //---------------------------------------------
                // Linear solver
                //---------------------------------------------
                t_solve.tic();
                vector<real> sol = spsolve(reducedmatrix,-reducedresidual);
                t_solve.toc(FBLU(BOLD("Solver time is: ")));
                //---------------------------------------------
                // Update unconstrained variables
                //---------------------------------------------
                fields.UpdateUnconstrainedSolutionFields(mesh,boundary_condition,sol);
                //---------------------------------------------
                // Call the assembly
                //---------------------------------------------
                std::tie(gstiffness,gtraction) = assembly.Assemble(mesh,formulation,material,fields);
                //---------------------------------------------
                // Update the residual
                //---------------------------------------------
                residual             =  gtraction - accumulated_load_factor*total_nodal_forces;
                //---------------------------------------------
                // Reduced residual and stiffness
                //---------------------------------------------
                reducedresidual      =  GetReducedVectors(residual,boundary_condition);
                reducedmatrix        =  GetReducedMatrix(gstiffness,boundary_condition);
                //---------------------------------------------
                // Check Newton-Raphson convergence
                //---------------------------------------------
                //---------------------------------------------
                // Print diagnostics for Newton-Raphson algorithm
                //---------------------------------------------
                println(FGRN(BOLD("Iteration")), NRiteration, FGRN(BOLD(" Norm of residual is: ")), 
                    norm(reducedresidual),"\n");
            }
        }


        //}
    vector<real> GetDirichletForceContribution(const BoundaryCondition &boundary_condition,
                                               integer n_dofs_formulation,const spmatrix<real> &gstiffness, const tensor<real,1> &diff_fields){
        //-----------------------------------------------------------
        // Contribution of Dirichlet boundary conditions on the right hand side
        //-----------------------------------------------------------
        // vector<real> dirichletforce = zeros(n_dofs_formulation);
        // vector<real> rowmatrix;
        // for (auto ifree=0; ifree<boundary_condition.n_free_dofs; ifree++){
        //     rowmatrix = static_cast<vector<real>>(gstiffness.row(boundary_condition.dirichlet_free_dofs(ifree)));
        //     for (auto idof=0; idof<n_dofs_formulation; idof++){
        //         dirichletforce(boundary_condition.dirichlet_free_dofs(ifree)) =+ rowmatrix(idof)*diff_fields(idof);
        //     }
        // }
        // return dirichletforce;

        // NEW APPROACH
        // -----------------------------------------------------------
        auto mapped_diff = Eigen::Map<const vector<real>>(diff_fields.data(),size(diff_fields));
        auto reducedmatrix = seqNC(gstiffness,boundary_condition.dirichlet_free_dofs.data(),boundary_condition.n_free_dofs);
        return reducedmatrix*mapped_diff;
        // -----------------------------------------------------------

    }
    SC_INLINE vector<real>  GetReducedVectors(const vector<real> &gtraction, const BoundaryCondition &boundary_condition){
        //-----------------------------------------------------------
        // Remove the fixed dofs from the traction vector
        //-----------------------------------------------------------
        vector<real> reducedtraction(boundary_condition.n_free_dofs);
        for (auto ifree=0; ifree<boundary_condition.n_free_dofs; ifree++){
            reducedtraction(ifree)  =  gtraction(boundary_condition.dirichlet_free_dofs(ifree));
        }
        return reducedtraction;
    }
    SC_INLINE spmatrix<real>  GetReducedMatrix(const spmatrix<real> &gstiffness, 
                const BoundaryCondition &boundary_condition){
        //-----------------------------------------------------------
        // Remove the fixed dofs from the stiffness
        //-----------------------------------------------------------

        // spmatrix<real> reducedmatrix(boundary_condition.n_free_dofs,boundary_condition.n_free_dofs);
        // auto b_ones    = block_sparse_extractor(boundary_condition.dirichlet_free_dofs.data(),size(boundary_condition.dirichlet_free_dofs),n_dofs_formulation);
        // reducedmatrix  =  b_ones*gstiffness*b_ones.transpose();
        // return reducedmatrix;

        // NEW APPROACH
        // -----------------------------------------------------------
        return seqNC(gstiffness,boundary_condition.dirichlet_free_dofs,boundary_condition.dirichlet_free_dofs);
    }
};

} // end of namespace SoftComp
