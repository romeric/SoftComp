#include "SoftComp.h"
#include "examples/examples.h"

using namespace SoftComp;

//---------------------------------------------------------------------
//---------------------------------------------------------------------
//
//  Boundary conditions for this example
//
//---------------------------------------------------------------------
//---------------------------------------------------------------------
matrix<real> user_defined_dirichlet_data(const Mesh &mesh) {
    //-----------------------------------------
    // Dirichlet boundary conditions
    //-----------------------------------------
    matrix<real> dirichlet_values(mesh.nnode,2);
    for(integer i=0; i<mesh.nnode; i++) {
        if (abs(mesh.points(i,1))<1e-12) {
            dirichlet_values(i,0) = 0;
            dirichlet_values(i,1) = 0;
        }
        else {
            dirichlet_values(i,0) = INFINITY;
            dirichlet_values(i,1) = INFINITY;
        }
    }
    return dirichlet_values;
}
matrix<real> user_defined_neumann_data(const Mesh &mesh) {
    //-----------------------------------------
    // Neumann boundary conditions
    //-----------------------------------------
    matrix<real> neumann_values =  zeros(mesh.nnode,2);
    for(integer i=0; i<mesh.nnode; i++) {
        if (abs(mesh.points(i,0) - 1)<1e-12){
            if (abs(mesh.points(i,1) - 10)<1e-12) {
               neumann_values(i,0) = 0.002;
               neumann_values(i,1) = 0.0;
            }
        }
        else {
            neumann_values(i,0) = 0;
            neumann_values(i,1) = 0;
        }
    }
    return neumann_values;
}


int main() {
    //---------------------------------------------------------------------
    // Mesh for this example
    //---------------------------------------------------------------------
    integer Nx  = 40; integer Ny  =  400;  real Lx  =  1;  real  Ly  =  10;  integer degree  =  1;
    auto mesh                  =  Mesh("quad");
    mesh.Rectangle(Lx,Ly,Nx,Ny,degree);
    //---------------------------------------------------------------------
    // Boundary conditions
    //---------------------------------------------------------------------
    BoundaryCondition boundary_condition  =  BoundaryCondition();
    boundary_condition.DirichletCriteria(user_defined_dirichlet_data,mesh);
    boundary_condition.NeumannCriteria(user_defined_neumann_data,mesh);
    //---------------------------------------------------------------------
    // Constitutive model (Composite)
    //---------------------------------------------------------------------
    auto material = MooneyRivlin(1.,1.,5.);
    //---------------------------------------------------------------------
    // Variational formulation
    //---------------------------------------------------------------------
    auto formulation   =  DisplacementFormulation<decltype(material)>("quad",degree);
    auto fields        =  SolutionFields<DisplacementFormulation<decltype(material)>>(); 
    fields.InitialisedSolutionFields(mesh);
    auto assembly      =  Assembly();
    //---------------------------------------------------------------------
    // Newton-Raphson solver
    //---------------------------------------------------------------------
    auto solver = NetwonRaphsonSolver(1e-10,1);
    solver.Solve(mesh,formulation,boundary_condition,material,fields,assembly);


    return 0;
}
