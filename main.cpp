#include "SoftComp.h"
#include "examples/examples.h"

using namespace SoftComp;


matrix<real> user_defined_dirichlet_data(Mesh mesh) {
    matrix<real> dirichlet_values(mesh.nnode,2);
    for(integer i=0; i<mesh.nnode; i++) {
        if (abs(mesh.points(i,1))<1e-12) {
            dirichlet_values(i,0) = 0;
            dirichlet_values(i,1) = 0;
        }
        else {
            dirichlet_values(i,0) = NAN;
            dirichlet_values(i,1) = NAN;
        }
    }
    return dirichlet_values;
}




//#include </home/roman/Dropbox/Fastor/Fastor.h>

int main() {

//    Fastor::Tensor<double,2,2,2> roger; roger.random();
//    print(roger);


//    return 0;


//    matrix<real> = readtxt();

// Mesh   //
    matrix<real> points = zeros(9,2);
    points << 0,0, 0.5,0, 1,0, 0,0.5, 0.5,0.5, 1.0,0.5, 0,1, 0.5,1, 1,1;
    matrix<integer> elements(4,4);
    elements << 0,1,4,3,1,2,5,4,3,4,7,6,4,5,8,7;


    auto mesh = Mesh("quad");
    mesh.elements = elements;
    mesh.points = points;
    mesh.nnode = size(mesh.points,0);
//    mesh.GetBoundaryEdges();
//    print(mesh.edges);

    // print(mesh.points,"\n", mesh.elements);

    // print("2222222222");
    auto quadbases = Quad(0,-0.5,0.5);
    print(quadbases.bases,quadbases.grad_bases);
    // auto quadbases = Line(7,-0.5);
    // Quad quadbases = Quad();
    // print(sum(quadbases.grad_bases));
    // auto dd = rand(100,100);
    // print(sum(dd,0),"\n",sum(dd,1)); 

    auto material = MooneyRivlin(2.,2.,5.);
    print(material.mu1,material.mu2,material.lambda);


//    quadbases.grad_bases = rand(8,4);
   tensor<real,3> grad_bases(4,2,8);
//   print(quadbases.grad_bases);
   matrix<real> X(4,2); X <<   0,0,
                               0.5,0,
                               0.5,0.5,
                               0,0.5;
   matrix<real> x(4,2); x <<   0,0,
                               0.5,0,
                               0.5,1,
                               0,1;

    auto kinematics = Kinematics(X,x,grad_bases);

    auto quad = QuadratureRule();
    quad.weights = rand(8);
    //quad1.Gauss(4);

//    auto quad= GaussLobatto(2);
//    print(quad0.points, quad0.weights);

    auto formulation = DisplacementFormulation();
    formulation.LocalResidual(quad,kinematics,material);
    print("WOW");
    return 0;

//    print(quadbases.bases, quadbases.grad_bases); 
   // quadbases.grad_bases = rand(4,2);


//    print(kinematics.F);
//    print(kinematics.H);
//    print(kinematics.J);


//    auto nodes =  seq(mesh.elements,0,{integer(0),size(mesh.elements,1)});
//    print(nodes);
//    vector<integer> nodes = mesh.elements.row(0);

//    vector<integer> dummy(2); dummy << 0,1;
//    auto X     =  seq(mesh.points,nodes,dummy);
//    print(nodes);



//// Boundary conditions //
////    BoundaryCondition boundary_condition = BoundaryCondition();
//    BoundaryCondition boundary_data = BoundaryCondition();
////    boundary_data.DirichletCriteria(vector<real>(*user_defined_dirichlet_data)(mesh));
////    boundary_data.DirichletCriteria(user_defined_dirichlet_data);
//    boundary_data.DirichletCriteria(user_defined_dirichlet_data,mesh);
//    print(boundary_data.dirichlet_values);
//    boundary_data.DirichletCriteria()

//    auto b = tensor(3,3,3);
//    b(0,1,2) = 1;

//    print(b);
//    integer a = 3;
//    matrix<integer> b = example_1(a);
//    print(b);


//    auto quadbases = Quad(2,-0.5,0.5);
//    print(quadbases.bases, quadbases.grad_bases);

    return 0;

    // Quadrature rule//
    QuadratureRule quad0 = QuadratureRule();  // Every class defines a type (Check equivalent line below!!)
    //auto quad0 = QuadratureRule();
    auto quad1 = QuadratureRule();
    //auto quad = QuadratureRule();
    quad0.GaussLobatto(4);
    //quad1.Gauss(4);

//    auto quad= GaussLobatto(2);
    print(quad0.points, quad0.weights);

    return 0;
}






//{
//    mesh
//    boundary_data
//    variation_formulation
//    function_space
//    linearsolver

//    solve(mesh,boundary);
//}































//    std::vector<int> vec  = {3, 2, 3, 3, 6, 5, 5, 6, 2, 6};
//    std::vector<int> uvec = {3, 2, 6, 5}; // vector of unique values
//    std::vector<int> ivec = {0, 1, 0, 0, 2, 3, 3, 2, 1, 2}; // vector of inverse
//    std::vector<int> orig_vec(ivec.size());
//    std::for_each(ivec.begin(), ivec.end(), [&uvec,&ivec,&orig_vec](int idx) {orig_vec[idx] = uvec[ivec[idx]];});
////    std::vector<int> orig_vec;
////    std::for_each(ivec.begin(), ivec.end(), [&uvec,&ivec,&orig_vec](int idx) {orig_vec.push_back(uvec[ivec[idx]]);});
////    std::for_each(ivec.begin(), ivec.end(), [&uvec,&orig_vec](int idx) {print(idx);});
//    print(orig_vec);

//    vector<int> vecc(10);
//#ifdef HAS_EIGEN
//    vecc << 3, 2, 3, 3, 6, 5, 5, 6, 2, 6;
//#else
//    vecc = {3, 2, 3, 3, 6, 5, 5, 6, 2, 6};
//    print(type_name<decltype(size(vecc))>());
//#endif

//    return 0;
//}
