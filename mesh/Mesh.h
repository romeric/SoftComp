#ifndef MESH_H
#define MESH_H


#include "commons/commons.h"
#include "nodal_arrangement.h"
#include "backends/matrix_backends/matrix_backend.h"

#include "commons/print.h"


namespace SoftComp {

enum mesh_reader {salome,gmsh,builtin};


class Mesh {
public:
    string element_type; // USE ENUM?
    matrix<real> points;
    matrix<integer> elements;
    matrix<integer> edges;
    matrix<integer> faces;
    matrix<integer> boundary_edges;
    matrix<integer> boundary_faces;
    matrix<integer> interior_edges;
    matrix<integer> interior_faces;
    integer nelem;
    integer nnode;
    integer degree;
    integer ndim;
    integer nodeperelem;
    integer edgeperelem;
    integer faceperelem;

    Mesh() = default;
    SC_INLINE Mesh(string etype) noexcept {
        this->element_type = etype;
        if (element_type=="tri" || element_type=="quad") {
            this->ndim = 2;
        }
        else if (element_type=="tet" || element_type=="pyr" || element_type=="pri" || element_type=="hex") {
            this->ndim = 3;
        }

        if (element_type=="tri") {this->edgeperelem = 3;}
        else if (element_type=="quad") {this->edgeperelem = 4;}
        else if (element_type=="tet") {this->edgeperelem = 6; this->faceperelem=4;}
        else if (element_type=="hex") {this->edgeperelem = 12; this->faceperelem=8;}
    }


    template<int rtype=mesh_reader::salome>
    SC_INLINE void Read(const string &fname) {
        std::cout << fname;
    }

    SC_INLINE void GetBoundaryEdges() {

        SC_ASSERT(!element_type.empty(),"ELEMENT TYPE NOT UNDERSTOOD");
        SC_ASSERT(!isempty(elements),"NODAL CONNECTIVITY IS EMPTY");
        SC_ASSERT(edgeperelem!=0,"NODAL CONNECTIVITY IS EMPTY");

        InferPolynomialDegree();
        nelem = rows(elements);
        auto arranger = NodalArrangement(element_type,degree-1);
        const auto &node_arranger = arranger.edge_arrangement;


//        print(arranger.edge_arrangement);
//        if (element_type=="quad") {
            edges = zeros<integer>(edgeperelem*nelem,2);
            for (integer i=0; i<nelem; ++i) {
                for (integer j=0; j<edgeperelem; ++j) {
                    edges(i+j*nelem,0) = elements(i,node_arranger(j,0));
                    edges(i+j*nelem,1) = elements(i,node_arranger(j,1));
                }

//                edges(i+nelem,0) = elements(i,node_arranger(1,0));
//                edges(i+nelem,1) = elements(i,node_arranger(1,1));

//                edges(i+2*nelem,0) = elements(i,node_arranger(2,0));
//                edges(i+2*nelem,1) = elements(i,node_arranger(2,1));

//                edges(i+3*nelem,0) = elements(i,node_arranger(3,0));
//                edges(i+3*nelem,1) = elements(i,node_arranger(3,1));
            }
//                    vstack((self.elements[:,node_arranger[0,:]],self.elements[:,node_arranger[1,:]],
//                                         self.elements[:,node_arranger[2,:]]));
//            print(edges);
//        }
    }

    SC_INLINE void GetBoundaryFaces() {
    }

    SC_INLINE integer InferPolynomialDegree() {

        SC_ASSERT(!element_type.empty(),"ELEMENT TYPE NOT UNDERSTOOD");

        integer p = 0;
        if (element_type == "tri") {
            for (auto i=0; i<105; ++i) {
                if ((i+1)*(i+2)/2==cols(elements)) {
                    p = i;
                    break;
                }
            }
        }
        else if (element_type == "tet") {
            for (auto i=0; i<105; ++i) {
                if ((i+1)*(i+2)*(i+3)/6==cols(elements)) {
                    p = i;
                    break;
                }
            }
        }
        else if (element_type == "quad") {
            for (auto i=0; i<105; ++i) {
                if ((i+1)*(i+1)==cols(elements)) {
                    p = i;
                    break;
                }
            }
        }
        else if (element_type == "hex") {
            for (auto i=0; i<105; ++i) {
                if ((i+1)*(i+1)*(i+1)==cols(elements)) {
                    p = i;
                    break;
                }
            }
        }

        degree = p;
        return p;
    }

    void Rectangle(real Lx, real Ly, integer Nx, integer Ny, integer deg, string etype="quad") {
        //--------------------------------------------
        // Generate the points of the mesh
        //--------------------------------------------

        if (etype!="quad") {
            println("Generating mesh for ", etype, " not yet implemented");
            return;
        }

        real Dx                        =  Lx/(Nx*deg);
        real Dy                        =  Ly/(Ny*deg);
        vector<real> XX                =  zeros(Nx*deg+1);
        vector<real> YY                =  zeros(Ny*deg+1);
        for (auto ix=1; ix<Nx*deg+1; ix++){
            XX(ix)                     =  XX(ix-1) + Dx;
        }
        for (auto iy=1; iy<Ny*deg+1; iy++){
            YY(iy)                     =  YY(iy-1) + Dy;
        }
        integer npoints                =  (Nx*deg + 1)*(Ny*deg + 1);
        integer ndim                   =  2;
        integer dummy                  =  0;
        points                         =  zeros(npoints,ndim);
        for (auto iy=0; iy<Ny*deg+1;  iy++){
            for (auto ix=0; ix<Nx*deg+1; ix++){
                points(dummy,0)        =  XX(ix);
                points(dummy,1)        =  YY(iy);
                dummy++;
            }
        }

        //--------------------------------------------
        // Generate the connectivity of the mesh
        //--------------------------------------------
        elements = matrix<integer>(Nx*Ny,(deg+1)*(deg+1));
        if (deg==0){
        }
        if (deg==1){
           integer ielem              =  0;
           for (auto iy=0; iy<Ny; iy++){
               for (auto ix=0; ix<Nx; ix++){
                   elements(ielem,0)  =  ix + 0 + iy*(Nx + 1);
                   elements(ielem,1)  =  ix + 1 + iy*(Nx + 1);
                   elements(ielem,3)  =  ix + 0 + (iy + 1)*(Nx + 1);
                   elements(ielem,2)  =  ix + 1 + (iy + 1)*(Nx + 1);
                   ielem++;
               }
           }
        }
        if (deg==2){
            integer ielem              =  0;
            for (auto iy=0; iy<Ny; iy++){
                for (auto ix=0; ix<Nx; ix++){
                    elements(ielem,0)  =  ix + 0 + (Nx*deg + 1)*0 + ix*(deg - 1) + iy*(Nx*deg+1)*deg;
                    elements(ielem,4)  =  ix + 1 + (Nx*deg + 1)*0 + ix*(deg - 1) + iy*(Nx*deg+1)*deg;
                    elements(ielem,1)  =  ix + 2 + (Nx*deg + 1)*0 + ix*(deg - 1) + iy*(Nx*deg+1)*deg;
                    elements(ielem,7)  =  ix + 0 + (Nx*deg + 1)*1 + ix*(deg - 1) + iy*(Nx*deg+1)*deg;
                    elements(ielem,8)  =  ix + 1 + (Nx*deg + 1)*1 + ix*(deg - 1) + iy*(Nx*deg+1)*deg;
                    elements(ielem,5)  =  ix + 2 + (Nx*deg + 1)*1 + ix*(deg - 1) + iy*(Nx*deg+1)*deg;
                    elements(ielem,3)  =  ix + 0 + (Nx*deg + 1)*2 + ix*(deg - 1) + iy*(Nx*deg+1)*deg;
                    elements(ielem,6)  =  ix + 1 + (Nx*deg + 1)*2 + ix*(deg - 1) + iy*(Nx*deg+1)*deg;
                    elements(ielem,2)  =  ix + 2 + (Nx*deg + 1)*2 + ix*(deg - 1) + iy*(Nx*deg+1)*deg;
                    ielem++;
                }
            }
        }

        element_type          =  etype;
        degree                =  deg;
        nnode                 =  size(points,0);
        nodeperelem           =  size(elements,1);
        nelem                 =  size(elements,0);

}

};












}

#endif // MESH_H
