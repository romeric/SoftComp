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

};












}

#endif // MESH_H
