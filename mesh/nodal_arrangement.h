#ifndef NODAL_ARRANGEMENT_H
#define NODAL_ARRANGEMENT_H

#include "backends/matrix_backends/matrix_backend.h"

namespace SoftComp {

class NodalArrangement {
public:
    matrix<integer> edge_arrangement;
    matrix<integer> face_arrangement;
    matrix<integer> element_arrangement;
    string element_type;
    integer degree;

    NodalArrangement(const string &etype, integer deg) {
        element_type = etype;
        degree = deg;
        auto C = degree - 1;
        const auto nsize = (C+2)*(C+2);

        if (element_type=="quad") {
            if (C>0) {
                vector<integer> edge0(1), edge1(1), edge2, edge3(1);
                edge0(0) = 4;
                edge1(0) = 2*C+5;
                edge3(0) = C + 4;

                for (auto i=1; i<C; ++i) {
                    edge0 = append(edge0,i+4);
                    edge1 = append(edge1, 2*C+5 + i*(C+2));
                    edge3 = append(edge3, C + 4 + i*(C+2));
                }
                std::cout << edge0 << "\n\n" << edge1 << std::endl;

                edge2 = arange<integer>(nsize-C,nsize);
                edge0 = vstack(arange<integer>(2),edge0);
                edge1 = vstack(arange<integer>(1,3),edge1);
                edge2 = vstack(arange<integer>(2,4),flipud(edge2));
    //            vector<integer> a3(2); a3(0)=3; a3(1)=0;
                vector<integer> a3(2); a3[0]=3; a3[1]=0;
                edge3 = vstack(a3,flipud(edge3));

                edge_arrangement = vstack(hstack(edge0,edge1),hstack(edge2,edge3));
            }
            else {
                edge_arrangement = zeros<integer>(2,4);
                edge_arrangement << 0, 1,
                                    1, 2,
                                    2, 3,
                                    3, 0;
            }


            // Node arrangement
            element_arrangement = zeros<integer>(nsize);
            element_arrangement(1) = C+1;
            element_arrangement(2) = nsize-1;
            element_arrangement(3) = nsize - (C+2);

            auto counter = 0;
            for (integer i=1; i<nsize; ++i) {
                if ( (i!=C+1) && (i!=nsize-1) && i!=(nsize-(C+2)) ) {
                    element_arrangement(counter+4) = i;
                    counter++;
                }
            }

        }
    }
};

}

#endif // NODAL_ARRANGEMENT_H
