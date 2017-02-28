#ifndef FUNCTIONALSPACESVARIATIONALFORMULATIONS_H
#define FUNCTIONALSPACESVARIATIONALFORMULATIONS_H

#include "mesh/Mesh.h"
#include "quadrarure_rules/QuadratureRule.h"

namespace SoftComp{


SC_INLINE void GetQuadratureRules() {
    if (is_quadrature_computed) return;
    is_quadrature_computed = true;
    mass_quadrature.Gauss(2);
    volume_quadrature.Gauss(1);
    surface_quadrature.Gauss(1); // this has to be ndim-1
}


}

#endif // FUNCTIONALSPACESVARIATIONALFORMULATIONS_H
