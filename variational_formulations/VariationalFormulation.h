#ifndef VARIATIONALFORMULATION_H
#define VARIATIONALFORMULATION_H

#include "commons/commons.h"
#include "constitutive_models/ConstitutiveModels.h"
#include "quadrarure_rules/QuadratureRule.h"

namespace SoftComp {

class VariationalFormulation {
public:
    integer nvar;

    VariationalFormulation() = default;

    SC_INLINE void LocalStiffness() {}
    SC_INLINE void LocalMass() {}
    SC_INLINE void LocalResidual() {}

};


class MixedFormulation : public VariationalFormulation {

};

}

#endif // VARIATIONALFORMULATION_H
