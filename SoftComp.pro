TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp

HEADERS += \
    mesh/Mesh.h \
    commons/commons.h \
    commons/print.h \
    commons/timeit.h \
    commons/types.h \
    commons/utils.h \
    commons/write.h \
    SoftComp.h \
    backends/matrix_backends/backend_eigen.h \
    backends/matrix_backends/matrix_backend.h \
    backends/matrix_backends/backend_arma.h \
    quadrarure_rules/QuadratureRule.h \
    variational_formulations/VariationalFormulation.h \
    variational_formulations/DisplacementFormulation.h \
    readers/text_reader.h \
    readers/parser.h \
    function_space/one_d/Jacobi.h \
    mesh/nodal_arrangement.h \
    function_space/one_d/Line.h \
    backends/matrix_backends/backend_blaze.h \
    backends/matrix_backends/backend_builtin.h \
    examples/examples.h \
    boundary_condition/BoundaryCondition.h \
    constitutive_models/ConstitutiveModels.h \
    kinematics/Kinematics.h \
    backends/matrix_backends/backend_tensor.h \
    backends/matrix_backends/backend_constitutivetensors.h \
    assembly/Assembly.h \
    variational_formulations/UPFormulation.h \
    variational_formulations/DisplacementPhasefieldFormulation.h \
    kinematics/GradientsPhaseField.h \
    function_space/FunctionSpace.h \
    function_space/FunctionalSpacesVariationalFormulations.h

# BUILTIN
#QMAKE_CXXFLAGS += -std=c++11


# EIGEN
# QMAKE_CXXFLAGS += -std=c++11 -DNDEBUG -DHAS_EIGEN -Wno-unused-parameter #-O2
QMAKE_CXXFLAGS += -std=c++11 -DHAS_EIGEN -Wno-unused-parameter #-O2
#QMAKE_CXXFLAGS += -std=c++11 -O3 -mavx -DHAS_EIGEN #-fopenmp
#LIBS += -lgomp

# ARMA
#QMAKE_CXXFLAGS += -std=c++11 -O3 -mavx -DHAS_ARMA -DARMA_DONT_USE_WRAPPER
#LIBS += -llapack -lblas # XEON HAS DEFAULT FAST MULTI-THREAD BLAS


# BLAZE
#QMAKE_CXXFLAGS += -std=c++14 -DHAS_BLAZE -O2


#INCLUDEPATH += /home/roman/Dropbox/Fastor/

#INCLUDEPATH += /media/MATLAB/eigen-3.3.2/
#INCLUDEPATH += /media/MATLAB/blaze-3.0/
#INCLUDEPATH += /media/MATLAB/armadillo/include/



INCLUDEPATH += /home/rogelio/Documents/eigen3.3.2/
#INCLUDEPATH += /home/roman/eigen-3.3.2
#INCLUDEPATH += /home/roman/blaze/
#INCLUDEPATH += /home/roman/armadillo/include/

#DISTFILES += \
#    SoftComp \
#    SoftComp.pro.user.3deb69e \
#    SoftComp.pro.user.b02c5e2
