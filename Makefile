# EIGEN_INCLUDE=/media/MATLAB/eigen-3.3.2/
# ARMA_INCLUDE=/media/MATLAB/armadillo/include/
# BLAZE_INCLUDE=/media/MATLAB/blaze-3.0/


EIGEN_INCLUDE=/home/roman/eigen-3.3.2/ 
ARMA_INCLUDE=/home/roman/armadillo/include/
BLAZE_INCLUDE=/home/roman/blaze/

# CXX_FLAGS = -std=c++11 -O3 -mavx -Wno-ignored-attributes -DNDEGUG
CXX_FLAGS = -std=c++11 -Wno-ignored-attributes -DNDEBUG

all:
	$(CXX) $(CXX_FLAGS) -I$(BLAZE_INCLUDE) -I$(EIGEN_INCLUDE) -I$(ARMA_INCLUDE) -I. main.cpp -o SoftComp -DHAS_EIGEN

eigen_par:
	$(CXX) $(CXX_FLAGS) -I$(BLAZE_INCLUDE) -I$(EIGEN_INCLUDE) -I$(ARMA_INCLUDE) -I. -DHAS_EIGEN main.cpp -o SoftComp -fopenmp


arma:
	$(CXX) $(CXX_FLAGS) -I$(BLAZE_INCLUDE) -I$(EIGEN_INCLUDE) -I$(ARMA_INCLUDE) -I. -DHAS_ARMA main.cpp -o SoftComp -lblas -llapack

blaze:
	$(CXX) $(CXX_FLAGS) -I$(BLAZE_INCLUDE) -I$(EIGEN_INCLUDE) -I$(ARMA_INCLUDE) -I. main.cpp -o SoftComp -DHAS_BLAZE -std=c++14