cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(spline_fitting)

# Eigen
set(EIGEN_DIR  c:/develop/Eigen)
#set(EIGEN_DIR  D:/develop/Eigen-3.2.7/eigen-lib)
include_directories(${EIGEN_DIR})

# ANN
set(ANN_DIR  c:/develop/ann)
#set(ANN_DIR  D:/develop/ann_1.1.2/ann-lib)
include_directories(${ANN_DIR}/include)
link_directories(${ANN_DIR}/lib)

#include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/core)



set(incs    
			core/cubic_b_spline.h
			core/spline_curve_fitting.h			
            )
			
set(srcs    core/cubic_b_spline.cpp
            core/spline_curve_fitting.cpp
			)
			

add_executable(spline_fitting_test ${incs} ${srcs} read_write_asc.h  read_write_asc.cpp main.cpp )


target_link_libraries(spline_fitting_test
			debug ANNd optimized ANN
			)