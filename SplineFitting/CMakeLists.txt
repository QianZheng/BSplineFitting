cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(spline_fitting)

# Eigen
set(EIGEN_DIR  C:/Develop/eigen)
include_directories(${EIGEN_DIR}/include)

# ANN
set(ANN_DIR  C:/Develop/ann)
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