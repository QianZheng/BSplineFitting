BSplineFitting
==============

Fitting cubic spline curve to 2d points

## Introduction ##

This is an implementation of paper "Fitting B-spline Curves to Point Clouds
by Curvature-Based Squared Distance Minimization".

Link to the paper: [http://www.geometrie.tuwien.ac.at/ig/sn/2006/wpl_curves_06/wpl_curves_06.html](http://www.geometrie.tuwien.ac.at/ig/sn/2006/wpl_curves_06/wpl_curves_06.html "link")

The input is a set of 2d points, the output are control points of a close cubic spline curve.

1. Input file: a file that contains n rows and each row reprents a point with x y positions.
2. Output files: a file with the control points of the output curve; and a file with sampling points of the output file. 

## Third-Party Dependencies ##

This project depends on two code librarie:

1. Eigen 3: all matrix and vector operations are based on this library
2. ANN: it is used to compute the nearest neighbor of a given point.

## Building and Running ##
BSplineFitting should be able to run in any environment, but it is **only** tested in windows enviroment.

I use cmake to configure and generate project files.

## Main Files##
drawResult.m: a simple .m file to visulize the input and output

core/cubic_b_spline.h: a class encode the cubic b spline

core/spline_curve_fitting.h: 

read_write_asc.h: a simple class that reads/writes files









