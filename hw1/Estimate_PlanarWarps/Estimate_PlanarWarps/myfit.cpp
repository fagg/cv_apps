//
//  myfit.cpp
//  Estimate_Homography
//
//  Created by Simon Lucey on 9/21/15.
//  Copyright (c) 2015 CMU_16432. All rights reserved.
//

#include "myfit.h"

// Use the Armadillo namespace
using namespace arma;
using namespace std;

//-----------------------------------------------------------------
// Function to return the affine warp between 3D points on a plane
//
// <in>
// X = concatenated matrix of 2D projected points in the image (2xN)
// W = concatenated matrix of 3D points on the plane (3XN)
//
// <out>
// A = 2x3 matrix of affine parameters
fmat myfit_affine(fmat &X, fmat &W) {

    // Fill in the answer here.....
    
    // create x', which is a 2n x 1 matrix of "3D" values
    fmat x_prime;
    
    x_prime << X(0,0) << X(1,0) << X(0,1) << X(1,1)
      << X(0,2) << X(1,2) << X(0,3) << X(1,3);
    x_prime = x_prime.t();
    
    // create matrix (A) of 2n x 6 values from 2D
    fmat x_new;
    uword num_points = X.n_cols;
    
    for (int i=0; i < num_points; i++) {
        double x_3d = W(0,i);  // u
        double y_3d = W(1,i);  // v
        
        fmat tmp_row;
        // insert first row
        tmp_row << x_3d << y_3d << 1 << 0 << 0 << 0;
        x_new.insert_rows(2*i, tmp_row);
        // insert second row
        tmp_row << 0 << 0 << 0 << x_3d << y_3d << 1;
        x_new.insert_rows(2*i+1, tmp_row);
    }
    
    fmat tmp_A = solve(x_new, x_prime); // returns 6 x 1
    
    // reshape so A is 2x3
    fmat A;
    A << tmp_A(0,0) << tmp_A(1,0) << tmp_A(2,0)<< endr
      << tmp_A(3,0) << tmp_A(4,0) << tmp_A(5,0);
    
//    A.print("A: ");
    
    return A;
}
//-----------------------------------------------------------------
// Function to project points using the affine transform
//
// <in>
// W = concatenated matrix of 3D points on the plane (3XN)
// A = 2x3 matrix of affine parameters
//
// <out>
// X = concatenated matrix of 2D projected points in the image (2xN)
fmat myproj_affine(fmat &W, fmat &A) {

    // Fill in the answer here.....
    // 3rd row of W is all zeros because it is on a plane
    // This essentially make 3rd ("Z") column of A zeros
    // Since we're not representing Z in the A matrix, we can simply
    // replace the 3rd row of 0s in W with 1s to represent the tx, ty offset
    
    fmat X;
    fmat tmp_W;
    
    tmp_W << W(0,0) << W(0,1) << W(0,2) << W(0,3) << endr
          << W(1,0) << W(1,1) << W(1,2) << W(1,3) << endr
          << 1 << 1 << 1 << 1;
    
    X = A * tmp_W;

    cout << X.n_rows << " " << X.n_cols << endl;
    
    return X;
}

//-----------------------------------------------------------------
// Function to return the affine warp between 3D points on a plane
//
// <in>
// X = concatenated matrix of 2D projected points in the image (2xN)
// W = concatenated matrix of 3D points on the plane (3XN)
//
// <out>
// H = 3x3 homography matrix
fmat myfit_homography(fmat &X, fmat &W) {
    
    // Fill in the answer here.....
    fmat H;
    
    fmat A_mat;
    uword num_points = X.n_cols;
    
    for (int i=0; i < num_points; i++) {
        double u = W(0,i);   // u
        double v = W(1,i);   // v
        double x = X(0,i);
        double y = X(1,i);
        
        fmat tmp_row;
        // insert first row
        tmp_row << 0 << 0 << 0 << -u << -v << -1 << y*u << y*v << y;
        A_mat.insert_rows(2*i, tmp_row);
        // insert second row
        tmp_row << u << v << 1 << 0 << 0 << 0 << -x*u << -x*v << -x;
        A_mat.insert_rows(2*i+1, tmp_row);
    }
    
    mat U;
    vec s;
    fmat V;
    
    mat _V = conv_to<mat>::from(V);
    mat _A_mat = conv_to<mat>::from(A_mat);
    
    svd(U,s,_V,_A_mat);
    
    V = conv_to<fmat>::from(_V);
    H = reshape(V.col(V.n_cols-1), 3, 3);
    H = H.t();
    H.print("H: ");
    
    return H;
}

//-----------------------------------------------------------------
// Function to project points using the affine transform
//
// <in>
// W = concatenated matrix of 3D points on the plane (3XN)
// H = 3x3 homography matrix
//
// <out>
// X = concatenated matrix of 2D projected points in the image (2xN)
fmat myproj_homography(fmat &W, fmat &H) {
    
    // Fill in the answer here.....
    fmat X;
    fmat tmp_W;
    
    tmp_W << W(0,0) << W(0,1) << W(0,2) << W(0,3) << endr
    << W(1,0) << W(1,1) << W(1,2) << W(1,3) << endr
    << 1 << 1 << 1 << 1;
    
    X = H * tmp_W;

    // divid by last row to get homogeneous coordinate
    X.row(0) = X.row(0) / X.row(2);
    X.row(1) = X.row(1) / X.row(2);
    
    X.shed_row(2);  // get rid of last row

    return X;
}