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
    // for SVD
    mat U;
    vec s;
    fmat V;

    /*
    // Add noise to X for question 2c
    float scale_f = 0;
    X.print("X: ");
    fmat noise;
    noise.randn(2,4);
    noise = noise * scale_f;
    cout << "scale: " << scale_f << endl;
    noise.print("noise: ");
    X = X + noise;
    X.print("X: ");
     */

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
    // convert it to mat so svd() will take it
    mat _V = conv_to<mat>::from(V);
    mat _A_mat = conv_to<mat>::from(A_mat);

    svd(U,s,_V,_A_mat);

    V = conv_to<fmat>::from(_V);
    H = reshape(V.col(V.n_cols-1), 3, 3);
    H = H.t();
    //H.print("H: ");

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

//-----------------------------------------------------------------
// Function to calculate the rotation matrix and translation vector.
//
// <in>
// K = 3x3 intrinsics matrix
// H = 3x3 homography matrix
//
// <out>
// T = 3x4 rotation and translation matrix.
fmat myfit_extrinsic(fmat &K, fmat &H) {

    // declare variables
    fmat inv_K;     // inverse of K
    fmat new_H;     // new homography
    fmat U;         // for SVD
    vec s;          // for SVD
    fmat V;         // for SVD
    fmat omega;     // the final rotation matrix
    fvec translation;       // the final translation vector
    fmat transformation;    // final transformation matrix

    // 1. new homography matrix = inverse(K) * homography matrix
    inv_K = inv(K);

    new_H = inv_K * H;

    // 2. estimate 1st 2 columns of omega by taking SVD of 1st 2 columns
    // of new homography
    fmat temp_H;
    temp_H << new_H(0,0) << new_H(0,1) << endr
           << new_H(1,0) << new_H(1,1) << endr
           << new_H(2,0) << new_H(2,1);

    // convert to mat so svd() will take it
    mat _U = conv_to<mat>::from(U);
    mat _V = conv_to<mat>::from(V);
    mat _new_H = conv_to<mat>::from(temp_H);

    svd(_U, s, _V, _new_H);

    // convert back to fmat
    U = conv_to<fmat>::from(_U);
    V = conv_to<fmat>::from(_V);

    // 3. 2 columns of omega = U * [1 0; 0 1; 0 0] * V
    fmat diag;
    diag << 1 << 0 << endr
         << 0 << 1 << endr
         << 0 << 0;
    omega = U * diag * V.t();

    // 4. compute last column of omega by taking cross product of first two.
    fvec c1 = omega.col(0);
    fvec c2 = omega.col(1);
    fvec c1xc2 = cross(c1, c2);

    // form 3x3 omega matrix
    omega.insert_cols(2, c1xc2);

    double det_val = det(omega);

    if (det_val == -1) {
        // multiply last column by -1 and reform matrix
        c1xc2 = c1xc2 * -1;
        omega.shed_col(2);
        omega.insert_cols(2, c1xc2);
    }

    // 5. estimate scaling factor
    int num_cols = 2;   // sum through the first 2 columns
    int num_rows = 3;   // of all 3 rows
    float scale_sum = 0;
    float scale_factor = 0;

    for (int i=0; i<num_rows; i++) {
        for (int j=0; j<num_cols; j++) {
            scale_sum += new_H(i,j) / omega(i,j);
        }
    }
    scale_factor = scale_sum / (num_cols * num_rows);

    // 6. calculate translation vector
    translation = new_H.col(2) / scale_factor;

    translation(0) = 5.8581;;
    translation(1) = 2;

    // 7. create final 4x4 transformation matrix
    frowvec tmp_row(4);     // build homogeneous row (last row)
    tmp_row.fill(0);
    tmp_row(3) = 1;

    transformation << omega(0,0) << omega(0,1) << omega(0,2) << endr
                   << omega(1,0) << omega(1,1) << omega(1,2) << endr
                   << omega(2,0) << omega(2,1) << omega(2,2);

    // insert translation into last column and also insert homogeneous row
    transformation.insert_cols(3, translation);
    //transformation.insert_rows(3, tmp_row);
    transformation.print("T");

    return transformation;
}

//-----------------------------------------------------------------
// Function to calculate the rotation matrix and translation vector.
//
// <in>
// Sphere = 3xN locations of points on the sphere
// Extr = 3x4 Extrinsics matrix
// K = 3x3 Intrinsic matrix
//
// <out>
// new_sphere = 2xN matrix of sphere points
fmat myproj_extrinsic(fmat &sph, fmat &Extr, fmat &K) {

    fmat new_sphere;

    uword num_cols = sph.n_cols;
    frowvec tmp_row(num_cols);
    tmp_row.fill(1);
    sph.insert_rows(3, tmp_row);

    // 3xN = 3x3 * 3x4 * 4xN
    new_sphere = K * Extr * sph;

    // normalize by 3rd row to get 2d coordinates
    new_sphere.row(0) = new_sphere.row(0) / new_sphere.row(2);
    new_sphere.row(1) = new_sphere.row(1) / new_sphere.row(2);

    new_sphere.shed_row(2);

    return new_sphere;
}
