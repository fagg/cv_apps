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
// Function to return the homography warp between 3D points on a plane
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
//    H.print("H: ");
    
    return H;
}

//-----------------------------------------------------------------
// Function to project points using the homography transform
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
    fmat tmp_W(W);
    
//    W.print("W:");

    // replace last row of 0s with 1s
    /*
    tmp_W << W(0,0) << W(0,1) << W(0,2) << W(0,3) << endr
    << W(1,0) << W(1,1) << W(1,2) << W(1,3) << endr
    << 1 << 1 << 1 << 1;
     */
    fmat tmp_row(1, W.n_cols, fill::ones);
    if (W.n_rows == 3) {
        tmp_W.row(2) = tmp_row;
    } else {
        tmp_W.insert_rows(2, tmp_row);
    }
    
//    tmp_W.print("tmp_W:");
    
    X = H * tmp_W;

    // divid by last row to get homogeneous coordinate
    X.row(0) = X.row(0) / X.row(2);
    X.row(1) = X.row(1) / X.row(2);
    
    X.shed_row(2);  // get rid of last row

    return X;
}

//-----------------------------------------------------------------
// Function to project points using the affine transform
//
// <in>
// W = concatenated matrix of 3D points on the plane (3XN)
// H = 3x3 homography matrix
//
// <out>
// best_h = 3x3 homography matrix with the most inliers
fmat my_ransac(fmat &W,
               vector<cv::KeyPoint> kp1,
               vector<cv::KeyPoint> kp2,
               vector<cv::DMatch> matches) {

    int min = 0;
    int max = matches.size();
    int num_points = 4;         // number of points to get from keypoints
    double r_thresh = 10;       // ransac threshold for point distance
    int ransac_iterations = 10000;  // loops of ransac to run
    int max_inliers = 0;        // max inliers we've found so far
    fmat best_h;                // best homography with most inliers
    
    // Uniform random number generator
    std::random_device rd;      // only used once to initialize (seed) engine
    std::mt19937 rng(rd());     // random-number engine (Mersenne-Twister)
    std::uniform_int_distribution<int> uni(min,max-1);
    
    // turn kp1 and kp2 to fmat so we can find homography easier
    fmat keypoints1;
    fmat keypoints2;
    for (int i=0; i<kp1.size(); i++) {
        fmat tmp_col;
        tmp_col << kp1[i].pt.x << endr << kp1[i].pt.y;
        keypoints1.insert_cols(i, tmp_col);
    }
    for (int i=0; i<kp2.size(); i++) {
        fmat tmp_col;
        tmp_col << kp2[i].pt.x << endr << kp2[i].pt.y;
        keypoints2.insert_cols(i, tmp_col);
    }
    // fill ones row so we can run homography with this
    fmat kp1_last_row(1, kp1.size(), fill::ones);
    keypoints1.insert_rows(2, kp1_last_row);

    // To reject points outside of the book cover
    // the projected 2D coordinates of the book cover
 //   arma::fmat homo_x = myproj_homography(W, H);
    
    for (int r_iter=0; r_iter<ransac_iterations; r_iter++) {
        // number of ransac inliers for this iteration
        int iter_inliers = 0;
        fmat new_kp1_points;      // help build the new homography
        fmat new_kp2_points;      // help build the new homography
        fmat new_h;             // new homography to test each point on
    
        // find 4 pairs of matching points and calculating the homography between them
        for (int i = 0; i<num_points; i++) {
            int idx = uni(rng);
            fmat tmp_col1;
            fmat tmp_col2;
            
            uword kp1_idx = matches[idx].queryIdx;
            uword kp2_idx = matches[idx].trainIdx;

            // create column vector and append to our new_x and new_w
            tmp_col1 << kp1[kp1_idx].pt.x << endr << kp1[kp1_idx].pt.y << endr << 1;
            new_kp1_points.insert_cols(i, tmp_col1);
            tmp_col2 << kp2[kp2_idx].pt.x << endr << kp2[kp2_idx].pt.y << endr << 1;
            new_kp2_points.insert_cols(i, tmp_col2);
        }
        new_h = myfit_homography(new_kp2_points, new_kp1_points);  // 3x3 homography
        
        // run through each correspondence and see if new_h gives good correspondence
        for (int i=0; i<matches.size(); i++) {
            uword kp1_idx = matches[i].queryIdx;
            uword kp2_idx = matches[i].trainIdx;
            
            // needs to be 3xN
            
            fmat tmp_kp1;
            tmp_kp1 << kp1[kp1_idx].pt.x << endr << kp1[kp1_idx].pt.y << endr << 1;
            fmat proj_kp2 = myproj_homography(tmp_kp1, new_h);
             
            double dist_x = proj_kp2(0,0) - kp2[kp2_idx].pt.x;
            double dist_y = proj_kp2(1,0) - kp2[kp2_idx].pt.y;
            double euc_dist = sqrt(dist_x * dist_x + dist_y * dist_y);
            
            if (euc_dist < r_thresh) {
                iter_inliers++;
            }
        }
        
        if (iter_inliers > max_inliers) {
            max_inliers = iter_inliers;
            best_h = new_h;
        }
        
    }
    cout << "max_inliers: " << max_inliers << endl;
    
    return best_h;
}