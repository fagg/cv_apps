//
//  ViewController.m
//  Estimate_PlanarWarps
//
//  Created by Simon Lucey on 9/21/15.
//  Copyright (c) 2015 CMU_16432. All rights reserved.
//

#import "ViewController.h"
#import "myfit.h"

#ifdef __cplusplus
#include <opencv2/opencv.hpp> // Includes the opencv library
#include <stdlib.h> // Include the standard library
#include "armadillo" // Includes the armadillo library
#include <opencv2/nonfree/features2d.hpp>
#include <opencv2/nonfree/nonfree.hpp>
#endif

using namespace std;

@interface ViewController () {
    // Setup the view
    UIImageView *imageView_1;
}
@end

@implementation ViewController

- (void)viewDidLoad {
    [super viewDidLoad];
    // Do any additional setup after loading the view, typically from a nib.
    
    // 3D planar points of the Prince Computer Vision textbook (in cm)
    arma::fmat W;
    W << 0.0 << 18.2 << 18.2 <<  0.0 << arma::endr
      << 0.0 <<  0.0 << 26.0 << 26.0 << arma::endr
      << 0.0 <<  0.0 <<  0.0 << 0.0;
    
    // Corresponding 2D projected points of the book in the image
    arma::fmat X;
    X << 483 << 1704 << 2175 <<  67 << arma::endr
      << 810 <<  781 << 2217 << 2286;
    
    // Intrinsics matrix for the device it was caputured from...
    arma::fmat K;
    K << 3043.72 <<       0 << 1196 << arma::endr
      <<       0 << 3043.72 << 1604 << arma::endr
      <<       0 <<    0    <<    1;
    
    // Load the 3D sphere points (dimensions of ball are in cm)
    NSString *str = [[NSBundle mainBundle] pathForResource:@"sphere" ofType:@"txt"];
    const char *SphereName = [str UTF8String]; // Convert to const char *
    arma::fmat sphere; sphere.load(SphereName); // Load the Sphere into memory should be 3xN
    
    // Read in the image
    UIImage *image1 = [UIImage imageNamed:@"prince_book.jpg"];
    if(image1 == nil) cout << "Cannot read in the file prince_book.jpg!!" << endl;
    // For Harris Corner
    UIImage *image2 = [UIImage imageNamed:@"new_view.jpg"];
    if(image2 == nil) cout << "Cannot read in the file new_view.jpg!!" << endl;
    
    // Setup the display
    // Setup the your imageView_1 view, so it takes up the entire App screen......
    imageView_1 = [[UIImageView alloc] initWithFrame:CGRectMake(0.0, 0.0, self.view.frame.size.width, self.view.frame.size.height)];
    
    // Another way to convert between cvMat and UIImage (using member functions)
    cv::Mat cvImage = [self cvMatFromUIImage:image1];
    cv::Mat gray; cv::cvtColor(cvImage, gray, CV_RGBA2GRAY); // Convert to grayscale
    cv::Mat display_im; cv::cvtColor(gray,display_im,CV_GRAY2BGR); // Get the display image
    
    
    // Important: add OpenCV_View as a subview
    [self.view addSubview:imageView_1];
    
    // Ensure aspect ratio looks correct
    imageView_1.contentMode = UIViewContentModeScaleAspectFit;

    // Run homography, Harris Corner, or SIFT
    int run_type = 2;
    
    if (run_type == 0) {
        //-------------------------------------------------------
        // Affine Warp / Homography
        
        const cv::Scalar RED = cv::Scalar(0,0,255); // Set the RED color
        const cv::Scalar GREEN = cv::Scalar(0,255,0); // Set the GREEN color
        const cv::Scalar BLUE = cv::Scalar(255,0,0); // Set the BLUE color
        
        int run_warp = 2;
        
        // Default homography warp given
        if (run_warp == 0) {
            display_im = DrawPts(display_im, X, RED);
            display_im = DrawLines(display_im, X, RED);
        } else if (run_warp == 1) {
            // Affine warp
            arma::fmat A = myfit_affine(X, W);
            arma::fmat affine_x = myproj_affine(W, A);
            display_im = DrawPts(display_im, affine_x, GREEN);
            display_im = DrawLines(display_im, affine_x, GREEN);
        } else if (run_warp == 2) {
            // Homography
            arma::fmat Hom = myfit_homography(X, W);
            arma::fmat homo_x = myproj_homography(W, Hom);
            display_im = DrawPts(display_im, homo_x, BLUE);
            display_im = DrawLines(display_im, homo_x, BLUE);
        }
        
        cv::cvtColor(display_im, display_im, CV_BGRA2RGBA);
    } else if (run_type == 1) {
        //--------------------------------------------------------------
        // Harris Corner
        
        // read in new image
        cv::Mat cvImage2 = [self cvMatFromUIImage:image2];
        cv::cvtColor(cvImage2, gray, CV_RGBA2GRAY); // Convert to grayscale
        cv::cvtColor(gray,display_im,CV_GRAY2BGR); // Get the display image
        
        cv::Mat dst, dst_norm_scaled;
        dst = cv::Mat::zeros( gray.size(), CV_32FC1 );
        
        dst_norm_scaled = RunHarrisCorner(dst_norm_scaled, gray, dst);
        
        // Switch colors to account for how UIImage and cv::Mat lay out their color channels differently
        cv::cvtColor(dst_norm_scaled, dst_norm_scaled, CV_GRAY2BGRA);
        cv::cvtColor(dst_norm_scaled, display_im, CV_BGRA2RGBA);
        
        std::cout << display_im.rows << "x" << display_im.cols << std::endl;
        std::cout << dst_norm_scaled.rows << "x" << dst_norm_scaled.cols << std::endl;
    } else if (run_type == 2) {
        //--------------------------------------------------------------
        // SIFT feature detector / descriptor extractor
        
        cv::Mat gray2, display_im2;
        
        // read in new image
        cv::Mat cvImage2 = [self cvMatFromUIImage:image2];
        cv::cvtColor(cvImage2, gray2, CV_RGBA2GRAY); // Convert to grayscale
        cv::cvtColor(gray2,display_im2,CV_GRAY2BGR); // Get the display image
        
        // 1. Detect keypoints using SIFT
        double threshold = 500;
        cv::SiftFeatureDetector detector(threshold);
        vector<cv::KeyPoint> keypoints_1, keypoints_2;
        detector.detect(display_im, keypoints_1);
        detector.detect(display_im2, keypoints_2);
        cout << "Detected " << (int) keypoints_1.size() << " keypoints" <<endl;
        cout << "Detected " << (int) keypoints_2.size() << " keypoints" <<endl;
        
        // 2. Compute feature description.
        cv::SiftDescriptorExtractor extractor;
        cv::Mat descriptors_1, descriptors_2;
        extractor.compute(display_im, keypoints_1, descriptors_1);
        extractor.compute(display_im2, keypoints_2, descriptors_2);

        // 3. Match descriptor vectors using brute force matcher
        cv::BFMatcher matcher(cv::NORM_L2, true);   // true for bidirectional cross_check
        std::vector<cv::DMatch> matches;
        matcher.match(descriptors_1, descriptors_2, matches);
        cout << "matches: " << matches.size() << endl;
        
        // Draw matches
        cv::Mat img_matches;
        cv::drawMatches(display_im, keypoints_1, display_im2, keypoints_2, matches, img_matches);
        
        // Find homography of book cover
        //arma::fmat H = myfit_homography(X, W);
        arma::fmat best_h = my_ransac(W, keypoints_1, keypoints_2, matches);
        
        best_h.print("best_h:");
        
        // display homography of book cover in new_img.jpg
        const cv::Scalar YELLOW = cv::Scalar(0,255,255); // Set the BLUE color
        arma::fmat homo_x = myproj_homography(X, best_h);
        homo_x.print("homo_x:");
        display_im2 = DrawPts(display_im2, homo_x, YELLOW);
        display_im2 = DrawLines(display_im2, homo_x, YELLOW);
        
        cv::cvtColor(display_im2, display_im, CV_BGRA2RGBA);
    }
    
    // Finally setup the view to display
    imageView_1.image = [self UIImageFromCVMat:display_im];

}

- (void)didReceiveMemoryWarning {
    [super didReceiveMemoryWarning];
    // Dispose of any resources that can be recreated.
}

cv::Mat RunHarrisCorner(cv::Mat &dst_norm_scaled, cv::Mat &gray, cv::Mat &dst) {
    // Detector parameters
    int blockSize = 3;      // 3 // 2
    int apertureSize = 19;  // 19 // 3
    double k = 0.065;        // 0.07 // 0.04
    int thresh = 200;   // 200
    
    cv::Mat dst_norm;
    
    /// Detecting corners
    cornerHarris( gray, dst, blockSize, apertureSize, k, cv::BORDER_DEFAULT );
    
    /// Normalizing
    normalize( dst, dst_norm, 0, 255, cv::NORM_MINMAX, CV_32FC1, cv::Mat() );
    convertScaleAbs( dst_norm, dst_norm_scaled );
    
    /// Drawing a circle around corners
    for( int j = 0; j < dst_norm.rows ; j++ ) {
        for( int i = 0; i < dst_norm.cols; i++ ) {
            if( (int) dst_norm.at<float>(j,i) > thresh) {
                circle( dst_norm_scaled, cv::Point(i,j), 5, cv::Scalar(0), 2, 8, 0 );
            }
        }
    }
    
    return dst_norm_scaled;
}

//---------------------------------------------------------------------------------------------------------------------
// You should not have to touch these functions below to complete the assignment!!!!
//---------------------------------------------------------------------------------------------------------------------
// Quick function to draw points on an UIImage
cv::Mat DrawPts(cv::Mat &display_im, arma::fmat &pts, const cv::Scalar &pts_clr)
{
   vector<cv::Point2f> cv_pts = Arma2Points2f(pts); // Convert to vector of Point2fs
   for(int i=0; i<cv_pts.size(); i++) {
       cv::circle(display_im, cv_pts[i], 5, pts_clr,5); // Draw the points
   }
    return display_im; // Return the display image
}
// Quick function to draw lines on an UIImage
cv::Mat DrawLines(cv::Mat &display_im, arma::fmat &pts, const cv::Scalar &pts_clr)
{
    vector<cv::Point2f> cv_pts = Arma2Points2f(pts); // Convert to vector of Point2fs
    for(int i=0; i<cv_pts.size(); i++) {
        int j = i + 1; if(j == cv_pts.size()) j = 0; // Go back to first point at the enbd
        cv::line(display_im, cv_pts[i], cv_pts[j], pts_clr, 3); // Draw the line
    }
    return display_im; // Return the display image
}
// Quick function to convert Armadillo to OpenCV Points
vector<cv::Point2f> Arma2Points2f(arma::fmat &pts)
{
 vector<cv::Point2f> cv_pts;
 for(int i=0; i<pts.n_cols; i++) {
 cv_pts.push_back(cv::Point2f(pts(0,i), pts(1,i))); // Add points
 }
 return cv_pts; // Return the vector of OpenCV points
}
// Member functions for converting from cvMat to UIImage
- (cv::Mat)cvMatFromUIImage:(UIImage *)image
{
    CGColorSpaceRef colorSpace = CGImageGetColorSpace(image.CGImage);
    CGFloat cols = image.size.width;
    CGFloat rows = image.size.height;
    
    cv::Mat cvMat(rows, cols, CV_8UC4); // 8 bits per component, 4 channels (color channels + alpha)
    
    CGContextRef contextRef = CGBitmapContextCreate(cvMat.data,                 // Pointer to  data
                                                    cols,                       // Width of bitmap
                                                    rows,                       // Height of bitmap
                                                    8,                          // Bits per component
                                                    cvMat.step[0],              // Bytes per row
                                                    colorSpace,                 // Colorspace
                                                    kCGImageAlphaNoneSkipLast |
                                                    kCGBitmapByteOrderDefault); // Bitmap info flags
    
    CGContextDrawImage(contextRef, CGRectMake(0, 0, cols, rows), image.CGImage);
    CGContextRelease(contextRef);
    
    return cvMat;
}
// Member functions for converting from UIImage to cvMat
-(UIImage *)UIImageFromCVMat:(cv::Mat)cvMat
{
    NSData *data = [NSData dataWithBytes:cvMat.data length:cvMat.elemSize()*cvMat.total()];
    CGColorSpaceRef colorSpace;
    
    if (cvMat.elemSize() == 1) {
        colorSpace = CGColorSpaceCreateDeviceGray();
    } else {
        colorSpace = CGColorSpaceCreateDeviceRGB();
    }
    
    CGDataProviderRef provider = CGDataProviderCreateWithCFData((__bridge CFDataRef)data);
    
    // Creating CGImage from cv::Mat
    CGImageRef imageRef = CGImageCreate(cvMat.cols,                                 //width
                                        cvMat.rows,                                 //height
                                        8,                                          //bits per component
                                        8 * cvMat.elemSize(),                       //bits per pixel
                                        cvMat.step[0],                            //bytesPerRow
                                        colorSpace,                                 //colorspace
                                        kCGImageAlphaNone|kCGBitmapByteOrderDefault,// bitmap info
                                        provider,                                   //CGDataProviderRef
                                        NULL,                                       //decode
                                        false,                                      //should interpolate
                                        kCGRenderingIntentDefault                   //intent
                                        );
    
    
    // Getting UIImage from CGImage
    UIImage *finalImage = [UIImage imageWithCGImage:imageRef];
    CGImageRelease(imageRef);
    CGDataProviderRelease(provider);
    CGColorSpaceRelease(colorSpace);
    
    return finalImage;
}

@end
