//
//  ViewController.m
//  Canny_Edge
//
//  Created by Johnny Wang on 10/9/15.
//  Copyright © 2015 CV_Apps. All rights reserved.
//

#import "ViewController.h"

#ifdef __cplusplus
#include <opencv2/opencv.hpp> // include opencv library
#include <stdlib.h> // include standard library
#endif

using namespace std;

@interface ViewController () {
    // Setup the view
    UIImageView *imageView_;
}
@end

@implementation ViewController

- (void)viewDidLoad {
    [super viewDidLoad];
    // Do any additional setup after loading the view, typically from a nib.
    
    // Load the prince book image
    UIImage *image = [UIImage imageNamed:@"prince_book.jpg"];
    if (image == nil) cout << "Cannot read in the file prince_book.jpg!!" << endl;
    
    // Setup the display
    // Setup the imageView_ view so it takes up the entire App screen...
    imageView_ = [[UIImageView alloc] initWithFrame:CGRectMake(0.0, 0.0, self.view.frame.size.width, self.view.frame.size.height)];
    // Important: add OpenCV_View as a subview
    [self.view addSubview:imageView_];
    
    // Ensure aspect ratio looks correct
    imageView_.contentMode = UIViewContentModeScaleAspectFit;
    
    cv::Mat cvImage = [self cvMatFromUIImage:image];
    cv::Mat gray;
    cv::cvtColor(cvImage, gray, CV_RGBA2GRAY);  // convert to grayscale
    cv::Mat display_im;
    
    // apply Canny filter
    int lowThreshold = 25;   // vary this
    int ratio = 3;
    int kernel_size = 3; // 3
    float f_scale = 0.10;
    float blur_val = 5.9;
    cv::Mat detected_edges;
    // Make gray image smaller using scale factors
    cv::resize(gray, gray, cv::Size(), f_scale, f_scale);
    cv::GaussianBlur(gray, detected_edges, cv::Size(blur_val,blur_val),1.5,1.5);  // 5.5 5.5 1.5
    cv::Canny(detected_edges, detected_edges, lowThreshold, lowThreshold*ratio, kernel_size);
    // Uisng Canny's output as a mask, display result
    cv::Mat dst;
    dst.create(gray.size(), gray.type());
    dst = cv::Scalar::all(0);
    gray.copyTo(dst, detected_edges);   // copy gray to dst while applying the mask
    
    // Display the view
    cv::cvtColor(dst, display_im, CV_GRAY2BGR);    // get the display image
    imageView_.image = [self UIImageFromCVMat:display_im];
    
    // Calculate the HoughLines
    cv::Mat cdst;
    vector<cv::Vec2f> lines;
    cv::cvtColor(dst, cdst, CV_GRAY2BGR);  // convert to color so we can draw red lines

    HoughLines(dst, lines, 1, CV_PI/180, 100, 0, 0 );
    
    for( size_t i = 0; i < lines.size(); i++ )
    {
        float rho = lines[i][0], theta = lines[i][1];
        cv::Point pt1, pt2;
        double a = cos(theta), b = sin(theta);
        double x0 = a*rho, y0 = b*rho;
        pt1.x = cvRound(x0 + 1000*(-b));
        pt1.y = cvRound(y0 + 1000*(a));
        pt2.x = cvRound(x0 - 1000*(-b));
        pt2.y = cvRound(y0 - 1000*(a));
        line( cdst, pt1, pt2, cv::Scalar(0,0,255), 3, CV_AA);   // draw red line
    }
    // Switch colors to account for how UIImage and cv::Mat lay out their color channels differently
    cv::cvtColor(cdst, cdst, CV_BGRA2RGBA);
    imageView_.image = [self UIImageFromCVMat:cdst];

}

- (void)didReceiveMemoryWarning {
    [super didReceiveMemoryWarning];
    // Dispose of any resources that can be recreated.
}

//----------------------------------------------------------------------------------------------------
// Additional functions
//----------------------------------------------------------------------------------------------------
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
