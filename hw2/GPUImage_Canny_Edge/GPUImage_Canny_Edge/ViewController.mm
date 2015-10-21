//
//  ViewController.m
//  GPUImage_Canny_Edge
//
//  Created by Johnny Wang on 10/9/15.
//  Copyright © 2015 CV_Apps. All rights reserved.
//

#import "ViewController.h"
#import <GPUImage/GPUImage.h>

#ifdef __cplusplus
#include <stdlib.h> // include the standard library
#endif

@interface ViewController () {
    // Setup the view (this time using GPUImageView)
    GPUImageView *imageView_;
}

@end

@implementation ViewController

- (void)viewDidLoad {
    [super viewDidLoad];
    // Do any additional setup after loading the view, typically from a nib.
    
    // Setup GPUImageView (note we are not using UIImageView here)...
    imageView_ = [[GPUImageView alloc] initWithFrame:CGRectMake(0.0, 0.0, self.view.frame.size.width, self.view.frame.size.height)];
    
    // Important: add as a subview
    [self.view addSubview:imageView_];
    
    // Read in the image (of the book)
    UIImage *inputImage = [UIImage imageNamed:@"prince_book.jpg"];
    
    // Initialize filters
    GPUImagePicture *stillImageSource = [[GPUImagePicture alloc] initWithImage:inputImage smoothlyScaleOutput:YES];
    GPUImageGrayscaleFilter *grayscaleFilter = [[GPUImageGrayscaleFilter alloc] init];
    GPUImageLanczosResamplingFilter *scaleFilter = [[GPUImageLanczosResamplingFilter alloc] init];
    GPUImageGaussianBlurFilter *gausFilter = [[GPUImageGaussianBlurFilter alloc] init];
    GPUImageCannyEdgeDetectionFilter *cannyEdgeFilter = [[GPUImageCannyEdgeDetectionFilter alloc] init];
    
    CGSize imgSize = inputImage.size;
    NSLog(@"width = %f, height = %f", imgSize.width, imgSize.height); // 2448 x 3264
    float down_scale = 0.17;
    float blur_radius = 1;  // in pixels
    float width_scale = imgSize.width * down_scale;
    float height_scale = imgSize.height * down_scale;
    float lower_thresh = 0.11;
    float upper_thresh = lower_thresh * 3;
    
    // Set filter variables
    // scale/resize image
    [scaleFilter forceProcessingAtSizeRespectingAspectRatio:CGSizeMake(width_scale,height_scale)];
    // blur
    [gausFilter setBlurRadiusInPixels:blur_radius];
    [cannyEdgeFilter setBlurRadiusInPixels:blur_radius];
    [cannyEdgeFilter setLowerThreshold:lower_thresh];
    [cannyEdgeFilter setUpperThreshold:upper_thresh];
    
    // Daisy chain the filters together (you can add as many filters as you like)
    [stillImageSource addTarget:grayscaleFilter];
    [grayscaleFilter addTarget:scaleFilter];
    [scaleFilter addTarget:gausFilter];
    [gausFilter addTarget:cannyEdgeFilter];
//    [scaleFilter addTarget:cannyEdgeFilter];
    [cannyEdgeFilter addTarget:imageView_];
    
    // Process the image
    [stillImageSource processImage];
}

- (void)didReceiveMemoryWarning {
    [super didReceiveMemoryWarning];
    // Dispose of any resources that can be recreated.
}

@end
