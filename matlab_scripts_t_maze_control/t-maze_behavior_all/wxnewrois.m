%% wxnewrois
 % This script allows the user to manually draw and see up to 6 ROIs
 % regions of interest
 % JL Alatorre Warren

%% Directory paths
pathRoot = 'I:';
pathData = [pathRoot, '/', '2017_m001-m006_backup'];
pathRois = [pathData, '/', 'rois'];

%% DIPimage and DIPlib
%{
  The function |dipstart| initializes DIPimage (MATLAB)
  DIPimage is a MATLAB toolbox for scientific image processing and
  analysis.
  Most image processing operations are deferred to DIPlib, a dedicated
  library written in C. In this way we extend MATLAB with new core
  functions, and thereby overcoming speed limitations of interpreted
  scripts.
  DIPlib is a platform independent scientific image processing library
  written in C. It consists of a large number of functions for processing
  and analysing multi-dimensional image data. The library provides
  functions for performing transforms, filter operations, object
  generation, local structure analysis, object measurements and
  statistical analysis of images.
%}
dipstart;
 
 %% Test Lumenera camera

% Assign camera number
% Display video preview
cameraNumber = 1; 
LucamShowPreview(cameraNumber);

%% Image acquisition, processing and analysis

% Capture and display a single monochrome frame
frameMonochrome = LucamCaptureMonochromeFrame(cameraNumber);
figure('Name','MonochromeColors');
imagesc(frameMonochrome, [0,255]);

% Draw and display regions of interest (ROIs)
% These ROIs are in DIPimage format
[roi01, roi02, roi03, roi04, roi05, roi06] = wxdrawrois(6,frameMonochrome);

% Save ROIs
writeim(roi01,[pathRois, '/', 'roi01.tif'],'TIFF','no')
writeim(roi02,[pathRois, '/', 'roi02.tif'],'TIFF','no')
writeim(roi03,[pathRois, '/', 'roi03.tif'],'TIFF','no')
writeim(roi04,[pathRois, '/', 'roi04.tif'],'TIFF','no')
writeim(roi05,[pathRois, '/', 'roi05.tif'],'TIFF','no')
writeim(roi06,[pathRois, '/', 'roi06.tif'],'TIFF','no')

% Union of all ROI (logical OR)
% Here, roi05 is not included on purpose. So, it's not an error.
roiAll = roi01 | roi02 | roi03 | roi04 | roi06;

% Convert |frameMonochrome| to DIPimage format
frameMonochrome = dip_image(frameMonochrome);

% Join channels: |frameMonochrome| and |roiAll|
frameMonochrome_roiAll = joinchannels('RGB', ... 
                                     roi05/max(roi05)*100, ...
                                     roiAll/max(roiAll)*100, ...
                                     frameMonochrome/max(frameMonochrome)*255);

% Display joint channels
dipshow(frameMonochrome_roiAll)