%% wxgetrois
 % Displays the location of all ROIs (regions of interest)
 % JL Alatorre-Warren

function [roi01, infoRoi01, ...
          roi02, infoRoi02, ...
          roi03, infoRoi03, ...
          roi04, infoRoi04, ...
          roi06, infoRoi06] = wxgetrois(cameraNumber, pathRois) 
 
% Capture a single monochrome frame
frameMonochrome = LucamCaptureMonochromeFrame(cameraNumber);

% Open ROIs
[roi01, infoRoi01] = readim([pathRois, '/', 'roi01.tif'],'TIFF');
[roi02, infoRoi02] = readim([pathRois, '/', 'roi02.tif'],'TIFF');
[roi03, infoRoi03] = readim([pathRois, '/', 'roi03.tif'],'TIFF');
[roi04, infoRoi04] = readim([pathRois, '/', 'roi04.tif'],'TIFF');
[roi06, infoRoi06] = readim([pathRois, '/', 'roi06.tif'],'TIFF');

% Union of all ROI (logical OR)
% Here, roi05 is not included on purpose. So, it's not an error.
roiAll = roi01 | roi02 | roi03 | roi04 | roi06;

% Convert |frameMonochrome| to DIPimage format
frameMonochrome = dip_image(frameMonochrome);

% Join channels: |frameMonochrome| and |roiAll|
frameMonochrome_roiAll = joinchannels('RGB', ... 
                                     roiAll*0, ...
                                     roiAll/max(roiAll)*100, ...
                                     frameMonochrome/max(frameMonochrome)*255);

% Display joint channels
% Figure setup
dipshow(frameMonochrome_roiAll)
figureRoiAll = gcf;
figureRoiAll.Position = [1289 49 853 635];