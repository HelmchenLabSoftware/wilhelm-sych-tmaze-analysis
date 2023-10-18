%% wxtworoisdetection
 % Detects running mice in either the right or the left arm of the T-maze
 % JL Alatorre-Warren

function [detectedSide, ...
          mainLogbook] = wxtworoisdetection(cameraNumber, ...
                                            mainLogbook, ...
                                            roiRight, ...
                                            roiLeft)

% Union of both ROIs (logical OR)
roiBoth = roiRight | roiLeft;

% Crop ROI and get ROI edges
[roiBothCropped,...
 columnIndexMin,...
 columnIndexMax,...
 rowIndexMin,...
 rowIndexMax] = wxcroproi(roiBoth);
                                          
% Capture a single monochrome frame
% Convert it to DIPimage format
% Crop the image with the same cropping parameters as in roiBothCropped
framePrevious = LucamCaptureMonochromeFrame(cameraNumber);
framePrevious = dip_image(framePrevious);
framePreviousCropped = framePrevious(columnIndexMin:columnIndexMax,rowIndexMin:rowIndexMax);
clear framePrevious

ii = 0;
endLoopCondition = 1000000;
while ii < endLoopCondition

  % Get a new frame and convert it to DIPimage format
  frameCurrent = LucamCaptureMonochromeFrame(cameraNumber);
  frameCurrent = dip_image(frameCurrent);
  frameCurrentCropped = frameCurrent(columnIndexMin:columnIndexMax,rowIndexMin:rowIndexMax);

  % Image masking
  previousFrameMaskedWithRoiBoth = framePreviousCropped*roiBothCropped;
  currentFrameMaskedWithRoiBoth  = frameCurrentCropped*roiBothCropped;

  % Compute differences between the current and previous frames
  % Compute an score of the sum of the absolute differences
  differenceInRoiBoth = previousFrameMaskedWithRoiBoth - currentFrameMaskedWithRoiBoth;
  sumOfAbsoluteDifferencesInRoiBoth = sum(abs(differenceInRoiBoth));

  % Detect changes in region of interest ROI
  % Increase the sensitivity factor to make the detection less sensitive
  sensitivityFactor = 7.0;
  detectionThreshold = sensitivityFactor*sum(roiBoth);

  % Actual detection and registration of timestamp
  if sumOfAbsoluteDifferencesInRoiBoth > detectionThreshold
    ii = endLoopCondition;
    timestampVector = clock;
  else
    framePreviousCropped = frameCurrentCropped;
  end
  
  % Ratio between the actual measurement and the detection threshold
  detectionScore = sumOfAbsoluteDifferencesInRoiBoth/detectionThreshold;
  
  % Display detection score and ROI name
  disp(['Detection score (2 ROIs): ' num2str(detectionScore)])

  ii = ii+1;
end

% Crop right and left ROIs using the same cropping parameters as in roiBothCropped
roiRightCropped = roiRight(columnIndexMin:columnIndexMax,rowIndexMin:rowIndexMax);
roiLeftCropped = roiLeft(columnIndexMin:columnIndexMax,rowIndexMin:rowIndexMax);

% Right and left ROIs separately
% Image masking
previousFrameMaskedWithRoiRight = framePreviousCropped*roiRightCropped;
previousFrameMaskedWithRoiLeft  = framePreviousCropped*roiLeftCropped;
currentFrameMaskedWithRoiRight  = frameCurrentCropped*roiRightCropped;
currentFrameMaskedWithRoiLeft   = frameCurrentCropped*roiLeftCropped;

% Right and left ROIs separately
% Compute differences between the current and previous frames
% Compute an score of the sum of the absolute differences
differenceInRoiRight = previousFrameMaskedWithRoiRight - currentFrameMaskedWithRoiRight;
differenceInRoiLeft  = previousFrameMaskedWithRoiLeft  - currentFrameMaskedWithRoiLeft;
sumOfAbsoluteDifferencesInRoiRight = sum(abs(differenceInRoiRight));
sumOfAbsoluteDifferencesInRoiLeft  = sum(abs(differenceInRoiLeft));

% Right and left ROIs separately
% Detect changes in region of interest ROI
% A sensitivity factor of 3 seems to offer excellent detection
% Increase the sensitivity factor to make the detection more sensitive
sensitivityFactor = 3;
detectionThresholdRight = sensitivityFactor*sum(roiRight);
detectionThresholdLeft  = sensitivityFactor*sum(roiLeft);

% Right and left ROIs separately
% Ratio between the actual measurement and the detection threshold
detectionScoreRight = sumOfAbsoluteDifferencesInRoiRight/detectionThresholdRight;
detectionScoreLeft  = sumOfAbsoluteDifferencesInRoiLeft/detectionThresholdLeft;

if detectionScoreRight > detectionScoreLeft
  detectedSide = 'right';
else
  detectedSide = 'left';
end

% Display detection score and ROI name
disp(['Detection score (right): ' num2str(detectionScoreRight)])
disp(['Detection score (left): ' num2str(detectionScoreLeft)])
disp(['The mouse went to the ' detectedSide])

% Update logbook
mainLogbook = wxupdatelogbook(mainLogbook, ['roi' '_' detectedSide], timestampVector);