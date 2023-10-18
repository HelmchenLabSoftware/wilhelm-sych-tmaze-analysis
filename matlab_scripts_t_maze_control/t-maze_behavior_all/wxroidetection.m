%% wxroidetection
 % Detects running mice in a given ROI (region of interest)
 % JL Alatorre-Warren

function mainLogbook = wxroidetection(cameraNumber,roi,roiName,mainLogbook)

% Select appropriate ROI
switch(roiName)
  case 'roi01'
    roiName = 'start';
  case 'roi02'
    roiName = 'corridor';
  case 'roi03'
    roiName = 'right';
  case 'roi04'
    roiName = 'left';
  case 'roi06'
    roiName = 'vertex';    
end

% Crop ROI and get ROI edges
[roiCropped,...
 columnIndexMin,...
 columnIndexMax,...
 rowIndexMin,...
 rowIndexMax] = wxcroproi(roi);

% Capture a single monochrome frame
% Convert it to DIPimage format
% Crop the image
framePrevious = LucamCaptureMonochromeFrame(cameraNumber);
framePrevious = dip_image(framePrevious);
framePreviousCropped = framePrevious(columnIndexMin:columnIndexMax,rowIndexMin:rowIndexMax);
clear framePrevious

ii = 0;
endLoopCondition = 1000000;
while ii < endLoopCondition

  % Get a new frame
  % Convert it to DIPimage format
  % Crop the image
  frameCurrent = LucamCaptureMonochromeFrame(cameraNumber);
  frameCurrent = dip_image(frameCurrent);
  frameCurrentCropped = frameCurrent(columnIndexMin:columnIndexMax,rowIndexMin:rowIndexMax);

  % Image masking
  previousFrameMaskedWithRoi = framePreviousCropped*roiCropped;
  currentFrameMaskedWithRoi  = frameCurrentCropped*roiCropped;

  % Compute differences between the current and previous frames
  % Compute an score of the sum of the absolute differences
  differenceInRoi = previousFrameMaskedWithRoi - currentFrameMaskedWithRoi;
  sumOfAbsoluteDifferencesInRoi = sum(abs(differenceInRoi));

  % Detect changes in region of interest ROI
  % Note: Increase the sensitivity factor to make the detection less sensitive
  % Create timestamp if the mouse is detected
  sensitivityFactor = 10.0;
  detectionThreshold = sensitivityFactor*sum(roiCropped);
  if sumOfAbsoluteDifferencesInRoi > detectionThreshold
    ii = endLoopCondition;
    timestampVector = clock;
    mainLogbook = wxupdatelogbook(mainLogbook, ['roi' '_' roiName], timestampVector);
    clear frameCurrent
  end
  
  % Ratio between the actual measurement and the detection threshold
  detectionScore = sumOfAbsoluteDifferencesInRoi/detectionThreshold;
  
  % Display detection score and ROI name
  disp(['Detection score (' roiName '): ' num2str(detectionScore)])

  % The current frame becomes the previous one
  framePreviousCropped = frameCurrentCropped;

  ii = ii+1;
end