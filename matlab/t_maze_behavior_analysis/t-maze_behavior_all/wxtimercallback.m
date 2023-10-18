%% wxtimercallback
 % Timer callback function to detect changes in subsequent images frames
 % JL Alatorre Warren

function wxtimercallback(obj, ...
                         event, ...
                         cameraNumber, ...
                         previousFrame, ...
                         roi01, ...
                         roi02)
                       
global tmr
global aa
global bb

% Capture a single monochrome frame
currentFrame = LucamCaptureMonochromeFrame(cameraNumber);

% Convert |currentFrame| to DIPimage format
currentFrame = dip_image(currentFrame);

% Image masking
previousFrameMaskedWithRoi01 = previousFrame*roi01;
previousFrameMaskedWithRoi02 = previousFrame*roi02;
currentFrameMaskedWithRoi01  = currentFrame *roi01;
currentFrameMaskedWithRoi02  = currentFrame *roi02;

differenceInRoi01 = previousFrameMaskedWithRoi01 - currentFrameMaskedWithRoi01;
differenceInRoi02 = previousFrameMaskedWithRoi02 - currentFrameMaskedWithRoi02;

differenceInRoi01Score = sum(abs(differenceInRoi01));
differenceInRoi02Score = sum(abs(differenceInRoi02));

display(differenceInRoi01Score)
display(differenceInRoi02Score)

if differenceInRoi01Score > 5000
  aa = 1;
  stop(tmr)
else
  aa = 0;
end

if differenceInRoi02Score > 5000
  bb = 1;
  stop(tmr)
else
  bb = 0;
end