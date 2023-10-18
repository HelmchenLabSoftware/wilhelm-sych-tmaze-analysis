%% wxtestservos
 % This script tests the Lumenera camera
 % JL Alatorre Warren


% Assign camera number
% Display video preview
cameraNumber = 1; 
LucamShowPreview(cameraNumber);

% Capture and display one single frame
frameStandard = LucamCaptureFrame(cameraNumber);
figure('Name','LucamCaptureFrame test');
image(frameStandard);

% Capture and display a single monochrome frame
frameMonochrome = LucamCaptureMonochromeFrame(cameraNumber);
figure('Name','LucamCaptureMonochromeFrame test');
imagesc(frameMonochrome, [0,255]); colormap(gray);

% Capture and display multiple frames
frameMultiple = LucamCaptureMultipleFrames(8, cameraNumber);
figure('Name','LucamCaptureMultipleFrame Test')
imaqmontage(frameMultiple)

% Display a specific (third) frame within multiple_frames
image(frameMultiple(:,:,:,3))

% Capture and display a single monochrome frame
frameRaw= LucamCaptureRawFrame(cameraNumber);
figure('Name','LucamCaptureRawFrame test');
imagesc(frameRaw,[0,255]);