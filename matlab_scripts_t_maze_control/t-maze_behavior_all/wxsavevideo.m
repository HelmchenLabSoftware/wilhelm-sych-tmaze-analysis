%% wxsavevideo
 % Captures multiple video frames using the camera Lumenera
 % JL Alatorre Warren

% Preparing DIPimage and the camera Lumenera
dipstart
cameraNumber = 1; 
LucamShowPreview(cameraNumber);

% Set frame rate and number of frames
LucamSetFrameRate(20, 1)
numberOfFrames = 1000;

% Actual video acquisition
tic
frames = LucamTakeVideo(numberOfFrames, 1);
toc

% Convert to DIPimage format
% Convert to monochrome image
% Reduce dimensionality
% Crop field of view
% Display
frames = dip_image(frames);
frames = frames(:,:,1,:);
framesSqueeze = squeeze(frames);
framesSqueezeCrop = framesSqueeze(300:860,0:765,:);
dipshow(framesSqueezeCrop)

% Convert back to MATLAT format
framesSqueezeCropUint8 = uint8(framesSqueezeCrop);

% Write video
videoPath = 'C:/Users/chernysheva/Documents/MATLAB/data';
videoName = 'run.avi';
videoObj = VideoWriter([videoPath '/' videoName]);
videoObj.FrameRate = 7.5;
open(videoObj);
for k = 1:numberOfFrames
  writeVideo(videoObj,framesSqueezeCropUint8(:,:,k));
end
close(videoObj);

toc