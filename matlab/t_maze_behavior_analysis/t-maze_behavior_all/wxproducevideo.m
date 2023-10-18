%% wxproducevideo
 % Produce video using the MAT-files acquired using the 
 % camera Logitech (USB Video Device)
 % JL Alatorre-Warren
 
function wxproducevideo(pathFrames)

% Create a list with all the frames of the experiment
structureOfFrames   = dir(pathFrames);
listOfFrames        = ({structureOfFrames(3:end).name})';
totalNumberOfFrames = size(listOfFrames);
totalNumberOfFrames = totalNumberOfFrames(1);

% A maximum number of frames per video is used because MATLAB does not
% allow to create a single large video. So, if there are too many frames,
% more than one video will be created.
maxNumberOfFramesPerVideo = 20000;

% Calculate the total number of videos
numberOfFramesInLastVideo = rem(totalNumberOfFrames,maxNumberOfFramesPerVideo);
if numberOfFramesInLastVideo ~= 0
  numberOfVideos = floor(totalNumberOfFrames/maxNumberOfFramesPerVideo)+1;
elseif numberOfFramesInLastVideo == 0
  numberOfVideos = floor(totalNumberOfFrames/maxNumberOfFramesPerVideo);  
end

% Load structure currentFrame
currentFrame = load([pathFrames '/' char(listOfFrames(1))], 'currentFrame');

% Convert the structure currentFrame to a variable with the same name (currentFrame)
currentFrame = currentFrame.currentFrame;

% Get the resolution of the frames
frameResolution = size(currentFrame);

% Do the following for each video
for counterOfVideos = 1:numberOfVideos
  
  % Calculate number of frames in current video
  if counterOfVideos < numberOfVideos
    numberOfFramesInCurrentVideo = maxNumberOfFramesPerVideo;
  elseif counterOfVideos == numberOfVideos
    numberOfFramesInCurrentVideo = numberOfFramesInLastVideo;
  end    

  % Preallocate memory space for setOfAllFrames
  setOfFramesOfCurrentVideo = zeros([frameResolution(1) 371 numberOfFramesInCurrentVideo]);
  
  % Do the following for each of the frames in the current video
  for counterOfFramesInCurrentVideo = 1:numberOfFramesInCurrentVideo

    % Transfer all image data into the variable setOfFramesOfCurrentVideo
    currentFrameIndex = counterOfFramesInCurrentVideo + ((counterOfVideos-1)*maxNumberOfFramesPerVideo);
    currentFrame = load([pathFrames '/' char(listOfFrames(currentFrameIndex))], 'currentFrame');
    currentFrame = currentFrame.currentFrame;

    % Convert RGB image to grayscale
    currentFrame = rgb2gray(currentFrame);

    % Crop image
    currentFrame = currentFrame(:,130:500);

    % Transfer current image to setOfFramesOfCurrentVideo
    setOfFramesOfCurrentVideo(:,:,counterOfFramesInCurrentVideo) = currentFrame;

    % Save the frames as image files with JPEG compression
    currentFileName = char(listOfFrames(currentFrameIndex));
    imwrite(currentFrame, [pathFrames '/' currentFileName(1:(end-4)) '.jpg']);

    % Delete the image MAT-files
    delete([pathFrames '/' currentFileName])

    % Display percentageStatus
    percentageStatus = 100*(currentFrameIndex/totalNumberOfFrames);
    disp(['Preprocessing images (video ' num2str(counterOfVideos) '): ' num2str(percentageStatus) '%'])

  end
  clear counterOfFramesInCurrentVideo
  clear currentFrame
  clear currentFileName

  % Convert matrix to grayscale image
  % This is performed to be able to produce the video using writeVideo
  setOfFramesOfCurrentVideo = mat2gray(setOfFramesOfCurrentVideo);

  % Name the video and set the directory in which it will be saved
  currentVideoName = ['video' num2str(counterOfVideos) '.avi'];
  videoObj = VideoWriter([pathFrames '/' currentVideoName]); %#ok<TNMLP>

  % Define the number of frames per second
  videoObj.FrameRate = 50;

  % Create and save video
  open(videoObj);
  for counterOfFramesInCurrentVideo = 1:numberOfFramesInCurrentVideo
    
    % Calculate the current frame index from the total amouunt of frames
    currentFrameIndex = counterOfFramesInCurrentVideo + ((counterOfVideos-1)*maxNumberOfFramesPerVideo);

    % Write video
    writeVideo(videoObj,setOfFramesOfCurrentVideo(:,:,counterOfFramesInCurrentVideo));

    % Display percentageStatus
    percentageStatus = 100*currentFrameIndex/totalNumberOfFrames;
    disp(['Producing video ' num2str(counterOfVideos) ': ' num2str(percentageStatus) '%']);

  end
  close(videoObj);
  clear counterOfFramesInCurrentVideo
  clear percentageStatus
  clear currentVideoName
  clear numberOfFramesInCurrentVideo
  
end
clear numberOfFramesInLastVideo
clear counterOfVideos