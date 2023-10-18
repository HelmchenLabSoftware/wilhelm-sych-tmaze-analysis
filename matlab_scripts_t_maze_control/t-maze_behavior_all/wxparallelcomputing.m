%% wxparallelcomputing
 % Root script that allows parallel computing of different routines
 % JL Alatorre Warren

%% Test Lumenera camera

% Assign camera number
% Display video preview
% cameraNumber = 1; 
% LucamShowPreview(cameraNumber);

dipstart;

%% Parallel Computing 

delete(gcp('nocreate'))
 
sInfoParpool = parpool(2);
parfor ii = 1:3
    if ii == 1
      wxfunction1(1,1)
    elseif ii == 2
      wxfunction2(1,1)
    elseif ii == 3
      tic
      disp('starts')
      cameraNumber = 1; 
      display(cameraNumber)
      LucamShowPreview(cameraNumber);
      LucamSetFrameRate(20, 1)
      toc
      tic
      frames = LucamTakeVideo(1000, 1);
      toc
      frames = dip_image(frames);
      frames = frames(:,:,1,:);
      framesSqueeze = squeeze(frames);
      framesSqueezeCrop = framesSqueeze(300:860,0:765,:);
      dipshow(framesSqueezeCrop)
      toc
    elseif ii == 4
      frames = wxsavevideo;
    elseif ii == 5
      for xx = 1:1000
        tic
        cameraNumber = 1; 
        LucamShowPreview(cameraNumber);
        frameMonochrome = LucamCaptureMonochromeFrame(cameraNumber);
        figure('Name','MonochromeColors');
        imagesc(frameMonochrome, [0,255]);
        frameMonochrome = dip_image(frameMonochrome)
        display(xx)
        toc
      end
    elseif ii == 6
      for jj = 1:1000
        pause(2)
        display(jj)
      end
    elseif ii == 7
      frames = LucamTakeVideo(2, 1);
      frames = dip_image(frames)
      frames = frames(:,:,1,:);
      dipshow(frames)
    end
end