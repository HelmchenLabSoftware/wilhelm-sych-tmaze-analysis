%% wxnotifyservoerror
 % Warns the user that the servomotor could not open and requests immediate action
 % JL Alatorre Warren
 
% Create warning image
warningImage = newim(1000,1000);
warningImage(:,:) = 1;
warningImage = joinchannels('RGB', warningImage*255, warningImage*0, warningImage*0);

% Display warning image
dipfig warningImage
dipshow(warningImage)
figureWarningImage = gcf;
figureWarningImage.Position = [1289 49 853 635];