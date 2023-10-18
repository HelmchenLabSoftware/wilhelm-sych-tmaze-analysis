%% wxdrawrois
 % This script allows the user to manually draw up to 6 ROIs (regions
 % of interest
 % JL Alatorre Warren

function [roi01, roi02, roi03, roi04, roi05, roi06] = wxdrawrois(numberOfRois, image)

% Prevent errors from the user
if isscalar(numberOfRois) == 0
  numberOfRois = 1;
end

% Set the maximum number of ROIs to 6
if numberOfRois > 6
  numberOfRois = 6;
end

% Draw the ROIs
for i = 1:numberOfRois
  roi = roipoly(image);
  switch(i)
    case 1
      roi01 = roi;
    case 2
      roi02 = roi;
    case 3
      roi03 = roi;
    case 4
      roi04 = roi;
    case 5
      roi05 = roi;
    case 6
      roi06 = roi;      
  end
end

% Set nonexisting ROIs as roi01
switch(numberOfRois)
  case 1
    roi02 = roi01;
    roi03 = roi01;
    roi04 = roi01;
    roi05 = roi01;
    roi06 = roi01;    
  case 2
    roi03 = roi01;
    roi04 = roi01;
    roi05 = roi01;
    roi06 = roi01;    
    case 3
    roi04 = roi01;
    roi05 = roi01;
    roi06 = roi01;    
  case 4
    roi05 = roi01;
    roi06 = roi01;
  case 5
    roi06 = roi01;        
end

% DIPimage format
roi01 = dip_image(roi01);
roi02 = dip_image(roi02);
roi03 = dip_image(roi03);
roi04 = dip_image(roi04);
roi05 = dip_image(roi05);
roi06 = dip_image(roi06);