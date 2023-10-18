%% wxcroproi
 % Crops ROI and gets ROI edges
 % JL Alatorre-Warren
 
function [roiCropped,...
          columnIndexMin,...
          columnIndexMax,...
          rowIndexMin,...
          rowIndexMax] = wxcroproi(roi)
         
% Convert from DIPimage to uint8 format
roiUint8 = uint8(roi);

% Find indices of nonzero elements
[rowIndices, columnIndices] = find(roiUint8);

% In DIPimage format, the index (0,0) is the origin instead of (1,1)
% In other words, DIPimage follows the C/C++ convention and not MATLAB
% convention; therefore, we substract 1 here to account for this.
rowIndexMin    = min(rowIndices) - 1;
rowIndexMax    = max(rowIndices) - 1;
columnIndexMin = min(columnIndices) - 1;
columnIndexMax = max(columnIndices) - 1;

% Get an additional set of rows and columns surrounding the cropped image
rowIndexMin    = rowIndexMin - 1;
rowIndexMax    = rowIndexMax + 1;
columnIndexMin = columnIndexMin - 1; 
columnIndexMax = columnIndexMax + 1;

% In DIPimage format, MATLAB rows are columns and vice versa
% Therefore, here the columns and row are swapped
roiCropped = roi(columnIndexMin:columnIndexMax,rowIndexMin:rowIndexMax);         