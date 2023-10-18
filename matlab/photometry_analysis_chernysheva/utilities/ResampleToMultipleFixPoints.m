function [out] = ResampleToMultipleFixPoints(in, inFP, alignFP, out_pnts) 

%   function - ResampleToMultipleFixPoints
% 
%   Fritjof Helmchen, Brain Research Institute (Hifo), University of Zurich
%   last updated: 26.5.2016
%
%   resample a time series to warp the trace segementwise such that certain fix points
%   (FP) are aligned to a common reference (e.g. taking average segment durations (as fraction of total trial time) as alignFP)
%
%   out         - resampled trace (1..out_pnts)
%   in          - original trace (in_pnts elements)  
%   inFP        - vector with list of fix points of in vector  (e.g. start, cue,
%                turn, .., end)
%   alignFP     - reference fix points in normalized units (fractions), e.g.   (0, ..0.2, 0.5 ,..., 1)    
%   out_npnts   - total number of output points

in_pnts = numel(in); % number of fix points have to match; otherwise return

if numel(inFP)~=numel(alignFP) return; end
    
num_seg = numel(inFP)-1;
out = [];
seg_start = 1;
for i = 1:num_seg;
    if i==1; inseg_pnts = inFP(i+1)-inFP(i)+1; else inseg_pnts = inFP(i+1)-inFP(i); end     % input segment number of pnts; for i>1 take the next point
    if i < num_seg; seg_end = round(out_pnts*alignFP(i+1)); else seg_end = out_pnts; end  % define segment end point, assure that it ends with out_pnts 
    outseg_pnts = seg_end-seg_start+1; 
          
    xin = 1:1:inseg_pnts;
    dx = (inseg_pnts-1)/(outseg_pnts-1);   % time step after warp
    xout = dx:dx:dx*outseg_pnts;
    if i==1; tmpin = in(inFP(i):inFP(i+1)); else tmpin = in(inFP(i)+1:inFP(i+1)); end    % for i>1 take the next point
    
    tmpout = interp1(xin,tmpin,xout,'spline');                                           % interpolation of segment number i
    
    if i==1; out = tmpout; else out = horzcat(out,tmpout); end   % concatenate all segments one after the other
    
    seg_start = seg_end+1;
end

