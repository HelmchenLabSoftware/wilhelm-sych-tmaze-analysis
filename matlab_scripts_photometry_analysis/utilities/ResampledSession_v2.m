function [SessionRes,dt_resample]=ResampledSession_v2(demodSig,Session,sm_dT, baseline_end_event_number)
%makes a matrix SessionRes composed from resampled trials, then makes a
%plot and saves it
% addpath('P:');
% dt=[];
% sm_dt=[];
% sd_dt=[];
% dT=[];
% i=1;
% % alignFP    - reference fix points in normalized units (fractions)

alignFP=[]; 
% Time_trial=1+sum(sm_dT)+1;
Time_trial=1+sum(sm_dT);
alignFP=[0,1/Time_trial];
i=1;
for i=1:length(sm_dT)
alignFP=horzcat(alignFP,alignFP(end)+(sm_dT(i)./Time_trial));
end 
% alignFP=horzcat(alignFP,1);
out_npnts=round(20*Time_trial);
dt_resample=50e-3;%Time_trial/(out_npnts-1);
% 
inFP=Session/50e-3; %number of point in demodSig vector

% input:
inFP_1=[]; %make the fix ponts start from "1" to apply it later in ResampleToMultipleFixPoints to the vectors of single trials
inFP_2=[];


for j=1:size(inFP,1)
inFP_1(j,:)=round(inFP(j,:)-inFP(j,1)+1/dt_resample);%approximately 20 points is 1 second, but it can be 19.928
inFP_1(j,1)=21;
% inFP_2(j,:)=[1,inFP_1(j,:),round(inFP_1(j,end)+1/dt_resample)];%add 1second before and after vector of points
inFP_2(j,:)=[1,inFP_1(j,:)];%add 1 second before
end

% deleete  trials when events are too close to each other, so no interpolation is possible
j=1;
TrialToDelete=[];
for k=1:size(inFP_2,1)
   for i=2:size(inFP_2,2)
      if inFP_2(k,i)-inFP_2(k,i-1)<2
          TrialToDelete(j)=k;
      else
      end
   end
   j=j+1;
end
 
TrialToDelete(TrialToDelete==0)=[];
inFP_2(TrialToDelete,:)=[];
inFP(TrialToDelete,:)=[];
% 

% analize the single trial and save only resampled OUT
SessionRes=[];
for i=1:1:size(inFP,1)
    
range1=inFP_2(i,baseline_end_event_number)-20;
range2=inFP_2(i,baseline_end_event_number);


F_in=zscore_dF_norm(demodSig(round(inFP(i,1)-20):round(inFP(i,end))),range1,range2);


[SessionRes(i,:)]=ResampleToMultipleFixPoints(F_in, inFP_2(i,:), alignFP, out_npnts);
end

%delete the points on the border of the resampled periods (it is mostly
%giving a "jump" in the signal
 PointsDel=round(alignFP*out_npnts);
 PointsDel(1)=[];
SessionRes(:,PointsDel)=[];

alignFT=alignFP.*dt_resample;

end
