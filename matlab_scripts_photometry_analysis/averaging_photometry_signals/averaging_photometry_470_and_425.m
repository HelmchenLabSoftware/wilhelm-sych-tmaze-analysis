%% Averaging 488 and 425 together
Sessions2Average{1}=New_lPr0911_m007;%
Sessions2Average{2}=New_lPr1010_m008;%%
Sessions2Average{3}=New_lPr1010_m009;%%
Sessions2Average{4}=New_lPr1011_m011;%

%% 
Trial_LWhole_Correct_Average=Sessions2Average{1}.Trial_LWhole_Correct;
Trial_RWhole_Correct_Average=Sessions2Average{1}.Trial_RWhole_Correct;
Trial_LWhole_Mistake_Average=Sessions2Average{1}.Trial_LWhole_Mistake;
Trial_RWhole_Mistake_Average=Sessions2Average{1}.Trial_RWhole_Mistake;

for i=2:size(Sessions2Average,2)
 
Trial_LWhole_Correct_Average=vertcat(Trial_LWhole_Correct_Average,Sessions2Average{1,i}.Trial_LWhole_Correct);
Trial_RWhole_Correct_Average=vertcat(Trial_RWhole_Correct_Average,Sessions2Average{1,i}.Trial_RWhole_Correct);

end
%% 
Trial_LWhole_Correct_Average_reduced=Trial_LWhole_Correct_Average;
Trial_LWhole_Correct_Average_reduced(:,[4,7,11,13,16])=[];

Trial_RWhole_Correct_Average_reduced=Trial_RWhole_Correct_Average;
Trial_RWhole_Correct_Average_reduced(:,[4,7,11,13,16])=[];


Trial_LWhole_Correct_Average_3Phases_reduced=Trial_LWhole_Correct_Average_reduced(:,(1:11));
Trial_RWhole_Correct_Average_3Phases_reduced=Trial_RWhole_Correct_Average_reduced(:,(1:11));


Trial_LWhole_Correct_Average_3Phases=Trial_LWhole_Correct_Average(:,(1:15));
Trial_RWhole_Correct_Average_3Phases=Trial_RWhole_Correct_Average(:,(1:15));

%% Correct trials

sm_dT_RL_Whole=[]; % vector of mean dt between "event" time points
sd_dT_RL_Whole=[];
dT_RL_Whole=[];% matrix, colums: delta T between neighboring timepoints, raws: different trials
dt=[]; % vector std of dt between "event" time points
% Trial_LWhole=Trial_LWhole_Correct_Average_3Phases;
% Trial_RWhole=Trial_LWhole_Correct_Average_3Phases;

Trial_LWhole=Trial_LWhole_Correct_Average_3Phases_reduced;
Trial_RWhole=Trial_RWhole_Correct_Average_3Phases_reduced;

i=1;
 for i=1:size(Trial_LWhole,2)-1
[dt,sm_dt,sd_dt]=NewDeltaT(vertcat(Trial_LWhole(:,i),Trial_RWhole(:,i)),vertcat(Trial_LWhole(:,i+1),Trial_RWhole(:,i+1))); %counting dT between i and i+1 time point
dT_RL_Whole(i,:)=dt;
sm_dT_RL_Whole(i)=sm_dt; % median time intervals between events
sd_dT_RL_Whole(i)=sd_dt;
 end 
 %% 
 baseline_end_event_number = 6;
 for i=1:size(Sessions2Average,2) 
     TrialR=Sessions2Average{i}.Trial_RWhole_Correct;
     TrialR(:,[4,7,11,13,16])=[];
     TrialR_reduced=TrialR(:,1:11);
     TrialL=Sessions2Average{i}.Trial_LWhole_Correct;
     TrialL(:,[4,7,11,13,16])=[];
     TrialL_reduced=TrialL(:,1:11);
    
     
[SessionResR_Wh_488{i},dTR_resample]=ResampledSession_v2(Sessions2Average{i}.demodSig488,TrialR_reduced,sm_dT_RL_Whole, baseline_end_event_number);
[SessionResL_Wh_488{i},dTL_resample]=ResampledSession_v2(Sessions2Average{i}.demodSig488,TrialL_reduced,sm_dT_RL_Whole, baseline_end_event_number);
 SessionRes{i}=vertcat(SessionResR_Wh_488{i},SessionResL_Wh_488{i});
 end
 
 TrialsCorrect=SessionRes{1};
 ii=1;
 for ii=2:size(Sessions2Average,2)
     TrialsCorrect=vertcat(TrialsCorrect,SessionRes{ii});
 end
 
 SessionResR_Wh_all_GCaMP_488=SessionResR_Wh_488{1};
 SessionResL_Wh_all_GCaMP_488=SessionResL_Wh_488{1};
 for i=2:size(Sessions2Average,2)
     SessionResR_Wh_all_GCaMP_488=vertcat(SessionResR_Wh_all_GCaMP_488,SessionResR_Wh_488{i});
     SessionResL_Wh_all_GCaMP_488=vertcat(SessionResL_Wh_all_GCaMP_488,SessionResL_Wh_488{i});
 end
 %% 
  baseline_end_event_number = 6;
 for i=1:size(Sessions2Average,2) 
     TrialR=Sessions2Average{i}.Trial_RWhole_Correct;
     TrialR(:,[4,7,11,13,16])=[];
     TrialR_reduced=TrialR(:,1:11);
     TrialL=Sessions2Average{i}.Trial_LWhole_Correct;
     TrialL(:,[4,7,11,13,16])=[];
     TrialL_reduced=TrialL(:,1:11);
    
     
[SessionResR_Wh_425{i},dTR_resample]=ResampledSession_v2(Sessions2Average{i}.demodSig425,TrialR_reduced,sm_dT_RL_Whole, baseline_end_event_number);
[SessionResL_Wh_425{i},dTL_resample]=ResampledSession_v2(Sessions2Average{i}.demodSig425,TrialL_reduced,sm_dT_RL_Whole, baseline_end_event_number);
 SessionRes{i}=vertcat(SessionResR_Wh_425{i},SessionResL_Wh_425{i});
 end
 
 TrialsCorrect=SessionRes{1};
 ii=1;
 for ii=2:size(Sessions2Average,2)
     TrialsCorrect=vertcat(TrialsCorrect,SessionRes{ii});
 end
 
 SessionResR_Wh_all_GCaMP_425=SessionResR_Wh_425{1};
 SessionResL_Wh_all_GCaMP_425=SessionResL_Wh_425{1};
 for i=2:size(Sessions2Average,2)
     SessionResR_Wh_all_GCaMP_425=vertcat(SessionResR_Wh_all_GCaMP_425,SessionResR_Wh_425{i});
     SessionResL_Wh_all_GCaMP_425=vertcat(SessionResL_Wh_all_GCaMP_425,SessionResL_Wh_425{i});
 end
 %% 
 NamePlot='GCaMP 488 vs 425 control';

color_input2=[128,128,128]/255; %grey
color_input1='k';
transparent=0;
y_axis1=-2;
y_axis2=2.5;


figure

SessionRes_GCaMP_488=vertcat(SessionResL_Wh_all_GCaMP_488,SessionResR_Wh_all_GCaMP_488);
MeanTrial_GCaMP_488=mean(SessionRes_GCaMP_488,1); 
sem_MeanTrial_GCaMP_488=std(SessionRes_GCaMP_488,1)/sqrt(size(SessionRes_GCaMP_488,1));% std(x)/sqrt(length(x));
dt_resample=0.0501;

t= -1:dt_resample:dt_resample*(length(MeanTrial_GCaMP_488)-1)-1;
x_axis1=-1;
x_axis2=max(t);%Time4Plot;
H=NewShadedErrorBar(t,smooth(MeanTrial_GCaMP_488),smooth(sem_MeanTrial_GCaMP_488),'k',transparent,21,y_axis1,y_axis2,-1,max(t),'','',size(SessionRes_GCaMP_488,1),color_input1) ;

hold on
SessionRes_GCaMP_425=vertcat(SessionResL_Wh_all_GCaMP_425,SessionResR_Wh_all_GCaMP_425);

MeanTrial_GCaMP_425=mean(SessionRes_GCaMP_425,1); 
sem_MeanTrial_GCaMP_425=std(SessionRes_GCaMP_425,1)/sqrt(size(SessionRes_GCaMP_425,1));% std(x)/sqrt(length(x));
dt_resample=dTR_resample; %(1+sum(sm_dT_RL_Whole))/(size(SessionRes,2)+length(sm_dT_RL_Whole)+1);
dt_resample=0.0501;%0018234865063;
t= -1:dt_resample:dt_resample*(length(MeanTrial_GCaMP_425)-1)-1;
x_axis1=-1;
x_axis2=max(t);%Time4Plot;
H=NewShadedErrorBar(t,smooth(MeanTrial_GCaMP_425)-1.5,smooth(sem_MeanTrial_GCaMP_425),'g',transparent,21,y_axis1,y_axis2,-1,max(t),'','',size(SessionRes_GCaMP_425,1),color_input1) ;
hold on


for qq=1:length(sm_dT_RL_Whole)
line([sum(sm_dT_RL_Whole(1:qq))-dt_resample*(qq+1) sum(sm_dT_RL_Whole(1:qq))-dt_resample*(qq+1)],[ y_axis1 y_axis2],'Color','k','LineWidth',1,'LineStyle','--')
end
i=1;

ylabel ''

axis([x_axis1 x_axis2 y_axis1 y_axis2])

hf=gcf;

axis('square')

title ({sprintf('%s',NamePlot);''})
 %% 
  excel_source_datafile = horzcat (t',...
                                   MeanTrial_GCaMP_488',...
                                   sem_MeanTrial_GCaMP_488',...
                                   MeanTrial_GCaMP_425',...
                                   sem_MeanTrial_GCaMP_425',...
                                   NaN(length(t),1),...
                                   SessionRes_GCaMP_488',...
                                   NaN(length(t),1),...
                                   SessionRes_GCaMP_425');

 %% calculate mean over period for each session
%Encoding:20: sum(sm_dT_RL_Whole(1:qq))-dt_resample*(qq+1)%qq=5
%Maintenance:sum(sm_dT_RL_Whole(1:qq1))-dt_resample*(qq1+1):sum(sm_dT_RL_Whole(1:qq2))-dt_resample*(qq2+1)
%Retrieval: sum(sm_dT_RL_Whole(1:qq3))-dt_resample*(qq3+1):end qq3=8
dt_resample=dTR_resample;
T1=1; %1 s before the start
tt=5; % 6 events belong to Encoding
tt3=8; %Retrieval starts after event 9
T2=sum(sm_dT_RL_Whole(1:tt))-dt_resample*(tt+1);
T3=sum(sm_dT_RL_Whole(1:tt3))-dt_resample*(tt3+1);
Fs=20;


SessionResL_Wh=SessionResL_Wh_488;
SessionResR_Wh=SessionResR_Wh_488;
 %by ----Trial----
 for i=1:size(Sessions2Average,2) 
     EncodingResL{i}=SessionResL_Wh{1,i}(:, Fs:round (Fs*T2)+Fs);
     MaintenanceResL{i}=SessionResL_Wh{1,i}(:, round (Fs*T2)+Fs:round (Fs*T3)+Fs);
     RetrievalResL{i}=SessionResL_Wh{1,i}(:, round (Fs*T3)+Fs:end);
     
     EncodingResR{i}=SessionResR_Wh{1,i}(:, Fs:round (Fs*T2)+Fs);
     MaintenanceResR{i}=SessionResR_Wh{1,i}(:, round (Fs*T2)+Fs:round (Fs*T3)+Fs);
     RetrievalResR{i}=SessionResR_Wh{1,i}(:, round (Fs*T3)+Fs:end);
     
%      EncodingResL_mean(i)=mean((EncodingResL{i},1));
%      MaintenanceResL_mean(i)=mean((MaintenanceResL{i},1));
%      RetrievalResL_mean(i)=mean((RetrievalResL{i},1));
%      
%      EncodingResR_mean(i)=mean(mean(EncodingResR{i},1));
%      MaintenanceResR_mean(i)=mean(mean(MaintenanceResR{i},1));
%      RetrievalResR_mean(i)=mean(mean(RetrievalResR{i},1));
%      
%      EncodingRes_meanAll(i)=mean(mean(vertcat(EncodingResR{i},EncodingResL{i}),1));
%      MaintenanceRes_meanAll(i)=mean(mean(vertcat(MaintenanceResR{i},MaintenanceResL{i}),1));
%      RetrievalRes_meanAll(i)=mean(mean(vertcat(RetrievalResR{i},RetrievalResL{i}),1));

 end
     EncodingResL_mean_Trials=[];
     MaintenanceResL_mean_Trials=[];
     RetrievalResL_mean_Trials=[];
for ii=1:size(Sessions2Average,2) 
     
    for kk=1:size(EncodingResL{ii},1)
     EncodingResL_mean_Trial(kk)=mean((EncodingResL{ii}(kk,1:end)));
     MaintenanceResL_mean_Trial(kk)=mean((MaintenanceResL{ii}(kk,1:end)));
     RetrievalResL_mean_Trial(kk)=mean((RetrievalResL{ii}(kk,1:end)));  
    end
    
     EncodingResL_mean_Trials=[EncodingResL_mean_Trials,EncodingResL_mean_Trial];
     MaintenanceResL_mean_Trials=[MaintenanceResL_mean_Trials,MaintenanceResL_mean_Trial];
     RetrievalResL_mean_Trials=[RetrievalResL_mean_Trials,RetrievalResL_mean_Trial];
      
clearvars EncodingResL_mean_Trial MaintenanceResL_mean_Trial RetrievalResL_mean_Trial   
end
     EncodingResR_mean_Trials=[];
     MaintenanceResR_mean_Trials=[];
     RetrievalResR_mean_Trials=[];
for ii=1:size(Sessions2Average,2) 
     
    for kk=1:size(EncodingResR{ii},1)
     EncodingResR_mean_Trial(kk)=mean((EncodingResR{ii}(kk,1:end)));
     MaintenanceResR_mean_Trial(kk)=mean((MaintenanceResR{ii}(kk,1:end)));
     RetrievalResR_mean_Trial(kk)=mean((RetrievalResR{ii}(kk,1:end)));  
    end
    
     EncodingResR_mean_Trials=[EncodingResR_mean_Trials,EncodingResR_mean_Trial];
     MaintenanceResR_mean_Trials=[MaintenanceResR_mean_Trials,MaintenanceResR_mean_Trial];
     RetrievalResR_mean_Trials=[RetrievalResR_mean_Trials,RetrievalResR_mean_Trial];
      
clearvars EncodingResR_mean_Trial MaintenanceResR_mean_Trial RetrievalResR_mean_Trial   
end
%  EncodingRes_mean_Trials_RL_425=[EncodingResL_mean_Trials,EncodingResR_mean_Trials];
%  MaintenanceRes_mean_Trials_RL_425=[MaintenanceResL_mean_Trials, MaintenanceResR_mean_Trials];
%  RetrievalRes_mean_Trials_RL_425=[RetrievalResL_mean_Trials,RetrievalResR_mean_Trials];
 
%    
 EncodingRes_mean_Trials_RL_488=[EncodingResL_mean_Trials,EncodingResR_mean_Trials];
 MaintenanceRes_mean_Trials_RL_488=[MaintenanceResL_mean_Trials, MaintenanceResR_mean_Trials];
  RetrievalRes_mean_Trials_RL_488=[RetrievalResL_mean_Trials,RetrievalResR_mean_Trials];
  %% 
  MaintenanceRes_mean_Trials_RL_488 = MaintenanceRes_mean_Trials_RL_488';
  MaintenanceRes_mean_Trials_RL_425 = MaintenanceRes_mean_Trials_RL_425';
  [p1,h] = signrank(MaintenanceRes_mean_Trials_RL_488, MaintenanceRes_mean_Trials_RL_425);

  