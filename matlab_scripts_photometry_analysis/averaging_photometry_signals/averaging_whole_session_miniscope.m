%This function allows to average miniscopr signals collected during T-maze
%alternation task


%Data should be pre-processed (e.g. ceells activity and trial structure extracted and sorted
%for each mouse). Since data for miniscope is recorded and stored
%differently, code is slightly adapted 

% Author: Maria Wilhelm

%Clean the previous matlab sesion

close all
clearvars
clc

folder_with_matlab_functions = 'J:\H50_Disk_D\WM_materials_for_the_paper\Matlab_Scripts_final'; %give the folder with the matlab functions:
addpath(folder_with_matlab_functions)

folder_with_data_for_averaging = 'J:\H50_Disk_D\WM_materials_for_the_paper\Data_for_Figures\Data_Sup_Figure8\sorted_data\behavior'; %give the folder with the raw data
cd (folder_with_data_for_averaging)


files= dir('*.mat');
[mice_to_pool, ~]=size(files);

for count_mice= 1:mice_to_pool
    
    current_file_name = files(count_mice,:).name;
    current_mouseID = files(count_mice,:).name(1:4);
    % load data by mouse
    load(current_file_name);
    
end
%% load all resorded behavior sessions

session_number = 1;
for count_mice = 1:mice_to_pool
    
    current_file_name = files(count_mice,:).name(1:end-4);
    current_mouse_data = eval(current_file_name);
    number_of_sessions(count_mice) = size(fieldnames(current_mouse_data),1);
    
    for current_session = 1:number_of_sessions(count_mice)
        
        session_names = fieldnames (current_mouse_data);
        current_session_name = session_names {current_session};
        
        sessions_to_average {session_number} = current_mouse_data.(current_session_name);
        session_number = session_number + 1;
    end
    
end
%% load miniscope recording data
folder_with_data_for_averaging = 'J:\H50_Disk_D\WM_materials_for_the_paper\Data_for_Figures\Data_Sup_Figure8\sorted_data\miniscope_recording'; %give the folder with the raw data
cd (folder_with_data_for_averaging)


files_miniscope= dir('*.mat');
[mice_to_pool, ~]=size(files_miniscope);

for count_mice= 1:mice_to_pool
    
    current_file_name = files_miniscope(count_mice,:).name;
    current_mouseID = files_miniscope(count_mice,:).name(1:4);
    % load data by mouse
    load(current_file_name);
end

%

session_number = 1;
for count_mice = 1:mice_to_pool
    
    current_file_name = files_miniscope(count_mice,:).name(1:end-4);
    current_miniscope_data = eval(current_file_name);
    number_of_sessions(count_mice) = size(fieldnames(current_miniscope_data),1);
    
    for current_session = 1:number_of_sessions(count_mice)
        
        session_names = fieldnames (current_miniscope_data);
        current_session_name = session_names {current_session};
        
        sessions_signal_to_average{session_number} = sum(current_miniscope_data.(current_session_name),2);
        session_number = session_number + 1;
    end
    
end

%% Averaging all sessions
N_TTL = sessions_to_average{1}.allSignals(:,1)>2;
FirstTTL = find(N_TTL,1);
index2delete = round(FirstTTL/100);
Time2Allign = index2delete*0.05;


Trial_LWhole_Correct_Average = sessions_to_average{1}.Trial_LWhole_Correct;
Trial_RWhole_Correct_Average = sessions_to_average{1}.Trial_RWhole_Correct;

Trial_LWhole_Correct_alligned = Trial_LWhole_Correct_Average - Time2Allign;
Trial_RWhole_Correct_alligned = Trial_RWhole_Correct_Average - Time2Allign;

for current_session = 2:size(sessions_to_average,2)

 current_trial_RWhole_Correct=[];
 current_trial_LWhole_Correct=[];
 
current_Trial_LWhole_Correct_Average = sessions_to_average{1,current_session}.Trial_LWhole_Correct;
current_Trial_RWhole_Correct_Average = sessions_to_average{1,current_session}.Trial_RWhole_Correct;

N_TTL = sessions_to_average{current_session}.allSignals(:,1)>2;
FirstTTL = find(N_TTL,1);
index2delete = round(FirstTTL/100);
Time2Allign = index2delete*0.05;

current_trial_LWhole_Correct = current_Trial_LWhole_Correct_Average-Time2Allign;
current_trial_RWhole_Correct = current_Trial_RWhole_Correct_Average-Time2Allign;

Trial_LWhole_Correct_alligned = vertcat(Trial_LWhole_Correct_alligned,current_trial_LWhole_Correct);
Trial_RWhole_Correct_alligned = vertcat(Trial_RWhole_Correct_alligned,current_trial_RWhole_Correct);


Trial_LWhole_Correct_Average = vertcat(Trial_LWhole_Correct_Average,sessions_to_average{1,current_session}.Trial_LWhole_Correct);
Trial_RWhole_Correct_Average = vertcat(Trial_RWhole_Correct_Average,sessions_to_average{1,current_session}.Trial_RWhole_Correct);

end

% 
Trial_LWhole_Correct_Average_downsampled = Trial_LWhole_Correct_Average;
Trial_LWhole_Correct_Average_downsampled(:,[4,7,11,13,16]) = [];

Trial_RWhole_Correct_Average_downsampled = Trial_RWhole_Correct_Average;
Trial_RWhole_Correct_Average_downsampled(:,[4,7,11,13,16]) = [];

Trial_LWhole_Correct_Average_3Phases_downsampled = Trial_LWhole_Correct_Average_downsampled(:,(1:11));
Trial_RWhole_Correct_Average_3Phases_downsampled = Trial_RWhole_Correct_Average_downsampled(:,(1:11));

sm_dT_RL_Whole=[]; % vector of mean dt between "event" time points
sd_dT_RL_Whole=[];
dT_RL_Whole=[];% matrix, colums: delta T between neighboring timepoints, raws: different trials
dt=[]; % vector std of dt between "event" time points


Trial_LWhole=Trial_LWhole_Correct_Average_3Phases_downsampled;
Trial_RWhole=Trial_RWhole_Correct_Average_3Phases_downsampled;

current_session=1;
 for current_session=1:size(Trial_LWhole,2)-1
[dt,sm_dt,sd_dt]=NewDeltaT(vertcat(Trial_LWhole(:,current_session),Trial_RWhole(:,current_session)),vertcat(Trial_LWhole(:,current_session+1),Trial_RWhole(:,current_session+1))); %counting dT between i and i+1 time point
dT_RL_Whole(current_session,:)=dt;
sm_dT_RL_Whole(current_session)=sm_dt;
sd_dT_RL_Whole(current_session)=sd_dt;
 end 
 %% resample all trials to the mean delta time between events
baseline_end_event_number = 6; %event after the t-junstion
 for current_session = 1:size(sessions_to_average,2) 
     
Trial_RWhole_Correct_alligned = [];
Trial_LWhole_Correct_alligned = [];
current_Trial_LWhole_Correct_Average = [];
current_Trial_RWhole_Correct_Average = [];

     
current_Trial_LWhole_Correct_Average = sessions_to_average{1,current_session}.Trial_LWhole_Correct;

current_Trial_LWhole_Correct_Average(:,[4,7,11,13,16]) = [];
TrialL_reduced = current_Trial_LWhole_Correct_Average(:,1:11);
     
current_Trial_RWhole_Correct_Average = sessions_to_average{1,current_session}.Trial_RWhole_Correct;
current_Trial_RWhole_Correct_Average(:,[4,7,11,13,16]) = [];
TrialR_reduced = current_Trial_RWhole_Correct_Average(:,1:11);


N_TTL = sessions_to_average{current_session}.allSignals(:,1)>2;
FirstTTL = find(N_TTL,1);
index2delete = round(FirstTTL/100);
Time2Allign = index2delete*0.05;

Trial_LWhole_Correct_alligned = TrialL_reduced - Time2Allign;
Trial_RWhole_Correct_alligned = TrialR_reduced - Time2Allign;
[SessionResR_Wh{current_session},dTR_resample] = ResampledSession_v2D_average(sessions_signal_to_average{current_session},Trial_RWhole_Correct_alligned,sm_dT_RL_Whole,baseline_end_event_number);
[SessionResL_Wh{current_session},dTL_resample] = ResampledSession_v2D_average(sessions_signal_to_average{current_session},Trial_LWhole_Correct_alligned,sm_dT_RL_Whole,baseline_end_event_number);

 end
%  
 SessionResR_Wh_all_Miniscope=SessionResR_Wh{1};
 SessionResL_Wh_all_Miniscope=SessionResL_Wh{1};
 
 for current_session=2:size(sessions_to_average,2)
     
     SessionResR_Wh_all_Miniscope=vertcat(SessionResR_Wh_all_Miniscope,SessionResR_Wh{current_session});
     SessionResL_Wh_all_Miniscope=vertcat(SessionResL_Wh_all_Miniscope,SessionResL_Wh{current_session});
     
 end
 %% plot average data
name_plot = 'Miniscope SumCells all mice RandL deconvolved';
color_input1=[0,128,0]/255;
transparent=0;
y_axis1=-2;
y_axis2=4;


NamesRWh={''};
NamesL={''};

hF=figure;
hold on
Session_to_plot=vertcat(SessionResL_Wh_all_Miniscope,SessionResR_Wh_all_Miniscope);
MeanTrial=mean(Session_to_plot,1); 
sem_MeanTrial=std(Session_to_plot,1)/sqrt(size(Session_to_plot,1));% std(x)/sqrt(length(x));

dt_resample = (1+sum(sm_dT_RL_Whole(1:end)))/length(MeanTrial); dTL_resample;
t= -1:dt_resample:dt_resample*(length(MeanTrial)-1)-1;
x_axis1=-1;
x_axis2=max(t);%Time4Plot;
H=NewShadedErrorBar(t,movmean(MeanTrial,5),movmean(sem_MeanTrial,5),'k',transparent,21,y_axis1,y_axis2,-1,max(t),'','',size(Session_to_plot,1),'k') ;
hold on


for qq=1:length(sm_dT_RL_Whole)
line([sum(sm_dT_RL_Whole(1:qq))-dt_resample*(qq+1) sum(sm_dT_RL_Whole(1:qq))-dt_resample*(qq+1)],[ y_axis1 y_axis2],'Color','k','LineWidth',1,'LineStyle','--')
end

hold on
x_axis1=-1;
x_axis2=sum(sm_dT_RL_Whole)-1;
% set(hF, 'Position', get(0, 'Screensize')) 
axis([x_axis1 x_axis2 y_axis1 y_axis2])
h=gcf;
set(h,'PaperOrientation','landscape');
h_title=title({name_plot;'';'';'';'';'';''});
axis('square')
set(gca,'Position',[0.1300    0.1100    0.7750    0.5455])
save_plots(name_plot)
 

%% collect source data to excel files

    individual_trials_z_score = Session_to_plot;
    average_z_score = movmean(MeanTrial,5);
    sem = movmean(sem_MeanTrial,5);
    
    
    excel_source_data = vertcat (t,...
                                 movmean(MeanTrial,5),...
                                 sem,...
                                 individual_trials_z_score);


save('variables_17082023','individual_trials_z_score','excel_source_data')

%% calculate mean values for each period and plot (SFig8b)
time_interval_encoding = [1,sum(sm_dT_RL_Whole(1:4))];
index_encoding_start = round(time_interval_encoding(1)/dt_resample);
index_encoding_end = round(time_interval_encoding(2)/dt_resample);

time_interval_maintenance = [sum(sm_dT_RL_Whole(1:4)),sum(sm_dT_RL_Whole(1:6))];   
index_maintenance_start = round(time_interval_maintenance(1)/dt_resample);
index_maintenance_end = round(time_interval_maintenance(2)/dt_resample);

time_interval_retrieval = [sum(sm_dT_RL_Whole(1:6)),sum(sm_dT_RL_Whole(1:10))];   
index_retrieval_start = round(time_interval_retrieval(1)/dt_resample);
index_retrieval_end = round(time_interval_retrieval(2)/dt_resample);

mean_per_session_resampled = [];

for current_session=1:size(sessions_to_average,2)
    
mean_signal_per_session = mean(vertcat(SessionResL_Wh{current_session},SessionResR_Wh{current_session}),1);

    
mean_per_session_resampled.encoding(current_session) = mean(mean_signal_per_session(index_encoding_start:index_encoding_end));
mean_per_session_resampled.maintenance(current_session) = mean(mean_signal_per_session(index_maintenance_start:index_maintenance_end));
mean_per_session_resampled.retrieval(current_session) = mean(mean_signal_per_session(index_retrieval_start:index_retrieval_end));

end

figure

input_data = [{mean_per_session_resampled.encoding},...
              {mean_per_session_resampled.maintenance},...
              {mean_per_session_resampled.retrieval}];
          
color_cohort = [{[0.59,0.53,0.82]},{[0.87,0.65,0.37]},{[0.53,0.82,0.82]}];
legend_cohort = [{'Enc'},{'M'},{'Ret'}];
plot_title = 'Mean per Period';
font_size = 20;
paired = 1;

SEM_Input = bar_plots_cohorts(input_data,legend_cohort,color_cohort, plot_title,'z-score',font_size, paired);

[p,t,stats] = anova1([input_data{1}', input_data{2}', input_data{3}']);
multicompare_results = multcompare(stats);


excel_source_data_SFig8_panel_b = vertcat ([input_data{1}', input_data{2}', input_data{3}'],...
                                           [mean(input_data{1}'), mean(input_data{2}'), mean(input_data{3}')],...
                                           [SEM_Input(1), SEM_Input(2), SEM_Input(3)],...
                                           [multicompare_results(1,6), multicompare_results(2,6), multicompare_results(3,6)]);
