%Clean the previous matlab sesion

close all
clearvars
clc

folder_with_matlab_functions = 'J:\H50_Disk_D\WM_materials_for_the_paper\Matlab_Scripts'; %give the folder with the matlab functions:
addpath(folder_with_matlab_functions)

folder_with_all_data_for_averaging = ('J:\H50_Disk_D\WM_materials_for_the_paper\Data_Figure4\sorted_data\individual_mice_data\Cmpd_and_Veh'); %give the folder with the raw data
[sessions_to_average] = extract_session_pharma (folder_with_all_data_for_averaging); 

% Averaging times between all events in all sessions

sm_dT_RL_Whole = average_time_between_events (sessions_to_average);

clearvars -except sm_dT_RL_Whole

% 
folder_for_figures = 'J:\H50_Disk_D\WM_materials_for_the_paper\Data_Figure4\figures';
folder_with_data_for_averaging.Cmpd = ('J:\H50_Disk_D\WM_materials_for_the_paper\Data_Figure4\sorted_data\individual_mice_data\Cmpd'); %give the folder with the raw data
[sessions_to_average_Cmpd] = extract_session_pharma (folder_with_data_for_averaging.Cmpd);


folder_with_data_for_averaging.Veh = ('J:\H50_Disk_D\WM_materials_for_the_paper\Data_Figure4\sorted_data\individual_mice_data\Veh'); %give the folder with the raw data
[sessions_to_average_Veh] = extract_session_pharma (folder_with_data_for_averaging.Veh);

%% 

baseline_end_event_number = 6;
[SessionResR_Wh_all_photometry_Cmpd,...
 SessionResL_Wh_all_photometry_Cmpd,...
 SessionResR_Wh.Cmpd,...
SessionResL_Wh.Cmpd] = calculate_resampled_signals (sessions_to_average_Cmpd, sm_dT_RL_Whole, baseline_end_event_number);

[SessionResR_Wh_all_photometry_Veh,...
 SessionResL_Wh_all_photometry_Veh,...
 SessionResR_Wh.Veh,...
SessionResL_Wh.Veh] = calculate_resampled_signals (sessions_to_average_Veh, sm_dT_RL_Whole, baseline_end_event_number);

%% 
name_plot = 'Cmpd vs Veh all mice RandL';
color_input1=[0,128,0]/255;
transparent=0;
y_axis1=-2;
y_axis2=4;


NamesRWh={''};%{'Start'  'Water' 'Finish Lick' 'T-j Back' '' 'End' 'Start' '' '' '' 'Water' 'Finish Lick' 'T-j Back' '' 'End' 'End+10s'};
NamesL={''};

hF=figure;
hold on
Session_Cmpd_to_plot = vertcat(SessionResR_Wh_all_photometry_Cmpd,SessionResL_Wh_all_photometry_Cmpd);
MeanTrial_Cmpd = mean(Session_Cmpd_to_plot,1); 
sem_MeanTrial_Cmpd = std(Session_Cmpd_to_plot,1)/sqrt(size(Session_Cmpd_to_plot,1));% std(x)/sqrt(length(x));

dt_resample = (1+sum(sm_dT_RL_Whole(1:end)))/length(MeanTrial_Cmpd);
t= -1:dt_resample:dt_resample*(length(MeanTrial_Cmpd)-1)-1;
x_axis1=-1;
x_axis2=max(t);%Time4Plot;
H=NewShadedErrorBar(t,movmean(MeanTrial_Cmpd,5),movmean(sem_MeanTrial_Cmpd,5),'k',transparent,21,y_axis1,y_axis2,-1,max(t),'','',size(Session_Cmpd_to_plot,1),'g') ;
hold on


Session_Veh_to_plot=vertcat(SessionResR_Wh_all_photometry_Veh,SessionResL_Wh_all_photometry_Veh);
MeanTrial_Veh=mean(Session_Veh_to_plot,1); 
sem_MeanTrial_Veh=std(Session_Veh_to_plot,1)/sqrt(size(Session_Veh_to_plot,1));% std(x)/sqrt(length(x));

dt_resample = (1+sum(sm_dT_RL_Whole(1:end)))/length(MeanTrial_Veh);
t= -1:dt_resample:dt_resample*(length(MeanTrial_Veh)-1)-1;
x_axis1=-1;
x_axis2=max(t);%Time4Plot;
H=NewShadedErrorBar(t,movmean(MeanTrial_Veh,5),movmean(sem_MeanTrial_Veh,5),'k',transparent,21,y_axis1,y_axis2,-1,max(t),'','',size(Session_Veh_to_plot,1),'k') ;
hold on
for qq=1:length(sm_dT_RL_Whole)
line([sum(sm_dT_RL_Whole(1:qq))-dt_resample*(qq+1) sum(sm_dT_RL_Whole(1:qq))-dt_resample*(qq+1)],[ y_axis1 y_axis2],'Color','k','LineWidth',1,'LineStyle','--')
end


hold on
h=text(0,y_axis2,'Start','FontSize',30,'FontWeight','bold');
set(h,'rotation',90)
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
cd(folder_for_figures)
save_plots(name_plot)

%% collect source data to excel files

%SFig 5a
    individual_trials_z_score = Session_Cmpd_to_plot;
    average_z_score = movmean(MeanTrial_Cmpd,5);
    sem = movmean(sem_MeanTrial_Cmpd,5);
    
    
    excel_source_data_Cmpd = vertcat (t,...
                                 average_z_score,...
                                 sem,...
                                 individual_trials_z_score);
                             
    individual_trials_z_score = Session_Veh_to_plot;
    average_z_score = movmean(MeanTrial_Veh,5);
    sem = movmean(sem_MeanTrial_Veh,5);
    
    
    excel_source_data_Veh = vertcat (t,...
                                 average_z_score,...
                                 sem,...
                                 individual_trials_z_score);
                             %% 
%SFig 5c
cd ('J:\H50_Disk_D\WM_materials_for_the_paper\Data_Figure4')
MK801_Cmpd_behavior = xlsread('MK801_and_Cmpd');
trial_duration_sfig15c = MK801_Cmpd_behavior(:,11:14);
mean_trial_duration_sfig15c = mean(trial_duration_sfig15c,1);
sem = std(trial_duration_sfig15c,1)/sqrt(size(trial_duration_sfig15c,1));


excel_source_data_sfig15c = vertcat (trial_duration_sfig15c,mean_trial_duration_sfig15c, sem);

[p,t,stats] = anova1(trial_duration_sfig15c);
multicompare_results = multcompare(stats);

%% 
trial_duration_sfig15d = MK801_Cmpd_behavior(:,17:20);
mean_trial_duration_sfig15d = mean(trial_duration_sfig15d,1);
sem = std(trial_duration_sfig15d,1)/sqrt(size(trial_duration_sfig15d,1));


excel_source_data_sfig15d_trial_duration = vertcat (trial_duration_sfig15d,mean_trial_duration_sfig15d, sem);

[p,t,stats] = anova1(trial_duration_sfig15d);
multicompare_results = multcompare(stats);
%% 
performance_sfig15d = MK801_Cmpd_behavior(:,22:25);
mean_trial_duration_sfig15d_performance = mean(performance_sfig15d,1);
sem = std(performance_sfig15d,1)/sqrt(size(performance_sfig15d,1));


excel_source_data_sfig15d_performance = vertcat (performance_sfig15d,mean_trial_duration_sfig15d_performance, sem);

[p,t,stats] = anova1(performance_sfig15d);
multicompare_results = multcompare(stats);
%% fig4 b

cd ('J:\H50_Disk_D\WM_materials_for_the_paper\Data_Figure4')
Cmpd_behavior_design1 = xlsread('Compound Testing Design 1');
Cmpd_performance_design1 = Cmpd_behavior_design1(:,1:4);
mean_trial_duration_fig4b = mean(Cmpd_performance_design1,1);
sem = std(Cmpd_performance_design1,1)/sqrt(size(Cmpd_performance_design1,1));


excel_source_data_fig4b = vertcat (Cmpd_performance_design1,mean_trial_duration_fig4b, sem);

[p,t,stats] = anova1(Cmpd_performance_design1);
multicompare_results = multcompare(stats);
%% fig4 b
cd ('J:\H50_Disk_D\WM_materials_for_the_paper\Data_Figure4')
MK801_Cmpd_behavior = xlsread('MK801_and_Cmpd');
performance_fig4c = MK801_Cmpd_behavior(:,1:4);
mean_trial_duration_fig4c = mean(performance_fig4c,1);
sem = std(performance_fig4c,1)/sqrt(size(performance_fig4c,1));


excel_source_data_fig4c = vertcat (performance_fig4c,mean_trial_duration_fig4c, sem);

[p,t,stats] = anova1(performance_fig4c);
multicompare_results = multcompare(stats);
