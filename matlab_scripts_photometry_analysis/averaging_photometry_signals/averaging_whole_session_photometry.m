%This function allows to average photometry signals collected during T-maze
%alternation task

%Data should be pre-processed (e.g. Trial structure extracted and sorted
%for each mouse

% Author: Maria Wilhelm

%Clean the previous matlab sesion

close all
clearvars
clc

folder_with_matlab_functions = 'D:\H50_Disk_D\WM_materials_for_the_paper\Matlab_Scripts_final'; %give the folder with the matlab functions:
addpath(folder_with_matlab_functions)

folder_with_all_data_for_averaging = ('D:\H50_Disk_D\WM_materials_for_the_paper\Data_for_Figures\Data_Figure2\sorted_data\individual_mice_data\GCaMP_and_GFP'); %give the folder with all raw data
[sessions_to_average, number_of_sessions] = extract_session (folder_with_all_data_for_averaging);

% Averaging delta times between all events in all sessions (including both GCaMP
% and GFP injected mice)


sm_dT_RL_Whole = average_time_between_events (sessions_to_average);

clearvars -except sm_dT_RL_Whole

% 
folder_for_figures = 'J:\H50_Disk_D\WM_materials_for_the_paper\Data_for_Figures\Data_Figure2\figures\figures_backup';
folder_with_data_for_averaging.GCaMP = ('D:\H50_Disk_D\WM_materials_for_the_paper\Data_for_Figures\Data_Figure2\sorted_data\individual_mice_data\GCaMP'); %give the folder with the raw data
[sessions_to_average_GCaMP, number_of_sessions.GCaMP] = extract_session (folder_with_data_for_averaging.GCaMP);

% 

folder_with_data_for_averaging.GFP = ('J:\H50_Disk_D\WM_materials_for_the_paper\Data_for_Figures\Data_Figure2\sorted_data\individual_mice_data\GFP'); %give the folder with the raw data
[sessions_to_average_GFP, number_of_sessions.GFP] = extract_session (folder_with_data_for_averaging.GFP);

%% calculate trial to trial signal resampled to average time between events
baseline_end_event_number = 6; % in the original trial structure, the even before the t-junction is 6 (correspond to event 5 in the paper)
[SessionResR_Wh_all_photometry_GCaMP,...
 SessionResL_Wh_all_photometry_GCaMP,...
 SessionResR_Wh.GCaMP,...
SessionResL_Wh.GCaMP] = calculate_resampled_signals (sessions_to_average_GCaMP, sm_dT_RL_Whole, baseline_end_event_number);

[SessionResR_Wh_all_photometry_GFP,...
 SessionResL_Wh_all_photometry_GFP,...
 SessionResR_Wh.GFP,...
SessionResL_Wh.GFP] = calculate_resampled_signals (sessions_to_average_GFP, sm_dT_RL_Whole, baseline_end_event_number);
 %% Plot averaged signal
name_plot = 'Photometry all mice RandL';

color_input1=[0,128,0]/255;
transparent=0;
y_axis1=-2;
y_axis2=4;

hF=figure;
hold on
Session_GCaMP_to_plot = vertcat(SessionResL_Wh_all_photometry_GCaMP,SessionResR_Wh_all_photometry_GCaMP);
MeanTrial_GCaMP = mean(Session_GCaMP_to_plot,1); 
sem_MeanTrial_GCaMP = std(Session_GCaMP_to_plot,1)/sqrt(size(Session_GCaMP_to_plot,1));% std(x)/sqrt(length(x));

dt_resample = (1+sum(sm_dT_RL_Whole(1:end)))/length(MeanTrial_GCaMP);
t= -1:dt_resample:dt_resample*(length(MeanTrial_GCaMP)-1)-1;
x_axis1=-1;
x_axis2=max(t);
H=NewShadedErrorBar(t,movmean(MeanTrial_GCaMP,5),movmean(sem_MeanTrial_GCaMP,5),'k',transparent,21,y_axis1,y_axis2,-1,max(t),'','',size(Session_GCaMP_to_plot,1),'k') ;
hold on

Session_GFP_to_plot=vertcat(SessionResL_Wh_all_photometry_GFP,SessionResR_Wh_all_photometry_GFP);
MeanTrial_GFP=mean(Session_GFP_to_plot,1); 
sem_MeanTrial_GFP=std(Session_GFP_to_plot,1)/sqrt(size(Session_GFP_to_plot,1));% std(x)/sqrt(length(x));

dt_resample = (1+sum(sm_dT_RL_Whole(1:end)))/length(MeanTrial_GFP);
t= -1:dt_resample:dt_resample*(length(MeanTrial_GFP)-1)-1;
x_axis1=-1;
x_axis2=max(t);
H=NewShadedErrorBar(t,movmean(MeanTrial_GFP,5),movmean(sem_MeanTrial_GFP,5),'k',transparent,21,y_axis1,y_axis2,-1,max(t),'','',size(Session_GFP_to_plot,1),'g') ;
hold on
for qq=1:length(sm_dT_RL_Whole)
line([sum(sm_dT_RL_Whole(1:qq))-dt_resample*(qq+1) sum(sm_dT_RL_Whole(1:qq))-dt_resample*(qq+1)],[ y_axis1 y_axis2],'Color','k','LineWidth',1,'LineStyle','--')
end

hold on
x_axis1=-1;
x_axis2=sum(sm_dT_RL_Whole)-1;
 
axis([x_axis1 x_axis2 y_axis1 y_axis2])
h=gcf;
set(h,'PaperOrientation','landscape');
h_title=title({name_plot;'';'';'';'';'';''});
axis('square')
set(gca,'Position',[0.1300    0.1100    0.7750    0.5455])
cd(folder_for_figures)
save_plots(name_plot)
 

%% collect source data to excel files

    individual_trials_z_score = Session_GCaMP_to_plot;
    average_z_score = movmean(MeanTrial_GCaMP,5);
    sem = movmean(sem_MeanTrial_GCaMP,5);
    
    
    excel_source_data_GCaMP = vertcat (t,...
                                 average_z_score,...
                                 sem,...
                                 individual_trials_z_score);
                             
    individual_trials_z_score = Session_GFP_to_plot;
    average_z_score = movmean(MeanTrial_GFP,5);
    sem = movmean(sem_MeanTrial_GFP,5);
    
    
    excel_source_data_GFP = vertcat (t,...
                                 average_z_score,...
                                 sem,...
                                 individual_trials_z_score);

%% COMPARISON WITHIN EACH MOUSE AND FISHER´S METHOD FOR META-ANALYSIS
sensor_name = 'GCaMP'; %change to 'GFP' for the GFP plot
sessions_to_quatify = sessions_to_average_GFP;

time_interval_encoding = [1,sum(sm_dT_RL_Whole(1:4))];
index_encoding_start = round(time_interval_encoding(1)/dt_resample);
index_encoding_end = round(time_interval_encoding(2)/dt_resample);

time_interval_maintenance = [sum(sm_dT_RL_Whole(1:4)),sum(sm_dT_RL_Whole(1:6))];   
index_maintenance_start = round(time_interval_maintenance(1)/dt_resample);
index_maintenance_end = round(time_interval_maintenance(2)/dt_resample);

time_interval_retrieval = [sum(sm_dT_RL_Whole(1:6)),sum(sm_dT_RL_Whole(1:10))];   
index_retrieval_start = round(time_interval_retrieval(1)/dt_resample);
index_retrieval_end = round(time_interval_retrieval(2)/dt_resample);


all_trials_per_mouse = [];  
cd (folder_with_data_for_averaging.(sensor_name))
files= dir('*.mat');
[mice_to_pool, ~] = size(files);

for current_mouse= 1:mice_to_pool
    
   sessions_range = [1+sum(number_of_sessions.(sensor_name)(1:current_mouse-1)):sum(number_of_sessions.(sensor_name)(1:current_mouse))];  
   all_trials_per_mouse{current_mouse} = vertcat(SessionResL_Wh.(sensor_name){sessions_range(1)},SessionResR_Wh.(sensor_name){sessions_range(1)});
   
   for current_session = sessions_range(2):sessions_range(end)
    
   all_trials_per_mouse{current_mouse} = vertcat(all_trials_per_mouse{current_mouse},...
                                                 SessionResL_Wh.(sensor_name){current_session},...
                                                 SessionResR_Wh.(sensor_name){current_session});
   end
    
mean_per_trial.(sensor_name).encoding{current_mouse} = mean(all_trials_per_mouse{current_mouse}(:,index_encoding_start:index_encoding_end),2);

mean_per_trial.(sensor_name).maintenance{current_mouse} = mean(all_trials_per_mouse{current_mouse}(:,index_maintenance_start:index_maintenance_end),2);
mean_per_trial.(sensor_name).retrieval{current_mouse} = mean(all_trials_per_mouse{current_mouse}(:,index_retrieval_start:index_retrieval_end),2);  
   
   
   
figure

input_data = [{mean_per_trial.(sensor_name).encoding{current_mouse}},...
              {mean_per_trial.(sensor_name).maintenance{current_mouse}},...
              {mean_per_trial.(sensor_name).retrieval{current_mouse}}];
          
color_cohort = [{[0.59,0.53,0.82]},{[0.87,0.65,0.37]},{[0.53,0.82,0.82]}];
legend_cohort = [{'Enc'},{'M'},{'Ret'}];
plot_title = ['Mean per Period individual Comparisons_', sensor_name, '_mouse',num2str(current_mouse)];
font_size = 20;
paired = 1;

SEM_input{current_mouse} = bar_plots_cohorts(input_data,legend_cohort,color_cohort, plot_title,'z-score',font_size, paired);
set(gca,'Position',[0.1300    0.1100    0.7750    0.7417])
save_plots(plot_title)

[p,t,stats] = anova1([input_data{1}, input_data{2}, input_data{3}]);
multicompare_results = multcompare(stats);

p_values_corrected_Tukey.(sensor_name).encoding_vs_maintenance{current_mouse} = multicompare_results(1,6);
p_values_corrected_Tukey.(sensor_name).encoding_vs_retrieval{current_mouse} = multicompare_results(2,6);
p_values_corrected_Tukey.(sensor_name).maintenance_vs_retrieval{current_mouse} = multicompare_results(3,6);
end

% meta analysis:grouped p-values calculated with Fisher´s method

group_pval.encoding_vs_maintenance = fisher_pvalue_meta_analysis([p_values_corrected_Tukey.(sensor_name).encoding_vs_maintenance{:}]);

group_pval.encoding_vs_retrieval = fisher_pvalue_meta_analysis([p_values_corrected_Tukey.(sensor_name).encoding_vs_retrieval{:}]);

group_pval.maintenance_vs_retrieval = fisher_pvalue_meta_analysis([p_values_corrected_Tukey.(sensor_name).maintenance_vs_retrieval{:}]);

%% source data for Fig 2 b
for current_mouse= 1:mice_to_pool
excel_source_data_Fig2_panel_b.(sensor_name){current_mouse} = vertcat([mean(mean_per_trial.(sensor_name).encoding{current_mouse}),mean(mean_per_trial.(sensor_name).maintenance{current_mouse}),mean(mean_per_trial.(sensor_name).retrieval{current_mouse})],...
                                                        SEM_input{current_mouse},...
                                                        [p_values_corrected_Tukey.(sensor_name).encoding_vs_maintenance{current_mouse}, p_values_corrected_Tukey.(sensor_name).maintenance_vs_retrieval{current_mouse}, p_values_corrected_Tukey.(sensor_name).encoding_vs_retrieval{current_mouse}],...
                                                        [mean_per_trial.(sensor_name).encoding{current_mouse}, mean_per_trial.(sensor_name).maintenance{current_mouse}, mean_per_trial.(sensor_name).retrieval{current_mouse}]);
end                                                        
%% Fig 2 panel b
if sensor_name == 'GCaMP'
input_1 = vertcat(mean(mean_per_trial.(sensor_name).encoding{1}), mean(mean_per_trial.(sensor_name).encoding{2}), mean(mean_per_trial.(sensor_name).encoding{3}), mean(mean_per_trial.(sensor_name).encoding{4}), mean(mean_per_trial.(sensor_name).encoding{5}), mean(mean_per_trial.(sensor_name).encoding{6}));
input_2 = vertcat(mean(mean_per_trial.(sensor_name).maintenance{1}), mean(mean_per_trial.(sensor_name).maintenance{2}), mean(mean_per_trial.(sensor_name).maintenance{3}), mean(mean_per_trial.(sensor_name).maintenance{4}), mean(mean_per_trial.(sensor_name).maintenance{5}), mean(mean_per_trial.(sensor_name).maintenance{6}));
input_3 = vertcat(mean(mean_per_trial.(sensor_name).retrieval{1}), mean(mean_per_trial.(sensor_name).retrieval{2}), mean(mean_per_trial.(sensor_name).retrieval{3}), mean(mean_per_trial.(sensor_name).retrieval{4}), mean(mean_per_trial.(sensor_name).retrieval{5}), mean(mean_per_trial.(sensor_name).retrieval{6}));
elseif sensor_name == 'GFP'
input_1 = vertcat(mean(mean_per_trial.(sensor_name).encoding{1}), mean(mean_per_trial.(sensor_name).encoding{2}), mean(mean_per_trial.(sensor_name).encoding{3}), mean(mean_per_trial.(sensor_name).encoding{4}));
input_2 = vertcat(mean(mean_per_trial.(sensor_name).maintenance{1}), mean(mean_per_trial.(sensor_name).maintenance{2}), mean(mean_per_trial.(sensor_name).maintenance{3}), mean(mean_per_trial.(sensor_name).maintenance{4}));
input_3 = vertcat(mean(mean_per_trial.(sensor_name).retrieval{1}), mean(mean_per_trial.(sensor_name).retrieval{2}), mean(mean_per_trial.(sensor_name).retrieval{3}), mean(mean_per_trial.(sensor_name).retrieval{4}));
end
figure

input_data = [{input_1},...
              {input_2},...
              {input_3}];
          
color_cohort = [{[0.59,0.53,0.82]},{[0.87,0.65,0.37]},{[0.53,0.82,0.82]}];
legend_cohort = [{'Enc'},{'M'},{'Ret'}];
plot_title = ['Mean per Period per Mouse', sensor_name];
font_size = 20;
paired = 1;
mean_values = [mean(input_1), mean(input_2),mean(input_3)];
SEM_Input = bar_plots_cohorts(input_data,legend_cohort,color_cohort, plot_title,'z-score',font_size, paired);
%% source data Fig 1

performance.saline = [93,78,85,88,88,74, 75, 82, 79, 92, 94];
performance.MK801 = [52,53,47,60,73,48,52,68,81,72,66];
[p.performance,h.performance] = signrank(performance.saline,performance.MK801);

trial_duration.saline = [109.0185
152.426
106.504
83.5315
96.513
118.573
85.895
173.9535
173.541
106.279
87.122];


trial_duration.MK801 = [64.357
56.113
59.8355
47.18
55.941
57.4885
47.1995
60.334
152.077
78.659
56.1945];

[p.trial_duration,h.trial_duration] = signrank(trial_duration.saline,trial_duration.MK801);

input_1 = performance.saline;
input_2 = performance.MK801;
input_data = [{input_1'},...
              {input_2'}];
mean_values = [mean(input_1), mean(input_2)];
SEM_Input = bar_plots_cohorts(input_data,[{'Veh'},{'MK801'}],[{'k'},{'b'}], 'Performance','',12, 1);
%% SFig1a
sensor_name = 'GCaMP';
all_trials_per_mouse = [];  
performance_vector = [];
performance_per_session =[];
cd (folder_with_data_for_averaging.(sensor_name))
files= dir('*.mat');
[mice_to_pool, ~] = size(files);
performance_per_mouse = [];
for current_mouse= 1:mice_to_pool
    
   sessions_range = [1+sum(number_of_sessions.(sensor_name)(1:current_mouse-1)):sum(number_of_sessions.(sensor_name)(1:current_mouse))];  
   
   session_number_per_mouse = 1;
   for current_session = sessions_range(1):sessions_range(end)
       
    performance_vector = vertcat(sessions_to_average_GCaMP{1, current_session}.sequenceOfSides{:, 6});
    performance_per_session (session_number_per_mouse) = 100*performance_vector(end); 
    session_number_per_mouse = session_number_per_mouse+1;  
   end
   performance_vector = [];
   performance_per_mouse {current_mouse} = performance_per_session';
   mean_performance_per_mouse(current_mouse) = mean (performance_per_mouse {current_mouse});
   sem_performance_per_mouse(current_mouse) = std(performance_per_mouse {current_mouse})/sqrt(length(performance_per_mouse {current_mouse}));
   performance_vector = [];
   performance_per_session =[];
end
dark_grey = [0.4 0.4 0.4];
figure
e=errorbar(1:1:6,mean_performance_per_mouse,sem_performance_per_mouse,'o','MarkerSize',10,'MarkerEdgeColor',dark_grey,'MarkerFaceColor',dark_grey,'LineWidth',2);
e.Color = dark_grey; 
hold on
axis([0 7 50 100])
hold on
color_grey = [0.8 0.8 0.8];

for current_mouse= 1:mice_to_pool
    scatter(current_mouse*ones(1,length(performance_per_mouse {current_mouse}))+0.2,performance_per_mouse {current_mouse},'MarkerEdgeColor',color_grey,'MarkerFaceColor',color_grey)
end
hold off
yticks(50:10:100);

title ('Mice Performance')
ylabel ('Percent Correct')


%% 

event_number_maintenance_start = 6; %6
event_number_retrieval_start = 9; %9
event_numbuer_retrieval_end = 11;

for i=1:size(sessions_to_average_GCaMP,2)
Trial_LWhole_Correct=sessions_to_average_GCaMP{1,i}.Trial_LWhole_Correct;
Trial_RWhole_Correct=sessions_to_average_GCaMP{1,i}.Trial_RWhole_Correct;

    for ii=1:size(Trial_LWhole_Correct,1)
Encoding_dt_LC(ii)=Trial_LWhole_Correct(ii,event_number_maintenance_start)-Trial_LWhole_Correct(ii,1);
Delay_dt_LC(ii)=Trial_LWhole_Correct(ii,event_number_retrieval_start)-Trial_LWhole_Correct(ii,event_number_maintenance_start);
Retrieval_dt_LC(ii)=Trial_LWhole_Correct(ii,15)-Trial_LWhole_Correct(ii,event_number_retrieval_start);
    end
    
Encoding_dt_LC_Session(i)=mean(Encoding_dt_LC);
Delay_dt_LC_Session(i)=mean(Delay_dt_LC);
Retrieval_dt_LC_Session(i)=mean(Retrieval_dt_LC);

    for ii=1:size(Trial_RWhole_Correct,1)       
Encoding_dt_RC(ii)=Trial_RWhole_Correct(ii,event_number_maintenance_start)-Trial_RWhole_Correct(ii,1);
Delay_dt_RC(ii)=Trial_RWhole_Correct(ii,event_number_retrieval_start)-Trial_RWhole_Correct(ii,event_number_maintenance_start);
Retrieval_dt_RC(ii)=Trial_RWhole_Correct(ii,15)-Trial_RWhole_Correct(ii,event_number_retrieval_start);
    end
Encoding_dt_RC_Session(i)=mean(Encoding_dt_RC);
Delay_dt_RC_Session(i)=mean(Delay_dt_RC);
Retrieval_dt_RC_Session(i)=mean(Retrieval_dt_RC);

Encoding_dt_AllC=([Encoding_dt_LC Encoding_dt_RC]);
Delay_dt_AllC=([Delay_dt_LC Delay_dt_RC]);
Retrieval_dt_AllC=([Retrieval_dt_LC Retrieval_dt_RC]);

Encoding_dt_AllC_Session(i)=mean(Encoding_dt_AllC);
Delay_dt_AllC_Session(i)=mean(Delay_dt_AllC);
Retrieval_dt_AllC_Session(i)=mean(Retrieval_dt_AllC);
clearvars Trial_LWhole_Correct Trial_RWhole_Correct Encoding_dt_LC Delay_dt_LC Retrieval_dt_LC Encoding_dt_RC Delay_dt_RC Retrieval_dt_RC
end

MedianEncodingCorrectSession=mean(Encoding_dt_AllC_Session);
MedianDelayCorrectSession=mean(Delay_dt_AllC_Session);
MedianRetrievalCorrectSession=mean(Retrieval_dt_AllC_Session);

var_EncodingCorrectSession=std(Encoding_dt_AllC_Session)/sqrt(length(Encoding_dt_AllC_Session));
var_DelayCorrectSession=std(Delay_dt_AllC_Session)/sqrt(length(Delay_dt_AllC_Session));
var_RetrievalCorrectSession=std(Retrieval_dt_AllC_Session)/sqrt(length(Retrieval_dt_AllC_Session));

meanMouse_C=horzcat(MedianEncodingCorrectSession,MedianDelayCorrectSession,MedianRetrievalCorrectSession);
varMouse_C=horzcat(var_EncodingCorrectSession,var_DelayCorrectSession,var_RetrievalCorrectSession);

%mistake
for i=1:size(sessions_to_average_GCaMP,2)
if any(ismember(fields(sessions_to_average_GCaMP{1, i}),'Trial_LWhole_Mistake'))==1
Trial_LWhole_Mistake=sessions_to_average_GCaMP{1,i}.Trial_LWhole_Mistake;


    for ii=1:size(Trial_LWhole_Mistake,1)
Encoding_dt_LM(ii)=Trial_LWhole_Mistake(ii,event_number_maintenance_start)-Trial_LWhole_Mistake(ii,1);
Delay_dt_LM(ii)=Trial_LWhole_Mistake(ii,event_number_retrieval_start)-Trial_LWhole_Mistake(ii,event_number_maintenance_start);
Retrieval_dt_LM(ii)=Trial_LWhole_Mistake(ii,12)-Trial_LWhole_Mistake(ii,event_number_retrieval_start);
    end
else
    Encoding_dt_LM=[];
    Delay_dt_LM=[];
    Retrieval_dt_LM=[];
end
Encoding_dt_LM_Session(i)=mean(Encoding_dt_LM);
Delay_dt_LM_Session(i)=mean(Delay_dt_LM);
Retrieval_dt_LM_Session(i)=mean(Retrieval_dt_LM);


if any(ismember(fields(sessions_to_average_GCaMP{1, i}),'Trial_RWhole_Mistake'))==1
Trial_RWhole_Mistake=sessions_to_average_GCaMP{1,i}.Trial_RWhole_Mistake;

    for ii=1:size(Trial_RWhole_Mistake,1)       
Encoding_dt_RM(ii)=Trial_RWhole_Mistake(ii,event_number_maintenance_start)-Trial_RWhole_Mistake(ii,1);
Delay_dt_RM(ii)=Trial_RWhole_Mistake(ii,event_number_retrieval_start)-Trial_RWhole_Mistake(ii,event_number_maintenance_start);
Retrieval_dt_RM(ii)=Trial_RWhole_Mistake(ii,12)-Trial_RWhole_Mistake(ii,event_number_retrieval_start);
    end
    else
    Encoding_dt_RM=[];
    Delay_dt_RM=[];
    Retrieval_dt_RM=[];
end
Encoding_dt_RM_Session(i)=mean(Encoding_dt_RM);
Delay_dt_RM_Session(i)=mean(Delay_dt_RM);
Retrieval_dt_RM_Session(i)=mean(Retrieval_dt_RM);

Encoding_dt_AllM=([Encoding_dt_LM Encoding_dt_RM]);
Delay_dt_AllM=([Delay_dt_LM Delay_dt_RM]);
Retrieval_dt_AllM=([Retrieval_dt_LM Retrieval_dt_RM]);

Encoding_dt_AllM_Session(i)=mean(Encoding_dt_AllM);
Delay_dt_AllM_Session(i)=mean(Delay_dt_AllM);
Retrieval_dt_AllM_Session(i)=mean(Retrieval_dt_AllM);
clearvars Trial_LWhole_Mistake Trial_RWhole_Mistake Encoding_dt_LM Delay_dt_LM Retrieval_dt_LM Encoding_dt_RM Delay_dt_RM Retrieval_dt_RM
end

MedianEncodingMistakeSession=mean(Encoding_dt_AllM_Session);
MedianDelayMistakeSession=mean(Delay_dt_AllM_Session);
MedianRetrievalMistakeSession=mean(Retrieval_dt_AllM_Session);

var_EncodingMistakeSession=std(Encoding_dt_AllM_Session)/sqrt(length(Encoding_dt_AllM_Session));
var_DelayMistakeSession=std(Delay_dt_AllM_Session)/sqrt(length(Delay_dt_AllM_Session));
var_RetrievalMistakeSession=std(Retrieval_dt_AllM_Session)/sqrt(length(Retrieval_dt_AllM_Session));

meanMouse_M=horzcat(MedianEncodingMistakeSession,MedianDelayMistakeSession,MedianRetrievalMistakeSession);
varMouse_M=horzcat(var_EncodingMistakeSession,var_DelayMistakeSession,var_RetrievalMistakeSession);
figure
bar(1:2:5,meanMouse_C(1:3),0.3)
%,'s','MarkerSize',20,'MarkerEdgeColor','green','MarkerFaceColor','green','LineWidth',1);

hold on
bar(2:2:6,meanMouse_M(1:3),0.3)

% axis('square')
hold on
%[ones(size(ChoiceDurationOFF_AllMice))-0.1+0.2*rand(size(ChoiceDurationOFF_AllMice)); 2*ones(size(ChoiceDurationON_AllMice))-0.1+0.2*rand(size(ChoiceDurationON_AllMice)) 
plot([ones(size(Encoding_dt_AllC_Session))-0.15+0.3*rand(size(Encoding_dt_AllC_Session)); 3*ones(size(Delay_dt_AllC_Session))-0.15+0.3*rand(size(Delay_dt_AllC_Session)); 5*ones(size(Retrieval_dt_AllC_Session))-0.15+0.3*rand(size(Retrieval_dt_AllC_Session))],[Encoding_dt_AllC_Session; Delay_dt_AllC_Session; Retrieval_dt_AllC_Session],'-o',...
    'MarkerSize',8,'LineStyle', 'none' ,...
    'MarkerEdgeColor','g','MarkerFaceColor',[1,1,1]);
hold on
errorbar(1:2:5,meanMouse_C(1:3) ,varMouse_C(1:3),'LineStyle', 'none');
hold on
plot([2*ones(size(Encoding_dt_AllM_Session))-0.15+0.3*rand(size(Encoding_dt_AllM_Session)); 4*ones(size(Delay_dt_AllM_Session))-0.15+0.3*rand(size(Delay_dt_AllM_Session)); 6*ones(size(Retrieval_dt_AllM_Session))-0.15+0.3*rand(size(Retrieval_dt_AllM_Session))],[Encoding_dt_AllM_Session; Delay_dt_AllM_Session; Retrieval_dt_AllM_Session],'-o',...
    'MarkerSize',8,'LineStyle', 'none' ,...
    'MarkerEdgeColor','k','MarkerFaceColor',[1,1,1]);
hold on
errorbar(2:2:6,meanMouse_M(1:3) ,varMouse_M(1:3),'LineStyle', 'none');%,'s','MarkerSize',20,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',1);


[p,t,stats] = anova1([Encoding_dt_AllC_Session',...
                      Encoding_dt_AllM_Session',...
                      Delay_dt_AllC_Session',...
                      Delay_dt_AllM_Session',...
                      Retrieval_dt_AllC_Session'...
                      Encoding_dt_AllM_Session']);
multicompare_results = multcompare(stats);


excel_source_data_SFig1b = horzcat(Encoding_dt_AllC_Session',...
                          Encoding_dt_AllM_Session',...
                          Delay_dt_AllC_Session',...
                          Delay_dt_AllM_Session',...
                          Retrieval_dt_AllC_Session',...
                          Encoding_dt_AllM_Session');
mean_all_perdids = [meanMouse_C(1),...
                    meanMouse_M(1),...
                    meanMouse_C(2),...
                    meanMouse_M(3),...
                    meanMouse_C(3),...
                    meanMouse_M(3)];
sem_all_perdids = [ varMouse_C(1),...
                    varMouse_M(1),...
                    varMouse_C(2),...
                    varMouse_M(3),...
                    varMouse_C(3),...
                    varMouse_M(3)];
                
cases_to_compare = {'EC', 'EM', 'DC', 'DM', 'RC', 'RM'};               
for each_p_value = 1:size (multicompare_results,1)
    
    p_value{each_p_value} = ...
                            ['p_value(',char(cases_to_compare(multicompare_results(each_p_value,1))),'_',...
                              char(cases_to_compare(multicompare_results(each_p_value,2))),...
                              ') = ',...
                              num2str(multicompare_results(each_p_value,6)), '; '];
end

all_p_values_string = [p_value{:}];