%Clean the previous matlab sesion

close all
clearvars
clc

folder_with_matlab_functions = 'H:\H50_Disk_D\WM_materials_for_the_paper\Matlab_Scripts_final'; %give the folder with the matlab functions:
addpath(folder_with_matlab_functions)

folder_with_data_for_averaging = uigetdir('H:\H50_Disk_D\WM_materials_for_the_paper\Data_for_Figures\Data_Sup_Figure8\sorted_data\behavior','Select a data folder for analysis'); %give the folder with the raw data
cd (folder_with_data_for_averaging)


files= dir('*.mat');
[mice_to_pool, ~]=size(files);

for count_mice= 1:mice_to_pool
    
    current_file_name = files(count_mice,:).name;
    current_mouseID = files(count_mice,:).name(1:4);
    % load data by mouse
    load(current_file_name);
    
end
%% 

session_number = 1;
for count_mice = 1:mice_to_pool
    
    current_file_name = files(count_mice,:).name(1:end-4);
    current_mouse_data = eval(current_file_name);
    number_of_sessions = size(fieldnames(current_mouse_data),1);
    
    for current_session = 1:number_of_sessions
        
        session_names = fieldnames (current_mouse_data);
        current_session_name = session_names {current_session};
        
        sessions_to_average {session_number} = current_mouse_data.(current_session_name);
        session_number = session_number + 1;
    end
    
end
%% load miniscope recording data
folder_with_data_for_averaging = uigetdir('H:\WM_materials_for_the_paper\Data_Sup_Figure8\sorted_data\miniscope_recording','Select a folder with miniscope data'); %give the folder with the raw data
cd (folder_with_data_for_averaging)


files_miniscope= dir('*.mat');
[mice_to_pool, ~]=size(files_miniscope);

for count_mice= 1:mice_to_pool
    
    current_file_name = files_miniscope(count_mice,:).name;
    current_mouseID = files_miniscope(count_mice,:).name(1:4);
    % load data by mouse
    load(current_file_name);
end
%% 
session_number = 1;
for count_mice = 1:mice_to_pool
    
    current_file_name = files_miniscope(count_mice,:).name(1:end-4);
    current_miniscope_data = eval(current_file_name);
    number_of_sessions = size(fieldnames(current_miniscope_data),1);
    
    for current_session = 1:number_of_sessions
        
        session_names = fieldnames (current_miniscope_data);
        current_session_name = session_names {current_session};
        
        sessions_signal_to_average{session_number} = sum(current_miniscope_data.(current_session_name),2);
        session_number = session_number + 1;
    end
    
end
%% 

%

trial_left_correct_average = sessions_to_average{1}.Trial_LWhole_Correct;
trial_right_correct_average = sessions_to_average{1}.Trial_RWhole_Correct;
trial_left_mistake_average = sessions_to_average{1}.Trial_LWhole_Mistake;
trial_right_mistake_average = sessions_to_average{1}.Trial_RWhole_Mistake;

for current_session = 2:size(sessions_to_average,2)
    
    trial_left_correct_average = vertcat(trial_left_correct_average,sessions_to_average{1,current_session}.Trial_LWhole_Correct);
    trial_right_correct_average = vertcat(trial_right_correct_average,sessions_to_average{1,current_session}.Trial_RWhole_Correct);
    
    
    if any(ismember(fields(sessions_to_average{1, current_session}),'Trial_LWhole_Mistake'))==1
        trial_left_mistake_average=vertcat(trial_left_mistake_average,sessions_to_average{1,current_session}.Trial_LWhole_Mistake);
    end
    
    if any(ismember(fields(sessions_to_average{1, current_session}),'Trial_RWhole_Mistake'))==1
        trial_right_mistake_average=vertcat(trial_right_mistake_average,sessions_to_average{1,current_session}.Trial_RWhole_Mistake);
    end
    
end

performanceR = 100*size(trial_right_correct_average,1)/(size(trial_right_correct_average,1)+size(trial_right_mistake_average,1));
performanceL = 100*(size(trial_left_correct_average,1))/(size(trial_left_correct_average,1)+size(trial_left_mistake_average,1));

%

names_of_events = {'Corridor' 'End' 'Start'};

trial_left_correct_average_maintanance = trial_left_correct_average(:,[6,7,8,9]);
trial_right_correct_average_maintanance = trial_right_correct_average(:,[6,7,8,9]);
trial_left_mistake_average_maintanance = trial_left_mistake_average(:,[6,7,8,9]);
trial_right_mistake_average_maintanance = trial_right_mistake_average(:,[6,7,8,9]);

sm_dT_RL_Whole = []; % vector of mean dt between "event" time points
sd_dT_RL_Whole = [];
dT_RL_Whole = [];% matrix, colums: delta T between neighboring timepoints, raws: different trials
dt=[]; % vector std of dt between "event" time points

current_event = 1;
for current_event = 1:size(trial_left_correct_average_maintanance,2)-1
    [dt,sm_dt,sd_dt] = NewDeltaT(vertcat(trial_left_correct_average_maintanance(:,current_event),...
        trial_right_correct_average_maintanance(:,current_event),...
        trial_left_mistake_average_maintanance(:,current_event),...
        trial_right_mistake_average_maintanance(:,current_event)),...
        vertcat(trial_left_correct_average_maintanance(:,current_event+1),...
        trial_right_correct_average_maintanance(:,current_event+1),...
        trial_left_mistake_average_maintanance(:,current_event+1),...
        trial_right_mistake_average_maintanance(:,current_event+1))); %counting dT between i and i+1 time point
    
    dT_RL_Whole(current_event,:) = dt;
    sm_dT_RL_Whole(current_event) = sm_dt;
    sd_dT_RL_Whole(current_event) = sd_dt;
end


%%  Correct trials
baseline_end_event_number = 2;
for current_session = 1:size(sessions_to_average,2)
    
    [session_resampled_right_correct_per_session{current_session},dt_resample_R] = ResampledSession_v2D_average(sessions_signal_to_average{current_session},...
                                                                                                       sessions_to_average{current_session}.Trial_RWhole_Correct(:,[6,7,8,9]),...
                                                                                                       sm_dT_RL_Whole,...
                                                                                                       baseline_end_event_number);
    [session_resampled_left_correct__per_session{current_session},dt_resample_L] = ResampledSession_v2D_average(sessions_signal_to_average{current_session},...
                                                                                                       sessions_to_average{current_session}.Trial_LWhole_Correct(:,[6,7,8,9]),...
                                                                                                       sm_dT_RL_Whole,...
                                                                                                       baseline_end_event_number);
    
end

maintanence_resampled_right_correct_all = session_resampled_right_correct_per_session{1};
maintanence_resampled_left_correct_all = session_resampled_left_correct__per_session{1};

for current_session = 2:size(session_resampled_left_correct__per_session,2)
    maintanence_resampled_left_correct_all = vertcat(maintanence_resampled_left_correct_all,session_resampled_left_correct__per_session{current_session});
end

for current_session = 2:size(session_resampled_right_correct_per_session,2)
    maintanence_resampled_right_correct_all = vertcat(maintanence_resampled_right_correct_all,session_resampled_right_correct_per_session{current_session});
end

dt_resample = sum(sm_dT_RL_Whole)/size(maintanence_resampled_left_correct_all,2); %0.0501;
time = -1:dt_resample:dt_resample*(size(maintanence_resampled_left_correct_all,2)-1)-1;

%% mistakes

session_resampled_right_mistake_per_session = [];
session_resampled_left_mistake_per_session = [];
for current_session = 1:size(sessions_to_average,2)-1
    
    if any(ismember(fields(sessions_to_average{1, current_session}),'Trial_RWhole_Mistake'))==1
        [session_resampled_right_mistake_per_session{current_session},dt_resample_R_Mistake] = ResampledSession_v2D_average(sessions_signal_to_average{current_session},...
                                                                                                                    sessions_to_average{current_session}.Trial_RWhole_Mistake(:,[6,7,8,9]),...
                                                                                                                    sm_dT_RL_Whole,...
                                                                                                                    baseline_end_event_number);
    end
    if any(ismember(fields(sessions_to_average{1, current_session}),'Trial_LWhole_Mistake'))==1
        [session_resampled_left_mistake_per_session{current_session},dt_resample_L_Mistake] = ResampledSession_v2D_average(sessions_signal_to_average{current_session},...
                                                                                                                    sessions_to_average{current_session}.Trial_LWhole_Mistake(:,[6,7,8,9]),...
                                                                                                                    sm_dT_RL_Whole,...
                                                                                                                    baseline_end_event_number);
    end
end

% delete sessions without mistakes in SessionResL_Wh_Mistake ( Left: 10, 23, Right: 26,32)

session_resampled_right_mistake_per_session = session_resampled_right_mistake_per_session(~cellfun('isempty',session_resampled_right_mistake_per_session));
session_resampled_left_mistake_per_session = session_resampled_left_mistake_per_session(~cellfun('isempty',session_resampled_left_mistake_per_session));


maintanence_resampled_right_mistake_all = session_resampled_right_mistake_per_session{1};
maintanence_resampled_left_mistake_all = session_resampled_left_mistake_per_session{1};

for current_session=2:size(session_resampled_right_mistake_per_session,2)
    
    maintanence_resampled_right_mistake_all = vertcat(maintanence_resampled_right_mistake_all,...
                                                  session_resampled_right_mistake_per_session{current_session});
end

for current_session=2:size(session_resampled_left_mistake_per_session,2)
    
    maintanence_resampled_left_mistake_all = vertcat(maintanence_resampled_left_mistake_all,...
                                                 session_resampled_left_mistake_per_session{current_session});
end

%% plot correct and mistake trials on the same graph
cases_to_compare = {maintanence_resampled_left_correct_all,...
                    maintanence_resampled_right_correct_all,...
                    maintanence_resampled_left_mistake_all,...
                    maintanence_resampled_right_mistake_all};
names_of_comparissons = {'LC', 'RC', 'LM', 'RM'};   
font_size = 20;
color_input = {'b','r','k',[128,128,128]/255};
cd ('H:\H50_Disk_D\WM_materials_for_the_paper\Data_for_Figures\Data_Sup_Figure8\figures')

for current_input_1 = 1:size (cases_to_compare,2)-1
    
input_1 = cases_to_compare{current_input_1};

for current_input_2 = current_input_1+1:size(cases_to_compare,2)
    
input_2 = cases_to_compare{current_input_2};

NamePlot = [names_of_comparissons{current_input_1}, '_vs_', names_of_comparissons{current_input_2}];

figure
y_axis1 = -2;
y_axis2 = 3;
[smoothed_mean_trial.(names_of_comparissons{current_input_1}), smoothed_sem_mean_trial.(names_of_comparissons{current_input_1})] = ...
                                                              plot_input_shaded_error (input_1, color_input{current_input_1}, time, y_axis1,y_axis2);
hold on
[smoothed_mean_trial.(names_of_comparissons{current_input_2}), smoothed_sem_mean_trial.(names_of_comparissons{current_input_2})] = ...
                                                              plot_input_shaded_error (input_2, color_input{current_input_2}, time, y_axis1,y_axis2);
hold on


for qq=1:length(sm_dT_RL_Whole)
    
    line([sum(sm_dT_RL_Whole(1:qq))-dt_resample*(qq+1) sum(sm_dT_RL_Whole(1:qq))-dt_resample*(qq+1)],[ y_axis1 y_axis2],'Color','k','LineWidth',1,'LineStyle','--')

end

for qq=1:length(names_of_events)
    
    text_handle = text(sum(sm_dT_RL_Whole(1:qq))-dt_resample*(qq+1) ,y_axis2,names_of_events(qq),'FontSize',font_size);
    set(text_handle,'rotation',90)

end

hold on
text_handle = text(0,y_axis2,'T-junction','FontSize',font_size);
set(text_handle,'rotation',90)
title({NamePlot;'';'';''},'Interpreter', 'none');
axis('square')

set (gcf,'Position', [442   181   773   707]);

hold off

print(gcf, NamePlot,'-dpdf','-fillpage')

number_of_draws = 20;
number_of_trials_sampled = 200; 

[AUC.(NamePlot), AUC_shuffle.(NamePlot),...
 h.(NamePlot),adj_p.(NamePlot),...
 smoothed_mean_trial.(NamePlot),...
 smoothed_sem_mean_trial.(NamePlot)] = compare_vectors_ROC (input_1,...
                                                            input_2,...
                                                            time, ...
                                                            number_of_draws, ...
                                                            number_of_trials_sampled,...
                                                            NamePlot);

end
end

%% collect source data to excel files
for current_cases_to_compare = 1:size (cases_to_compare,2)
    
    individual_trials_z_score = cases_to_compare{current_input_1};
    average_z_score = smoothed_mean_trial.(names_of_comparissons{current_cases_to_compare});
    sem = smoothed_sem_mean_trial.(names_of_comparissons{current_input_1});
    
    
    excel_source_data.(names_of_comparissons{current_cases_to_compare}) = vertcat (time,...
                                                                                   average_z_score,...
                                                                                   sem,...
                                                                                   individual_trials_z_score);
end

%collect AUC source data
all_comparissons = fieldnames(AUC);
for current_comparisson = 1:size(all_comparissons)
    
    
%     individual_trials_z_score = cases_to_compare{current_input_1};
    mean_AUC = smoothed_mean_trial.(all_comparissons{current_comparisson}).AUC;
    sem_AUC = smoothed_sem_mean_trial.(all_comparissons{current_comparisson}).AUC;
    
    mean_AUC_shuffle = smoothed_mean_trial.(all_comparissons{current_comparisson}).AUC_shuffle;
    sem_AUC_shuffle = smoothed_sem_mean_trial.(all_comparissons{current_comparisson}).AUC_shuffle;
    adjusted_p_value = adj_p.(all_comparissons{current_comparisson});
    h_significance_vector = h.(all_comparissons{current_comparisson});
    individual_draws_AUC = AUC.(all_comparissons{current_comparisson})';
    individual_draws_AUC_shuffle = AUC_shuffle.(all_comparissons{current_comparisson})';
    
    excel_source_data_AUC.(all_comparissons{current_comparisson}) = vertcat (time,...
                                                                             mean_AUC,...
                                                                             sem_AUC,...
                                                                             mean_AUC_shuffle,...
                                                                             sem_AUC_shuffle,...
                                                                             adjusted_p_value,...
                                                                             h_significance_vector,...
                                                                             individual_draws_AUC,...
                                                                             individual_draws_AUC_shuffle);
    
end

save('variables_15082023','excel_source_data_AUC','excel_source_data','AUC','AUC_shuffle')
