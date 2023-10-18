function sm_dT_RL_Whole = average_time_between_events (sessions_to_average)

Trial_LWhole_Correct_Average = sessions_to_average{1}.Trial_LWhole_Correct;
Trial_RWhole_Correct_Average = sessions_to_average{1}.Trial_RWhole_Correct;


for current_session = 2:size(sessions_to_average,2)

 current_Trial_RWhole_Correct_Average=[];
 current_Trial_LWhole_Correct_Average=[];
 
current_Trial_LWhole_Correct_Average = sessions_to_average{1,current_session}.Trial_LWhole_Correct;
current_Trial_RWhole_Correct_Average = sessions_to_average{1,current_session}.Trial_RWhole_Correct;

Trial_LWhole_Correct_Average = vertcat(Trial_LWhole_Correct_Average,current_Trial_LWhole_Correct_Average);
Trial_RWhole_Correct_Average = vertcat(Trial_RWhole_Correct_Average,current_Trial_RWhole_Correct_Average);

end


Trial_LWhole_Correct_Average(:,[4,7,11,13,16]) = [];

Trial_RWhole_Correct_Average(:,[4,7,11,13,16]) = [];

Trial_LWhole_Correct_Average_3Phases_downsampled = Trial_LWhole_Correct_Average(:,(1:11));
Trial_RWhole_Correct_Average_3Phases_downsampled = Trial_RWhole_Correct_Average(:,(1:11));

sm_dT_RL_Whole=[]; % vector of mean dt between "event" time points
sd_dT_RL_Whole=[];
dT_RL_Whole=[];% matrix, colums: delta T between neighboring timepoints, raws: different trials
dt=[]; % vector std of dt between "event" time points


Trial_LWhole=Trial_LWhole_Correct_Average_3Phases_downsampled;
Trial_RWhole=Trial_RWhole_Correct_Average_3Phases_downsampled;

current_point=1;
 for current_point=1:size(Trial_LWhole,2)-1
[dt,sm_dt,sd_dt]=NewDeltaT(vertcat(Trial_LWhole(:,current_point),Trial_RWhole(:,current_point)),vertcat(Trial_LWhole(:,current_point+1),Trial_RWhole(:,current_point+1))); %counting dT between i and i+1 time point
dT_RL_Whole(current_point,:)=dt;
sm_dT_RL_Whole(current_point)=sm_dt;
sd_dT_RL_Whole(current_point)=sd_dt;
 end 
 
end