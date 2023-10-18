function [SessionResR_Wh_all_photometry,...
         SessionResL_Wh_all_photometry,...
         SessionResR_Wh,...
         SessionResL_Wh] = calculate_resampled_signals (sessions_to_average, sm_dT_RL_Whole, baseline_end_event_number)

 for current_session = 1:size(sessions_to_average,2) 
     
Trial_RWhole_Correct_Average = [];
Trial_LWhole_Correct_Average = [];
current_Trial_LWhole_Correct_Average = [];
current_Trial_RWhole_Correct_Average = [];
     
current_Trial_LWhole_Correct_Average = sessions_to_average{1,current_session}.Trial_LWhole_Correct;

current_Trial_LWhole_Correct_Average(:,[4,7,11,13,16]) = [];
Trial_LWhole_Correct_Average = current_Trial_LWhole_Correct_Average(:,1:11);
     
current_Trial_RWhole_Correct_Average = sessions_to_average{1,current_session}.Trial_RWhole_Correct;
current_Trial_RWhole_Correct_Average(:,[4,7,11,13,16]) = [];
Trial_RWhole_Correct_Average = current_Trial_RWhole_Correct_Average(:,1:11);

[SessionResR_Wh{current_session},dTR_resample] = ResampledSession_v2(sessions_to_average{current_session}.demodSig488,Trial_RWhole_Correct_Average,sm_dT_RL_Whole,baseline_end_event_number);
[SessionResL_Wh{current_session},dTL_resample] = ResampledSession_v2(sessions_to_average{current_session}.demodSig488,Trial_LWhole_Correct_Average,sm_dT_RL_Whole,baseline_end_event_number);

 end
%  
 SessionResR_Wh_all_photometry = SessionResR_Wh{1};
 SessionResL_Wh_all_photometry = SessionResL_Wh{1};
 
 for current_session=2:size(sessions_to_average,2)
     
     SessionResR_Wh_all_photometry = vertcat(SessionResR_Wh_all_photometry,SessionResR_Wh{current_session});
     SessionResL_Wh_all_photometry = vertcat(SessionResL_Wh_all_photometry,SessionResL_Wh{current_session});
     
 end
 
end