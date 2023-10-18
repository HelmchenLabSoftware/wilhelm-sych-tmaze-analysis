function [sessions_to_average, number_of_sessions] = extract_session (folder_with_data_for_averaging)

cd (folder_with_data_for_averaging)


files= dir('*.mat');
[mice_to_pool, ~]=size(files);

for count_mice= 1:mice_to_pool
    
    current_file_name = files(count_mice,:).name;
    current_mouseID = files(count_mice,:).name(1:4);
    % load data by mouse
    load(current_file_name);
    
end

% 

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
end