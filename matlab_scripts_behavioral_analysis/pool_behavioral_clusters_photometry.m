cd 'D:\Documents\data\wm_t_maze\photometry\ca_data';

load('leftCorrect.mat')

% create categories out of mice
temp_name =[];
for count_names = 1:numel(all_files)

    temp_name(:,count_names) = str2num(all_files(count_names).name(3:4));
    %temp_name(:,count_names) = str2num(all_trials{1, count_names}.coord.name(3:4));

end
% create a categorical variable
mice_cat.left =categorical(temp_name);
all_cat = categories(mice_cat.left);


clusters.temp_left=cluster_temp_input;
clusters.mean_left=cluster_input;
clusters.trial_size_left= median_trials_size;

% save to test angles
test.pool_angle_tj_LEFT = pool_angle_tj;
test.pool_angle_start_LEFT = pool_angle_start;

test.pool_turn_tj_LEFT = pool_turn_tj;
test.pool_turn_start_LEFT = pool_turn_start;

load('rightCorrect.mat')

test.pool_angle_tj_RIGHT = pool_angle_tj;
test.pool_angle_start_RIGHT = pool_angle_start;

test.pool_turn_tj_RIGHT = pool_turn_tj;
test.pool_turn_start_RIGHT = pool_turn_start;

% calculate frequency of direction to
frequency_of_l = and(test.pool_angle_start_LEFT<30, test.pool_angle_start_LEFT>0);
frequency_of_r = and(test.pool_angle_start_RIGHT<30, test.pool_angle_start_RIGHT>0);

[h,p] =ztest(frequency_of_l,mean(frequency_of_r), std(frequency_of_r));

% calculate frequency of turns to
% Orientation towards T-junction
% "+" sign is clockwise (RIGH Turn) "-" sign is anti-clockwise (LEFT Turn)
% Orientation towards Startbox
% "-" sign is clockwise (RIGH Turn) "+" sign is anti-clockwise (LEFT Turn)
frequency_of_l = test.pool_turn_start_LEFT>0;
frequency_of_r = test.pool_turn_start_RIGHT>0;

[h,p] =ztest(frequency_of_l,mean(frequency_of_r), std(frequency_of_r));

% ---create categories out of mice---
temp_name =[];
for count_names = 1:numel(all_files)

    temp_name(:,count_names) = str2num(all_files(count_names).name(3:4));
    %temp_name(:,count_names) = str2num(all_trials{1, count_names}.coord.name(3:4));
end
% create a categorical variable
mice_cat.right =categorical(temp_name);

clusters.temp_right=cluster_temp_input;
clusters.mean_right=cluster_input;
clusters.trial_size_right= median_trials_size;

% ---resize clusters.temp_right to clusters.temp_left to use it in UMAP---
% pool behavioral variables into clusters
input = [];
input_cat = [];
temp_size = clusters.trial_size_right;

temp_start = 0;
for count_behav = 1:9 % 9 behavioral vectors, defined as inputs# above

    % create vecotor
    temp_stop = temp_start+temp_size(count_behav);
    % size of clusters Trials x Time (time concatenated across all behav vars)
    input= resample(clusters.temp_right(:,temp_start+1:temp_stop)',clusters.trial_size_left(count_behav),temp_size(count_behav));
    input_cat = [input_cat input'];
    temp_start = temp_stop;
    
end
clusters.temp_right_resampled =input_cat;

input_to_cluster = [clusters.temp_left; clusters.temp_right_resampled];
input_trial_type = [zeros(1,size(clusters.temp_left,1)) ones(1,size(clusters.temp_right,1))];
[reduction, umap, clusterIdentifiers, extras]=run_umap(input_to_cluster);
[idx,C] =kmeans(reduction,3);

X = reduction;

figure;
subplot(1,2,1)
gscatter(X(:,1),X(:,2),input_trial_type,'br')
hold on
plot(C(:,1),C(:,2),'kx')
legend('Cluster 1','Cluster 2','Cluster 3','Cluster Centroid')
axis('square')
title('Clusters Left vs Right')

subplot(1,2,2)
input_category = [mice_cat.left mice_cat.right];
gscatter(X(:,1),X(:,2),input_category)
axis('square')
title('Clusters across mice')

%% plot Left vs Right for one mouse

subplot(1,2,1)
ind = find(input_category == '11');
input_to_plot = X(ind,:);
gscatter(input_to_plot(:,1),input_to_plot(:,2),input_trial_type(ind),'br')

axis('square')
title('Clusters Left vs Right')
