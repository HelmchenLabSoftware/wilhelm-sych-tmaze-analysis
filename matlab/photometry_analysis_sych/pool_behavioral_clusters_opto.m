clear
cd 'D:\Documents\data\wm_t_maze\opto_activation\b_mdl_results';

%load('leftCorrect.mat')
load('leftCorrect.mat')
all_trial_names= [];
all_trial_names= [all_trial_names; all_files];

% create categories out of mice
temp_name =[];
all_trial_durations = [];

for count_names = 1:numel(all_files)

    if ~isempty(b_mdl{1, count_names})
        temp_name(:,count_names) = str2num(all_files(count_names).name(3:4));
        %temp_name(:,count_names) = str2num(all_trials{1, count_names}.coord.name(3:4));
        trial_durations(count_names) = numel(b_mdl{1, count_names}.position_y);  
    end

end
all_trial_durations = [all_trial_durations trial_durations];

% create a categorical variable
mice_cat.left =categorical(temp_name);
all_cat = categories(mice_cat.left);


clusters.temp_left=cluster_temp_input;
clusters.mean_left=cluster_input;
clusters.trial_size_left= [median_trials_size 1];

% save to test angles
test.pool_angle_tj_LEFT = pool_angle_tj;
test.pool_angle_start_LEFT = pool_angle_start;

test.pool_turn_tj_LEFT = pool_turn_tj;
test.pool_turn_start_LEFT = pool_turn_start;

cd 'D:\Documents\data\wm_t_maze\opto_activation\b_mdl_results';
load('optoleftCorrect.mat')

all_trial_names= [all_trial_names; all_files];

test.pool_angle_tj_RIGHT = pool_angle_tj;
test.pool_angle_start_RIGHT = pool_angle_start;

test.pool_turn_tj_RIGHT = pool_turn_tj;
test.pool_turn_start_RIGHT = pool_turn_start;

% calculate frequency of direction to
frequency_of_l = and(test.pool_angle_tj_LEFT<30, test.pool_angle_tj_LEFT>0);
frequency_of_r = and(test.pool_angle_tj_RIGHT<30, test.pool_angle_tj_RIGHT>0);

[h,p] =ztest(frequency_of_l,mean(frequency_of_r), std(frequency_of_r));

% calculate frequency of turns to
% Orientation towards T-junction
% "+" sign is clockwise (RIGH Turn) "-" sign is anti-clockwise (LEFT Turn)
% Orientation towards Startbox
% "-" sign is clockwise (RIGH Turn) "+" sign is anti-clockwise (LEFT Turn)
frequency_of_l = test.pool_turn_tj_LEFT>0;
frequency_of_r = test.pool_turn_tj_RIGHT>0;

[h,p] =ztest(frequency_of_l,mean(frequency_of_r), std(frequency_of_r));

% ---create categories out of mice---
temp_name =[];
for count_names = 1:numel(all_files)

    if ~isempty(b_mdl{1, count_names})
        temp_name(:,count_names) = str2num(all_files(count_names).name(3:4));
        %temp_name(:,count_names) = str2num(all_trials{1, count_names}.coord.name(3:4));
        trial_durations(count_names) = numel(b_mdl{1, count_names}.position_y);
    end
end
all_trial_durations = [all_trial_durations trial_durations];

% create a categorical variable
mice_cat.right =categorical(temp_name);

clusters.temp_right=cluster_temp_input;
clusters.mean_right=cluster_input;
clusters.trial_size_right= [ median_trials_size 1];

% ---resize clusters.temp_right to clusters.temp_left to use it in UMAP---
% pool behavioral variables into clusters
input = [];
input_cat = [];
temp_size = clusters.trial_size_right;

% find and remove NaN from Time of Turn clusters
clusters.mean_right( find(isnan(clusters.mean_right(:,9))) ,9) = median(clusters.mean_right(:,9),"omitnan");
clusters.mean_left( find(isnan(clusters.mean_left(:,9))) ,9) = median(clusters.mean_left(:,9),"omitnan");

clusters.temp_right( find(isnan(clusters.temp_right(:,end))) ,end) = median(clusters.mean_right(:,end),"omitnan");
clusters.temp_left( find(isnan(clusters.temp_left(:,end))) ,end) = median(clusters.mean_left(:,end),"omitnan");

clusters.temp_right( find(isnan(clusters.temp_right(:,end))) ,end) = median(clusters.mean_right(:,end),"omitnan");
clusters.temp_left( find(isnan(clusters.temp_left(:,end))) ,end) = median(clusters.mean_left(:,end),"omitnan");

temp_start = 0;
for count_behav = 1:10 % 9 behavioral vectors, defined as inputs# above

    % create vecotor
    temp_stop = temp_start+temp_size(count_behav);
    % size of clusters Trials x Time (time concatenated across all behav vars)
    input= resample(clusters.temp_right(:,temp_start+1:temp_stop)',clusters.trial_size_left(count_behav),temp_size(count_behav));
    input_cat = [input_cat input'];
    temp_start = temp_stop;
    
end
clusters.temp_right_resampled =input_cat;

input_to_cluster = [clusters.temp_left; clusters.temp_right_resampled];
%input_to_cluster = [clusters.mean_left(:,[2:4 7:10]); clusters.mean_right(:,[2:4 7:10])];

input_trial_type = [zeros(1,size(clusters.temp_left,1)) ones(1,size(clusters.temp_right,1))];
[reduction, umap, clusterIdentifiers, extras]=run_umap(input_to_cluster);
[idx,C] =kmeans(reduction,2);

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

%% Identify clusters and behaviors within

figure;
subplot(1,2,1)
gscatter(X(:,1),X(:,2),idx)
axis('square')
hold on
plot(C(:,1),C(:,2),'kx')
legend('Cluster 1','Cluster 2','Cluster 3','Cluster Centroid')

subplot(1,2,2)
silhouette(X,idx)


% pool behavioral variables into clusters
input_to_plot = cell(1,10);
temp_size = clusters.trial_size_left;
input_to_cluster = [clusters.temp_left; clusters.temp_right_resampled];

temp_start = 0;
for count_behav = 1:10 % 9 behavioral vectors, defined as inputs# above

    % create vecotor
    temp_stop = temp_start+temp_size(count_behav);
    input_to_plot{1,count_behav}.c1 = input_to_cluster(idx==1, temp_start+1:temp_stop);
    input_to_plot{1,count_behav}.c2 = input_to_cluster(idx==2, temp_start+1:temp_stop);
    input_to_plot{1,count_behav}.c3 = input_to_cluster(idx==3, temp_start+1:temp_stop);

    input_to_plot{1,count_behav}.all = input_to_cluster(:, temp_start+1:temp_stop);

    % partition into clusters by the trial type
    input_to_plot{1,count_behav}.cl = input_to_cluster(input_trial_type==0, temp_start+1:temp_stop);
    input_to_plot{1,count_behav}.cr = input_to_cluster(input_trial_type==1, temp_start+1:temp_stop);
    
    temp_start = temp_stop;
    
end

%% choose NOT resampled movement to compare across clusters

% --------------------------------------------------------------------------------------
% NOTE: the movement seems to be higher in Light ON trials in the plots
% below this cell is only because movement resampled to the median duration across trials 
% and then is averaged. This means that the major difference between Laser
% On - Off comes from the overall lower duration during Light On trials
% --------------------------------------------------------------------------------------
load('leftCorrect.mat')
all_trial_names= [];
all_trial_names= [all_trial_names; all_files];

% create categories out of mice
temp_name =[];
all_trial_movement = [];
n_leftCorrect = numel(all_files);
for count_names = 1:numel(all_files)

    if ~isempty(b_mdl{1, count_names})
        temp_name(:,count_names) = str2num(all_files(count_names).name(3:4));
        %temp_name(:,count_names) = str2num(all_trials{1, count_names}.coord.name(3:4));
        trial_durations(count_names) = numel(b_mdl{1, count_names}.position_y);  
        trial_movement(count_names) = mean(b_mdl{1, count_names}.speed.full_trial,'omitnan'); 
    end
end
all_trial_movement = [all_trial_movement trial_movement];

load('optoleftCorrect.mat')
for count_names = 1:numel(all_files)

    if ~isempty(b_mdl{1, count_names})
        temp_name(:,count_names) = str2num(all_files(count_names).name(3:4));
        %temp_name(:,count_names) = str2num(all_trials{1, count_names}.coord.name(3:4));
        trial_durations(count_names) = numel(b_mdl{1, count_names}.position_y);  
        trial_movement(count_names) = mean(b_mdl{1, count_names}.speed.full_trial,'omitnan'); 
    end
end
all_trial_movement = [all_trial_movement trial_movement];

subplot(2,1,1)
c1 = all_trial_movement(idx == 1);
c2 = all_trial_movement(idx == 2);
c3 = all_trial_movement(idx == 3);

barwitherr([std(c1, 'omitnan')/sqrt(numel(c1)) std(c2, 'omitnan')/sqrt(numel(c2)) ...
     ], ...
    [mean(c1, 'omitnan') mean(c2, 'omitnan')  ]);
set(gca,'XTickLabel',{'Cluster1', 'Cluster2', 'Cluster3', 'Cluster4'})
ylim([0 50])
p12 = ranksum(c1,c2);
p23 = ranksum(c3,c2);

title(['Movement Cluster 2 p_{12}= ' num2str(p12) ' p_{23}= ' num2str(p23)])

% plot by trial type
subplot(2,1,2)
c1 = all_trial_movement(1:n_leftCorrect);
c2 = all_trial_movement(n_leftCorrect+1:end);

barwitherr([std(c1, 'omitnan')/sqrt(numel(c1)) std(c2, 'omitnan')/sqrt(numel(c2)) ], ...
    [mean(c1, 'omitnan') mean(c2, 'omitnan')]);
set(gca,'XTickLabel',{'Cluster1', 'Cluster2'})
ylim([0 50])
p12 = ranksum(c1,c2);

title(['Movement Light On-Off p_{12}= ' num2str(p12) ])
%% plot averaged behavior in trial types

figure;

% remove zeros
input_to_plot{1,7}.cl(input_to_plot{1,7}.cl == 0) = [];
input_to_plot{1,7}.cr(input_to_plot{1,7}.cr == 0) = [];

input_to_plot{1,8}.cl(input_to_plot{1,8}.cl == 0) = [];
input_to_plot{1,8}.cr(input_to_plot{1,8}.cr == 0) = [];

input_to_plot{1,10}.cl(input_to_plot{1,10}.cl == 0) = [];
input_to_plot{1,10}.cr(input_to_plot{1,10}.cr == 0) = [];

cl = [];
cr = [];

current_lim =[0 500; 0 60; 0 70; 0 60; -1 1; -15 15; -0.3 0.3; -0.4 0.3; 0 10; 0 5];
for count_behav = 1:10

subplot(10,1,count_behav)

cl = median(input_to_plot{1,count_behav}.cl,2, 'omitnan') ;
cr = median(input_to_plot{1,count_behav}.cr,2, 'omitnan') ;

Names_of_behaviors ={'Maze Y position', 'Speed full trial' , 'Speed towards Start' , 'Speed towards Tjunction', ...
    'Angle Towards Start' ,'Angle Towards Tjunction', 'Direction of turn Start' ,'Direction of turn Tjunction', 'Time of Turn', 'Number of Turns'  };

barwitherr([std(cl, 'omitnan')/sqrt(numel(cl)) std(cr, 'omitnan')/sqrt(numel(cr)) ], ...
    [mean(cl, 'omitnan') mean(cr, 'omitnan')]);
set(gca,'XTickLabel',{'ClusterL', 'ClusterR'})
ylim(current_lim(count_behav,:));

p12 = ranksum(cl,cr);
title([Names_of_behaviors{1,count_behav} ' p_{12}= ' num2str(p12)])

end

%%
cl = [];
cr = [];

figure;
for count_behav = 1:6

subplot(6,1,count_behav)

cl.median = median(input_to_plot{1,count_behav}.cl,1) ;
cr.median = median(input_to_plot{1,count_behav}.cr,1) ;

cl.std = std(input_to_plot{1,count_behav}.cl,[],1)/sqrt(size(input_to_plot{1,count_behav}.cl,1)) ;
cr.std = std(input_to_plot{1,count_behav}.cr,[],1)/sqrt(size(input_to_plot{1,count_behav}.cr,1)) ; 

t = 0:1/all_trials{1, count_trials}.angle_fps: (numel(cl.median)-1)/all_trials{1, count_trials}.angle_fps;
shadedErrorBar(t,cl.median,cl.std,'b');
hold on
shadedErrorBar(t,cr.median,cr.std,'r');
hold off

%set(gca,'XTickLabel',{'Cluster1', 'Cluster2', 'Cluster3'})
title(Names_of_behaviors{1,count_behav})
end
%% plot averaged behavior in clusters

figure;

% remove zeros
input_to_plot{1,7}.c1(input_to_plot{1,7}.c1 == 0) = [];
input_to_plot{1,7}.c2(input_to_plot{1,7}.c2 == 0) = [];
input_to_plot{1,7}.c3(input_to_plot{1,7}.c3 == 0) = [];

input_to_plot{1,8}.c1(input_to_plot{1,8}.c1 == 0) = [];
input_to_plot{1,8}.c2(input_to_plot{1,8}.c2 == 0) = [];
input_to_plot{1,8}.c3(input_to_plot{1,8}.c3 == 0) = [];

input_to_plot{1,10}.c1(input_to_plot{1,10}.c1 == 0) = [];
input_to_plot{1,10}.c2(input_to_plot{1,10}.c2 == 0) = [];
input_to_plot{1,10}.c3(input_to_plot{1,10}.c3 == 0) = [];

c1 = [];
c2 = [];
c3 = [];
c4 = [];
current_lim =[0 500; 0 10; 0 70; 0 60; -1 1; -15 15; -0.7 0.7; -0.7 0.7; 0 10; 0 5];

Names_of_behaviors ={'Maze Y position', 'Speed full trial' , 'Speed towards Start' , 'Speed towards Tjunction', ...
    'Angle Towards Start' ,'Angle Towards Tjunction', 'Direction of turn Start' ,'Direction of turn Tjunction', 'Time of Turn', 'Number of Turns' };

for count_behav = 1:10

subplot(10,1,count_behav)

c1 = median(input_to_plot{1,count_behav}.c1,2, 'omitnan') ;
c2 = median(input_to_plot{1,count_behav}.c2,2, 'omitnan') ;
c3 = median(input_to_plot{1,count_behav}.c3,2, 'omitnan') ;
c4 = median(input_to_plot{1,count_behav}.c4,2, 'omitnan') ;

barwitherr([std(c1, 'omitnan')/sqrt(numel(c1)) std(c2, 'omitnan')/sqrt(numel(c2)) std(c3, 'omitnan')/sqrt(numel(c3)) std(c4, 'omitnan')/sqrt(numel(c4))], ...
    [mean(c1, 'omitnan') mean(c2, 'omitnan') mean(c3, 'omitnan') mean(c4, 'omitnan')]);
set(gca,'XTickLabel',{'Cluster1', 'Cluster2', 'Cluster3'})
ylim(current_lim(count_behav,:));

p12 = ranksum(c1,c2);
% p13 = ranksum(c1,c3);
% p23 = ranksum(c2,c3);
title([Names_of_behaviors{1,count_behav} ' p_{12}= ' num2str(p12) ' p_{13}= ' num2str(p13) ' p_{23}= ' num2str(p23)])


end

%% plot averaged behavior in clusters
c1 = [];
c2 = [];
c3 = [];

figure;
for count_behav = 1:6

subplot(6,1,count_behav)

c1.median = median(input_to_plot{1,count_behav}.c1,1) ;
c2.median = median(input_to_plot{1,count_behav}.c2,1) ;
c3.median = median(input_to_plot{1,count_behav}.c3,1) ;

c1.std = std(input_to_plot{1,count_behav}.c1,[],1)/sqrt(size(input_to_plot{1,count_behav}.c1,1)) ;
c2.std = std(input_to_plot{1,count_behav}.c2,[],1)/sqrt(size(input_to_plot{1,count_behav}.c2,1)) ; 
c3.std = std(input_to_plot{1,count_behav}.c3,[],1)/sqrt(size(input_to_plot{1,count_behav}.c3,1)) ; 

t = 0:1/all_trials{1, count_trials}.angle_fps: (numel(c1.median)-1)/all_trials{1, count_trials}.angle_fps;
shadedErrorBar(t,c1.median,c1.std,'b');
hold on
shadedErrorBar(t,c2.median,c2.std,'r');
hold on
shadedErrorBar(t,c3.median,c3.std,'g');
hold off

%set(gca,'XTickLabel',{'Cluster1', 'Cluster2', 'Cluster3'})
title(Names_of_behaviors{1,count_behav})
end