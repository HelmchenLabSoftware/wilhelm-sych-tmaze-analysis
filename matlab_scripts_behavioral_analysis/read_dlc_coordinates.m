clear;
% ---load behavioral data from DLC tracked csv files---
path_to_data = 'D:\Documents\data\wm_t_maze\photometry\video\LeftCorrect';
path_to_ca_data = 'D:\Documents\data\wm_t_maze\photometry\ca_data\LeftCorrect';

cd(path_to_data)
% Get a list of all files and folders in the mouse folder
all_files = dir('*.csv');

p_threshold = 0.99; % (corresponds to p= 0.05 )

% coordinate (0,0)
zero_x =  1;
zero_y =  50;

for count_trials = 1:numel(all_files)

    temp_file_name = all_files(count_trials).name;
    temp_data = importdata(temp_file_name);

    % save support variables, name etc
    data.name = all_files(count_trials).name;
    data.trial_length = numel(temp_data.data(:,1));

    if data.trial_length>1

    % read and filter by likelihood coordinates
    % ---------- nose ----------
    temp_pval = temp_data.data(:,4) > p_threshold;
    data.nose.x = temp_data.data(:,2);
    data.nose.x (temp_pval == 0 ,1) = NaN;
    data.nose.x_int = inpaint_nans(data.nose.x,3);

    data.nose.y = temp_data.data(:,3);
    data.nose.y (temp_pval == 0 ,1) = NaN;
    data.nose.y_int = inpaint_nans(data.nose.y,3);

%     plot(data.nose.y);
%     hold on
%     plot(data.nose.y_int+1);
%     hold off

    % ---------- leftEar ----------
    temp_pval = temp_data.data(:,7) > p_threshold;
    data.leftear.x = temp_data.data(:,5);
    data.leftear.x (temp_pval == 0 ,1) = NaN;
    data.leftear.x_int =  inpaint_nans(data.leftear.x,3);

    data.leftear.y = temp_data.data(:,6);
    data.leftear.y (temp_pval == 0 ,1) = NaN;
    data.leftear.y_int =  inpaint_nans(data.leftear.y,3);

    % ---------- rightEar ----------
    temp_pval = temp_data.data(:,10) > p_threshold;
    data.rightear.x = temp_data.data(:,8);
    data.rightear.x (temp_pval == 0 ,1) = NaN;
    data.rightear.x_int =  inpaint_nans(data.rightear.x,3);

    data.rightear.y = temp_data.data(:,9);
    data.rightear.y (temp_pval == 0 ,1) = NaN;
    data.rightear.y_int =  inpaint_nans(data.rightear.y,3);

    % ---------- tailBase ----------
    temp_pval = temp_data.data(:,13) > p_threshold;
    data.tail.x = temp_data.data(:,11);
    data.tail.x (temp_pval == 0 ,1) = NaN;
    data.tail.x_int = inpaint_nans(data.tail.x,3);

    data.tail.y = temp_data.data(:,12);
    data.tail.y (temp_pval == 0 ,1) = NaN;
    data.tail.y_int = inpaint_nans(data.tail.y,3);

    % ---------- laser Spot(for optogenetics)  ----------
    temp_pval = temp_data.data(:,16) > p_threshold;
    data.opto.x = temp_data.data(:,14);
    data.opto.x (temp_pval == 0 ,1) = NaN;
    data.opto.x_int = inpaint_nans(data.opto.x,3);

    data.opto.y = temp_data.data(:,15);
    data.opto.y (temp_pval == 0 ,1) = NaN;
    data.opto.y_int = inpaint_nans(data.opto.y,3);


    % ---- calculate behavioral parameters ----
    forward_x = mean([data.rightear.x_int data.leftear.x_int data.nose.x_int],2,'omitnan');
    forward_y = mean([data.rightear.y_int data.leftear.y_int data.nose.y_int],2,'omitnan');

    % bin data in steps of 1s and calculate difference
    forward_x_bin = [];
    forward_y_bin = [];

    angle = [];
    angle_sign = [];

    delta = 2;
    for count= 1 : floor(numel(forward_y)/delta)
        range = (count-1)*delta +1: count*delta;

        % forward coordinates binned
        fincrement_x = mean(forward_x(range),1,'omitnan');
        fincrement_y = mean(forward_y(range),1,'omitnan');

        % back/tale coordinates binned
        bincrement_x = mean(data.tail.x_int(range),1,'omitnan');
        bincrement_y = mean(data.tail.y_int(range),1,'omitnan');

        %calculate the angle
        % u relative to the T-maze
        %u = [zero_x  zero_y 0];
        % u relative to front-tale
        u = [fincrement_x-bincrement_x  zero_y 0];
        v = [fincrement_x-bincrement_x fincrement_y-bincrement_y 0];

        % check if mouse changes x direction
        count_direction = 1;
        count_direction((fincrement_x-bincrement_x)<0) = -1;
        angle.x_direction(count) = count_direction;

        % check if mouse changes y direction
        count_direction = 1;
        count_direction((fincrement_y-bincrement_y)<0) = -1;
        angle.y_direction(count) = count_direction;

        angle.rad(count) = atan2(norm(cross(u,v)),dot(u,v));
        angle.deg(count) = rad2deg(angle.rad(count));

        temp_angle = asin(cross(u,v)/(norm(u)*norm(v)));
        angle_sign.rad(count) = temp_angle(3);

% angle between twovectors from (0,0) to front and to tale- always small!
%         a = subspace([bincrement_x  bincrement_y]',[fincrement_x fincrement_y]');
%         angle.deg(count) = rad2deg(a);
%         angle.rad(count) = a;

    end

    % rotation
    all_trials{1, count_trials}.coord = data;
    all_trials{1, count_trials}.angle_sm = smooth(angle.deg,'moving');
    all_trials{1, count_trials}.angle_raw = angle.deg;
    all_trials{1, count_trials}.angle_rad = angle.rad;

    all_trials{1, count_trials}.angle_sign =rad2deg(smooth(angle_sign.rad,'moving'));

    all_trials{1, count_trials}.direction_x = angle.x_direction;
    all_trials{1, count_trials}.direction_y = angle.y_direction;


    all_trials{1, count_trials}.angle_fps = 30/delta;

% plot(angle.deg)
% hold on
% plot(test)
% hold off
    data = [];
    angle= [];
    angle_sign= [];

    end % if trial length is >1

end % count_trials

%% ------------------------------------
%  ----- Create Behavioral Model ------
%  ------------------------------------

b_mdl = cell(1,size(all_trials,2));
bm = [];
p_val=[];

pool_angle_start = [];
pool_angle_tj = [];
pool_turn_tj = [];
pool_turn_start = [];
pool_turnFreq_tj = [];
pool_turnFreq_start = [];

% use start sample to avoid beggining of the trial when the view on the
% mouse is incomplete
start_sample = 5;
stop_sample = 5;
gradient_range = 2;

for count_trials = 1:size(all_trials,2)

    % check if not empty
    if ~isempty(all_trials{1, count_trials})

    % --------------------- Y Position in the maze ------------------------

    temp_front= mean([all_trials{1, count_trials}.coord.nose.y_int ...
        all_trials{1, count_trials}.coord.leftear.y_int ...
        all_trials{1, count_trials}.coord.rightear.y_int],2,'omitnan');

    bm.position_y = temp_front;

    % ----------------------------- Speed -----------------------------
    temp_speed = smooth(abs(diff(temp_front))/(1/all_trials{1, count_trials}.angle_fps)) ;

    bm.speed.full_trial = temp_speed;

    % ind1 is the direction towards Startbox
    ind1= smooth(all_trials{1, count_trials}.direction_y) ==1 ;
    samples= ind1(start_sample:end);
    bm.speed.towards_start= temp_speed(samples);

    % ind2 is the direction towards the T-junction
    ind2= smooth(all_trials{1, count_trials}.direction_y) == -1;
    samples= ind2(start_sample:end);
    bm.speed.towards_tj= temp_speed(samples);

    % -----------------Test direction, angle, likelihood ----------------

    % angle in each quadrant
    % ind1 is Y the direction towards the startbox
    ind1= smooth(all_trials{1, count_trials}.direction_y) ==1 ;

    % ind3 is X the direction towards the startbox
    ind3= smooth(all_trials{1, count_trials}.direction_x) ;
    orientation(ind3==1) = 1;
    orientation(ind3==-1) = -1;

    % use same samples for angle orientation in X
    % and maze direction in Y
    samples= ind1(start_sample:end);
    bm.angle.towards_start= orientation(samples)'.*all_trials{1, count_trials}.angle_sm(samples);
    pool_angle_start = [pool_angle_start bm.angle.towards_start'];

    % ind2 is the direction towards the T-junction
    ind2= smooth(all_trials{1, count_trials}.direction_y) == -1;

    samples=ind2(start_sample:end);
    bm.angle.towards_tj= orientation(samples)'.*(all_trials{1, count_trials}.angle_sm(samples)-180);
    pool_angle_tj = [pool_angle_tj bm.angle.towards_tj'];

    if numel(bm.angle.towards_tj)>10
    % test Likelihood that mice are biased in a T-junction
    if or(isempty(bm.angle.towards_tj(bm.angle.towards_tj>0)) , ...
            isempty(bm.angle.towards_tj(bm.angle.towards_tj<0)))
        p_val.tj(count_trials) = nan;
        p_val.mean_tj(count_trials) = nan;
    else
        % do z test for proportions for binary variables
        % h = ztest(x,m,sigma)
        [h,p] =ztest(bm.angle.towards_tj<0,...
            mean(bm.angle.towards_tj>0),std(bm.angle.towards_tj>0));

        p_val.tj(count_trials) = p;
        % calculate mean angle
        p_val.mean_tj(count_trials) = mean(bm.angle.towards_tj,'omitnan');

    end % end if not empty vector
    end % end if more than 4 sample points are present

    if numel(bm.angle.towards_start)>10
    % test Likelihood that mice are biased in a start-box
    if or(isempty(bm.angle.towards_start(bm.angle.towards_start<0)) , ...
            isempty(bm.angle.towards_start(bm.angle.towards_start>0)))
        p_val.start(count_trials) = nan;
        p_val.mean_start(count_trials) = nan;
    else
        % do z test for proportions for binary variables
        % h = ztest(x,m,sigma)
        [h,p] =ztest(bm.angle.towards_start<0,...
            mean(bm.angle.towards_start>0),std(bm.angle.towards_start>0));
        p_val.start(count_trials) = p;
        % calculate mean angle
        p_val.mean_start(count_trials) = mean(bm.angle.towards_start,'omitnan');
    end % end if not empty vector
    end % end if more than 4 sample points are present

    % ---------------------------------------------------------------------
    %  --------------- Time of the turn towards t-junction ----------------
    
    % --- 1. Find the instances of the Turn from the Y-direction ---
    % input = smoothed vector containing Y direction change:
    % -1 towards T-junction; +1 towards start box
    input = round(smooth(all_trials{1, count_trials}.direction_y(start_sample:end-stop_sample),10));
    bm.time_of_turn =  find(input==-1, 1, 'first') / all_trials{1, count_trials}.angle_fps;

    % ---------------------------------------------------------------------
    %  ----------------- Number and direction of turns --------------------
    ind_direction_tj =[];
    ind_direction_start =[];

    % --- 1. Find the instances of the Turn from the Y-direction ---
    for count= 2:numel(input)
        if and(input(count)==-1, input(count-1)>-1)
            ind_direction_tj= [ind_direction_tj count];
        elseif and(input(count)==1, input(count-1)<1)
            ind_direction_start= [ind_direction_start count];
        end
    end

    % input = smoothed vector containing X direction change:
    input2 = all_trials{1, count_trials}.angle_sign(start_sample:end);

    if isempty(ind_direction_tj)
        bm.direction_of_turn.tj = [];
    else

    % find gradients/direction for ind_direction_tj
    for count = 1:numel(ind_direction_tj)

        % find patch of X to do gradient
        patch_X= (input2(ind_direction_tj(count)-gradient_range:ind_direction_tj(count)+gradient_range));

        % use regression on patch to identify the sign of the slope

%         mdl = fitlm(1:numel(patch_X),patch_X);
%         coefs = mdl.Coefficients.Estimate;
%         p = mdl.Coefficients.pValue;
%         if p(2,1)<0.05
%             % on the video -coef is the counter-clockwise turn 
%             % on the video +coef is the clockwise turn 
%             bm.direction_of_turn.tj(count) = coefs(2,1);
%         else
%             bm.direction_of_turn.tj(count) = nan;
%         end

    % use sign of the mean angle directly
    bm.direction_of_turn.tj(count) = sign(mean(patch_X));

    end

    end % if isempty ind

    pool_turn_tj = [pool_turn_tj bm.direction_of_turn.tj];
    pool_turnFreq_tj = [pool_turnFreq_tj numel(bm.direction_of_turn.tj)];


    % find gradients/direction for ind_direction_start
    if isempty(ind_direction_start)
        bm.direction_of_turn.start = [];
    else

    for count = 1:numel(ind_direction_start)

        % find patch of X to do gradient
        patch_X= (input2(ind_direction_start(count)-gradient_range:ind_direction_start(count)+gradient_range));

        % use sign of the mean angle directly
        bm.direction_of_turn.start(count) = sign(mean(patch_X));

%         mdl = fitlm(1:numel(patch_X),patch_X);
%         coefs = mdl.Coefficients.Estimate;
%         p_val = mdl.Coefficients.pValue;
%         if p_val(2,1)<0.05
%             bm.direction_of_turn.start(count) = coefs(2,1);
%         else
%             bm.direction_of_turn.start(count) = [];
%         end
    end

    end % if isempty ind

    pool_turn_start = [pool_turn_start bm.direction_of_turn.start];
    pool_turnFreq_start = [pool_turnFreq_start numel(bm.direction_of_turn.start)];

    b_mdl{1,count_trials} = bm;
    bm = [];

    end % end if empty

end % count_trials

%% -------------------------------------------------------------------
%  ----- plot pooled behavioral data with Orientation and Turns ------
%  -------------------------------------------------------------------

subplot(3,2,1)
polarhistogram(deg2rad(pool_angle_tj),24)
thetalim([-90 90])
%p= ranksum(pool_angle_tj(pool_angle_tj<0), pool_angle_tj(pool_angle_tj>0));
title(['Angle towards T-junction'])

subplot(3,2,2)
polarhistogram(deg2rad(pool_angle_start),24)
thetalim([-90 90])

%test the angle in bins of [0,+30] and [-30,0] degrees
% input1 = and(pool_angle_start>0, pool_angle_start<30);
% input2 = and(pool_angle_start<0, pool_angle_start>-30);
% [h,p] =ztest(input1, mean(input2), std(input2));
title('Angle towards Startbox')

subplot(3,2,3)
% Orientation towards T-junction
% "+" sign is clockwise (RIGH Turn) "-" sign is anti-clockwise (LEFT Turn)
barwitherr([std(pool_turn_tj<0)/sqrt(size(all_trials,2)) std(pool_turn_tj>0)/sqrt(size(all_trials,2))],...
    [mean(pool_turn_tj<0) mean(pool_turn_tj>0)])
set(gca,'XTickLabel',{'LEFT','RIGH'});

[h,p] =ztest(pool_turn_tj<0, mean(pool_turn_tj>0), std(pool_turn_tj>0));
title(['Turns to T-Junction ' num2str(p)])

subplot(3,2,4)

% Orientation towards Startbox
% "-" sign is clockwise (RIGH Turn) "+" sign is anti-clockwise (LEFT Turn)
barwitherr([std(pool_turn_start<0)/sqrt(size(all_trials,2)) std(pool_turn_start>0)/sqrt(size(all_trials,2))],...
    [mean(pool_turn_start<0) mean(pool_turn_start>0)])
set(gca,'XTickLabel',{'RIGH','LEFT'});

[h,p] =ztest(pool_turn_start<0, mean(pool_turn_start>0), std(pool_turn_start>0));
title(['Turns to Startbox ' num2str(p)])

subplot(3,2,5)
% Orientation towards T-junction
% "+" sign is clockwise (RIGH Turn) "-" sign is anti-clockwise (LEFT Turn)
barwitherr([0 0], [sum(pool_turnFreq_tj==0)/size(all_trials,2) sum(pool_turnFreq_tj~=0)/size(all_trials,2)])

set(gca,'XTickLabel',{'undef','def'});
ylim([0 1])
title(['Frequency of Turns to T-Junction '])

subplot(3,2,6)

% Orientation towards Startbox
% "-" sign is clockwise (RIGH Turn) "+" sign is anti-clockwise (LEFT Turn)
barwitherr([0 0], [sum(pool_turnFreq_start==0)/size(all_trials,2) sum(pool_turnFreq_start~=0)/size(all_trials,2)])
set(gca,'XTickLabel',{'undef','def'});
ylim([0 1])
title(['Frequency of Turns to Startbox '])

% scatter plot to test if within each trial there are biases towards certain direction
% scatter(p_val.mean_tj,p_val.tj);
% set(gca,'yscale','log')

%% ------------------------------------------------------------------
% --- do UMAP clustering on all trials based on Behavioral model ---
% -------------------------------------------------------------------
cluster_input =[];
trials_size =[];

% -------------------------------------------------------------------------
% ---- Create Cluster Input of Behavioral Variables with Temporal Info ----
% -------------------------------------------------------------------------

% ----- 1. Find median trial size across all data -----
for count_trials = 1:size(b_mdl,2)

    temp = b_mdl{1,count_trials};
    if ~isempty(temp) 

        size_vector= [size(temp.position_y,1) size(temp.speed.full_trial,1) ...
                size(temp.speed.towards_start,1) size(temp.speed.towards_tj,1) ...
                size(temp.angle.towards_start,1) size(temp.angle.towards_tj,1) ...
                size(temp.direction_of_turn.start,2) size(temp.direction_of_turn.tj,2) ...
                1];
    end % if full trial is empty

trials_size =cat(1,trials_size, size_vector);

end % count_trials

max_trials_size = max(trials_size,[],1);
median_trials_size = round(median(trials_size,1));

% ----- 2. Resample all to median -----
cluster_temp_input =[];
cluster_input =[];

backup_mdl_speed_tj = [];
backup_mdl_angle_tj = [];
backup_mdl_time_of_turn = [];

for count_trials = 1:size(b_mdl,2)

    temp = b_mdl{1,count_trials};
    if ~isempty(temp)

    if ~isempty(temp.speed.towards_tj)

    input1 = resample(temp.position_y', median_trials_size(1), numel(temp.position_y));
    input2 = resample(temp.speed.full_trial', median_trials_size(2), numel(temp.speed.full_trial));
    input3 = resample(temp.speed.towards_start', median_trials_size(3), numel(temp.speed.towards_start));
    input4 = resample(temp.speed.towards_tj', median_trials_size(4), numel(temp.speed.towards_tj));
    input5 = resample(temp.angle.towards_start', median_trials_size(5), numel(temp.angle.towards_start));
    input6 = resample(temp.angle.towards_tj', median_trials_size(6), numel(temp.angle.towards_tj));
    input7 = mean(temp.direction_of_turn.start);
    input7(isnan(input7)) = 0;
    input8 = mean(temp.direction_of_turn.tj);
    input8(isnan(input8)) = 0;

    input9 = temp.time_of_turn;
    input9(isempty(input9)) = median(backup_mdl_time_of_turn);

    else

    input1 = resample(temp.position_y', median_trials_size(1), numel(temp.position_y));
    input2 = resample(temp.speed.full_trial', median_trials_size(2), numel(temp.speed.full_trial));
    input3 = resample(temp.speed.towards_start', median_trials_size(3), numel(temp.speed.towards_start));
    input4 = median(backup_mdl_speed_tj);
    input5 = resample(temp.angle.towards_start', median_trials_size(5), numel(temp.angle.towards_start));
    input6 = median(backup_mdl_angle_tj);
    input7 = mean(temp.direction_of_turn.start);
    input7(isnan(input7)) = 0;
    input8 = mean(temp.direction_of_turn.tj);
    input8(isnan(input8)) = 0;

    input9 = temp.time_of_turn;
    input9(isempty(input9)) = median(backup_mdl_time_of_turn);

    end % if the vector to T-junction is not empty

    temp_measures=[input1 input2 input3 input4 input5 input6 input7 input8 input9];
    % use measures for reduced behavioral vector - no temporal information
    measures=[mean(input1) mean(input2) mean(input3) mean(input4) mean(input5) mean(input6) input7 input8 input9];

    % use backup models containing data from past trials to approximate
    % with median current trial is the model variables are NaN
    backup_mdl_speed_tj = cat(1,backup_mdl_speed_tj,input4);
    backup_mdl_angle_tj = cat(1,backup_mdl_angle_tj,input6);
    backup_mdl_time_of_turn = cat(1,backup_mdl_time_of_turn,input9);

    end % if full trial is empty

cluster_temp_input =cat(1,cluster_temp_input,temp_measures);
cluster_input =cat(1,cluster_input,measures);

end % count_trials

% ----- 3. Do Dimentionality Reduction -----

% do tSNE on conctenated behavioral vectors with temporal structure
% [Y, loss] = tsne(cluster_temp_input,'Algorithm','exact','Distance','cosine');
% gscatter(Y(:,1),Y(:,2));

[reduction, umap, clusterIdentifiers, extras]=run_umap(cluster_temp_input);
[idx,C] =kmeans(reduction,3);

% idx(idx == 2) =0;
% idx(idx == 1) =2;
% idx(idx == 2) =1;

X = reduction;

figure;
subplot(1,2,1)
gscatter(X(:,1),X(:,2),idx)
axis('square')
hold on
plot(C(:,1),C(:,2),'kx')
legend('Cluster 1','Cluster 2','Cluster 3','Cluster Centroid')

subplot(1,2,2)
silhouette(X,idx)

% % do tSNE on means of behavioral vectors
% [Y, loss] = tsne(cluster_input,'Algorithm','exact','Distance','cosine');
% gscatter(Y(:,1),Y(:,2));




% pool behavioral variables into clusters
input_to_plot = cell(1,9);
temp_size = [median_trials_size(1:6) 1 1 1];

temp_start = 0;
for count_behav = 1:9 % 9 behavioral vectors, defined as inputs# above

    % create vecotor
    temp_stop = temp_start+temp_size(count_behav);
    input_to_plot{1,count_behav}.c1 = cluster_temp_input(idx==1, temp_start+1:temp_stop);
    input_to_plot{1,count_behav}.c2 = cluster_temp_input(idx==2, temp_start+1:temp_stop);
    input_to_plot{1,count_behav}.c3 = cluster_temp_input(idx==3, temp_start+1:temp_stop);
    input_to_plot{1,count_behav}.all = cluster_temp_input(:, temp_start+1:temp_stop);
    temp_start = temp_stop;
    
end

%% plot averaged behavior in clusters

figure;

% remove zeros
input_to_plot{1,7}.c1(input_to_plot{1,7}.c1 == 0) = NaN;
input_to_plot{1,7}.c2(input_to_plot{1,7}.c2 == 0) = NaN;
input_to_plot{1,7}.c3(input_to_plot{1,7}.c3 == 0) = NaN;

input_to_plot{1,8}.c1(input_to_plot{1,8}.c1 == 0) = NaN;
input_to_plot{1,8}.c2(input_to_plot{1,8}.c2 == 0) = NaN;
input_to_plot{1,8}.c3(input_to_plot{1,8}.c3 == 0) = NaN;

c1 = [];
c2 = [];
c3 = [];
current_lim =[0 500; 0 60; 0 70; 0 60; -1 1; -15 15; -0.3 0.3; -0.4 0.3; 0 8];
for count_behav = 1:9

subplot(9,1,count_behav)

c1 = median(input_to_plot{1,count_behav}.c1,2, 'omitnan') ;
c2 = median(input_to_plot{1,count_behav}.c2,2, 'omitnan') ;
c3 = median(input_to_plot{1,count_behav}.c3,2, 'omitnan') ;

Names_of_behaviors ={'Maze Y position', 'Speed full trial' , 'Speed towards Start' , 'Speed towards Tjunction', ...
    'Angle Towards Start' ,'Angle Towards Tjunction', 'Direction of turn Start' ,'Direction of turn Tjunction', 'Time of Turn' };

barwitherr([std(c1, 'omitnan')/sqrt(numel(c1)) std(c2, 'omitnan')/sqrt(numel(c2)) std(c3, 'omitnan')/sqrt(numel(c3))], ...
    [mean(c1, 'omitnan') mean(c2, 'omitnan') mean(c3, 'omitnan')]);
set(gca,'XTickLabel',{'Cluster1', 'Cluster2', 'Cluster3'})
ylim(current_lim(count_behav,:));

p12 = ranksum(c1,c2);
p13 = ranksum(c1,c3);
p23 = ranksum(c2,c3);
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

%% ----- check how Ca dynamics changes across clusters -----

cd(path_to_ca_data)
% Get a list of all files and folders in the mouse folder
all_cafiles = dir('*.mat');
count_all_trials = 1;
count_ind = 1;
cadata = [];

for count_sessions = 1:numel(all_cafiles)

    temp_file_name = all_cafiles(count_sessions).name;
    temp_data = importdata(temp_file_name);

    % save support variables, name etc
    cadata.name = all_cafiles(count_sessions).name;

    for count_trials = 1:size(temp_data,2)

        cadata.trial_names{1,count_all_trials}= [cadata.name(1:end-4) '_trial_' num2str(count_trials) ...
            'DLC_mobnet_100_WM_T-maze2Dec2shuffle1_120000.csv'];
        
        
        % find a matching name in all_files
        for count_names= 1: numel(all_files)
            %strfind(all_files(:).name, cadata.trial_names{1,count_all_trials});

            if matches( all_files(count_names).name, cadata.trial_names{1,count_all_trials})
                index_name_to_behavior = count_names;
            end
        end % count_names

        cadata.trial_ind(count_all_trials) = index_name_to_behavior;
        % alternative for the bar plot
        %cadata.trial_ca(count_all_trials,:) = mean(temp_data{1, count_trials},2,'omitnan');  

        cadata.trial_ca(count_all_trials,:) = resample(temp_data{1, count_trials},500,numel(temp_data{1, count_trials})); 
        
        count_all_trials = count_all_trials+1;

    end % count each trial within session

    temp_file_name = [];
    temp_data =[];

end

%%
% align ca data to trials as in behavioral pooling
input_ca.all = cadata.trial_ca(cadata.trial_ind,:);

% split into behavioral clusters
input_ca.c1 = input_ca.all(idx==1,:);
input_ca.c2 = input_ca.all(idx==2,:);
input_ca.c3 = input_ca.all(idx==3,:);


c1.median = mean(input_ca.c1,1) ;
c2.median = mean(input_ca.c2,1) ;
c3.median = mean(input_ca.c3,1) ;

c1.std = std(input_ca.c1,[],1)/sqrt(size(input_ca.c1,1)) ;
c2.std = std(input_ca.c2,[],1)/sqrt(size(input_ca.c2,1)) ; 
c3.std = std(input_ca.c3,[],1)/sqrt(size(input_ca.c3,1)) ; 

subplot(3,1,1)
t = 0:499;
shadedErrorBar(t,c1.median,c1.std,'b');
hold on
shadedErrorBar(t,c2.median,c2.std,'r');
hold on
shadedErrorBar(t,c3.median,c3.std,'g');
hold off

subplot(3,1,2)
for count = 1:500
    p = ranksum(input_ca.c1(:,count),input_ca.c2(:,count));
    p_ts(count) = p;
end

p_ts_fdr = fdr_bh(p_ts);
%fdr_bh(pval)
%h=fdr_bh(pval(~isnan(pval/10)))

imagesc(p_ts_fdr>0)
map =[1 1 1
    0 0 0]
colormap(map)

subplot(3,1,3)
for count_trials = 1:numel(cadata.trial_ind)

    if ~isempty(b_mdl{1, count_trials})
    percent_of_tj(count_trials)= 100- 100*numel(b_mdl{1, count_trials}.speed.towards_tj)...
        /numel(b_mdl{1, count_trials}.speed.full_trial);
    end
end
histogram(percent_of_tj)
%%
input_ca.c1 = mean(input_ca.all(idx==1,:),2);
input_ca.c2 = mean(input_ca.all(idx==3,:),2);
input_ca.c3 = mean(input_ca.all(idx==2,:),2);

barwitherr([std(input_ca.c1, 'omitnan')/sqrt(numel(input_ca.c1)) std(input_ca.c2, 'omitnan')/sqrt(numel(input_ca.c2)) ...
    std(input_ca.c3, 'omitnan')/sqrt(numel(input_ca.c3))], ...
    [mean(input_ca.c1, 'omitnan') mean(input_ca.c2, 'omitnan') mean(input_ca.c3, 'omitnan')]);
set(gca,'XTickLabel',{'Cluster1', 'Cluster2', 'Cluster3'})
[p12,h] = ranksum(input_ca.c1,input_ca.c2);
[p13,h] = ranksum(input_ca.c1,input_ca.c3);
[p23,h] = ranksum(input_ca.c3,input_ca.c2);
title([' p12 =' num2str(p12) ' p23 =' num2str(p23) ' p13 =' num2str(p13)])

% save workspace to compare left and right trials
%save('D:\Documents\data\wm_t_maze\photometry\ca_data\rightMistake.mat')
