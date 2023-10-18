%% plot mice turns
clear
%cd 'D:\Documents\data\wm_t_maze\opto_inhibition\video\behavior_pool_t-junction'
cd 'D:\Documents\data\wm_t_maze\photometry\video\behavior_pool_start-box'
input_Left = load('mice_turns_angle_LeftCorrect') ;
input_Right = load('mice_turns_angle_RightCorrect') ;

% for PHOTOMETRY comparison mice_turns_angle_LeftCorrect vs mice_turns_angle_optoLeftMistake
unique_mice_cat = [7 8 9 11 17 41];
% for OPTO comparison mice_turns_angle_LeftMistake vs mice_turns_angle_optoLeftMistake
%unique_mice_cat = [7 8 9 11 17 41 42];

figure(1)
for count_mice= 1:numel(unique_mice_cat)

    input1 = input_Left.mice_turns{1,count_mice};
    input2 = input_Right.mice_turns{1,count_mice};

    subplot(numel(unique_mice_cat), 1, count_mice)
    barwitherr([std(input1,[],'omitnan')/sqrt(numel(input1)) ...
        std(input2,[],'omitnan')/sqrt(numel(input2)) ], ...
        [ mean(input1,'omitnan') mean(input2,'omitnan') ]);

    p = ranksum(input1,input2);
    title(['Turns Left vs Right p ' num2str(p)])
    ylim([-0.30 0.30])

end

figure(2)
for count_mice= 1:numel(unique_mice_cat)

     input1 = input_Left.mice_angle{1,count_mice};
    input2 = input_Right.mice_angle{1,count_mice};

    subplot(numel(unique_mice_cat), 1, count_mice)
    barwitherr([std(input1,[],'omitnan')/sqrt(numel(input1)) ...
        std(input2,[],'omitnan')/sqrt(numel(input2)) ], ...
        [ mean(input1,'omitnan') mean(input2,'omitnan') ]);

    p = ranksum(input1,input2);
    title(['Angle Left vs Right p ' num2str(p)])
    ylim([-30 30])

end    

%% plot pooled data and mean across mice on one plot

input1_to_plot = [];
input2_to_plot = [];

input1_mean = [];
input2_mean = [];

% collect trial averaged and mouse averaged data separately
for count_mice= 1:numel(unique_mice_cat)

    input1 = input_Left.mice_turns{1,count_mice};
    input2 = input_Right.mice_turns{1,count_mice};

    [p,h,stats] = ranksum(input1,input2);
    p_val(count_mice) = p;
    %z_val(count_mice) = stats.zval;

    input1_to_plot = [input1_to_plot input1];
    input2_to_plot = [input2_to_plot input2];

    input1_mean = [input1_mean mean(input1,'omitnan')];
    input2_mean = [input2_mean mean(input2,'omitnan')];

end

subplot(2, 3, 1)
    barwitherr([std(input1_to_plot,[],'omitnan')/sqrt(numel(input1_to_plot)) ...
        std(input2_to_plot,[],'omitnan')/sqrt(numel(input2_to_plot)) ], ...
        [ mean(input1_to_plot,'omitnan') mean(input2_to_plot,'omitnan') ]);

hold on

    scale = 0.1;
    scatter([1*ones(1,numel(input1_mean))+ scale*rand(1,numel(input1_mean)) 2*ones(1, numel(input2_mean)) + scale*rand(1,numel(input1_mean))], ...
        [input1_mean input2_mean] )

hold on

    plot([ones(1,numel(input1_mean)); 2*ones(1,numel(input1_mean))] ,[input1_mean; input2_mean] ,'LineWidth',2,'Color',[0.7 0.7 0.7]);
    title(['Turns Left vs Right '])
    ylim([-0.90 0.50])
hold off

subplot(2, 3, 3)
scatter(ones(1,numel(input1_mean)) , log(p_val))
hold on
hline(log(0.05))
hold off

% ------------------------
% ------ plot angles -----
% ------------------------

input1_to_plot = [];
input2_to_plot = [];

input1_mean = [];
input2_mean = [];

% collect trial averaged and mouse averaged data separately
for count_mice= 1:numel(unique_mice_cat)

    input1 = input_Left.mice_angle{1,count_mice};
    input2 = input_Right.mice_angle{1,count_mice};

    [p,h,stats] = ranksum(input1,input2);
    p_val(count_mice) = p;
    %z_val(count_mice) = stats.zval;

    input1_to_plot = [input1_to_plot input1];
    input2_to_plot = [input2_to_plot input2];

    input1_mean = [input1_mean mean(input1,'omitnan')];
    input2_mean = [input2_mean mean(input2,'omitnan')];

end

[p,h,stats] = ranksum(input1_to_plot,input2_to_plot)

subplot(2, 3, 4)

polarhistogram(deg2rad(input1_to_plot),25,'FaceColor','blue','FaceAlpha',.3);
hold on
polarhistogram(deg2rad(input2_to_plot),25,'FaceColor','red','FaceAlpha',.3);
hold on
thetalim([-90 90])


subplot(2, 3, 5)
polarhistogram(deg2rad(input1_mean),25,'FaceColor','blue','FaceAlpha',.3);
hold on
polarhistogram(deg2rad(input2_mean),25,'FaceColor','red','FaceAlpha',.3);
hold on
thetalim([-90 90])


subplot(2, 3, 6)
scatter(ones(1,numel(input1_mean)) , log(p_val))
hold on
hline(log(0.05))
hold off
