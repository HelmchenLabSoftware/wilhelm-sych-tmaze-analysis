function [AUC, AUC_shuffle,...
          h,adj_p,...
          smoothed_mean_trial,...
          smoothed_sem_mean_trial] = compare_vectors_ROC (input_1,...
                                                            input_2r,...
                                                            time,...
                                                            number_of_draws,...
                                                            number_of_trials_sampled,...
                                                            name_figure)

predictor = [];
label = [];
label_shuffle = [];
AUC = [];
AUC_shuffle = [];

[n_trials1,length1] = size(input_1);

input_2 = [];

[n_trials2,length2] = size(input_2r);

% reshape test_1 and test_2 to a common length
for i = 1:n_trials2
    input_2(i,:) = resample(input_2r(i,:),length1,length2);
end


% --- Input EVENT name into the structure field as a predictor---
predictor = [];
label = [];
label_shuffle = [];
AUC = [];
AUC_shuffle = [];


% --- bin time series =1 if Yes else =0? ---
bin_ON = 0;

if bin_ON == 1
    [n,length1] = size(input_1);
    test_1 = [];
    for i = 1:n
        test_1(i,:) = bin_time_series(input_1(i,:), bin_size);
    end
    
    [n,length2] = size(input_2);
    test_2 = [];
    for i = 1:n
        test_2(i,:) = bin_time_series(input_2(i,:), bin_size);
    end
else
    test_1 = input_1;
    test_2 = input_2;
end

[n1,all_steps] = size(test_1);
[n2,all_steps] = size(test_2);

label = [];
label_shuffle = [];
label=[ones(n1,1); zeros(n2,1)];
% labels are shuffled after random sampling
%label_shuffle= label(randperm(length(label)));


% use interval of samples to predict trial types
start_interval = 0;
stop_interval = 0;

for count_draws = 1:number_of_draws
    
    for count_time_step = start_interval+1:all_steps-stop_interval
        
        % use a vector of 3 samples as a predictor (not one time step)
        pdt_interval = count_time_step-start_interval:count_time_step+stop_interval;
        predictor = [test_1(:,pdt_interval); test_2(:,pdt_interval)];
        [predictor_sampled, index_sampled] = datasample(predictor,number_of_trials_sampled,'Replace',false);
        label_sampled = label(index_sampled);
        label_shuffle = label_sampled(randperm(length(label_sampled)));
        % fit General Linear Model
        mdl = fitglm(predictor_sampled,label_sampled,'Distribution','binomial','Link','logit');
        score_log = mdl.Fitted.Probability; % Probability estimates
        % fit General Linear Model for shuffled labels
        mdl_shuffle = fitglm(predictor_sampled,label_shuffle,'Distribution','binomial','Link','logit');
        score_shuffle = mdl_shuffle.Fitted.Probability; % Probability estimates
        
        % Compute the standard ROC curve using the probabilities for scores.
        [Xlog,Ylog,Tlog,AUClog] = perfcurve(label_sampled,score_log,1);
        [Xshf,Yshf,Tshf,AUCshf] = perfcurve(label_shuffle,score_shuffle,1);
        
        AUC(count_time_step,count_draws)=AUClog;
        AUC_shuffle(count_time_step,count_draws)=AUCshf;
        
    end
end


color_grey = [128,128,128]/255;

figure

[time_bin,~] = size(AUC);

subplot(2,1,1)


y_axis1 = 0.5;
y_axis2 = 0.65;

[smoothed_mean_trial.AUC_shuffle, smoothed_sem_mean_trial.AUC_shuffle] = plot_input_shaded_error (AUC_shuffle', color_grey, time, y_axis1,y_axis2);
hold on
[smoothed_mean_trial.AUC, smoothed_sem_mean_trial.AUC] = plot_input_shaded_error (AUC', 'r', time, y_axis1,y_axis2);

axis ('square')
title(['AUC_',name_figure],'Interpreter', 'none');
hold off


p_val=[];

for count_time_step= start_interval+1:all_steps-stop_interval
    
    p = ranksum(AUC(count_time_step,:), AUC_shuffle(count_time_step,:));
    p_val(count_time_step)=p;
    
end

%
p_value_threshold = 0.01;
[h, crit_p, adj_ci_cvrg, adj_p] = false_discovery_rate_control(p_val,p_value_threshold,'pdep','yes');

subplot(2,1,2)

Z=zeros(1,length(time));
for ii=1:length(time)
    plot ([time(ii) time(ii)],[Z(ii) h(ii)],'k')
    hold on
end
axis([-1 max(time) 0 1])
axis ('square')
title (['p_adj < ',num2str(p_value_threshold)],'Interpreter', 'none')
set (gcf,'Position', [442   165   786   808]);
set(findall(gcf,'-property','FontSize'),'FontSize',20)
save_plots(['AUC_',name_figure])
end