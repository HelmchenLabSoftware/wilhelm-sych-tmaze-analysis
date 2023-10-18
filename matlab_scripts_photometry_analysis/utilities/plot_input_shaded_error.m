function [smoothed_MeanTrial, smoothed_sem_MeanTrial] = plot_input_shaded_error (input_to_plot, color_to_plot, t, y_axis1,y_axis2)


transparent = 0;
MeanTrial=mean(input_to_plot,1);
sem_MeanTrial=std(input_to_plot,1)/sqrt(size(input_to_plot,1));

smoothed_MeanTrial = movmean(MeanTrial,5);
smoothed_sem_MeanTrial = movmean(sem_MeanTrial,5);

H = NewShadedErrorBar(t,movmean(MeanTrial,5),movmean(sem_MeanTrial,5),'k',transparent,21,y_axis1,y_axis2,-1,max(t),'','',size(input_to_plot,1),color_to_plot) ;

end