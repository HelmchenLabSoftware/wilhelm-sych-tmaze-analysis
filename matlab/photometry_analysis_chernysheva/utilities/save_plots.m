%just to save Plots and Figures in different formats

function save_plots(name_plot)

print(gcf, name_plot,'-dpdf') 
print(gcf, name_plot,'-dpng') 
savefig([name_plot,'.fig']) %saves figure as matlab file
% saveas(gcf,name_plot,'epsc') %saves figure as epc for Affinity
end