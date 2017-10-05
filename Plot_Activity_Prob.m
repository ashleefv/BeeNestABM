clear; clc; close all
%
% Change the COLONY ID to select which colony data you want to process
Colony_ID = 4;
%
% Loads the relevant data file
% The data file should be inside "..\data" folder
% This relative location should be maintained for proper functioning
file2load = strcat('..\data\Essential_Info_Col_',num2str(Colony_ID),'.mat');
load(file2load);
%
%
figure('units','normalized','outerposition',[0 0 1 1])
%
subplot(1,2,1)
bar1 = bar(Activity_Prob_Dist_Pre(:,3:5)');
ylim([0,1]);
legend('Active','Inactive');
set(gca,'XTickLabel',{'Group-I','Group-II','Group-III'},'fontsize',18);
title_TEXT = strcat('Pre-Exposure (Colony ID: ',num2str(Colony_ID),')');
title(title_TEXT);
%
subplot(1,2,2)
bar2 = bar(Activity_Prob_Dist_Post(:,3:5)');
ylim([0,1]);
legend('Inactive','Active');
set(gca,'XTickLabel',{'Group-I','Group-II','Group-III'},'fontsize',18);
title_TEXT = strcat('Pre-Exposure (Colony ID: ',num2str(Colony_ID),')');
title(title_TEXT);
%
image_2_save = strcat('..\results\Activity_Level_Col_',num2str(Colony_ID),'.png');
saveas(gcf,image_2_save)