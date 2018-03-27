function Speed_PDF_Object = compute_speed_dist_allColonies(colony_2_include,Pre_Post_Indctr)
%% 
% Inputs:
%    ---colony_2_include---
%    An array of Colony-IDs
%    Size: m x 1
%          where m is the number of colonies being combined to
%          compute the estimates of speed distributions
%    ---Pre_Post_Indctr---
%    A binary 0/1 variable
%    0: Consider Pre-Exposure Data
%    1: Consider Post-Exposure Data
%
% Outputs:
%    ---Speed_PDF_Object---
%    An object defining the PDF of the log10 values of speed
%    * Use "random(Speed_PDF_Object,*,*)" to generate samples from this
%    distribution.
%
data_FrameRate = 2; % in fps
Speed_across_Colonies = [];
%
%%
%
for colonyNumber = colony_2_include
    file2load = strcat('data/allDataCol',num2str(colonyNumber),'.mat');
    %
    if (Pre_Post_Indctr == 0)
        load(file2load,'preNest');
        nestData = preNest;
    elseif  (Pre_Post_Indctr == 1)
        load(file2load,'postNest');
        nestData = postNest;
    end
    %
    preDiffVel = abs(diff(nestData(:,:,1:2)));
    preVels = sqrt(preDiffVel(:,:,1).^2 + preDiffVel(:,:,2).^2);
    % Correction for frame rate
    preVels = preVels*data_FrameRate; 
    speedMat = preVels;
    %
    Speed_across_Colonies = [Speed_across_Colonies; speedMat(:)];
    %
end
%
Speed_across_Colonies(isnan(Speed_across_Colonies)==1) = [];
%
Speed_across_Colonies__LOG = log10(Speed_across_Colonies);
Speed_across_Colonies__LOG(Speed_across_Colonies__LOG<-4) = [];
%
% figure('units','normalized','outerposition',[0 0 1 1])
histfit(Speed_across_Colonies__LOG,200,'kernel');
title('Distribution of the log of Speed (across Colonies)');
saveas(gcf,'Distr_Speed_log10_acrss_Colonies.png');
%
Speed_PDF_Object = fitdist(Speed_across_Colonies__LOG,'kernel');