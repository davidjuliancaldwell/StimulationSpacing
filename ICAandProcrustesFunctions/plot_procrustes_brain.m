function [] = plot_procrustes_brain(subjid,mat,stims)
% USAGE: This is a function to plot procrustes weights on the brain
% Input includes a subject ID, and a matrix of values to plot on the brain
% subjid = subject ID
% mat = data matrix to average across
% stims = stimulatio channels 


sid = subjid;
load(fullfile(getSubjDir(subjid),'trodes.mat'))

ave_mat = mean(mat,2);

w = nan(size(Grid, 1), 1);
for i = 1:64
    if (i~=stims)
        w(i) = ave_mat(i);
    end
end
%     w = zeros(size(Grid,1),1);
%     w(stims) = -1;
%     w(beta) = 1;


clims = [0 max(w)];

CT = cbrewer('seq','YlOrRd',11);

figure
PlotDotsDirect(subjid, Grid, w, determineHemisphereOfCoverage(subjid), clims, 20, CT, 1:size(Grid, 1), true);

% very often, after plotting the brain and dots, I add a colorbar for
% reference as to what the dot colors mean
% needs to be the same as what was used in the function call above
colormap(CT);
h = colorbar;
ylabel(h,'Procrustes Distance metric')
title({'Procrustes Distance Metric'})
set(gca,'fontsize', 14)

% plot stim channels 
w = nan(size(Grid, 1), 1);
w(stims) = 1;
color  = [0 1 0];
clims = [0 max(w)];
plot3(Grid(stims,1),Grid(stims,2),Grid(stims,3),'o','MarkerFaceColor',color,'MarkerSize',20,'MarkerEdgeColor','k')


end

