% fit a plane to the electrodes - DJC 4-6-2017

load('20f8a3_temporary_BIS.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
locs = [Grid1; Grid2];

locs(:,1) = locs(:,1) - mean(locs(:,1));
locs(:,2) = locs(:,2) - mean(locs(:,2));
locs(:,3) = locs(:,3) - mean(locs(:,3));

sf = fit([locs(:,1),locs(:,2)],locs(:,3),'cubicinterp');
figure
H = plot(sf);
H.FaceAlpha = 0.4;
H.EdgeColor = [0.25 0.25 0.25];
H.FaceColor = 'none';
hold on
%%%%%%%%%%%%%%%%%%%%%%%
g1_size = length(Grid1);

scatter3(locs(1:g1_size,1),locs(1:g1_size,2),locs(1:g1_size,3),[150],'MarkerFaceColor',[0 .8 .8]);

total = g1_size;
trodeLabels = [1:total];
for chan = 1:total
    txt = num2str((chan));
    t = text(locs(chan,1),locs(chan,2),locs(chan,3),txt,'FontSize',10,'HorizontalAlignment','center','VerticalAlignment','middle');
    set(t,'clipping','on');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g2_size = length(Grid2);

scatter3(locs(g1_size+1:end,1),locs(g1_size+1:end,2),locs(g1_size+1:end,3),[150],'MarkerFaceColor',[1 .75 .75]);

total = g2_size;
trodeLabels = [1:total];
for chan = 1:total
    txt = num2str((chan));
    t = text(locs(chan+g1_size,1),locs(chan+g1_size,2),locs(chan+g1_size,3),txt,'FontSize',10,'HorizontalAlignment','center','VerticalAlignment','middle');
    set(t,'clipping','on');
end


xlabel('L/R')
ylabel('A/P')
zlabel('S/I')

