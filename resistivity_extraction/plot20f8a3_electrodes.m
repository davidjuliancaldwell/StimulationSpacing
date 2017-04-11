locs = [Grid1; Grid2];

figure
h = scatter3(locs(:,1),locs(:,2),locs(:,3),[100],'filled');
% PlotCortex('20f8a3','r',[],0.25)
hold on

total = length(locs);
trodeLabels = [1:total];
for chan = 1:total
    txt = num2str((chan));
    t = text(locs(chan,1),locs(chan,2),locs(chan,3),txt,'FontSize',10,'HorizontalAlignment','center','VerticalAlignment','middle');
    set(t,'clipping','on');
end

xlabel('L/R')
ylabel('A/P')
zlabel('S/I')

figure
h = plot3(locs(:,1),locs(:,2),locs(:,3),'o','Color',[0 0 0], 'MarkerSize', 15);
trodeLabels = [1:total];
for chan = 1:total
    txt = num2str(chan);
    t = text(locs(chan,1),locs(chan,2),locs(chan,3),txt,'FontSize',10,'HorizontalAlignment','center','VerticalAlignment','middle');
    set(t,'clipping','on');
end

xlabel('L/R')
ylabel('A/P')
zlabel('S/I')