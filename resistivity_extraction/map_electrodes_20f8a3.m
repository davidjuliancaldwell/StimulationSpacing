%% DJC 3-24-2017 - script to map the clinical electrodes to the ones we recorded for subject 20f8a3

sid = '20f8a3';
SUB_DIR = fullfile(myGetenv('subject_dir'));
montageFilepath = strcat(SUB_DIR,'\',sid,'\',sid,'_Montage.mat');
load(montageFilepath);

%%
%interestNames = {'Grid1','Grid2','LTP','LAT','LMT','LPT','LOF','LAHD','LMHD','LAID','LMID','LPID'};
%interestElecs = {[1:32],[1:32],[1,3:6],[1:6],[1:6],[1,3,5],[1,3,5,7],[1,3,5],[1,3,5],[1:8,10,12,14],[2,4:13],[1:12]};

%interestNames = {'Grid2','LTP','LMT','LPT','LOF','LAHD','LMHD','LAID','LMID','LPID'};
%interestElecs = {[1:32],[1,3:6],[1:6],[1,3,5],[1,3,5,7],[1,3,5],[1,3,5],[1:8,10,12,14],[2,4:13],[1:12]};

interestNames = {'Grid1','Grid2','Depth4','Depth1','Depth3','OL','MTL','Depth2'};
interestElecs = {[1:32],[1:16],[1:8],[1:8],[1:8],[1:6],[1:6],[1:8]};

convertBis = 0;

if convertBis
    montageBISFilepath = strcat(SUB_DIR,'\',sid,'\','bis_trodes.mat');
    load(montageBISFilepath);
    Montage.MontageTrodes = AllTrodes; % convert to BIS trodes
end

[locs,identifier] = splitMontage(Montage,sid,interestNames,interestElecs);

% ok - now we have to remap it 

locs = [locs(1:48,:);zeros(16,3);locs(49:78,:);zeros(2,3);locs(79:end,:)];
identifier = [identifier(1:48,:);zeros(16,3);identifier(49:78,:);zeros(2,3);identifier(79:end,:)];

saveIt = 0;
plotIt = 1;
%%
if saveIt
    save('20f8a3_mappedElectrodes_NEW_BIS.mat','locs','identifier','interestNames','interestElecs')
end

%%
if plotIt
    
    figure
    n_elems = max(max(identifier));
    j = 0;
    colors = distinguishable_colors(n_elems+2);
    figure
    [h1,h2] = PlotCortex(sid,'b',[],0.5);
    hold on
    
    for i = 1:n_elems
        tempInd = identifier==i;
        tempInd = tempInd(:,1);
        scatter3(locs(tempInd,1),locs(tempInd,2),locs(tempInd,3),[100],colors(i+1,:),'filled')
        
        
    end
    
    total = 48;
    
    trodeLabels = [1:total];
    for chan = 1:total
        txt = num2str(trodeLabels(chan));
        t = text(locs(chan,1),locs(chan,2),locs(chan,3),txt,'FontSize',10,'HorizontalAlignment','center','VerticalAlignment','middle');
        set(t,'clipping','on');
    end
    
    set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    legend(interestNames)
    %legend({'','',interestNames{:}})
    
    
    %% compare to all electrodes
    
    n_elems = length(Montage.Montage);
    j = 0;
    colors = distinguishable_colors(n_elems);
    figure
    hold on
    
    [h1,h2] = PlotCortex(sid,'r',[],0.5);
    
    
    for i = 1:n_elems
        
        combined_info = split(Montage.MontageTokenized{i},["(",")"]);
        name = combined_info{1};
        elecs = str2num(combined_info{2});
        total = length(elecs);
        
        sub_sel = Montage.MontageTrodes((elecs+j),:);
        locsS = sub_sel;
        j = j + total;
        
        h = scatter3(locsS(:,1),locsS(:,2),locsS(:,3),[100],colors(i,:),'filled');
        
        trodeLabels = [1:total];
        for chan = 1:total
            txt = num2str(trodeLabels(chan));
            t = text(locsS(chan,1),locsS(chan,2),locsS(chan,3),txt,'FontSize',10,'HorizontalAlignment','center','VerticalAlignment','middle');
            set(t,'clipping','on');
        end
    end
    set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
    legend(Montage.MontageTokenized)
    %legend({'','',Montage.MontageTokenized{:}})
end

