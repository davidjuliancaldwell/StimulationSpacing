%% DJC 3-24-2017 - script to map the clinical electrodes to the ones we recorded for subject 2fd831

sid = '2fd831';
SUB_DIR = fullfile(myGetenv('subject_dir'));
montageFilepath = strcat(SUB_DIR,'\',sid,'\',sid,'_Montage.mat');
load(montageFilepath);

interestNames = {'Grid1','Grid2','LTP','LAT','LMT','LPT','LOF','LAHD','LMHD','LAID','LMID','LPID'};
interestElecs = {[1:32],[1:32],[1,3:6],[1:6],[1:6],[1,3,5],[1,3,5,7],[1,3,5],[1,3,5],[1:8,10,12,14],[2,4:13],[1:12]};

%interestNames = {'Grid2','LTP','LMT','LPT','LOF','LAHD','LMHD','LAID','LMID','LPID'};
%interestElecs = {[1:32],[1,3:6],[1:6],[1,3,5],[1,3,5,7],[1,3,5],[1,3,5],[1:8,10,12,14],[2,4:13],[1:12]};

%interestNames = {'LPT'};
%interestElecs = {1};

convertBis = 0;

if convertBis
    montageBISFilepath = strcat(SUB_DIR,'\',sid,'\','bis_trodes.mat');
    load(montageBISFilepath);
    Montage.MontageTrodes = AllTrodes; % convert to BIS trodes
end

%%
[locs,identifier] = splitMontage(Montage,sid,interestNames,interestElecs);

saveIt = 0;
plotIt = 1;
%%
if saveIt && convertBis
    save('2fd831_BIS_mappedElectrodes.mat','locs','identifier','interestNames','interestElecs')
elseif saveIt && ~convertBis
    save('2fd831_mappedElectrodes.mat','locs','identifier','interestNames','interestElecs')
end

%%
if plotIt
    
    figure
    n_elems = max(max(identifier));
    j = 0;
    colors = distinguishable_colors(n_elems);
    figure
    [h1,h2] = PlotCortex(sid,'b',[],0.5);
    hold on
    
    for i = 1:n_elems
        tempInd = identifier==i;
        tempInd = tempInd(:,1);
        scatter3(locs(tempInd,1),locs(tempInd,2),locs(tempInd,3),[100],colors(i,:),'filled')
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
    
    [h1,h2] = PlotCortex(sid,'b',[],0.5);
    
    
    for i = 1:n_elems
        
        combined_info = split(Montage.MontageTokenized{i},["(",")"]);
        name = combined_info{1};
        elecs = str2num(combined_info{2});
        total = length(elecs);
        
        sub_sel = Montage.MontageTrodes((elecs+j),:);
        locsS = sub_sel;
        j = j + total;
        
        h = scatter3(locsS(:,1),locsS(:,2),locsS(:,3),[100],colors(i,:),'filled');
        
        
    end
    set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
    legend(Montage.MontageTokenized)
    %legend({'','',Montage.MontageTokenized{:}})
end

