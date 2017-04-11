%% DJC - 3/27/2017 - Function to split montage in desired parts

function [locsTotal,identifier] = splitMontage(Montage,sid,interestNames,interestElecs,allElecs)
% need montage loaded

if (~exist('allElecs','var'))
    allElecs = 0;
end

if ~allElecs
    locsTotal = [];
    names = {};
    identifier = [];
    
    % keep track of where in the total trodes you are
    for i = 1:length(interestNames)
        
        ind = cellfun(@(x)( ~isempty(x) ), regexp(Montage.MontageTokenized, interestNames{i}));
        
        % keep track of where in the total trodes you are
        if find(ind) == 1
            startInd = 1;
        else
            startInd = sum(Montage.Montage(1:(find(ind)-1)))+1;
        end
        
        locsOfInterest = Montage.MontageTrodes(startInd+interestElecs{i}-1,:);
        locsTotal = [locsTotal; locsOfInterest];
        
        identifierTemp = ones(size(locsOfInterest))*i;
        
        identifier = [identifier; identifierTemp];
    end
    
end



if allElecs
    
    n_elems = length(Montage.Montage);
    j = 0;
    colors = distinguishable_colors(n_elems);
    figure
    %PlotCortex(sid,'b',[],0.5)
    hold on
    
    
    for i = 1:n_elems
        
        combined_info = split(Montage.MontageTokenized{i},["(",")"]);
        name = combined_info{1};
        elecs = str2num(combined_info{2});
        total = length(elecs);
        
        sub_sel = Montage.MontageTrodes((elecs+j),:);
        locsS = sub_sel;
        j = j + total;
        
        scatter3(locsS(:,1),locsS(:,2),locsS(:,3),[100],colors(i,:),'filled');
   
        
        trodeLabels = [1:total];
        for chan = 1:total
            txt = num2str(trodeLabels(chan));
            t = text(locsS(chan,1),locsS(chan,2),locsS(chan,3),txt,'FontSize',10,'HorizontalAlignment','center','VerticalAlignment','middle');
            set(t,'clipping','on');
        end
        
    end
    
end

end