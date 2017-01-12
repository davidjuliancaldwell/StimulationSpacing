function [sigCCEPs] =  ica_sigCCEPextract(signal,t,pre_begin,pre_end,post_begin,post_end,fs_data,filter_it,plotIt)

if (~exist('plotIt','var'))
    plotIt = false;
end

prompt = {'whats the zscore threshold? e.g. 15'};
dlg_title = 'zthresh';
num_lines = 1;
defaultans = {'10'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
zThresh = str2num(answer{1});
% get the data

zMat = {};
magMat = {};
latencyMat = {};

for idx = 1:size(signal,2)
    
    sig = mean(signal(:,idx,:),3);
    
    if strcmp(filter_it,'y')
        
        sig_pre = notch(sig(t<pre_end & t>pre_begin),[60 120 180 240],fs_data);
        sig_post = notch(sig(t>post_begin & t<post_end),[60 120 180 240],fs_data);
    else
        sig_pre = sig(t>pre_begin & t<pre_end );
        sig_post = sig(t>post_begin & t<post_end);
    end
    
    %[z_ave,mag_ave,latency_ave,w_ave,p_ave,zI,magI,latencyI,wI,pI] = zscoreStimSpacing(signal,signal,t,pre_begin,pre_end,...
    %post_begin,post_end,plotIt);
    
    [z_ave,mag_ave,latency_ave,w_ave,p_ave] = zscoreStimSpacing(sig_pre,sig_post,t,pre_begin,pre_end,...
        post_begin,post_end,plotIt);
    
    zMat{idx} = z_ave;
    magMat{idx} = mag_ave;
    latencyMat{idx} = latency_ave;
    
end


% plot significant CCEPs
zConverted = cell2mat(zMat);

sigCCEPs = find(zConverted>zThresh);

end