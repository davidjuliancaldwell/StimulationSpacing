%% make the wavelet matrix for SVD 

function [ncPost_m,t_post,fw] = make_waveletMatrix(signal,t,pre_end,pre_begin,post_begin,post_end,fs_data,filter_it,plotIt)

if (~exist('plotIt','var'))
    plotIt = false;
end

fw = [1:3:200]; %frequency bins 

numFreqs = length(fw);

t_post = t(t<post_end & t>post_begin);
t_pre = t(t<pre_end & t>pre_begin);

chans = 1:size(signal,2);

ncPost_m = zeros(numFreqs,length(t_post),length(chans),size(signal,3)); % initialize matrix 

for i = chans
    
    sig = squeeze(signal(:,i,:));
    
    for j = 1:size(sig,2)
        
        if strcmp(filter_it,'y')
            sig_pre = notch(sig((t<pre_end & t>pre_begin),j),[60 120 180 240],fs_data);
            sig_postL = notch(sig((t>post_begin & t<post_end),j),[60 120 180 240],fs_data);
        else
            
            sig_pre = sig((t<pre_end & t>pre_begin),j);
            sig_postL = (sig((t>post_begin & t<post_end),j));
            
        end
        
        if plotIt == true
            figure
            [f_pre,P1_pre] = spectralAnalysis(fs_data,t_pre,sig_pre);
            [f_postL,P1_postL] = spectralAnalysis(fs_data,t_post,sig_postL);
            
            legend({'pre','high'})
        end
        % do some time frequency analysis
        
        if plotIt == true
            figure
        end
        
        [t_post,fw,ncPost] = timeFrequencyAnalWavelet(sig_pre,sig_postL,t_pre,t_post,fs_data,plotIt);
        ncPost_m(:,:,i,j) = ncPost;
        
    end
end

end
