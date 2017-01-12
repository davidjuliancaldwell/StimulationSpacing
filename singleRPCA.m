function [ R1r, R2r ] = singleRPCA( data, lambda, tol, maxIter, dispTTL)
%This function applies a single robust PCA iteration to real data with a
%given lambda
% The parameter lambda used in the inexact_alm_rpca algorithm can be tuned
% to best separate the sparse from low-rank matrices

%% Apply robust PCA
% the corrupt data matrix is decomposed into real and imaginary parts
% before applying the inexact_alm_rpca algorithm
ur=real(data);

% low-rank (R1) and sparse (R2) portions
[R1r, R2r] = inexact_alm_rpca(real(ur.'),lambda, tol, maxIter, dispTTL);
end

