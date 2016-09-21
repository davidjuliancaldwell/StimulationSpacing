%% RPCA example from Nathan Kutz's book
% For this example, the low-rank matrix L has two modes of significance
% while the sparse matrix S has 60 nonzero entries (and needs 26 of them to
% represent 99% of the original data.
close all, clear all, clc

addpath('C:\Users\jcronin\Code\Matlab\Experiments\stimulation_spacing\inexact_alm_rpca')
addpath('C:\Users\jcronin\Code\Matlab\Experiments\stimulation_spacing\inexact_alm_rpca\PROPACK')


n=200;
x=linspace(-10,10,n);
t=linspace(0,10,30);
[X,T]=meshgrid(x,t);
usol=sech(X).*(1-0.5*cos(2*T))+(sech(X).*tanh(X)).*(1-0.5*sin(2*T));
figure(1)
subplot(2,2,1), waterfall(x,t,abs(usol)); colormap([0 0 0]); title('ideal data')

% randintrlv is used to produce noise spikes (60 spikes in total) in
% certain matrix/pixel locations, thus corrupting the data matrix
sam=60;
Atest2=zeros(length(t),n);
Arand1=rand(length(t),n);
Arand2=rand(length(t),n);
r1 = randintrlv([1:length(t)*n],793);
r1k= r1(1:sam);
for j=1:sam
    Atest2(r1k(j))=-1;
end
Anoise=Atest2.*(Arand1+i*Arand2);
unoise=usol+2*Anoise;
subplot(2,2,2), waterfall(x,t,abs(unoise)); colormap([0 0 0]); title('noisy data')

%% SVD
A1=usol; A2=unoise;
[U1,S1,V1]=svd(A1, 'econ');
[U2,S2,V2]=svd(A2, 'econ');

%% Plot the singular values
figure(2)
subplot(3,1,1)
plot(diag(S1),'ko','Linewidth',[2])
title('ideal data')
subplot(3,1,2)
plot(diag(S2),'ko','Linewidth',[2])
title('noisy data')

%% Plot the temporal portions of the first 4 modes
figure(3)
subplot(3,1,1)
plot(x,V1(:,1:4))
title('ideal data'), xlabel('x'), ylabel('V')
legend('mode 1', 'mode 2', 'mode 3', 'mode 4')
subplot(3,1,2)
plot(x,real(V2(:,1:4)))
title('noisy data'), xlabel('x'), ylabel('V')

%% Apply robust PCA
% the corrupt data matrix is decomposed into real and imaginary parts
% before applying the inexact_alm_rpca algorithm
ur=real(unoise);
ui=imag(unoise);

% The parameter lambda used in the inexact_alm_rpca algorithm can be tuned
% to best separate the sparse from low-rank matrices
lambda=0.2; 
% See below for examples of changing lambda

% low-rank (R1) and sparse (R2) portions
[R1r, R2r] = inexact_alm_rpca(real(ur.'),lambda);
[R1i, R2i] = inexact_alm_rpca(real(ui.'),lambda);

% the real and imaginary portions are put back together here to SVD the
% low-rank portion:
R1 = R1r + 1i*R1i;
R2 = R2r + 1i*R2i;

[U3,S3,V3]=svd(R1.');

%% Add the subplots for the RPCA results
figure(1)
subplot(2,2,3), waterfall(x,t,abs(R1)'), colormap([0 0 0]), title('RPCA: low-rank')
subplot(2,2,4), waterfall(x,t,abs(R2)'), colormap([0 0 0]), title('RPCA: sparse')

figure(2)
subplot(3,1,3)
plot(diag(S3),'ko','Linewidth',[2])
title('RPCA: low-rank')

figure(3)
subplot(3,1,3)
plot(x,V3(:,1:4))
title('RPCA: low-rank'), xlabel('x'), ylabel('V')

%% Low-rank approximations using the SVD of the corrupt data matrix vs. the RPCA
% low-rank approximation of corrupt data matrix
figure(4)
for j=2:5
    ff=U2(:,1:j)*S2(1:j,1:j)*V2(:,1:j).';
    subplot(2,2,j-1), waterfall(x,t,abs(ff)), colormap([0 0 0])
    title(['rank= ', num2str(j)])
    xlabel('x'), ylabel('t'), zlabel('|u|')
end
subplot(2,2,1), title(['low-rank approximation of corrupt data matrix;' ' rank= ' num2str(2)])

% approximation of low-rank matrix from robust PCA
figure(5)
for j=2:5
    ff=U3(:,1:j)*S3(1:j,1:j)*V3(:,1:j).';
    subplot(2,2,j-1), waterfall(x,t,abs(ff)), colormap([0 0 0])
    title(['rank= ', num2str(j)])
    xlabel('x'), ylabel('t'), zlabel('|u|')
end
subplot(2,2,1), title(['low-rank approximation  from robust PCA;' ' rank= ' num2str(2)])

%% What happens when you change lambda
lambda=0.5; 

% low-rank (R1) and sparse (R2) portions
[R1r, R2r] = inexact_alm_rpca(real(ur.'),lambda);
[R1i, R2i] = inexact_alm_rpca(real(ui.'),lambda);

% the real and imaginary portions are put back together here to SVD the
% low-rank portion:
R1 = R1r + 1i*R1i;
R2 = R2r + 1i*R2i;

figure(6)
subplot(2,2,1), waterfall(x,t,abs(R1)'), colormap([0 0 0]), title(['RPCA: low-rank; lambda = ', num2str(lambda)])
subplot(2,2,2), waterfall(x,t,abs(R2)'), colormap([0 0 0]), title(['RPCA: sparse; lambda = ', num2str(lambda)])

lambda=0.8; 

% low-rank (R1) and sparse (R2) portions
[R1r, R2r] = inexact_alm_rpca(real(ur.'),lambda);
[R1i, R2i] = inexact_alm_rpca(real(ui.'),lambda);

% the real and imaginary portions are put back together here to SVD the
% low-rank portion:
R1 = R1r + 1i*R1i;
R2 = R2r + 1i*R2i;

subplot(2,2,3), waterfall(x,t,abs(R1)'), colormap([0 0 0]), title(['RPCA: low-rank; lambda = ', num2str(lambda)])
subplot(2,2,4), waterfall(x,t,abs(R2)'), colormap([0 0 0]), title(['RPCA: sparse; lambda = ', num2str(lambda)])

