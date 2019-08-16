% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ningning Zhao @ UCLA, department of radiation oncology
% Email: nzhaonzhao@g.ucla.edu
%        buaazhaonn@gmail.com
% some parts are based on code from Salman Asif and Javier Royuela del Val
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mname = mfilename;
mpath = mfilename('fullpath');

addpath(genpath('utils_wavelet'));
addpath(genpath('data'));
addpath('operators');
clear;

%  Parse parameters and options
opts = [];

%%% Reconstruction options
opts.out_maxIter = 3;
opts.maxIter = 600;   
opts.stepsize = 10;
opts.tol = 1e-3;
opts.showTrigger = 20;
opts.stopCriterion = 3;
opts.flag_stopCriterion = 2;

%  regularization parameters
opts.priorType = 'l1+tv3d';
opts.rpSP = 1e-3;  % sparse regularization parameter
opts.rpLR = 1e-3;  % low rank regularization parameter
opts.rpOF = 1e-3;  % optical flow regularization parameter
opts.rpTVx = 1e-5;
opts.rpTV = 1e-5; % TV regularization of the motion vectors
opts.rpVB = 2;
opts.ma = 1e-5;

opts.T=TempFFT(3); % @The CLASS function must be called from a class constructor
opts.scaleSet = [5,4,3]; % multi-scale resulution for motion estimation
opts.periodic = true;
opts.MEmethod = 'Liu';
opts.reduction_factor = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data read    
% readData;
load('cardiac.mat');
I_cube = seq/(max(seq(:)));
I_cube = abs(I_cube);
opts.roi_ver = [81:168];
opts.roi_hor = [56:155];

global ROW COL T_frames
[ROW,COL,T_frames] = size(I_cube);
DS_Type = 'GoldenAngle';
switch DS_Type
    case 'GoldenAngle'
        lines = 30; % Specify the no. of radial lines for under-sampling the k space;
        [mask] = double(goldenratio_samp(ROW,COL,T_frames,lines,0));
        mask = fftshift(fftshift(mask,1),2);
    case 'RadialMask' % too slow when numFras is big
        mask = createRadialMask(ROW,COL,T_frames,100);
        mask = ifftshift(mask);        
    case 'CartesianMask'
        factor = 2;
        [mask,Tsize] = createCartesianMask(ROW,COL,T_frames,factor); 
%         mask = ifftshift(mask);
end
idxList=find(mask~=0);
percent = sum(mask(:))/numel(mask(:));
disp(['Keeping ',num2str(100*percent),' percent of Fourier coefficients.'])

A = @(x)A_fhp3D(x,idxList,ROW,COL,T_frames);
At = @(x)At_fhp3Dc(x,idxList,ROW,COL,T_frames);
b = A(I_cube);

% opts.reduction_factor = 10;
opts.observation = b;
opts.xtrue = I_cube;
opts.xinit = At(b);
opts.A = A;
opts.At = At;
opts.mname = 'Cardiac';
opts.errFcn = @(x) norm( abs(x(:)) - abs(I_cube(:))) / norm(abs(I_cube(:)));
opts.outFcn = @(x) [norm( abs(x(:)) - abs(I_cube(:)), 'inf' ), norm( abs(x(:)) - abs(I_cube(:))) / norm(I_cube(:))];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Output options
opts.save_results = false;
opts.mc_option = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SSIM options
opts.ssim.K       = [0.01 0.03];
opts.ssim.window  = fspecial('gaussian', 11, 1.5);
opts.ssim.maps    = false;
%%% RMSE options

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Wavelets parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
J = floor(log2(min(ROW,COL)))-4; % Number of scales in wavelet transform

% Complex redundant wavelets % (Fhh +/- Fgg, Fhg +/- Fgh)
red = 4;
SYM = 3;
[Faf, Fsf, af, sf] = BiOrthDualFilt_mod;
C2D = 0;    % 1 -- complex coefficients % 0 -- real and imaginary separate
psiT = @(z) CWT2D_op(SD2seq(z), Faf, af, Fsf, sf, J, SYM, C2D, ROW, COL);
psi = @(z) seq2SD(adj_CWT2D_op(z, Faf, af, Fsf, sf, J, SYM, C2D, ROW, COL));

U_sp = @(z) psiT(z);
Ut_sp = @(z) psi(z);
len_U = red*numel(I_cube);
L1reg = 1;
U_kt = @(z) U_sp(z);
Ut_kt = @(z) Ut_sp(z);  
opts.U_h = U_kt;
opts.Ut_h = Ut_kt;

% frame = 1;
% subbox = [81,41,128,128];
% figure(1), set(gca,'FontSize',15,'FontWeight','bold')
% imshow(abs(I_cube(:,:,frame))); 
% hold on, rectangle('Position',subbox, 'LineWidth',2,'LineStyle','--','EdgeColor','b');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Launch reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
method ='PDAL';
opts.rpSP = 5e-4;  % sparse regularization parameter
opts.rpLR = 5e-4;  % low rank regularization parameter
opts.rpOF = 5e-4;  % optical flow regularization parameter
opts.rpTV = 1e-6; 
opts.tol = 5e-4;
opts.tol_mc = 1e-4;
opts.ka = 0.1;
opts.ma = 3e-4;
opts.priorMotion = 'l2';

opts.rpOF = 1e-5;
opts.tol = 1e-5;
opts.ma = 5e-5;
[Ir_mc, VOF, costs] = JPDAL_function_motionPrior(opts);

   

