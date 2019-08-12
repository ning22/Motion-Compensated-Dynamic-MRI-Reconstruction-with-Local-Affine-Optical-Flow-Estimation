%------------------------------------------------------------------------
% Copyright or © or Copr. CREATIS laboratory, Lyon, France.
% 
% Contributor: Martino Alessandrini, Post Doctoral Fellow at the 
% Centre de Recherche en Acquisition et Traitement de l'Image pour la Sant�
% CREATIS (CNRS 5220, INSERM U630, INSA, Claude Bernard Lyon 1 University) 
% in France (Lyon).
% 
% Date of creation: March 26th 2012
% 
% E-mail of the author: martino.alessandrini@creatis.insa-lyon.fr
% 
% This folder provides a MATLAB implementation of an Optical Flow
% estimation algorithm based on the monogenic phase. Given two input images
% the algorithm compute the displacement field between the two by assuming
% the conservation of the monogenic phase. This feature is much less
% sensitive to changes in the illumination conditions as compared to the
% traditional pixel intensity. To reduce dependency on the size of the
% windowing function, the computation is carried out at different scales in
% a coarse-to-fine fashion. The estimation is then refined iteratively in a
% pyramidal scheme.
% 
% The algorithm herein implemented is described in:
% M. Alessandrini, A. Basarab, H. Liebgott and O. Bernard, "Multiscale 
% Optical Flow Computation from the Monogenic Signal", submitted fot
% buplication to IEEE Transactions on Image Processing
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Built filter set for multiresolution moments computation. Each output is
% a cell with a number of entries equal to the number of scales.
% 
% The implementation follows the paper:
% S�hling M, Arigovindan M, Hunziker P, Unser M., "Multiresolution moment 
% filters: theory and applications." EEE Trans Image Process. 2004 Apr;
% 13(4):484-95. 
%------------------------------------------------------------------------

function [Bpyramyd, Wfilt, W2filt] = build_b_pyramid(J,bspline)

% bspline coefficients at scale 0
N = str2double(bspline(2)); % bspline order

if mod(N,2) == 0
    error('use odd order bsplines');
end

x0 = -(N+1)/2:(N+1)/2; 

bNjj = bsplineN(x0,N);
Bpyramyd = cell(J+1,1);
Wfilt = cell(J+1,1); % filters for 0 moments
W2filt = cell(J+1,1); % filters for 0 moments

Bpyramyd{1} = bNjj;
Wfilt{1} = fliplr(x0).*bNjj;
W2filt{1} = x0.^2.*bNjj;

% first scale
[u2N, c] = u2N_FIR_coefs(N);

for jj = 1:J
    upfac = 2^(jj-1);
    ujjN = upsample(u2N,upfac); 
    ujjN = ujjN(1:end-upfac+1); % remove final zeros
    bNjj = c*conv(ujjN,bNjj);
    Bpyramyd{jj+1} = bNjj;
    
    n = numel(bNjj);
    x0 = -(n-1)/2:(n-1)/2;
    
    Wfilt{jj+1} = bNjj.*fliplr(x0);
    W2filt{jj+1} = bNjj.*x0.^2;
end



