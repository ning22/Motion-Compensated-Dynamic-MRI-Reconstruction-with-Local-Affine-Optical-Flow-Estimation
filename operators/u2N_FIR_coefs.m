%------------------------------------------------------------------------
% Copyright or Â© or Copr. CREATIS laboratory, Lyon, France.
% 
% Contributor: Martino Alessandrini, Post Doctoral Fellow at the 
% Centre de Recherche en Acquisition et Traitement de l'Image pour la Santé
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
% M. Alessandrini, A. Basarab, H. Liebgott and O. Bernard, "Myocardial 
% Motion Estimation from Medical Images Using the Monogenic Signal, 
% accepted for buplication to IEEE Transactions on Image Processing
%------------------------------------------------------------------------

function [u2N, c] = u2N_FIR_coefs(N)
% u2n_FIR_coefs - filter coefficients for B-spline scale relation
%
% [u2N, c] = u2N_FIR_coefs(N) generates the FIR filter u2N which relates two
% scalings of the b-spline basis:
% 
% b_N_2 = c*conv(u2N, b_N)
%
% where b_N and b_N_2 are two discretized versions of the basic B-spline
% of order N, b_N_2 is dilated with a factor of two:
%
%       b_N_2 = beta_N(k/2) for k = -Inf, Inf
%       b_N = beta_N(k) for k = -Inf, Inf
%
% c is an additional scale factor.
u2N = zeros(1,N+2);
c = 2^(-N);

for k = (-(N+1)/2):((N+1)/2);
    u2N(k+(N+1)/2+1) = nchoosek((N+1),k+(N+1)/2);
end
