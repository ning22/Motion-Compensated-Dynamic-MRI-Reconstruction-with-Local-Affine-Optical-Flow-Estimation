% BiOrthDualFilt.m 

% First stage filters are Farras filters organized for the dual-tree
% complex DWT.
%
% Later stages have two seperate filter types as explained in paper by Yu
% and Ozkaramanli ("Hilbert Transform Pairs of Biorthogonal Wavelet Bases",
% IEEE Trans. on SP, vol. 54, No. 6, June 2006)
%
% Tree h has 9/7 Biothogonal filters
% Tree g has 10/10 Biorthogonal filters designed by Yu and Ozkaramanli
%
% First stage Farras filter taken from WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/
%
% Usage: [Faf, Fsf, af, sf] = BiOrthDualFilt();
% Outputs: Faf, Fsf, af, sf
%
% Faf{i}: First stage analysis filter coefficients for tree i
% Faf{i}(:,j) : j = 1: low pass filter
%               j = 2: high pass filter
%
% Fsf : First stage synthesis filter bank coefficients,
% af, sf : analysis and synthesis filter bank coeffs. for later stages
% 
% Filter coefficients for the first stage in Complex wavelet Filter bank
%
% Faf : First stage analysis filter for tree 1 (Faf{1}) & 2 (Faf{2})
% Fsf : First stage synthesis filter for tree 1 & 2 
%
% First column is low-pass / scaling filter
%

function [Faf, Fsf, af, sf] = BiOrthDualFilt_mod()

% Filter bank for Tree h
% daub79.m
%
% Returns the filter coefficients (for lowpass and highpass) for the
% Daubechies 7,9 biorthogonal wavelet set (9/7 Biorthogonal also known as
% Antonini filters)
% h0 - lowpass analysis
% h1 - highpass analysis
% ht0 - lowpass synthesis
% ht1 - highpass synthesis
%

b = sqrt(2)*[0.6029490182363579 0.2668641184428723 -0.07822326652898785 ...
  -0.01686411844287495 0.02674875741080976];
c = sqrt(2)/2*[1.115087052456994 -0.5912717631142470 -0.05754352622849957 ...
      0.09127176311424948]; 

h0 = [0 fliplr(b) b(2:5)];
h1 = [0 -fliplr(c) -c(2:4) 0 0];
h0t = ((-1).^(1:10)).*h1;
h1t = ((-1).^(0:9)).*h0;

% Filter bank for Tree g
%
% g0 - lowpass analysis
% g1 - highpass analysis
% gt0 - lowpass synthesis
% gt1 - highpass synthesis
% 
% Frist half of coefficients for decomposition scaling filter
g0n = [.00509532607545 .01801021269583 -.06233231520576 .03683791155347 .50238886488101 ];
% Frist half of coefficients for reconstruction scaling filter
g0tn = [.00271860510468 -.00960932734160 -.06326036368884 .08639550120371 .48375558472205];

% Here is the even filter delayed by 0.5 samples from odd Daub79 filter
g0 = sqrt(2)*[g0n fliplr(g0n)]; % low pass filter on analysis side
g0t = sqrt(2)*[g0tn fliplr(g0tn)]; % low pass filter on synthesis side
g1 = -(-1).^(1:10).*fliplr(g0t); % high pass filter on analysis side
g1t = -(-1).^(0:9).*fliplr(g0); % high pass filter on synthesis side

% Analysis and synthesis filters for later stages
af{1} = [h0' h1'];
af{2} = [g0' g1'];
sf{1} = [h0t' h1t'];
sf{2} = [g0t' g1t'];

% % First stage with filter shifted back
% We have to use this one because we have to flip the filters after shift,
% so it means a delay of 1 in h(n) is a shift of -1 in h(-n);

h0 = [0 fliplr(b) b(2:5)];
h1 = [0 -fliplr(c) -c(2:4) 0 0];
h0t = ((-1).^(1:10)).*h1;
h1t = ((-1).^(0:9)).*h0;

h0_p = [fliplr(b) b(2:5) 0 ];
h1_p = [0 0 -fliplr(c) -c(2:4) 0];
h0t_p = ((-1).^(0:9)).*h1_p;
h1t_p = ((-1).^(1:10)).*h0_p;

Faf{1} = [h0' h1'];
Faf{2} = [h0_p' h1_p'];
Fsf{1} = [h0t' h1t'];
Fsf{2} = [h0t_p' h1t_p'];

