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

function [Faf, Fsf, af, sf] = BiOrthDualFilt();

Faf{1} = [
                  0                  0
  -0.08838834764832  -0.01122679215254
   0.08838834764832   0.01122679215254
   0.69587998903400   0.08838834764832
   0.69587998903400   0.08838834764832
   0.08838834764832  -0.69587998903400
  -0.08838834764832   0.69587998903400
   0.01122679215254  -0.08838834764832
   0.01122679215254  -0.08838834764832
                  0                  0
 ];
Fsf{1} = Faf{1}(end:-1:1, :);

Faf{2} = [
   0.01122679215254                  0
   0.01122679215254                  0
  -0.08838834764832  -0.08838834764832
   0.08838834764832  -0.08838834764832
   0.69587998903400   0.69587998903400
   0.69587998903400  -0.69587998903400
   0.08838834764832   0.08838834764832
  -0.08838834764832   0.08838834764832
                  0   0.01122679215254
                  0  -0.01122679215254
];
Fsf{2} = Faf{2}(end:-1:1, :);


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
g0 = sqrt(2)*[0 g0n fliplr(g0n) 0]; % low pass filter on analysis side
g0t = sqrt(2)*[0 g0tn fliplr(g0tn) 0 ]; % low pass filter on synthesis side
g1 = -(-1).^(0:11).*fliplr(g0t); % high pass filter on analysis side
g1t = -(-1).^(1:12).*fliplr(g0); % high pass filter on synthesis side

% % According to the paper (this gives 0.5 delay w.r.t 'h' but we can swap 
% 'g' and 'h'. Use 'g' as the filter for real part so that we can use 
% symmetric extension in both branches. 

% The filter is described in the form of g(-n) in the paper. 
% So to use g(-n) on analysis side we flip it. 
% For the first stage when we use delayed filter flipped around 
% zero it means filter is actually advanced by 1). 
% So output of g will be imaginary.
g0_p = sqrt(2)*[0 0 g0n fliplr(g0n)]; % low pass filter on analysis side
g0t_p = sqrt(2)*[g0tn fliplr(g0tn) 0 0]; % low pass filter on synthesis side
g1_p = (-1).^(1:12).*(g0t_p); % high pass filter on analysis side
g1t_p = (-1).^(0:11).*(g0_p); % high pass filter on synthesis side


% Analysis and synthesis filters for later stages
af{1} = [h0' h1'];
af{2} = [g0' g1'];
sf{1} = [h0t' h1t'];
sf{2} = [g0t' g1t'];

% af{1} = [h0' h1'];
% af{2} = [g0_p' g1_p'];
% sf{1} = [h0t' h1t'];
% sf{2} = [g0t_p' g1t_p'];

% First stage filters as Biorthogonal odd 9/7 Daubechies filter
% First stage with filter delayed by 1
h0_p = [0  0 0 fliplr(b) b(2:5)];
h1_p = [0 -fliplr(c) -c(2:4) 0 0 0 0];
h0t_p = ((-1).^(1:12)).*h1_p;
h1t_p = ((-1).^(0:11)).*h0_p;

% % First stage with filter shifted back
% We have to use this one because we have to flip the filters after shift,
% so it means a delay of 1 in h(n) is a shift of -1 in h(-n);
h0_p = [0  fliplr(b) b(2:5) 0 0 ];
h1_p = [0 0 0 -fliplr(c) -c(2:4) 0 0];
h0t_p = ((-1).^(1:12)).*h1_p;
h1t_p = ((-1).^(0:11)).*h0_p;
% % 
Faf{1} = [h0' h1'];
Faf{2} = [h0_p' h1_p'];
Fsf{1} = [h0t' h1t'];
Fsf{2} = [h0t_p' h1t_p'];

% % Filter coeffs. from Selensick Codes
% % % 
% af{1} = [
%    0.03516384000000                  0
%                   0                  0
%   -0.08832942000000  -0.11430184000000
%    0.23389032000000                  0
%    0.76027237000000   0.58751830000000
%    0.58751830000000  -0.76027237000000
%                   0   0.23389032000000
%   -0.11430184000000   0.08832942000000
%                   0                  0
%                   0  -0.03516384000000
%  ];
%  
% af{2} = [
%                   0  -0.03516384000000
%                   0                  0
%   -0.11430184000000   0.08832942000000
%                   0   0.23389032000000
%    0.58751830000000  -0.76027237000000
%    0.76027237000000   0.58751830000000
%    0.23389032000000                  0
%   -0.08832942000000  -0.11430184000000
%                   0                  0
%    0.03516384000000                  0
% ];
%  
% sf{1} = af{1}(end:-1:1, :);
%  
% sf{2} = af{2}(end:-1:1, :);
