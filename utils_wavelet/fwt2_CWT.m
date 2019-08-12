% fwt2_CWT.m
%
% Fast 2D wavelet transform using Complex Wavelets for one pair of tree out
% of four. 
% In CWT2.m we can calculate Whh, Wgg, Whg and Wgh seperately and
% then add/subtract accordingly to get real and imaginary wavelet coeffs. 
%
% Usage: w = fwt2_CWT(x, h0_r, h1_r, h0_c, h1_c, J, symr, symc);
% x - nxn input image (n must be a power of 2)
% h0_r - lowpass decomposition filter on rows
% h1_r - highpass decomposition filter on rows
%      h0 and h1 should be zero padded appropriately so that they have the
%      same length, and this length is even.
% h0_c - lowpass decomposition filter on columns
% h1_c - highpass decomposition filter on columns

% J - number of levels in the filter bank.
%     Right now, it must be chose so that n*2^(-J) >= length(h0)
% symr, symc - How to extend the input for row or column filters
%     1: type-I symmetric extension ( [... x(2) x(1) x(2) x(3) ...])
%            The wavelet filters must have type-I even symmetry 
%            (e.g. daub79)
%     2: type-II symmetric extension ( [... x(2) x(1) x(1) x(2) x(3) ...])
%            The lowpass filter must have type-II even symmetry, 
%            The highpass filter must have type-II even symmetry.
%            (e.g. daub1018)
%
% Written by: Salman Asif
% Created: February 2008