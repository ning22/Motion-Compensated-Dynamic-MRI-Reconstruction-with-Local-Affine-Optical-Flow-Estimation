% ICWT2.m
% Fast Inverse 2D wavelet transform using Complex Wavelets
%
% function x = iCWT2_separate(Fsep, Fsf, sf, J, symF, symh, symg);
%
% Usage: x = iCWT2_separate(Fsep, Fsf, sf, J, symF, symh, symg);
%
% Fsep = [Whh Wgg Whg Wgh]; 
% x - nxn input image (n must be a power of 2)
% Wxx - wavelet coeffs. from Dual tree CWT2

% Fsf{i}: First stage synthesis filter coefficients for tree i
% Fsf{i}(:,j) : j = 1: low pass filter
%               j = 2: high pass filter
%
% sf : synthesis filter bank coeffs. for later stages
%
% J - number of levels in the filter bank.
%     Right now, it must be chose so that n*2^(-J) >= length(h0)
% symF, symh, symg - How to extend the input for First stage, Tree h and
% Tree g
%     0: Periodic Extension (circular convolution)
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

function x = iCWT2_separate(Fsep, Fsf, sf, J, symF, symh, symg);

ROW = size(Fsep,1);
COL = size(Fsep,2)/4;

Whh = Fsep(:,1:COL);
Wgh = Fsep(:,COL+1:2*COL);
Whg = Fsep(:,2*COL+1:3*COL);
Wgg = Fsep(:,3*COL+1:end);

if J>1
    % Filter banks for later stages
    h0t = sf{1}(:,1);
    h1t = sf{1}(:,2);
    g0t = sf{2}(:,1);
    g1t = sf{2}(:,2);

    % Later stages outputs
    iFhh = ifwt2_CWT(Whh(1:ROW/2, 1:COL/2), h0t, h1t, h0t, h1t, J-1, symh, symh, 0,0);
    iFgg = ifwt2_CWT(Wgg(1:ROW/2, 1:COL/2), g0t, g1t, g0t, g1t, J-1, symg, symg, 0,0);
    iFhg = ifwt2_CWT(Whg(1:ROW/2, 1:COL/2), h0t, h1t, g0t, g1t, J-1, symh, symg, 0,0);
    iFgh = ifwt2_CWT(Wgh(1:ROW/2, 1:COL/2), g0t, g1t, h0t, h1t, J-1, symg, symh, 0,0);

%     iFhh = ifwt2_CWT(Whh(1:ROW/2, 1:COL/2), h0t, h1t, h0t, h1t, J-1, symh, symh);
%     iFgg = ifwt2_CWT(Wgg(1:ROW/2, 1:COL/2), g0t, g1t, g0t, g1t, J-1, symg, symg);
%     iFhg = ifwt2_CWT(Whg(1:ROW/2, 1:COL/2), h0t, h1t, g0t, g1t, J-1, symh, symg);
%     iFgh = ifwt2_CWT(Wgh(1:ROW/2, 1:COL/2), g0t, g1t, h0t, h1t, J-1, symg, symh);
else
    iFhh = Whh(1:ROW/2, 1:COL/2);
    iFgg = Wgg(1:ROW/2, 1:COL/2);
    iFhg = Whg(1:ROW/2, 1:COL/2);
    iFgh = Wgh(1:ROW/2, 1:COL/2);
end

% First stage filter banks
Fh0t = Fsf{1}(:,1);
Fh1t = Fsf{1}(:,2);
Fg0t = Fsf{2}(:,1);
Fg1t = Fsf{2}(:,2);

Thh = Whh;
Thh(1:ROW/2, 1:COL/2) = iFhh;

Tgg = Wgg;
Tgg(1:ROW/2, 1:COL/2) = iFgg;

Thg = Whg;
Thg(1:ROW/2, 1:COL/2) = iFhg;

Tgh = Wgh;
Tgh(1:ROW/2, 1:COL/2) = iFgh;

% First stage filter bank outputs
% since our first stage filter in tree g is shifted from center. 
% we have shift of -1 so there is no need to use delay (nu+1)
% But for the inverse we need Es(2,1) extension on scaling and Es(1,2)
% extenstion on wavelet coeffs. (which is opposite to what we do in simple
% zero centered Daub79 filters).
% So I use rosym_flip and cosym_flip flags to tell if we need to flip the
% odd symmetry for row or column. The last two parameters are rosym_flip
% and cosym_flip

xhh = ifwt2_CWT(Thh, Fh0t, Fh1t, Fh0t, Fh1t, 1, symF, symF, 0,0);
xgg = ifwt2_CWT(Tgg, Fg0t, Fg1t, Fg0t, Fg1t, 1, symF, symF, 1,1);
xhg = ifwt2_CWT(Thg, Fh0t, Fh1t, Fg0t, Fg1t, 1, symF, symF, 0,1);
xgh = ifwt2_CWT(Tgh, Fg0t, Fg1t, Fh0t, Fh1t, 1, symF, symF, 1,0);

% xhh = ifwt2_CWT(Thh, Fh0t, Fh1t, Fh0t, Fh1t, 1, symF, symF);
% xgg = ifwt2_CWT(Tgg, Fg0t, Fg1t, Fg0t, Fg1t, 1, symF, symF);
% xhg = ifwt2_CWT(Thg, Fh0t, Fh1t, Fg0t, Fg1t, 1, symF, symF);
% xgh = ifwt2_CWT(Tgh, Fg0t, Fg1t, Fh0t, Fh1t, 1, symF, symF);


x = (xhh+xgg+xhg+xgh)/2;