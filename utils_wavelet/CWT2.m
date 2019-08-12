% CWT2.m
% Fast 2D wavelet transform using Complex Wavelets
%
% function FO2D = CWT2(x, Faf, af, J, symF, symh, symg);
%
% Usage: FO2D = CWT2(x, Faf, af, J, symF, symh, symg);
% x - nxn input image (n must be a power of 2)
%
% FO2D - 4nxn output containing real and imaginary parts of CWT
%
% Faf{i}: First stage analysis filter coefficients for tree i
% Faf{i}(:,j) : j = 1: low pass filter
%               j = 2: high pass filter
%
% af : analysis filter bank coeffs. for later stages
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

function FO2D = CWT2(x, Faf, af, J, symF, symh, symg);

% if (size(x,1) ~= size(x,2))
%     error('This method Works only with square images');
% end

[ROW COL] = size(x);

% First stage filter banks
Fh0 = Faf{1}(:,1);
Fh1 = Faf{1}(:,2);
Fg0 = Faf{2}(:,1);
Fg1 = Faf{2}(:,2);

% First stage filter bank outputs
Fhh = fwt2_CWT(x, Fh0, Fh1, Fh0, Fh1, 1, symF, symF);
Fgg = fwt2_CWT(x, Fg0, Fg1, Fg0, Fg1, 1, symF, symF);
Fhg = fwt2_CWT(x, Fh0, Fh1, Fg0, Fg1, 1, symF, symF);
Fgh = fwt2_CWT(x, Fg0, Fg1, Fh0, Fh1, 1, symF, symF);

if (J > 1)
    % Filter banks for later stages
    h0 = af{1}(:,1);
    h1 = af{1}(:,2);
    g0 = af{2}(:,1);
    g1 = af{2}(:,2);

    % Later stages outputs
    Thh = fwt2_CWT(Fhh(1:ROW/2,1:COL/2), h0, h1, h0, h1, J-1, symh, symh);
    Tgg = fwt2_CWT(Fgg(1:ROW/2,1:COL/2), g0, g1, g0, g1, J-1, symg, symg);
    Thg = fwt2_CWT(Fhg(1:ROW/2,1:COL/2), h0, h1, g0, g1, J-1, symh, symg);
    Tgh = fwt2_CWT(Fgh(1:ROW/2,1:COL/2), g0, g1, h0, h1, J-1, symg, symh);

    % Put all the coeffs. into output images.
    Whh = Fhh;
    Whh(1:ROW/2,1:COL/2) = Thh;
    
    Wgg = Fgg;
    Wgg(1:ROW/2,1:COL/2) = Tgg;

    Whg = Fhg;
    Whg(1:ROW/2,1:COL/2) = Thg;

    Wgh = Fgh;
    Wgh(1:ROW/2,1:COL/2) = Tgh;
    
else
    Whh = Fhh;
    Wgg = Fgg;
    Whg = Fhg;
    Wgh = Fgh;
end
FO2D = [Whh-Wgg, Whh+Wgg, Wgh+Whg, Wgh-Whg]/sqrt(8);