% adj_CWT2.m
% Adjoint of CWT2
%
% function u = adj_CWT2([Whh, Wgg, Whg, Wgh], Faf, af, J, symF, symh, symg);
%
% Usage: u = adj_CWT2([Whh, Wgg, Whg, Wgh], Faf, af, J, symF, symh, symg);
% u - nxn image (n must be a power of 2)
% FO2D = [Whh-Wgg; Whh+Wgg; Wgh+Whg; Wgh-Whg]/sqrt(4);

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

function u = adj_CWT2(FO2D, Faf, af, J, symF, symh, symg);

ROW = size(FO2D,1);
COL = size(FO2D,2)/4;

Whh = FO2D(:,1:COL)+FO2D(:,COL+1:2*COL);
Wgg = -FO2D(:,1:COL)+FO2D(:,COL+1:2*COL);
Wgh = FO2D(:,2*COL+1:3*COL)+FO2D(:,3*COL+1:end);
Whg = FO2D(:,2*COL+1:3*COL)-FO2D(:,3*COL+1:end);

if J>1
    % Filter banks for later stages
    h0 = af{1}(:,1);
    h1 = af{1}(:,2);
    g0 = af{2}(:,1);
    g1 = af{2}(:,2);

    % Later stages outputs
    aFhh = afwt2_CWT(Whh(1:ROW/2,1:COL/2), h0, h1, h0, h1, J-1, symh, symh);
    aFgg = afwt2_CWT(Wgg(1:ROW/2,1:COL/2), g0, g1, g0, g1, J-1, symg, symg);
    aFhg = afwt2_CWT(Whg(1:ROW/2,1:COL/2), h0, h1, g0, g1, J-1, symh, symg);
    aFgh = afwt2_CWT(Wgh(1:ROW/2,1:COL/2), g0, g1, h0, h1, J-1, symg, symh);

else
    aFhh = Whh(1:ROW/2,1:COL/2);
    aFgg = Wgg(1:ROW/2,1:COL/2);
    aFhg = Whg(1:ROW/2,1:COL/2);
    aFgh = Wgh(1:ROW/2,1:COL/2);
end

% First stage filter banks
Fh0 = Faf{1}(:,1);
Fh1 = Faf{1}(:,2);
Fg0 = Faf{2}(:,1);
Fg1 = Faf{2}(:,2);

Thh = Whh;
Thh(1:ROW/2,1:COL/2) = aFhh;

Tgg = Wgg;
Tgg(1:ROW/2,1:COL/2) = aFgg;

Thg = Whg;
Thg(1:ROW/2,1:COL/2) = aFhg;

Tgh = Wgh;
Tgh(1:ROW/2,1:COL/2) = aFgh;

% First stage filter bank outputs

uhh = afwt2_CWT(Thh, Fh0, Fh1, Fh0, Fh1, 1, symF, symF);
ugg = afwt2_CWT(Tgg, Fg0, Fg1, Fg0, Fg1, 1, symF, symF);
uhg = afwt2_CWT(Thg, Fh0, Fh1, Fg0, Fg1, 1, symF, symF);
ugh = afwt2_CWT(Tgh, Fg0, Fg1, Fh0, Fh1, 1, symF, symF);

u = (uhh+ugg+uhg+ugh)/sqrt(8);
