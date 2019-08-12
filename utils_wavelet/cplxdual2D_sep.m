function w = cplxdual2D_sep(x, J, Faf, af)

% Dual-Tree Complex 2D Discrete Wavelet Transform
%
% USAGE:
%   w = cplxdual2D(x, J, Faf, af)
% INPUT:
%   x - 2-D array
%   J - number of stages
%   Faf{i}: first stage filters for tree i
%   af{i}:  filters for remaining stages on tree i
% OUTPUT:
%   w{j}{i}{d1}{d2} - wavelet coefficients
%       j = 1..J (scale)
%       i = 1 (real part); i = 2 (imag part)
%       d1 = 1,2; d2 = 1,2,3 (orientations)
%   w{J+1}{m}{n} - lowpass coefficients
%       d1 = 1,2; d2 = 1,2
%
%%%%%%
% My changes, my notations:
% w{j}{i}{d1}{d2} - wavelet coefficients
%       j = 1..J (scale)
%       i = 1 (real part); i = 2 (imag part)
%       d1 = 1,2; (angle) 1 - positive, 2 - negative.
%       d2 = 1,2,3 (orientations, lo_hi (vertical), hi_lo (horizontal),
%       hi_hi (diagonal) respectively in the order of
%       column_row (first filter on columns (in vertical direction) second
%       on rows (in horizontal direction))
%%%%%%
%
% EXAMPLE:
%   x = rand(256);
%   J = 5;
%   [Faf, Fsf] = FSfarras;
%   [af, sf] = dualfilt1;
%   w = cplxdual2D(x, J, Faf, af);
%   y = icplxdual2D(w, J, Fsf, sf);
%   err = x - y;
%   max(max(abs(err)))
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

% normalization
x = x/2;

for m = 1:2
    for n = 1:2
        [lo w{1}{m}{n}] = afb2D(x, Faf{m}, Faf{n});
        for j = 2:J
            [lo w{j}{m}{n}] = afb2D(lo, af{m}, af{n});
        end
        w{J+1}{m}{n} = lo;
    end
end
%
% {m}{n} in the above set of wavelets correspond to {col}{row}
% where 1 - regular wavelet, 2 - Hilbert transformed wavelet.
%

% Rest of it in seperate2cw.m