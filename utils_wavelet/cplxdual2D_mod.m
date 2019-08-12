function w = cplxdual2D_mod(x, J, Faf, af)

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
%       hi_hi (diagonal) respectively in the order of column_row (first filter on columns second on rows))
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

for j = 1:J+1
    if j == J+1
        [wp_r wm_r] = pm(w{j}{1}{1},w{j}{2}{2});
        [wp_i wm_i] = pm(w{j}{2}{1},w{j}{1}{2});
        w{j}{1}{1} = wm_r;
        w{j}{1}{2} = wp_r;
        w{j}{2}{1} = -wp_i;
        w{j}{2}{2} = -wm_i;
    else
        for m = 1:3
            %         [w{j}{1}{1}{m} w{j}{2}{2}{m}] = pm(w{j}{1}{1}{m},w{j}{2}{2}{m});
            %         [w{j}{1}{2}{m} w{j}{2}{1}{m}] = pm(w{j}{1}{2}{m},w{j}{2}{1}{m});
            
            % I think it should have been like this to ensure that we use
            % upper half of the Fourier plane (+ve vertical frequency).
            % The way low pass filter (phi_h + j phi_g) in dualfilt1 gives negative
            % frequency band. The following changes use (phi_h - j phi_g)
            % instead.
            switch m
                case 1
                    [wp_r wm_r] = pm(w{j}{1}{1}{m},w{j}{2}{2}{m});
                    [wp_i wm_i] = pm(w{j}{2}{1}{m},w{j}{1}{2}{m});
                    
                    w{j}{1}{1}{m} = wp_r;
                    w{j}{1}{2}{m} = wm_r;
                    w{j}{2}{1}{m} = -wm_i;
                    w{j}{2}{2}{m} = -wp_i;
                case 2
                    [wp_r wm_r] = pm(w{j}{1}{1}{m},w{j}{2}{2}{m});
                    [wp_i wm_i] = pm(w{j}{2}{1}{m},w{j}{1}{2}{m});
                    
                    w{j}{1}{1}{m} = wp_r;
                    w{j}{1}{2}{m} = wm_r;
                    w{j}{2}{1}{m} = wm_i;
                    w{j}{2}{2}{m} = wp_i;
                case 3
                    [wp_r wm_r] = pm(w{j}{1}{1}{m},w{j}{2}{2}{m});
                    [wp_i wm_i] = pm(w{j}{2}{1}{m},w{j}{1}{2}{m});
                    
                    w{j}{1}{1}{m} = wm_r;
                    w{j}{1}{2}{m} = wp_r;
                    w{j}{2}{1}{m} = wp_i;
                    w{j}{2}{2}{m} = wm_i;
                otherwise
                    disp('na re na');
            end
        end
    end
end