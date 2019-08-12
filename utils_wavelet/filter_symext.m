function y = filter_symext(x,h,sym, shift)

% x: input signal
% h: filter coeffs.
% sym: type of symmetric extension 1 or 2 for (1,1) or (2,2).
% shift: shift of center from zero. (only working for (1,1) symmetric ext.)
% y: output of a single subband

% This function returns downsampled value after passing the 
% specified symmetric extension through the filter. 
%
% For odd biorthogonal filter sym = 1 (1,1) sym. extension
% For even biorthogonal filter sym = 2 (2,2) sym. extension
% 
% shift takes care of the symmetry when filter is not centered at the zero
%
% On analysis side if scaling filter is delayed by 'p' samples, wavelet 
% filter will be advanced by 'p' samples, and converse will hold on the 
% synthesis side.
%
% References: 
% Christopher Brislawn - Classification of nonexpansive
% symmetric extension transforms for multirate filter banks
% "Preservation of subband symmetry in multirate signal coding" IEEE TSP 95
% Satyabrata Rout - MS Thesis (Virginia Tech) "Orthogonal vs. Biorthogonal
% wavelets for image compression"
%
% Written by Salman Asif - Georgia Tech
% Created: February 2008

l = length(h);
L = l/2;
ln = l-L-1;
n = length(x);
y = zeros(1,n/2);

if sym == 2
    for m = 0:n/2-1
        taps = [];
        indces = [];
        for k = 0:l-1
            % kx = k-L+1+m;
            kx = k+2*m-ln+shift;
            if(kx<0)
                a = mod((-kx-1),2*n);
                if a > n-1
                    kxp = 2*n-a-1;
                else
                    kxp = a;
                end

            else
                if (kx>n-1)
                    a = mod(kx,2*n);
                    if a > n-1
                        kxp = 2*n-a-1;
                    else
                        kxp = a;
                    end

                else
                    kxp = kx;

                end
            end
            taps = [taps l-k-1];
            indces = [indces kxp];
            y(m+1) = y(m+1)+x(kxp+1)*h(l-k);
        end
    end
end


if sym == 1
    for m = 0:n/2-1
        taps = [];
        indces = [];
        for k = 0:l-1
            % kx = k-L+1+m;
            kx = 2*m - ln + k + 2*shift;
            if(kx<0)
                a = mod((-kx),2*(n-1));
                if a > n-1
                    kxp = 2*(n-1)-a;
                else
                    kxp = a;
                end

            else
                if (kx>n-1)
                    a = mod(kx,2*(n-1));
                    if a > n-1
                        kxp = 2*(n-1)-a;
                    else
                        kxp = a;
                    end

                else
                    kxp = kx;

                end
            end
            taps = [taps l-k-1];
            indces = [indces kxp];
            y(m+1) = y(m+1)+x(kxp+1)*h(l-k);
        end
%         m
%         taps
%         indces
    end
end
