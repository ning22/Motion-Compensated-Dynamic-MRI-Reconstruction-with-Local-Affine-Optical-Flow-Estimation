function W_tot = CWTstruct2mat(w, J)

% w{j}{m}{n} 
% {m}{n} in the above set of wavelets correspond to {col}{row}
% where 1 - regular wavelet, 2 - Hilbert transformed wavelet.
%
% Need to apply separate2cw to achieve the final CWT form.
% 
% For example with separate wavelet coefficients 
% W11 W12 W21 W22 - matrix forms of wavelet coefficients denote
% Wmn: {m,n} -- 1 or 2, corresponding to h or g filter on {col,row}
% OR with combined coefficients 
% Wmn: m - 1 or 2, real or imaginary
%      n - 1 or 2, positive or negative orientation
%      
% 	      LL | LH
%  Wmn =  --   --
%         HL | HH
%

COL = size(w{1}{1}{1}{1},2)*2;
ROW = size(w{1}{1}{1}{1},1)*2;
W = zeros(ROW,COL);
W11 = W;
W12 = W;
W21 = W;
W22 = W;

for i = 1:2
    for d1 = 1:2
        for j = 1:J+1
            if j==J+1
                switch (i-1)*2+d1
                    case 1
                        W11(1:ROW/2^J,1:COL/2^J) = w{J+1}{i}{d1};
                    case 2
                        W12(1:ROW/2^J,1:COL/2^J) = w{J+1}{i}{d1};
                    case 3
                        W21(1:ROW/2^J,1:COL/2^J) = w{J+1}{i}{d1};
                    case 4
                        W22(1:ROW/2^J,1:COL/2^J) = w{J+1}{i}{d1};
                    otherwise
                        disp('na na re na');
                end
            else
                switch (i-1)*2+d1
                    
                    case 1
                        W11(1:ROW/2^j,COL/2^j+1:COL/2^(j-1)) = w{j}{i}{d1}{1};
                        W11(ROW/2^j+1:ROW/2^(j-1),1:COL/2^j) = w{j}{i}{d1}{2};
                        W11(ROW/2^j+1:ROW/2^(j-1),COL/2^j+1:COL/2^(j-1)) = w{j}{i}{d1}{3};
                    case 2
                        W12(1:ROW/2^j,COL/2^j+1:COL/2^(j-1)) = w{j}{i}{d1}{1};
                        W12(ROW/2^j+1:ROW/2^(j-1),1:COL/2^j) = w{j}{i}{d1}{2};
                        W12(ROW/2^j+1:ROW/2^(j-1),COL/2^j+1:COL/2^(j-1)) = w{j}{i}{d1}{3};
                    case 3
                        W21(1:ROW/2^j,COL/2^j+1:COL/2^(j-1)) = w{j}{i}{d1}{1};
                        W21(ROW/2^j+1:ROW/2^(j-1),1:COL/2^j) = w{j}{i}{d1}{2};
                        W21(ROW/2^j+1:ROW/2^(j-1),COL/2^j+1:COL/2^(j-1)) = w{j}{i}{d1}{3};
                    case 4
                        W22(1:ROW/2^j,COL/2^j+1:COL/2^(j-1)) = w{j}{i}{d1}{1};
                        W22(ROW/2^j+1:ROW/2^(j-1),1:COL/2^j) = w{j}{i}{d1}{2};
                        W22(ROW/2^j+1:ROW/2^(j-1),COL/2^j+1:COL/2^(j-1)) = w{j}{i}{d1}{3};
                    otherwise
                        disp('chal oye chal');
                end
            end
        end
    end
end
W_tot = [W11, W12, W21, W22];