function w = CWTmat2struct(W_tot, J)

% w{j}{m}{n}{d} 
% {m}{n} in the above set of wavelets correspond to {col}{row}
% where 1 - regular wavelet, 2 - Hilbert transformed wavelet.       
%  d = 1,2,3 (orientations, lo_hi (vertical), hi_lo (horizontal),
%  hi_hi (diagonal) respectively in the order of column_row (first filter on columns second on rows))
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

ROW = size(W_tot,1);
COL = size(W_tot,2)/4;

W11 = W_tot(:,1:COL);
W12 = W_tot(:,COL+1:2*COL);
W21 = W_tot(:,2*COL+1:3*COL);
W22 = W_tot(:,3*COL+1:4*COL);

for i = 1:2
    for d1 = 1:2
        switch (i-1)*2+d1
            case 1
                w{J+1}{i}{d1} = W11(1:ROW/2^J,1:COL/2^J);
            case 2
                w{J+1}{i}{d1} = W12(1:ROW/2^J,1:COL/2^J);
            case 3
                w{J+1}{i}{d1} = W21(1:ROW/2^J,1:COL/2^J);
            case 4
                w{J+1}{i}{d1} = W22(1:ROW/2^J,1:COL/2^J);
            otherwise
                disp('na na re na');
        end
    end
end
for i = 1:2
    for d1 = 1:2
        for j = 1:J
            switch (i-1)*2+d1
                
                case 1
                    w{j}{i}{d1}{1} = W11(1:ROW/2^j,COL/2^j+1:COL/2^(j-1));
                    w{j}{i}{d1}{2} = W11(ROW/2^j+1:ROW/2^(j-1),1:COL/2^j);
                    w{j}{i}{d1}{3} = W11(ROW/2^j+1:ROW/2^(j-1),COL/2^j+1:COL/2^(j-1));
                case 2
                    w{j}{i}{d1}{1} = W12(1:ROW/2^j,COL/2^j+1:COL/2^(j-1));
                    w{j}{i}{d1}{2} = W12(ROW/2^j+1:ROW/2^(j-1),1:COL/2^j);
                    w{j}{i}{d1}{3} = W12(ROW/2^j+1:ROW/2^(j-1),COL/2^j+1:COL/2^(j-1));
                case 3
                    w{j}{i}{d1}{1} = W21(1:ROW/2^j,COL/2^j+1:COL/2^(j-1));
                    w{j}{i}{d1}{2} = W21(ROW/2^j+1:ROW/2^(j-1),1:COL/2^j);
                    w{j}{i}{d1}{3} = W21(ROW/2^j+1:ROW/2^(j-1),COL/2^j+1:COL/2^(j-1));
                case 4
                    w{j}{i}{d1}{1} = W22(1:ROW/2^j,COL/2^j+1:COL/2^(j-1));
                    w{j}{i}{d1}{2} = W22(ROW/2^j+1:ROW/2^(j-1),1:COL/2^j);
                    w{j}{i}{d1}{3} = W22(ROW/2^j+1:ROW/2^(j-1),COL/2^j+1:COL/2^(j-1));
                otherwise
                    disp('chal oye chal');
            end
        end
    end
end