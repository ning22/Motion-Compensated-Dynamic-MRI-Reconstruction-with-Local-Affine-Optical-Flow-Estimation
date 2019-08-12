function W11 = DWTstruct2mat(w, J)

% w{j}{d}
% Need to apply separate2cw to achieve the final CWT form.
%      
% 	      LL | LH
%  Wmn =  --   --
%         HL | HH
%

COL = size(w{1}{1},2)*2;
ROW = size(w{1}{1},1)*2;
W11 = zeros(ROW,COL);

for j = 1:J+1
    if j==J+1
        W11(1:ROW/2^J,1:COL/2^J) = w{J+1};
    else
        W11(1:ROW/2^j,COL/2^j+1:COL/2^(j-1)) = w{j}{1};
        W11(ROW/2^j+1:ROW/2^(j-1),1:COL/2^j) = w{j}{2};
        W11(ROW/2^j+1:ROW/2^(j-1),COL/2^j+1:COL/2^(j-1)) = w{j}{3};
    end
end
