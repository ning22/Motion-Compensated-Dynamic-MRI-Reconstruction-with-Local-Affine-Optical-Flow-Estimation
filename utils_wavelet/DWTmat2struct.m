function w = DWTmat2struct(W11, J)

% w{j}{d} 
%  d = 1,2,3 (orientations, lo_hi (vertical), hi_lo (horizontal),
%  hi_hi (diagonal) respectively in the order of column_row (first filter on columns second on rows)) 
%
% 	      LL | LH
%  Wmn =  --   --
%         HL | HH
%

ROW = size(W11,1);
COL = size(W11,2);

w{J+1} = W11(1:ROW/2^J,1:COL/2^J);

for j = 1:J
    w{j}{1} = W11(1:ROW/2^j,COL/2^j+1:COL/2^(j-1));
    w{j}{2} = W11(ROW/2^j+1:ROW/2^(j-1),1:COL/2^j);
    w{j}{3} = W11(ROW/2^j+1:ROW/2^(j-1),COL/2^j+1:COL/2^(j-1));
end