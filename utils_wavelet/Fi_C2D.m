function W_tot = Fi_C2D(W_c, ROW,COL); 
% W_tot: N x 4*N matrix with 
% W11 W12 W21 W22 - matrix forms of wavelet coefficients.
% Wmn: m - 1 or 2, real or imaginary
%      n - 1 or 2, positive or negative orientation
% W_c: complex for of CWT

W_c = reshape(W_c,ROW,4*COL);

W1 = W_c(:,1:COL);
W2 = W_c(:,COL+1:2*COL);
W3 = W_c(:,2*COL+1:3*COL);
W4 = W_c(:,3*COL+1:4*COL);

W11 = (W1+W3);
W21 = (-1i*W1+1i*W3);
W12 = (W2+W4);
W22 = (-1i*W2+1i*W4);

W_tot = [W11 W12 W21 W22]/sqrt(2);