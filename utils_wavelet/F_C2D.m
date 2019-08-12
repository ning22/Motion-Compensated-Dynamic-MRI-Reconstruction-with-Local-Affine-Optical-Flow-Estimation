function W_c = F_C2D(W_tot,ROW,COL); 
% W_tot: N x 4*N matrix with 
% W11 W12 W21 W22 - matrix forms of wavelet coefficients.
% Wmn: m - 1 or 2, real or imaginary
%      n - 1 or 2, positive or negative orientation
% W_c: complex for of CWT

W_tot = reshape(W_tot,ROW,4*COL);

W11 = W_tot(:,1:COL);
W12 = W_tot(:,COL+1:2*COL);
W21 = W_tot(:,2*COL+1:3*COL);
W22 = W_tot(:,3*COL+1:4*COL);

W1 = W11+1i*W21;
W2 = W12+1i*W22;
W3 = W11-1i*W21;%conj(W1);
W4 = W12-1i*W22;%conj(W2);

W_c = [W1 W2 W3 W4]/sqrt(2);
