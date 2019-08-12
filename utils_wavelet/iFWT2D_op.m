function x_vec = iFWT2D_op(w_vec, g0, g1, J, sym, ROW, COL)
% Takes ifwt of multiple frames in given block
% Usage: x_vec = iWT2D_op(w_vec, g0, g1, J, sym, ROW, COL)
% x_vec - T_frames*ROW*COL vector
% w_vec - T_frames*ROW*COL vector

N = ROW*COL;
T_count = length(w_vec)/N;

x_vec = zeros(ROW*COL*T_count,1);
for frame = 1:T_count
    wf = w_vec((frame-1)*N+1:frame*N);
    % xt = IWT2_PO(reshape(wf,n,n), L, QMF);
    if isreal(wf)
        xt = ifwt2(reshape(wf,ROW,COL), g0, g1, J, sym);
    else
        xt_r = ifwt2(reshape(real(wf),ROW,COL), g0, g1, J, sym);
        xt_i = ifwt2(reshape(imag(wf),ROW,COL), g0, g1, J, sym);
        xt = xt_r+1i*xt_i;
    end
    x_vec((frame-1)*ROW*COL+1:frame*ROW*COL) = xt(:);
end
