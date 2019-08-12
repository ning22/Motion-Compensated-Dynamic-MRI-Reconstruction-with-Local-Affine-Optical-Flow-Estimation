function x_vec = adj_FWT2D_op(w_vec, h0, h1, J, sym, ROW, COL)
% Takes adjoint for FWT2D_op 
% Usage: x_vec = adj_FWT2D_op(w_vec, J, h0, h1, sym, ROW, COL)
% x_vec - T_frames*ROW*COL vector
% w_vec - T_frames*ROW*COL vector

N = ROW*COL;
T_count = length(w_vec)/N;

x_vec = zeros(ROW*COL*T_count,1);
for frame = 1:T_count
    wf = w_vec((frame-1)*N+1:frame*N);
    if isreal(wf)
        xt = afwt2(reshape(wf,ROW,COL), h0, h1, J, sym);
    else
        xt_r = afwt2(reshape(real(wf),ROW,COL), h0, h1, J, sym);
        xt_i= afwt2(reshape(imag(wf),ROW,COL), h0, h1, J, sym);
        xt = xt_r+1i*xt_i;
    end
    
    x_vec((frame-1)*ROW*COL+1:frame*ROW*COL) = xt(:);
end

