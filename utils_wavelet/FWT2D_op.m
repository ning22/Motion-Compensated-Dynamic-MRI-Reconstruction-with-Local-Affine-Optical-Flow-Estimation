function w_vec = FWT2D_op(x_vec, h0, h1, J, sym, ROW, COL) 

% Takes fwt2 of multiple frames in given block
% Usage: w_vec = FWT2D_op(x_vec, J, h0, h1, sym, ROW, COL)
% x_vec - T_frames*ROW*COL vector
% w_vec - T_frames*ROW*COL vector

T_count = length(x_vec(:))/(ROW*COL);
N = ROW*COL;

w_vec = zeros(ROW*COL*T_count,1);
for frame = 1:T_count
    xf = x_vec((frame-1)*N+1:frame*N);
    if isreal(xf)
        wt = fwt2(reshape(xf,ROW,COL), h0, h1, J, sym);
    else
        wt_r = fwt2(reshape(real(xf),ROW,COL), h0, h1, J, sym);
        wt_i = fwt2(reshape(imag(xf),ROW,COL), h0, h1, J, sym);
        wt = wt_r+1i*wt_i;
    end    
    w_vec((frame-1)*N+1:frame*N) = wt(:);
end
