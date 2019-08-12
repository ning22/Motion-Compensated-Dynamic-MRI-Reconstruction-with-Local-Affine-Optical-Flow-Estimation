function w_vec = DWT2D_op(x_vec, J, af, ROW, COL)
% Takes DWT of multiple frames in given block
% Usage: w_vec = DWT2D_op(x_vec, J, af, ROW, COL)
% x_vec - T_frames*ROW*COL
% w_vec - T_frames*ROW*COL vector

T_count = length(x_vec(:))/(ROW*COL);
N = ROW*COL;

w_vec = zeros(N*T_count,1);
for frame = 1:T_count
    xf = reshape(x_vec((frame-1)*N+1:frame*N),ROW,COL);
    w = dwt2D(xf, J, af);
    wt = DWTstruct2mat(w,J);
    w_vec((frame-1)*N+1:frame*N) = wt(:);
end
