function x_vec = iDWT2D_op(w_vec, J, sf, ROW, COL)
% Takes IWT of multiple frames in given block
% Usage: x_vec = IWT2D_op(w_vec, J, sf, ROW, COL)
% x_vec - T_frames*ROW*COL
% w_vec - T_frames*ROW*COL vector

N = ROW*COL;
T_count = length(w_vec(:))/(ROW*COL);

x_vec = [];
x_vec = zeros(T_count*N,1);
for frame = 1:T_count
    wf = reshape(w_vec((frame-1)*N+1:frame*N),ROW,COL);
    w = DWTmat2struct(wf,J);
    x = idwt2D(w,J,sf);
    x_vec((frame-1)*N+1:frame*N) = x(:);
end
