function x_SD = seq2SD(x_seq, desired_frames)
% Divides a sequence of images to Static and Dynamic regions

global ROW COL T_frames STATIC_MASK

if nargin == 1
    desired_frames = 1:T_frames;
end

if (nnz(STATIC_MASK)>0) || (nargin>1) && sum(desired_frames(:)-[1:T_frames]')
    static_index = find(STATIC_MASK);
    dynamic_index = ~STATIC_MASK;
    N_S = length(static_index);
    N_D = ROW*COL-N_S;
    N = ROW*COL;
    
    fvec_s = zeros(N_S,1);
    fvec_d = zeros(N_D*T_frames,1);
    y_ind = 0;
    for frame = desired_frames
        f_frame = x_seq(y_ind+1:y_ind+N);
        y_ind = y_ind+N;
        f_s = f_frame(static_index);
        f_d = f_frame(dynamic_index);
        
        fvec_s = fvec_s+f_s;
        fvec_d((frame-1)*N_D+1:frame*N_D,:) = f_d;
    end
    x_SD = [fvec_s; fvec_d];
else
    x_SD = x_seq;
end