function x_seq = SD2seq(x_SD, desired_frames)
% Merge Static and Dynamic regions to create image sequence

global ROW COL T_frames STATIC_MASK

if nargin == 1
    desired_frames = 1:T_frames;
end

if (nnz(STATIC_MASK)>0) || (nargin>1) && sum(desired_frames(:)-[1:T_frames]')
    x_seq = zeros(ROW*COL*length(desired_frames(:)),1);
    
    static_index = find(STATIC_MASK);
    dynamic_index = find(~STATIC_MASK);
    N_S = length(static_index);
    N_D = ROW*COL-N_S;
    
    fvec_s = x_SD(1:N_S);
    fvec_d = x_SD(N_S+1:end);
    y_ind = 0;
    for frame = desired_frames
        f_frame = zeros(ROW, COL);
        f_frame(static_index) = fvec_s;
        f_frame(dynamic_index) = fvec_d((frame-1)*N_D+1:frame*N_D);
        x_seq(y_ind+1:y_ind+ROW*COL) = f_frame(:);
        y_ind = y_ind+ROW*COL;
    end
else
    x_seq = x_SD;
end