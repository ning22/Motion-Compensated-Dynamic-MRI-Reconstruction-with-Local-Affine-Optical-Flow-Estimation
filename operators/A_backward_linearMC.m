function BI = A_backward_linearMC(x_vec, opts)
% This function computes backward motion-compensated images using bilinear
% interpolation
% suppose desired_frames = [st:end];
% x_mc = [x_st; x_st+1; ... x_end] = [B_st+1 x_st+1; F_st+2 x_st+2; ... F_end+1 x_end+1]
%
% Sticking with the notation that x_i-1 = B_i x_i + e_i-1

global MOTION_FIELD_BACKWARD T_frames ROW COL

% Desired frames
if isfield(opts,'desired_frames')
    desired_frames = opts.desired_frames;
else
    desired_frames = 1:T_frames-1;
end
% Motion type
if isfield(opts,'mtype')
    mtype = opts.mtype;
else
    mtype = 'CWT';
end
% ROW = opts.ROW;
% COL = opts.COL;
% T_frames = opts.T_frames;
% MOTION_FIELD_BACKWARD = opts.MOTION_FIELD_BACKWARD;

N = ROW*COL;
BI = zeros(ROW*COL*(length(desired_frames(:))),1);
switch mtype
   case 'OBMC'
        hor_OB = opts.OBMC.hor_OB;
        ver_OB = opts.OBMC.ver_OB;
        MASK = opts.OBMC.MASK;
        x_ind = 0;
        for frame = desired_frames+1;
            if frame == T_frames+1 % assuming periodicity of image seq, last one related to the first
                frame = 1;
            end
            Ic = reshape(x_vec((frame-1)*N+1:frame*N),ROW,COL);
            
            hor_MC = MOTION_FIELD_BACKWARD(1,:,frame)';
            ver_MC = MOTION_FIELD_BACKWARD(2,:,frame)';
            B_Ic = operator_OBMC(Ic, hor_OB, ver_OB, hor_MC, ver_MC, MASK);
            
            BI(x_ind+1:x_ind+ROW*COL) = B_Ic(:);
            x_ind = x_ind+ROW*COL;
        end
    otherwise
        x_ind = 0;
        for frame = desired_frames+1;
            if frame == T_frames+1 % assuming periodicity of image seq, last one related to the first
                frame = 1;
            end
            Ic = reshape(x_vec((frame-1)*N+1:frame*N),ROW,COL);
            
            hor_ind = MOTION_FIELD_BACKWARD{frame}{1};
            ver_ind = MOTION_FIELD_BACKWARD{frame}{2};
            B_Ic = linear_mc(Ic,hor_ind, ver_ind, ROW, COL);
            
            BI(x_ind+1:x_ind+N) = B_Ic(:);
            x_ind = x_ind+N;
        end
    
end