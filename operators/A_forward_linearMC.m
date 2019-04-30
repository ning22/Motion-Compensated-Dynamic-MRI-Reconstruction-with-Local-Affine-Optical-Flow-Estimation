function FI = A_forward_linearMC(x_vec,opts)
% This function computes forward motion-compensated images using bilinear
% interpolation
% suppose desired_frames = [st:end];
% x_mc = [x_st; x_st+1; ... x_end] = [F_st-1 x_st-1; F_st x_st; ... F_end-1 x_end-1]
%
% Sticking with the notation that x_i+1 = F_i x_i + e_i+1

global MOTION_FIELD_FORWARD T_frames ROW COL

% Desired frames
if isfield(opts,'desired_frames')
    desired_frames = opts.desired_frames;
else
    desired_frames = 2:T_frames;
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
% MOTION_FIELD_FORWARD = opts.MOTION_FIELD_FORWARD;

N = ROW*COL;
FI = zeros(ROW*COL*(length(desired_frames(:))),1);
switch mtype
    case 'OBMC'
        hor_OB = opts.OBMC.hor_OB;
        ver_OB = opts.OBMC.ver_OB;
        MASK = opts.OBMC.MASK;
        x_ind = 0;
        for frame = desired_frames-1
            if frame == 0
                frame = T_frames;
            end
            Ic = reshape(x_vec((frame-1)*N+1:(frame)*N),ROW,COL);
            
            hor_MC = MOTION_FIELD_FORWARD(1,:,frame)';
            ver_MC = MOTION_FIELD_FORWARD(2,:,frame)';
            F_Ic = operator_OBMC(Ic, hor_OB, ver_OB, hor_MC, ver_MC, MASK);
            
            FI(x_ind+1:x_ind+ROW*COL) = F_Ic(:);
            x_ind = x_ind+ROW*COL;
        end
    otherwise
        x_ind = 0;
        for frame = desired_frames-1
            if frame == 0
                frame = T_frames;
            end
            Ic = reshape(x_vec((frame-1)*N+1:(frame)*N),ROW,COL);
            
%             [xI, yI] = meshgrid(1:COL,1:ROW);
%             locs = xI+opts.u(:,:,frame);
%             locs(locs<1) = 1; 
%             locs(locs>COL) = COL;
%             hor_ind = locs; 
% 
%             locs = yI+opts.v(:,:,frame);
%             locs(locs<1) = 1; 
%             locs(locs>ROW) = ROW;
%             ver_ind = locs;
            hor_ind = MOTION_FIELD_FORWARD{frame}{1};
            ver_ind = MOTION_FIELD_FORWARD{frame}{2};
            
            F_Ic = linear_mc(Ic,hor_ind, ver_ind, ROW, COL);
            
            FI(x_ind+1:x_ind+N) = F_Ic(:);
            x_ind = x_ind+ROW*COL;
        end    
end