function x_vec = At_forward_linearMC(y_vec, opts)

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

x_vec = zeros(ROW*COL*T_frames,1);
y_mr = y_vec; % This should be (T_frames-1)*n^2
switch mtype
    case 'OBMC'
        hor_OB = opts.OBMC.hor_OB;
        ver_OB = opts.OBMC.ver_OB;
        MASK = opts.OBMC.MASK;
        y_ind = 0;
        for frame = desired_frames-1
            if frame == 0
                frame = T_frames;
            end
            hor_MC = MOTION_FIELD_FORWARD(1,:,frame)';
            ver_MC = MOTION_FIELD_FORWARD(2,:,frame)';

            y_t = reshape(y_mr(y_ind+1:y_ind+N),ROW,COL);
            y_ind = y_ind+N;
            f_temp = adjoint_OBMC(y_t, hor_OB, ver_OB, hor_MC, ver_MC, MASK);
            
            x_vec((frame-1)*N+1:(frame)*N) = x_vec((frame-1)*N+1:(frame)*N)+f_temp(:);
        end
    otherwise
         y_ind = 0;
        for frame = desired_frames-1
            if frame == 0
                frame = T_frames;
            end
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
            
            y_t = reshape(y_mr(y_ind+1:y_ind+N),ROW,COL);
            y_ind = y_ind+N;
            f_temp = linear_adjoint(y_t,hor_ind, ver_ind, ROW, COL);
            
            x_vec((frame-1)*N+1:(frame)*N) = x_vec((frame-1)*N+1:(frame)*N)+f_temp(:);
        end
end
