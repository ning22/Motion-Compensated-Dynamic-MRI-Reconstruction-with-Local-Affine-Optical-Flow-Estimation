function x_vec = At_backward_linearMC(y_vec,opts)

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
x_vec = zeros(ROW*COL*T_frames,1);
y_mr = y_vec; % This should be (T_frames-1)*n^2
switch mtype
    case 'OBMC'
        hor_OB = opts.OBMC.hor_OB;
        ver_OB = opts.OBMC.ver_OB;
        MASK = opts.OBMC.MASK;
        y_ind = 0;
        for frame = desired_frames+1
            if frame == T_frames+1 % assuming periodicity of image seq, last one related to the first
                frame = 1;
            end
            hor_MC = MOTION_FIELD_BACKWARD(1,:,frame)';
            ver_MC = MOTION_FIELD_BACKWARD(2,:,frame)';
            
            y_t = reshape(y_mr(y_ind+1:y_ind+N),ROW,COL);
            y_ind = y_ind + N;
            f_temp = adjoint_OBMC(y_t, hor_OB, ver_OB, hor_MC, ver_MC, MASK);
            
            x_vec((frame-1)*N+1:(frame)*N) = x_vec((frame-1)*N+1:(frame)*N)+f_temp(:);
        end
    otherwise
        y_ind = 0;
        for frame = desired_frames+1
            if frame == T_frames+1 % assuming periodicity of image seq, last one related to the first
                frame = 1;
            end
            hor_ind = MOTION_FIELD_BACKWARD{frame}{1};
            ver_ind = MOTION_FIELD_BACKWARD{frame}{2};
            
            y_t = reshape(y_mr(y_ind+1:y_ind+N),ROW,COL);
            y_ind = y_ind + N;
            f_temp = linear_adjoint(y_t,hor_ind, ver_ind, ROW, COL);
            
            x_vec((frame-1)*N+1:(frame)*N) = x_vec((frame-1)*N+1:(frame)*N)+f_temp(:);
        end
end
