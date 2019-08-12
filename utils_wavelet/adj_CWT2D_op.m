function x_vec = adj_CWT2D_op(w_vec, Faf, af, Fsf, sf, J, SYM, C2D, ROW, COL) 
% Takes adjoint of CWT with multiple frames in given block
% Usage: x_vec = adj_CWT2D_op(w_vec)
% x_vec - T_frames*N^2 vector
% w_vec - T_frames*N^2 vector

% global  Faf Fsf af sf J SYM ROW COL

% J     - number of scales
% Faf   - first stage analysis filter coeffs. 
% Fsf   - first stage synthesis filter coeffs.
% af    - analysis filter coeffs. for next stages
% sf    - synthesis filter coeffs. for next stages
% SYM   - symmetric extension options 

%n = sqrt(length(w_vec)/T_frames/4);

N = ROW*COL;
T_count = length(w_vec)/(N)/4;

switch SYM
    case 0
        x_vec = zeros(T_count*N,1);
        for frame = 1:T_count
            wf = w_vec((frame-1)*4*N+1:frame*4*N);
            if C2D
                wf = Fi_C2D(wf,ROW,COL);
            end
            xt = icplxdual2D_mod(CWTmat2struct(reshape(wf,ROW,4*COL),J),J,Fsf,sf);
            x_vec((frame-1)*N+1:frame*N) = xt(:);
        end
    case 1
        symF = 0; symh = 0; symg = 0; % Adjoint NOT WORKING PROPERLY
        x_vec = zeros(T_count*N,1);
        for frame = 1:T_count
            w_temp = w_vec((frame-1)*4*N+1:frame*4*N);
            wsym_temp = CWTmat2struct(reshape(w_temp,ROW,4*COL),J);
            wsym = cw2separate(wsym_temp,J);
            wf = CWTstruct2mat(wsym,J);
            xt = iCWT2_separate(reshape(wf,ROW,4*COL), Fsf, sf, J, symF, symh, symg);
            x_vec((frame-1)*N+1:frame*N) = xt(:);
        end
    case 2
        symF = 0; symh = 0; symg = 0;
        x_vec = zeros(T_count*N,1);
        for frame = 1:T_count
            wf = w_vec((frame-1)*4*N+1:frame*4*N);
            xt = adj_CWT2(reshape(wf,ROW,4*COL), Faf, af, J, symF, symh, symg);
            % xt = iCWT2(reshape(wf,ROW,4*COL), Fsf, sf, J, symF, symh, symg);
            x_vec((frame-1)*N+1:frame*N) = xt(:);
        end
    case 3
        symF = 1; symh = 1; symg = 2;
        x_vec = zeros(T_count*N,1);
        for frame = 1:T_count
            wf = w_vec((frame-1)*4*N+1:frame*4*N);
            if C2D
                wf = Fi_C2D(wf,ROW,COL);
            end
            if isreal(wf)
                xt = adj_CWT2(reshape(wf,ROW,4*COL), Faf, af, J, symF, symh, symg);
            else
                xt_r = adj_CWT2(reshape(real(wf),ROW,4*COL), Faf, af, J, symF, symh, symg);
                xt_i= adj_CWT2(reshape(imag(wf),ROW,4*COL), Faf, af, J, symF, symh, symg);
                xt = xt_r+1i*xt_i;
            end
            % x = xt(:);
            % x_vec = [x_vec; xt(:)];
            x_vec((frame-1)*N+1:frame*N) = xt(:);
        end
    otherwise
        disp('cant do it sire');
end
