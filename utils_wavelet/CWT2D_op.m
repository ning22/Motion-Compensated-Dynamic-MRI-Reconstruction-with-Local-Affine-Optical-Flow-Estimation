function w_vec = CWT2D_op(x_vec, Faf, af, Fsf, sf, J, SYM, C2D, ROW, COL) 
% Takes CWT of multiple frames in given block
% Usage: w_vec = CWT2D_op(x_vec)
% x_vec - T_frames*N^2 vector
% w_vec - T_frames*N^2 vector

% global ROW COL Faf Fsf af sf J SYM  C2D

% J     - number of scales
% Faf   - first stage analysis filter coeffs. 
% Fsf   - first stage synthesis filter coeffs.
% af    - analysis filter coeffs. for next stages
% sf    - synthesis filter coeffs. for next stages
% SYM   - symmetric extension options 

% A friendly way to understand/implement symmetric wavelets and adjoints 
% is to treat symmetric extension and filtering as two separate operations. 

N = ROW*COL;
T_count = length(x_vec(:))/(N);
switch SYM
    case 0
        % WORKS ONLY WITH ORTHOGONAL WAVELETS... 
        % for SYM = 0, min(ROW,COL)/2^J should not get smaller than filter length
        w_vec = zeros(N*4*T_count,1);
        for frame = 1:T_count
            xf = x_vec((frame-1)*N+1:frame*N);
            wt = CWTstruct2mat(cplxdual2D_mod(reshape(xf,ROW,COL), J,Faf, af),J);
            if C2D
                wt = F_C2D(wt,ROW,COL);
            end
            w_vec((frame-1)*N*4+1:frame*N*4) = wt(:);
        end
    case 1 
        % WORKS ONLY WITH ORTHOGONAL WAVELETS
        % case 1 is equivalent to case 0 and maintains the proper pm
        % combination, which is required in motion estimation
        symF = 0; symh = 0; symg = 0;
        w_vec = zeros(N*4*T_count,1);
        
        for frame = 1:T_count
            xf = x_vec((frame-1)*N+1:frame*N);
            w_temp = CWT2_separate(reshape(xf,ROW,COL), Faf, af, J, symF, symh, symg);
            wsym_temp = CWTmat2struct(w_temp,J);
            wsym = separate2cw(wsym_temp,J);
            wt = CWTstruct2mat(wsym, J);
            w_vec((frame-1)*N*4+1:frame*N*4) = wt(:);
        end
    case 2
        symF = 0; symh = 0; symg = 0; 
        % Use iCWT2 in place of adjoint only with orthogonal wavelets 
        w_vec = zeros(N*4*T_count,1);
        for frame = 1:T_count
            xf = x_vec((frame-1)*N+1:frame*N);
            wt = CWT2(reshape(xf,ROW,COL), Faf, af, J, symF, symh, symg);
            % w = wt(:)/sqrt(2);
            w_vec((frame-1)*N*4+1:frame*N*4) = wt(:);
        end
    case 3
        symF = 1; symh = 1; symg = 2;
        w_vec = zeros(N*4*T_count,1);
        for frame = 1:T_count
            xf = x_vec((frame-1)*N+1:frame*N);
            if isreal(xf)
                wt = CWT2(reshape(xf,ROW,COL), Faf, af, J, symF, symh, symg);
            else
                wt_r = CWT2(reshape(real(xf),ROW,COL), Faf, af, J, symF, symh, symg);
                wt_i = CWT2(reshape(imag(xf),ROW,COL), Faf, af, J, symF, symh, symg);
                wt = wt_r+1i*wt_i;
            end
            if C2D
                wt = F_C2D(wt,ROW,COL);
            end
            % w = wt(:);
            % w_vec = [w_vec; wt(:)];
            w_vec((frame-1)*N*4+1:frame*N*4) = wt(:);
        end
    otherwise
        disp('cant do it sire');
end

