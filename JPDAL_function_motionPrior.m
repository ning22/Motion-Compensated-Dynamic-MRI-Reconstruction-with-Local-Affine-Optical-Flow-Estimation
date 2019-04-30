function [Ir_cube, VOF, costs] = JPDAL_function_motionPrior(opts)

% operators and initialization
b = opts.observation;
A = opts.A;
At = opts.At;
xinit = opts.xinit;
xtrue = opts.xtrue;
periodic = opts.periodic;
if periodic 
    I0 = circshift(xinit,1,3); % periodic
else
    I0 = xinit;
    I0(:,:,2:end) = opts.xinit(:,:,1:end-1);
end
% save_results = opts.save_results;
% mc_option = opts.mc_option; 
global MOTION_FIELD_FORWARD ROW COL T_frames

% algorithm tuning parameters
maxIter = opts.maxIter;
t = opts.stepsize;
tol = opts.tol;
showTrigger = opts.showTrigger;

delta = 0.99;
alpha = 1e-1; 
reductionFactor = 0.5;
theta = 1;

% regularization parameters
tau = opts.rpOF;
eta = opts.rpSP;
etas = opts.rpLR;
gammap = opts.rpTV;
prior_x = opts.priorType;
prior_motion = opts.priorMotion;
T = opts.T;
U = @(x) vec(opts.U_h(x));
Ut = @(x) cube(opts.Ut_h(x));

scaleSet = sort(opts.scaleSet,'descend');
numscale = numel(scaleSet);
J_coarsest = scaleSet(1);
bspline = 'b5';
[bfilt,wfilt1] = build_b_pyramid(J_coarsest,bspline);
w0 = bfilt{1}; 
w1 = wfilt1{1};
N = str2double(bspline(2)); % bspline order
H = two_scale_filters_computation(J_coarsest,N); % filters for two scale filters computation
opts.H = H;
opts.w0 = w0;
opts.w1 = w1;
[ROW,COL,T_frames] = size(I0);
Ntotal = numel(I0);
[xI, yI, zI] = meshgrid(1:COL,1:ROW,1:T_frames);
uJ = zeros(ROW,COL,T_frames); % 
vJ = uJ;
interp_type = 'linear';

cube = @(x) reshape(x,ROW,COL,T_frames);
Ir_set = cell(3,1);
opts.method = 'PDAL';
costs = [];
costsLL = [];  costsOF = []; costsRGX = [];  costsRGVOF=[];
%% Reconstruction...
tic
for iscale= 1:numscale
    nscale = scaleSet(iscale);
    dd = 2^(nscale-1);
    nr = ceil(ROW/dd);
    nc = ceil(COL/dd);   
    
    if iscale == 1
        x_km1 = xinit;
        u_km1 = zeros(nr,nc,T_frames);
        v_km1 = zeros(nr,nc,T_frames);
        ux_km1 = zeros(nr,nc,T_frames);
        vx_km1 = zeros(nr,nc,T_frames);
        uy_km1 = zeros(nr,nc,T_frames);
        vy_km1 = zeros(nr,nc,T_frames);       
    else
        x_km1 = x_k;
        u_km1 = uJ(1:dd:end,1:dd:end,:);    
        v_km1 = vJ(1:dd:end,1:dd:end,:);
        ux_km1 = uxJ(1:dd:end,1:dd:end,:);    
        vx_km1 = vxJ(1:dd:end,1:dd:end,:);
        uy_km1 = uyJ(1:dd:end,1:dd:end,:);    
        vy_km1 = vyJ(1:dd:end,1:dd:end,:);

%         tol = tol; % lungcor
        tol = tol*0.5;  %cardiac
%          tol = tol*0.05; %wrap
        tau = tau*opts.ka; %comment lungcor
        eta = eta*opts.ka;
        etas = etas*opts.ka;
%         opts.ma = opts.ma*0.8; % cardiac
        opts.ma = opts.ma*0.5; % wrap
    end
    fprintf('No. scale: %3d ~ Regularization of OF: %.2e  \n', nscale, tau);
    tkm1 = t;
    z1_km1 = zeros(size(b));
    z2_km1 = zeros(nr,nc,T_frames);
    switch prior_x
        case 'l1_w'
            z3_km1 = zeros(4*Ntotal,1);
        case 'l1'
            z3_km1 = zeros(ROW,COL,T_frames);
        case {'TV','tv'}
            z3_km1 = zeros(ROW,COL,2*T_frames);
        case {'tv3d','TV3d'}
            z3_km1 = zeros(ROW,COL,3*T_frames);
        case {'nuclear','low_rank'}
            z3_km1 = zeros(ROW,COL,T_frames);
        case 'l+s'
            z3_km1 = zeros(ROW,COL,T_frames);
            z3s_km1 = zeros(ROW,COL,T_frames);
        case 'l+tv'
            z3_km1 = zeros(ROW,COL,2*T_frames);
            z3s_km1 = zeros(ROW,COL,T_frames);        
        case 'l+tv3d'
            z3_km1 = zeros(ROW,COL,3*T_frames);
            z3s_km1 = zeros(ROW,COL,T_frames);
        case 'l1+tv3d'
            z3_km1 = zeros(ROW,COL,3*T_frames);
            z3s_km1 = zeros(ROW,COL,T_frames);            
    end
    % prior of the motion fields
    switch prior_motion
        case {'TV','tv'}
            z4_km1 = zeros(nr,nc,2*T_frames);
            z5_km1 = zeros(nr,nc,2*T_frames);
            z4x_km1 = zeros(nr,nc,2*T_frames);
            z5x_km1 = zeros(nr,nc,2*T_frames);
            z4y_km1 = zeros(nr,nc,2*T_frames);
            z5y_km1 = zeros(nr,nc,2*T_frames);
        case 'l2'
            z4_km1 = zeros(nr,nc,T_frames);
            z5_km1 = zeros(nr,nc,T_frames);
            z4x_km1 = zeros(nr,nc,T_frames);
            z5x_km1 = zeros(nr,nc,T_frames);
            z4y_km1 = zeros(nr,nc,T_frames);
            z5y_km1 = zeros(nr,nc,T_frames);            
        case 'l1'
            z4_km1 = zeros(nr,nc,T_frames);
            z5_km1 = zeros(nr,nc,T_frames);
            z4x_km1 = zeros(nr,nc,T_frames);
            z5x_km1 = zeros(nr,nc,T_frames);
            z4y_km1 = zeros(nr,nc,T_frames);
            z5y_km1 = zeros(nr,nc,T_frames);               
    end
    
    [D1I0,D2I0] = gradient_op(abs(I0));
    [WD1I0,MW] = ApplyW(D1I0,nscale,H,w0,w1,0);
    [WD2I0] = ApplyW(D2I0,nscale,H,w0,w1,0);
    WI0 = ApplyW(abs(I0),nscale,H,w0,w1,1);
    
    applyI0hat = @(u0,ux,uy,v0,vx,vy) WD1I0.m_00.*u0 + WD1I0.m_01.*ux + WD1I0.m_10.*uy + ...
        WD2I0.m_00.*v0 + WD2I0.m_01.*vx +WD2I0.m_10.*vy;
    applyI0hatTrans_u0 = @(x) WD1I0.m_00.*x;
    applyI0hatTrans_ux = @(x) WD1I0.m_01.*x;
    applyI0hatTrans_uy = @(x) WD1I0.m_10.*x;
    applyI0hatTrans_v0 = @(x) WD2I0.m_00.*x;
    applyI0hatTrans_vx = @(x) WD2I0.m_01.*x;
    applyI0hatTrans_vy = @(x) WD2I0.m_10.*x;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for k = 0:maxIter % The Chambolle-Pock iteration begins here.
        switch prior_x
            case 'l1_w'
                Ktz1 = At(z1_km1) - ApplyWTrans(z2_km1,nscale,H,w0,MW) + Ut(z3_km1);
            case 'l1'
                Ktz1 = At(z1_km1) - ApplyWTrans(z2_km1,nscale,H,w0,MW) + T'*z3_km1;
            case {'TV','tv'}
                Ktz1 = At(z1_km1) - ApplyWTrans(z2_km1,nscale,H,w0,MW) + applyDTrans(z3_km1,T_frames);
            case {'tv3d','TV3d'}
                Ktz1 = At(z1_km1) - ApplyWTrans(z2_km1,nscale,H,w0,MW) + applyD3dTrans(z3_km1,T_frames);
            case {'nuclear','low_rank'}
                Ktz1 = At(z1_km1) - ApplyWTrans(z2_km1,nscale,H,w0,MW) + z3_km1;
            case 'l+s'
                Ktz1 = At(z1_km1) - ApplyWTrans(z2_km1,nscale,H,w0,MW) + T'*z3_km1 + z3s_km1;
            case 'l+tv'
                Ktz1 = At(z1_km1) - ApplyWTrans(z2_km1,nscale,H,w0,MW) + applyDTrans(z3_km1,T_frames)+ z3s_km1;
            case 'l+tv3d'
                Ktz1 = At(z1_km1) - ApplyWTrans(z2_km1,nscale,H,w0,MW) + applyD3dTrans(z3_km1,T_frames) + z3s_km1;
            case 'l1+tv3d'
                Ktz1 = At(z1_km1) - ApplyWTrans(z2_km1,nscale,H,w0,MW) + applyD3dTrans(z3_km1,T_frames) + T'*z3s_km1;
        end
        
        switch prior_motion
            case 'TV'
                Ktz2 = applyI0hatTrans_u0(z2_km1) + applyDTrans(z4_km1,T_frames);
                Ktz3 = applyI0hatTrans_v0(z2_km1) + applyDTrans(z5_km1,T_frames);
                Ktz2x = applyI0hatTrans_ux(z2_km1) + applyDTrans(z4x_km1,T_frames);
                Ktz3x = applyI0hatTrans_vx(z2_km1) + applyDTrans(z5x_km1,T_frames);
                Ktz2y = applyI0hatTrans_uy(z2_km1) + applyDTrans(z4y_km1,T_frames);
                Ktz3y = applyI0hatTrans_vy(z2_km1) + applyDTrans(z5y_km1,T_frames);
            case {'l2', 'l1'}
                Ktz2 = applyI0hatTrans_u0(z2_km1) + z4_km1;
                Ktz3 = applyI0hatTrans_v0(z2_km1) + z5_km1;
                Ktz2x = applyI0hatTrans_ux(z2_km1) + z4x_km1;
                Ktz3x = applyI0hatTrans_vx(z2_km1) + z5x_km1;
                Ktz2y = applyI0hatTrans_uy(z2_km1) + z4y_km1;
                Ktz3y = applyI0hatTrans_vy(z2_km1) + z5y_km1;                               
        end
        x_k = x_km1 - tkm1*Ktz1;
        u_k = u_km1 - tkm1*Ktz2;
        v_k = v_km1 - tkm1*Ktz3;
        ux_k = ux_km1 - tkm1*Ktz2x;
        vx_k = vx_km1 - tkm1*Ktz3x;
        uy_k = uy_km1 - tkm1*Ktz2y;
        vy_k = vy_km1 - tkm1*Ktz3y;    

        % calculate the cost function
        WIt = ApplyW(I0-abs(x_k),nscale,H,w0,w1,1);
        [cost,costLL,costOF,costRGX,costRGVOF] = computeCost...
            (A,b,x_k,u_k,v_k,ux_k,vx_k,uy_k,vy_k,applyI0hat,WIt,T,tau,eta,etas,gammap,prior_x,prior_motion,ROW,COL,T_frames);
        costs = [costs,cost];
        costsLL = [costsLL,costLL];
        costsOF = [costsOF,costOF];
        costsRGX = [costsRGX,costRGX];
        costsRGVOF = [costsRGVOF,costRGVOF];
               
        t = tkm1*sqrt(1+theta);
        accept_t = 0;      

        if k>2
            if mod(k,showTrigger) == 0           
                if numel(opts.xtrue)>1
                    rmse = calRMSE(x_k(:),xtrue(:));
                    errX = opts.errFcn(x_k);
                    fprintf('Iter: %3d  ~ Cost: %.3e ~ t: %.2e ~ ~ RMSE:  %.2e ~ ~ Error: %.2e, stopping rule: %.2e\n',k,cost,t,rmse,errX, abs(costs(k)-costs(k-1))); 
                else
                    fprintf(' ite: %d , cost: %f3, update: %f3\n', k,cost,norm(x_k(:)-x_km1(:))/norm(x_km1(:))); 
                end   
            end
            
            
            if abs(costs(k)-costs(k-1))<tol
                break
            end      
        end

        % linesearch
        while accept_t ==0
            theta = t/tkm1;
            s = alpha*t;

            % proximal operator of g1* (data fidelity term A(x0)-b_kspace)    
            term1 = x_k+theta*(x_k-x_km1);
            in_1 = z1_km1 + s*A(term1);
            z1_k = (in_1 - s*b)/(1 + s);

            % proximal operator of g2* (data fidelity term I-func(I0))                          
            term2 = u_k+theta*(u_k-u_km1);
            term3 = v_k+theta*(v_k-v_km1);
            term2x = ux_k+theta*(ux_k-ux_km1);
            term3x = vx_k+theta*(vx_k-vx_km1);
            term2y = uy_k+theta*(uy_k-uy_km1);
            term3y = vy_k+theta*(vy_k-vy_km1);        
            in_2 = z2_km1 + s*(-ApplyW(abs(term1),nscale,H,w0,w1,1) + applyI0hat(term2,term2x,term2y,term3,term3x,term3y)+WI0);   
%             z2_k = in_2./max(1,abs(in_2)/tau); % l1 term
            z2_k = tau*in_2/(tau + s); % l2 term

            % proximal operator of g3* (regularization of image x)
            switch prior_x
                case 'l1'
                    in_3 = z3_km1 + s*(T*term1);
                    z3_k = in_3./max(1,abs(in_3)/eta);
                case {'TV','tv'}
                    in_3 = z3_km1 + s*applyD(term1);
                    z3_k = in_3./max(1,abs(in_3)/eta);
                case {'tv3d','TV3d'}
                    in_3 = z3_km1 + s*applyD3d(term1);
                    z3_k = in_3./max(1,abs(in_3)/eta);            
                case {'nuclear','low_rank'}
                    in_3 = z3_km1 + s*term1;
                    [temp, St] = proxNuclearnorm(in_3/s,ROW,COL,T_frames,eta/s);
                    z3_k = in_3 - s*temp;
                case 'l+s'
                    in_3 = z3_km1 + s*(T*term1);
                    z3_k = in_3./max(1,abs(in_3)/eta);

                    in_3s = z3s_km1 + s*term1;
                    [temp,St] = proxNuclearnorm(in_3s/s,ROW,COL,T_frames,etas/s);
                    z3s_k = in_3s - s*temp;
                case 'l+tv'
                    in_3 = z3_km1 + s*applyD(term1);
                    z3_k = in_3./max(1,abs(in_3)/eta);

                    in_3s = z3s_km1 + s*term1;
                    [temp,St] = proxNuclearnorm(in_3s/s,ROW,COL,T_frames,etas/s);
                    z3s_k = in_3s - s*temp;                  
                case 'l+tv3d'
                    in_3 = z3_km1 + s*applyD3d(term1);
                    z3_k = in_3./max(1,abs(in_3)/eta);

                    in_3s = z3s_km1 + s*term1;
                    [temp,St] = proxNuclearnorm(in_3s/s,ROW,COL,T_frames,etas/s);
                    z3s_k = in_3s - s*temp;    
                case 'l1+tv3d'
                    in_3 = z3_km1 + s*applyD3d(term1);
                    z3_k = in_3./max(1,abs(in_3)/eta);

                    in_3s = z3s_km1 + s*(T*term1);
                    z3s_k = in_3s./max(1,abs(in_3s)/etas);
            end

            % proximal operator of g4* -- g10* (regularization of the motion vectors)
            switch prior_motion
                case 'TV'
                    in_4 = z4_km1 +s*applyD(term2);
                    in_5 = z5_km1 +s*applyD(term3);
                    in_4x = z4x_km1 +s*applyD(term2x);
                    in_5x = z5x_km1 +s*applyD(term3x);
                    in_4y = z4y_km1 +s*applyD(term2y);
                    in_5y = z5y_km1 +s*applyD(term3y);        
                    z4_k = in_4./max(1,abs(in_4)/gammap);
                    z5_k = in_5./max(1,abs(in_5)/gammap);
                    z4x_k = in_4x./max(1,abs(in_4x)/gammap);
                    z5x_k = in_5x./max(1,abs(in_5x)/gammap);
                    z4y_k = in_4y./max(1,abs(in_4y)/gammap);
                    z5y_k = in_5y./max(1,abs(in_5y)/gammap);
                case 'l2'            
                    z4_k = gammap*(z4_km1 + s*term2)/(gammap+s);
                    z5_k = gammap*(z5_km1 + s*term3)/(gammap+s);
                    z4x_k = gammap*(z4x_km1 + s*term2x)/(gammap+s);
                    z4y_k = gammap*(z4y_km1 + s*term2y)/(gammap+s);
                    z5x_k = gammap*(z5x_km1 + s*term3x)/(gammap+s);
                    z5y_k = gammap*(z5y_km1 + s*term3y)/(gammap+s);
                case 'l1'
                    in_4 = z4_km1 + s*term2;
                    in_5 = z5_km1 + s*term3;
                    in_4x = z4x_km1 + s*term2x;
                    in_4y = z4y_km1 + s*term2y;
                    in_5x = z5x_km1 + s*term3x;
                    in_5y = z5y_km1 + s*term3y;
                    z4_k = in_4./max(1,abs(in_4)/gammap);                   
                    z5_k = in_5./max(1,abs(in_5)/gammap);                    
                    z4x_k = in_4x./max(1,abs(in_4x)/gammap);                    
                    z4y_k = in_4y./max(1,abs(in_4y)/gammap);                   
                    z5x_k = in_5x./max(1,abs(in_5x)/gammap);                    
                    z5y_k = in_5y./max(1,abs(in_5y)/gammap);
            end
            
            switch prior_x
                case 'l1_w'
                    Ktz1_new = At(z1_k) - ApplyWTrans(z2_k,nscale,H,w0,MW) + Ut(z3_k);
                case 'l1'
                    Ktz1_new = At(z1_k) - ApplyWTrans(z2_k,nscale,H,w0,MW) + T'*z3_k;
                case {'TV','tv'}
                    Ktz1_new = At(z1_k) - ApplyWTrans(z2_k,nscale,H,w0,MW) + applyDTrans(z3_k,T_frames); 
                case {'tv3d','TV3d'}
                    Ktz1_new = At(z1_k) - ApplyWTrans(z2_k,nscale,H,w0,MW) + applyD3dTrans(z3_k,T_frames); 
                case {'nuclear','low_rank'}
                    Ktz1_new = At(z1_k) - ApplyWTrans(z2_k,nscale,H,w0,MW) + z3_k; 
                case 'l+s'
                    Ktz1_new = At(z1_k) - ApplyWTrans(z2_k,nscale,H,w0,MW) + T'*z3_k + z3s_k;
                case 'l+tv'
                    Ktz1_new = At(z1_k) - ApplyWTrans(z2_k,nscale,H,w0,MW) + applyDTrans(z3_k,T_frames) + z3s_k;
                case 'l+tv3d'
                    Ktz1_new = At(z1_k) - ApplyWTrans(z2_k,nscale,H,w0,MW) + applyD3dTrans(z3_k,T_frames) + z3s_k;
                case 'l1+tv3d'
                    Ktz1_new = At(z1_k) - ApplyWTrans(z2_k,nscale,H,w0,MW) + applyD3dTrans(z3_k,T_frames) + T'*z3s_k;
            end
            
            switch prior_motion
                case 'TV'
                    Ktz2_new = applyI0hatTrans_u0(z2_k) + applyDTrans(z4_k,T_frames);
                    Ktz3_new = applyI0hatTrans_v0(z2_k) + applyDTrans(z5_k,T_frames);
                    Ktz2x_new = applyI0hatTrans_ux(z2_k) + applyDTrans(z4x_k,T_frames);
                    Ktz3x_new = applyI0hatTrans_vx(z2_k) + applyDTrans(z5x_k,T_frames);   
                    Ktz2y_new = applyI0hatTrans_uy(z2_k) + applyDTrans(z4y_k,T_frames);
                    Ktz3y_new = applyI0hatTrans_vy(z2_k) + applyDTrans(z5y_k,T_frames);        
                case {'l2', 'l1'}
                    Ktz2_new = applyI0hatTrans_u0(z2_k) + z4_k;
                    Ktz3_new = applyI0hatTrans_v0(z2_k) + z5_k;
                    Ktz2x_new = applyI0hatTrans_ux(z2_k) + z4x_k;
                    Ktz3x_new = applyI0hatTrans_vx(z2_k) + z5x_k;   
                    Ktz2y_new = applyI0hatTrans_uy(z2_k) + z4y_k;
                    Ktz3y_new = applyI0hatTrans_vy(z2_k) + z5y_k;                      
            end
            Ktz_diff1 = Ktz1_new - Ktz1;
            Ktz_diff2 = Ktz2_new - Ktz2;
            Ktz_diff3 = Ktz3_new - Ktz3;
            Ktz_diff2x = Ktz2x_new - Ktz2x;
            Ktz_diff3x = Ktz3x_new - Ktz3x;
            Ktz_diff2y = Ktz2y_new - Ktz2y;
            Ktz_diff3y = Ktz3y_new - Ktz3y;
            z_diff1 = z1_k - z1_km1;
            z_diff2 = z2_k - z2_km1;
            z_diff3 = z3_k - z3_km1;
            z_diff4 = z4_k - z4_km1;
            z_diff5 = z5_k - z5_km1;
            z_diff4x = z4x_k - z4x_km1;
            z_diff5x = z5x_k - z5x_km1;
            z_diff4y = z4y_k - z4y_km1;
            z_diff5y = z5y_k - z5y_km1;

            lhs = t*s*(norm(Ktz_diff1(:))^2 + norm(Ktz_diff2(:))^2 + norm(Ktz_diff3(:))^2+ norm(Ktz_diff2x(:))^2 + norm(Ktz_diff3x(:))^2+ ...
                norm(Ktz_diff2y(:))^2 + norm(Ktz_diff3y(:))^2);
            if strcmp(prior_x,'l+s') || strcmp(prior_x,'l+tv')|| strcmp(prior_x,'l+tv3d')|| strcmp(prior_x,'l1+tv3d')
                z_diff3s= z3s_k -z3s_km1;
                rhs = delta^2*(norm(z_diff1(:))^2 + norm(z_diff2(:))^2 + norm(z_diff3(:))^2 +norm(z_diff3s(:))^2 + norm(z_diff4(:))^2 + norm(z_diff5(:))^2+...
                    norm(z_diff4x(:))^2 + norm(z_diff5x(:))^2+ norm(z_diff4y(:))^2 + norm(z_diff5y(:))^2);
            else
                rhs = delta^2*(norm(z_diff1(:))^2 + norm(z_diff2(:))^2 + norm(z_diff3(:))^2 + norm(z_diff4(:))^2 + norm(z_diff5(:))^2 + norm(z_diff4x(:))^2 ...
                    + norm(z_diff5x(:))^2+ norm(z_diff4y(:))^2 + norm(z_diff5y(:))^2);    
            end
            if lhs <= rhs
                tkm1 = t;
                accept_t = 1;
            else
                t = reductionFactor*t;
                tkm1 = t;
            end 
        end

%         if iscale >1 && k>=1
%             if errX<errXp
%                 if periodic 
%                     I0 = circshift(x_k,1,3); % periodic
%                 else
%                     I0 = x_k;
%                     I0(:,:,2:end) = x_k(:,:,1:end-1);
%                 end 
%             end
%         else
            if periodic 
                I0 = circshift(x_k,1,3); % periodic
            else
                I0 = x_k;
                I0(:,:,2:end) = x_k(:,:,1:end-1);
            end     

%         end
        
        x_km1 = x_k;
        u_km1 = u_k;
        v_km1 = v_k;
        ux_km1 = ux_k;
        vx_km1 = vx_k;
        uy_km1 = uy_k;
        vy_km1 = vy_k;
        z1_km1 = z1_k;
        z2_km1 = z2_k;
        z3_km1 = z3_k;
        z4_km1 = z4_k;
        z5_km1 = z5_k;
        z4x_km1 = z4x_k;
        z5x_km1 = z5x_k;
        z4y_km1 = z4y_k;
        z5y_km1 = z5y_k;
        if strcmp(prior_x,'l+s')|| strcmp(prior_x,'l+tv3d') || strcmp(prior_x,'l1+tv3d')
            z3s_km1 = z3s_k;       
        end 
    end      
    if nscale ==1
        uJ = u_k;
        vJ = v_k;
        uxJ = ux_k;
        vxJ = vx_k;
        uyJ = uy_k;
        vyJ = vy_k;
    else
        xx = xI(1:dd:end,1:dd:end,:);
        yy = yI(1:dd:end,1:dd:end,:);
        zz = zI(1:dd:end,1:dd:end,:);

        uJ = interp3(xx,yy,zz,u_k,xI,yI,zI,interp_type,0);
        vJ = interp3(xx,yy,zz,v_k,xI,yI,zI,interp_type,0);
        uxJ = interp3(xx,yy,zz,ux_k,xI,yI,zI,interp_type,0);
        vxJ = interp3(xx,yy,zz,vx_k,xI,yI,zI,interp_type,0);
        uyJ = interp3(xx,yy,zz,uy_k,xI,yI,zI,interp_type,0);
        vyJ = interp3(xx,yy,zz,vy_k,xI,yI,zI,interp_type,0);   
    end 
    
    Ir_set{iscale} = x_k;
    VOF.uJ{iscale} = uJ;
    VOF.vJ{iscale} = vJ;
    
    if numel(opts.xtrue) ==1
        xinit_roi = opts.xinit(opts.roi_hor,opts.roi_ver,:);
        xk_roi = x_k(opts.roi_hor,opts.roi_ver,:);
        rg{iscale} = Calc_RG(xinit_roi,xk_roi);
        fprintf('RG Min: %.4f, RG Max: %.4f\n', min(rg{iscale}), max(rg{iscale}));
    end
    
    %% motion compensation
%     % Liu method for motion estimation
%     alpha = 1;
%     ratio = 0.5;
%     minWidth = 40;
%     nOuterFPIterations = 3;
%     nInnerFPIterations = 1;
%     nSORIterations = 20;
%     para = [alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations];
%     opts.para = para;    
% 
%     Ir_cube_me = x_k;
%     uF = zeros(size(x_k));
%     vF = uF;
% 
%     for frame = 1:T_frames
%         if frame ==1
%             Ip_m = Ir_cube_me(:,:,T_frames);
%         else
%             Ip_m = Ir_cube_me(:,:,frame-1);
%         end
%         Ie = Ir_cube_me(:,:,frame); % TARGET FRAME
% 
%         Ip_m = abs(Ip_m); 
%         Ie = abs(Ie);
% 
%         [u,v] = handle_OF(Ip_m,Ie,'Liu',opts);
%         uF(:,:,frame) = u;
%         vF(:,:,frame) = v;
%     end
%     opts.uJ = uF;
%     opts.vJ = vF;  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     load('./results/Ref_motion_LC.mat');
%     opts.uJ = uR;
%     opts.vJ = vR;

    opts.stepsize = 10;  
    opts.uJ = uJ;
    opts.vJ = vJ;
    opts_motion = [];
    
    opts_motion.mtype = 'OF';
    opts_motion.desired_frames = 1:T_frames;
    opts.U = @(z) cube(A_forward_linearMC(SD2seq(z),opts_motion)-SD2seq(z,opts_motion.desired_frames));   
    opts.Ut = @(z) cube(seq2SD(At_forward_linearMC(z(:),opts_motion))-seq2SD(z(:),opts_motion.desired_frames));   
    opts.maxIter_mc = 200;
    if strcmp(opts.mname,'Cardiac')
        opts.maxIter_mc = 50;
    end
    [Ir_cube_mc,costs_mc] = MC_function(A,At,x_k,b,opts);
    costs = [costs,costs_mc];
    if numel(opts.xtrue) ==1
        Ir_mc_roi = Ir_cube_mc(opts.roi_hor,opts.roi_ver,:);
        rg_mc = Calc_RG(xinit_roi,Ir_mc_roi);
        fprintf('RG Min: %.4f, RG Max: %.4f\n', min(rg_mc), max(rg_mc));
        figure(10), plot(rg{iscale}); hold on, plot(rg_mc,'r'); hold off;
    end
%     
    x_k = Ir_cube_mc;

end
timeTake = toc;
fprintf('Computational time: %.2f, prior is %s\n', timeTake/60, prior_x);
% figure, semilogy(costs);
Ir_cube = x_k;
end

function [cost,cost1,cost2,cost3,cost4] = computeCost(A,b,x,u,v,ux,vx,uy,vy,ApplyI0hat,WIt,T,tau,eta,etas,gamma, prior_x,prior_motion,nr,nc,nf) 
term1 = A(x)-b;
cost1 = .5*sum(abs(term1(:)).^2);

term2 = ApplyI0hat(u,ux,uy,v,vx,vy)+WIt;
% cost2 = tau*sum(abs(term2(:))); %l1 term
cost2 = tau*sum(abs(term2(:)).^2); %l2 term

if eta~=0 || etas~=0
switch prior_x
    case 'l1'
        term3 = T*x;
        cost3 = eta*sum(abs(term3(:)));
    case {'TV','tv'}
        term3 = applyD(x);
        cost3 = eta*sum(abs(term3(:)));
    case {'tv3d','TV3d'}
        term3 = applyD3d(x);
        cost3 = eta*sum(abs(term3(:)));        
    case {'nuclear','low_rank'}
        [Ut,St,Vt]=svd(reshape(x,nr*nc,nf),0);
        cost3 = eta*sum(diag(St));
    case 'l+s'
        term3 = T*x;
        cost3 = eta*sum(abs(term3(:)));   
        
        [Ut,St,Vt]=svd(reshape(x,nr*nc,nf),0);
        cost3s = etas*sum(diag(St));
        cost3 = cost3 + cost3s;
    case 'l+tv'
        term3 = applyD(x);
        cost3 = eta*sum(abs(term3(:)));
        
        [Ut,St,Vt]=svd(reshape(x,nr*nc,nf),0);
        cost3s = etas*sum(diag(St));
        cost3 = cost3 + cost3s;        
    case 'l+tv3d'
        term3 = applyD3d(x);
        cost3 = eta*sum(abs(term3(:)));
        
        [Ut,St,Vt]=svd(reshape(x,nr*nc,nf),0);
        cost3s = etas*sum(diag(St));
        cost3 = cost3 + cost3s;
    case 'l1+tv3d'
        term3 = applyD3d(x);
        cost3 = eta*sum(abs(term3(:)));
        
        term3s = T*x;
        cost3s = etas*sum(abs(term3s(:)));
        cost3 = cost3 + cost3s;       
end
else 
    cost3 = 0;
end
switch prior_motion
    case 'TV'
        [Du] = applyD(u);
        [Dv] = applyD(v);
        [Dux] = applyD(ux);
        [Dvx] = applyD(vx);
        [Duy] = applyD(uy);
        [Dvy] = applyD(vy);
        cost4 = gamma*(sum(abs(Du(:))) + sum(abs(Dv(:))) + sum(abs(Dux(:))) + sum(abs(Dvx(:)))+ sum(abs(Duy(:))) + sum(abs(Dvy(:))));
    case 'l2'
        cost4 = gamma*(sum(abs(u(:)).^2) + sum(abs(v(:)).^2) + sum(abs(ux(:)).^2) + sum(abs(uy(:)).^2) + sum(abs(vx(:)).^2) + sum(abs(vy(:)).^2));
    case 'l1'
        cost4 = gamma*(sum(abs(u(:))) + sum(abs(v(:))) + sum(abs(ux(:))) + sum(abs(uy(:))) + sum(abs(vx(:))) + sum(abs(vy(:))));
end

cost = cost1+cost2+cost3+cost4;
end


function [Du,D1x,D2x] = applyD(u)
[D1x,D2x] = gradient_op(u);
Du = cat(3,D1x,D2x);
end

function Dtu = applyDTrans(u,nf) % Input: u is of size nr*nc*2*nf
u1 = u(:,:,1:nf);
u2 = u(:,:,1+nf:end);
Dtu = div_op(u1, u2);
Dtu = -Dtu;
end

function [Du,D1x,D2x] = applyD3d(u)
[D1x,D2x,D3x] = gradient_op3d(u);
Du = cat(3,D1x,D2x,D3x);
end

function Dtu = applyD3dTrans(u,nf) % Input: u is of size nr*nc*2*nf
u1 = u(:,:,1:nf);
u2 = u(:,:,1+nf:2*nf);
u3 = u(:,:,1+2*nf:end);
Dtu = div_op3d(u1, u2,u3);
Dtu = -Dtu;
end

function [out,St] = proxNuclearnorm(in,nr,nc,nf,T)
[Ut,St,Vt]=svd(reshape(in,nr*nc,nf),0);
St=diag(prox1Norm(diag(St),St(1)*T));
x=Ut*St*Vt';
out = reshape(x,nr,nc,nf);
end

function out = prox1Norm(in,T)
out =  max(abs(in)-T,0)./(max(abs(in)-T,0)+T).*in;
end

function H = two_scale_filters_computation(J,N)

l = -(N+1)/2:(N+1)/2;
[u2N, c] = u2N_FIR_coefs(N);
h = c*u2N; % h00

h = fliplr(h);
l = fliplr(l);

H = cell(1,J+1);

H{1}.h_00 = h;
H{1}.h_10 = l.*h;
H{1}.h_11 = h;
H{1}.h_20 = 1*l.^2.*h;
H{1}.h_21 = 2*l.*h;
H{1}.h_22 = h;
for kk = 2:J+1
    H{kk}.h_00 = H{1}.h_00;
    H{kk}.h_10 = 2^(kk-1)*H{1}.h_10;
    H{kk}.h_11 = H{1}.h_11;
    H{kk}.h_20 = 2^(2*(kk-1))*H{1}.h_20;
    H{kk}.h_21 = 2^(kk-1)*H{1}.h_21;
    H{kk}.h_22 = H{1}.h_22;
end
end

function [Wx,Mx] = ApplyW(x,nscale,H,w0,w1,outCase)
[Wx,Mx] = bspline_mom(x,H,w0,w1,nscale,outCase);
end

function [output,M] = bspline_mom(I,H,wfilt0,wfilt1,nscale,outCase)
M = cell(nscale,1);

% y - filtering
m_0_x0 = convolution2D(wfilt0,1,I);
m_0_x1 = convolution2D(wfilt1,1,I);
clear I
% x - filtering
m_0_00 = convolution2D(1,wfilt0,m_0_x0);
m_0_10 = convolution2D(1,wfilt1,m_0_x0);
m_0_01 = convolution2D(1,wfilt0,m_0_x1);   
clear m_0_x0 m_0_x1 wfilt0 wfilt1

M{1} = struct('m_00',m_0_00,'m_01',m_0_01,'m_10',m_0_10);
M{1}.msize = size(m_0_00);
clear m_0_00 m_0_01 m_0_10 

J = numel(H)-1;
if nscale<=J
    for jj = 1:nscale-1
        M{jj+1} = moments_pyramid(M{jj},H{jj});
    end
else
    error('Invalid value of nscale');
end
if outCase==0
    output = M{nscale};
else
    output = M{nscale}.m_00;
end
end

function [M_coarser] = moments_pyramid(M_finer,H_finer)
m_00_x = filter_and_downsample_x(M_finer.m_00,H_finer.h_00);
m_01_x = filter_and_downsample_x(M_finer.m_01,H_finer.h_00);
M_coarser.m_00 = filter_and_downsample_y(m_00_x,H_finer.h_00);
M_coarser.m_01 = filter_and_downsample_y(m_00_x,H_finer.h_10) + filter_and_downsample_y(m_01_x,H_finer.h_11);
clear m_00_x m_01_x    
m_00_y = filter_and_downsample_y(M_finer.m_00,H_finer.h_00);
m_10_y = filter_and_downsample_y(M_finer.m_10,H_finer.h_00);

M_coarser.m_10 = filter_and_downsample_x(m_00_y,H_finer.h_10) + filter_and_downsample_x(m_10_y,H_finer.h_11);
clear m_00_y m_10_y
M_coarser.msize = size(M_coarser.m_00);
end

function [WTx,MTx] = ApplyWTrans(x,nscale,H,w0,MW)                                                                                                                                                                                                            
[WTx,MTx] = bspline_mom_trans(x,H,w0,nscale,MW);
end

function [output,M] = bspline_mom_trans(I,H,wfilt0,nscale,MW)
M = cell(nscale,1);
M{nscale} = struct('m_00',I);
if nscale ~=1
    for jj=nscale-1:-1:1
    M{jj} = filter_and_upsample_xy(M{jj+1}.m_00,H{jj+1}.h_00,MW{jj}.msize);
    end
end
out = M{1}.m_00;
out = convolution2D(1,rot90(wfilt0,2),out);
output = convolution2D(rot90(wfilt0,2),1,out);
end

function Ifilt = filter_and_downsample_x(I,v)
Ifilt = convolution2D(1,v,I);
Ifilt = Ifilt(:,1:2:end,:);
end

function Ifilt = filter_and_downsample_y(I,v)
Ifilt = convolution2D(v,1,I);
Ifilt = Ifilt(1:2:end,:,:);
end

function output = convolution2D(h1,h2,signal)
nt = size(signal,3);
output = zeros(size(signal));
for i = 1:nt
    output(:,:,i) = conv2(h1,h2,signal(:,:,i),'same');
end
end

function out = filter_and_upsample_xy(I,h,msize)
Ifilt = filter_and_upsample_y(I,h,msize);
Ifilt = filter_and_upsample_x(Ifilt,h,msize);
out.m_00 = Ifilt;
end

function Ifilt = filter_and_upsample_x(I,v,msize) 
[nr,nc,nt] = size(I);
Ifilt = zeros(nr,msize(2),nt);
Ifilt(:,1:2:end,:) = I;
Ifilt = convolution2D(1,v,Ifilt);
end

function Ifilt = filter_and_upsample_y(I,v,msize)
[nr,nc,nt] = size(I);
Ifilt = zeros(msize(1),nc,nt);
Ifilt(1:2:end,:,:) = I;
Ifilt = convolution2D(v,1,Ifilt);
end


