function [x_k,costs] = MC_function(A,At,x,b,opts)
x_km1 = x;
if opts.periodic 
    xp = circshift(x,1,3); % periodic
    xpb = circshift(x,-1,3);
else
    xp = x;
    xp(:,:,2:end) = x(:,:,1:end-1);
end 
z_km1 = zeros(size(b));
t = opts.stepsize;
delta = 0.99;
alpha = 1e-1; 
reductionFactor = 0.5;
maxIter = opts.maxIter_mc;
tol = opts.tol_mc;
showTrigger = opts.showTrigger;
errFcn = opts.errFcn;
tkm1 = t;
ma = opts.ma;
theta = 1;
costs = [];
opts.Fx_case = 2;
% forwardmotion = @(z) reshape(opts.U(z(:)),size(x));

for k = 0:maxIter   
    if k== 0
        errXp = errFcn(x_km1);
    else
        errXp = errX;
    end
    Ktz1 = At(z_km1);
    in_1 = x_km1 - tkm1*Ktz1;
    x0 = forwardmotion(xp,opts);
    x_k = x0 + prox1Norm(in_1 - x0,ma);
    
    if opts.periodic 
        xp = circshift(x_k,1,3); % periodic
    else
        xp = x_k;
        xp(:,:,2:end) = x_k(:,:,1:end-1);
    end 
    t = tkm1*sqrt(1+theta);
    accept_t = 0;
    cost = Cost_MC(A,b,x_k,x0,ma);
    costs = [costs,cost];
    errX = errFcn(x_k);
    if mod(k,showTrigger) == 0
        fprintf('Motion compensation starts... ');
        fprintf('Iter: %3d  ~ Cost: %.3e ~ t: %.2e ~ ~ Error: %.2e\n',k,cost,t,errX); 
    end

    if k>2
        if abs(costs(k)-costs(k-1))<tol
            costs = costs(1:k);
            break
        end  
%         if errXp<=errX
%             break
%         end
    end
        
    while accept_t ==0
        theta = t/tkm1;
        s = alpha*t;    
        term = x_k+theta*(x_k-x_km1);
        in_2 = z_km1 + s*A(term);
        z_k = (in_2 - s*b)/(1+s);
        Ktz1_new = At(z_k);
        Ktz_diff = Ktz1_new - Ktz1;
        z_diff = z_k - z_km1;
        lhs = t*s*(norm(Ktz_diff(:))^2);
        rhs = delta^2*(norm(z_diff(:))^2);    

        if lhs <= rhs
            tkm1 = t;
            accept_t = 1;
        else
            t = reductionFactor*t;
            tkm1 = t;
        end 
    end
    x_km1 = x_k;
    z_km1 = z_k;
end
end

function cost = Cost_MC(A,b,x,xp,alpha)
term1 = A(x)-b;
cost1 = .5*sum(abs(term1(:)).^2);
cost2 = alpha*sum(abs(x(:)-xp(:)));
cost = cost1+cost2;
end


function x1 = forwardmotion(x0,opts)
u = opts.uJ;
v = opts.vJ;

[nr,nc,nf] = size(x0);
[xI, yI] = meshgrid(1:nc,1:nr);
x1 = x0;
for frame = 1:nf
    locs = xI+u(:,:,frame);
    locs(locs<1) = 1; 
    locs(locs>nc) = nc;
    hor_ind = locs; 

    locs = yI+v(:,:,frame);
    locs(locs<1) = 1; 
    locs(locs>nr) = nr;
    ver_ind = locs; 

    Ip_m = x0(:,:,frame); % REFERENCE FRAME
    F_I0 = linear_mc(Ip_m,hor_ind, ver_ind, nr, nc);
    x1(:,:,frame) = F_I0;
      
end
end

function x1 = backwardmotion(x0,opts)
u = opts.ubJ;
v = opts.vbJ;

[nr,nc,nf] = size(x0);
[xI, yI] = meshgrid(1:nc,1:nr);
x1 = x0;
for frame = 1:nf
    locs = xI+u(:,:,frame);
    locs(locs<1) = 1; 
    locs(locs>nc) = nc;
    hor_ind = locs; 

    locs = yI+v(:,:,frame);
    locs(locs<1) = 1; 
    locs(locs>nr) = nr;
    ver_ind = locs; 

    Ip_m = x0(:,:,frame); % REFERENCE FRAME
    B_I0 = linear_mc(Ip_m,hor_ind, ver_ind, nr, nc);
    x1(:,:,frame) = B_I0;
      
end
end
