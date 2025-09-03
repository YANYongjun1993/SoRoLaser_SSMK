%  Function: [x,fval,it,X] = seq_quad_prog (f, gradf, hessf,
%                                    A, b, G, r, x0, itmax, tol)
%
% min.  1/2 z'Hz + f'z
% s.t.  Gz = h
%       Az <= b
%
%  using the sequential quadratic programming method.
%
%  f     : The name of the objective function f : R^n -> R.
%  gradf : The name of the function that returns the gradient of f.
%  hessf : The name of the function that returns the hessian matrix of f.
%  A     : A matrix with dimension m x n for the linear equality contraints.
%  b     : A vector with m elements for the linear equality contraints.
%  G     : A matrix with dimension p x n for the linear inequality constraints.
%  r     : A vector with p elements for the linear inequality constraints.
%  x0    : The start point for the algorithm.
%  itmax : The maximal number of iterations allowed.
%  tol   : The bound needed for the stop criteria.

function [x,fval,flag] = SeqQuadProg(xhat,km,opt,size_x,size_u,p,N,der)
%% Initialization
coder.varsize('DuDp','DxDc','DuDg','DxDJ','DuDJ','G','H','DuuDL','DxxDL','DxuDL','He','f','A','b','G','h','eqlin','ineqlin');
stop = false;
it = 1;
flag = 1;
pos = p;
size_Cbar = 4;
size_f = 2;

x0 = km.z;
eqlin = km.l;
ineqlin = km.v;
    while( ~stop )
        xhat = singularity_check(xhat);
        uu = reshape(km.z(1:size_u*N,1),[size_u,N]);
        xx = reshape(km.z(size_u*N+1:size_u*N+size_x*N,1),[size_x,N]);
        ss = km.z(size_u*N+size_x*N + 1:size_u*N+size_x*N + N, 1);
        ref = [pos(2,:);pos(3,:)];
       %% compute gradients of the constraints 
        DuDp = [];
        DxDc = [];
        for j = 1:N
            der.Cu(:,:,j) = DCDu(xx(:,j)',uu(:,j)');
            DuDp = blkdiag(DuDp,der.Cu(:,:,j));
            der.Cbarx(:,:,j) = DCbarDx(xx(:,j)');
            DxDc = blkdiag(DxDc,der.Cbarx(:,:,j));
        end
        Gamma = kron(speye(N), ones(size_Cbar,1));
        DzDh = blkdiag(DuDp, DxDc);
        DzDh = [DzDh,[zeros(size(DuDp,1),N);-Gamma];zeros(N,size(DzDh,2)),-eye(N)];
        % compute gradients of the dynamic equation, referring to the page 38 of modele6 slides
        DuDg = [];
        for j = 1:N
            if j == 1
                der.fu(:,:,j) = DfDu(xhat',uu(:,j)');
                DuDg = blkdiag(DuDg,-der.fu(:,:,j));
            else
                der.fu(:,:,j) = DfDu(xx(:,j-1)',uu(:,j)');
                DuDg = blkdiag(DuDg,-der.fu(:,:,j));
            end
        end
        for j = 1:N-1
            der.fx(:,:,j) = DfDx(xx(:,j)',uu(:,j+1)');                
        end
        DxDg = zeros(size_x*N,size_x*N);
        for j = 1:N-1
            E = zeros(N,N);
            E(j+1,j) = -1;
            DxDg = DxDg + kron(E,der.fx(:,:,j));
        end
        E = eye(N,N);
        DxDg = DxDg + kron(E,eye(size_x));
        DzDg = [DuDg,DxDg,zeros(N*size_f,N)];
        % compute gradients of the cost function, referring to the page 37 of modele6 slides
        DxDJ = [];
        DuDJ = [];
        for j = 1:N
            if j  < N
                der.Lx(:,:,j) = DLDx(xx(:,j)',uu(:,j+1)',ref(:,j)');
                DxDJ = [DxDJ,der.Lx(:,:,j)];
            else
                der.Phix(:,:,j) = DPhiDx(xx(:,j)',uu(:,j)',ref(:,j)');
                DxDJ = [DxDJ,der.Phix(:,:,j)];
            end
            if j == 1
                der.Lu(:,:,j) = DLDu(xhat',uu(:,j)',ref(:,j)');
                DuDJ = [DuDJ,der.Lu(:,:,j)];
            else
                der.Lu(:,:,j) = DLDu(xx(:,j-1)',uu(:,j)',ref(:,j)');
                DuDJ = [DuDJ,der.Lu(:,:,j)];
            end
        end
        DzDJ = [DuDJ,DxDJ,opt.gamma*ones(1,N)];
        % compute gradient of the Lagrangian, referring to the page 41 of modele6 slides
        DzDL = DzDJ' + DzDh'*ineqlin + DzDg'*eqlin;
        % natural residual
        G = [ ];
        G = [xx(:,1)- valf(xhat', uu(:,1)')];
        % G = [xx(:,1)- valf(xhat', u_0')];
        for j=1:N-1                   % generate dynamic constraints
             G = [G;xx(:,j+1)- valf(xx(:,j)', uu(:,j+1)')];
        end
        H = [];
        for j = 1:N
            H = [H; valC(xx(:,j)',uu(:,j)')];
        end
        for j = 1:N
            H = [H; valCbar(xx(:,j)')];
        end
        H = [H; -ss];
        F_NR = [DzDL; G; min(-H,ineqlin)];  % |-H|?
        % compute Hessian of the cost function, referring to the page 37 of modele6 slides
        DuuDL = [];
        DxxDL = [];
        DxuDL = [];
        for j = 1:N
            if j == 1
                % der.Luu(:,:,j) = D2LDuu(xhat',u_0',ref(:,j)');
                der.Luu(:,:,j) = D2LDuu(xhat',uu(:,j)',ref(:,j)');
                DuuDL = blkdiag(DuuDL,der.Luu(:,:,j));
            else
                der.Luu(:,:,j) = D2LDuu(xx(:,j-1)',uu(:,j)',ref(:,j)');
                DuuDL = blkdiag(DuuDL,der.Luu(:,:,j));
            end
            if j == 1
                der.Lxu(:,:,j) = D2LDxu(xx(:,j)',uu(:,j)',ref(:,j)');
                DxuDL = blkdiag(DxuDL,der.Lxu(:,:,j));
            else
                der.Lxu(:,:,j) = D2LDxu(xx(:,j-1)',uu(:,j)',ref(:,j)');
                DxuDL = blkdiag(DxuDL,der.Lxu(:,:,j));
            end
            if j == 1
                der.Lxx(:,:,j) = D2LDxx(xhat',uu(:,j)',ref(:,j)');
                % der.Lxx(:,:,j) = D2LDxx(xhat',u_0',ref(:,j)');
                DxxDL = blkdiag(DxxDL,der.Lxx(:,:,j));
            elseif j < N
                der.Lxx(:,:,j) = D2LDxx(xhat',uu(:,j)',ref(:,j)');
                % der.Lxx(:,:,j) = D2LDxx(xx(:,j)',uu(:,j)',ref(:,j)');
                DxxDL = blkdiag(DxxDL,der.Lxx(:,:,j));
            else
                der.Phixx(:,:,j) = D2PhiDxx(xx(:,j)',uu(:,j)',ref(:,j)');
                DxxDL = blkdiag(DxxDL,der.Phixx(:,:,j));
            end
        end
        Hessian = [DuuDL,DxuDL,zeros(size(DuuDL,1),N); DxuDL',DxxDL,zeros(size(DxxDL,1),N); zeros(N,size(DuuDL,1)),zeros(N,size(DxxDL,1)),zeros(N,N)];
        He = real(Hessian); 
        eigenvalues = eig(He);  
        weight = eps;
        while all(eigenvalues(:) >= 0) == false
            He = He + eye(size(Hessian))*weight;
            weight = weight*10;
            eigenvalues = eig(He);
        end   
        f = DzDJ';        
        G_aeq = DzDg; 
        h = -G;  
        A = full(DzDh);
        b = -H;
        % Add the constraints for changing the output sampling rate of the
        % control variables
        % CrtSam = zeros((N-NumCrt)*size_u,size(G_aeq,2));
        % for raw = 1:(N-NumCrt)*size_u
        %     CrtSam(raw,raw+fix((raw-1)/16)*4) = 1;
        %     CrtSam(raw,raw+4+fix((raw-1)/16)*4) = -1;
        % end
        % G_aeq = [G_aeq;CrtSam];
        % h_CrtSam = zeros(size(CrtSam,1),1);
        % h = [h;h_CrtSam];
        options = optimoptions('quadprog','Algorithm','active-set','MaxIterations',3000,'Display','iter');
        % options = optimoptions(options,'MaxIterations',30);
        [x,fval,exitflag,output,lambda] = quadprog(He,f,A,b,G_aeq,h,[],[],x0,options);
        x0 = x;
        if( norm(F_NR) < opt.tol )
			stop = true; % => x is the solution
        else
            km.z = km.z + x;
            eqlin = lambda.eqlin;
            ineqlin = lambda.ineqlin;
			it = it + 1;
        end

        if (it >= opt.itmax)
            flag = 0;
			stop = true;
        end
    end  
    x = km.z(1:(size_u+size_x+1)*N,1);
    fval =0;
end