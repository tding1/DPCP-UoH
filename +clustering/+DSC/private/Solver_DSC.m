% -------------------------------------------------------------------------
% This function solves:
% ||A^T*X||_1,p + gamma*||Z||_1
% Subject to  A=X*Z    AND    diag(A^T * Q) = 1

% X: The given data matrix
% mu: The parameter mu in the ADMM algorithm
% gamma: The parameter gamma in the algorithm
% itr: The maximum number of iterations allowed for the ADMM algorithm
% nrm: It determines if we use ell_1 norm or ell_2 norm.
% -------------------------------------------------------------------------
% Copyright @ Mostafa Rahmani, 2017
% -------------------------------------------------------------------------

function t = Solver_DSC(X,Q , mu , gamma, itr , nrm)

[N1,N2] = size(X) ;
[~,N2q] = size(Q) ;

a = normailize_column(Q) ;

mi = mu^-1 ;

y1  = zeros(N1,N2q) ;
y2 = zeros(N2q,1) ;
y3 = zeros(N2,N2q) ;
y4 = zeros(N2,N2q) ;
u = randn(N2,N2q)/3 ;



G = mi*((eye(N1) + Q*Q' + X*X')^(-1)) ;
G2 = mi*((X'*X + eye(N2))^(-1)) ;




for k = 1:itr
    Aux_1 = X'*a ;
    t = Aux_1 - mi*y3 ;
    
    if (nrm == 1)
        t = sign(t).*( max(abs(t) - mi , 0) ) ;
        
    else
        t = shrinkage_12(t,mi) ;
        
    end
    
    h = u - mi*y4 ;
    z = sign(h).*( max(abs(h) - mi*gamma , 0) ) ;
    
    Aux_2 = X*u ;
    a = G*(mu*Aux_2 + mu*Q + mu*X*t - y1 - Q*diag((y2)) + X*y3) ;
    
    u = G2*(mu*Aux_1 + mu*z + X'*y1 + y4) ;
    
    y1 = y1 + mu*(a - Aux_2) ;
    y2 = y2 + mu*(diag(a'*Q) - ones(N2q,1)) ;
    y3 = y3 + mu*(t - Aux_1) ;
    y4 = y4 +  mu*(z - u);
    
    %norm(vec(t - Aux_1))/norm(vec(t))
    if( norm(t - Aux_1) <= 0.001)
        break
    end
    
    solver_iteration_number = k ;
    if(rem(solver_iteration_number,10) == 0)
        solver_iteration_number
    end
end

