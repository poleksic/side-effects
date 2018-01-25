% This is Algorithm 2 from Witten's paper 
% (Witten et al. Biostatistics 10, no. 3 (2009): 515-534.)
% X is nxp matrix 
% v is px1 vector
% u is nx1 vector
% reasonable is to use maxiter ~ 200 and precision ~ 1e-7

%     INPUT:
% C1 (reg on u) should range from 1 and sqrt(n)
% C1 (reg on v) should range from 1 and sqrt(p)

%     OUTPUT:
% columns of u are u_1, ..., u_k (each one is n-dimensional)
% columns of v are v_1, ..., v_k (each one is p-dimensional)
% columns of d are d_1, ..., d_k (each one is 1-dimensional)
function  [u, v, d, convu, convv, lastPossK] = K_PMD(K, X, C1, C2, initv, maxiter, precision)
    [n, p] = size(X);
    copyX = X;
    u = zeros(n,K);
    v = zeros(p,K);
    d = zeros(1,K);
    convu = zeros(1,K);
    convv = zeros(1,K);
    lastPossK = K;
    for k=1:K
       [ u(:,k), v(:,k), d(:,k) convu(:,k) convv(:,k) ] = PMD(copyX, C1, C2, initv, k, maxiter, precision);
       if convu(:,k) * convv(:,k) < 1
           warning('last cannonical component is %d', k-1)
           lastPossK = k-1;
           break
       end
       copyX = copyX - d(:,k) * u(:,k) * v(:,k)';
    end
end