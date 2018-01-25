%    INPUT:
% Y - mx881 matrix of 881 fingerprints for each one of m drugs
% R - the incomplete mxn matrix of drug-ADR associations (to be
% reconstracted)
% K - the number of canonical components 
% (see Pauwels E. et al. (2011) BMC Bioinformatics, 12:169)

%    OUTPUT:


function [ u, v, d, convu, convv, lastPossK ] = ComputeCCA( Y, R, K )

    X = full(R')*full(Y); % n=1385 se vs p=881 sub
    
    [n,p] = size(X);
    C1 = sqrt(n)/2;
    C2 = sqrt(p)/2;
    maxiter = 400;
    precision = 1e-7;

%     Mn = mean(X);
%     Sd = std(X);
%     Sd(Sd==0) = 1;
% 
%     Xn = bsxfun(@minus,X,Mn);
%     X = bsxfun(@rdivide,Xn,Sd);
    
    % initialize series of K v's
    [U d V]=svd(X);
    initv = V(:,1:K);
    % end initialize
    
    % note u is nxK, v is pxk, d is 1xk
    [u, v, d, convu, convv, lastPossK] = K_PMD(K, X, C1, C2, initv, maxiter, precision);
end

