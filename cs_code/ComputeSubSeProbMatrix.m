function [ SUB_SE_PROB ] = ComputeSubSeProbMatrix(TEST_MTX, M, N, CFEAT, J, lR, lM, lN, iter, rnk, cutoff)
    % these will be count matrix: how many times a SUB (row) is associated with a SE (column) 
    
    SUB_SE = CFEAT' * TEST_MTX;
    
    A = logical(SUB_SE);
    sum_first = sum(A(:));
    sum_first
    
    REF_SUB_SE = RefineSubSe(SUB_SE, CFEAT, cutoff);
    
    sum_second = sum(REF_SUB_SE(:));
    sum_second
    
    EXCL_COLS = sum(REF_SUB_SE,1) <= 0;
    EXCL_ROWS = sum(REF_SUB_SE,2) <= 0;
    
    INCL_COLS = sum(REF_SUB_SE,1) > 0;
    INCL_ROWS = sum(REF_SUB_SE,2) > 0;
    
    [m, n] = size(REF_SUB_SE);

    RND_M = 0.5 + 0.001 * rand(m,m);
    RND_N = 0.5 + 0.00001 * rand(n,n);

    M = max(eye(m),RND_M);
%     N = max(eye(n),RND_N);

    [DM nM]= GetDiag(M,J);
    [DN nN]= GetDiag(N,J);
    DMM = DM-nM;
    DNN = DN-nN;

    W = max(1, 6 * REF_SUB_SE); Q = zeros(m,n);

    [F G] = WeightImputeLogFactorization(REF_SUB_SE,DMM,DNN,W,Q,lR,lM,lN,iter,rnk);
    SUB_SE_PROB = GetP(F*G');

end
