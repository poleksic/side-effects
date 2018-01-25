function [ AVG_AUC_REF AVG_AUPR_REF AVG_AUC_RND AVG_AUPR_RND AVG_AUC_COS AVG_AUPR_COS AVG_AUC_ML AVG_AUPR_ML AVG_AUC_CCA AVG_AUPR_CCA ] = ParallelCrossVal(LOWER_DRUGS_PER_SE, UPPER_DRUGS_PER_SE, mode, workers, folds, round, ADR_FOLDER, outfileID)
    
    
    % we will take advantage of parallel computing environment
    myCluster = parcluster('local');
    myCluster.NumWorkers = workers;
    parpool(myCluster,workers);
    
    % specify parameters: X
    % K - number of canonical components for CCA;
    % J - number of nearest neighbors in ML
    % rnk, iter and lambdas for CS
    K = 60; J = 5; rnk = 100; iter = 100; lR = 1.0; lM = 0.1; lN = 10;

    % specify the input matrices
    RmQ = dlmread(strcat(ADR_FOLDER,'INTERACTION_MATRIX'),' ',1,1); % SIDER associatiowithout postmarketing
    Q = dlmread(strcat(ADR_FOLDER,'IMPUTE_MATRIX'),' ',1,1); % Q is only postmarketing
    M = dlmread(strcat(ADR_FOLDER,'TANIMOTO_MATRIX'),' ',0,0); % Tanimoto similarity
    N = dlmread(strcat(ADR_FOLDER,'ADR_PATH_LESK_MATRIX'),' ',1,1); % UMLS similarity of ADRs
    CFEAT = dlmread(strcat(ADR_FOLDER,'FEATURE_VECTORS'),' ',1,1); % 881 PubChem feature
    R = RmQ + Q;
    
    % set degree matrices to be used in the CS method
    [DM nM]= GetDiag(M,J);
    [DN, nN]= GetDiag(N,J);
    DMM = DM-nM;
    DNN = DN-nN;
    
    [m,n] = size(R);


    if mode == 2
        fprintf('mode = 2\n');
        % select collections of columns (ADRs) to run CV segregated by ADRs
        indices = ColdStartTestIndicesSet( R, folds, 'CSCOLS' ); 
	all_indices = indices;
    else 
        SCOL = sum(R,1);
        INCLUDE_COLS = find(SCOL >= LOWER_DRUGS_PER_SE & SCOL <= UPPER_DRUGS_PER_SE);
        TMP = zeros(m,n);
        TMP(:,INCLUDE_COLS) = 1;

        % these are indices of all elements where the adr has specified number
        % of drugs
        INCLUDE_COLS_IND = find(TMP > 0.5);

        if mode == 1
            fprintf('mode = 1\n');
            % select collections of rows (drugs) to run CV segregated by ADRs
            all_indices = ColdStartTestIndicesSet( R, folds, 'CSROWS' );
            for t = 1:folds
                indices{t} = intersect(all_indices{t}, INCLUDE_COLS_IND);
            end 
        else
            fprintf('mode = 0\n');
            all_indices = OffTargetTestIndicesSet( R, folds );
            for t = 1:folds
                indices{t} = intersect(all_indices{t}, INCLUDE_COLS_IND);
            end 
        end
    end
        
    AUC_REF = zeros(1, folds);
    AUPR_REF = zeros(1, folds);    
    AUC_RND = zeros(1, folds);
    AUPR_RND = zeros(1, folds);
    AUC_COS = zeros(1, folds);
    AUPR_COS = zeros(1, folds);
    AUC_ML = zeros(1, folds);
    AUPR_ML = zeros(1, folds);
    AUC_CCA = zeros(1, folds);
    AUPR_CCA = zeros(1, folds);
    
    AVG_AUC_REF = 0; AVG_AUPR_REF = 0;
    AVG_AUC_RND = 0; AVG_AUPR_RND = 0;
    AVG_AUC_COS = 0; AVG_AUPR_COS = 0;
    AVG_AUC_ML = 0; AVG_AUPR_ML = 0;
    AVG_AUC_CCA = 0; AVG_AUPR_CCA = 0;

    parfor t = 1:folds
        if(length(indices{t}) <= 0) 
            continue
        end
        
        % these are the true values (that are beingmasked out)
        GOLDMTX = R(indices{t});
        GOLD = GOLDMTX(:);
        
        TEST_MTX = R;
        % hide (mask out) the set of rows
        TEST_MTX(all_indices{t}) = 0; 

        % REF method (see the paper)
        FREQ = sum(TEST_MTX,1)/m;
        EXC_REF = TEST_MTX;        
        [ii, jj] = ind2sub(size(TEST_MTX), indices{t});
        for x = 1:length(ii)
            EXC_REF(ii(x),jj(x)) = FREQ(jj(x));
        end
        
        % RANDOM CLASSIFIER
        EXC_RND = TEST_MTX;
        EXC_RND(indices{t}) = rand(1, length(indices{t}));
        
        % our own Compressed Sensing method (CS)
        MULT_LAB = MLKNN_TEST(R, M, indices{t}, J, 1);
        raw_normalization = 0;
        [ unimportant, ExcludedColumns ] = find(sum(TEST_MTX,1)==0); 
        [ ExcludedRows unimportant ] = find(sum(TEST_MTX,2)==0); 
        W = max(1, 6 * TEST_MTX);
        IMPUTE = zeros(m,n);
        [F G] = WeightImputeLogFactorization(TEST_MTX,DMM,DNN,W,IMPUTE,lR,lM,lN,iter,rnk);
        [F G] = WeightedProfile(F, G, M, N, ExcludedRows, ExcludedColumns, J, raw_normalization);
        EXC_COS = GetP(F*G');
        if mode == 1
            EXC_COS(ExcludedRows,:) = EXC_COS(ExcludedRows,:) .* MULT_LAB(ExcludedRows,:);
        end
        
        % multi-label learning method (ML)
        EXC_ML = MULT_LAB;
        
        % Canonical Correlation Analysis (CCA)
        [ u, v, d, convu, convv, lastPossK ] = ComputeCCA(CFEAT, TEST_MTX, K);
        SUB_SE_MTX = GenerateCCAMatrix( u, v, d, 1, K, lastPossK );
        EXC_CCA = CFEAT * SUB_SE_MTX'; 
        
        PREDMTX_REF = EXC_REF(indices{t});
        PRED_REF = PREDMTX_REF(:);
        PREDMTX_RND = EXC_RND(indices{t});
        PRED_RND = PREDMTX_RND(:);
        PREDMTX_COS = EXC_COS(indices{t});
        PRED_COS = PREDMTX_COS(:);
        PREDMTX_ML = EXC_ML(indices{t});
        PRED_ML = PREDMTX_ML(:);
        PREDMTX_CCA = EXC_CCA(indices{t});
        PRED_CCA = PREDMTX_CCA(:);
        
        [ ~, ~, ~, AUC_REF(t) ] = perfcurve(GOLD,PRED_REF,1);
        [ ~, ~, ~, AUC_RND(t) ] = perfcurve(GOLD,PRED_RND,1);
        [ ~, ~, ~, AUC_COS(t) ] = perfcurve(GOLD,PRED_COS,1);        
        [ ~, ~, ~, AUC_ML(t) ] = perfcurve(GOLD,PRED_ML,1);
        [ ~, ~, ~, AUC_CCA(t) ] = perfcurve(GOLD,PRED_CCA,1);

        
        [ ~, ~, ~, AUPR_REF(t) ] = perfcurve(GOLD, PRED_REF, 1, 'xCrit', 'reca', 'yCrit', 'prec');
        [ ~, ~, ~, AUPR_RND(t) ] = perfcurve(GOLD, PRED_RND, 1, 'xCrit', 'reca', 'yCrit', 'prec');
        [ ~, ~, ~, AUPR_COS(t) ] = perfcurve(GOLD, PRED_COS, 1, 'xCrit', 'reca', 'yCrit', 'prec');
        [ ~, ~, ~, AUPR_ML(t) ] = perfcurve(GOLD, PRED_ML, 1, 'xCrit', 'reca', 'yCrit', 'prec');
        [ ~, ~, ~, AUPR_CCA(t) ] = perfcurve(GOLD, PRED_CCA, 1, 'xCrit', 'reca', 'yCrit', 'prec');

        fprintf('AUC -- cos: %f  ml: %f  cca: %f freq: %f  rnd: %f\n', AUC_COS(t), AUC_ML(t), AUC_CCA(t), AUC_REF(t), AUC_RND(t));
        fprintf('AUPR - cos: %f  ml: %f  cca: %f freq: %f  rnd: %f\n', AUPR_COS(t), AUPR_ML(t), AUPR_CCA(t), AUPR_REF(t), AUPR_RND(t));

    end
    
    for t = 1:folds
        AVG_AUC_REF = AVG_AUC_REF + AUC_REF(t);
        AVG_AUPR_REF = AVG_AUPR_REF + AUPR_REF(t);
        AVG_AUC_RND = AVG_AUC_RND + AUC_RND(t);
        AVG_AUPR_RND = AVG_AUPR_RND + AUPR_RND(t);
        AVG_AUC_COS = AVG_AUC_COS + AUC_COS(t);
        AVG_AUPR_COS = AVG_AUPR_COS + AUPR_COS(t);
        AVG_AUC_ML = AVG_AUC_ML + AUC_ML(t);
        AVG_AUPR_ML = AVG_AUPR_ML + AUPR_ML(t);
        AVG_AUC_CCA = AVG_AUC_CCA + AUC_CCA(t);
        AVG_AUPR_CCA = AVG_AUPR_CCA + AUPR_CCA(t);
    end
    
    AVG_AUC_REF = AVG_AUC_REF / folds;
    AVG_AUPR_REF = AVG_AUPR_REF / folds;    
    AVG_AUC_RND = AVG_AUC_RND / folds;
    AVG_AUPR_RND = AVG_AUPR_RND / folds;
    AVG_AUC_COS = AVG_AUC_COS / folds;
    AVG_AUPR_COS = AVG_AUPR_COS / folds;
    AVG_AUC_ML = AVG_AUC_ML / folds;
    AVG_AUPR_ML = AVG_AUPR_ML / folds;
    AVG_AUC_CCA = AVG_AUC_CCA / folds;
    AVG_AUPR_CCA = AVG_AUPR_CCA / folds;

    fprintf(outfileID,'round: %d\n', round);
    fprintf(outfileID,'%f  %f  %f  %f  %f\n', AVG_AUC_REF, AVG_AUC_RND, AVG_AUC_COS, AVG_AUC_ML, AVG_AUC_CCA);
    fprintf(outfileID,'%f  %f  %f  %f  %f\n\n', AVG_AUPR_REF, AVG_AUPR_RND, AVG_AUPR_COS, AVG_AUPR_ML, AVG_AUPR_CCA);
    
    % shut down parallel cluster
    poolobj = gcp('nocreate');
    delete(poolobj);
end
