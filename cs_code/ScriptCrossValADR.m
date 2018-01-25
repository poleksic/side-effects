clear;
workers = 10; % specify number of available workers (typically # of cores) for MATLAB Paraller Computing Toolbox
folds = 10; % number of CV folds
tot_rounds = 2; % number of CV rounds
mode = 2;

ADR_FOLDER = '/home/aleksandar/cs_code_data/cs_data/'; % default folder for output; also the folder where data matrices are stored
OUT_FILE = 'out_file_ADR';
outfileID = fopen(strcat(ADR_FOLDER,OUT_FILE),'w');

fprintf(outfileID,'AUC_MAX  AUC_RND  AUC_COS  AUC_ML  AUC_CCA\n');

% run tot_rounds of CV experiments
for round=1:tot_rounds
    ParallelCrossVal(0,0,mode, workers, folds, round, ADR_FOLDER, outfileID);
end

fclose(outfileID);
