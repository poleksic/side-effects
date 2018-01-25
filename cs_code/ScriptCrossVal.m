clear;
workers = 10; % specify number of available workers (typically # of cores) for MATLAB Paraller Computing Toolbox
folds = 10; % number of CV folds
tot_rounds = 2; % number of CV rounds
mode = 0;

ADR_FOLDER = '/home/aleksandar/cs_code_data/cs_data/';

OUT_FILE_1_12 = 'out_file_1_12';
fileID_1_12 = fopen(strcat(ADR_FOLDER,OUT_FILE_1_12),'w');
OUT_FILE_1_25 = 'out_file_1_25';
fileID_1_25 = fopen(strcat(ADR_FOLDER,OUT_FILE_1_25),'w');
OUT_FILE_1_50 = 'out_file_1_50';
fileID_1_50 = fopen(strcat(ADR_FOLDER,OUT_FILE_1_50),'w');
OUT_FILE_1_100 = 'out_file_1_100';
fileID_1_100 = fopen(strcat(ADR_FOLDER,OUT_FILE_1_100),'w');
OUT_FILE_1_200 = 'out_file_1_200';
fileID_1_200 = fopen(strcat(ADR_FOLDER,OUT_FILE_1_200),'w');
OUT_FILE_1_400 = 'out_file_1_400';
fileID_1_400 = fopen(strcat(ADR_FOLDER,OUT_FILE_1_400),'w');
OUT_FILE_1_800 = 'out_file_1_800';
fileID_1_800 = fopen(strcat(ADR_FOLDER,OUT_FILE_1_800),'w');
OUT_FILE_1_1500 = 'out_file_1_1500';
fileID_1_1500 = fopen(strcat(ADR_FOLDER,OUT_FILE_1_1500),'w');


fprintf(fileID_1_12,'MAX  RND  COS  ML  CCA\n');
fprintf(fileID_1_25,'MAX  RND  COS  ML  CCA\n');
fprintf(fileID_1_50,'MAX  RND  COS  ML  CCA\n');
fprintf(fileID_1_100,'MAX  RND  COS  ML  CCA\n');
fprintf(fileID_1_200,'MAX  RND  COS  ML  CCA\n');
fprintf(fileID_1_400,'MAX  RND  COS  ML  CCA\n');
fprintf(fileID_1_800,'MAX  RND  COS  ML  CCA\n');
fprintf(fileID_1_1500,'MAX  RND  COS  ML  CCA\n');

for round=1:tot_rounds
    ParallelCrossVal(1, 12, mode, workers, folds, round, ADR_FOLDER, fileID_1_12);
    ParallelCrossVal(1, 25, mode, workers, folds, round, ADR_FOLDER, fileID_1_25);
    ParallelCrossVal(1, 50, mode, workers, folds, round, ADR_FOLDER, fileID_1_50);
    ParallelCrossVal(1, 100, mode, workers, folds, round, ADR_FOLDER, fileID_1_100);
    ParallelCrossVal(1, 200, mode, workers, folds, round, ADR_FOLDER, fileID_1_200);
    ParallelCrossVal(1, 400, mode, workers, folds, round, ADR_FOLDER, fileID_1_400);
    ParallelCrossVal(1, 800, mode, workers, folds, round, ADR_FOLDER, fileID_1_800);
    ParallelCrossVal(1, 1500, mode, workers, folds, round, ADR_FOLDER, fileID_1_1500);
end

fclose(fileID_1_12); fclose(fileID_1_25); fclose(fileID_1_50); fclose(fileID_1_100);
fclose(fileID_1_200); fclose(fileID_1_400); fclose(fileID_1_800); fclose(fileID_1_1500);
