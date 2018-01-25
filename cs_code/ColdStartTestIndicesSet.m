function [ indices ] = ColdStartTestIndicesSet( MTX, nfolds, MODE )
    rowsorcols = cell(1, nfolds);
    indices = cell(1, nfolds);

    % MODE = 'CSROWS' means cold start on rows (e.g. new drugs)
    % MODE = 'CSCOLS' means cold start on cols (e.g. new proteins or new side effects)
    [nrows, ncols] = size(MTX)
    
    if strcmp(MODE, 'CSROWS')
        maindim = nrows;
        otherdim = ncols;
    elseif strcmp(MODE, 'CSCOLS')
        maindim = ncols;
        otherdim = nrows;
    else
        msg = 'ColdStartTestIndicesSet: Incorrect mode specified!';
        error(msg)
        return;
    end 
    
    if mod(maindim,nfolds) == 0
        binsize = maindim/nfolds - 1;  
    else
        binsize = floor(maindim/nfolds);
    end

    Perm = randperm(maindim);    
    for t = 1:nfolds
        rng shuffle;
        rowsorcols{t} = sort(Perm((t-1)*binsize+1:t*binsize));
        if t <= maindim-binsize*nfolds
            rowsorcols{t} = sort([rowsorcols{t} Perm(t+binsize*nfolds)]);
        end
        i = 1;
        for j=1:length(rowsorcols{t})
            for o = 1:otherdim
                if strcmp(MODE, 'CSROWS')
                    indices{t}(i)= sub2ind(size(MTX),rowsorcols{t}(j),o);
                elseif strcmp(MODE, 'CSCOLS')
                    indices{t}(i)= sub2ind(size(MTX), o, rowsorcols{t}(j));
                else
                    msg = 'ColdStartTestIndicesSet: Incorrect mode specified!';
                    error(msg)
                    return;
                end
                i = i + 1;
            end
        end
    end
end