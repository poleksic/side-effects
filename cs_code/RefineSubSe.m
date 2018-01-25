function [ REF_SUB_SE ] = RefineSubSe( SUB_SE, CFEAT, cutoff )

    [m, n] = size(SUB_SE);
    REF_SUB_SE = zeros(m,n);
    
    for i=1:m
        sub_freq = sum(CFEAT(:,i));
        for j=1:n
            sum_horiz = sum(SUB_SE(i,:));
            sum_vert = sum(SUB_SE(:,j));
            prob_horiz = SUB_SE(i,j) / sum_horiz;
            prob_vert = SUB_SE(i,j) / sum_vert;
            prob_freq = SUB_SE(i,j) / sub_freq;
            prob_f = SUB_SE(i,j) / (sub_freq * sum_horiz);
            % note, if SUB_SE(i,j) = 0, naturally, the probability will be zero
            REF_SUB_SE(i,j) = prob_horiz;
        end
    end
    REF_SUB_SE = (REF_SUB_SE > cutoff);
end