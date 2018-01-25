

function [ SUB_SE_MTX ] = GenerateCCAMatrix( u, v, d, firstK, lastK, lastPossK )
    if lastPossK < firstK
        SUB_SE_MTX = u(:,1:lastPossK) * diag(d(:,1:lastPossK)) * v(:,1:lastPossK)';
    else
        last = min(lastK, lastPossK);
        SUB_SE_MTX = u(:,firstK:last) * diag(d(:,firstK:last)) * v(:,firstK:last)';
    end
end

