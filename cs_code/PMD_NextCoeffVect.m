% Compute solution u (along with corresponding Delta) given in Witten's Lemma 2.2; 
% use maxiter ~ 200 and precision ~ 1e-6
% normS is normalized function S(a,c) from the paper
function [ u ] = PMD_NextCoeffVect(a, c, maxiter, precision)
    u = normS(a,0);
    if norm(u,1) <= c + precision
        return
    end

    Delta1 = 0;
    Delta2 = max(abs(a));
    iter = 1;
    while iter < maxiter
        Delta = (Delta1+Delta2)/2;
        u = normS(a,Delta);
        L1 = norm(u,1);
        if abs(L1 - c) < precision
            return
        elseif L1 < c
            Delta2 = (Delta1+Delta2)/2;
        else 
            Delta1= (Delta1+Delta2)/2;
        end
        if (Delta2 - Delta1) < precision
            Delta=(Delta1+Delta2)/2;
            u = normS(a,Delta);
            return
        end
        iter = iter + 1;
    end
    Delta=(Delta1+Delta2)/2;
    u = normS(a,Delta);
    return
end

function [normS] = normS(a,delta)
    S = sign(a) .* max(0, abs(a)-delta);
    if sum(S) ~= 0
        normS = S ./ norm(S,2);
    else
        warning('normS: Zero vector given in input')
        normS = S;
    end
end
