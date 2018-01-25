% This is Algorithm 1 (with details specified in Algorithm 3) from Witten
% X is nxp matrix 
% v is px1 vector
% u is nx1 vector
function  [u, v, d, convu, convv] = PMD(X, C1, C2, initv, k, maxiter, precision)
    [n, p] = size(X);
    v = initv(:, k);
    [~, K] = size(initv);
    cntr = 0;
    while cntr < K
        v1 = zeros(p,1);
        for i= 1:maxiter
            if norm(v-v1,1)> precision
                v1 = v;
                [ u ] = PMD_NextCoeffVect(X*v, C1, maxiter, precision);
                [ v ] = PMD_NextCoeffVect(X'*u, C2, maxiter, precision);
            else
                break
            end
        end
        nu1 = norm(u,1);
        nv1 = norm(v,1);
        if nu1 < C1 + precision
            convu = 1;
        else
            convu = 0;
        end
        if nv1 < C2 + precision
            convv = 1;
        else
            convv = 0;
        end
        d =u'*X*v;
        if convu * convv > 0
            return
        end
        cntr = cntr + 1;
        if cntr ~= k
            v = initv(:, cntr);
        else
            cntr = cntr + 1;
        end
    end
end


