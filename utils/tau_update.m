function [G, S] = tau_update(Sigma, radius, tau, n)
% For Kullback Leibler, set tau = 0
    P = Sigma(1:n, 1:n) - Sigma(1:n, n+1:end)/(Sigma(n+1:end, n+1:end))*Sigma(n+1:end, 1:n);
    L = chol(P)';
    value = 1;
    t1 = 0;
    if tau == 1
        t2 = 10 / max(eig(P));
    else
        e = eig(P);
        r = max(abs(e));
        t2 = (1-10^-5)*((1-tau)*r)^-1;
    end
    while abs(value)>=10^-9
        th = 0.5 * (t1 + t2);
        if tau == 0
            
            value = trace(inv(eye(n)-th*P)-eye(n)) + log(det(eye(n)-th*P)) - radius;
        end
        if tau > 0 && tau < 1
            value = trace(-1/(tau*(1-tau))*(eye(n)-(1-tau)*th*L'*L)^(tau/(tau-1)) + 1/(1-tau)*(eye(n)-(1-tau)*th*L'*L)^(1/(tau-1))+1/tau*eye(n)) - radius;
        end
        if tau == 1
            value = trace(th*L'*L*expm(th*L'*L) - expm(th*L'*L)+eye(n)) - radius;
        end
        if value > 0
            t2=th;
        else
            t1 = th;
        end
    end
    if tau == 1
        V = L*expm(th*L'*L)*L';
    else
        V = L*(eye(n)-(1-tau)*th*L'*L)^(1/(tau-1))*L';
    end
    S_xx = V + Sigma(1:n, n+1:end)/(Sigma(n+1:end, n+1:end))*Sigma(n+1:end, 1:n);
    S = [S_xx, Sigma(1:n, n+1:end); Sigma(n+1:end, 1:n), Sigma(n+1:end, n+1:end)];
    G = S(1:n, n+1:end)/(S(n+1:end, n+1:end));
end