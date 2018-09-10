function xhat = apply_kalman_gain(sys, G, y, x_0)
    A = sys.A;
    C = sys.C;
    [n,~,T] = size(G);
    xhat = zeros(n,T);
    x_prev = x_0;
    for t = 1 : T

        mu_t = [A; C*A] * x_prev;
        xhat(:,t) = G(:,:,t)*(y(:,t) - mu_t(n+1:end)) + mu_t(1:n);
        x_prev = xhat(:,t);

    end
end