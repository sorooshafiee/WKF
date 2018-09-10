function [] = plotErrorEllipse(mu, Sigma, p, txt)
    s = -2 * log(1 - p);

    [V, D] = eig(Sigma * s);

    t = linspace(0, 2 * pi);
    a = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];
    %disp(['plot(a(1, :) + mu(1), a(2, :) + mu(2)', txt, ');']);
    eval(['plot(a(1, :) + mu(1), a(2, :) + mu(2)', txt, ');']);
end

