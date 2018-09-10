function [x,y,y0] = generate_data(sys, x_0, T, coeff, is_TV)

    [m,n] = size(sys.C);
    d = m + n;

    y0 = sys.C * x_0 + sys.D * randn(d,1);
    x_prev = x_0;
    x = zeros(n,T);
    y = zeros(m,T);
    Delta = 2 * rand - 1;
    A_purt = sys.A + [0, coeff * Delta; 0, 0];
    for t = 1 : T
        x(:,t) = A_purt * x_prev + sys.B * randn(d,1);
        y(:,t) = sys.C * x(:,t) + sys.D * randn(d,1);
        x_prev = x(:,t);
        if is_TV
            Delta = 2 * rand - 1;
            A_purt = sys.A + [0, coeff * Delta; 0, 0];
        end
    end

end