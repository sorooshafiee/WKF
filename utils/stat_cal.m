function [xx, yy, xy, xw, ww] = stat_cal(S)
% For the system y = x + w
% The joint least-favorable covariance S can be decomposed as  S = [xx, xy; yx, yy]
% Variance of signal: xx
% Variance of observation: yy
% Variance of noise: ww
% Covariance between noise and signal: xw

    xx = S(1, 1);
    xy = S(1, 2);
    yy = S(2, 2);
    
    xw = xy - xx;
    ww = yy -xx - 2*xw;
    
    if ww < 0
        display('Error: negative noise variance');
    end
end

