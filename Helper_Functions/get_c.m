function [c, c_min, c_max] = get_c(N, K, df, f_param, P_min, P_max)
    c_min = zeros(N, 1); 
    c_max = zeros(N, 1);
    for n=1:N
        c_min(n) = min(df(n, 1:K, P_min(:,1:K), f_param));
        c_max(n) = max(df(n, 1:K, P_max(:,1:K), f_param));
    end
    c = (c_min+c_max)/2;
    c_min = c_min - 0.1;
    c_max = c_max + 0.1;