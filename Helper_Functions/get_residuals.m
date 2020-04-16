function residuals = get_residuals(N, K, df, Q, f_param, P_max, P_min, dw, c)

residuals = zeros(1, N);
for n=1:N      
    u = min( df(n, 1:K, Q, f_param) - c(n) , 0);
    v = max( df(n, 1:K, Q, f_param) - c(n) , 0);
        
    residuals(n) = (Q(n,:)-P_max(n,:))*u'*dw + (Q(n,:)-P_min(n,:))*v'*dw;
end