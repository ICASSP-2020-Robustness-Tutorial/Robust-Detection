% Example for least favorable densities under band uncertainty

% Sample Space
wmin = -10;
wmax = 10;
dx = 0.01;
x = wmin : dx : wmax;

% Density Bands
N = 3;
p1 = normpdf(x, -2, 3);
p2 = normpdf(x, 2, 2);
p3 = normpdf(x, 0, 1);
P = [p1; p2; p3];
Pmin = 0.8*P;
Pmax = 1.2*P;

% Objective Function parameter
a = [0.5 0.5];

% The objective function itself and its derivative are defined at the end of the file


% LFDs under density band uncertainty

% Denisty band model using regular algorithm
Q = multi_lfds_density_band(@f, @df, a, Pmin, Pmax, dx, true);

% Plot lfds
figure; plot(x,Q)
legend('q_0', 'q_1', 'q_2')
title('Density band uncertainty - Regular')

% Denisty band model using proximal algorithm
Q = multi_lfds_density_band_proximal(@f, @df, a, Pmin, Pmax, dx, true);

% Plot lfds
figure; plot(x,Q)
legend('q_0', 'q_1', 'q_2')
title('Density band uncertainty - Proximal')


% LFDs under 10%  contamination (outliers)

% outlier model using regular algorithm
Q = multi_lfds_outlier(@f, @df, a, P, dx, 0.1, true, P);

% Plot lfds
figure; plot(x,Q)
legend('q_0', 'q_1', 'q_2')
title('Epsilon contamination - Regular')

% outlier model using proximal algorithm
Q = multi_lfds_outlier_proximal(@f, @df, a, P, dx, 0.1, true, P);

% Plot lfds
figure; plot(x,Q)
legend('q_0', 'q_1', 'q_2')
title('Epsilon contamination - Proximal')


% Objective function
function val = f(k, x, f_param)
    a = f_param;
    val = 0;
    for n=1:2
        val = val + a(n)*log(x(3,k)./x(n,k)).*x(3,k);
    end
end


% Partial derivatives
function val = df(n, k, x, f_param)
    a = f_param;
    if n == 1 || n == 2
        val = -a(n)*x(3,k)./x(n,k);
    elseif n == 3
        val = a(1)*log(x(3,k)./x(1,k)) + a(2)*log(x(3,k)./x(2,k));
    end
end