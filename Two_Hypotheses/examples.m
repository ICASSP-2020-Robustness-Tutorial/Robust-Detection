% Example for least favorable densities under band uncertainty

% Sample space
dx = 0.01;
w = -8:dx:8;

% nominal densities
p0 = normpdf(w, -2, 2);
p1 = normpdf(w, 1, 1);

% LFDs under density band uncertainty

% bands
p_min(1, :) = 0.8 * p0;
p_min(2, :) = 0.8 * p1;
p_max(1, :) = 1.2 * p0;
p_max(2, :) = 1.2 * p1;

% solve for LFDs
[q0, q1, llr, c, nit] = lfds_density_band(p_min, p_max, dx);

% plot lfds
figure;
plot(w, q0, w, q1)
legend('q_0', 'q_1')
title('Density band uncertainty')

% plot log-likelihood ratio
figure;
plot(w, llr)
legend('log-likelihood ratio')
title('Density band uncertainty')

% LFDs under 10% and 15% contamination (outliers)
 
eps = [0.1, 0.15];
 
% solve for LFDs
[q0, q1, llr, c] = lfds_outliers(p0, p1, dx, eps);

% plot lfds
figure;
plot(w, q0, w, q1)
legend('q_0', 'q_1')
title('Epsilon Contamination')

% plot log-likelihood ratio
figure;
plot(w, llr)
legend('log-likelihood ratio')
title('Epsilon Contamination')
