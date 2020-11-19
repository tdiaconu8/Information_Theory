% NO CSIT

% VARIABLES
M = 1;
N = 1;
K = 10;
H_los = ones(N,M);
I_N = eye(N);
I_M = eye(M);
P = 10;
Q = (P/M)*I_M;

L = []; % List containing the values of capacity
L_inf = []; % List containing the values of capacity when K is infinite

for k = 1:1:10000 % Computing the Monte Carlo simulation (10.000 iterations)
    
    % GENERATING A CHANNEL
    H_r = randn(N,M) + i*randn(N,M);
    H = sqrt(K/(K+1))*H_los + sqrt(1/(K+1))*H_r;
    H_H = ctranspose(H);
    % Sometimes, we need this code line to correctly compute the matrix product
    %C = log2(det(I_N+H.*Q.*H_H)); 
    C = log2(det(I_N+mtimes(H,mtimes(Q,H_H))));
    % We compute the capacity for K = infinite, H = H_los
    C_inf = log2(det(I_N+mtimes(H_los,mtimes(Q,ctranspose(H_los)))));
    L = [L C];
    L_inf = [L_inf C_inf];
end


cdfplot(L);
hold on;
cdfplot(L_inf);
title('SIMO fading case');
legend('Rayleigh fading channel', 'AWGN');
xlabel('C (bits/channel)');
ylabel('CDF');