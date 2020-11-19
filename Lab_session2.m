% WITH CSIT: COMPUTING THE WATER FILLING ALGORITHM

% VARIABLES
M = 3;
N = 3;
K = 10;
H_los = ones(N,M);
I_N = eye(N);
I_M = eye(M);
P = 10;
Q = P/M*I_M;

L = []; % List containing the values of capacity with no CSIT
L_CSIT = []; % List containing the values of capacity with CSIT

for k = 1:1:10000 % Computing the Monte Carlo simulation (10.000 iterations)
    
    % GENERATING A CHANNEL
    H_r = randn(N,M) + i*randn(N,M);
    H = sqrt(K/(K+1))*H_los + sqrt(1/(K+1))*H_r;
    H_H = ctranspose(H);
    
    % WITH NO CSIT
    C = log2(det(I_N+mtimes(H,mtimes(Q,H_H))));
    % Sometimes, we need this code line to correctly compute the matrix product
    %C = log2(det(I_N+H.*Q.*H_H));
    L = [L C];
    
    % WITH CSIT
    epsilon = 0.01;
    Q_CSIT = zeros(M,M);
    u=0;
    A = zeros(M,M);
    [U,S,V]=svd(H);
    while trace(Q_CSIT)<P-epsilon % WATER FILLING ALGORITHM
        A(1,1)=max(0,u-1/S(1,1)^2);
        A(2,2)=max(0,u-1/S(2,2)^2);
        A(3,3)=max(0,u-1/S(2,2)^2);
        Q_CSIT=mtimes(V,mtimes(A,ctranspose(V)));
        u=u+0.01;
    end
    C_CSIT = log2(det(I_N+mtimes(H,mtimes(Q_CSIT,H_H))));
    %C_CSIT = log2(det(I_N+H.*Q_CSIT.*H_H);
    L_CSIT = [L_CSIT C_CSIT];
end

cdfplot(L_CSIT);
hold on;
cdfplot(L_CSIT);
title('MIMO Ricean case');
legend('Without CSIT', 'With CSIT');
xlabel('C (bits/channel)');
ylabel('CDF');