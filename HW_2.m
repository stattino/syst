%% HW 2
% Exercise 1 - Marginal allocation

lambda = 1/1000 * [55 43 36 70 29 45 111]; % Intensity of arrivals
c = [5 18 14 17 16 24 70]; % Cost of spare parts
T = [8 4 14 3 14 9 25]; % Repair times
C_max = 500;

% Marginal allocation

% R_j(s_j) = cdf(pd, s_j +1 , 'upper')
m = 10;
n = 7;
A = [];
for i = 0:m
    for j = 1:n
        A(i+1,j) = poisscdf( i+1 , lambda(j)*T(j), 'upper') / c(j);
    end
end

s = zeros(m,7);
k = 0;
C = 0;

EBO = [lambda*T'];
C = [0];
while (C <= C_max)
    k = k+1;
    [val, ind] = max(max(A));
    
    A( s(k,ind) + 1 , ind) = -1;
    
    s(k+1,:) = s(k,:);
    s(k+1,ind) = s(k,ind) +1;
    
    
    R_ind = poisscdf( s(k, ind) , lambda(ind)*T(ind), 'upper') / c(ind);
    
    C(k+1) = C(k) + c(ind);
    
    EBO(k+1) = EBO(k) - R_ind;
    
end

plot(1:k+1, EBO)
figure()
plot(1:k+1, C)

%% Task 2
% Dynamic programming
% Define problem: stages k, states s_k, and decision x_k, feasibility set
% F_k(s_k), next state fctn S_k+1 = h_k(s_k, x_k)

% The stage k is the k:th bought spare part

% The state s_k is a vector with entries for how many spare parts of each
% component is bought at the k:th bought spare part.

% x_k is a bit vector where all entries are zero except for one which is
% the part we decide to buy.

% The feasibility set is all the avail. components s.t the cost is less
% than C_max

% The function s_k+1 = h_k(s_k, x_k) = s_k + x_k
N = 0; C = 0;
EBO = [lambda*T'];
s = zeros(1,7);
k = 0;
f_n = -EBO;
while (C <= C_max)
    k = k+1;
    
    f_s_z = zeros(1,7);
    x = zeros(1,7);
    
    for z = 1:7
        
        if (C + c(z) <= C_max) % Check feasibility
            G_n = -poisscdf( s(k, z) , lambda(z)*T(z), 'upper') / c(z);
            f_s_z(z) = G_n + f_n ;
        end
        
    end
    
    [val, ind] = min(f_s_z);
    f_n = f_n + val;
    x(ind) = 1;
    s(k+1, :) = s(k, :) + x;
    
    C = C + c(ind);
    
end

