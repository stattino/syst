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

while (C <= C_max)
    k = k+1;
    [val, ind] = max(max(A));
    
    A( s(ind) +1 , ind) = -1;
    
    s(k+1,:) = s(k,:);
    s(k+1,ind) = s(k,ind) +1;
    
    C = C + c(ind);
    
    R_ind = poisscdf( s(k, ind) , lambda(ind)*T(ind), 'upper') / c(ind);
    EBO(k+1) = EBO(k) - R_ind;
    
end