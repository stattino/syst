%% HW 2
% Exercise 1 - Marginal allocation

lambda = 1/1000 * [55 43 36 70 29 45 111]; % Intensity of arrivals
c = [5 18 14 17 16 24 70]; % Cost of spare parts
T = [8 4 14 3 14 9 25]; % Repair times
C_max = 500;
m = 10;
n = 7;
A = [];

for i = 0:m
    for j = 1:n
        A(i+1,j) = poisscdf( i , lambda(j)*T(j), 'upper') / c(j); %R_j(0)/c_j
    end
end

s = zeros(m,7);
k = 0; C = 0;

EBO = [lambda*T'];
C = [0];
while (C <= C_max)
    k = k+1;
    [val, ind] = max(max(A));
    
    A( s(k,ind) + 1 , ind) = -1;
    
    s(k+1,:) = s(k,:);
    s(k+1,ind) = s(k,ind) +1;
    
    R_ind = poisscdf( s(k, ind) , lambda(ind)*T(ind), 'upper');
    C(k+1) = C(k) + c(ind);
    EBO(k+1) = EBO(k) - R_ind;
end
%%
plot(C, EBO, 'rd')
hold on
plot([500, 500], [0 5], 'r--')
plot([C], EBO)
legend('Efficient points', 'Max-budget')
xlabel('Cost')
ylabel('EBO')
title('Efficient Curve for EBO and Cost')

set(gca,'FontSize',18,'Fontname','Helvetica','Box','off','Tickdir','out','Ticklength',[.02 .02],'xcolor',0.5*[1 1 1],'ycolor',0.5*[1 1 1]);
grid on
%% Task 2
% Dynamic programming
% KNAPSACK PROBLEM W INFINITE ITEMS?
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
lambda = 1/1000 * [55 43 36 70 29 45 111]; % Intensity of arrivals
c = [5 18 14 17 16 24 70]; % Cost of spare parts
T = [8 4 14 3 14 9 25]; % Repair times
C_max = 500;

% EBO and s start at 0 money spent
EBO = zeros(1,501);
EBO(1) = [lambda*T'];
s = zeros(501,7);

for C = 1:500
    
    f_s_z = Inf(1,7);
    x = zeros(1,7);
    spent = s(C,:)*c'; % Money spent
    
    for z = 1:7
        
        if (spent+c(z) <= C_max) % Check feasibility
            if (C-c(z)>=0) % Check feasibility
                G_n = -poisscdf( s(C-c(z)+1, z)+1 , lambda(z)*T(z), 'upper');
                f_s_z(z) = G_n + EBO(C-c(z)+1);
            end
        end
        
    end
    
    [val, ind] = min(f_s_z);
    
    if(val<EBO(C))
        EBO(C+1) = val;
        x(ind) = 1;
        s(C+1, :) = s(C-c(ind)+1, :) + x;
    else
        EBO(C+1) = EBO(C);
        s(C+1, :) = s(C, :);
    end
    
end




