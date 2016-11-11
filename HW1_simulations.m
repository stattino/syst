%% Homework 1. Simulation of a MC
clc
clear all
%% Ex. 2
% Q - intensity matrix, T - months simulated
% mu - parameters of exponential distributions, 1/q_ij
% t initialized at 0, state initialized at 0
Q = [-1/3 1/3; 2 -2];
T = 1000;
mu = [3, 1/2];
t = 0; state = 0;
timevec = [t]; statevec = [state];
while (t < T)
    t_inc = exprnd(mu(state+1));
    state = mod(state+1,2);
    t = t + t_inc;
    timevec = [timevec; t];
    statevec = [statevec; state];
end

timedif = timevec(2:end) - timevec(1:end-1);

time_working = timedif(1:2:end);
time_broken = timedif(2:2:end);

MTTR = mean(time_broken)
MTTF = mean(time_working)
MTBF = MTTF+ MTTR


%% Plots ex. 1
statevec_2 = zeros(size(statevec,1)*2,1);
for i = 1:size(statevec,1)
    statevec_2(2*i) = statevec(i);
    statevec_2(2*(i-1)+1) = statevec(i);
end


time_2 = zeros(size(timevec,1)*2,1);
for i = 1:size(timevec,1)
    time_2(2*i) = timevec(i);
    time_2(2*(i-1)+1) = timevec(i);
end

plot(time_2(2:end), statevec_2(1:end-1))


%% Ex 3. Discrete time 
% States 0 and 1.
Q = [-1/3 1/3; 2 -2];
n = 1/100;
N = 100000;
P =    [ 1+n*Q(1,1), n*Q(1,2); ...
        n*Q(2,1), 1+n*Q(2,2)];

state = 0;
statevec = [state];
sw = zeros(2,2);
for t = 1:N
    x = rand;
    switch statevec(t)
        case 0
            if x > P(1,1)
                state = mod(state+1,2);
                sw(1,1) = sw(1,1) +1;
            else
                sw(1,2) = sw(1,2) +1;
            end
        case 1
            if x > P(2,2)
                state = mod(state+1,2);
                sw(2,1) = sw(2,1) +1;
                
            else
                sw(2,2) = sw(2,2) +1;
            end
    end
    
    statevec = [statevec state];
end

plot(1:N+1, statevec)

% sw-matrix contains in (1,1) how many times it stays in state 0, (1,2) how
% many times it leaves it. Simliarily (2,1) how many times stays in state
% 1, and (2,2) leaving it. 

MTTF = sw(1,2) / sw(1,1) * n
MTTR = sw(2,2) / sw(2,1) * n
%MTBF = MTTF+ MTTR

%% Ex 4 version 2.0
% Definierar två states: 0 och 1, för trasig och helt

mu = [3, 1/2];
t = 0; state = 0;
timevec = [t]; statevec = [state];
while (t < T)
    switch statevec(end)
        case 0
            components = exprnd(3, 3, 1); %  simulates 3 components
            t_inc = min(components) ;
            state = state + 1;
            t = t+ t_inc;
        case 1
            repair = exprnd(1/2);
            t = t + repair;
            state = state - 1;
    end
    timevec = [timevec; t];
    statevec = [statevec; state];
end

timedif = timevec(2:end) - timevec(1:end-1);

time_working = timedif(1:2:end);
time_broken = timedif(2:2:end);

MTTR = mean(time_broken)
MTTF = mean(time_working)
%MTBF = MTTF+ MTTR

%% Plots series

%% Ex. 5 n=2 components in parallel
% Two components in parallel. 
% 3 states -> state 0: all working, state 1: one broken, state 2: two
% broken. Then this simplifies as state (1 U 2): working, state (3): broken.
Q = [   -1/3, 1/3, 0, 0; ...
        2, -7/3, 1/3, 0; ...
        0, 0, 2, -2     ];

T = 10000; t = 0; state = 0;

timevec = [t];
statevec = [state];
no_failures = 0;
while (t < T)
    switch (statevec(end))
        case 0
            up = exprnd(3);
            t = t + up;
            state = state + 1;
        case 1
            up = exprnd(3);
            down = exprnd(1/2);
            if up < down
                state = state + 1;
                t = t + up;
            else
                state = state - 1;
                t = t + down;
            end
        case 2 
            down = exprnd(1/2);
            t = t + down;
            state = state - 1;
    end

    statevec = [statevec state];
    timevec = [timevec t];
end
% Broken when in state 2, 

%

statevec_broken = (statevec==2); 
statevec_working = (statevec~=2);

timedif = timevec(2:end) - timevec(1:end-1);

timedif_broken = statevec_broken(2:end).*timedif;
timedif_working = statevec_working(2:end).*timedif;

total_time_broken = sum(timedif_broken);
number_time_broken = sum(timedif_broken~=0); 
MTTR = total_time_broken/number_time_broken

total_time_working = sum(timedif_working);
MTTF = total_time_working/(number_time_broken+1)
%MTBF = MTTF+ MTTR
%% Plots parallel 
plot(timevec, statevec_broken)
figure()
plot(timevec, statevec)


