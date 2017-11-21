function [time_series, time99] = func_q1p3(param_a, x0)

%input the initial condition x0
% and the a parameter

%output the timecourse of X and the time
% required for the population to reach 99% of its maximum value

soly = 0;
solx = 0;

while soly < 0.99
N = 100;
span = 0:0.1:100;
rhs = @(t, x) (param_a*x*(1-x/N));
sol = ode23(rhs, span, x0);
soly = sol.y(end);
solx = sol.x;
end

time99 = length(sol.y);
time_series = length(solx);
