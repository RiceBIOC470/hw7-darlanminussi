%HW7
%GB comments
1a 90 Almost enterily correct. 1 signifies the carrying capacity 
1b 100
1c 100
1d 80 There is more to the graph. If you extend a, you will notice that the system continues to bifurcate which generates a growing number of fixed points. 
2a 30 equation is not correct. This models activation and not a toggle switch. Use [V/(1+x(2)^4)-x(1); V/(1+x(1)^4)-x(2)];
2b. 100 correctly generated plots, but the equations used are wrong. I will give full credit
2c  100 same as 2b. 
overall 86


% Problem 1: Modeling population growth
% The simplest model for a growing population assumes that each current
% individual has equal likelihood to divide, which yields a differential
% equation dx/dt = a*x where a is the division rate. It is easy to see that
% this will grow exponentially without bound. A simple modification is to
% assume that the growth rate slows done as the population reaches some
% maximum value N so that dx/dt = a*x*(1-x/N). Defining X = x/N, we have 
% dX/dt = a*X*(1-X).

a = 1;
N = 100;
span = 0:0.1:100;
rhs = @(t, x) (a*x*(1-x/N));
sol = ode23(rhs, span, 1);
plot(sol.x, sol.y,'r-','LineWidth', 3, 'MarkerSize', 18);
xlabel('Time'); ylabel('Population');
set(gca,'FontSize', 24);


% Part 1. This equation has two fixed points at 0 and 1. Explain the
% meaning of these two points.

% At 0 there is no population making it a stable point
% At 1 there is only one individual with enough resources from one making
% the population stable


% Part 2: Evaluate the stability of these fixed points. Does it depend on
% the value of the parameter a? 

colors = ['r-' 'g-' 'b-' 'm-' 'c-' 'k-' 'm-'];

figure;
hold on;
for aa = 1:7
    
N = 100;
span = 0:0.1:100;
rhs = @(t, x) (aa*x*(1-x/N));
sol = ode23(rhs, span, 1);
plot(sol.x, sol.y,colors(aa),'LineWidth', 3, 'MarkerSize', 18);
xlabel('Time'); ylabel('Population');
set(gca,'FontSize', 24);
legendInfo{aa} = ['a: ' num2str(aa)];
end

legend(legendInfo);
hold off;

% as parameter a increases the population reaches its maximum at an earlier
% time though it's possible to see that it's a very small increase

% Part 3: Write a function that takes two inputs - the initial condition x0
% and the a parameter and integrates the equation forward in time. Make
% your code return two variables - the timecourse of X and the time
% required for the population to reach 99% of its maximum value. 

% [time_series, time99] = func_q1p3(param_a, x0)

[time_series, time99] = func_q1p3(2, 5);


% Part 4: Another possible model is to consider discrete generations
% instead allowing the population to vary continuously. e.g. X(t+1) = a*
% X(t)*(1-X(t)). Consider this model and vary the a parameter in the range 0
% < a <= 4. For each value of a choose 200 random starting points  0 < x0 < 1 
% and iterate the equation forward to steady state. For each final
% value Xf, plot the point in the plane (a,Xf) so that at the end you will
% have produced a bifucation diagram showing all possible final values of
% Xf at each value of a. Explain your results. 

figure;

for aa = 0.1:0.2:4
    for ii = 1:200
        x = rand();
        Xf = x;
        counter = 1;
        while counter < 10
            Xf = aa*Xf*(1-Xf);
            counter = counter + 1;
        end        
        plot(aa, Xf, 'r.')
        hold on;
    end
    xlabel('a'); ylabel('Xf');
end
hold off;

% the function converges at 2.5 where it bifurcates

% Problem 2. Genetic toggle switches. 
% Consider a genetic system of two genes A and B in which each gene
% product represses the expression of the other. Make the following
% assumptions: 
% a. Repression is cooperative:  each promotor region of one gene has 4
% binding sites for the other protein and that all of these need to be
% occupied for the gene to be repressed. 
% b. You can ignore the intermediate mRNA production so that the product of
% the synthesis of one gene can be assumed to directly repress the other
% c. the system is prefectly symmetric so that the degradation
% times, binding strengths etc are the same for both genes. 
% d. You can choose time and concentration scales so that all Michaelis
% binding constants and degradation times are equal to 1. 
%
% Part 1. Write down a two equation model (one for each gene product) for
% this system. Your model should have one free parameter corresponding to the
% maximum rate of expression of the gene, call it V. 
%

% dz/dt = V*x(2).^4./(1+x(2).^4)-x(1));

eq = @(t,x) [x(2); (5*x(2).^4)./(1+x(2).^4)-x(1)];


% Part 2. Write code to integrate your model in time and plot the results for V = 5 for two cases, 
% one in which A0 > B0 and one in which B0 > A0. 


% A0> B0
V = 5;
sol = ode23(eq, [0 1], [1; 3]);
plot(sol.x, sol.y(1,:),'g-','LineWidth', 3, 'MarkerSize', 18); hold on;
plot(sol.x, sol.y(2,:),'r-','LineWidth', 3, 'MarkerSize', 18);
xlabel('Time'); ylabel('Expression');

%B0 > A0
V = 5;
sol = ode23(eq, [0 1], [3; 1]);
plot(sol.x, sol.y(1,:),'g-','LineWidth', 3, 'MarkerSize', 18); hold on;
plot(sol.x, sol.y(2,:),'r-','LineWidth', 3, 'MarkerSize', 18);
xlabel('Time'); ylabel('Expression');

%
% Part 3. By any means you want, write code to produce a bifurcation diagram showing all
% fixed points of the system as a function of the V parameter. 

ku = 0;
kb = 3;

figure;
hold on;
for V = 0:0.05:5
    polycoeff = [1 V*(-kb) 1 V*(-ku)];
    rts = roots(polycoeff);
    rts = rts(imag(rts) == 0);
    plot(V*ones(length(rts),1), rts, 'r.');
end
hold off;

