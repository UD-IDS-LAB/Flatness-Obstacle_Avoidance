clear; close all; clc;
global D total_time P0 V0 Pf;

% system setup
D = 1.25;         % obstacle + agent radius
total_time = 10;  % trajectory time
P0 = [-2; -2];    % initial position
Pf = [1; 2];      % final position
V0 = [0; 0];      % initial speed

% PLOT = true shows the generated trajectory and prints data about the
% initial/final states and optimality. Set it to false to decrease runtime,
% then set it to true after finding a solution to view it.
global PLOT;
PLOT = false;

%create 100 values of t for the plot 
%  this is used to plot and to check constraint violations, the solver o
%  only needs to use t(1) and t(end).
t = linspace(0, total_time);
%make an initial plot of the unconstrained trajectory
[p, ~, ~, ~] = Unconstrained(P0, Pf, V0, t);
figure(10); clf; hold on;
plot(p(1,:), p(2,:), 'linewidth', 3);
rectangle('Position', [-D, -D, 2*D, 2*D], 'curvature', 1);
axis equal;


tic % start timing!
%solve unconstrained
[p, v, u, udot] = Unconstrained(P0, Pf, V0, t);
%determine the distance to the obstacle
dist = dot(p,p);
%find constraint violations along the trajectory
idx = find(dist <= D^2);
%if there are no constraint violations end early
if isempty(idx)
    toc
    plot(p(1,idx(1)), p(2,idx(1)), 'or');
    plot(p(1,idx(end)), p(2,idx(end)), 'or');
    return
end


%find the time of first constraint violation, that will be our guess for t1
t1_unc = t(idx(1));

%find the angle that we violated the constraint, that will be our initial
%      guess for theta
thIn = atan2(p(2,idx(1)), p(1,idx(1)));

%make theta_in positive
if thIn < 0
    thIn = thIn + 2*pi;
end

%the solver assumes d = [-D^2, 0], so we want angle = pi - theta.
angle = pi - thIn; %pi - thIn;
%we set t_1 = 0 to cancel out the nonlinearities in the constraint
%activation. As a consequence t0 becomes negative.
t0 = -t1_unc;

%here our unknowns are (pi-theta) and t^f. The condition u(t^f) = 0 is a
%     bit easier to solve than the equations for t_1 and t^0.
x0 = [pi - thIn;   total_time - t1_unc];
%we also need to lower and upper bound tf to ensure 
%   that t^0 < t_1 = 0 < t^f
lb = [-inf, 0];
ub = [ inf, total_time];

%use the default algorithm and set display to off to speed it up a bit
opts = optimoptions(@fmincon,'Algorithm','interior-point', ... 
    'Display','off');%,'StepTolerance',eps);

%we are solving a constrained system of equations using fmincon, see the
%docs: https://www.mathworks.com/help/optim/ug/nonlinear-systems-with-constraints.html
sln = fmincon(@(x) 0, x0, [], [], [],[], lb, ub, @Instantaneous, opts);

%now we evaluate the solution and check for constraint violations
[c_neq, c_eq] = Instantaneous(sln);

% if there are no constraint violations end early
if max(c_neq) == 0
    toc
    PLOT = true;
    [c_neq, c_eq] = Instantaneous(sln);
    PLOT = false;
    return;
end



%finally, we select t2 and re-use the values of theta and t1 from before
t2 = 1e-3;
%in this case we use t0 and t_2 as our unknown variables and constrain them
%accordingly to maintain t0 < t1 = 0 < t2 < tf
x0 = [angle; t0; t2];
lb = [-inf, -total_time, 0];
ub = [ inf, 0, total_time];


opts = optimoptions(@fmincon,'Algorithm','interior-point', ... 
    'Display','off');%,'StepTolerance',eps);
sln = fmincon(@(x) 0, x0, [], [], [],[], lb, ub, @SystemFcn, opts);

%if we reach this far, print the time. This is the best solution we could
%find.
toc


angle = sln(1);
t0 = sln(2);
t2 = sln(3);
sln;

PLOT = true;
[~, ~] = SystemFcn(sln);
PLOT = false;



