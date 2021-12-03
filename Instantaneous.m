function [c_neq, c_eq] = Instantaneous(x)
%parameters
global D total_time P0 V0 Pf;

global PLOT;

%extract state values and calculate times
angle = x(1);
tf = x(2);
t1 = 0; %we assume the obstacle is contacted at t1 = 0.
t0 = tf - total_time;

%rotate the entire coordinate system by -angle = -(pi-theta)
R = [cos(-angle), -sin(-angle); sin(-angle), cos(-angle)];

p0 = R*P0;
pf = R*Pf;
v0 = R*V0;


%solve for the coefficients (via InstantaneousEquations.mlx)
d = [-D; 0];
a1 = -(3*d*t0^2 + 2*d*tf^2 - 2*p0*tf^2 - 3*pf*t0^2 + 2*t0*tf^2*v0 - 3*t0^2*tf*v0 - 6*d*t0*tf + 6*p0*t0*tf)/(t0^3*tf*(3*t0 - 4*tf));
a2 = (3*d*tf - 2*d*t0 - 3*p0*tf + 2*pf*t0 + t0*tf*v0)/(t0*tf^2*(3*t0 - 4*tf));
b = -(3*(3*d*tf - 2*d*t0 - 3*p0*tf + 2*pf*t0 + t0*tf*v0))/(t0*tf*(3*t0 - 4*tf));
c = (6*d*tf^2 - 3*d*t0^2 - 6*p0*tf^2 + 3*pf*t0^2 + 2*t0*tf^2*v0)/(t0*tf*(3*t0 - 4*tf));


%calculate the first and second unconstrained arcs
T1 = linspace(t0, t1);
P1 = a1*T1.^3 + b*T1.^2 + c*T1 + d;
V1 = 3*a1*T1.^2 + 2*b*T1 + c;
U1 = 6*a1*T1 + 2*b;

T2 = linspace(t1, tf);
P2 = a2*T2.^3 + b*T2.^2 + c*T2 + d;
%V2 = 3*a2*T2.^2 + 2*b*T2 + c;
U2 = 6*a2*T2 + 2*b;


%print some info if PLOT has been set to true in Main.m
if PLOT
    figure(1); clf; hold on;
    plot(P1(1,:), P1(2,:), 'b', 'linewidth', 3)
    plot(P2(1,:), P2(2,:), 'r', 'linewidth', 3)
    rectangle('position', [-D, -D, 2*D, 2*D], 'curvature', 1);
    axis equal
    

    fprintf('P1 is (%g, %g)\n', P1(:,1));
    fprintf('V1 is (%g, %g)\n', V1(:,1));
    fprintf('Pf is (%g, %g)\n', P2(:,end));
    fprintf('Uf is (%g, %g)\n', U2(:,end));

   cost_1 = trapz(T1, sum(U1.*U1));
   cost_2 = trapz(T2, sum(U1.*U2));
   
   fprintf('Total Cost: %g\n', (cost_1 + cost_2) / 2);
    
  %plot it in the original coordinates
    figure(2); clf; hold on;
    P = [P1, P2]; %aggregate Ps
    %rotate back
    for i = 1:length(P)
        P(:,i) = R'*P(:,i);
    end
   
    pTouch = R'*P1(:,end);
    
    %linewidth
    lw = 3;
    
    %plot it!
    rectangle('position', [-D, -D, 2*D, 2*D], 'curvature', 1, ...
       'linewidth', lw, 'edgecolor', 'r', 'facecolor', [0.95, 0.85, 0.85]);
   
   
    plot(P(1,:), P(2,:), '-k', 'linewidth', lw)
    
    plot(P(1,1), P(2,1), '^k', 'markerfacecolor', 'k', 'markersize', 10);
    plot(P(1,end), P(2,end), 'vk', 'markerfacecolor', 'k', 'markersize', 10);

        
    thetaLine = pTouch / sqrt(sum(pTouch.^2)) * (D+0.5);
    
    
    plot([D+0.5, 0, thetaLine(1)], [0, 0, thetaLine(2)], ':b', ...
        'linewidth', lw - 0.5);
    
    text(0.02, 0.25, '\theta', 'FontSize', 18);
    
    th_arc = linspace( 0, atan2(thetaLine(2), thetaLine(1)));
    r_arc = 0.1;  %or whatever radius you want
    x = r_arc*cos(th_arc);
    y = r_arc*sin(th_arc);
    plot(x, y, 'b', 'linewidth', lw - 0.5);
    %text(pTouch(1) - 0.1, pTouch(2) + 0.3, 't_1', 'FontSize', 16);
    
    grid on; box on; axis equal;
    xlabel('x (m)');
    ylabel('y (m)');
    set(gca,'FontSize', 12, 'FontName', 'Times')

   
end

% these are the 2 equations that our two unknowns must satisfy
%make sure a1_y = a2_y (from \dot{u}^+\cdot v = \dot{u}^- \cdot v)
c1 = a1(2) - a2(2);
%make sure dot(p,v) = 0
c2 = dot(P1(:,end),V1(:,end));

%make sure sol'n is feasible
p1_const = trapz(T1, max(D*D - dot(P1,P1), 0));
p2_const = trapz(T2, max(D*D - dot(P2,P2), 0));

%inequality and equality constraints
c_neq = (p1_const + p2_const);
c_eq = [c1; c2];


end