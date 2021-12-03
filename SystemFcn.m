function [c_neq, c_eq] = SystemFcn(x)
%this solves for the trajectory of the system when the constrained motion
%primitive is active over a non-zero interval of time

%parameters
global D C total_time P0 V0 Pf;

global PLOT;

%inputs are (pi-theta, t0, t2)
angle = x(1);
t0 = x(2);
t2 = x(3);

t1 = 0; %we assume the obstacle is contacted at t1 = 0.
tf = total_time + t0;

%rotate coordinate system by -angle
R = [cos(-angle), -sin(-angle); sin(-angle), cos(-angle)];
p0 = R*P0;
pf = R*Pf;

v0 = R*V0;

%in the rotated coordinate system p(t1) = p(0) = d = [-D, 0], and
%  v \cdot \hat{p} = c_x d_x + c_y d_y = 0.
d = [-D; 0];
cx = 0;

%first solve bx
bx = 3/(t0^2)*(p0(1)+D) - (v0(1)/t0);
%cy from bx
cy = sqrt(2*D*bx);
c = [cx; cy];
%by from cy
by = 3/(t0^2)*p0(2) - v0(2)/t0 - 2/t0*cy;
b = [bx; by];

a = (v0 - c)/(3*t0^2) - 2*b/(3*t0);

%calculate the trajectory
T1 = linspace(t0, t1);
P1 = a*T1.^3 + b*T1.^2 + c*T1 + d;
%V1 = 3*a*T1.^2 + 2*b*T1;
U1 = 6*a*T1 + 2*b;

if PLOT
    figure(1); clf; hold on;
    plot(P1(1,:), P1(2,:), 'b', 'linewidth', 3)
    rectangle('position', [-D, -D, 2*D, 2*D], 'curvature', 1);
    axis equal
end

%calculate C from first unconstrained arc
C = 6*dot(a, c) - 2*dot(b, b);

%calculate t for t2
tspan = [t1, t2];


%px, py, vx, vy, ux, uy
p1 = a*t1^3 + b*t1^2 + c*t1 + d;
v1 = 3*a*t1^2 + 2*b*t1 + c;
u1 = 6*a*t1 + 2*b;

%IC for the constrained ODE
y0 = [p1; v1; u1];

% solve for arc along circle
[T2, Y2] = ode45(@odefun, tspan, y0);
if PLOT
    plot(Y2(:,1), Y2(:,2), 'r', 'linewidth', 3)
end


%solve for second unconstrained trajectory using the exit states from the
%constrained arc and the final conditions
p2 = Y2(end,1:2)';
v2 = Y2(end,3:4)';
%u2 = Y2(end,5:6)';

T3 = linspace(T2(end), tf);

[P3, V3, U3, Udot3] = Unconstrained(p2, pf, v2, T3);


if PLOT
    plot(P3(1,:), P3(2,:), 'g', 'linewidth', 3)
    camroll(rad2deg(angle))

    cv1 = max(max(D*D - dot(P1,P1), 0));
    
    PC = Y2(:,1:2)';
    
    cv2 = max(max(D*D - dot(PC, PC), 0));
    cv3 = max(max(D*D - dot(P3,P3), 0));
    
    T = [T1'; T2; T3'];
    U = [U1'; Y2(:,5:6); U3'];
    UX = U(:,1);
    UY = U(:,2);
    
    UM = UX.^2 + UY.^2;
    
    figure(100);
    plot(T, UM)
    title('Cost vs Time');
    
    cost = trapz(T, UM);
    
    fprintf('Max violation at grid points: %g\nCost: %g\n', ...
                max([cv1,cv2,cv3]), cost);
    
end

% we have 3 unknowns, these are the 3 corresponding equations that we have 
% not used yet - continuity in U and continuity in \dot{U} \dot V.
c1 = Y2(end,5) - U3(1,1);
c2 = Y2(end,6) - U3(2,1);
c3 = C + 1/2*dot(Y2(end,5:6),Y2(end,5:6)) - dot( Udot3(:,1), V3(:,1) );


%make sure we don't violate safety or time constraints
p1_const = trapz(T1, max(D*D - dot(P1,P1), 0));
p2_const = trapz(T3, max(D*D - dot(P3,P3), 0));
t_constraint = t2 - tf;


c_neq = [p1_const, p2_const, t_constraint];
c_eq = [c1; c2; c3];

end