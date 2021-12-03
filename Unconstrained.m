function [p, v, u, udot] = Unconstrained(p0, pf, v0, t)
%this solves for the unconstrained trajectory given p0, pf, v0, and uf=0

%get t0 and tf
t0 = t(1);
tf = t(end);

%solve the linear equations A*C = b
bx = [p0(1); v0(1); pf(1); 0];
by = [p0(2); v0(2); pf(2); 0];

%build the matrix A
A = [t0^3,   t0^2,    t0, 1; ... %p0
     3*t0^2, 2*t0,     1, 0; ... %v0
     tf^3,   tf^2,    tf, 1; ... %pf
     6*tf,      2,     0, 0; ... %uf
    ];

%do linear algebra to find the coefficients in x and y
cx = A\bx;
cy = A\by;


%generate trajectory
px = cx(1)*t.^3 + cx(2)*t.^2 + cx(3)*t + cx(4);
py = cy(1)*t.^3 + cy(2)*t.^2 + cy(3)*t + cy(4);

p = [px; py];

vx = 3*cx(1)*t.^2 + 2*cx(2)*t + cx(3);
vy = 3*cy(1)*t.^2 + 2*cy(2)*t + cy(3);

v = [vx; vy];

ux = 6*cx(1)*t + 2*cx(2);
uy = 6*cy(1)*t + 2*cy(2);

u = [ux; uy];

udotx = 6*cx(1)*ones(size(t));
udoty = 6*cy(1)*ones(size(t));

udot = [udotx; udoty];

end

