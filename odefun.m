function dydt = odefun(t,y)
%C is the known constant of integration
%D is the size of the obstacle + agent
global C D;

P = y(1:2);
V = y(3:4);
U = y(5:6);

%udot p and v components
udotP = -3*dot(U, V) / D;
udotV = (C + dot(U,U) / 2) / norm(V);


phat = P  / norm(P);
if norm(V) > 0
    vhat = V / norm(V);
else
    vhat = [0; 0];
    udotV = 0;
end

dydt = zeros(size(y));

dydt(1:2) = V; %pdot = v
dydt(3:4) = U; %vdot = u
dydt(5:6) = (udotP * phat) + (udotV * vhat);



end

