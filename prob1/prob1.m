clearvars
close all

% Part (a)
syms x;
Psi2= @(x) [(x-2).*(x-3)/2;-(x-1).*(x-3);(x-1).*(x-2)/2];
p = 1.3;
psi2P = Psi2(p);
fprintf('Part (a)\n')
fprintf('Psi 2 of elem 2 at p = %f: %.4e\n',p,psi2P(2))

% Part (b)
Psi1 = @(x) [1-x; x];
A = @(x) 2*(2-x);
E = 1;
K1 = int(E*A(x)*diff(Psi1(x),x)*diff(Psi1(x),x)',0,1);
fprintf('Part (b)\n')
fprintf('K1(2,1) = %.4e\n',K1(2,1))

% Part (c)
c = 2;
b = 1;
omega = 3;
A2 = c*b;
f2 = -omega*A2; %Weight per unit length at the second element
F1 = [-5.0;-4.0];
F = zeros(4,1);
F2 = int(f2*Psi2(x),1,3); %Local load vector of the second element

F(1:2) = F1;
F(2:4) = F(2:4) + F2;
fprintf('Part (c)\n')
fprintf('F(2) = %.4e\n',F(2))

% Part (d)
K2 = A2*E*int(diff(Psi2(x),x)*diff(Psi2(x),x)',1,3); %Local stiffness matrix 
                                                     %of element Omega2
K = zeros(4);
Q = zeros(4,1);
K(1:2,1:2) = K1;
K(2:4,2:4) = K(2:4,2:4) + K2; %global stifness matrix                                     
u = zeros(4,1);

fixedNods = 1;
freeNods = setdiff(1:4,fixedNods);

%set the boundary conditons, B.C.
%Natural B.C.

Q(freeNods)=0; 

%Essential B.C.
u(1)=0;

fixedNods = [fixedNods, 3, 4];
freeNods = setdiff(1:4,fixedNods);
u(fixedNods)=[u(1);-59/6;-34/3];

%Reduced system
Fm = F(freeNods) + Q(freeNods) - K(freeNods,fixedNods)*u(fixedNods);
Km = K(freeNods,freeNods);
um = Km\Fm;
u(freeNods)=um;

fprintf('Part (d)\n')
fprintf('u(2) = %.4e\n',u(2))

%Part (e)
%Post process
Q = K*u - F; %Q (= reaction forces at node 1) 
fprintf('Part (e)\n')
fprintf('Q(1) = %.4e\n',Q(1))