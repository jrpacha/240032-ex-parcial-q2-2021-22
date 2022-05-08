clearvars
close all

x = 0.337;
dudx = 1.0;
a = -1.0;
alpha = 2.0;

numDiv = 100;

[minU,interpU] = funcioProb3(dudx,a,alpha);

%Part A
fprintf('Part A\n')
fprintf('For a = -1, alpha = 2.0, the minumum of the nodal solution\n')
fprintf('    is min {u(i), i = 1,...,%d} = %.4e\n',numDiv+1,minU)

minUHint = funcioProb3(dudx,0,alpha);
fprintf('Hint. For a = 0, alpha = 2.0, min {u(i), i = 1,...,%d}\n',...
    numDiv+1)
fprintf('    = %.4e\n\n',minUHint)

%Part B
fprintf('Part B\n')
fprintf('For a = -1, alpha = 2.0, the interpolarted  value of u, U,\n')
fprintf('    at x = %f is U = %.4e\n\n',x,interpU)

%Part C
[alpha,u] = funcioProb3(dudx,a);
fprintf('Part C\n')
fprintf('For a = -1, the value of alpha s.t. u(27)-u(70) = 2u(52) is\n')
fprintf('    alpha = %.4e\n',alpha)
fprintf('Check: u(27) - 2u(52) - u(70) = %.4e\n',u(27)-2*u(52)-u(70))