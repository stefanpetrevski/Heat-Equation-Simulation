function [errorM] = heat_eqn1error(T,M,lambda)

%set up grid points
h = 2/M; % grid spacing for splitting the x-axis
x = (linspace(-1,1,M+1))'; % splitting the interval [-1,1]
k = lambda*h^2; % grid spacing for going ahead in time
N = ceil(T/k);

%initial conditions
u = cos(pi*x/2) + (sin(pi*x))/2; % ux,0)
t = 0; % starting at 0 time

%finite difference operator
e = ones(M-1,1); % generating an (M-1)x1 column vector of 1s only
L = spdiags([e  -2*e  e], [-1 0 1], M-1, M-1); 
% tridiagonal and sparse, with main diag of -2s, adjacent of -1s
I = speye(M-1); % sparse identity matrix
A = I + lambda*L; % A as defined in the problem

%main loop
while (t<T)
    u_inner = u(2:M,1); % gets the M-1 inner components of U^n
	u_inner = A*u_inner; % calculating solution forward in time
    u = [0; u_inner; 0]; % adding the boundary conditions
    t = t + k; % going ahead in time by the specified amount
end

u_exact = exp(-((pi^2)*t)/4)*cos(pi*x/2)+(exp(-(pi^2)*t)*(sin(pi*x)))/2;
% actual solution at endpoint as after loop, t is the final time

errorM = norm(u - u_exact, inf); % the infinite norm = error

end