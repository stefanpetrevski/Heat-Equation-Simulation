function u = allen_cahn(T,M,lambda,beta)

%set up grid points
h = 2/M; % grid spacing for splitting the x-axis
x = (linspace(-1,1,M+1))'; % splitting the interval [-1,1]
k = lambda*h^2; % grid spacing for going ahead in time
N = ceil(T/k);

%initial conditions
u = cos(pi*x/2) + (sin(pi*x))/2; % u(x,0)
t = 0; 

%finite difference operator
e = ones(M-1,1); % generating an (M-1)x1 column vector of 1s only
L = spdiags([e  -2*e  e], [-1 0 1], M-1, M-1); 
I = speye(M-1); % sparse identity matrix
C = inv(I - (lambda/2)*L);
D = ((lambda/2)*L + (beta*k + 1)*I);

%main loop
while (t<T)
    u_inner = u(2:M,1); % gets the M-1 inner components of U^n
	u_inner = C*(D*u_inner -  beta*k*u_inner.^3);
    u = [0; u_inner; 0];
    t = t + k; % going ahead in time by the specified amount
end