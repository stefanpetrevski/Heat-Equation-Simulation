function u = heat_eqn2(T,M,lambda)

%set up grid points
h = 2/M; % grid spacing for splitting the x-axis
x = (linspace(-1,1,M+1))'; % splitting the interval [-1,1]
k = lambda*h^2; % grid spacing for going ahead in time
N = ceil(T/k);

%initial conditions
u = cos(pi*x/2) + (sin(pi*x))/2; % u(x,0)
t = 0;

%plot u(x,0)
plot(x,u,'b','LineWidth',2); % plotting the function u(x,0)
axis([-1 1 0 1.3]); % setting the min/max values of the axes displayed
title('t = 0'); % adding a title that displays the initial time
xlabel('x'); % labeling x-axis
ylabel('u'); % labeling y-axis

%finite difference operator
e = ones(M-1,1); % generating an (M-1)x1 column vector of 1s only
L = spdiags([e  -2*e  e], [-1 0 1], M-1, M-1); 
% tridiagonal and sparse, with main diag of -2s, adjacent of -1s
I = speye(M-1); % sparse identity matrix
B = inv(I - (lambda/2)*L)*(I + (lambda/2)*L); % B as defined 

%main loop
while (t<T)
    u_inner = u(2:M,1); % gets the M-1 inner components of U^n
	u_inner = B*u_inner; % calculating solution forward in time
    u = [0; u_inner; 0];
    %plot current solution
    hold off
    plot(x,u,'LineWidth',2); % plotting u(x,t)
    axis([-1 1 0 1.3]); % setting the min/max of the axes
	title(['t = ' num2str(t)]); % converts t to character to display
    xlabel('x'); % labeling x-axis
    ylabel('u(x,t)'); % labeling y-axis
	pause(0.02); % pausing to display solution sufficiently long
    t = t + k; % going ahead in time by the specified amount
end