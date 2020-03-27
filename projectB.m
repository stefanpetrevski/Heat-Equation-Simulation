%% Project B Numerical solution of the heat equation

%% Exercise 4.1

%%
% See comments next to each line in the code of the following function for
% intepretation/explanation of what specific commands do.

type heat_eqn.m

%% Part a)

heat_eqn(1,20,0.3)

%% Part b)

heat_eqn(1,20,1)

%%
% For $$ \lambda = 0.3 $$, the solution is bounded and stable: the heat is 
% at first mostly concentrated just to the left of $$ x = 0.4 $$. Then, as
% time progresses, the heat is being distributed more and more evenly 
% across the interval $$ -1 \leq x \leq 1 $$, and ultimately it is 
% in the shape reminiscent of a low-arcing, symmetric 
% parabola with the peak heat having shifted to where $$ x = 0 $$. The
% physical intepretation of this problem is that it is a Dirichlet
% problem with a zero steady-state situation, whereby the ends are not
% thermally insulated and the heat may escape the object, which is
% what happens as we go further in time. From the Fourier Series and PDEs
% course, we would expect the coefficients for the distribution of heat 
% with time to decay exponentially, solidifying this claim further.
% Meanwhile, for $$ \lambda = 1 $$, our solution initially 
% behaves almost identically to the previous one until just after
% $$ t = 0.35 $$ Indeed, upon reaching a shape similar to that previously
% obtained in a), it then becomes unstable, shooting off to $$ \pm \infty $$  
% at some points and remaining $$ 0 $$ at others. It continues to exhibit 
% this extreme behaviour by alternating the positions of these vertical
% lines all the way until $$ t = 1 $$.

%% Exercise 4.2

%%
% We use the stable_test function to create a $$ 1 \times 10 $$ vector whose
% entries are a logical true (i.e. 1) if a stable solution occurs, or a 
% logical false (i.e. 0) otherwise. The $$ (1,k) $$-th entry corresponds
% to the type of solution for $$ \lambda = 0.k $$, for the first 9 
% positive integers $$ k $$ and the last entry corresponds to 
% $$ \lambda = 1 $$.

type stable_test.m
lambdaresults = [];
for lambda = 0.1 : 0.1 : 1
    lambdaresults = [lambdaresults stable_test(lambda)];
end
lambdaresults

%%
% Hence, if we run the stable_test.m code for $$ M = 20 $$, we see that 
% our solutions will be stable as long as $$ \lambda $$ is strictly
% smaller than $$ 0.6 $$, i.e. it enables a stable solution for the first
% 5 values of $$ \lambda $$ and an unstable one from this point onwards.

%% Exercise 4.3

type heat_eqn1.m

%% Exercise 4.4

[u u_exact] = heat_eqn1(1,20,0.5);
error = norm(u - u_exact, inf)

%% Exercise 4.5

%%
% In this case, the heat_eqn1 code has been altered notably to
% obtain heat_eqn1error: instead of outputting the exact and the 
% numerical approximation of the solution, we now directly get the value
% of the infinite norm (our convention for the error). Furthermore,
% the code is vectorized, and operations including logarithming and
% applying the heat_eqn1error function are carried out in a 
% component-wise manner as shown below. 

type heat_eqn1error.m


mvector = []; % starting off with an empty vector
for M = 20:1:100 % looping over desired values of M
    mvector = [mvector 2/M]; % generating [2/20 2/21 2/22 ... 2/100]
end        

x = log(mvector); % x = ln(epsilon) as defined in the question
y = log(arrayfun(@(M) heat_eqn1error(1, M, 0.5), linspace(20,100,81)));
% applying the function component-wise for M = {20, 21, ..., 100}
p = polyfit(x, y, 1) % linear regression for the data above
f = polyval(p,x); % returns a linear function based on polyfit
plot(x,y,'*',x,f,'-') 
xlabel('ln(h)'); % labeling x-axis
ylabel('ln(epsilon)'); % labeling y-axis
legend('Actual data','Linear fit','Location','southwest')
title('Linear fit: ln(epsilon) vs. ln(h) for M = 20,21,...,100')

%%
% From the linear regression model, we can see that the slope is 
% approximately $$ 1.9986 $$, hence after getting rid of the logarithms,
% we obtain approx. that $$ \varepsilon \propto h^2 $$, with constant
% of proportionality equal to around $$ e^{-2.4580} $$. Therefore 
% $$ \varepsilon $$ and $$ h $$ satisfy approximately a positive quadratic
% relationship. This makes sense intuitively, as when we increase $$ M $$
% (i.e. generate smaller steps for the $$ x $$ values) then $$ h $$
% decreases. Hence, when using the finite difference approximations,
% we use points increasingly closer to the ones we were given before,
% thus (taking into account that our function is continuous and behaves
% well) the calculated estimates for the derivatives are better, so we 
% should expect to see smaller errors.
 
%% Exercise 4.6
 
type heat_eqn2.m

%% Exercise 4.7

%% Part a)

heat_eqn2(1,20,0.3)

%% Part b)

heat_eqn2(1,20,1)

%%
% For $$ \lambda = 0.3 $$, the solution is again stable: the heat is being 
% distributed more and more evenly across the interval, and ultimately 
% it is in the shape reminiscent of a low-arcing, symmetric 
% parabola with the peak heat having shifted to where $$ x = 0 $$. In other
% words, it behaves almost identically to the numerical approximation 
% initially for the same value of $$ \lambda $$. Furthermore, for 
% $$ \lambda = 1 $$, the difference from initially is notable: the odd 
% vertical line behaviour no longer exists, and instead we have an 
% analogous behaviour to that for $$ \lambda = 0.3 $$, i.e. our solution 
% still is bounded and stable. This is consistent with the observation in
% the project booklet that solutions will stay bounded when using the more
% elaborate method regardless of the value of $$ \lambda $$.

%% Exercise 4.8

%%
% Generalizing the component-wise expressions to obtain a matrix form 
% for the equations, we obtain the following:
%
% $$(I- \frac{\lambda}{2}L)U^{n+1}=((\beta k+1)I+
% \frac{\lambda}{2}L)U^{n}-\beta k(U^{n})^3 $$
% 
% where the last vector is cubed component-wise.
% Hence, we can set:
% 
% $$ C = (I- \frac{\lambda}{2}L)^{-1} $$
%
% $$ D = (\beta k+1)I+\frac{\lambda}{2}L $$.
%
% Our problem then becomes:
%
% $$ U^{n+1} = C(DU^n-\beta k(U^{n})^3) $$
%
% where again we have component-wise cubing. This is reflected in the 
% code of the allen.cahn function below, where this step is repeated in
% the iterative while loop.

type allen_cahn.m

%% Exercise 4.9

%%
% For aesthetic and spatial purposes, we again place the vector solutions
% in a matrix (as we did previously for $$ \lambda $$). The $$ k $$-th
% column of the matrix represents the solution for $$ \beta = k $$, where
% $$ 1 \leq k \leq 5 $$. Generate also a solution matrix for $$ t = 2000 $$
% using the allen_cahn_infinite code, which only suppresses the animation.
% The first matrix is called acmatrix, whereas the latter is infmatrix.

acmatrix = [];
for beta  = 1:5
   acmatrix = [acmatrix allen_cahn(1,20,0.5,beta)];
   % placing the solutions in a 21x5 matrix
end    
acmatrix

infmatrix = [];
for beta  = 1:5
   infmatrix = [infmatrix allen_cahn_infinite(2000,20,0.5,beta)];
   % placing the solutions in a 21x5 matrix
end    
infmatrix

%%
% The Allen-Cahn equation provides us with a trend that is quite a bit
% different from the ones previously observed, whereby the heat tended 
% to zero symmetrically. Consider as previously the case when 
% $$ \lambda = 0.5 $$. In the animations obtained now, the heat initially
% behaves analogously, as it decreases, and tends to a parabolic-like shape
% symmetric about the line $$ x = 0 $$. For the initial values of
% $$ \beta = 1, 2 $$, the final solutions seems to be tending very 
% closely to $$ 0 $$, although in the latter case this convergence is 
% somewhat slower. From $$ \beta = 3 $$ and onwards, the solutions tend to
% a parabolic-likerelationship that is now quite taller than in the
% previous cases, and are clearly distinct from the zero solutions.
% The peaks of these "parabolae" are greater as $$ \beta $$ increases.
% Upon a conversation with my tutor, I was advised to conceptually think 
% of the Allen-Cahn equation as being governed by a (chemical) reaction
% in addition to the flow of the heat. After heat flows away from the peak 
% at the right of $$ x = 0 $$, the value of $$ u(x,t) $$ drops below
% $$ 1 $$ everywhere. As soon as $$ u < 1 $$, $$ \beta(u-u^3) > 0 $$,
% as the cubic term has less of an effect than the linear one. Note that
% further in time, as the solution "stabilizes" around a certain value 
% of $$ u(x,t) $$, then $$ u_{t} $$ approximately vanishes. Hence, the total 
% contribution of the right-hand side must be approximately zero, hence
% $$ u_{xx} $$ is of about the same magnitude as $$ \beta(u-u^3) $$ but
% negative. This is consistent with the observation that as $$ \beta $$
% increases so does the net effect of the reaction, and there is an 
% increase in the magnitude of $$ u_{xx} $$ so that each "parabola" is 
% steeper (and taller) than the one for the previous smaller value of 
% $$ \beta $$. For completeness, I have included an allen_cahn_infinite.m
% code which takes the same inputs as the allen_cahn function, yet skips 
% the animation and only outputs the solution vector. I repeatedly used 
% this in the command window to consolidate the behaviour as t becomes very
% large, say $$ 10000 $$ or $$ 20000 $$, so I concluded that the negligible
% change at these great times should support my conjecture.
