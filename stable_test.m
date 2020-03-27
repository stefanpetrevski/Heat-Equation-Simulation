function result = stable_test(lambda)
M = 20; % taking M = 20 as previously in the problem
e = ones(M-1,1);
L = spdiags([e  -2*e  e], [-1 0 1], M-1, M-1);
I = speye(M-1);
A = I + lambda*L;
eigsvector = eigs(A); % returns the 6 largest abs of eigenvalues
        if max(abs(eigsvector)) <= 1
            result = true; % when every abs(eigenvalue)<=1
        else
            result = false; % if there exists an abs(eigenvalue)>1
        end
end
        
