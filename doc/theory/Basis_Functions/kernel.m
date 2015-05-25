function psi = kernel(km2, x)

% Get the order 
k = km2 + 2;
sigma = 1/sqrt(2/(2*k - 1));

psi = -2*sigma/(k - 1).*Jacobi(km2, 1, 1, x);
end