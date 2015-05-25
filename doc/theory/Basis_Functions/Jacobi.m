function P = jacobi(N,alpha,beta,x)


%P0 = ones(size(x,1),size(x,2));
P0=ones(size(x));
P1 = 1/2*(alpha - beta + (alpha + beta + 2).*x);
if( N == 0) 
    P = P0;
elseif( N == 1) 
    P = P1;
else
   for i = 2:N
       n = i - 1;
       an1 = 2*(n + 1)*(n + alpha + beta + 1)*(2*n + alpha + beta);
       an2 = (2*n + alpha + beta + 1)*(alpha^2 - beta^2);
       an3 = (2*n + alpha + beta)*(2*n + alpha + beta + 1)*(2*n + alpha + beta + 2);
       an4 = 2*(n + alpha)*(n + beta)*(2*n + alpha + beta + 2);
       
       P = 1./an1*((an2.*ones(size(x,1),size(x,2)) + an3.*x).*P1 - an4.*P0);
       P0 = P1;
       P1 = P;
   end 
end