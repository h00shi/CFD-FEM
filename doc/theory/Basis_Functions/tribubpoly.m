function bubp = tribubpoly(ibub, p, L1, L2, L3)
n1 = ibub;
n2 = (p - 1) - n1;
bubp = L1.*L2.*L3.*kernel(n1 - 1 ,L2 - L1).*...
    kernel(n2 - 1, L1 - L3);
end