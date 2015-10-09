function edgep = edgepoly(n, L1, L2)
    edgep = L1.*L2.*kernel(n - 2 ,L2 - L1);
end