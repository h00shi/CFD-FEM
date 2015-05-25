function bubp = tetbubpoly(ibub, p, L1, L2, L3, L4 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    np = p - 4 + 1;
    k = 0;
    n1 = 0;
    n2 = 0;
    n3 = 0;
    
    %---> Based on value of ibub get the indicies n1, n2, n3;  
    %---> Ok this is a bit strange...but the index triples are layed out 
    %  in a triangular pattern of "points".  This triangular patter is such 
    %  that the max value of any n1,n2,n3 is np with is p - 3.  We can actually 
    %  represent these values exactly using linear affine(barycentric 
    %  coodrinates).  We can also rep phyiscal coorindates by i,j and using 
    % all this gives us n1, n2, n3.  */
    
    k = 0;
    for j = 0:np
        for i = 0:np - j
            k =k + 1;
            if( k == ibub)
                n1 = np - i - j;
                n2 = 1 + i;
                n3 = 1 + j;
            end
        end
    end
    bubp = L1.*L2.*L3.*L4.*...
        kernel(n1 - 1, L3 - L1).*...
        kernel(n2 - 1, L2 - L1).*...
        kernel(n3 - 1, L4 - L1);

end

