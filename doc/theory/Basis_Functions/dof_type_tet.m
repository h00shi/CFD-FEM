function dtype = dof_type_tet(k, p )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    ndofpm1 = (p - 1 + 1)*(p - 1 + 2)*(p - 1 + 3)/6;
    
    li = k - ndofpm1;
    
    ne = 6; 
    nf = (p - 2)*4; 
    nb = (p - 2)*(p - 3)/2; 
   
    if( li - ne < 0 )
      dtype = 1;
    elseif( li - ne - nf < 0) 
      dtype = 2;
    elseif( li - ne - nf - nb < 0 ) 
      dtype = 3; 
    else 
      dtype = -1;
    end     
    
end

