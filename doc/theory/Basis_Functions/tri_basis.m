clc
clear all
close all
p = input('Enter the edge Function Order: ') ;

% A quick script to play with polynomals in heirachcal
% sense for triangles, we will use the reference triangle
%[-1,-1],[-1,1],[1,-1]

xi = linspace(-1,1,20);
eta = xi;
for i = 1:size(xi,2)
    for j = 1:size(xi,2)
        XI(i,j) = xi(i);
        ETA(i,j) = (1 + eta(j))*(1 - xi(i))/2-1;
    end
end

% We want to make and plot the vertex edge and bubble functions of a
% certain order

%Vertex Functions
L1 = -(XI + ETA)/2;
L2 = (XI + 1)/2;
L3 = (ETA + 1)/2;

phi(:,:,1) = L1;
phi(:,:,2) = L2;
phi(:,:,3) = L3;

for ip = 2:p
    nmode = (ip + 1)*(ip + 2)/2;
    nmode_pm1 = ip*(ip + 1)/2;
    for i = nmode_pm1+1:nmode
        li = i - nmode_pm1 - 1;
       
       
        if( li - 3 < 0)
            dtype = 1
        elseif( (li - 3) - (ip - 2) < 0 )
            dtype = 2
        else
            dtype = -1;
        end
        
        if( dtype == 1)
            % Edge Functions
            if(li == 0)
                phi(:,:,i) = edgepoly(p, L1, L2);
            elseif(li == 1)
                phi(:,:,i) = edgepoly(p, L2, L3);
            elseif(li == 2)
                phi(:,:,i) = edgepoly(p, L3, L1);
            end
            % Bubble Functions
            elseif(dtype == 2)
                ibub = (i - nmode_pm1 - 1) - 3 + 1;
                phi(:,:,i) = tribubpoly(ibub, ip, L1, L2, L3);  
            end
        end
    end

nmode = (p + 1)*(p + 2)/2;
for i = 1:nmode
    figure(i)
    surf(XI, ETA, phi(:,:,i),'FaceColor','interp' )
    xlabel('\xi','Fontsize',14)
    ylabel('\eta','Fontsize',14);
end
