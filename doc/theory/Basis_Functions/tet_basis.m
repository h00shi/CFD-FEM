clc
clear all
close all
p = input('Enter the edge Function Order: ') ;

% A quick script to play with polynomals in heirachcal
% sense for triangles, we will use the reference triangle
%[-1,-1],[-1,1],[1,-1]
r = linspace(-1,1,20);
s = r;
u = r;
v = (1-u)/2.0;
w = (1-u)/2.0;
for i = 1:size(r,2)
    for j = 1:size(r,2)
        R(i,j) = r(i);
        S(i,j) = (1.0 + s(j))*(1.0 - r(i))/2.0 - 1.0;
    end
end

% for i = 1:size(u,2)
%     for j = 1:size(u,2)
%         for k = 1:size(u,2)
%             U(i,j,k) = u(i);
%             V(i,j,k) = u(j);
%             W(i,j,k) = u(k);
%         end
%     end
% end
[U,V,W]=meshgrid(u,u,u);
XI1 = (1 + (1 + U).*(1 - W)/2.0 - 1).*(1 - V)/2.0 - 1;
XI2 = (1 + V).*(1 - W)/2.0 - 1;
XI3 = W;

%[XI1,XI2,XI3]=meshgrid(u,v,w);


% FACE POINTS
% Face 1
xif(1:3,1) = [-1,1,-1];
etaf(1:3,1) = [-1,-1,-1];
zetaf(1:3,1) = [-1,-1,1];
% Face 2
xif(1:3,2) = [1,-1,-1];
etaf(1:3,2) = [-1,1,-1];
zetaf(1:3,2) = [-1,-1,1];

% Face 3
xif(1:3,3) = [-1,-1,-1];
etaf(1:3,3) = [1,-1,-1];
zetaf(1:3,3) = [-1,-1,1];

% Face 4
xif(1:3,4) = [-1, -1, 1];
etaf(1:3,4) = [-1, 1, -1];
zetaf(1:3,4) = [-1, -1, -1];

l1 = -(R + S)/2;
l2 = (R + 1)/2;
l3 = (S + 1)/2;
phi=zeros(size(r,2), size(r,2), 4, (p + 1)*(p + 2)*(p + 3)/6);

% Vertex Functions
for i = 1:4
    XI(:,:,i)   = l1*xif(1,i)   + l2*xif(2,i)   + l3*xif(3,i);
    ETA(:,:,i)  = l1*etaf(1,i)  + l2*etaf(2,i)  + l3*etaf(3,i);
    ZETA(:,:,i) = l1*zetaf(1,i) + l2*zetaf(2,i) + l3*zetaf(3,i);
    
    L1(:,:) = -(1+ XI(:,:,i) + ETA(:,:,i) + ZETA(:,:,i))/2;
    L2(:,:) = (XI(:,:,i) + 1)/2;
    L3(:,:) = (ETA(:,:,i) + 1)/2;
    L4(:,:) = (ZETA(:,:,i) + 1.0)/2;
      
    phi(:,:,i,1) = L1;
    phi(:,:,i,2) = L2;
    phi(:,:,i,3) = L3;
    phi(:,:,i,4) = L4;
    
    for ip = 2:p
        nmode = (ip + 1)*(ip + 2)*(ip + 3)/6;
        nmode_pm1 = ip*(ip + 1)*(ip + 2)/6;
        for dof = nmode_pm1+1:nmode
            dtype = dof_type_tet(dof - 1, ip);
            
            switch dtype
                case 1
                    eindex = dof - nmode_pm1;
                    
                    switch eindex
                        case 1
                            phi(:,:,i,dof) = edgepoly(ip, L1, L2);
                        case 2
                            phi(:,:,i,dof) = edgepoly(ip, L2, L3);
                        case 3
                            phi(:,:,i,dof) = edgepoly(ip, L3, L1);
                        case 4
                            phi(:,:,i,dof) = edgepoly(ip, L1, L4);
                        case 5
                            phi(:,:,i,dof) = edgepoly(ip, L2, L4);
                        case 6
                            phi(:,:,i,dof) = edgepoly(ip, L3, L4);
                    end
                case 2
                    li = (dof - 1) - nmode_pm1 - 6;
                    face = floor(li/(ip - 2));
                    fbub = li - face*(ip - 2) + 1;
                                        
                    switch face
                        case 0
                            phi(:,:,i,dof) = tribubpoly(fbub, ip, L1, L2, L4);
                        case 1
                            phi(:,:,i,dof) = tribubpoly(fbub, ip, L2, L3, L4);
                        case 2
                            phi(:,:,i,dof) = tribubpoly(fbub, ip, L3, L1, L4);
                        case 3
                            phi(:,:,i,dof) = tribubpoly(fbub, ip, L1, L3, L2);
                    end
                case 3
                    ibub = dof - 1 - nmode_pm1 - 6 - (ip - 2)*4  + 1;
                    phi(:,:,i,dof) = tetbubpoly(ibub, ip, L1, L2, L3, L4);
                    
                          
            end
            
        end
    end
    
end   
clear L1;
clear L2;
clear L3;
clear L4;
L1(:,:,:) = -(1.0 + XI1 + XI2 + XI3)/2.0;
L2(:,:,:) = (XI1 + 1.0)/2.0;
L3(:,:,:) = (XI2 + 1.0)/2.0;
L4(:,:,:) = (XI3 + 1.0)/2.0;
phiV=zeros(20,20,20,max(1,(p - 3)*(p - 2)*(p - 1)/6));
for ip=4:p
    nbub = (ip - 3)*(ip - 2)*(ip - 1)/6;
    nbubpm1 = (ip - 1 - 3)*(ip - 1 - 2)*(ip - 1 - 1)/6;
    for dof = nbubpm1 + 1:nbub;
        phiV(:,:,:,dof) = tetbubpoly(dof, ip, L1, L2, L3, L4);
        %tetbubpoly(dof, ip, L1, L2, L3, L4);
    end
end


for i = 1:(p + 1)*(p + 2)*(p + 3)/6
    figure(i);
    hold on
    surf(XI(:,:,1),ETA(:,:,1), ZETA(:,:,1),phi(:,:,1,i),'FaceColor','Interp','FaceLighting','phong');
    surf(XI(:,:,2),ETA(:,:,2), ZETA(:,:,2),phi(:,:,2,i),'FaceColor','Interp','FaceLighting','phong');
    surf(XI(:,:,3),ETA(:,:,3), ZETA(:,:,3),phi(:,:,3,i),'FaceColor','Interp','FaceLighting','phong');
    surf(XI(:,:,4),ETA(:,:,4), ZETA(:,:,4),phi(:,:,4,i),'FaceColor','Interp','FaceLighting','phong');
    xlabel('\xi','Fontsize',14);
    ylabel('\eta','Fontsize',14);
    zlabel('\zeta','Fontsize',14);
    view(-30,30);
    colorbar;
    hold off
end
 XI_const_xi=zeros(size(u,2));
[XI_const_eta, XI_const_zeta]=meshgrid(u,u);



figure;
slice(U,V,W, phiV(:,:,:,1),0,0,0);
xlabel('\xi','Fontsize',14);
ylabel('\eta','Fontsize',14);
zlabel('\zeta','Fontsize',14);
colorbar;
   