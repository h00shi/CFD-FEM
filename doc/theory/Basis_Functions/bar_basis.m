clc
clear all
close all
p = input('Enter the edge Function Order: ') ;

% A quick script to play with polynomals in heirachcal
% sense for triangles, we will use the reference triangle
%[-1,-1],[-1,1],[1,-1]
xi = linspace(-1,1,50);
X = xi;

% We want to make and plot the vertex edge and bubble functions of a
% certain order

%Vertex Functions
L1 = (1 - X)/2.0;
L2 = (X + 1)/2.0;

phi(:,1) = L1;
phi(:,2) = L2;

for i = 2:p
    phi(:,i+1) = edgepoly(i, L1, L2);
end

for i = 1:p+1
    figure(i)
    plot(X,phi(:,i),'Linewidth',2.0);
    xlabel('$\xi$','interpreter','latex','Fontsize',14)
    ylabel('$\phi(\xi)$','interpreter','latex','Fontsize',14)
end
