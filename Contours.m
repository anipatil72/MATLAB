close all
clear 
clc
 
% Parameters provided
L = 1;
E = 200*10^9;
a = 0.8;
b = 0.05;
h = 0.155;
I = (1/12)*(b*h^3);
f = 25000;
M1 = f*a;
M2 = f*(L-a);
 
x = 0:0.01:0.8;
y = 0:0.01:0.8;
 
delta_1 = - (f*x.^2)./(6.*E.*I).*(3.*a - x);  % for 0 <= x <= a
delta_2 = - (f*a.^2)./(6.*E.*I).*(3.*x - a);  % for a <= x <= L
 
sigma_1 = (M1.*delta_1)./I;
 
x = rand(100,1)*4-2;
y = rand(100,1)*4-2;
z = (M1.*((f*x.^2)./(6.*E.*I).*(3.*a - x)))./I;
z1 = (M2.*((f*a.^2)./(6.*E.*I).*(3.*x - a)))./I;
 
F = TriScatteredInterp(x,y,z);
ti = -1:.25:1;
[qx,qy] = meshgrid(ti,ti);
qz = F(qx,qy);
 
figure(1)
scatter3(x,y,z)
title('Mesh Plot for 0 <= x <= a')
hold on
mesh(qx,qy,qz)
 
figure(2)
contour(qx,qy,qz)
title('Contour Plot for 0 <= x <= a')
 
F = TriScatteredInterp(x,y,z1);
ti = -1:0.25:1;
[qx,qy] = meshgrid(ti,ti);
qz = F(qx,qy);
 
figure(3)
scatter3(x,y,z)
title('Mesh Plot for a <= x <= L')
hold on
mesh(qx,qy,qz)
 
figure(4)
contour(qx,qy,qz)
title('Contour Plot for a <= x <= L')