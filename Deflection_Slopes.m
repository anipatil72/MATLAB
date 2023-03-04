close all
clear 
clc

%Program for Calculating Deflection and slopes throughout the beam from
%point loading only
 
% Parameters from the question
E = 200*10^9;
a = 0.8;
b = 0.05;
h = 0.155;
I = (1/12)*(b*h^3);
f = 25000;
 
x = 0:0.01:0.8;
delta_1 = - (f*x.^2)./(6.*E.*I).*(3.*a - x);  % for 0 <= x <= a
delta_2 = - (f*a.^2)./(6.*E.*I).*(3.*x - a);  % for a <= x <= L
 
figure(1)
plot(x, delta_1)
xlabel('Length (m)')
ylabel('Beam Deflection (m)')
title('Cantilever Beam Deflection for 0 <= x <= a')
 
figure(2)
plot(x, delta_2)
xlabel('Length (m)')
ylabel('Beam Deflection (m)')
title('Cantilever Beam Deflection for a <= x <= L')