close all
clear
clc

%Program for Calculating Shear Force and Bending Moments throughout the beam from
%point loading only

% For simplicity length is taken as unity 
 
L = 1; % Unit in meters
n = 2; % Number of point loads on the beam
 
for i = 1:n
fprintf('Enter the point load and distance where load acts for %d node\n',i)
Wc(i)=input('Enter the load in Newton\n');
Lc(i) = input('Enter the distance of the point load from the fixed end\n');
            
if Lc(i)>L || Lc(i)<0                
    error('Please check the Length')            
end
end
 
NL = zeros(1,n);        
NW = zeros(1,n);
       
for i=1:n
    [a,b]=max(Lc);            
    NL(i)=a;           
    NW(i)=Wc(b);          
    Lc(b)=[];           
    Wc(b)=[];        
end
 
NL(n+1) = 0;
        
% Shear force diagram
figure(1)
Ra = sum(NW);
if NL(1)==L
X1 = L;
F1 = NW(1);
elseif i==1
X1 = [L NL(1) NL(1) 0];
F1 = [0 0 Ra Ra];
else
X1 = [L NL(1)];
F1 = [0 0];
end
S=[]; X=[]; F=[];
        if n>=2
            for i=2:n+1
                x = [NL(i-1) NL(i)];
                S= [S NW(i-1)];
                f = sum(S);
                X= [X x];
                F = [F f f];
            end
        end
subplot(2,1,1)
        
plot([X1 X],[F1 F],'b')
xlim([0 L]);
hline = refline(0,0);
hline.Color = 'k';
legend('Shear Force','Reference')          
title('Shear Force Diagram of Cantilever Beam')
xlabel('Length of the Beam in Meter')
ylabel('Shear Force in Newton')
  
hold on
grid on
fprintf('Reaction Force =%d N\n',Ra)
    
% Bending moment diagram
subplot(2,1,2)
if NL(1)<L
X1=NL(1): L;
M1 = zeros(size(X1));
plot(X1,M1,'r')
hold on 
grid on
end
 
Xm=[]; 
Mm=[];
    
for i=1:n
X = NL(i+1):0.1:NL(i);
M = -(NW(1,1:i)*((NL(1,1:i))'-newX(X,i))) ;
Xm =[X Xm];
Mm =[M Mm];
end    
    
plot(Xm,Mm,'r') 
xlim([0 L]);
hline = refline(0,0);
hline.Color = 'k';
legend('Bending Moment','Reference')  
    
title('Bending Moment Diagram of Cantilever Beam')
xlabel('Length of the Beam in Meter')
ylabel('Bending Moment in Newton-Meter')
hold off
grid on
 
function x = newX(X,i)                                                
[~,d] = size(X);                                                 
x = zeros(i,d);                                                
for j=1:i                                                        
x(j,1:d) = X;                                                   
end                                                             
end    