%Ekansh Roy
%190104030


syms d(x)
syms g(x)
syms a(x)

%------------------------------------------------------------------------------

syms x
Funct_Shape_1=x*(1-x);
Funct_Shape_2=(x^2)*(1/2-x)*((1-x)^2);
Funct_Shape_3=(x^3)*(1/3-x)*(2/3-x)*((1-x)^3);
Funct_Shape_4=(x^4)*(1/4-x)*(1/2-x)*(3/4-x)*((1-x)^4);
kMatrix=1+x;
f=x^2;

%------------------------------------------------------------------------------

pMat=[Funct_Shape_1 Funct_Shape_2 Funct_Shape_3 Funct_Shape_4 ];

K=zeros(4,4);

for loopi=1:4
    for loopj=1:4
        K(loopi,loopj)=eip(pMat(loopi),pMat(loopj),kMatrix);
    end
end

%------------------------------------------------------------------------------
%------------------------------------------------------------------------------


syms y(x)

ode=(1+x)*diff(y,x,2)+diff(y,x)-x^2==0;

boundcond1=y(0)==0;
boundcond2=y(1)==0;
boundcond=[boundcond1 boundcond2];

ySol(x)=dsolve(ode,boundcond);


%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

F=zeros(4,1);

for loopi=1:4
    F(loopi)=fip(pMat(loopi),f);
end

c=K\F;
syms P(x)
P=0;

difference=ySol-P;


for loopj=1:4
    P=P+c(loopj)*pMat(loopj);
end


fplot(P,'r',[0 1])
hold on
fplot(difference,'b',[0 1])
fplot(ySol,'g', [0 1])


function I=eip(d,g,a)
    d=diff(d);
    g=diff(g);
    I=int(d*g*a,[0 1]);
end

function Zeka=fip(d,g)
    Zeka=int(d*g,[0 1]);
end






