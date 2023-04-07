clc;
clear;
disp("Code by Patil Aniruddha Ramesh")
disp("Roll No. : 200104072")
disp("For Trigonometric Basis : ")
disp("Basis are : ")
disp("Trigonometric Matrix : ")

disp(" ")

    syms x ;


    %p1=sin(pi*x);
    %p2=sin(2*pi*x);
    %p3=sin(3*pi*x);
    %p4=sin(4*pi*x);
    %p5=sin(5*pi*x);

    p1=x*(1-x) ;
    p2=x*(1/2-x)*(1-x) ;
    p3=x*(1/3-x)*(2/3-x)*(1-x) ;
    p4=x*(1/4-x)*(1/2-x)*(3/4-x)*(1-x) ;
    p5=x*(1/5-x)*(2/5-x)*(3/5-x)*(4/5-x)*(1-x) ;


    pMat=[p1 p2 p3 p4 p5];
    
    disp(pMat)
    disp(" ")


for n_basis=1:5 

    formatStr = "For Number of Basis as %d\n";
    fprintf(formatStr,n_basis);

    kMatrix=1+x;
    K=zeros(n_basis,n_basis);




    for var=1:n_basis
        for let=1:n_basis
            K(var,let)= eip(pMat(var),pMat(let),kMatrix);
        end
    end

    disp('K=');
    disp(K);

    f = x^2;
    F=zeros(n_basis,1);

    for i=1:n_basis
        F(i)=ip(pMat(i),f);
    end
    disp('F=') ;
    disp (F)

    c=K\F;
    disp('c=') ;
    disp(c) ;
    
    
    syms P
    
    P =0;
    
    for g=1:n_basis
        P = P + c(g)*pMat(g);
    end
        


    %P=c(1)*p1+c(2)*p2+c(3)*p3+c(4)*p4+p5*c(5);
    %disp(P);


    %syms P(x)
    %P=0;
    %for loopj=1:5
    %   P=P+c(loopj)*pMat(loopj);
    %end

    syms y(x)
    ode=(1+x)*diff(y,x,2)+diff(y,x)+x^2==0;
    syms resi
    residual = simplify((1+x)*diff(P,x,2)+diff(P,x)+x^2);
    boundcond1=y(0)==0;
    boundcond2=y(1)==0;
    boundcond=[boundcond1 boundcond2];

    ySol(x)=dsolve(ode,boundcond);
    
    fprintf("The Exact Solution is : ySol(x)=");
    
    disp(ySol);
    
    fprintf("The Approximate Solution is : P(x)=");
    
    disp(P);
    
    
    %syms lhs
    %syms rhs
    %syms residual
    
    
   
    %lhs = simplify(subs(ode, y(x), P)); % Substitute P(x) for y(x) in the left-hand side of the equation
    %rhs = simplify(subs(lhs, diff(P, x), diff(y(x), x))); % Substitute the first derivative of P(x) for the first derivative of y(x) in the right-hand side of the equation
    %residual = simplify(lhs - rhs); % Compute the residual as the difference between the left-hand side and the right-hand side


    %difference=ySol-P;
    
    fprintf("The residual is : R(x)=");
    
    disp(residual);
    
    fprintf("The Error is : Error =");
    
    disp(simplify(ySol -P));

    divs=100; % Number of divisions for domain for plotting
    X = linspace(0,1,divs); % Define the domain
    
    figure(n_basis)

    fplot(P,'g',[0 1], "Linewidth",3)
    hold on
    fplot(ySol,'r',[0 1], 'LineStyle', '--', "Linewidth",1.5)
    hold on
    fplot(residual,'b',[0 1], "Linewidth",1.5)
    ylabel("u")
    xlabel("x")
    form = "Plot of u for number of basis ";
    title(form,n_basis,'color','black')
    legend("u_{actual}","u_{trial}","Residual")


end
function I= eip (f, g, k)
syms x  ;
I= int (k*diff(f, x)*diff(g, x),x ,0, 1);
end 

function I= ip(g,f)
syms x  ;
I= int (f*g,x,0,1);
end