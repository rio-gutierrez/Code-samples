%Code to solve BVP u''(x) = u + 2/3*e^x    on 0 <= x <= 1
%BCs: u(0) = 0; u(1) = 1/3 *e

%Set constants:
a = 0;
b = 1;
alph = 0;
bet = (1/3)*exp(1);
n = 69;
h = (b-a)/(n+1);

%Generate the matrix A (size nxn):
A = zeros(n);  %initialize nxn matrix
for i = 1:n
    for j = 1:n
        if i == j
            A(i,i) = -(2+h^2)/h^2;
        elseif (j == i+1) || (i == j+1)
            A(i,j) = 1/h^2;
        end
    end
end

%define the x-grid:
x = linspace(a,b,n+2);  %size n+2 (n interior pts + 2 BCs)

%define function vector f:
f = zeros(n,1);     %initialize nx1 vector
for i = 1:n
    if i==n
        f(i) = (2/3)*exp(x(i)) - exp(1)/(3*h^2);
    else
        f(i) = (2/3)*exp(x(i));
    end
end

%Use linear solver to solve Au=f for u:
u  = linsolve(A,f);
usol = [alph; u; bet]; %extend solution to include BCs


%generate vector of errors and vector of exact solution:
abserror_vec =  zeros(n,1);
funct_vec =  zeros(n,1);
funct = @(t) (1/3)*t*exp(t); %closed-form solution
for i = 1:n
    funct_vec(i) = funct (x(i));
    abserror_vec(i) = abs( usol(i) - funct_vec(i) );
end

%extend exact solution vector to include BCs
funct_vec = [alph; funct_vec; bet]; 

%generate vector of errors
abserror_vec =  zeros(n,1);
funct = @(t) (1/3)*t*exp(t); %closed-form solution
for i = 1:n
    abserror_vec(i) = abs( usol(i) - funct(x(i)) );
end

%extend abs error vector to include BCs
abserror_vec = [0; abserror_vec; 0]; 

%Semilog plot of the absolute error vs x:
semilogy(x,abserror_vec, "m+")
ylabel('Error')
xlabel('x')
exportgraphics(gcf,'abserror_Prob4.pdf')
close 


%Plot results:
plot(x,usol, "r*")
hold on
plot(x,funct_vec, "g+-", "LineWidth",2)
ylabel('u(x)')
xlabel('x')
legend("Numerical Solution", "Exact Solution", 'Location','northwest')
exportgraphics(gcf,'BVP_Prob4.pdf')
close 
