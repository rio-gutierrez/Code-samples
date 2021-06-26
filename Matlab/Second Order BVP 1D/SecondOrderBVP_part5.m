%Code to solve BVP u''(x) = (2 + 4x^2)u = 0  on 0 <= x <= 1
%BCs: u(0) = 1; u(1) = e

%Set constants:
a = 0;
b = 1;
alph = 1;
bet = exp(1);
n = 69;
h = (b-a)/(n+1);

%define the x-grid
%size n (n interior pts only for now; 2 bd pts will be added later):
x = linspace(a+h,b-h,n);  

% %define the Psi vector
 Psi = zeros(n,1);     %initialize nx1 vector
 for i = 1:n
     Psi(i) = h^2 * (2 + 4*x(i)^2);
 end
 
 

%Generate the matrix A (size nxn):
A = zeros(n);  %initialize nxn matrix
for i = 1:n
    for j = 1:n
        if i == j
            A(i,i) = - ( 2 + Psi(i) );
        elseif (j == i+1) || (i == j+1)
            A(i,j) = 1;
        end
    end
end


%define function vector f:
f = zeros(n,1);     %initialize nx1 vector
for i = 1:n
    if i==1
        f(i) = - 1;
    elseif i == n
        f(i) = - exp(1);
    end
end

 
%Use linear solver to solve Au=f for u:
u  = linsolve(A,f);
usol = [alph; u; bet]; %extend solution to include BCs


%redefine x to include bd pts
x = linspace(a,b,n+2);  %size n+2 (n interior pts + 2 BCs)


%generate vector of errors and vector of exact solution:
abserror_vec =  zeros(n,1);
funct_vec =  zeros(n,1);
funct = @(t) exp(t^2); %closed-form solution
for i = 1:n
    funct_vec(i) = funct (x(i));
    abserror_vec(i) = abs( usol(i) - funct_vec(i) );
end

%extend abs error vector to include BCs
abserror_vec = [0; abserror_vec; 0]; 

%Semilog plot of the absolute error vs x:
semilogy(x,abserror_vec, "m+")
ylabel('Error')
xlabel('x')
exportgraphics(gcf,'abserror_Prob5.pdf')
close 


%extend exact solution vector to include BCs
funct_vec = [alph; funct_vec; bet]; 

%Plot results:
plot(x,usol, "b*")
hold on
plot(x,funct_vec, "m-o", "LineWidth",2)
ylabel('u(x)')
xlabel('x')
legend("Numerical Solution", "Exact Solution", 'Location','northwest')
exportgraphics(gcf,'BVP_Prob5.pdf')
close 