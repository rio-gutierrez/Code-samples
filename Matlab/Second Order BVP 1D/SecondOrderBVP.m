%Code to solve BVP u_{xx}(x) = 1 - |x| on -1 < x < 1
%BCs: u(-1) = 5; u(1) = 7

%Set constants:
a = -1;
b = 1;
alph = 5;
bet = 7;
n = 49;
h = (b-a)/(n+1);

%Generate the matrix A (size nxn):
A = zeros(n);  %initialize nxn matrix
for i = 1:n
    for j = 1:n
        if i == j
            A(i,i) = -2/h^2;
        elseif (j == i+1) || (i == j+1)
            A(i,j) = 1/h^2;
        end
    end
end

%define the x grid in either of the two ways:

%Method 1:
% for i = 0:n+1
%     x(i+1) = -1+i*h;
% end

%or Method 2:
x = linspace(a,b,n+2);  %size n+2 (n interior pts + 2 BCs)

%define function vector f:
f = zeros(n,1);     %initialize nx1 vector
for i = 1:n
    if i==1
        f(i) = 1 - abs(x(i+1)) - alph/h^2;
    elseif i == n
        f(i) = 1 - abs(x(i+1)) - bet/h^2;
    else
        f(i) = 1 - abs(x(i+1));
    end
end

%Use linear solver to solve Au=f for u:
u  = linsolve(A,f);
usol = [alph, u', bet]; %extend solution to include BCs

%Plot results:
plot(x,usol', "r--x")
hold on
funct = @(t) -(abs(t)^3)/6 + t^2/2 + t + 17/3; %closed-form solution
fplot(funct, [-1,1], "g--o")
ylabel('u(x)')
xlabel('x')
legend("Numerical Solution", "Exact Solution", 'Location','northwest')
exportgraphics(gcf,'BVP_1.pdf')
close 
