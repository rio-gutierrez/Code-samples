%------------------------------------------------------------------------------
% FEM code to solve the BVP (-ku')' + pu = f  w/ vanishing Dirichlet BCs
%------------------------------------------------------------------------------


%Interval endpoints and number of subintervals
a = 0;
b = 1;
n = 20;

x = linspace(a,b, n+1);    %uniform mesh

%functions to be called
k_funct  = @(x) 1+x;
p_funct  = @(x) 5 .*x .* exp(x);
f_funct   = @(x) -5 .* (x.^3) .* exp(x) + 5 .* exp(x) .*x.^2 + 4 .*x + 1;

M = MassMatD0(x, p_funct);      %call mass matrix
K = StiffMatD0(x, k_funct);     %call stiffness matrix
F = LoadVecD0(x, f_funct);      %call load vector

U = (M+K)\F;           %Solve (M+K)U = F

U_full = [0; U; 0];    %extend solution to include BCs

%Plot results:
plot(x, U_full, "r--x")
hold on
funct = @(x) x - x.^2;     %closed-form solution
fplot(funct, [0,1], "g--o")
ylabel('U(x)')
xlabel('x')
legend("Numerical Solution", "Exact Solution", 'Location','northwest')
exportgraphics(gcf,'FEM_Prob1.pdf')
close
