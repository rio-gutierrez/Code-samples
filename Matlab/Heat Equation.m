%Backward Euler solution to the Heat Eq using only three hat functions
%(result expected to look rough!)

m = 3;
h = 1/4;
dt = h^2/2;    %(can be larger since we are using an implicit scheme)

%Generate arrays:
M = zeros(m);
K = zeros(m);
U_0 = ones(m,1);    %initial conditions
rhs = zeros(m,1);   %initialize rhs of Eq


for i = 1:m
    M(i,i) = 1/6;
    K(i,i) = 8;
    if i ~= m
        M(i,i+1) = 1/24;
        M(i+1,i) = M(i,i+1);
        K(i,i+1) = -4;
        K(i+1,i) = K(i,i+1);
    end
end

Mat = M + dt*K;

f_vec = [17/384; 23/384; 17/384];
f = @(t) cos(t);



%------------------------------------------
%        BACKWARD EULER CODE
%------------------------------------------

it_max = 500;      %max number of iterations allowed
tol = 1e-5;        %tolerance allowed
it = 0;

for n = 1 : it_max
    it = it +1;
    rhs = M* U_0 + dt * f((n+1)*dt) * f_vec;
    U = Mat\rhs;

    if norm(U - U_0) <= tol
        disp(['It took ', num2str(it),
              ' iterations for the solution to converge.'])
        break
    elseif  it == it_max
      disp('No convergence; max number of iterations reached.')
    end

    disp(num2str( norm(U - U_0) ));

    U_0 = U;  %update U_0 value for next iteration
end

%----------------------------------------------
%      END OF BACKWARD EULER CODE
%----------------------------------------------


%extend solution to include boundaries
U = [0; U; 0];

x = linspace(0,1,m+2);

%Plot results:
plot(x,U, "m--x")
ylabel('U(x)')
xlabel('x')
legend("Numerical solution of Heat Equation with three hats",'Location','north')
exportgraphics(gcf,'BE_Heateq_S_3.pdf')
close
