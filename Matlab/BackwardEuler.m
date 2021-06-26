%Backward Euler function
function err = BackwardEuler(a, b, h, k, t_0, T, sgm)

    m = (b-a)/h - 1;
    Psi = (k*sgm)/(h^2);

    %Generate the matrix A (size mxm) to be applied in BTCS:
    A = zeros(m);  %initialize mxm matrix
%     A = spalloc(m,m,3*m);   % initialize to zero matrix
    A(1,1) = 1+Psi;
    A(m,m) = 1+Psi;

    for i = 1:m
        for j = 1:m
            if (i == j) &&  ( (i ~= 1) && (i ~= m) )
                A(i,i) = 1 + 2*Psi;
            elseif (j == i+1) || (i == j+1)
                A(i,j) = -Psi;
            end
        end
    end

    %function f & exact solution u_exact
    f  = @(t,x) cos(t) * sin(x.^2) - 2 * sgm * sin(t) * ...
         ( cos(x.^2) - 2 * x.^2 * sin(x.^2) );
    u_true  = @(x) sin(T) .* sin(x.^2);
    u_true_vec = zeros(m,1);    %vectorize true solution...initializing
    rhs = zeros(m,1);           %initialize rhs of Eq 12e) from the problem set


    %------------------------------------------
    %        BACKWARD EULER CODE
    %------------------------------------------

    u_0 = zeros(m,1);                  %Initial condition
    it_max = ceil( (T - t_0)/k );

    for n = 0 : it_max

        Xi = 4*h*Psi*cos(4) * sin( t_0+(n+1)*k ); %\Xi function from Neumann BC

        for i = 1 : m
            if i == m
                rhs(i) = u_0(i) + Xi + k * f(t_0 + (n+1)*k, a + i*h);
            else
                rhs(i) = u_0(i) + k * f(t_0 + (n+1)*k, a + i*h);
            end

            if n == it_max
                u_true_vec(i) = u_true(a + i*h);
            end
        end

        u = A\rhs;          %solve Au = rhs

        if n == it_max
            %extend solution to include boundaries
            u_true_vec_full = [u_true(a); u_true_vec; u_true(b)];
            u_full = [u(1); u; 4*h*cos(4) * sin(t_0+n*k) + u(m)];

            err = norm(u_full - u_true_vec_full, Inf);    %global error (output of function)

            x = linspace(a,b,m+2);     %x-grid (for plotting purposes)

            %Plot results:
            plot(x,u_full, "r--x")
            hold on
            plot(x, u_true_vec_full, "g--o")
            ylabel('u(x)')
            xlabel('x')
            legend("Numerical Solution", "Exact Solution", 'Location','northwest')
            hold off
    %             exportgraphics(gcf,'BE_Heateq.pdf')
            shg
        end

        u_0 = u;  %update u_0 value for next iteration

    end

    %----------------------------------------------
    %      END OF BACKWARD EULER CODE
    %----------------------------------------------

end
