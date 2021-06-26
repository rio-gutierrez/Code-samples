%Forward Euler function
function err = ForwardEuler(a, b, h, k, t_0, T, sgm)

  m = (b-a)/h - 1;

  %function f & exact solution u_exact
  f  = @(t,x) cos(t) * sin(x.^2) - 2 * sgm * sin(t) * ...
       ( cos(x.^2) - 2 * x.^2 * sin(x.^2) );
  u_true  = @(x) sin(T) .* sin(x.^2);
  u_true_vec = zeros(m,1);            %vectorize true solution...initializing


  %------------------------------------------
  %        FORWARD EULER CODE
  %------------------------------------------

  u_0 = zeros(m,1);          %Initial condition
  u = u_0;                   %initialize solution vector

  Psi = (k*sgm)/(h^2);
  it_max = ceil( (T - t_0)/k );  %use ceiling in case t_0+nk
                                 %never equals T (it_max) exactly

  for n = 1 : it_max

      Xi = 4*h*cos(4) * sin( t_0+n*k );      %\Xi function from Neumann BC

      for i = 1 : m
          if i == 1
              u(i) = (1 - Psi) * u_0(i) + Psi * u_0(i+1) + ...
                     k * f(t_0 + n*k, a + i*h);
          elseif i == m
              u(i) = Psi * u_0(i-1) +  (1 - Psi) * u_0(i) + Xi*Psi + ...
                     k *  f(t_0 + n*k, a + i*h);
          else
              u(i) = Psi * u_0(i-1) +  (1 - 2*Psi) * u_0(i) + ...
                     Psi * u_0(i+1) + k *  f(t_0 + n*k, a + i*h);
          end

          if n == it_max
              u_true_vec(i) = u_true(a + i*h);
          end
      end

      if n == it_max
          %extend solution to include boundaries
          u_true_vec_full = [ u_true(a); u_true_vec; u_true(b)  ];
          u_full = [u(1); u; Xi + u(m)];

          %global error (output of function)
          err = norm(u_full - u_true_vec_full, Inf);

          x = linspace(a,b,m+2);     %x-grid (for plotting purposes)

          %Plot results:
          plot(x,u_full, "r--x")
          hold on
          plot(x, u_true_vec_full, "g--o")
          ylabel('u(x)')
          xlabel('x')
          legend("Numerical Solution", "Exact Solution", 'Location','northwest')
          hold off
%             exportgraphics(gcf,'FE_Heateq.pdf')
          shg
       end

      u_0 = u;  %update u_0 value for next iteration

  end

  %----------------------------------------------
  %     END OF FORWARD EULER CODE
  %----------------------------------------------

end
