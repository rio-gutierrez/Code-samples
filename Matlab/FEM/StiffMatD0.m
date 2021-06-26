function K = StiffMatD0(x, k)
%input mesh vector x and function k(x) to  StiffMatD0
%output Stiffness Matrix K


  n = length(x)-1;         %number of subintervals
  K = zeros(n-1, n-1);     %allocate stiffness matrix
  %No need for half-hats due to vanishing BCs;
  %otherwise M would have dim (n+1)x(n+1)

  for i = 1:n-1
      h_minus = x(i+1) - x(i);                   %h_i       (index offset)
      xmid = (x(i+1) + x(i))/2;                  %m_i       (index offset)
      h_plus = x(i+2) - x(i+1);                  %h_{i+1}   (index offset)
      xmid_plus = (x(i+2) + x(i+1))/2;           %m_{i+1}   (index offset)

      K(i,i) = (1/(6*h_minus)) * ( k(x(i)) + 4*k(xmid) + k(x(i+1)) )  + ...
               (1/(6*h_plus)) * ( k(x(i+1)) + 4*k(xmid_plus) + k(x(i+2)) ) ;

      if i ~= n-1
          K(i+1,i)  =  - (1/(6*h_plus)) * ( k(x(i+1)) + ...
                       4*k(xmid_plus) + k(x(i+2)) );
          K(i,i+1) = K(i+1,i);
      end
  end

end
