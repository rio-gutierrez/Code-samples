function F = LoadVecD0(x, f)
%input mesh vector x and function f to  LoadVecD0
%output Load Vector F


  n = length(x)-1;        %number of subintervals
  F = zeros(n-1, 1);      %allocate load vector
  %No need for half-hats due to vanishing BCs;
  %otherwise F would have dim n+1

  for i = 1:n-1
      h_minus = x(i+1) - x(i);                  %h_i        (index offset)
      xmid = (x(i+1) + x(i))/2;                 %m_i        (index offset)
      h_plus = x(i+2) - x(i+1);                 %h_{i+1}    (index offset)
      xmid_plus = (x(i+2) + x(i+1))/2;          %m_{i+1}    (index offset)

      F(i) = (h_minus/6) * ( f(x(i+1)) + 2*f(xmid) )  + ...
             (h_plus/6) * ( f(x(i+1)) + 2*f(xmid_plus) );
  end

end
