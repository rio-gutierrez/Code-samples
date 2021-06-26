%Rayleigh Quotient Iteration
%Input: matrix A, initial (nonzero) vector x, number of steps k
%Output: eigenvalue lam and

function [lam,u] = RQI(A,x,k)

    for j = 1:k
      u = x/norm(x);                       %normalize
      lam = u'* A * u;                     %Rayleigh Quotient
      x = (A - lam * eye(size(A))) \ u;    %Inverse power iteration
    end

    u = x/norm(x);
    lam = u' * A * u;                      %Rayleigh quotient

end
