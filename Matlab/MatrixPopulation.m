%Periodic matrix population model
%(Exercise from Problem set 4, Math Modeling 2)

%Import parameters
A1 = 0.125;
A2 = 0.2458;
A3 = 0.586;
A4 = 0.85;
F = 45;
L = 0.03426;
C = 0.85;

n = 3;

%Initialize matrices
M1 = zeros(n);
M2 = zeros(n);
M3 = zeros(n);
M4 = zeros(n);


%Complete matrices
M1(1,n) = F;
M1(2,1) = L;
M1(n,2) = C;
M1(n,n) = A1;

M2(1,n) = F;
M2(2,1) = L;
M2(n,2) = C;
M2(n,n) = A2;

M3(n,n) = A3;
M4(n,n) = A4;


%Build matrix M  (for Problem 1)
M = M4 * M3 * (M2^6) * (M1^6);

%Build matrix \tilde{M} (for Problem 2)
Mt =  (M1^6) * M4 * M3 * (M2^6);

%Build matrix \hat{M} (for Problem 4)
h_l = 0.1;
h_c = 0.5;
h_a = 0.5;
h = [h_l; h_c; h_a];
H = diag(h);
Mh =  M4 * M3 * (M2^6) * (M1^4) * H * (M1^2);



%------------------------------------------------------------------
%     CODE FOR SENSITIVITY & ELASTICITY MATRICES
%------------------------------------------------------------------

% %Initialize Sensitivity & Elasticity matrices for M and \tilde{M}
% S = zeros(n);
% E = zeros(n);
% St = zeros(n);
% Et = zeros(n);
%
% %eigenvalues and left & right eigenvectors for  M and \tilde{M}
% [V,D,W] = eig(M);
% [Vt,Dt,Wt] = eig(Mt);
%
% %find index where dominant eigenvalue lies
% index = find( diag(D) == max(diag(D)) );
% index_t = find( diag(Dt) == max(diag(Dt)) );
%
% %define the dominant eigenvalue
% lmb = D(index ,index);
% lmb_t = D(index_t, index_t);
%
% %define the dominant right (v) and left (w) eigenvectors
% v = V(:,index);
% w = W(:,index);
% vt = Vt(:,index_t);
% wt = Wt(:,index_t);
%
%
% for i = 1:n
%     for j = 1:n
%         if M(i,j) ~= 0
%             S(i,j) = ( w(i) * v(j) )/( w' * v );
%             E(i,j) = ( M(i,j)/lmb ) * S(i,j);
%         end
%          if Mt(i,j) ~= 0
%             St(i,j) = ( wt(i) * vt(j) )/( wt' * vt );
%             Et(i,j) = ( Mt(i,j)/lmb_t ) * St(i,j);
%         end
%     end
% end
%
%
% disp('The sensitivity matrix for M is ')
% S
% disp('The elasticity matrix for M is ')
% E
% disp('The sensitivity matrix for \tilde{M} is ')
% St
% disp('The elasticity matrix for \tilde{M} is ')
% Et

%---------------------------------------------------------------------------
%    END OF CODE FOR SENSITIVITY & ELASTICITY MATRICES
%----------------------------------------------------------------------------



% %---------------------------------------------------
% %    POWER ITERATION CODE
% %---------------------------------------------------
%
% x0 = [5000; 5000; 5000];    %initial guess
% x0 = x0/norm(x0, Inf) ;         %normalize
% x   = zeros(n,1);                    %initialize vector
% tol = 1e-6;                            %accepted tolerance
% it = 0;                                   %initialize number of iterations
% it_max = 100;                       %max number of iterations allowed
%
% for i = 1:it_max
%
%     it = it + 1;
%     x = Mt * x0;
%     x = x/norm(x, Inf);    %Normalize
%
%     if norm(x-x0) <= tol
%         lmb = (x' * Mt * x)/(x' * x);        %Rayleigh quotient
%         disp(['The dominant eigenvalue is ', num2str(lmb), '. It took ', num2str(it), ' iterations to converge.'])
%         disp('The dominant eigenvector is ')
%         disp(x)
%         break
%     end
%
%     x0 = x;  %update x0 value for next iteration
%
%     if it == it_max
%         disp('No convergence; max number of iterations reached.')
%     end
% end
%
% % ---------------------------------------------------
% %    END OF POWER ITERATION CODE
% % ---------------------------------------------------





% %--------------------------------------------------------------
% %    POWER ITERATION CODE  (APPLIED TO \HAT{M}
% %--------------------------------------------------------------

x0 = [5000; 5000; 5000];    %initial guess
x0 = x0/norm(x0, Inf) ;         %normalize
x   = zeros(n,1);                    %initialize vector
tol = 1e-6;                            %accepted tolerance
it = 0;                                   %initialize number of iterations
it_max = 100;                       %max number of iterations allowed

for i = 1:it_max

    it = it + 1;
    x = Mh * x0;
    x = x/norm(x, Inf);    %Normalize

    if norm(x-x0) <= tol
        lmb = (x' * Mh * x)/(x' * x);        %Rayleigh quotient
        disp(['The dominant eigenvalue is ', num2str(lmb), '. It took ', num2str(it), ' iterations to converge.'])
        disp('The dominant eigenvector is ')
        disp(x)
        break
    end

    x0 = x;  %update x0 value for next iteration

    if it == it_max
        disp('No convergence; max number of iterations reached.')
    end
end


% % ---------------------------------------------------
% %    END OF POWER ITERATION CODE
% % ---------------------------------------------------
