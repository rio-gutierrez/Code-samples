%Set constants:
a = -1;
b = 1;
alph = 5;
bet = 7;
N = [24,49,99,199];

%initialize vectors of norms to be used later:
norml1_vec =  zeros(size(N,2),1);
norml2_vec =  zeros(size(N,2),1);
normlinf_vec =  zeros(size(N,2),1);

for n = N
      
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

    %define the x grid:
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

    %closed-form solution
    funct = @(t) -(abs(t)^3)/6 + t^2/2 + t + 17/3; 

    
    %generate vector of errors
    error_vec =  zeros(n,1);
    for i = 1:n
        error_vec(i) = usol(i) - funct(x(i));
    end

    l1 = h* norm(error_vec,1);
    l2 = sqrt(h)* norm(error_vec,2);
    l_inf = norm(error_vec,Inf);
    

    %fill in vectors of norms:
    it =  find(N==n);      %get the n-index of the tuple N
    norml1_vec(it) = l1;
    norml2_vec(it) = l2;
    normlinf_vec(it) = l_inf;
    
%     %displays the norms, for each n
%     norm_display_1 = ['The L1 norm for n= ',num2str(n), ' is ', num2str(l1)];
%     norm_display_2 = ['The L2 norm for n= ',num2str(n), ' is ', num2str(l2)];
%     norm_display_inf = ['The L^inf norm for n= ',num2str(n), ' is ', num2str(l_inf)];
%     disp(norm_display_1)
%     disp(norm_display_2)
%     disp(norm_display_inf)
end


% %Get the slopes
slope_l1 = ( norml1_vec(end-1) -  norml1_vec(end-2) )/ ( N(end-1) - N(end-2));
slope_l2 = ( norml2_vec(end-1) -  norml2_vec(end-2) )/ ( N(end-1) - N(end-2));
slope_linf = ( normlinf_vec(end-1) -  normlinf_vec(end-2) )/ ( N(end-1) - N(end-2));


% %Plot results:
% loglog(N, norml1_vec, "ms")
% hold on
% loglog(N, norml2_vec, "bd")
% loglog(N, normlinf_vec, "ro")
% ylabel('Error')
% xlabel('Number of grid points (n)')
% legend("1-norm","2-norm", "\infty-norm")
% hold off
% exportgraphics(gcf,'norms_loglog.pdf')
% close