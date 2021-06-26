%Set constants:
a = -1;
b = 1;
alph = 5;
bet = 7;
N = [9, 49, 99, 999, 4999, 9999];

for n = N
    
    %Start CPU clock
    tStart = cputime;
    
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

    tEnd = cputime - tStart;
    %displays the cpu time
    cpu_display = ['For n= ',num2str(n), '  the CPU time was ', num2str(tEnd), ' seconds.'];
    disp(cpu_display)
end