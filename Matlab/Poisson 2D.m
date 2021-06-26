
% poisson2.m  -- solve the Poisson problem u_{xx} + u_{yy} = f(x,y)
% on [a,b] x [a,b].
%
% The 5-point Laplacian is used at interior grid points.
% This system of equations is then solved using backslash.
%

% a = 0;
% b = 1;
% m = 20;


%Rectangular (i.e., non-square) domain ( x \in [1,2]; y \in [0,1] )
ax = 1;
bx = 2;
ay = 0;
by = 1;

Mx = [24, 49, 99];   %corresponds to hx = 1/25, 1/50, and 1/100, respectively
My = [19, 39, 79];   %corresponds to hy = 1/20, 1/5=40, and 1/80, respectively

for mx = Mx
    for my = My
        hx = (bx-ax)/(mx+1);
        hy = (by-ay)/(my+1);
        x = linspace(ax,bx,mx+2);   % grid points x including boundaries
        y = linspace(ay,by,my+2);   % grid points y including boundaries


        [X,Y] = meshgrid(x,y);      % 2d arrays of x,y values
        X = X';                     % transpose so that X(i,j),Y(i,j) are
        Y = Y';                     % coordinates of (i,j) point

        Iint = 2:mx+1;              % indices of interior points in x
        Jint = 2:my+1;              % indices of interior points in y
        Xint = X(Iint,Jint);       % interior points
        Yint = Y(Iint,Jint);

        f = @(x,y) 1.25*exp(x+y/2);         % f(x,y) function

         rhs = f(Xint,Yint); %evaluate f at interior points for right hand side
                             %rhs is modified below for boundary conditions.

        utrue = exp(X+Y/2);        % true solution for test problem

        % set boundary conditions around edges of usoln array:

        usoln = utrue;        % use true solution for this test problem
                              % This sets full array, but only boundary values
                              % are used below.  For a problem where utrue
                              % is not known, would have to set each edge of
                              % usoln to the desired Dirichlet boundary values.


        % adjust the rhs to include boundary terms:
        rhs(:,1) = rhs(:,1) - usoln(Iint,1)/hy^2;
        rhs(:,my) = rhs(:,my) - usoln(Iint,my+2)/hy^2;
        rhs(1,:) = rhs(1,:) - usoln(1,Jint)/hx^2;
        rhs(mx,:) = rhs(mx,:) - usoln(mx+2,Jint)/hx^2;

        %convert the 2d grid function rhs into a column vector
        %for rhs of system:
        F = reshape(rhs,mx*my,1);

%         % form matrix A:
%         I = speye(mx,my);
%         e = ones(my,1);
%         T = spdiags([e -4*e e],[-1 0 1],mx,my);
%         S = spdiags([e e],[-1 1],mx,my);
%         A = (kron(I,T) + kron(S,I)) /(hx*hy);

        % form matrix A:
        Ix = speye(my);
        ex = ones(mx,1);
        Tx = spdiags([ex -2*ex ex],[-1 0 1],mx,mx) * (1/hx^2);

        Iy = speye(mx);
        ey = ones(my,1);
        Ty = spdiags([ey -2*ey ey],[-1 0 1],my,my) * (1/hy^2);

        A =  kron(Ix,Tx) + kron(Ty,Iy) ;

        % Solve the linear system:
        uvec = A\F;


        % reshape vector solution uvec as a grid function and
        % insert this interior solution into usoln for plotting purposes:
        % (recall boundary conditions in usoln are already set)

        usoln(Iint,Jint) = reshape(uvec,mx,my);

        usoln = reshape(usoln, (my+2)*(mx+2), 1);
        utrue = reshape(utrue, (my+2)*(mx+2), 1);


        %disp(['For hx = ', num2str(hx), ' and hy =  ', num2str(hy),
        %' the L^inf norm is ', num2str( max(abs(usoln-utrue)) )])
        disp(['For hx = ', num2str(hx), ' and hy =  ', num2str(hy),
        ' the L^2 norm is ', num2str(   norm(usoln-utrue,2)  * sqrt(hx*hy) )])



        % assuming true solution is known and stored in utrue:
    %     err = max(max(abs(usoln-utrue)));
    %     err = max(max(abs(usoln-utrue)));
    %     fprintf('Error relative to true solution of PDE = %10.3e \n',err)
    end
end

% % plot results as a surface
% figure(2)
% surf(X,Y,usoln)
