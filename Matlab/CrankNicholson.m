function [h,k,error] = CNmod(m)
% Solve u_t = kappa * u_{xx} on [ax,bx] with Dirichlet boundary conditions,
% using the Crank-Nicolson method with m interior points.

clf              %clear graphics
hold on          %Put all plots on the same graph (comment out if desired)

ax = -1;
bx = 1;
kappa = .02;               %heat conduction coefficient
tfinal = 1;                %final time

h = (bx-ax)/(m+1);             % h = delta x
x = linspace(ax,bx,m+2)';   %note x(1)=-1 and x(m+2)=1
                            %u(1)=g0 and u(m+2)=g1 are known from BC's
% k = 4*h;                  %k = Delta t
% k = h^2/kappa;
k = h;

nsteps = round(tfinal/k);    % number of time steps
% nplot = 1;      % plot solution every nplot time steps
                         % (set nplot=2 to plot every 2 time steps, etc.)
nplot = nsteps;  % only plot at final time

if abs(k*nsteps - tfinal) > 1e-5
   % The last step won't go exactly to tfinal.
   disp(' ')
   display(['WARNING *** k does not divide tfinal, k = ', num2str(k)]);
   disp(' ')
end


% True solution
utrue = @(t,x) 0.5 * erfc( x/(sqrt(4*kappa*t)) );

syms z;
y = piecewise(z<0, 1, z>=0, 0) ;
f = symfun(y,z);
u_0 = f(x);
u_0 = double(u_0);

% Each time step we solve MOL system U' = AU + g using the Trapezoidal method

% set up matrices:
r = (1/2) * kappa* k/(h^2);
e = ones(m,1);
A = spdiags([e -2*e e], [-1 0 1], m, m);
A1 = eye(m) - r * A;
A2 = eye(m) + r * A;

% initial data on fine grid for plotting:
xfine = linspace(ax,bx,1001);

% initialize u and plot:
tn = 0;
u = u_0;

% plot(x,u,'b.-', xfine,ufine,'r')
% legend('computed','true')
% title('Initial data at time = 0')

plot(x,u,'g.-')
legend('Initial Data')
% title('Initial data at time = 0')

axis padded

% exportgraphics(gcf,'CN_Heateq_Init.pdf')

input('Hit <return> to continue  ');

% main time-stepping loop:

for n = 1:nsteps
     tnp = tn + k;   % = t_{n+1}

     % boundary values u(0,t) and u(1,t) at times tn and tnp:

     g0n = u(1);
     g1n = u(m+2);
     g0np = utrue(tnp,ax);
     g1np = utrue(tnp,bx);

     % compute right hand side for linear system:
     uint = u(2:(m+1));   % interior points (unknowns)
     rhs = A2*uint;
     % fix-up right hand side using BC's (i.e. add vector g to A2*uint)
     rhs(1) = rhs(1) + r*(g0n + g0np);
     rhs(m) = rhs(m) + r*(g1n + g1np);

     % solve linear system:
     uint = A1\rhs;

     % augment with boundary values:
     u = [g0np; uint; g1np];

     % plot results at desired times:
     if mod(n,nplot)==0 || n==nsteps
        ufine = utrue(tnp, xfine);
        plot(x,u,'b.-', xfine,ufine,'r')
        legend('Initial Data','Numerical','Exact')
        title(['t = ', num2str(tnp), ' after ', num2str(n),
         ' time steps with ', num2str(m+2), ' grid points']);
        error = max(abs(u-utrue(tnp,x)));
        display(['at time t = ', num2str(tnp), ' max error = ', num2str(error)])
        if n<nsteps
            input('Hit <return> to continue  ');
        end
     end

     tn = tnp;   % for next time step

     axis padded

%      exportgraphics(gcf,'CN_Heateq_final38.pdf')

end
