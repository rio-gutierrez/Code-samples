function advection_up_Prob3(m)

global a
a = 1;        %advection velocity

clf           %clear graphics
% set(gcf,'PaperType','A4')

ax = 0;
bx = 2;
tfinal = 0.5;                 %final time

h = (bx-ax)/(m+1);            %h = delta x (m=99 for this exercise)
k = 0.5*h;                    %time step

nu = a*k/h;                 %Courant number
x = linspace(ax,bx,m+2)';   %note x(1)=0 and x(m+2)=1
                            %With periodic BC's there are m+1 unknowns u(2:m+2)

I = 2:(m+2);                %indices of unknowns

nsteps = round(tfinal / k);   %number of time steps
% nplot = nsteps;             %only plot at final time

if abs(k*nsteps - tfinal) > 1e-5
   % The last step won't go exactly to tfinal.
   disp(' ')
   display(['  WARNING *** k does not divide tfinal, k = ', num2str(k)])
   disp(' ')
end

% initial conditions:
tn = 0;
u0 = eta(x);
u = u0;

% periodic boundary conditions:
u(1) = u(m+2);   % copy value from rightmost unknown to ghost cell on left
u(m+3) = u(2);   % copy value from leftmost unknown to ghost cell on right

% plot initial data:
plot(x,u0,'rx')
hold on


%create movie for the plot from t_0 to t_final
v = VideoWriter('Upwind_movie.avi');
open(v);

% main time-stepping loop:
for n = 1:nsteps
     tnp = tn + k;   % = t_{n+1}

     % Upwind:
     u(I) = u(I) - nu*(f(u(I)) - f(u(I-1))) ;

     % periodic boundary conditions:
     u(1) = u(m+2);   % copy value from rightmost unknown to ghost cell on left
     u(m+3) = u(2);   % copy value from leftmost unknown to ghost cell on right

      uint = u(1:m+2);      %points on the interval (drop ghost cell on right)

      plot(x,uint, 'c.-');
      title(['Time =  ', num2str(tnp,6), ' after ', ...
      num2str(n), ' time steps with ', num2str(m+1), ' grid points']);
      legend('Initial Data', 'Numerical Solution')
      axis
      pause(0.1);
      frame = getframe(gcf);
      writeVideo(v,frame);

     tn = tnp;   % for next time step

end

exportgraphics(gcf,'UW_Prob3.pdf');

close(v);



function f = f(x)
        f = 0.5 .* x.^2;
return



function eta = eta(x)                %Initial data

    eta = zeros(length(x), 1);
    for i = 1 : length(x)
        if ( (x(i) >= 0) && (x(i) <= 1) )
            eta(i) = 1 - cos(2*pi.*x(i));
        end
    end

return
