%Solve an elliptic equation of form -div(c*grad(u))+a*u=f
%on Omega = [-1,1]^2 with u=0 on the boundary

clear all;
close all;

c = 1;

%function a (a = 4 + xy)
a = @(location,state) 4 + (location.x .* location.y);

%function f (f = 5 e^y cos((3/2)x)
f = @(location,state) 5 .* exp(location.y) .* cos( (3/2) .* location.x );

model = createpde(1);
geometryFromEdges(model,@squareg);

% apply zero Dirichlet boundary conditions to all edges
applyBoundaryCondition(model,'dirichlet','Edge',...
                            1:model.Geometry.NumEdges,'u',0);

specifyCoefficients(model,'m', 0, ...
                          'd', 0, ...
                          'c', c,...
                          'a', a,...
                          'f', f);

% generate a finite element mesh (hmax specifies maximum mesh 'size')
% to get a finer mesh choose a smaller 'Hmax'
generateMesh(model,'GeometricOrder','linear','Hmax',0.03);

figure(1);
% plot the mesh (show element labels and node labels)
% pdeplot(model,'ElementLabels','on', 'NodeLabels','on');
% plot the mesh without the element labels
pdeplot(model);
axis('tight');
% exportgraphics(gcf,'mesh.pdf')



% solve the problem and get the nodal solution
results = solvepde(model);
u = results.NodalSolution;

% or alternatively
% obtain finite element matrices
FEM = assembleFEMatrices(model, 'nullspace');
% solve the problem
u1 = FEM.B*(FEM.Kc\FEM.Fc) + FEM.ud;

% surface (3D) plot of the solution
figure(2);
pdeplot(model,'XYData',u,'ZData',u,'Colormap','jet', 'Mesh','on');
% exportgraphics(gcf,'3dplot.pdf')




% 2D plot
figure(3);
pdeplot(model,'XYData',u,'Colormap','jet', 'Mesh','on');
% exportgraphics(gcf,'contour.pdf')
