% test example
% solving Poisson problem with Dirichlet boundary data
% -\Delta u = 0                               in \Omega,
%    u(r,alpha) = r^(2/3)* sin(2/3*alpha)     on \Gamma=\partial\Omega
% with u_ex = r^(2/3)* sin(2/3*alpha).

%Remark:
%
%    This program is a supplement to the paper 
%    >> Adaptive Mesh Refinement in 2D - An Efficient Implementation in Matlab <<
%    by S. Funken, and A. Schmidt. The reader should 
%    consult that paper for more information.   
%
%Authors:
% 
%    S. Funken, A.Schmidt  21-08-18
clear all;close all; clf;
format compact

% restoredefaultpath
addpath('../')
addpath('../Refinement/')

%*** Geometry
elements   = [3,4,7,6; ...
              1,2,4,3; ...
              4,5,8,7];
coordinates = [-1,-1; 0,-1; -1,0; 0,0; 1,0; -1,1; 0,1; 1,1];
dirichlet   = [1,2;2,4;4,5;5,8;8,7;7,6;6,3;3,1];  neumann     = [];               
constrains  = zeros(0,3);  % zeros(0,3) for QrefineR or zeros(0,4)for QrefineR2; 

%*** Funtions (volume force, Neuman and Dirichlet data)
uD     = @(x)   uD_exact(x);  %x(:,1).^3;                 
g      = @(x)   zeros(size(x,1),1); 
f      = @(x)   zeros(size(x,1),1);                 

%*** Compute solution
nEmax = 5000;
rho = 0.5;
[x,coordinates,elements,constrains,indicators,data] = adaptiveAlgorithm(...
       coordinates,elements,constrains,dirichlet,neumann,f,g,uD,nEmax,rho,uD); 

%*** Plot convergence
figure(1)
loglog(data(:,1),data(:,2),'*-', ...
       data(:,1),data(:,5),'o:', ...
       data(:,1),data(:,6),'o:', ...
       data(:,1),.5*data(:,1).^(-1/2), ...
       data(:,1),0.15*data(:,1).^(-1))
legend('estimator','H1-error','L2-error') 
xlabel('number of elements')
ylabel('error')

%*** Plot solution
if size(coordinates,1) < 1e5
  figure(2)
  x = x./max(abs(x));
  patch('Faces',elements,'Vertices',coordinates,'FaceVertexCData',...
        x,'Facecolor','interp')
end  
