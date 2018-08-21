% test example
% solving Poisson problem with Dirichlet boundary data
% -\Delta u = 1     in \Omega,
%    u(x,y) = 0     on \Gamma=\partial\Omega

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

clear all, close all, clc
format compact

% restoredefaultpath
addpath '../' 
addpath '../../refinement/'

%*** Geometry
elements   = [3,4,7,6; ...
              1,2,4,3; ...
              4,5,8,7];
coordinates = [-1,-1; 0,-1; -1,0; 0,0; 1,0; -1,1; 0,1; 1,1];
dirichlet   = [1,2;2,4;4,5;5,8;8,7;7,6;6,3;3,1];  neumann     = [];        E = 0.21407587e-0;            
%dirichlet   = [2,4;4,5;8,7;7,6;6,3;3,1];          neumann     = [5,8;1,2]; E = 0.2681285e-0;           
constrains  = zeros(0,4);

%*** Functions (volume force, Neuman and Dirichlet data)
uD = @(x)     zeros(size(x,1),1);
g  = @(x)     zeros(size(x,1),1);
f  = @(x)     ones(size(x,1),1);

%*** Compute solution
nEmax = 5000;
rho = 0.5;
[x,coordinates,elements,constrains,indicators,data] = adaptiveAlgorithm(...
       coordinates,elements,constrains,dirichlet,neumann,f,g,uD,nEmax,rho); 

%*** Plot convergence
figure(1)
loglog(data(:,1),data(:,2),'*-', ...
       data(:,1),sqrt(E-data(:,3)),'o:', ...
       data(:,1),2.5*data(:,1).^(-1/2))
%loglogTriangle(100,60000,10^(-1.2),0.5,'l')
%loglogTriangle(100,60000,10^(0),0.5,'u')
legend('estimator','error')
legend('estimator','error')
xlabel('number of elements')
ylabel('energy error')    

%*** Plot solution
if size(coordinates,1) < 1e5
  figure(2)
  x = x./max(abs(x));
patch('Faces',elements,'Vertices',coordinates,'FaceVertexCData',...
        x,'Facecolor','interp')
end  

