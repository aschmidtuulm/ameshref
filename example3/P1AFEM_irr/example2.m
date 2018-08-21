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

clear all;close all;clf
format compact

addpath '../../refinement/'
addpath '../'
%*** Geometry
elements = [1,2,3;1,3,4;1,4,5;1,5,6;1,6,7;1,7,8];
coordinates = [0,0; 1,0; 1,1; 0,1; -1,1;-1,0;-1,-1;0,-1];
dirichlet = [1,2;2,3;3,4;4,5;5,6;6,7;7,8;8,1];
constrains = zeros(0,3);
neumann = []; 

uD     = @(x)   uD_exact(x);        
g      = @(x)   zeros(size(x,1),1); 
f      = @(x)   zeros(size(x,1),1);                 

nEmax = 5000; 
rho = 0.5;
[x,coordinates,elements,constrains,indicators,data] = ...
                           adaptiveAlgorithm(coordinates,elements, ...
                           constrains,dirichlet,neumann,f,g,uD,nEmax,rho,uD);
                         
%*** Plot convergence
figure(1)
loglog(data(:,1),data(:,2),'*-', ...
       data(:,1),data(:,5),'o:', ...
       data(:,1),data(:,6),'o:', ...
       data(:,1),data(:,7),'o:', ...
       data(:,1),0.8*data(:,1).^(-1/2), ...
       data(:,1),0.35*data(:,1).^(-1))
legend('estimator','H^1-error','L^2-error','L^{oo}-error')
xlabel('number of elements')
ylabel('error')
%*** Plot solution
if size(coordinates,1) < 1e5
  figure(2)
    patch('Faces',elements,'Vertices',coordinates,'FaceVertexCData',...
        x,'Facecolor','interp')
end 

