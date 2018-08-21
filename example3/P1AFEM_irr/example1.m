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
clear all
format compact

addpath '../../refinement/'

% define initial mesh
elements = [1,2,3;1,3,4;1,4,5;1,5,6;1,6,7;1,7,8];
coordinates = [0,0; 1,0; 1,1; 0,1; -1,1;-1,0;-1,-1;0,-1];
dirichlet = [1,2;2,3;3,4;4,5;5,6;6,7;7,8;8,1];
constrains = zeros(0,3);
neumann = [];

uD     = @(x)   zeros(size(x,1),1);
g      = @(x)   zeros(size(x,1),1);
f      = @(x)   ones(size(x,1),1);
E = 0.21407587e-0;

nEmax = 5000; % number of maximal elements
rho = 0.5; % paramter in Dörfler marking
[x,coordinates,elements,constrains,indicators,data] = ...
    adaptiveAlgorithm(coordinates,elements, ...
    constrains,dirichlet,neumann,f,g,uD,nEmax,rho);

%*** Plot convergence
figure(1)
loglog(data(:,1),data(:,2),'*-', ...
    data(:,1),sqrt(E-data(:,3)),'o:')
% loglogTriangle(100,60000,10^(-1.2),0.5,'l')
% loglogTriangle(100,60000,10^(0),0.5,'u')
legend('estimator','error')
xlabel('number of elements')
ylabel('energy error')

%*** Plot solution
if size(coordinates,1) < 1e5
    figure(2)
    patch('Faces',elements,'Vertices',coordinates,'FaceVertexCData',...
        x,'Facecolor','interp')
end


