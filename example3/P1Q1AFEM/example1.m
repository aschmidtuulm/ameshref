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
clear all, close all, clf
format compact

% restoredefaultpath
addpath '../'
addpath '../../refinement/'

%*** Geometry
ref = 'TrefineNVB';   % 'QrefineRG2' 'QrefineRG' 'QrefineRB' 'TrefineRG' 'TrefineRGB' 'TrefineNVB'

if strcmp(ref,'TrefineRG') || strcmp(ref,'TrefineNVB') || strcmp(ref,'TrefineRGB')
    elements3   = [3,1,2;1,3,4;5,1,4;1,5,6;7,1,6;1,7,8];
    elements4   = zeros(0,4);
    coordinates = [0,0; 1,0; 1,1; 0,1; -1,1;-1,0;-1,-1;0,-1];
    dirichlet = [1,2;2,3;3,4;4,5;5,6;6,7;7,8;8,1];
    neumann = [];
else
    elements3   = zeros(0,3);
    elements4   = [3,4,7,6; ...
        1,2,4,3; ...
        4,5,8,7];
    coordinates = [-1,-1; 0,-1; -1,0; 0,0; 1,0; -1,1; 0,1; 1,1];
    dirichlet   = [1,2;2,4;4,5;5,8;8,7;7,6;6,3;3,1]; neumann     = [];
    
    
end
%*** Functions (volume force, Neumann and Dirichlet data)
uD = @(x)     0*(x(:,1)+x(:,2));
g  = @(x)     zeros(size(x,1),1);
f  = @(x)     ones(size(x,1),1);
E = 0.21407587;

nEmax = 5000;
rho = 0.5;
[x,coordinates,elements3,elements4,indicators,data] = adaptiveAlgorithm(...
    coordinates,elements3,elements4,dirichlet,neumann,f,g,uD,nEmax,rho,[],ref);

%*** Plot convergence
figure(1)
loglog(data(:,1),data(:,2),'*-', ...
    data(:,1),sqrt(E-data(:,3)),'o:')
%loglogTriangle(100,60000,10^(-1.2),0.5,'l')
%loglogTriangle(100,60000,10^(0),0.5,'u')
legend('estimator','error')
legend('estimator','error')
xlabel('number of elements')
ylabel('energy error')

%*** Plot solution
if size(coordinates,1) < 1e5
    figure(2)
    patch('Faces',elements4,'Vertices',coordinates,'FaceVertexCData',...
        x,'Facecolor','interp')
    hold on
    patch('Faces',elements3,'Vertices',coordinates,'FaceVertexCData',...
        x,'Facecolor','interp')
end

