function [x,energy] = solveLaplace(coordinates,elements3,elements4, ...
                                   dirichlet,neumann,f,g,uD)
%solveLaplace: computes P1Q1-finite element solution for the two dimensional
%              Laplace equation with mixed Dirichlet-Neumann boundary 
%              condition
%
%    solveLaplace solves Laplace equation 
%      - div(grad(u)) = f                   in Omega
%                   u = u_D                 on the Dirichlet boundary
%              d/dn u = g                   on the Neumann boundary
%    on a geometry described by regular meshes consisting of triangles or
%    quadrilaterals
%
%Usage:
%
%[x,energy] = solveLaplace(coordinates,elements,dirichlet,neumann,f,g,ud)
%
%Comments:
%
%    solveLaplace expects as input a finite element mesh described by the 
%    fields coordinates, elements3 ,elements4 , dirichlet and neumann. The volume
%    force f, the Neumann data g, and the (inhomogeneous) Dirichlet data
%    uD are given as M-files <f.m>, <g.m>, and <uD.m>. Either of these 
%    M-files is assumed to take n evaluation points as (n x 2) matrix and to
%    return an (n x 1) column vector. 
%
%    solveLaplace assembles the Galerkin data and solves the resulting 
%    linear system of equations to obtain the P1-Q1 finite element solution 
%    of the Laplace problem. The function returns a column vector X which
%    contains the nodal values of the FEM solution. Additionally, solveLaplace 
%    provides the energy of the discrete solution uh, i.e. 
%    energy = || grad(uh) ||_{L2(Omega)}^2.
%
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
                                 
nC = size(coordinates,1);
nE3 = size(elements3,1);
nE4 = size(elements4,1);
%*** Assembly of stiffness matrix for triangular elements
%*** First vertex of elements and corresponding edge vectors 
c1 = coordinates(elements3(:,1),:);
d21 = coordinates(elements3(:,2),:) - c1;
d31 = coordinates(elements3(:,3),:) - c1;
%*** Vector of element areas 4*|T|
area4 = 2*(d21(:,1).*d31(:,2)-d21(:,2).*d31(:,1));
%*** Assembly of stiffness matrix
I3 = reshape(elements3(:,[1 2 3 1 2 3 1 2 3])',9*nE3,1);
J3 = reshape(elements3(:,[1 1 1 2 2 2 3 3 3])',9*nE3,1);
a = (sum(d21.*d31,2)./area4)';
b = (sum(d31.*d31,2)./area4)';
c = (sum(d21.*d21,2)./area4)';
A3 = [-2*a+b+c;a-b;a-c;a-b;b;-a;a-c;-a;c];
%*** Tensor Gauss quadrature
[w,s] = gauss(4);
[xg,yg] = meshgrid(s,s); wg = w*w'; 
xg = xg(:); yg = yg(:); wg = wg(:);
%*** Assembly of stiffness matrix for quadrilateral elements
phi  = [(1-xg).*(1-yg),(1+xg).*(1-yg),(1+xg).*(1+yg),(1-xg).*(1+yg)]/4;
I4 = reshape(elements4(:,repmat(1:4,4,1) ),16*nE4,1);
J4 = reshape(elements4(:,repmat(1:4,4,1)'),16*nE4,1);
xx = reshape(coordinates(elements4',1),4,[]);
yy = reshape(coordinates(elements4',2),4,[]);
A4 = zeros(nE4,16);
L4 = zeros(nE4,4);
for k = 1:size(xg,1)
  dphi_ds = [yg(k)-1,  1-yg(k), 1+yg(k), -(1+yg(k))]/4;
  dphi_dt = [xg(k)-1,-(1+xg(k)),1+xg(k),   1-xg(k) ]/4;
  DF_ds = [dphi_ds * xx; dphi_ds * yy]';
  DF_dt = [dphi_dt * xx; dphi_dt * yy]';
  det_DF = sum((DF_ds *[0,1;-1,0 ]).*DF_dt,2);
  a_1 = wg(k) * sum(DF_dt.^2,2)    ./det_DF;
  b_1 = wg(k) * sum(DF_ds.*DF_dt,2)./det_DF;
  c_1 = wg(k) * sum(DF_ds.^2,2)    ./det_DF;
  L4 = L4 + (det_DF.*wg(k).*f([phi(k,:)*xx;phi(k,:)*yy]'))*phi(k,:);
  A4 = A4 + ( a_1 * reshape(dphi_ds' * dphi_ds,1,[]) ...
            + c_1 * reshape(dphi_dt' * dphi_dt,1,[]) ...
            - b_1 * reshape(dphi_ds' * dphi_dt + dphi_dt' * dphi_ds,1,[]));
end
A = sparse([I3;I4],[J3;J4],[A3(:);A4(:)]);
%*** Prescribe values at Dirichlet nodes
x = zeros(nC,1);
dirichlet = unique(dirichlet);
x(dirichlet) = feval(uD,coordinates(dirichlet,:));
%*** Assembly of right-hand side
fsT = feval(f,c1+(d21+d31)/3);
b = accumarray(elements4(:),L4(:),[nC 1]) ...
  + accumarray(elements3(:),repmat(12\area4.*fsT,3,1),[nC 1]) - A * x;
%*** Assembly of Neumann load
if ~isempty(neumann)
    cn1 = coordinates(neumann(:,1),:);
    cn2 = coordinates(neumann(:,2),:);
    gmE = feval(g,(cn1+cn2)/2);
    b = b + accumarray(neumann(:),...
                      repmat(2\sqrt(sum((cn2-cn1).^2,2)).*gmE,2,1),[nC 1]);
end
% Computation of P1-Q1-FEM approximation
freenodes = setdiff(1:nC, dirichlet);
x(freenodes) = A(freenodes,freenodes) \ b(freenodes); 
%*** Compute energy || grad(uh) ||^2 of discrete solution
energy = x(1:nC)'*A(1:nC,1:nC)*x(1:nC);
