function [x,energy] = solveLaplace(coordinates,elements,constrains, ...
                                   dirichlet,neumann,f,g,uD)
%solveLaplace: computes Q1-finite element solution for the two dimensional
%              Laplace equation with mixed Dirichlet-Neumann boundary 
%              condition
%
%    solveLaplace solves Laplace equation 
%      - div(grad(u)) = f                   in Omega
%                   u = u_D                 on the Dirichlet boundary
%              d/dn u = g                   on the Neumann boundary
%    on a geometry described by quadrilaterals with hanging nodes. 
%
%Usage:
%
%[x,energy] = solveLaplace(coordinates,elements,dirichlet,neumann,f,g,ud)
%
%Comments:
%
%    solveLaplace expects as input a finite element mesh described by the 
%    fields coordinates, elements,constrains, dirichlet and neumann. The volume
%    force f, the Neumann data g, and the (inhomogeneous) Dirichlet data
%    uD are given as M-files <f.m>, <g.m>, and <uD.m>. Either of these 
%    M-files is assumed to take n evaluation points as (n x 2) matrix and to
%    return an (n x 1) column vector. 
%
%    solveLaplace assembles the Galerkin data and solves the resulting 
%    linear system of equations to obtain the P1 finite element solution 
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
nE = size(elements,1);
%*** Tensor Gauss quadrature
[w,s] = gauss(4);
[xg,yg] = meshgrid(s,s); wg = w*w'; 
xg = xg(:); yg = yg(:); wg = wg(:);
%*** Assembly of stiffness matrix and volume force of right-hand side 
%    for quadrilateral elements
phi  = [(1-xg).*(1-yg),(1+xg).*(1-yg),(1+xg).*(1+yg),(1-xg).*(1+yg)]/4;
I = reshape(elements(:,repmat(1:4,4,1) ),16*nE,1);
J = reshape(elements(:,repmat(1:4,4,1)'),16*nE,1);
xx = reshape(coordinates(elements',1),4,[]);
yy = reshape(coordinates(elements',2),4,[]);
A = zeros(nE,16);
L = zeros(nE,4);
for k = 1:size(xg,1)
  dphi_ds = [yg(k)-1,  1-yg(k), 1+yg(k), -(1+yg(k))]/4;
  dphi_dt = [xg(k)-1,-(1+xg(k)),1+xg(k),   1-xg(k) ]/4;
  DF_ds = [dphi_ds * xx; dphi_ds * yy]';
  DF_dt = [dphi_dt * xx; dphi_dt * yy]';
  det_DF = sum((DF_ds *[0,1;-1,0 ]).*DF_dt,2);
  a_1 = wg(k) * sum(DF_dt.^2,2)    ./det_DF;
  b_1 = wg(k) * sum(DF_ds.*DF_dt,2)./det_DF;
  c_1 = wg(k) * sum(DF_ds.^2,2)    ./det_DF;
  L = L + (det_DF.*wg(k).*f([phi(k,:)*xx;phi(k,:)*yy]'))*phi(k,:);
  A = A + ( a_1 * reshape(dphi_ds' * dphi_ds,1,[]) ...
          + c_1 * reshape(dphi_dt' * dphi_dt,1,[]) ...
          - b_1 * reshape(dphi_ds' * dphi_dt + dphi_dt' * dphi_ds,1,[]));
end
A = sparse(I,J,A(:));
b = accumarray(elements(:),L(:),[nC 1]); 
%*** Prescribe values at Dirichlet nodes
x = zeros(nC,1);
fixed = unique(dirichlet);
x(fixed) = feval(uD,coordinates(fixed,:));
%*** Assembly of Neumann load
if ~isempty(neumann)
    cn1 = coordinates(neumann(:,1),:);
    cn2 = coordinates(neumann(:,2),:);
    gmE = feval(g,(cn1+cn2)/2);
    b = b + accumarray(neumann(:),...
                      repmat(2\sqrt(sum((cn2-cn1).^2,2)).*gmE,2,1),[nC 1]);
end
%*** Remove constrained nodes
if ~isempty(constrains)
  tconstrains = constrains;
  I = speye(size(A));
  while ~isempty(tconstrains)
    idx = logical(ones(nC,1)); 
    idx(tconstrains(:,3:end)) = 0;
    jdx = idx(tconstrains(:,1)) & idx(tconstrains(:,2)); 
    if size(constrains,2) == 3
      I = (speye(size(A))+sparse([tconstrains(jdx,3);tconstrains(jdx,3)], ...
                                 [tconstrains(jdx,1);tconstrains(jdx,2)], ...
                              ones(2*sum(jdx),1)/2,size(A,1),size(A,2)))*I ;
    else
      I = (speye(size(A))+sparse([tconstrains(jdx,3);tconstrains(jdx,3)], ...
                                 [tconstrains(jdx,1);tconstrains(jdx,2)], ...
           [ones(sum(jdx),1)/3,2*ones(sum(jdx),1)/3],size(A,1),size(A,2)) ...
                         +sparse([tconstrains(jdx,4);tconstrains(jdx,4)], ...
                                 [tconstrains(jdx,1);tconstrains(jdx,2)], ...
           [2*ones(sum(jdx),1)/3,ones(sum(jdx),1)/3],size(A,1),size(A,2)))*I ;
    end    
    tconstrains = tconstrains(~jdx,:);                      
  end
  A = I'*A*I;
  b = I'*b;
  fixed = [fixed;reshape(constrains(:,3:end),[],1)];
end
%*** Computation of Q1-FEM approximation
freenodes = setdiff(1:size(A,1), fixed);
b = b - A*x;
x(freenodes) = A(freenodes,freenodes)\b(freenodes);
if ~isempty(constrains)
    x(constrains(:,3:end)) = 0;
    x = I * x;
end
%*** Compute energy || grad(uh) ||^2 of discrete solution
energy = full(x(freenodes)'*A(freenodes,freenodes)*x(freenodes));


