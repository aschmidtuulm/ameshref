function [x,coordinates,elements,constrains,indicators,data] ...
    = adaptiveAlgorithm(coordinates,elements,constrains, ...
                        dirichlet,neumann,f,g,uD,nEmax,rho,u_ex)
%adaptiveAlgorithm  adaptive finite element algorithm for two-dimensional 
%                   Laplace equation
%
%    adaptiveAlgorithm computes the P1-AFEM solution of the Laplace 
%    equation for irregular meshes
%      - div(grad(u)) = f                   in Omega
%                   u = u_D                 on the Dirichlet boundary
%              d/dn u = g                   on the Neumann boundary
%    on a geometry described by triangles. 
%
%    adaptiveAlgorithm is the implementation of the following iterative
%    process for a given initial finite element mesh:
%
%          1) compute the discrete solution (via solveLaplace.m)
%          2) compute the elementwise error indicators (via computeEtaR.m)
%          3) mark the elements for refinement (via the Doerfler criterion)
%          4) refine the finite element mesh (via TrefineR.m)
%          5) return to 1)
%
%Usage:
%
%[x,coordinates,elements,constrains,indicators,data] ...
%    = adaptiveAlgorithm(coordinates,elements,constrains,dirichlet,neumann,f,g,uD,nemax,rho,u_ex)
%
%Comments:
%
%    adaptiveAlgorithm expects as input an initial finite element mesh 
%    described by the fields coordinates, elements, constrains, dirichlet
%    and neumann.
%
%    Volume force and boundary data are given as M-files <f.m>, <g.m>, and 
%    <uD.m>. Either of these M-files is assumed to take N evaluation
%    points as (N x 2) matrix and to return a (N x 1) column vector.
%
%    The stopping criterion is realized via the maximal number of elements 
%    NEMAX of the finite element mesh, and the computation is stopped a
%    soon as the current number of elements is larger than NEMAX.
%
%    The parameter RHO in (0,1) corresponds to the marking of elements by
%    use of the Doerfler marking with respect to the residual-based error
%    estimator. RHO = 1 means that (essentially) all elements are marked
%    for refinement, whereas small RHO leads to highly adapted meshes.
%
%    The function returns the adaptively generated mesh in terms of 
%    the matrices coordinates and elements, the vector x of the nodal
%    values of the P1-FEM solution as well the corresponding
%    refinement indicators, i.e., the value of the error estimator is
%    given by sqrt(sum(indicators)).
%
%    If known, the exact solution u_ex can be handed over as well to
%    compute the error directly.
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

data = []; H1 = NaN; L2 = NaN; H1s = NaN; Linf = NaN; k=1;
while 1
  %*** Compute discrete solution
  [x,energy] =  solveLaplace(coordinates,elements,constrains, ...
                               dirichlet,neumann,f,g,uD);
  if (nargin == 11)   
    [H1s,H1,L2,Linf] = estimateError(coordinates,elements,zeros(0,4),u_ex,x);
  end
  %*** Compute refinement indicators
  indicators = computeEtaR(x,coordinates,elements,constrains, ...
                           dirichlet,neumann,f,g);
  data = [data;size(coordinates,1),sqrt(sum(indicators)),energy,H1s,H1,L2,Linf];
  %*** Stopping criterion
  if size(elements,1) >= nEmax
      break
  end
  %*** Mark elements for refinement
  [indicators,idx] = sort(indicators,'descend');
  sumeta = cumsum(indicators);
  ell = find(sumeta>=sumeta(end)*rho,1);
  marked = idx(1:ell);
   
  %*** Refine mesh
  [coordinates,elements,constrains,dirichlet,neumann] ...
      = TrefineR(coordinates,elements,constrains,dirichlet,neumann,marked);
end
