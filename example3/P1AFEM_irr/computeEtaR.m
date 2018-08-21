function etaR = computeEtaR(x,coordinates,elements,constrains, ...
                            dirichlet,neumann,f,g)
% Residual Error Estimator
% compute 
%         eta_R =  int_T        h_T^2 | f + div(A grad u_h ... |^2 dT 
%                + int_(bdry T) h_E   [ [ J(u_h) ]]^2  dE
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

nE = size(elements,1);
[edge2nodes,element2edges,foo{1},dirichlet2edges,neumann2edges] ...
     = provideGeometricData([elements;constrains],zeros(0,4),dirichlet,neumann);
%*** First vertex of elements and corresponding edge vectors
c1  = coordinates(elements(:,1),:);
d21 = coordinates(elements(:,2),:) - c1;
d31 = coordinates(elements(:,3),:) - c1;
%*** Vector of element volumes 2*|T|
area2 = d21(:,1).*d31(:,2)-d21(:,2).*d31(:,1);
%*** Compute curl(uh) = (-duh/dy, duh/dx)
u21 = repmat(x(elements(:,2))-x(elements(:,1)), 1,2);
u31 = repmat(x(elements(:,3))-x(elements(:,1)), 1,2);
curl = (d31.*u21 - d21.*u31)./repmat(area2,1,2);
%*** Compute edge terms hE*(duh/dn) for uh
dudn21 =  sum(d21.*curl,2);
dudn13 = -sum(d31.*curl,2);
dudn32 =  -(dudn13+dudn21);
%element2edges(1:nE,:)
etaR = accumarray(reshape(element2edges(1:nE,:),[],1),...
                 [dudn21;dudn32;dudn13],[size(edge2nodes,1),1]);
%*** Incorporate Neumann data
if ~isempty(neumann)
  cn1 = coordinates(neumann(:,1),:);
  cn2 = coordinates(neumann(:,2),:);
  gmE = feval(g,(cn1+cn2)/2);
  etaR(neumann2edges) = etaR(neumann2edges)-sqrt(sum((cn2-cn1).^2,2)).*gmE;
end
%*** Incorporate constrains
if ~isempty(constrains)
  etaR(element2edges(nE+1:end,2)) = etaR(element2edges(nE+1:end,2))...
    + etaR(element2edges(nE+1:end,1))/2;
  etaR(element2edges(nE+1:end,3)) = etaR(element2edges(nE+1:end,3))...
    + etaR(element2edges(nE+1:end,1))/2;
%   etaR(element2edges(nE+1:end,1)) = ...
%                          sqrt(2*(etaR(element2edges(nE+1:end,2))...
%                                 +etaR(element2edges(nE+1:end,3))).^2);
  etaR(element2edges(nE+1:end,1)) = ...
                         sqrt(2*(etaR(element2edges(nE+1:end,2)).^2 ...
                                +etaR(element2edges(nE+1:end,3)).^2));
end
%*** Incorporate Dirichlet data
etaR(dirichlet2edges) = 0;
%*** Assemble edge contributions of indicators
etaR = sum(etaR(element2edges(1:size(elements,1),:)).^2,2);
%*** Add volume residual to indicators
fsT = feval(f,(c1+(d21+d31)/3));
etaR = etaR + (0.5*area2.*fsT).^2;
