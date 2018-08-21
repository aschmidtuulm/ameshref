function [etaR3,etaR4] = computeEtaR(x,coordinates,elements3,elements4, ...
                            dirichlet,neumann,f,g)
% Residual Error Estimator
% compute 
%         eta_R =  int_T        h_T^2 | f + div(A grad u_h ... |^2 dT 
%                + int_(bdry T) h_E   [ [ J(u_h) ]]^2  dE
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

%*** Obtain geometric information on edges
[edge2nodes,element3edges,element4edges,dirichlet2edges,neumann2edges] ...
             = provideGeometricData(elements3,elements4,dirichlet,neumann);
nE = size(edge2nodes,1);
%*** quadrature rule on edges
[bary,wg] = quad1d(8); 
nQ = size(bary,1);
%*** Compute jump term on triangles
%    First vertex of elements and corresponding edge vectors
c1  = coordinates(elements3(:,1),:);
d21 = coordinates(elements3(:,2),:) - c1;
d31 = coordinates(elements3(:,3),:) - c1;
%*** Vector of element volumes 2*|T|
area2 = d21(:,1).*d31(:,2)-d21(:,2).*d31(:,1);
%*** Compute curl(uh) = (-duh/dy, duh/dx)
u21 = repmat(x(elements3(:,2))-x(elements3(:,1)), 1,2);
u31 = repmat(x(elements3(:,3))-x(elements3(:,1)), 1,2);
curl = (d31.*u21 - d21.*u31)./repmat(area2,1,2);
%*** Compute edge terms hE*(duh/dn) for uh
dudn21 =  sum(d21.*curl,2);  
dudn13 = -sum(d31.*curl,2);  
dudn32 =  -(dudn13+dudn21);  
etaR = accumarray(element3edges(:),[dudn21;dudn32;dudn13],[nE 1]);
%eta=0*etaR;
%*** Incorporate Neumann data
if ~isempty(neumann)
  cn1 = coordinates(neumann(:,1),:);
  cn2 = coordinates(neumann(:,2),:);
  gmE = feval(g,(cn1+cn2)/2);
  etaR(neumann2edges) = etaR(neumann2edges) - sqrt(sum((cn2-cn1).^2,2)).*gmE;
end
%*** Compute jump term on quadrilaterals
etaR = etaR(:,ones(1,nQ));
xx = reshape(coordinates(elements4',1),4,[]);
yy = reshape(coordinates(elements4',2),4,[]);
dphi_ds = repmat([-1,1,0,0],nQ,1);
dphi_dt = [-bary(:,2),-bary(:,1),bary(:,1),bary(:,2)];
%*** loop over edges
s = [2,3,4,1];
for j = 1:4
  t = coordinates(elements4(:,s(j)),:)-coordinates(elements4(:,j),:);
  ind = elements4(:,j) > elements4(:,s(j));
  x2 = reshape(x(elements4),[],4);
  du_ds = x2 * dphi_ds'; 
  du_dt = x2 * dphi_dt'; 
  dF_ds = [dphi_ds * xx; dphi_ds * yy]';
  dF_dt = [dphi_dt * xx; dphi_dt * yy]';
  det_DF = dF_ds(:,1:nQ).*dF_dt(:,nQ+1:end) ...
         - dF_ds(:,nQ+1:end).*dF_dt(:,1:nQ);
  du_dn = (t(:,2*ones(1,nQ)).*( dF_dt(:,nQ+1:end).*du_ds ...
                               -dF_ds(:,nQ+1:end).*du_dt) ...
          -t(:,1*ones(1,nQ)).*(-dF_dt(:,   1:nQ ).*du_ds ...
                               +dF_ds(:,   1:nQ ).*du_dt))./det_DF;
  etaR(element4edges(ind,j),:) = ...
                    etaR(element4edges(ind,j),:) + du_dn(ind,:);
  etaR(element4edges(~ind,j),end:-1:1) = ...
                    etaR(element4edges(~ind,j),end:-1:1) + du_dn(~ind,:);            
  temp = dphi_dt;
  dphi_dt(:,s) = dphi_ds;
  dphi_ds(:,s) = -temp;
end
%*** Incorporate Dirichlet data
etaR(dirichlet2edges,:) = 0;
etaR = etaR.^2*wg;
etaR3 = sum(reshape(etaR(element3edges),[],3),2);
etaR4 = sum(reshape(etaR(element4edges),[],4),2);
%*** Add volume residual to indicators
fsT = feval(f,(c1+(d21+d31)/3));
etaR3 = etaR3 + (0.5*area2.*fsT).^2;
%*** Tensor Gauss quadrature on quadrilaterals
[w,s] = gauss(4);
[xg,yg] = meshgrid(s,s); wg = w * w'; 
xg = xg(:); yg = yg(:); wg = wg(:); nQ = size(wg,1);
%*** Compute volume term
phi  = [(1-xg).*(1-yg),(1+xg).*(1-yg),(1+xg).*(1+yg),(1-xg).*(1+yg)]/4;
xx = reshape(coordinates(elements4',1),4,[]);
yy = reshape(coordinates(elements4',2),4,[]);
etaR_V = zeros(size(elements4,1),1);
for k = 1:nQ
  dphi_ds = [yg(k)-1,  1-yg(k), 1+yg(k), -(1+yg(k))]/4;
  dphi_dt = [xg(k)-1,-(1+xg(k)),1+xg(k),   1-xg(k) ]/4;
  d2phi_dsdt = [1,-1,1,-1]/4;
  dF_ds = [dphi_ds * xx; dphi_ds * yy]';
  dF_dt = [dphi_dt * xx; dphi_dt * yy]';
  d2F_dsdt = [d2phi_dsdt * xx; d2phi_dsdt * yy]';
  det_DF = sum((dF_ds *[0,1;-1,0 ]).*dF_dt,2);
  du_ds = x(elements4) * dphi_ds';
  du_dt = x(elements4) * dphi_dt';
  d2u_dsdt = x(elements4) * d2phi_dsdt';
  aa = [dF_dt(:,2).*du_ds-dF_ds(:,2).*du_dt, ...
       -dF_dt(:,1).*du_ds+dF_ds(:,1).*du_dt];
  bb = -2*sum(dF_dt.*dF_ds,2)./det_DF.^2;      
  D2u = sum(bb.*(d2u_dsdt-sum(d2F_dsdt.*aa,2)./det_DF),2);
  etaR_V = etaR_V + det_DF.*wg(k) .* (D2u+f([phi(k,:)*xx;phi(k,:)*yy]')).^2;
end
etaR4 = etaR4 + sum((coordinates(elements4(:,2),:)...
                    -coordinates(elements4(:,1),:)).^2,2) .* etaR_V;


