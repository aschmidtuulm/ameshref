function etaR = computeEtaR(x,coordinates,elements,constrains, ...
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
nT = size(elements,1);
[edge2nodes,constrains2edges,element2edges,dirichlet2edges,neumann2edges] ...
  =provideGeometricData(constrains,elements,dirichlet,neumann);
nE = size(edge2nodes,1);
etaR = zeros(nE,1);
%*** Compute jump term on Neumann boundaries
if ~isempty(neumann)
  cn1 = coordinates(neumann(:,1),:);
  cn2 = coordinates(neumann(:,2),:);
  gmE = feval(g,(cn1+cn2)/2);
  etaR(neumann2edges) = etaR(neumann2edges) - sqrt(sum((cn2-cn1).^2,2)).*gmE;
end
%*** Compute jump term on quadrilateral elements
%[bary,wg] = quad1d(8); 
bary = eye(2); wg = [1;1]/2;
nQ = size(bary,1);
etaR = etaR(:,ones(1,nQ));
xx = reshape(coordinates(elements',1),4,[]);
yy = reshape(coordinates(elements',2),4,[]);
dphi_ds = repmat([-1,1,0,0],nQ,1);
dphi_dt = [-bary(:,2),-bary(:,1),bary(:,1),bary(:,2)];
%*** loop over edges
s = [2,3,4,1];
for j = 1:4
  t = coordinates(elements(:,s(j)),:)-coordinates(elements(:,j),:);
  ind = elements(:,j) > elements(:,s(j));
  x2 = reshape(x(elements),[],4);
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
%  etaR(element2edges(:,j)) = etaR(element2edges(:,j)) + du_dn*wg;
  etaR(element2edges(ind,j),:) = ...
                    etaR(element2edges(ind,j),:) + du_dn(ind,:);
  etaR(element2edges(~ind,j),end:-1:1) = ...
                    etaR(element2edges(~ind,j),end:-1:1) + du_dn(~ind,:);            

  temp = dphi_dt;
  dphi_dt(:,s) = dphi_ds;
  dphi_ds(:,s) = -temp;
end
%*** Incorporate constrains
if ~isempty(constrains)
  tmp = zeros(size(constrains));
  
  ind = constrains(:,1) < constrains(:,2);
  tmp( ind,[1,3]) = etaR(constrains2edges( ind,1),:)/2;
  tmp(~ind,[3,1]) = etaR(constrains2edges(~ind,1),:)/2;
  tmp(:,2) = sum(tmp(:,[1,3]),2)/2;
  
  ind = constrains(:,2) < constrains(:,3);
  tmp( ind,[3,2]) = tmp( ind,[3,2]) + etaR(constrains2edges( ind,2),:);
  tmp(~ind,[2,3]) = tmp(~ind,[2,3]) + etaR(constrains2edges(~ind,2),:);
  
  ind = constrains(:,3) < constrains(:,1);
  tmp( ind,1) = tmp( ind,1) + etaR(constrains2edges( ind,3),2);
  tmp(~ind,1) = tmp(~ind,1) + etaR(constrains2edges(~ind,3),1);

  etaR(constrains2edges(:,1),:) = 0;
  etaR(constrains2edges(:,2),:) = tmp(:,[2,3]);
  etaR(constrains2edges(:,3),:) = tmp(:,[1,2]);

end
%*** Incorporate Dirichlet data
etaR(dirichlet2edges,:) = 0;
etaR = 1/4*(etaR(:,1)+etaR(:,2)).^2+1/12*(etaR(:,2)-etaR(:,1)).^2;
%*** Incorporate constrains
if ~isempty(constrains)
  etaR(element2edges(nT+1:end,1)) = 2*sum(reshape(etaR(element2edges(nT+1:end,2:3)),[],2),2);
end
etaR = sum(reshape(etaR(element2edges(1:nT,:)),[],4),2);

%*** Add volume residual to indicators
[w,s] = gauss(4);
[xg,yg] = meshgrid(s,s); wg = w * w'; 
xg = xg(:); yg = yg(:); wg = wg(:); nQ = size(wg,1);
%*** Compute volume term
phi  = [(1-xg).*(1-yg),(1+xg).*(1-yg),(1+xg).*(1+yg),(1-xg).*(1+yg)]/4;
xx = reshape(coordinates(elements',1),4,[]);
yy = reshape(coordinates(elements',2),4,[]);
etaR_V = zeros(size(elements,1),1);
for k = 1:nQ
  dphi_ds = [yg(k)-1,  1-yg(k), 1+yg(k), -(1+yg(k))]/4;
  dphi_dt = [xg(k)-1,-(1+xg(k)),1+xg(k),   1-xg(k) ]/4;
  d2phi_dsdt = [1,-1,1,-1]/4;
  dF_ds = [dphi_ds * xx; dphi_ds * yy]';
  dF_dt = [dphi_dt * xx; dphi_dt * yy]';
  d2F_dsdt = [d2phi_dsdt * xx; d2phi_dsdt * yy]';
  det_DF = sum((dF_ds *[0,1;-1,0 ]).*dF_dt,2);
  du_ds = x(elements) * dphi_ds';
  du_dt = x(elements) * dphi_dt';
  d2u_dsdt = x(elements) * d2phi_dsdt';
  aa = [dF_dt(:,2).*du_ds-dF_ds(:,2).*du_dt, ...
       -dF_dt(:,1).*du_ds+dF_ds(:,1).*du_dt];
  bb = -2*sum(dF_dt.*dF_ds,2)./det_DF.^2;      
  D2u = sum(bb.*(d2u_dsdt-sum(d2F_dsdt.*aa,2)./det_DF),2);
  etaR_V = etaR_V + det_DF.*wg(k) .* (D2u+f([phi(k,:)*xx;phi(k,:)*yy]')).^2;
end
etaR = etaR + sum((coordinates(elements(:,2),:)...
                  -coordinates(elements(:,1),:)).^2,2) .* etaR_V;


