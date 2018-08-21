function [H1s,H1,L2,Linf] = estimateError(coordinates,elements3,elements4,u,x)
% estimates Error

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
if ~isempty(elements4) 
  idx = nC+(1:size(elements4,1))';
  coordinates = [coordinates;[sum(reshape(coordinates(elements4',1),4,[]))'/4, ...
                              sum(reshape(coordinates(elements4',2),4,[]))'/4]];
  x = [x(1:nC);sum(reshape(x(elements4',1),4,[]))'/4];     
  elements3 = [elements3;elements4(:,[1,2]),idx;elements4(:,[2,3]),idx; ...
               elements4(:,[3,4]),idx;elements4(:,[4,1]),idx];    
end
eC = x - u(coordinates);
err3 = [reshape(eC(elements3),[],3),(x(elements3(:,[2,3,1]))+x(elements3))/2 - ...
       reshape(u((coordinates(elements3(:,[2,3,1]),:)+coordinates(elements3,:))/2),[],3)];
Linf = max(abs(err3));     
%*** First vertex of elements and corresponding edge vectors
c1  = coordinates(elements3(:,1),:);
d21 = coordinates(elements3(:,2),:) - c1;
d31 = coordinates(elements3(:,3),:) - c1;
%*** Vector of element volumes 2*|T|
area4 = 2*(d21(:,1).*d31(:,2)-d21(:,2).*d31(:,1));
%***
a = (sum(d21.*d31,2)./area4);
b = (sum(d31.*d31,2)./area4);
c = (sum(d21.*d21,2)./area4);
Sa = -[6 1 1 -4 0 -4; 1 0 -1 -4 4 0; 1 -1 0 0 4 -4; ...
      -4 -4 0 8 -8 8; 0 4 4 -8 8 -8; -4 0 -4 8 -8 8]/3;
Sb =  [3 1 0 -4 0  0; 1 3  0 -4 0 0; 0  0 0 0 0  0; ...
      -4 -4 0 8  0 0; 0 0 0  0 8 -8;  0 0  0 0 -8 8]/3;
Sc =  [3 0 1  0 0 -4; 0 0  0  0 0 0; 1  0 3 0 0 -4;  ...
       0  0 0 8 -8 0; 0 0 0 -8 8  0; -4 0 -4 0  0 8]/3;
H1s = zeros(size(elements3,1),1);
for i = 1:6
  for j = 1:6
    H1s = H1s + (Sa(i,j)*a+Sb(i,j)*b+Sc(i,j)*c).*err3(:,i).*err3(:,j);
  end
end
H1 = sum(H1s);
H1s = sqrt(H1);
if nargout<2
  return
end

Mloc = [6, -1, -1,  0, -4,  0; -1, 6, -1,  0,  0, -4; -1, -1, 6, -4,  0,  0; ...
        0,  0, -4, 32, 16, 16; -4, 0,  0, 16, 32, 16;  0, -4, 0, 16, 16, 32]/720;
L2 = sum(sum( (err3*Mloc).*err3,2).*area4);
H1 = sqrt(H1+L2);
L2 = sqrt(L2);
