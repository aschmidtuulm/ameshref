%point2element: indicate elements containing points
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
%    S. Funken, A. Schmidt 25-04-19

% Example:
%   coordinates = [0,0;1,0;1,1;0,1]; elements = [1,2,3;1,3,4]; 
%   point2element(coordinates,elements,[0.3,0.2;0.2,0.3])
% 

function ind = point2element(coordinates,elements,pts)
% pts is a list of N points and has dimension N x 2
X = reshape(coordinates(elements,1),[],3);
Y = reshape(coordinates(elements,2),[],3);
ind = zeros(size(pts,1),1);
for k = 1 : size(pts,1)
  L = cross(X-pts(k,1),Y-pts(k,2),2);
  idx = find (L(:,1)>=0 & L(:,2)>=0 & L(:,3)>=0);
  if ~isempty(idx)
    ind(k) = idx(1);
  end
end
ind = setdiff(unique(ind),0);