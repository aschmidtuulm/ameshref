function [mark3,mark4] = markCircle(coordinates,elements3,elements4,C,R,h)
%markCircle:Marks all triangular or quadrilateral elements
%           with diameter > h which intersect the circle with
%           radius R and center C.
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
%    S. Funken, A. Schmidt  21-08-18

if ~isempty(elements3)
  hmax  = zeros(size(elements3,1),1);
  mark3 = zeros(size(elements3,1),1);
  for j=1:3
    k = mod(j,3)+1;
    mark3 = max(mark3,intersectSegmentCircle(...
                       coordinates(elements3(:,j),:),...
                       coordinates(elements3(:,k),:),C,R));
    hmax = max(hmax,sqrt(sum((coordinates(elements3(:,j),:) ...
                             -coordinates(elements3(:,k),:)).^2,2)));
  end  
  mark3 = mark3 & (hmax >= h);
else
  mark3 = [];
end

if ~isempty(elements4)
  hmax  = zeros(size(elements4,1),1);
  mark4 = zeros(size(elements4,1),1);
  for j=1:4
    k = mod(j,4)+1;
    mark4 = max(mark4,intersectSegmentCircle(...
                       coordinates(elements4(:,j),:), ...
                       coordinates(elements4(:,k),:),C,R));
    hmax = max(hmax,sqrt(sum((coordinates(elements4(:,j),:) ...
                             -coordinates(elements4(:,k),:)).^2,2)));
  end  
  hmax = max(hmax,sqrt(sum((coordinates(elements4(:,1),:) ...
                           -coordinates(elements4(:,3),:)).^2,2)));
  hmax = max(hmax,sqrt(sum((coordinates(elements4(:,2),:) ...
                           -coordinates(elements4(:,4),:)).^2,2)));
  mark4 = mark4 & (hmax >= h);
else
  mark4 = [];
end