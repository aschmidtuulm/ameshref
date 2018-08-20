function [coordinates,newElements,newIrregular,varargout] ...
    = QrefineR2(coordinates,elements,irregular,varargin)

%QrefineR2: local refinement of quadrilateral mesh by red2 refinement, 
%          where marked elements are refined by bisecting all edges 
%          of the element
%
%Usage:
%
% [coordinates,elements4,irregular,dirichlet,neumann] ...
%    = QrefineR2(coordinates,elements4,irregular,dirichlet,neumann,marked)
% or
%
% [coordinates,elements4,irregular] ...
%    = QrefineR2(coordinates,elements4,irregular,marked)
%
%Comments:
%
%    QrefineR2 expects as input a mesh described by the 
%    fields coordinates, elements4, irregular, dirichlet (optional) and neumann 
%    (optional). The vector marked contains the indices of elements which
%    are refined by refining all edges of the element.
%    Further elements will be refined by a red2 refinement to obtain. 
%    2-Irregularity of the mesh is ensured by the 2-Irregular Rule.
% 
%    The function returns the refined mesh in terms of the same data as
%    for the input.
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
%    S. Funken, A. Schmidt  20-08-18


nE = size(elements,1);
markedElements = varargin{end};
%*** Obtain geometric information on edges
[edge2nodes,~,element2edges,boundary2edges{1:nargin-4}] ...
    = provideGeometricData(zeros(0,3),[elements;irregular],...
    varargin{1:end-1});
irregular2edges = element2edges(nE+1:end,:);
element2edges = element2edges(1:nE,:);
orientation = elements(:,[2,3,4,1])<elements;
%*** Mark edges for refinement and existing hanging nodes
edge2newNode = zeros(1,size(edge2nodes,1));
edge2newNode(element2edges(markedElements,:)) = 1;
edge2newNode(irregular2edges(:,1)) = 1;
kdx = 1;
while ~isempty(kdx) || ~isempty(swap)
    markedEdge = edge2newNode(element2edges);
    %*** Change flags for elements
    kdx = find(sum(abs(markedEdge),2)<4 & ...
        (sum(abs(markedEdge),2)>2 | min(markedEdge,[],2)<0));
    [idx,jdx] = find(~markedEdge(kdx,:));
    edge2newNode(element2edges( kdx(idx)+(jdx-1)*nE)) = 1;
    %*** Change flags for irregular marker elements
    markedEdge = edge2newNode(irregular2edges);
    flag = irregular2edges(any(markedEdge(:,2:end),2),1);
    swap = find(edge2newNode(flag)~=-1);
    edge2newNode(flag(swap)) = -1;
end
%*** Generate new nodes on edges
edge2newNode = edge2newNode([1,1],:)';
edge2newNode(irregular2edges(:,1),:) = -1;
idx = edge2newNode(:,1)>0;
edge2newNode(idx,:) = size(coordinates,1)+reshape(1:2*nnz(idx),2,[])';
coordinates(edge2newNode(idx,1),:) ...
    = (2*coordinates(edge2nodes(idx,1),:)+...
    coordinates(edge2nodes(idx,2),:))/3;
coordinates(edge2newNode(idx,2),:) ...
    = (coordinates(edge2nodes(idx,1),:)+...
    2*coordinates(edge2nodes(idx,2),:))/3;
%*** Refine boundary conditions
varargout = cell(nargout-3,1);
for j = 1:nargout-3
    boundary = varargin{j};
    if ~isempty(boundary)
        newNodes = edge2newNode(boundary2edges{j},:);
        markedEdges = find(newNodes(:,1));
        if ~isempty(markedEdges)
            ind = boundary(markedEdges,1) > boundary(markedEdges,2);
            newNodes(markedEdges(ind),:) ...
                = newNodes(markedEdges(ind),[2,1]);
            boundary = [boundary(~newNodes(:,1),:); ...
                boundary(markedEdges,1),newNodes(markedEdges,1); ...
                newNodes(markedEdges,1),newNodes(markedEdges,2); ...
                newNodes(markedEdges,2),boundary(markedEdges,2)];
        end
    end
    varargout{j} = boundary;
end
%*** Provide new nodes for refinement of elements
edge2newNode(irregular2edges(:,1),:) = irregular(:,[4,3]);
newNodes = reshape(edge2newNode(element2edges,:),[],8);
%*** Determine type of refinement for each element
reftyp = (newNodes(:,1:4)~=0)*2.^(0:3)';
none   = reftyp < 15;
red2    = reftyp == 15;
%*** Generate new irregularity data
kdx = find(reftyp >0 & reftyp < 15);
if ~isempty(kdx)
    [idx,jdx] = find( newNodes(kdx,1:4));
    edx = element2edges(kdx(idx)+(jdx-1)*nE);  
    newIrregular = [edge2nodes(edx(:),:),...
        newNodes(kdx(idx(:))+(jdx(:)+3)*nE),...
        newNodes(kdx(idx(:))+(jdx(:)-1)*nE)];
else
    newIrregular = zeros(0,4);
end
newEdgeNodes = reshape(edge2newNode(irregular2edges(:,2:4),1),[],3);
kdx = find(sum(newEdgeNodes,2)~=0);
[idx,jdx,val] = find( newEdgeNodes(kdx,:));
edx = irregular2edges(kdx(idx)+(jdx-1+1)*size(irregular2edges,1));
newIrregular = [newIrregular;[edge2nodes(edx(:),:),val(:)+1,val(:)]];
%*** Generate new interior nodes if red elements are refined
sdx = find(red2);
midNodes = zeros(nE,4);
midNodes(sdx,:) = size(coordinates,1)+reshape(1:4*length(sdx),[],4);
%*** Generate element numbering for refined mesh
idx = zeros(nE,1);
idx(none)    = 1;
idx(red2)     = 9;
idx = [1;1+cumsum(idx)];
%*** Incorporate orientation of edges
[idx0,jdx0] = find(orientation);
tmp = newNodes(idx0 + (jdx0-1)*nE);
newNodes(idx0 + (jdx0-1)*nE) = newNodes(idx0 + (jdx0+3)*nE);
newNodes(idx0 + (jdx0+3)*nE) = tmp;
%*** Generate new elements
newElements = zeros(idx(end)-1,4);
newElements(idx(none),:) = elements(none,:);
s = [2,3,4,1]; p=[8,5,6,7];
c = [4,2,1,2;2,4,2,1;1,2,4,2;2,1,2,4];
for j=1:4
    newElements([2*(j-1)+idx(red2),2*j-1+idx(red2)],:) = ...
        [elements(red2,j),newNodes(red2,j),midNodes(red2,j),...
        newNodes(red2,p(j));newNodes(red2,j),newNodes(red2,j+4),...
        midNodes(red2,s(j)),midNodes(red2,j)];
    coordinates = [coordinates;...
        (c(1,j)*coordinates(elements(sdx,1),:)...
        + c(2,j)*coordinates(elements(sdx,2),:)...
        + c(3,j)*coordinates(elements(sdx,3),:)...
        + c(4,j)*coordinates(elements(sdx,4),:))/9];
end
newElements(8+idx(red2),:) = midNodes(red2,:);
             


