function [coordinates,newElements,varargout] ...
    = TrefineRGB(coordinates,elements,varargin)

%TrefineRGB: local refinement of triangular mesh by red-green-blue refinement, 
%            where marked elements are refined by bisecting all edges 
%            of the element
%
%Usage:
%
% [coordinates,elements3,dirichlet,neumann] ...
%    = TrefineRGB(coordinates,elements3,dirichlet,neumann,marked)
% or
%
% [coordinates,elements3] ...
%    = TrefineRGB(coordinates,elements3,marked)
%
%Comments:
%
%    TrefineRGB expects as input a mesh described by the 
%    fields coordinates, elements3, dirichlet (optional) and neumann 
%    (optional). The vector marked contains the indices of elements which
%    are refined by refining all edges of the element.
%    Further elements will be refined by a red-green-blue refinement to obtain
%    a regular triangulation. 
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

markedElements = varargin{end};
nE = size(elements,1);
%*** Obtain geometric information on edges
[edge2nodes,element2edges,~,boundary2edges{1:nargin-3}] ...
    = provideGeometricData(elements,zeros(0,4),varargin{1:end-1});
%*** Mark edges for refinement
edge2newNode = zeros(1,max(max(element2edges)));
edge2newNode(element2edges(markedElements,:)) = 1;
hash = [1,0,0;1,1,0;1,0,1;1,1,1];
[map,value] = hash2map((0:7)',hash); 
%*** Change flags for elements
swap = 1;
while ~isempty(swap)
    markedEdge = edge2newNode(element2edges);
    dec = sum(markedEdge(1:nE,:).*(ones(nE,1)*2.^(0:2)),2);
    val = value(dec+1);
    [idx,jdx] = find(~markedEdge(1:nE,:) & map(dec+1,:));
    swap = idx +(jdx-1)*nE;
    edge2newNode(element2edges(swap,1)) = 1;
end
%*** Generate new nodes
edge2newNode(edge2newNode~=0) = size(coordinates,1) + (1:nnz(edge2newNode));
idx = find(edge2newNode);
coordinates(edge2newNode(idx),:) ...
    = (coordinates(edge2nodes(idx,1),:)+coordinates(edge2nodes(idx,2),:))/2;
%*** Refine boundary conditions
varargout = cell(nargout-2,1);
for j = 1:nargout-2
    boundary = varargin{j};
    if ~isempty(boundary)
        newNodes = edge2newNode(boundary2edges{j})';
        markedEdges = find(newNodes);
        if ~isempty(markedEdges)
            boundary = [boundary(~newNodes,:); ...
                boundary(markedEdges,1),newNodes(markedEdges); ...
                newNodes(markedEdges),boundary(markedEdges,2)];
        end
    end
    varargout{j} = boundary;
end
%*** Provide new nodes for refinement of elements
newNodes = edge2newNode(element2edges);
%*** Determine type of refinement for each element
none = find(val==0) ;
green   = find(val==1);
bluer  = find(val==2);
bluel  = find(val==3);
red = find(val==4);
%*** Generate element numbering for refined mesh
idx = ones(nE,1);
idx(none) = 1;
idx(green)   = 2;
idx([bluer, bluel])  = 3;
idx(red) = 4;
idx = [1;1+cumsum(idx)];
%*** Generate new elements
newElements = zeros(idx(end)-1,3);
newElements(idx(none),:) = elements(none,:);
tmp = [elements(:,:), newNodes(:,:)];
newElements([idx(green),1+idx(green)],:) ...
    = [tmp(green,[3,1,4]); tmp(green,[2,3,4])];
newElements([idx(bluer),1+idx(bluer),2+idx(bluer)],:) ...
    = [tmp(bluer,[3,1,4]);tmp(bluer,[4,2,5]);tmp(bluer,[3,4,5])];
newElements([idx(bluel),1+idx(bluel),2+idx(bluel)],:) ...
    = [tmp(bluel,[4,3,6]);tmp(bluel,[1,4,6]);tmp(bluel,[2,3,4])];
newElements([idx(red),1+idx(red),2+idx(red),...
    3+idx(red)],:) = [tmp(red,[1,4,6]);tmp(red,...
    [4,2,5]);tmp(red,[6,5,3]);tmp(red,[5,6,4])];
