function [coordinates,newElements,varargout] ...
    = TrefineRG(coordinates,elements,varargin)

%TrefineRG: local refinement of triangular mesh by red-green refinement, 
%           where marked elements are refined by bisecting all edges 
%           of the element
%
%Usage:
%
% [coordinates,elements3,dirichlet,neumann] ...
%    = TrefineRG(coordinates,elements3,dirichlet,neumann,marked)
% or
%
% [coordinates,elements3] ...
%    = TrefineRG(coordinates,elements3,marked)
%
%Comments:
%
%    TrefineRG expects as input a mesh described by the 
%    fields coordinates, elements3, dirichlet (optional) and neumann 
%    (optional). The vector marked contains the indices of elements which
%    are refined by refining all edges of the element.
%    Further elements will be refined by a red-green refinement to obtain
%    a regular triangulation. To ensure shape regularity of the mesh,
%    green elements are coarsened before further refined.
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

persistent nG
nE = size(elements,1);
markedElements = varargin{end};
%*** Obtain geometric information on edges
[edge2nodes,element2edges,~,boundary2edges{1:nargin-3}] ...
    = provideGeometricData(elements,zeros(0,4),varargin{1:end-1});
%*** Count number of green sibling elements;
if isempty(nG)
    nG = 0;
end
nR = nE-nG;
%*** Mark edges for refinement
edge2newNode = zeros(1,size(edge2nodes,1));
marked = markedElements(markedElements<=nR);
edge2newNode(element2edges(marked,:)) = 1;
marked = ceil((markedElements(markedElements>nR)-nR)/2);
edge2newNode(element2edges(nR+[2*marked-1,nE+2*marked])) = 1;
hashR = logical([1,1,1;1,0,0;0,1,0;0,0,1]);
[mapR,valueR] = hash2map((0:7)',hashR); 
hashG = logical([1,0,0,1;1,1,0,1;1,0,1,1;1,1,1,1]);
[mapG,valueG] = hash2map((0:15)',hashG); 
swap = 1;
while ~isempty(swap) || any(flags(:))
    markedEdge = edge2newNode(element2edges);
    %*** Change flags for red elements
    bit = markedEdge(1:nR,:);
    dec = sum(bit.*(ones(nR,1)*2.^(0:2)),2);
    valR = valueR(dec+1);
    [idx,jdx] = find(~bit & mapR(dec+1,:));
    swap = idx +(jdx-1)*nE;
    edge2newNode(element2edges(swap)) = 1;
    %*** Change flags for green elements
    bit = [markedEdge(nR+1:2:end,1:2), markedEdge(nR+2:2:end,1:2)];
    dec = sum(bit.*(ones(nG/2,1)*2.^(0:3)),2);
    valG = valueG(dec+1);
    gdx = find(valG)';
    flags = ~bit & mapG(dec+1,:);
    edge2newNode(element2edges(nR+2*gdx(flags(gdx,1))-1,1)) = 1;
    edge2newNode(element2edges(nR+2*gdx(flags(gdx,2))-1,2)) = 1;
    edge2newNode(element2edges(nR+2*gdx(flags(gdx,3)),1)) = 1;
    edge2newNode(element2edges(nR+2*gdx(flags(gdx,4)),2)) = 1;
end
%*** Generate new nodes
edge2newNode(edge2newNode~=0) = size(coordinates,1) ...
                                    + (1:nnz(edge2newNode));
idx = find(edge2newNode);
coordinates(edge2newNode(idx),:) = (coordinates(edge2nodes(idx,1),:)...
                            +coordinates(edge2nodes(idx,2),:))/2;
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
%*** Determine type of refinement for each red element
none = find(valR == 0);
r2red = find(valR == 1);
r2green1 = find(valR == 2);
r2green2 = find(valR == 3);
r2green3 = find(valR == 4);
%*** Determine type of refinement for each green element
g2green = nR + find(valG == 0);
g2red = nR + find(valG == 1);
g2red1 = nR + find(valG == 2);
g2red2 = nR + find(valG == 3);
g2red12 = nR + find(valG == 4);
%*** Generate element numbering for refined mesh
rdx = zeros(nR+nG/2,1);
rdx(none)    = 1;
rdx([r2red,g2red])   = 4;
rdx([g2red1,g2red2])  = 3;
rdx(g2red12) = 2;
rdx = [1;1+cumsum(rdx)];
gdx = zeros(size(rdx));
gdx([r2green1, r2green2, r2green3, g2green, g2red1, g2red2]) = 2;
gdx(g2red12) = 4;
gdx = rdx(end)+[0;0+cumsum(gdx)];
%*** Generate new red elements
newElements = 1+zeros(gdx(end)-1,3);
newElements(rdx(none),:) = elements(none,:);
tmp = [elements(1:nR,:),newNodes(1:nR,:),zeros(nR,2);...
    elements(nR+1:2:end,2),elements(nR+2:2:end,[2,3,1])...
    newNodes(nR+2:2:end,2), newNodes(nR+1:2:end,1:2),...
    newNodes(nR+2:2:end,1)];
newElements([rdx(r2red),1+rdx(r2red),2+rdx(r2red),3+rdx(r2red)],:) ...
    = [tmp(r2red,[4,5,6]);tmp(r2red,[1,4,6]);tmp(r2red,[2,5,4]);...
    tmp(r2red,[3,6,5])];
newElements([rdx(g2red),1+rdx(g2red),2+rdx(g2red),3+rdx(g2red)],:) ...
    = [tmp(g2red,[4,5,6]);tmp(g2red,[1,4,6]);...
    tmp(g2red,[4,2,5]);tmp(g2red,[3,6,5])];
newElements([rdx(g2red1),1+rdx(g2red1),2+rdx(g2red1)],:) ...
    = [tmp(g2red1,[4,5,6]);tmp(g2red1,[4,2,5]);tmp(g2red1,[3,6,5])];
newElements([rdx(g2red2),1+rdx(g2red2),2+rdx(g2red2)],:) ...
    = [tmp(g2red2,[4,5,6]);tmp(g2red2,[1,4,6]);tmp(g2red2,[3,6,5])];
newElements([rdx(g2red12),1+rdx(g2red12)],:) ...
    = [tmp(g2red12,[4,5,6]);tmp(g2red12,[3,6,5])];
%*** New green elements
newElements([gdx(r2green1),1+gdx(r2green1)],:) ...
    = [tmp(r2green1,[3,1,4]);tmp(r2green1,[4,2,3])];
newElements([gdx(r2green2),1+gdx(r2green2)],:) ...
    = [tmp(r2green2,[1,2,5]);tmp(r2green2,[5,3,1])];
newElements([gdx(r2green3),1+gdx(r2green3)],:) ...
    = [tmp(r2green3,[2,3,6]);tmp(r2green3,[6,1,2])];
newElements([gdx(g2green),1+gdx(g2green)],:) ...
    = [tmp(g2green,[3,1,4]);tmp(g2green,[4,2,3]);];
newElements([gdx(g2red1),1+gdx(g2red1)],:) ...
    = [tmp(g2red1,[6,1,7]);tmp(g2red1,[7,4,6])];
newElements([gdx(g2red2),1+gdx(g2red2)],:) ...
    = [tmp(g2red2,[5,4,8]);tmp(g2red2,[8,2,5])];
newElements([gdx(g2red12),1+gdx(g2red12),2+gdx(g2red12),...
    3+gdx(g2red12)],:) = [tmp(g2red12,[6,1,7]);tmp(g2red12,[7,4,6]);...
    tmp(g2red12,[5,4,8]);tmp(g2red12,[8,2,5])];
nG = size(newElements,1)-rdx(end)+1;
