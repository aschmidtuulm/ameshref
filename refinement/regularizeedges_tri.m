function [coordinates,newElements4, newElements3] ...
    = regularizeedges_tri(coordinates,elements,irregular)

%regularizeedges_tri: regularization of irregular quadrilateral mesh  
%            where three different green refinement types (triangular and 
%            quadrilateral) come into use
%
%Usage:
%
% [coordinates,elements4,elements3] = 
%           regularizeedges_tri(coordinates,elements4,irregular)
%
%Comments:
%
%    regularizeedges_tri expects as input a mesh described by the 
%    fields coordinates, elements4, and irregular (virtual triangles).
%    It regularizes irregular edges by matching appropriate green
%    refinements.
%    The function returns the irregular mesh as regular mesh consisting 
%    of quadrilaterals and triangles with output parameters coordinates, 
%    elements3, elements4.
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

nE= size(elements,1);
%*** Provide geometric data
[~,irregular2edges,element2edges] ...
    = provideGeometricData(irregular,elements);
edge2element = zeros(nE,2);
for k = 1:4
    edge2element(element2edges(:,k),1) = 1:nE;
    edge2element(element2edges(:,k),2) = k;
end
%*** Provide irregularity data
IrregularEdges=zeros(nE,4);
for i = 1:length(irregular(:,1))
     id=edge2element(irregular2edges(i),:);
     IrregularEdges(id(:,1),id(:,2))=irregular(i,3);
 end
%*** Determine type of refinement
markedEdges=IrregularEdges(:,:);
markedEdges(markedEdges~=0)=1;
marks = sum(markedEdges,2);
none = marks==0;
green1 = marks==1;
green2 = min((markedEdges(:,1:2)+markedEdges(:,3:4)),[],2)==1;
green3 = max((markedEdges(:,1:2)+markedEdges(:,3:4)),[],2)==2;
% Determine orientation of refinement patterns
orientation_green2=markedEdges(green2,[4,1,2,3])*[0,1,1,4]';
orientation_green2(orientation_green2==5)=3;
[orientation_green3,~,~]=find(markedEdges(green3,1:2)'==1);
[orientation_green1,~,~]=find(markedEdges(green1(:,1),:)'==1);
% Standardizing elements such that z1<min{z2,z3,z4} for regular elements
elements_tmp=elements(none,:);
ind=elements_tmp(:,2)==min(elements_tmp,[],2);
elements_tmp(ind,:)=elements_tmp(ind,[2,3,4,1]);
newElements4 = [elements_tmp; zeros(length(orientation_green3)*2,4)];
newElements3 = zeros(length(orientation_green1)*3 ...
    +length(orientation_green2)*4,3);
%*** Rotate data until orientation of refinement pattern is at 1
elements_gr1=elements(green1,:);
elements_gr2=elements(green2,:);
elements_gr3=elements(green3,:);
edge2newNodegr2=IrregularEdges(green2(:,1),[4,1,2,3]);
edge2newNodegr3=IrregularEdges(green3(:,1),1:4);
A = [1,2,3,4;4,1,2,3;3,4,1,2;2,3,4,1];
for i =2:4
    idx=orientation_green1==i;
    jdx=orientation_green2==i;
    elements_gr1(idx,A(i,:))=elements_gr1(idx,:);
    elements_gr2(jdx,A(i,:))=elements_gr2(jdx,:);
    edge2newNodegr2(jdx,A(i,:))=edge2newNodegr2(jdx,:);
end
jdx=orientation_green3==2;
elements_gr3(jdx,A(2,:))=elements_gr3(jdx,:);
edge2newNodegr3(jdx,A(2,:))=edge2newNodegr3(jdx,:);
%*** Generate new nodes and elements
tmp_green1=[elements_gr1,sum(IrregularEdges(green1(:,1),:),2)];
tmp_green2=[elements_gr2,edge2newNodegr2(:,1:2)];
tmp_green3=[elements_gr3,edge2newNodegr3(:,[1,3])];
if ~isempty(tmp_green1)
    gdx1=1:3:length(orientation_green1)*3;
    newElements3([gdx1,1+gdx1,2+gdx1],:)=...
        [tmp_green1(:,[1,5,4]); tmp_green1(:,[3,4,5]);...
        tmp_green1(:,[2,3,5])];
end
if ~isempty(tmp_green2)
    gdx2 = length(orientation_green1)*3+...
        (1:4:length(orientation_green2)*4);
    newElements3([gdx2,1+gdx2,2+gdx2,3+gdx2],:)=...
        [tmp_green2(:,[4,5,3]);tmp_green2(:,[5,6,3]);...
        tmp_green2(:,[6,2,3]);tmp_green2(:,[5,1,6])];
end
if ~isempty(tmp_green3)
    gdx3=nnz(newElements4(:,1))+1:2:length(newElements4);
    newElements4([gdx3,1+gdx3],:)=[tmp_green3(:,[5,2,3,6]);...
        tmp_green3(:,[6,4,1,5])];
end
