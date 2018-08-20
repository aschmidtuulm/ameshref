function [coordinates,newElements] ...
    = regularizeedges(coordinates,elements,irregular)

%regularizeedges: regularization of irregular quadrilateral mesh  
%            where three different green refinement types come into use
%
%Usage:
%
% [coordinates,elements4] = regularizeedges(coordinates,elements4,irregular)
%
%Comments:
%
%    regularizeedges expects as input a mesh described by the 
%    fields coordinates, elements4, and irregular (virtual quadrilaterals).
%    It regularizes irregular edges by matching appropriate green
%    refinements.
%    The function returns the irregular mesh as regular mesh with output
%    parameters coordinates, elements4.
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
nC = size(coordinates,1);
nI = size(irregular,1);
%*** provide geometric data
[~,~,element2edges] ...
    = provideGeometricData(zeros(0,3),[elements;irregular]);
irregular2edges = element2edges(end-(nI)+1:end,1);
element2edges = element2edges(1:end-nI,:);
edge2element = zeros(nE,2);
for k = 1:4
    edge2element(element2edges(:,k),1) = 1:nE;
    edge2element(element2edges(:,k),2) = k;
end
%*** Provide irregularity data
IrregularEdges = zeros(nE,4,2);
for i = 1:length(irregular(:,1))
     id=edge2element(irregular2edges(i),:);
     IrregularEdges(id(:,1),id(:,2),:)=irregular(i,3:4);
 end
%*** Check orientation of hanging nodes
orientation_irr=elements(:,[2,3,4,1])>elements;
changeOr = [logical(zeros(size(orientation_irr))),orientation_irr];
tmp = IrregularEdges;
IrregularEdges(orientation_irr)=IrregularEdges(changeOr);
IrregularEdges(changeOr)=tmp(orientation_irr);
%*** Determine type of refinement
markedEdges=IrregularEdges(:,:,1);
markedEdges(markedEdges~=0)=1;
marks = sum(markedEdges,2);
none = marks==0;
green1 = marks==1;
green2 = min((markedEdges(:,1:2)+markedEdges(:,3:4)),[],2)==1;
green3=max((markedEdges(:,1:2)+markedEdges(:,3:4)),[],2)==2;
%*** Determine orientation of refinement patterns
orientation_green2=markedEdges(green2,[4,1,2,3])*[0,1,1,4]';
orientation_green2(orientation_green2==5)=3;
[orientation_green3,~,~]=find(markedEdges(green3,1:2)'==1);
[orientation_green1,~,~]=find(markedEdges(green1(:,1),:)'==1);
% standardize elements such that z1 < min{z2,z3,z4} for regular elements
elements_tmp=elements(none,:);
ind=elements_tmp(:,2)==min(elements_tmp,[],2);
elements_tmp(ind,:)=elements_tmp(ind,[2,3,4,1]);
newElements = [elements_tmp;zeros(length(orientation_green1)*4 ...
    +length(orientation_green2)*5+length(orientation_green3)*3,4)];
%*** Generate new interior nodes for green1 and green2 elements
midnodes1=nC+reshape(1:2*length(orientation_green1),[],2);
midnodes2=nC+2*max(size(midnodes1(:,1),1))+...
    reshape(1:2*length(orientation_green2),[],2);
%*** Rotate data until orientation of refinement patterns is at 1
elements_gr1=elements(green1,:);
elements_gr2=elements(green2,:);
elements_gr3=elements(green3,:);
edge2newNodegr2_1=IrregularEdges(green2(:,1),[4,1,2,3],1);
edge2newNodegr2_2=IrregularEdges(green2(:,1),[4,1,2,3],2);
edge2newNodegr3_1=IrregularEdges(green3(:,1),1:4,1);
edge2newNodegr3_2=IrregularEdges(green3(:,1),1:4,2);
A = [1,2,3,4;4,1,2,3;3,4,1,2;2,3,4,1];
c = [4,2,1,2;2,4,2,1;1,2,4,2;2,1,2,4];
for i =2:4
    idx=orientation_green1==i;
    jdx=orientation_green2==i;
    elements_gr1(idx,A(i,:))=elements_gr1(idx,:);
    elements_gr2(jdx,A(i,:))=elements_gr2(jdx,:);
    edge2newNodegr2_1(jdx,A(i,:))=edge2newNodegr2_1(jdx,:);
    edge2newNodegr2_2(jdx,A(i,:))=edge2newNodegr2_2(jdx,:);
end
jdx=orientation_green3==2;
elements_gr3(jdx,A(2,:))=elements_gr3(jdx,:);
edge2newNodegr3_1(jdx,A(2,:))=edge2newNodegr3_1(jdx,:);
edge2newNodegr3_2(jdx,A(2,:))=edge2newNodegr3_2(jdx,:);
%*** Generate new nodes and elements
tmp_green1=[elements_gr1,sum(IrregularEdges(green1(:,1),:,1),2),...
    sum(IrregularEdges(green1(:,1),:,2),2),midnodes1];
tmp_green2=[elements_gr2,edge2newNodegr2_1(:,1:2),...
    edge2newNodegr2_2(:,1:2),midnodes2];
tmp_green3=[elements_gr3,edge2newNodegr3_1(:,[1,3]),...
    edge2newNodegr3_2(:,[1,3])];
if ~isempty(tmp_green1)
    coordinates=[coordinates;...
        (c(1,3)*coordinates((tmp_green1(:,1)),:) ...
        + c(2,3)*coordinates((tmp_green1(:,2)),:) ...
        + c(3,3)*coordinates((tmp_green1(:,3)),:) ...
        + c(4,3)*coordinates((tmp_green1(:,4)),:))/9];
    coordinates=[coordinates;...
        (c(1,4)*coordinates((tmp_green1(:,1)),:) ...
        + c(2,4)*coordinates((tmp_green1(:,2)),:) ...
        + c(3,4)*coordinates((tmp_green1(:,3)),:) ...
        + c(4,4)*coordinates((tmp_green1(:,4)),:))/9];
    gdx1=nnz(newElements(:,1))+(1:4:(length(orientation_green1)*4));
    newElements([gdx1,gdx1+1,gdx1+2,gdx1+3],:)=...
        [tmp_green1(:,[7,8,5,6]);tmp_green1(:,[7,6,2,3]);...
        tmp_green1(:,[8,7,3,4]);tmp_green1(:,[5,8,4,1])];    
end
if ~isempty(tmp_green2)
    coordinates=[coordinates;...
        (c(1,1)*coordinates((tmp_green2(:,1)),:) ...
        + c(2,1)*coordinates((tmp_green2(:,2)),:) ...
        + c(3,1)*coordinates((tmp_green2(:,3)),:) ...
        + c(4,1)*coordinates((tmp_green2(:,4)),:))/9];
    coordinates=[coordinates;...
        (c(1,3)*coordinates((tmp_green2(:,1)),:) ...
        + c(2,3)*coordinates((tmp_green2(:,2)),:) ...
        + c(3,3)*coordinates((tmp_green2(:,3)),:) ...
        + c(4,3)*coordinates((tmp_green2(:,4)),:))/9];
    gdx2=nnz(newElements(:,1))+(1:5:(length(orientation_green2)*5));   
    newElements([gdx2,gdx2+1,gdx2+2,gdx2+3,gdx2+4],:)=...
        [tmp_green2(:,[10,9,6,8]);tmp_green2(:,[10,8,2,3]);...
        tmp_green2(:,[5,10,3,4]);tmp_green2(:,[10,5,7,9]);...
        tmp_green2(:,[1,6,9,7])];
end
if ~isempty(tmp_green3)
    gdx3=nnz(newElements(:,1))+(1:3:(length(orientation_green3)*3));
    newElements([gdx3,gdx3+1,gdx3+2],:)=...
        [tmp_green3(:,[6,7,2,3]);tmp_green3(:,[5,8,4,1]);...
        tmp_green3(:,[6,8,5,7])];
end
