function [coordinates,newElements,marked,irregular]...
    =recoarseedges(coordinates,elements,marked)

%recoarseedges: coarsening of regularity edges to obtain an irregular 
%                   quadrilateral mesh  
%
%Usage:
%
% [coordinates,elements4,marked,irregular] = 
%           recoarseedges(coordinates,elements4,marked4)
%
%Comments:
%
%    recoarseedges expects as input a mesh described by the 
%    fields coordinates, elements4 and marked4.
%    It recoarses green refinements to obtain an irregular mesh.
%    The vector marked4 contains the indices of elements that need to be
%    refined. Marked elements that are coarsened pass
%    their marking to their father elements.
%    The function returns the irregular mesh obtained by coarsening the green
%    refinements and translation of the marking. The output parameters are
%    coordinates, elements4, marked4, irregular.
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

if nnz(elements(:,1)~=min(elements,[],2))>0
    %*** Determine type of refinement
    rdx=find(elements(:,1)~=min(elements,[],2));
    green_elements=elements(rdx(1):end,:);
    tmp=flipud(green_elements(:,2));
    green3=tmp==[tmp(2:end,:);0];
    green3=find(green3((1:length(green3)))~=0); 
    green_elements=green_elements(1:end-length(green3)*3,:);
    green2=find(green_elements(:,1)==min(green_elements,[],2))-5;
    if ~isempty(green2)
        green1=(1:4:green2(1)-1)-1;
    elseif ~isempty(green3)
        green1=(1:4:green3(1)-1-rdx(1))-1;
    else
        green1=(1:4:length(green_elements(:,1))-1)-1;
    end
    gdx1=(green1+rdx(1))';
    nG1=length(green1);
    gdx2=green2+rdx(1);
    nG2=length(green2);
    gdx3=length(elements)-flipud(green3)-1;
    %*** Generate red elements
    newElements=elements(1:rdx(1)-1,:);
    irregular=[];
    %*** Mark father element for marked green elements
    markedelements=zeros(length(elements(:,1)),1);
    markedelements(marked)=1;
    newmarks=markedelements(1:rdx(1)-1);
    if ~isempty(gdx1)
        marks1 = ...
            max(reshape(markedelements(gdx1(1):gdx1(end)+3)',4,[]))';
        newmarks=[newmarks;marks1];
    end
    if ~isempty(gdx2)
        marks2 = ...
            max(reshape(markedelements(gdx2(1):gdx2(end)+4)',5,[]))';
        newmarks=[newmarks;marks2];
    end
    if ~isempty(gdx3)
        marks3 = ...
            max(reshape(markedelements(gdx3(1):gdx3(end)+2)',3,[]))';
        newmarks=[newmarks;marks3];
    end
    marked=find(newmarks==1);
    %*** Recoarse green elements and generate irregularity data
    if ~isempty(gdx1)
        tmp_green1=reshape(elements(gdx1(1):gdx1(end)+3,:)',16,[])';
        newElements=[newElements;tmp_green1(:,[7,8,12,16])];
        irregular=[irregular;tmp_green1(:,[7,16,3,4])];
    end
    if ~isempty(gdx2)
        tmp_green2=reshape(elements(gdx2(1):gdx2(end)+4,:)',20,[])';
        newElements=[newElements;tmp_green2(:,[17,7,8,12])];
        irregular=[irregular;tmp_green2(:,[7,17,3,4]);...
            tmp_green2(:,[12,17,15,14])];
    end
    if ~isempty(gdx3)
        tmp_green3=reshape(elements(gdx3(1):gdx3(end)+2,:)',12,[])';
        newElements=[newElements;tmp_green3(:,[3,4,7,8])];
        irregular=[irregular;tmp_green3(:,[3,8,5,2]);...
            tmp_green3(:,[4,7,10,9])];
    end
    %*** Delete midpoints of green elements
    nMidpoints=nG1*2+nG2*2;
    coordinates=coordinates(1:end-nMidpoints,:);
    %*** Ensure z1 < z2 for irregularity data
    change_ind=min(irregular,[],2)==(irregular(:,2));
    irregular(change_ind,:)=irregular(change_ind,[2,1,4,3]);
    %*** Ensure z1 < min(z2,z3,z3) for all elements
    A=[1,2,3,4;4,1,2,3;3,4,1,2;2,3,4,1];
    for i = 1:4
        change_ind=min(newElements,[],2)==newElements(:,i);
        newElements(change_ind,A(i,:))=newElements(change_ind,:);
    end
else
    newElements=elements;
    irregular=zeros(0,4);
end
end
