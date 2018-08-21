function [coordinates,newElements,marked4,irregular]...
    =recoarseedges_tri(coordinates,elements3,elements4,marked3,marked4)

%recoarseedges_tri: coarsening of regularity edges to obtain an irregular 
%                   quadrilateral mesh  
%
%Usage:
%
% [coordinates,elements4,marked4,irregular] = 
%           recoarseedges_tri(coordinates,elements3,elements4,marked3,marked4)
%
%Comments:
%
%    recoarseedges_tri expects as input a mesh described by the 
%    fields coordinates, elements3, elements4, marked3 and marked4.
%    It recoarses green refinements to obtain an irregular mesh.
%    The vector marked3 contains the indices of triangles
%    which have to be refined, and marked4 contains the indices of 
%    quadrilaterals, respectively. Marked elements that are coarsened pass
%    their marking to their father elements.
%    The function returns the irregular mesh obtained by coarsening the green
%    refinements and translation of the marking. The 
%    output parameters are coordinates, elements4, marked4, irregular.
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

if nnz(find(elements4(:,1)...
        ~= min(elements4,[],2))) || ~isempty(elements3)
    %*** Determine type of refinement for quadrilaterals
    nE=length(elements4(:,1));
    rdx=find(elements4(:,1)~=min(elements4,[],2));
    gdx3=rdx:2:nE;
    %*** Determine type of refinement for triangles
    nG = length(elements3);
    gdx2 = nG-3;
    while gdx2>0
        if elements3(gdx2,3) == elements3(gdx2+1,3) && ...
                elements3(gdx2+1,3) == elements3(gdx2+2,3)
            gdx2 = gdx2 -4;
            if gdx2 ==1
                break;
            end
        else
            gdx2 = gdx2+4;
            break;
        end
    end
    if gdx2 ==0
        gdx2 = 4:4:nG; 
        gdx1 = 1;
    elseif gdx2==-3 
        gdx2 = []; 
        gdx1 = []; 
    else
        gdx1 = 1:3:gdx2-1;
        gdx2 = gdx2: 4: nG;
    end
    %*** Initial triangulation is not sorted
    if isempty(rdx)
        rdx = nE+1;
    end
    %*** Generate red elements
    newElements=elements4(1:rdx(1)-1,:);
    irregular=[];
    % Mark father element for marked green elements
    markedelements3=zeros(length(elements3(:,1)),1);
    markedelements4=zeros(length(elements4(:,1)),1);
    markedelements3(marked3)=1;
    markedelements4(marked4)=1;
    newmarks=markedelements4(1:rdx(1)-1);
    if ~isempty(gdx1)
        marks1 = ...
        max(reshape(markedelements3(gdx1(1):gdx1(end)+2)',3,[]))';
        newmarks=[newmarks;marks1];
    end
    if ~isempty(gdx2)
        marks2 = ...
            max(reshape(markedelements3(gdx2(1):gdx2(end)+3)',4,[]))';
        newmarks=[newmarks;marks2];
    end
    if ~isempty(gdx3)
        marks3 = ...
            max(reshape(markedelements4(gdx3(1):gdx3(end)+1)',2,[]))';
        newmarks=[newmarks;marks3];
    end
    marked4=find(newmarks==1);
    %*** Recoarse green elements and generate irregularity data
    if ~isempty(gdx1)
        tmp_green1=reshape(elements3(gdx1(1):gdx1(end)+2,:)',9,[])';
        newElements = [newElements;tmp_green1(:,[1,7,4,3])];
        irregular = [irregular;tmp_green1(:,[1,7,2])];
    end
    if ~isempty(gdx2)
        tmp_green2 =reshape(elements3(gdx2(1):gdx2(end)+3,:)',12,[])';
        newElements = [newElements;tmp_green2(:,[11,8,3,1])];
        irregular = [irregular;tmp_green2(:,[11,8,5]);...
            tmp_green2(:,[1,11,2])];
    end
    if ~isempty(gdx3)
        tmpgreen3=reshape(elements4(gdx3(1):gdx3(end)+1,:)',8,[])';
        newElements = [newElements;tmpgreen3(:,[7,2,3,6])];
        irregular = [irregular;tmpgreen3(:,[7,2,1]);...
            tmpgreen3(:,[3,6,4])];
    end
    %*** Ensure z1 < z2 for irregularity data
    change_ind=min(irregular,[],2)==(irregular(:,2));
    irregular(change_ind,:)=irregular(change_ind,[2,1,3]);
    %*** Ensure z1 < min(z2,z3,z3) for all elements
    A=[1,2,3,4;4,1,2,3;3,4,1,2;2,3,4,1];
    for i = 1:4
        change_ind=min(newElements,[],2)==newElements(:,i);
        newElements(change_ind,A(i,:))=newElements(change_ind,:);
    end
else
    newElements=elements4;
    irregular=zeros(0,3);
end
end

