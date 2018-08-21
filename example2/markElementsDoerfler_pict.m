function marked = markElementsDoerfler_pict(eta,theta)
%markElementsDoerfler_pict: marks elements via Doerfler marking
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

[S,I] = sort(eta,'descend');
i = find(cumsum(S) > theta*sum(eta),1);
marked = I(1:i);

end