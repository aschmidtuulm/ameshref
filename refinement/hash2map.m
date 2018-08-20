function [map,val] = hash2map(dec,hash)
% hash2map: maps markings given in decimal number to given hashs
%
%Usage:
% [map,val] = hash2map(dec,hash)
%
%Comments:
%
%    hash2map expects as input a decimal number which determines the marking 
%    uniquely (e.g. 110 corresponds to dec = 3) and a hash where all possible
%    refinement patterns are described. The function maps the marking to one 
%    of the edges by finding the pattern for which at least further edges have
%    to be marked.
%
%    The function returns the mapping and an assigned val for the hash.
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

n = size(hash,2);
bin=rem(floor(dec*pow2(1-n:0)),2); 
bin = fliplr(bin);
map = zeros(size(bin));
val = zeros(1,size(bin,1));
[idx,jdx] = find(bin);
for i = 1:size(bin,1)
    if dec(i)
        kdx = idx==i;
        % already marked edges stay marked edges
        [pos,~] = find(hash(:,jdx(kdx)));
        % find corresponding hash
        [~,mdx] = min(sum(abs(hash(pos,:)...
                              -bin(i*ones(length(pos),1),:)),2));
        map(i,:)= hash(pos(mdx),:);
        val(i) = pos(mdx);
    end
end
end

