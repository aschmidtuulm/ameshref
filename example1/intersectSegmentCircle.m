function flag = intersectSegmentCircle(A,B,X,R)
%intersectSegmentCircle:Flag is one if the segment AB (including endpoints) and 
%                      the circle with radius R and origin X have a common point,
%                      otherwise flag is zero. It works for multiple segments and the same
%                      circle. The first and second coordinate of each point
%                      should be given by A(:,1) and A(:,2) resp. B(:,1) and B(:,2).
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

X = X(:)';
nP = size(A,1);
BmA = B-A;
ApBm2X = A+B-2*X(ones(nP,1),:);
a = dot(BmA,BmA,2);
b = dot(BmA,ApBm2X,2);
c = dot(ApBm2X,ApBm2X,2)-4*R^2;
b2mac = b.^2-a.*c;
SQR = sqrt(b2mac);
SQR(b2mac<0) = realmax;
lambda = [-b+SQR,-b-SQR]./a(:,[1,1]);
flag = min(abs(lambda),[],2)<=1;



