function boolean_value =  PointInTriangle2(P,elements,coordinates)
%PointInTriangle2: computes if a point is in a triangle
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

% Compute vectors 
A = coordinates(elements(:,1),:)';
B = coordinates(elements(:,2),:)';
C = coordinates(elements(:,3),:)';

n_p = size(P,2);
m_e = size(elements,1);

i_e = mod(0:n_p*m_e-1,m_e)+1;
i_p = ((1:n_p*m_e) - i_e)/m_e + 1;

v0 = C(:,i_e) - A(:,i_e);
v1 = B(:,i_e) - A(:,i_e);
v2 = P(:,i_p) - A(:,i_e);

% Compute dot products
dot00 = dot(v0, v0);
dot01 = dot(v0, v1);
dot02 = dot(v0, v2);
dot11 = dot(v1, v1);
dot12 = dot(v1, v2);

% Compute barycentric coordinates
invDenom = 1 ./ (dot00 .* dot11 - dot01 .* dot01);
u = (dot11 .* dot02 - dot01 .* dot12) .* invDenom;
v = (dot00 .* dot12 - dot01 .* dot02) .* invDenom;

boolean_value = zeros(1,m_e*n_p);
boolean_value( (u >= 0) & (v >= 0) & (u + v <= 1)) = 1;
boolean_value = reshape(boolean_value,m_e,n_p)';
end

