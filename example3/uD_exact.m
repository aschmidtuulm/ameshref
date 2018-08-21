function val = uD_exact(x)
%u_ex = r^(2/3)* sin(2/3*alpha)

%Remark:
%
%    This program is a supplement to the paper
%    >> Adaptive Mesh Refinement in 2D - An Efficient Implementation in Matlab <<
%    by S. Funken, and A. Schmidt. The reader should
%    consult that paper for more information.
%
%Authors:
%
%    S. Funken, A.Schmidt  21-08-18

  beta = 3 * pi/4;
  x = x * [cos(beta),-sin(beta); sin(beta), cos(beta)];
  [phi,r] = cart2pol(x(:,1),x(:,2));
  val = r.^(2/3) .* sin( 2/3 * (phi + beta ) );
