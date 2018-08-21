function [lambda, weights] = quad1d(order)
%*** Table of quadrature formulae of different order
npoints = ceil((order+1)/2);
switch npoints
  case 1
    lambda = 0.5;
    weights = 1;
  case 2
    lambda =  2.1132486540518712e-1;
    weights = 0.5;
  case 3
    lambda = [1.1270166537925831e-1;0.5];
    weights = [5/18;4/9];
  case 4
    lambda = [6.9431844202973712e-2; 3.3000947820757187e-1];
    weights = [1.7392742256872693e-1; 3.2607257743127307e-1];
  case 5
    lambda = [4.6910077030668004e-2; 2.3076534494715845e-1; 0.5];
    weights = [1.1846344252809454e-1; 2.3931433524968323e-1; 128/450];
  otherwise
      error(['Gauss rule of order ',num2str(order),' not supported.'])  
end
%*** Generate quadrature points and weights from data above
weights(npoints:-1:npoints-size(lambda,1)+1,1) = weights;
lambda(npoints:-1:npoints-size(lambda,1)+1,1) = 1-lambda;
lambda = [lambda(end:-1:1),lambda];
