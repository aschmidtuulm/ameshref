%example4: refine for a given set of discrete points all elements containing these 
%          points with different mesh refinement strategies
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
%    S. Funken, A. Schmidt 25-04-19

clear all
format compact
clf

addpath ../refinement

%% *** define discrete points (e.g. through a prescribed curve)
x = -pi:pi/100:pi;
y = sin(x); %tan(sin(x)) - sin(tan(x));
points = [x',y'];

%% *** define initial triangulation
c = [-3.5,-3;3.5,-3;3.5,3;-3.5,3];
e3  = [2,4,1;4,2,3];
e4 = [1,2,3,4];
maxE = 1e4; % max number of elements

%% *** Do TrefineNVB
coordinates = c; elements3 = e3;
patch('Faces', elements3, 'Vertices', coordinates, 'Facecolor','none')
hold on
plot(x,y,'or')
print -depsc2 'MeshEmptyT'

while 1
  mark3 = point2element(coordinates,elements3,points);
    [coordinates,elements3] = TrefineNVB(coordinates,elements3,mark3);
  if size(elements3,1)>maxE
    break
  end
end
clf
patch('Faces', elements3, 'Vertices', coordinates, 'Facecolor','none')
print -depsc2 'MeshTrefineNVB'
pause


%% *** Do TrefineRGB
coordinates = c; elements3 = e3; 

while 1
  mark3 = point2element(coordinates,elements3,points);
  [coordinates,elements3] = TrefineRGB(coordinates,elements3,mark3);
  if size(elements3,1)>maxE
    break
  end
end
clf
patch('Faces', elements3, 'Vertices', coordinates, 'Facecolor','none')
print -depsc2 'MeshTrefineRGB'
pause

%% *** Do TrefineRG
coordinates = c; elements3 = e3; 

while 1
mark3 = point2element(coordinates,elements3,points);
  [coordinates,elements3] = TrefineRG(coordinates,elements3,mark3);
  if size(elements3,1)>maxE
    break
  end
end
clf
patch('Faces', elements3, 'Vertices', coordinates, 'Facecolor','none')
print -depsc2 'MeshTrefineRG'
pause

%% *** Do TrefineR
coordinates = c; elements3 = e3; constrains = zeros(0,3);

while 1
    mark3 = point2element(coordinates,elements3,points);
   [coordinates,elements3,constrains] = ...
      TrefineR(coordinates,elements3,constrains,mark3);
   if size(elements3,1)>maxE
    break
  end
end
clf
patch('Faces', elements3, 'Vertices', coordinates, 'Facecolor','none')
print -depsc2 'MeshTrefineR'
pause
 
%% *** Do QrefineR
coordinates = c; elements4 = e4; irregular = zeros(0,3);
patch('Faces', elements4, 'Vertices', coordinates, 'Facecolor','none')
print -depsc2 'MeshEmptyQ'
while 1
  mark4 = unique([point2element(coordinates,elements4(:,1:3),points);point2element(coordinates,elements4(:,[1,3,4]),points)]);
  [coordinates,elements4,irregular] = QrefineR(coordinates,elements4,irregular,mark4);
  if size(elements4,1)>maxE
    break
  end
end
clf
patch('Faces', elements4, 'Vertices', coordinates, 'Facecolor','none')
print -depsc2 'MeshQrefineR'
pause

%% *** Do QrefineRB
coordinates = c; elements4 = e4;
while 1
  mark4 = unique([point2element(coordinates,elements4(:,1:3),points);point2element(coordinates,elements4(:,[1,3,4]),points)]);
  [coordinates,elements4] = QrefineRB(coordinates,elements4,mark4);
  if size(elements4,1)>maxE
    break
  end
end
clf
patch('Faces', elements4, 'Vertices', coordinates, 'Facecolor','none')
print -depsc2 'MeshQrefineRB'
pause

%% *** Do QrefineR2
coordinates = c; elements4 = e4;irregular = zeros(0,4);
while 1
  mark4 = unique([point2element(coordinates,elements4(:,1:3),points);point2element(coordinates,elements4(:,[1,3,4]),points)]);
  [coordinates,elements4,irregular] = QrefineR2(coordinates,elements4,irregular,mark4);
  if size(elements4,1)>maxE
    break
  end
end
clf
patch('Faces', elements4, 'Vertices', coordinates, 'Facecolor','none')
print -depsc2 'MeshQrefineR2'
pause

%% *** Do QrefineRG2
coordinates = c; elements4 = e4;

while 1
  mark4 = unique([point2element(coordinates,elements4(:,1:3),points);point2element(coordinates,elements4(:,[1,3,4]),points)]);
  %*** Refine mesh
  [coordinates,elements4,marked,irregular]...
      = recoarseedges(coordinates,elements4,mark4);
  [coordinates,elements4,irregular] ...
      = QrefineR2(coordinates,elements4,irregular,marked);
  [coordinates,elements4]...
      =regularizeedges(coordinates,elements4,irregular);

  clf
  if size(elements4,1)>maxE
    break
  end
end
clf
patch('Faces', elements4, 'Vertices', coordinates, 'Facecolor','none')
print -depsc2 'MeshQrefineRG2'
pause

%% *** Do QrefineRG
coordinates = c; elements3 = zeros(0,3); elements4 = e4;
while 1
  mark4 = unique([point2element(coordinates,elements4(:,1:3),points);point2element(coordinates,elements4(:,[1,3,4]),points)]);
  mark3 = point2element(coordinates,elements3,points);
  %*** Refine mesh
  [coordinates,elements4,marked,elements3]...
      =recoarseedges_tri(coordinates,elements3,elements4,mark3,mark4);
  
  [coordinates,elements4,irregular] ...
      = QrefineR(coordinates,elements4,elements3,marked);
  
  [coordinates,elements4,elements3]...
      =regularizeedges_tri(coordinates,elements4,irregular);

  clf
  if size(elements4,1)+size(elements3,1)>maxE
    break
  end
end
clf
patch('Faces', elements4, 'Vertices', coordinates, 'Facecolor','none')
hold on;
patch('Faces', elements3, 'Vertices', coordinates, 'Facecolor','none')
print -depsc2 'MeshQrefineRGtri2'
