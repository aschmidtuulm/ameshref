%test_refinement: refine given mesh along a circle with different
%                 mesh refinement strategies
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
clear all
format compact
clf

addpath ../refinement
%% Define initial mesh and circle
c = [0,0;1,0;1,1;0,1;2,0;2,1];
e3 = [3,1,2;1,3,4;2,6,3;6,2,5];
e4 = [1,2,3,4;2,5,6,3];

C = [0.8,0.7]; 
R = 0.5;
h = 1e-2; 

%% *** Do TrefineNVB
coordinates = c; elements3 = e3;
patch('Faces', elements3, 'Vertices', coordinates, 'Facecolor','none')
print -depsc2 'MeshEmptyT'

while 1
  mark3 = markCircle(coordinates,elements3,[],C,R,h);
  mark3 = find(mark3);
  [coordinates,elements3] = TrefineNVB(coordinates,elements3,mark3);
  if isempty(mark3)
    break
  end
end
patch('Faces', elements3, 'Vertices', coordinates, 'Facecolor','none')
print -depsc2 'MeshTrefineNVB'
pause


%% *** Do TrefineRGB
coordinates = c; elements3 = e3; 

while 1
  mark3 = markCircle(coordinates,elements3,[],C,R,h);
  mark3 = find(mark3);
  [coordinates,elements3] = TrefineRGB(coordinates,elements3,mark3);
  if isempty(mark3)
    break
  end
end
patch('Faces', elements3, 'Vertices', coordinates, 'Facecolor','none')
print -depsc2 'MeshTrefineRGB'
pause

%% *** Do TrefineRG
coordinates = c; elements3 = e3; 

while 1
  mark3 = markCircle(coordinates,elements3,[],C,R,h);
  mark3 = find(mark3);
  [coordinates,elements3] = TrefineRG(coordinates,elements3,mark3);
  if isempty(mark3)
    break
  end
end
patch('Faces', elements3, 'Vertices', coordinates, 'Facecolor','none')
print -depsc2 'MeshTrefineRG'
pause

%% *** Do TrefineR
coordinates = c; elements3 = e3; constrains = zeros(0,3);

while 1
  mark3 = markCircle(coordinates,elements3,[],C,R,h);
  mark3 = find(mark3);
   [coordinates,elements3,constrains] = ...
      TrefineR(coordinates,elements3,constrains,mark3);
   if isempty(mark3)
    break
  end
end
patch('Faces', elements3, 'Vertices', coordinates, 'Facecolor','none')
print -depsc2 'MeshTrefineR'
pause
 
%% *** Do QrefineR
coordinates = c; elements4 = e4; irregular = zeros(0,3);
patch('Faces', elements4, 'Vertices', coordinates, 'Facecolor','none')
print -depsc2 'MeshEmptyQ'
while 1
  [~,mark4] = markCircle(coordinates,[],elements4,C,R,h);
  mark4 = find(mark4);
  [coordinates,elements4,irregular] = QrefineR(coordinates,elements4,irregular,mark4);
  if isempty(mark4)
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
  [~,mark4] = markCircle(coordinates,[],elements4,C,R,h);
  mark4 = find(mark4);
  [coordinates,elements4] = QrefineRB(coordinates,elements4,mark4);
  if isempty(mark4)
    break
  end
end
clf
patch('Faces', elements4, 'Vertices', coordinates, 'Facecolor','none')
print -depsc2 'MeshQrefineRGquad'
pause

%% *** Do QrefineR2
coordinates = c; elements4 = e4;irregular = zeros(0,4);
while 1
  [~,mark4] = markCircle(coordinates,[],elements4,C,R,h);
  mark4 = find(mark4);
  [coordinates,elements4,irregular] = QrefineR2(coordinates,elements4,irregular,mark4);
  if isempty(mark4)
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
  [~,mark4] = markCircle(coordinates,[],elements4,C,R,h);
  mark4 = find(mark4)';
  %*** Refine mesh
  [coordinates,elements4,marked,irregular]...
      = recoarseedges(coordinates,elements4,mark4);
  [coordinates,elements4,irregular] ...
      = QrefineR2(coordinates,elements4,irregular,marked);
  [coordinates,elements4]...
      =regularizeedges(coordinates,elements4,irregular);

  clf
  if isempty(mark4)
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
  [mark3,mark4] = markCircle(coordinates,elements3,elements4,C,R,h);
  mark4 = find(mark4)';
  mark3 = find(mark3)';
  %*** Refine mesh
  [coordinates,elements4,marked,elements3]...
      =recoarseedges_tri(coordinates,elements3,elements4,mark3,mark4);
  
  [coordinates,elements4,irregular] ...
      = QrefineR(coordinates,elements4,elements3,marked);
  
  [coordinates,elements4,elements3]...
      =regularizeedges_tri(coordinates,elements4,irregular);

  clf
  if isempty(mark4)
    break
  end
end
clf
patch('Faces', elements4, 'Vertices', coordinates, 'Facecolor','none')
hold on;
patch('Faces', elements3, 'Vertices', coordinates, 'Facecolor','none')
print -depsc2 'MeshQrefineRGtri2'
