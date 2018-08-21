%mainAdaptiveScript: generates a triangulation of a picture via a
%                    gradient based error estimator using Dörfler marking
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

clear all, clc,close all

addpath ../refinement
%*** load picture
pic = imread('./pics/muenster.jpg');
pic = rgb2gray(pic);
[m,n] = size(pic);

%*** plot picture
figure(1);
imshow(pic);
pic = pic(m:-1:1,:);

%*** initial coordinates, elements, boundary
c = [1,1;n,1;n,m;1,m];
e3 = [2,4,1;4,2,3];
e4 = [1,2,3,4];
b = [1,2;2,3;3,4;4,1];
nit = 7; % number of iterations

%% *** Do TrefineNVB
coordinates = c; elements = e3; boundary = b;
for i =1:nit
    fprintf('******** TrefineNVB ******** STEP  %d **********\n',i);
    %*** compute refinement indicator
    etaR = computeEtaR_pict(pic,elements,coordinates);
    %*** mark elements
    marked = markElementsDoerfler_pict(etaR,0.8);
    %*** refine mesh
    [coordinates,elements,boundary]= TrefineNVB(coordinates,elements,boundary,marked);
    clf
    patch('Faces', elements, 'Vertices', coordinates, 'Facecolor','none')
    str = ['MeshTrefineNVB' num2str(i)];
    print(figure(1),str,'-depsc2')
end

%% *** Do TrefineRGB
coordinates = c; elements = e3; boundary = b;
for i =1:nit
    fprintf('******** TrefineRGB ******** STEP  %d **********\n',i);
    %*** compute refinement indicator
    etaR = computeEtaR_pict(pic,elements,coordinates);
    %*** mark elements
    marked = markElementsDoerfler_pict(etaR,0.8);
    %*** refine mesh
    [coordinates,elements,boundary]= TrefineRGB(coordinates,elements,boundary,marked);
    clf
    patch('Faces', elements, 'Vertices', coordinates, 'Facecolor','none')
    str = ['MeshTrefineRGB' num2str(i)];
    print(figure(1),str,'-depsc2')
end

%% *** Do TrefineRG
coordinates = c; elements = e3; boundary = b;
for i =1:nit
    fprintf('******** TrefineRG ******** STEP  %d **********\n',i);
    %*** compute refinement indicator
    etaR = computeEtaR_pict(pic,elements,coordinates);
    %*** mark elements
    marked = markElementsDoerfler_pict(etaR,0.8);
    %*** refine mesh
    [coordinates,elements,boundary]= TrefineRG(coordinates,elements,boundary,marked);
    clf
    patch('Faces', elements, 'Vertices', coordinates, 'Facecolor','none')
    str = ['MeshTrefineRG' num2str(i)];
    print(figure(1),str,'-depsc2')
end


%% *** Do TrefineR
coordinates = c; elements = e3; boundary = b; irregular = zeros(0,3);
%*** refine mesh
for i =1:nit
    fprintf('******** TrefineR ******** STEP  %d **********\n',i);
    %*** compute refinement indicator
    etaR = computeEtaR_pict(pic,elements,coordinates);
    %*** mark elements
    marked = markElementsDoerfler_pict(etaR,0.8);
    %*** refine mesh
    [coordinates,elements,irregular,boundary]= TrefineR(coordinates,elements,irregular,boundary,marked);
    clf
    patch('Faces', elements, 'Vertices', coordinates, 'Facecolor','none')
    str = ['MeshTrefineR' num2str(i)];
    print(figure(1),str,'-depsc2')
end

%% *** Do QrefineR
coordinates = c; elements = e4; boundary = b; irregular = zeros(0,3);
%*** refine mesh
for i =1:nit
    fprintf('******** QrefineR ******** STEP  %d **********\n',i);
    %*** compute refinement indicator
    etaR = computeEtaR_pict(pic,elements,coordinates);
    %*** mark elements
    marked = markElementsDoerfler_pict(etaR,0.8);
    %*** refine mesh
    [coordinates,elements,irregular,boundary]= QrefineR(coordinates,elements,irregular,boundary,marked);
    clf
    patch('Faces', elements, 'Vertices', coordinates, 'Facecolor','none')
    str =  ['MeshQrefineR' num2str(i)];
    print(figure(1),str,'-depsc2')
end

%% *** Do QrefineRB
coordinates = c; elements = e4; boundary = b;
%*** refine mesh
for i =1:nit
    fprintf('******** QrefineRB ******** STEP  %d **********\n',i);
    %*** compute refinement indicator
    etaR = computeEtaR_pict(pic,elements,coordinates);
    %*** mark elements
    marked = markElementsDoerfler_pict(etaR,0.8);
    %*** refine mesh
    [coordinates,elements,boundary]= QrefineRB(coordinates,elements,boundary,marked);
    clf
    patch('Faces', elements, 'Vertices', coordinates, 'Facecolor','none')
    str = ['MeshQrefineRGquad' num2str(i)];
    print(figure(1),str,'-depsc2')
end


%% *** Do QrefineR2
coordinates = c; elements = e4; boundary = b; irregular = zeros(0,4);
%*** refine mesh
for i =1:nit-2
    fprintf('******** QrefineR2 ******** STEP  %d **********\n',i);
    %*** compute refinement indicator
    etaR = computeEtaR_pict(pic,elements,coordinates);
    %*** mark elements
    marked = markElementsDoerfler_pict(etaR,0.8);
    %*** refine mesh
    [coordinates,elements,irregular,boundary]= QrefineR2(coordinates,elements,irregular,boundary,marked);
    clf
    patch('Faces', elements, 'Vertices', coordinates, 'Facecolor','none')
    str = ['MeshQrefineR2' num2str(i)];
    print(figure(1),str,'-depsc2')
end

%% *** Do QrefineRG
coordinates = c; elements4 = e4; elements3 = zeros(0,3); boundary = b;
%*** refine mesh
for i =1:nit
    fprintf('******** QrefineRG ******** STEP  %d **********\n',i);
    %*** compute refinement indicator
    etaR3 = computeEtaR_pict(pic,elements3,coordinates);
    etaR4 = computeEtaR_pict(pic,elements4,coordinates);
    %*** mark elements
    mark3 = markElementsDoerfler_pict(etaR3,0.8);
    mark4 = markElementsDoerfler_pict(etaR4,0.8);
    %*** refine mesh
    [coordinates,elements4,marked,elements3]...
        =recoarseedges_tri(coordinates,elements3,elements4,mark3,mark4);
    
    [coordinates,elements4,irregular] ...
        = QrefineR(coordinates,elements4,elements3,marked);
    
    [coordinates,elements4,elements3]...
        =regularizeedges_tri(coordinates,elements4,irregular);
    clf
    patch('Faces', elements4, 'Vertices', coordinates, 'Facecolor','none')
    hold on;
    patch('Faces', elements3, 'Vertices', coordinates, 'Facecolor','none')
    str = ['MeshQrefineRG2' num2str(i)];
    print(figure(1),str,'-depsc2')
end


%% *** Do QrefineRG2
coordinates = c; elements = e4; boundary = b;
%*** refine mesh
for i =1:nit-2
    fprintf('******** QrefineRG2 ******** STEP  %d **********\n',i);
    %*** compute refinement indicator
    etaR = computeEtaR_pict(pic,elements,coordinates);
    %*** mark elements
    marked = markElementsDoerfler_pict(etaR,0.8);
    %*** refine mesh
    [coordinates,elements,marked,irregular]...
        = recoarseedges(coordinates,elements,marked);
    [coordinates,elements,irregular] ...
        = QrefineR2(coordinates,elements,irregular,marked);
    [coordinates,elements]...
        =regularizeedges(coordinates,elements,irregular);
    clf
    patch('Faces', elements, 'Vertices', coordinates, 'Facecolor','none')
    str = ['MeshQrefineRG2' num2str(i)];
    print(figure(1),str,'-depsc2')
end
