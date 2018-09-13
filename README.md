# ameshref 
## Efficient Matlab Implementation of Adaptive Mesh Refinement in 2D

This project provides an efficient implementation of various adaptive mesh refinement strategies in two dimensions. The documentation provides a detailed investigation of the implemented mesh refinement methods and also presents how the methods are realized by utilization of reasonable data structure, use of Matlab built-in functions and vectorization. In the documentation, the focus is set on the deployment of ameshref in the context of the adaptive finite element method. To this end, also codes to realize the adaptive finite element method (AFEM) is additionally provided here. However, our code is kept general such that this package is not restricted to the AFEM context as our other test examples show. This package can be used to adaptively generate meshes in two dimensions in different ways. Moreover, it can also serve educational purposes on how to realize mesh refinement methods in an accessible, short but efficient way. The data structure used is very simple because only the coordinates of the vertices and the element-connectivities are needed. Marking of elements are interpreted as edge-based marking, i.e., an element is marked by marking each edge for bisection. In contrast to other mesh refinement packages, we do not make use of a recursive approach but prefer the realization in terms of vectorization.

## Getting Started

For a use of the refinement methods download the complete repository _/ameshref_ with the test examples and run them on your computer with Matlab. To use it within your own examples you only need the repository _/refinement_. Please see the following instructions for a correct use.

### Prerequisites

MATLAB

## Running the test examples

/example1: run _test_refinement.m_ (refinement along a circle)

/example2: run _mainAdaptiveScript.m_ (triangulation of a picture)

/example3: run _example1.m_ or _example2.m_ in each of the repositories (mesh refinement in the context of AFEM)


### Deployment for your own examples -  Data structure

To be able to deploy this package within your framework you need the following data strucure of the mesh:

For triangles you need to define _coordinates_, _elements3_, and optionally boundary data _dirichlet_ or _neumann_. 

![alt text](https://github.com/aschmidtuulm/ameshref/blob/master/TriangulationWithQuadrilaterals.png?raw=true)


Similarly, for quadrilaterals you need to define _coordinates, elements4_, and optionally boundary data _dirichlet_ or _neumann_. 

![alt text](https://github.com/aschmidtuulm/ameshref/blob/master/TriangulationWithTriangles.png?raw=true)

For meshes with hanging nodes an additional data vector named irregular is needed, where _irregular(l,1)_ and _irregular(l,2_) are the starting and end point of the lth-irregular edge with hanging node stored in _irregular(l,3)_. 

### How to call the functions

We abbreviate the data structure as _coordinates (C), elements3 (E3), elements4 (E4), irregular (I), dirichlet (D)_, and _neumann (N)_. Furthermore, marked elements are stored with the corresponding index in the variable _marked_.

The different mesh refinement methods are then called by (remember that _D_ and _N_ are optional arguments)

#### TrefineR

```
[C,E3,I,D,N] = TrefineR(C,E3,I,D,N,marked)
```
#### QrefineR

```
[C,E4,I,D,N] = QrefineR(C,E4,I,D,N,marked)
```
#### QrefineR2

```
[C,E4,I,D,N] = QrefineR2(C,E4,I,D,N,marked)
```
#### TrefineNVB

```
[C,E3,D,N] = TrefineNVB(C,E3,D,N,marked)
```
#### TrefineRGB

```
[C,E3,D,N] = TrefineRGB(C,E3,D,N,marked)
```
#### TrefineRG

```
[C,E3,D,N] = TrefineRG(C,E3,D,N,marked)
```
#### QrefineRB

```
[C,E4,D,N] = QrefineRB(C,E4,D,N,marked)
```
#### QrefineRG

3-step call
```
[C,E4,marked,I] = recoarseedges_tri(C,E3,E4,marked3,marked4)
[C,E4,D,N] = QrefineR(C,E4,D,N,marked)
[C,E4,E3] = regularizeedges_tri(C,E4,I)
```
#### QrefineRG2
3-step call
```
[C,E4,marked,I] = recoarseedges(C,E4,marked4)
[C,E4,I,D,N] = QrefineR2(C,E4,I,D,N,marked)
[C,E4] = regularizeedges_tri(C,E4,I)
```

#### Minimal example
The used functions are provided in the repository _/example1_ and _/refinement_:
```
C = [0.8,0.7]; % center of circle
r = 0.5; % radius of circle
min_d = 5e-6; % minimal diameter of an element allowed

%% exemplary use of a quadrilateral mesh refinement strategy
coordinates = [0,0;1,0;1,1;0,1;2,0;2,1];
elements4 = [1,2,3,4;2,5,6,3];
while 1
    mark4 = markCircle(coordinates,[],elements4,C,r,min_d);
    mark4 = find(mark4);  
    [coordinates,elements4] = QrefineRB(coordinates,elements4,mark4);
    if isempty(mark4)
        break
    end
end
```
#### Plot your mesh with Matlab
triangular mesh:
```
patch('Faces',elements3,'Vertices',coordinates,'Facecolor','none')
```
quadrilateral mesh:
```
patch('Faces',elements4,'Vertices',coordinates,'Facecolor','none')
```

## Authors

* **Stefan A. Funken** - **Anja Schmidt** Institute for Numerical Mathematics, Ulm University, Germany

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

