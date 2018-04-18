# Feedback Vertex Set

Comparison of parametrized algorithms for solving feedback vertex set problem. (Finding minimum set of vertices making the graph acyclic)

The research was inspired by results of PACE Challange 2016 (Parametrized algorithms and computational experiments).

Branching algorithms compared during the research (and their subversions depending on different reductions used):

- O*(8^k) by Yixin Cao (Can be found f.e. http://arxiv.org/abs/1707.08684)
- O*(5^k)  - (Can be f.e http://www.ii.uib.no/~yngvev/publications/konf/ChenFLLV07.pdf)
- O*((2+fi)^k) - By Kociumaka and Pilipczuk (Can be found f.e. https://www.mimuw.edu.pl/~kociumaka/files/ipl2.pdf)
- Branching based on Half Integral Relaxation by Imanishi and Iwata algorithm from Pace (https://arxiv.org/abs/1608.01463 , https://arxiv.org/abs/1310.2841)


To make things clear - O* notation is used in parametrized complexity to express complexity in terms of some parameter k, which 
is the expected size of solution. It disregards the parts of complexity that do not contain this parameter or (in most cases) are not exponential.


### Prerequisites

For compilation purposes we use CMake system, you will need to install it in order to compile.
Definetely works under versions >= 3.7.2, but should also work under older ones.

### Compilation

```
cd src/
mkdir build
cd build/
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```


### Usage (from /src/ level):

`./build/FeedbackVertexSet <sol_type> [Optimizations]`

### Parameters - sol_type

Possible parameters:

* 8 - Default (if no parameter), Runs Cao solution
* 5 - Iterative Compression solution
* 3.62 - Kociumaka Pilipczuk 
* PD - Cao 8^k solution, branching first on vertices with biggest deg including double neighbours
* PU - Cao 8^k, branching first on vertices with maximum number of neighbours in Undeletable set
* Full - Runs Cao solution with all possible optimizations
* LP - branching based on Iwata paper
* Appx - Runs approximate solution

### Parameters - Optimizations

Possible parameters:

* noCC - doesn't use breaking Graph into Connected components reduction
* noDeg3 - Doesn't use deg3 solver in solution
* deg3 - Uses deg3 solver in solution (parameter only for LP branching)
* IC - Uses iterative compression routine while branching (avaialable for 8, PD, PU)
* Full - uses all of the Optimizations available for chosen solution

### Input and Output

Program takes as an input list of edges (as pairs u v, separated by a whitespace), one pair per row. Vertices names can be any string, with only exception they can't start from # - it is used for commented out lines. 
The format is borrowed from PACE 2016. (To find out more look onto website https://pacechallenge.wordpress.com/pace-2016/track-b-feedback-vertex-set/)

Example Input:

```
a b
b c
c d
d a
```

It corresponds to a cycle graph made of 4 vertices.

On the standard output the solution returns two lines:
- In the first line a single integer N denoting size of the solution found
- In the second line N strings corresponding to the names of vertices that were put in solution

Example output (for the Input above):
```
1
c
```


During the main part of our experiments we have used 100 tests that were open during PACE 2016 + 130 that were hidden (during contest they were used for final scoring).


### Output Verificator (from /src/ level):

`./build/Verifier InFile OutFile [Optimal]`

Used for verifying, whether the output returned by solution is actually a correct Feedback Vertex Set and how does it compare optional optimal file.


## Authors

* **Krzysztof Kiljan** - (https://github.com/karek) 
* **Marcin Pilipczuk** - (https://github.com/malcin541) 

We hereby invite you to visit our project site: http://kernelization-experiments.mimuw.edu.pl/

This repository contains parts of the code by Marcin Pilipczuk from https://bitbucket.org/marcin_pilipczuk/fvs-pace-challenge

This code is licensed under simplified BSD license (see the LICENSE file).
