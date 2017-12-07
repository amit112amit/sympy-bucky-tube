# sympy-bucky-tube
Some symbolic calculations to create a flat triangular lattice and roll it up into a zig-zag bucky tube

Assuming a unit lattice spacing, we create a planar triangular lattice in the xy plane. We then roll it up into a cylinder.
The coordinates of the points in the plane and the cylinder are calculated symbolically. There are two functions that lambdify
the symbolic calculations and generate Legacy VTK files that can be viewed in Paraview for inspection of the algorithm.
