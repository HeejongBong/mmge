##Research Code for Fitting Multiple Gaussian Graphical Model


###Compilation

You need to have [R][R-website] installed before compiling the code.
You have to install R from source on a linux machine.

To compile the C/C++ code, just type `make` in the root directory. To
recompile, first type `make clean`, then type `make`. 



###Example
To run the 
provided example, just excute `example.R` file in R  

~~~ R
source("src/R/example.R")
~~~ 

This returns a list object `sol_path` storing the solution path for convex and non-convex method. Specifically, the returning list `sol_path` consists of two arrays --- `sol_convex` and `sol_nonconvex` --- with `sol_convex` refering to the solution path
 obtained by using $L_1$ penalty and `sol_nonconvex` refering to that using 
 truncated $L_1$ penalty. The array indices are arranged in a way such that 
 `sol_nonconvex[,,j,i]` is refering to the solution corresponding 
 to grid point `(Lambda1.vec[i], Lambda2.vec[j])`.   

[R-website]: http://www.r-project.org/
