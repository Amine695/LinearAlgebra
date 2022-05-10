# Linear Algebra project

C project with the implementations of several algebraic algorithmic methods such as : <br/>

- The LU decomposition for a square matrix <br/>
- Solving a linear system using the LU decomposition <br/>
- The "naive" matrix inversion using the LU decomposition <br/>
- The inverse matrix using the Strassen's inversion algorithm with naive product <br/>
- Compute the inverse matrix but with the Strassen's multiplication algorithm <br/>
- A benchmark to compute time measurements over the different algorithms. 



## Compilation
A makefile is available to make the compilation easier. To use it, please run `make` in the project folder. A **build** folder will be created with all the objects files inside as well as a **benchmarks** folder that will contains all the data needed for benchmarks.

## Execution 
### Algorithmes
To run all the algorithms, please run `./project --n --p` with **n** the size of the matrix (power of 2) and **p** a prime number.

### Benchmark
To run the benchmarks, please run `./project benchmarks n` with **n** the limit size in power of 2, so for exemple n = 11 would run with matrices of size 2^1 to 2^11. <br/>
Once executed, a menu prints out asking you which algorithm you want to compare.<br/>
Finally, at the end, a **plot** of the results is automatically generated and should appear on the screen.<br/> However, if nothing happens or an error is display, use the command `make plot` directly.<br/>
Compile with **Linux** to avoid problems. If you use **python** instead of **python3** in the temrinal, kindly change the python interpreter's path at the top of `plot.py`. <br/>.

## Documentation
A **Doxygen** documentation is available with clear explanations. Open the `index.html` which is in the *doc* folder.

Please refer to the  **compilation** section of the report for more details.
<br/><br/>

### Authors: <br/>
Amine Berbagui <br/>
Ghassen Hachani