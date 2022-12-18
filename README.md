# Fortran90 Examples of Dynamic Systems

## Table of contents

1. [Introduction](#Introduction)
2. [Usage](#Usage)
3. [Examples](#Examples)
	1. [Oscillator with non-linear forcing components](#Oscillator-with-non-linear-forcing-components)
	2. [Logistic map](#Logistic-map)
	3. [Henon map](#Henon-map)
	4. [Lorenz model](#Lorenz-model)
	5. [Multifractal measure](#Multifractal-measure)
4. [Online resources](#Online-resources)
5. [carepackage.f90](#carepackagef90)
6. [System requirements](#System-requirements)

## Introduction 

Fortran90 examples of simple dynamic systems, for an introductory university class. Aside from the code, all examples include a:

1. `Makefile` for ease of cleaning/compilation
2. `*_plotter.plt` (sub the `*` with the program name) Gnuplot script for the plots required by the class assignments 
3. `*_runner.sh` (see above for the `*`) bash helper script to create data folders and run both the compiled example and gnuplot

In addition, examples [3](https://github.com/Squar3wave/fortran90_dynsys/tree/master/3_henon_map) and [4](https://github.com/Squar3wave/fortran90_dynsys/tree/master/4_lorenz_model) include a `parameters.txt` file from which the `*_runner.sh` script takes the needed input values. To use these examples stand-alone, just uncomment the values i left in the code and comment all the `get_command_argument` related lines.  

Feel free to alter the code  to fit your needs, just don't keep it to yourself, fortran90 is already hard enough to use and/or get accustomed to, people need all the help they can get! I sure did.

## Usage

All examples come already compiled and with data generated.  
In order to run them from scratch:

1. Open your terminal of choice and navigate to the example folder
2. Either use the `make clean` comand to remove the compiled example and all data and plots, or do so manually
3. Either use the `make` comand to compile, or do so with preferred/required method of choice
4. Type `./*_runner.sh` to run the helper script (you may need to run `chmod +x *_runner.sh` to make it executable) in order to execute both the example and Gnuplot, or comment the gnuplot line if other plotting/fitting methods are preferred/required.

Important: should you choose not to use the included `Makefile` and/or `*_runner.sh` script, make sure to create the needed data folders and to feed the examples the right input values (like for examples [3](https://github.com/Squar3wave/fortran90_dynsys/tree/master/3_henon_map) and [4](https://github.com/Squar3wave/fortran90_dynsys/tree/master/4_lorenz_model)).


## Examples

### Oscillator with non-linear forcing components

$$\ddot{x} + 2\lambda\epsilon\dot{x} + \omega_0^2x +\epsilon gx^3 = \epsilon f \cos  [(\omega_0 + \epsilon \sigma) \cdot t]$$

1. Integration using [Runge-Kutta](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods) 4th order algorithm in two ways:  
	* Increasing values of $\sigma$. Starting from the second $\sigma$ value onward, for each evolution the starting conditions are the last values of the previous one.  
	* Decreasing values of $\sigma$. Same procedure as above.  
2. Study of maxima behavior in function of $\sigma$ once the forcing effects have subsided for the two kinds of evolution above mentioned  
3. Comparison between the maxima behavior in function of $\sigma$ for the increasing and decreasing $\sigma$ evolutions  	

### [Logistic map](https://en.wikipedia.org/wiki/Logistic_map)  

$$x(n+1) = \mu x(n)\cdot [1-x(n)]$$

1. Integration using [Euler method](https://en.wikipedia.org/wiki/Euler_method)  
2. Calculation of all possible [Lyapunov coefficients](https://en.wikipedia.org/wiki/Lyapunov_exponent) for different $\mu$ values  

### [Henon map](https://en.wikipedia.org/wiki/Henon_map)
	
$$\begin{align} 
x(n+1) &= y(n) + 1 - 1.4\cdot x^2(n)\\
y(n+1) &= 0.3 \cdot x(n) 
\end{align}$$

1. Integration using [Euler method](https://en.wikipedia.org/wiki/Euler_method)  
2. [Lyapunov coefficient](https://en.wikipedia.org/wiki/Lyapunov_exponent) calculation with [Benettin algorithm](https://cds.cern.ch/record/1453295/files/978-3-642-23666-2_BookBackMatter.pdf), for single and double trajectory  
3. [Fractal dimension](https://en.wikipedia.org/wiki/Fractal_dimension) calculation with [Grassberger-Procaccia algorithm](http://www.scholarpedia.org/article/Grassberger-Procaccia_algorithm) and comparison with [Kaplan-Yorke conjecture](https://en.wikipedia.org/wiki/Kaplan%E2%80%93Yorke_conjecture)  
	
### [Lorenz model](https://en.wikipedia.org/wiki/Lorenz_system)

$$\begin{align}
\dot{x} &= \sigma\cdot (y - x)\\
\dot{y} &= rx - y - xz\\
\dot{z} &= xy - bz
\end{align}$$
	
1. Integration with [Runge-Kutta](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods) 4th order  
2. [Lyapunov coefficient](https://en.wikipedia.org/wiki/Lyapunov_exponent) calculation with [Benettin algorithm](https://cds.cern.ch/record/1453295/files/978-3-642-23666-2_BookBackMatter.pdf), for single and triple trajectory  
3. [Fractal dimension](https://en.wikipedia.org/wiki/Fractal_dimension) calculation with [Grassberger-Procaccia algorithm](http://www.scholarpedia.org/article/Grassberger-Procaccia_algorithm) and comparison with [Kaplan-Yorke conjecture](https://en.wikipedia.org/wiki/Kaplan%E2%80%93Yorke_conjecture)  

### Multifractal measure

Generation and study of a multifractal measure with the multiplicative process 

## Online resources

1. [Runge-Kutta algorithm (Wikipedia)](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)
2. [Euler method (Wikipedia)](https://en.wikipedia.org/wiki/Euler_method)
3. [Lyapunov exponent (Wikipedia)](https://en.wikipedia.org/wiki/Lyapunov_exponent)
4. [Fractal dimension (Wikipedia)](https://en.wikipedia.org/wiki/Fractal_dimension)
5. [Benettin algorithm (CERN)](https://cds.cern.ch/record/1453295/files/978-3-642-23666-2_BookBackMatter.pdf)
6. [Kaplan-Yorke conjecture (Wikipedia)](https://en.wikipedia.org/wiki/Kaplan%E2%80%93Yorke_conjecture)
7. [Grassberger-Procaccia algorithm (Scholarpedia)](http://www.scholarpedia.org/article/Grassberger-Procaccia_algorithm)
8. [Logistic map (Wikipedia)](https://en.wikipedia.org/wiki/Logistic_map)
9. [Henon map (Wikipedia)](https://en.wikipedia.org/wiki/Henon_map)
10. [Lorenz model (Wikipedia)](https://en.wikipedia.org/wiki/Lorenz_system)

## carepackage.f90

The file `carepackage.f90` is a homemade compilation of useful subroutines not present by default in fortran90. Right now it contains

1. Ports of `arange` and `linspace` from python
2. A subroutine for linear fitting (without error calculation) with efficiency comparable with Gnuplot's fitting function, developed with [antonio-evangelista](https://github.com/antonio-evangelista)
3. A subroutine for searching the max value of a 1D vector, can be easily manipulated to find the minimum

## System requirements

All Work was done on [Ubuntu](https://ubuntu.com/) 20.04.5 LTS and [Arch Linux](https://archlinux.org/) using the following tools, available for all linux platforms (to my knowledge). Gnuplot and Make are optional, depending on user preference.  

1. [gfortran](https://gcc.gnu.org/wiki/GFortran) (gcc fortran compiler)  

   * for Ubuntu/Debian based systems  
  
         sudo apt install gfortran  
  
   * For Arch based systems  
    
         sudo pacman -S gfortran

2. [Gnuplot](http://www.gnuplot.info/) (powerful plotting software, also useful for data fitting)  

   * for Ubuntu/Debian based systems  
  
         sudo apt install gnuplot  
  
   * For Arch based systems  
    
         sudo pacman -S gnuplot

3. [Make](https://www.gnu.org/software/make/) (build system, usually already present on system)  

   * for Ubuntu/Debian based systems  
  
         sudo apt install make  
  
   * For Arch based systems  
    
         sudo pacman -S make

