<!-- Your Project title, make it sound catchy! -->

# Turing Patterns and Partial Differential Equations

<!-- Provide a short description to your project -->

## Description

This code is part of the Research Computing and Data Science Examples (ReCoDE) project. 
It is a non-linear partial differential equation (PDE) solver written in Fortran, designed to tackle both boundary value problems (BVP) and initial boundary value problems (IBVP) with temporal marching.
The solver is versatile, capable of handling problems in one or two dimensions and can solve single equations or pairs of coupled equations. This exemplar highlights the code's useful features, including the use of the Fortran Package Manager (FPM), use of LAPACK and the codes modular nature.
As a practical example, the code addresses a PDE from a predator-prey model, known for producing solutions that feature Turing patterns - the patterns that feature on many animals such as pufferfish. While knowledge of PDEs and predator-prey models is not required, all necessary mathematical concepts will be explained below.

<!-- What should the students going through your exemplar learn -->

## Learning Outcomes

- Complied Codes
- Fortran Package Manager (FPM)
- Modular Codes
- Multipurpose Codes
- Solving Mathematical Problems (PDEs)
- Generalising Problems
- Discretisation in multiple dimensions
- Use of external libraries (LAPACK and BLAS)
- Testing Fortran code


<!-- How long should they spend reading and practising using your Code.
Provide your best estimate -->

| Task       | Time    |
| ---------- | ------- |
| Reading    | 7 hours |
| Practising | 7 hours |

## Requirements

This exemplar is for entry-level researchers with basic knowledge of Fortran syntax. An [Introduction to Fortran Course](https://www.imperial.ac.uk/students/academic-support/graduate-school/professional-development/doctoral-students/research-computing-data-science/courses/intro-to-fortran/) 
is available to help you get started.


### Academic

A basic overview of differential equations is required:
- A discussion on ordinary differential equaitons can be found [here[(https://math.libretexts.org/Courses/Monroe_Community_College/MTH_211_Calculus_II/Chapter_8%3A_Introduction_to_Differential_Equations/8.1%3A_Basics_of_Differential_Equations#:~:text=A%20differential%20equation%20is%20an,are%20substituted%20into%20the%20equation.)
- We will be looking at partial differential equations, which are equations dependent on more than one variable. The numerical techniques to solve these are very similar to the ones for ordinary differential equations. A discussion on partial differential equations can be found [here](https://en.wikipedia.org/wiki/Partial_differential_equation).
- We will solve our equations using [finite differences](https://en.wikipedia.org/wiki/Finite_difference_method).
  
### System

- A Fortran compiler, such as `gfortran`: see [here](https://fortran-lang.org/learn/os_setup/install_gfortran/) for an installation guide.
- Fortran Package Manager (FPM): see [here](https://fpm.fortran-lang.org/install/index.html) for an installation guide.
- `BLAS` external library: see here for a `BLAS` installationa guide [here](https://coral.ise.lehigh.edu/jild13/2016/07/27/install-lapack-and-blas-on-linux-based-systems/). Mac users can install with homebrew. 
- `LAPACK` external library: see the `LAPACK` documentaion [here](https://www.netlib.org/lapack/) and [here](https://coral.ise.lehigh.edu/jild13/2016/07/27/install-lapack-and-blas-on-linux-based-systems/) for an intallation guide. Mac users can install with homebrew. `BLAS` must be installed first.
- Optional: For visualisation of solutions we have made use of `MATLAB_R2023a`. Other visulation software can be used.


## Getting Started

- Start by going through the notebook [`00.fortran_package_manager.md`](notebooks/00.fortran_package_manager.md). This will give you an introduction into FPM by going through and running the code in the repository `fortran_fibonacci`. 
- Then read the [`01.theory.md`](https://github.com/ImperialCollegeLondon/ReCoDE-Solving-Singular-PDEs-in-Fortran/blob/main/notebooks/01.theory.md). This will give a general overview of the problem we are trying to solve, and discuss the basic techniques will will use.
- Reading through notebooks [`03-1.solver_settings.md`](notebooks/03-1.solver_settings.md), [`03-2.equation_definitions.md`](notebooks/03-2.equation_definitions.md) and [`03-3.fpm_settings.md`](main/notebooks/03-3.fpm_settings.md) will enable you to run and use the code.
- Following these, notebook [`09.example_equations.md`](notebooks/09.example_equations.md) will detail how to solve particular problems with the code.
- The remaining notebooks explain specific techniques used in the coding process.

## License

This project is licensed under the [BSD-3-Clause license](LICENSE.md)
