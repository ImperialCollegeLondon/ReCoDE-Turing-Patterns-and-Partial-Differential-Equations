<!-- Your Project title, make it sound catchy! -->

# Turing Patterns and Partial Differential Equations

<!-- Provide a short description to your project -->

## Description

This code is a component of the Research Computing and Data Science Examples (ReCoDE) project. It comprises a non-linear partial differential equation (PDE) solver implemented in Fortran, designed to address both boundary value problems (BVP) and initial boundary value problems (IBVP) with temporal progression.
The solver's versatility allows it to handle problems in one or two dimensions and can accommodate single equations or pairs of coupled equations. This exemplar showcases several key features of the code, including:

- Integration with the Fortran Package Manager (FPM)
- Utilization of LAPACK libraries
- Modular architecture for enhanced maintainability and extensibility

To demonstrate its practical application, the code solves a PDE derived from a predator-prey model. This model is renowned for generating solutions that exhibit Turing patterns, which are observed in various biological systems, such as the skin patterns of pufferfish. While prior knowledge of PDEs and predator-prey models is not prerequisite for understanding this example, all relevant mathematical concepts will be elucidated in the subsequent sections.

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


## Academic

- Ordinary Differential Equations (ODEs): Ordinary differential equations are equations that involve functions of a single independent variable and their derivatives. These equations are fundamental in modeling various physical, biological, and economic phenomena where the rate of change of a quantity is related to the quantity itself. A discussion on ordinary differential equations can be found [here](https://math.libretexts.org/Courses/Monroe_Community_College/MTH_211_Calculus_II/Chapter_8%3A_Introduction_to_Differential_Equations/8.1%3A_Basics_of_Differential_Equations)
- Partial Differential Equations (PDEs): Partial differential equations are more complex than ODEs, as they involve functions of multiple independent variables and their partial derivatives. PDEs are crucial in describing many natural phenomena, including heat transfer, fluid dynamics, and quantum mechanics. They allow us to model systems that change with respect to multiple variables simultaneously, such as time and space. A discussion on partial differential equations can be found [here](https://en.wikipedia.org/wiki/Partial_differential_equation).
- Numerical Techniques: While some differential equations can be solved analytically, many real-world problems require numerical methods for approximation. The numerical techniques used for solving PDEs are often similar to those used for ODEs, but they must account for the additional complexity introduced by multiple variables.
- [Finite Difference Method](https://en.wikipedia.org/wiki/Finite_difference_method): One common numerical approach for solving PDEs is the finite difference method. This technique approximates derivatives by differences over small intervals. It discretises the continuous domain of the PDE into a grid or mesh, and the solution is computed at discrete points. This method is particularly useful for problems with regular geometries and is relatively straightforward to implement.

Understanding these concepts is crucial for working with the non-linear PDE solver included in this exemplar. The solver's ability to handle both BVPs and IVBPs makes it a versatile tool for a wide range of applications in scientific computing and mathematical modeling.

  
## System

- A Fortran compiler, such as `gfortran`: see [here](https://fortran-lang.org/learn/os_setup/install_gfortran/) for an installation guide.
- Fortran Package Manager (FPM): see [here](https://fpm.fortran-lang.org/install/index.html) for an installation guide.
- `BLAS` external library: see here for a `BLAS` installation guide [here](https://coral.ise.lehigh.edu/jild13/2016/07/27/install-lapack-and-blas-on-linux-based-systems/). Mac users can install with homebrew. 
- `LAPACK` external library: see the `LAPACK` documentation [here](https://www.netlib.org/lapack/) and [here](https://coral.ise.lehigh.edu/jild13/2016/07/27/install-lapack-and-blas-on-linux-based-systems/) for an installation guide. Mac users can install with homebrew. `BLAS` must be installed first.
- Optional: For visualisation of solutions we have made use of `MATLAB_R2023a` code. These are given in the `solver/examples` directory. Other visualisation software can be used.


## Getting Started
To begin working with this PDE solver, please follow these steps:

* **Introduction to Fortran Package Manager (FPM)**: Start by reading the review of Fortran Package Manager in Section 00 which is required for compiling this project. This will provide an introduction to FPM using the `fortran_fibonacci` repository as an example.

* **Theoretical Background**: Section 01 provides a general overview of the problem and the basic techniques employed in the solution.

* **Solver Configuration and Usage**: The three parts of Section 03 describes how to run and use the code.

* **Example Problems**: Consult Section 9 for detailed instructions on solving specific problems using this code.

* **Advanced Techniques**: The remaining sections explain specific coding techniques used in the development of this solver.

This structured learning sequence provides a comprehensive overview of the PDE solver, encompassing its theoretical foundations, practical implementation, operational use, and critically, the coding techniques employed in its development.

## License

This project is licensed under the [BSD-3-Clause license](LICENSE.md)
