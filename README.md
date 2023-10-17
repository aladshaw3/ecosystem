[![Checks](https://github.com/aladshaw3/ecosystem/actions/workflows/build-ecosystem.yml/badge.svg)](https://github.com/aladshaw3/ecosystem/actions/workflows/build-ecosystem.yml)

# README #

This README file outlines how to get setup using the ecosystem codes and executables. Please see below for guidelines on how to cite this repository if you are using any of these modules for a publication. Thanks. 

### What is this repository for? ###

* This repository is for the ecosystem tool set designed and developed by Austin Ladshaw at the Georgia Institute of Technology and Oak Ridge National Laboratory. Tools and algorithms developed here are designed for simulation and data anaylsis associated with adsorption phenomena in air and water chemistry. However, several sets of tools here are also generalized and could be used in nearly any realm of science and engineering. Below is a brief description of some of the modules...

	+ MACAW - templated C++ matrix object for simple matrix math
	+ LARK - C-style library of linear and non-linear solvers
	+ DOVE - C++ object for creating and solving systems of ODEs
	+ FINCH - C-style library for simulating 1D conservation laws with a MUSCL scheme
	+ CROW - C++ object for creating and simulating reaction mechanisms 
	+ IBIS - C++ object for solving radioactive decay chains
	+ EEL - C++ object to store and register specific atomic information
	+ MOLA - C++ object to store and register specific molecular information
	+ SHARK - C++ object for creating and simulating aqueous chemistry and adsorption
	+ YAML_WRAPPER - C++ object for reading and writing yaml formatted files 
	+ FAIRY - C++ object to utilize fission yield databases to predict fission products
	+ MAGPIE - C-style object to simulate multi-species gas phase adsorption thermodynamics
	+ GSTA_OPT - C-style optimization routine to process gas phase adsorption data 

* Version: 1.0.0

* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* Requirements: 
	+ gcc/g++ version 4.7 or newer (or llvm version 4.2 or newer)
	+ python 3.7 or newer (recommended for running python scripts)

* Summary of set up

	+ clone the ecosystem project to your machine
	+ type "make" in the ecosystem directory to build the project
	+ type "make install" to install the executable to your /usr/local/bin
	+ type "make lib" to build a shared library (used with some python scripts)
	+ type "make all" to build the executable and a shared library

* Configuration

	+ For Mac and Linux, no special or additional configuration required
	+ For Windows, I recommend using Cygwin to install the bash shell and bash terminal to your machine, then building the executable from that terminal. You will need gcc, g++, gdb, and make (minimum requirements). This has been shown to work for Windows operating systems. 

* Dependencies: No outside dependencies

* Database configuration: 
	+ Necessary database files include Fission Yields and a Nuclide Library
	+ All database files are distributed with the software source code
	+ Optional database files include atmospheric information (see input_files/CARDINAL)

* How to run tests

	+ The executable will have built-in tests you can run to check for runtime errors. To run these tests, first build the project from source, then type "eco -t" to open a test menu. You can then choose a test to run.
	+ Tests will report errors if any are present (no error messages = good to go)

* Deployment instructions: No special instructions

### Contribution guidelines ###

* Writing tests: Test your code before committing anything to the repository

* Code review: Stay organized and format code to match existing code format

* Other guidelines: 

	+ If you are changing any of the provided source code, it may be better to create your own branch off the master branch so that you do not alter any existing algorithms. 
	+ Please contact Austin Ladshaw (aladshaw3@outlook.com) to request any source code changes and do not try to change the source code yourself. 

### Who do I talk to? ###

* Lead Developer: Austin Ladshaw (austin.ladshaw@span.io, aladshaw3@outlook.com)
* Trajectory Contributor: Alex Wiechert (awiechert3@gatech.edu)

### How to cite this repository ###

* Please refer to the citation guidelines for your particular journal or publisher for citing software. If no guidelines are available, you may cite this work as the following:
* Ladshaw, A., Wiechert, A., Tsouris, C., Yiacoumi, S., "Ecosystem Software - A C/C++ library for environmental chemistry and adsorption," Version <(version #)>, Available at https://bitbucket.org/gitecosystem/ecosystem, Accessed on (Month) (Day), (Year). 
