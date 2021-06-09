Harmonic Inversion
==================

This repository contains a program for the harmonic inversion analysis
as described in Fuchs *et al.* 2014 *J. Phys. A: Math. Theor.* **47**
125304.

The source code is in the directory `src` together with a Makefile for
generating the program.  You can compile the program by simply running
`make` or `make prog` which generates the executable `prog`.  Note that
the `LAPACK` and `FFTW3` libraries need to be installed on your computer
in order to compile successfully.  Run the program as `prog -h` to see
the help message with further explanations.
