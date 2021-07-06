kSEMAW is a workspace for the analysis of Spectrophotometric (SP), Ellipsometric (ELI) and
Photothermal Deflection Spectroscopy (PDS) measurements. The letter “k” indicates the use of the
Qt libraries.
Features
••••simulate SP, ELI and PDS measurements of a multilayer structure, being known the
thicknesses and the complex refractive indexes of each material composing the different
layers
calculate the complex refractive index and the thickness of a given layer (if ``thin'' ) from
experimental measurements (SP, ELI, PDS), being known the thicknesses and the complex
refractive indexes of all the other layers composing the structure
evaluate the mean value of physical quantities, weighted on a given international standard
spectrum (such as ASTM G173-03) or on own customized reference spectrum
predict the angular trend by using a realistic model or the equivalent model algorithm
Data and software availability
Code source files, user manual as well as a sample of working directories populated with assorted
files can be freely downloaded from SourceForge at the https://sourceforge.net/projects/ksemaw/
The FORTRAN executable requires the two libraries MINPACK and PGPLOT.
In order to harmonize the MINPACK sources with the GNU FORTRAN compiler (gfortran), a
slightly modified version is added to kSEMAW sources; that complies with the MINPACK
disclaimer (https://www.netlib.org/minpack/disclaimer).
The library PGPLOT is offered as a binary package in most of LINUX distros. As declared in the
copyright: ”PGPLOT is not public-domain software. However, it is freely available for non-
commercial use'”, like in the case of kSEMAW.
The C++ executable is devoted to control the Graphical User Interface (GUI), which is based on the
open source version of the Qt library, which is offered as binary packages in any LINUX distro.
License
The code source files are distributed as open source software under the GNU General Public
License as published by the Free Software Foundation version 3.
