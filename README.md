# ksemaw_v2.6
kSEMAWc (Spectro-Ellipsometric Measurement Analysis Workbench) is a workspace for the analysis of Spectrophotometric (SP), Ellipsometric (ELI) and Photothermal Deflection Spectroscopy (PDS) measurements. The letter “k” indicates the use of the Qt libraries to generate the graphical user interface, typical of the desktop Linux KDE environment. The last letter “c” indicates the full software rewriting in C++;  considering this major improvement, the name was changed to kSEMAWc since version 1.0.0.

Features: 
1) simulate SP, ELI and PDS measurements of a multilayer structure, being known the
thicknesses and the complex refractive indexes of each material composing the different
layers;
2) calculate the complex refractive index and the thickness of a given layer (if thin ) from
experimental measurements (SP, ELI, PDS), being known the thicknesses and the complex
refractive indexes of all the other layers composing the structure;
3) optimize the parameters of the multilayer model by best-fitting of experimental measurements;
4) evaluate the mean value of physical quantities, weighted on a given international standard
spectrum (such as ASTM G173-03) or on own customized reference spectrum;
5) predict the angular trend by using a realistic model or the equivalent model algorithm.

Data and software availability:
kSEMAWc is written in C++ and is provided of a graphical user interface (GUI) for data input-output based on the Qt libraries, which are cross-platforms.

kSEMAWc generates 2D plots based on the Qwt library which is a graphics extension to the Qt GUI application framework.

Code source files, compiled Windows version, user manual as well as a sample of working directories populated with assorted files can be freely downloaded from

\texttt{https://github.com/mmonty1960/ksemaw}

kSEMAWc can be installed in different ways:  copying the 64bit Windows10 executable version folder or compiling the source files in Linux or Windows OS.

The installation of the precompiled Windows executable simply consists in:

• download the zipped file from https://github.com/mmonty1960/ksemaw

• extract the folder Workspace to any position on the user hard disk

• create a shortcut to the ksemawc.exe file contained in Workspace\qtSource\kSEMAWc\bin

License:
The code source files are distributed as open source software under the GNU General Public
License as published by the Free Software Foundation version 3.
