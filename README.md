Unfolding code for the neutron spectrometry with activation foils
========

This is the unfolding code for the neutron spectrometry with activation foils. It was developped for [the TRIUMF UCN experiment](http://www.triumf.ca/ucn).

1. External library and data
-----------

### ROOT
[ROOT](https://root.cern.ch/) is an object-oriented program and library developed by CERN. ROOT::Minuit2 is used for the unfolding. ROOT v5.34.32 was tested.

### Cross section data
Cross section data for the foils and covers are on xsec_foil and xsec_cover, respectively.
They were taken from [Evaluated Nuclear Data File (ENDF)](http://www.nndc.bnl.gov/exfor/endf.htm).
The first column is the energy (MeV) and the second column is the cross section (barn).
You can take new cross section data from ENDF.

2. Installation
------------------

### Installation of ROOT
- Download the codes from [ROOT website](https://root.cern.ch/downloading-root) and decompress it.
- Set environment variables as follows.
export ROOTSYS=(path to ROOT installation directory)
export PATH=${ROOTSYS}/bin:${PATH}
export LD_LIBRARY_PATH=${ROOTSYS}/lib:${LD_LIBRARY_PATH}

### Installation of this unfolding code
- make

3. Run the simulation
------------------

### Usage
- Unfold {input card file}

4. Input card file
------------------------

Input file of this code is the text file containing a sequence of option lines (so-called "cards").
"input.card" is a example of the input card file.

### GLOBAL
This card specifies the general information of the program.
- Spectrum:    Draw energy spectrum. 2:draw with error, 1:draw without error, 0: no.
- Covariance:  Draw covariance matrix of errors. 1:yes, 0:no.
- Correlation: Draw correlation matrix of errors. 1:yes, 0:no.

### COVERS
This card specifies the information about the covers surrouding the activation foils.
Properties of the covers are listed in the following style.
Name  Cross_section_file  Density[g/cm3]  Mass_number  Abundance  Thickness[cm]

### FOILS
This card specifies the information about the activation foils.
Properties of the foils are listed in the following style.
Name  Cross_section_file  Density[g/cm3]  Mass_number  Abundance  Thickness[cm]  Produced_RI  Error_on_produced_RI  Covers...
More than one cover can be used for one foil.
