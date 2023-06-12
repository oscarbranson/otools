# Oscar's Python Tools

Use at your own risk!

## Contents

### Geochem
Various functions for calculating (bio)mineral geochemistry from solution chemistry, given precipitation conditions. For example, implementations of DePaolo's Surface Kinetic Model (SKM), Rayleigh Fractionation, and Trans-Membrane-Transport models. 

### Chemistry
Functions for importing the periodic table of elements (scraped from webelements.com), and calculating the molecular mass of compounds.

### phreeqc
Functions for generating PHREEQC input strings, running them with phreeqpy and parsing the outputs.

### Peakshapes
Various peak shapes. Largely redundant... but hey!

### Uncertainties
A couple of helper wrappers to work with the `uncerainties` library.

### Plotting
Convenience functions for making plots. `rangecalc`, `spreadm`, `interval` and `unitpicker` are highlights.

## Installing

Git clone this repo, then run `pip install -e .` from the root directory. This will install the package in editable mode, so you can make changes to the code and they will be reflected in your environment.
