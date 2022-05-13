# Quick Plotting Library

Library for quickly plotting with just 1 line of code by parameterizing Matplotlib commands. In addition, plots will always have consistent formatting.

## Examples:

<img src="https://raw.githubusercontent.com/jbrillon/quickplotlib/master/examples/figures/example_01.png" width="60%"></img>

### Dependencies:

This library requires the installation of the LaTeX font libraries locally. 

On a clean install of Ubuntu 20.04, you may run the following script:
https://github.com/jbrillon/quickplotlib/blob/doc/install_ubuntu2004.sh

If you encounter compile errors, there are likely related to the LaTeX dependencies; as a solution, comment the following lines in `lib/quickplotlib.py`:\
`matplotlibrc('text.latex', preamble='\usepackage{color}')`\
`matplotlibrc('text', usetex=True)`\
`matplotlibrc('font', family='serif')`
