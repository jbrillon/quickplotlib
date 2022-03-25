# Quick Plotting Library

Library for quickly plotting with just 1 line of code by parameterizing Matplotlib commands. In addition, plots will always have consistent formatting.

## Examples:

<img src="https://raw.githubusercontent.com/jbrillon/quickplotlib/master/examples/figures/example_01.png" width="60%"></img>

### Dependencies:

This library requires the installation of the LaTeX font libraries locally. If you encounter compile errors, there are likely related to this; comment the following lines in `lib/quickplotlib.py` as a solution:
`matplotlibrc('text.latex', preamble='\usepackage{color}')
matplotlibrc('text', usetex=True)
matplotlibrc('font', family='serif')`