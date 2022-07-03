#-----------------------------------------------------
# import libraries
import numpy as np # NumPy: contains basic numerical routines
import sys
sys.path.append('../lib/')
from quickplotlib import plotfxn
#-----------------------------------------------------
# Example 1: Trigonometric functions
#-----------------------------------------------------
N = 400
x = np.linspace(0, 2*np.pi, N)
y1 = np.sin(x)
y2 = np.cos(x)
y3 = 1.0/np.sin(x) # csc(x)
y4 = 1.0/np.cos(x) # sec(x)
plotfxn(xdata=[x,x,x,x],
        ydata=[y1,y2,y3,y4],
        ylabel='$y$',
        xlabel='$x$',
        figure_filename='example_01',
        figure_filetype="png",
        title_label='Trigonometric Functions',
        legend_labels_tex=['$\\sin(x)$','$\\cos(x)$','$\\csc(x)$','$\\sec(x)$'],
        black_lines=False,
        xlimits=[0,2*np.pi],ylimits=[-2,2],
        which_lines_dashed=[2,3],
        remove_vertical_asymptotes_on_curve_number=[2,3])
#-----------------------------------------------------
# Example 2: Trigonometric function (single curve)
#-----------------------------------------------------
N = 400
x = np.linspace(0, 2*np.pi, N)
y1 = np.sin(x)
plotfxn(xdata=[x],
        ydata=[y1],
        ylabel='$y$',
        xlabel='$x$',
        figure_filename='example_02',
        figure_filetype="png",
        title_label='Trigonometric Function',
        legend_labels_tex=['$\\sin(x)$'],
        black_lines=True,
        legend_on=False,
        xlimits=[0,2*np.pi],ylimits=[-1,1])
#-----------------------------------------------------