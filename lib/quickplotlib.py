#-----------------------------------------------------
# Python code for easy post-processing
# - Written by Julien Brillon
# - Email: julien.brillon@mail.mcgill.ca
# - McGill University
#-----------------------------------------------------
# Import libraries
import numpy as np # NumPy: contains basic numerical routines
import scipy # SciPy: contains additional numerical routines to numpy
import matplotlib.pyplot as plt # Matlab-like plotting
import matplotlib.colors as colors # needed for sparsity pattern
import matplotlib # needed for the sparsity pattern code
#-----------------------------------------------------
# plotting essentials
#-----------------------------------------------------
from matplotlib.lines import Line2D
from matplotlib import rc as matplotlibrc
matplotlibrc('text.latex', preamble='\\usepackage{color}')
matplotlibrc('text', usetex=True)
matplotlibrc('font', family='serif')
clr = ['tab:blue','tab:red','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
mrkr = ['o','s','^','d','v','>','<']
lnstl = ['solid','dashed','dashdot','dotted']
# Font sizes
axisTitle_FontSize = 16
axisTickLabel_FontSize = 14
legend_fontSize = 16
fig_directory = "figures"
#-----------------------------------------------------
# define functions
#-----------------------------------------------------
def plot_any_axes(func,x,y,lc,mk,ls):
    # abstract plotting function to handle any kind of axes
    return func(x, y, color=lc, marker=mk, markersize=6, mfc='None', linestyle=ls)
#-----------------------------------------------------
def plotfxn(xdata=[],ydata=[],ylabel="ydata",xlabel="xdata",
            figure_filename="myfig",
            figure_filetype="pdf",
            title_label=" ",
            markers=False,
            legend_labels_tex=[],
            black_lines=False,
            xlimits=[],ylimits=[],
            log_axes=None,
            error_bars_on_curve_number=[],
            yerr_above=[],yerr_below=[],
            which_lines_black=[],
            which_lines_dashed=[],
            nlegendcols=1,legend_on=True,legend_inside=True,
            remove_vertical_asymptotes_on_curve_number=[],
            which_lines_only_markers=[]):
    print("---------------------------------------------")
    #-----------------------------------------------------
    # determine number of curves
    #-----------------------------------------------------
    if(hasattr(xdata[0],'__len__')):
        # then multiple curves
        ndata = int(len(xdata))
    #-----------------------------------------------------
    # pre-plotting data manipulation:
    #-----------------------------------------------------
    if(remove_vertical_asymptotes_on_curve_number!=[]):
        tolerance_for_vertical_asymptote_value_difference = 1.0
        for i in remove_vertical_asymptotes_on_curve_number:
            N = int(len(ydata[i]))
            for j in range(0,N-1):
                condition1 = ydata[i][j+1]*ydata[i][j]<0
                condition2 = np.abs(ydata[i][j+1]-ydata[i][j]) > tolerance_for_vertical_asymptote_value_difference
                if(condition1 and condition2):
                    ydata[i][j+1]=np.nan; ydata[i][j]=np.nan;
    #-----------------------------------------------------
    # plotting:
    #-----------------------------------------------------
    print('Plotting: ' + fig_directory + "/" + figure_filename + "." + figure_filetype)
    if(legend_inside):
        fig, ax = plt.subplots(figsize=(6,6))
    else:
        fig, ax = plt.subplots(figsize=(9,6))
    plt.grid()
    ax.set_xlabel(xlabel,fontsize=axisTitle_FontSize)
    ax.set_ylabel(ylabel,rotation=90,fontsize=axisTitle_FontSize)
    if(title_label!=" "):
        plt.title(title_label,fontsize=axisTitle_FontSize)
    plt.setp(ax.get_xticklabels(),fontsize=axisTickLabel_FontSize); plt.setp(ax.get_yticklabels(),fontsize=axisTickLabel_FontSize);
    if(xlimits!=[]):
        ax.set_xlim(xlimits)
    if(ylimits!=[]):
        ax.set_ylim(ylimits)
    ls = lnstl[0]  # default line style (ls)
    lc = clr[0] # default line color (lc)
    mk = 'None'
    if(black_lines):
        lc = 'k' # set color to black
    
    color_index_shift = 0
    leg_elements = []

    for i in range(0,ndata):
        if(black_lines):
            ls = lnstl[i]
            mk = 'None' # reset to default
            if(i in which_lines_only_markers):
                ls = 'None'
                mk = mrkr[i]
        else:
            if(i in which_lines_black):
                lc = 'k'
                color_index_shift += 1
            else:
                lc = clr[i-color_index_shift]
                # warning: this color_index_shift could be buggy for when there's multiple desired black lines with colour ones
                #          -- no issues when only one black line is specified
            mk = 'None' # reset to default
            if(i in which_lines_dashed):
                ls = "dashed"
            elif(i in which_lines_only_markers):
                ls = 'None'
                mk = mrkr[i]
            else:
                ls = lnstl[0]
        x = xdata[i]; y = ydata[i];
        if(markers):
            mk = mrkr[i]
            if(log_axes==None):
                if(error_bars_on_curve_number!=[] and i==error_bars_on_curve_number):
                    yerr = np.array([yerr_below,yerr_above])
                    fmt_string = lc+mrkr[i]
                    plt.errorbar(x, y, yerr, fmt=fmt_string, mfc='None')
                    # plt.errorbar(x, y, yerr, color=lc, marker=mrkr[i], markersize=6, mfc='None', linestyle='None')
                else:
                    plot_any_axes(plt.plot,x,y,lc,mk,ls)
        
        # add legend element 
        leg_elements.append(Line2D([0],[0], label=legend_labels_tex[i], color=lc, marker=mk, markersize=6, mfc='None', linestyle=ls))
        
        # plot command
        if(log_axes==None):
            plot_any_axes(plt.plot,x,y,lc,mk,ls)
        elif(log_axes=="both"):
            plot_any_axes(plt.loglog,x,y,lc,mk,ls)
        elif(log_axes=="x"):
            plot_any_axes(plt.semilogx,x,y,lc,mk,ls)
        elif(log_axes=="y"):
            plot_any_axes(plt.semilogy,x,y,lc,mk,ls)
    if(legend_on):
        if(legend_inside):
            leg = plt.legend(handles=leg_elements, loc="best", ncol=nlegendcols, shadow=False, fancybox=True, fontsize=legend_fontSize, framealpha=1.0,edgecolor='inherit')
        else:
            leg = plt.legend(handles=leg_elements, loc="upper center", bbox_to_anchor=(1.2, 1.0), ncol=nlegendcols, shadow=False, fancybox=True, fontsize=legend_fontSize, framealpha=1.0,edgecolor='inherit')
    plt.tight_layout()
    print('\t ... Saving figure ...')
    plt.savefig(fig_directory+"/"+figure_filename+'.'+figure_filetype,format=figure_filetype,dpi=500)
    plt.close()
    print('\t     Saved.')
    print("---------------------------------------------")
#-----------------------------------------------------
#=====================================================
# set the colormap and centre the colorbar
class MidpointNormalize(colors.Normalize):
    """
    Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

    e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
    """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))
#-----------------------------------------------------
def plot_matrix_sparsity_pattern(A=np.random.random((6,6)),
                                 colour_toggle='y',
                                 cutOff=0.01,
                                 figure_filename = "spy",
                                 figure_filetype = 'pdf'
                                 ):
    print('-----------------------------------------------------')
    # get matrix A number of rows and columns
    nRow, nCol = A.shape

    # Plot
    figure_title = "Sparsity Pattern, Matrix Size: %ix%i" % (nRow,nCol)
    figure_title_print = "Sparsity Pattern, Matrix Size: %ix%i" % (nRow,nCol)
    print('Plotting: ' + figure_title_print)
    fig = plt.figure(figure_title)
    plt.title(figure_title)

    if colour_toggle == 'y': # Plot in colour
        elev_min=np.min(A)
        elev_max=np.max(A)
        mid_val=0
        cmap=matplotlib.cm.RdBu_r
        plt.imshow(A, cmap=cmap, clim=(elev_min, elev_max), norm=MidpointNormalize(midpoint=mid_val,vmin=elev_min, vmax=elev_max))
        plt.colorbar(label='Operator Weights')
    elif colour_toggle == 'n': # Plot black and white without color bar
        # plt.spy(A,precision=0.001,markersize=5,color='k')
        plt.spy(A,precision=cutOff,markersize=int(25.0*9.0/min(nRow,nCol)),color='k')

    plt.xticks(np.arange(-0.5,nCol),0*list(range(0, nCol+1)))
    plt.yticks(np.arange(-0.5,nRow),0*list(range(0, nRow+1)))
    plt.grid(linewidth=0.001)
    plt.tight_layout()
    print('\t ... Saving figure ...')
    plt.savefig(fig_directory+"/"+figure_filename+'.'+figure_filetype,format=figure_filetype,dpi=500)
    plt.close()
    print('\t     Saved.')
    print('-----------------------------------------------------')
#=====================================================
