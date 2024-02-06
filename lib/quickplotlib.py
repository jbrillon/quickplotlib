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
#-----------------------------------------------------
# define functions
#-----------------------------------------------------
def plot_any_axes(func,x,y,lc,mk,ls):
    # abstract plotting function to handle any kind of axes
    return func(x, y, color=lc, marker=mk, markersize=6, mfc='None', linestyle=ls)
#-----------------------------------------------------
def plot_lines(axes,xdata,ydata,ndata,lnstl,clr,mrkr,black_lines,
    which_lines_black,which_lines_only_markers,which_lines_markers,which_lines_dashed,
    markers,log_axes,lnstl_input,clr_input,mrkr_input,
    error_bars_on_curve_number,yerr_below,yerr_above,
    legend_labels_tex,leg_elements_input):
    
    leg_elements = []
    ls = lnstl[0]  # default line style (ls)
    lc = clr[0] # default line color (lc)
    mk = 'None'
    if(black_lines):
        lc = 'k' # set color to black
    
    color_index_shift = 0 

    for i in range(0,ndata):
        if(black_lines):
            ls = lnstl[i]
            mk = 'None' # reset to default
            if(i in which_lines_only_markers):
                ls = 'None'
                mk = mrkr[i]
            elif(i in which_lines_markers):
                mk = mrkr[i]
        else:
            # colour
            if(clr_input!=[]):
                lc = clr[i]
            elif(i in which_lines_black):
                lc = 'k'
                color_index_shift += 1
            else:
                lc = clr[i-color_index_shift]
                # warning: this color_index_shift could be buggy for when there's multiple desired black lines with colour ones
                #          -- no issues when only one black line is specified
            
            # linestyle
            if(lnstl_input!=[]):
                ls = lnstl[i]
            elif(i in which_lines_dashed):
                ls = "dashed"
            elif(i in which_lines_only_markers):
                ls = 'None'
            else:
                ls = 'solid' # default
            
            # markers
            if(mrkr_input!=[]):
                mk = mrkr[i]
            elif(i in which_lines_only_markers):
                mk = mrkr[i]
            elif(i in which_lines_markers):
                mk = mrkr[i]
            else:
                mk = 'None' # reset to default

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
                    plot_any_axes(axes.plot,x,y,lc,mk,ls)

        # add legend element 
        if(legend_labels_tex!=[] and leg_elements_input==[]):
            leg_elements.append(Line2D([0],[0], label=legend_labels_tex[i], color=lc, marker=mk, markersize=6, mfc='None', linestyle=ls))
        
        # plot command
        if(log_axes==None):
            plot_any_axes(axes.plot,x,y,lc,mk,ls)
        elif(log_axes=="both"):
            plot_any_axes(axes.loglog,x,y,lc,mk,ls)
        elif(log_axes=="x"):
            plot_any_axes(axes.semilogx,x,y,lc,mk,ls)
        elif(log_axes=="y"):
            plot_any_axes(axes.semilogy,x,y,lc,mk,ls)
    return leg_elements
#-----------------------------------------------------
def plotfxn(xdata=[],ydata=[],ylabel="ydata",xlabel="xdata",
            figure_filename="myfig",
            figure_filetype="pdf",
            title_label=" ",
            markers=False,
            legend_labels_tex=[],
            leg_elements_input=[],
            black_lines=False,
            xlimits=[],ylimits=[],
            log_axes=None,
            error_bars_on_curve_number=[],
            yerr_above=[],yerr_below=[],
            which_lines_black=[],
            which_lines_dashed=[],
            nlegendcols=1,legend_on=True,legend_inside=True,
            remove_vertical_asymptotes_on_curve_number=[],
            which_lines_only_markers=[],
            which_lines_markers=[],
            figure_size=(6,6),
            transparent_legend=False,
            legend_border_on=True,
            grid_lines_on=True,
            fig_directory = "figures",
            clr_input=[],mrkr_input=[],lnstl_input=[],
            axisTitle_FontSize=16,
            axisTickLabel_FontSize=14,
            legend_fontSize=16,
            legend_location="best",
            legend_anchor=[],
            second_leg_elements_input=[],
            second_leg_anchor=[],
            plot_zoomed_section=False,
            x_limits_zoom=[],y_limits_zoom=[],
            zoom_box_origin_and_extent=[]):
    print("---------------------------------------------")
    #-----------------------------------------------------
    # Safeguard for when empty data is passed
    #-----------------------------------------------------
    if(xdata==[] or ydata==[]):
        print("quickplotlib error: x or y data is empty")
        print("aborting...")
        return
    #-----------------------------------------------------
    # determine number of curves
    #-----------------------------------------------------
    if(hasattr(xdata[0],'__len__')):
        # then lists were passed, single or multiple curves
        ndata = int(len(xdata))
    else:
        # then numpy array was passed, single curve
        ndata=1; xdata=[xdata]; ydata=[ydata]
    #-----------------------------------------------------
    # colors, markers, and linestyles
    #-----------------------------------------------------
    # color
    if(clr_input!=[]):
        clr = clr_input
    else:
        clr = ['tab:blue','tab:red','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    # markers
    if(mrkr_input!=[]):
        mrkr = mrkr_input
    else:
        mrkr = ['o','s','^','d','v','>','<']
    # linestyles
    if(lnstl_input!=[]):
        lnstl = lnstl_input
    else:
        lnstl = ['solid','dashed','dashdot','dotted']
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
        fig, ax = plt.subplots(figsize=figure_size)
        # fig, ax = plt.subplots(figsize=(8,6))
    else:
        fig, ax = plt.subplots(figsize=(9,6))
    if(grid_lines_on):
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
    
    leg_elements_plt = plot_lines(ax,xdata,ydata,ndata,lnstl,clr,mrkr,black_lines,
                    which_lines_black,which_lines_only_markers,which_lines_markers,which_lines_dashed,
                    markers,log_axes,lnstl_input,clr_input,mrkr_input,
                    error_bars_on_curve_number,yerr_below,yerr_above,
                    legend_labels_tex,leg_elements_input)

    if(legend_on):
        if(leg_elements_input!=[]):
            leg_elements=leg_elements_input
        else:
            leg_elements=leg_elements_plt
        if(legend_inside):
            if(legend_anchor!=[]):
                leg = plt.legend(handles=leg_elements, loc=legend_location, bbox_to_anchor=(legend_anchor[0], legend_anchor[1]), ncol=nlegendcols, shadow=False, fancybox=True, fontsize=legend_fontSize, framealpha=1.0,edgecolor='inherit')
            else:
                leg = plt.legend(handles=leg_elements, loc=legend_location, ncol=nlegendcols, shadow=False, fancybox=True, fontsize=legend_fontSize, framealpha=1.0,edgecolor='inherit')
        else:
            leg = plt.legend(handles=leg_elements, loc="upper center", bbox_to_anchor=(1.2, 1.0), ncol=nlegendcols, shadow=False, fancybox=True, fontsize=legend_fontSize, framealpha=1.0,edgecolor='inherit')
        if(transparent_legend):
            leg.get_frame().set_facecolor('None')
        else:
            leg.get_frame().set_facecolor('w')
        if(not legend_border_on):
            leg.get_frame().set_edgecolor('w')
            leg.get_frame().set_linewidth(0.0)
        else:
            leg.get_frame().set_edgecolor('k')
            # leg.get_frame().set_linewidth(1.0)
        ax.add_artist(leg)
        if(second_leg_elements_input!=[]):
            second_leg_elements=second_leg_elements_input
            # legend_inside==True; legend_anchor==[]
            if(legend_location=="upper left"):
                second_legend_location="upper right"
            elif(legend_location=="upper left"):
                second_legend_location="upper right"
            else:
                second_legend_location="best"
            # adding the legend to the plot
            if(second_leg_anchor!=[]):
                second_legend_location=legend_location # for same side
                second_leg = plt.legend(handles=second_leg_elements, bbox_to_anchor=(second_leg_anchor[0], second_leg_anchor[1]), loc=second_legend_location, ncol=nlegendcols, shadow=False, fancybox=True, fontsize=legend_fontSize, framealpha=1.0,edgecolor='inherit')
            else:
                second_leg = plt.legend(handles=second_leg_elements, loc=second_legend_location, ncol=nlegendcols, shadow=False, fancybox=True, fontsize=legend_fontSize, framealpha=1.0,edgecolor='inherit')
            if(transparent_legend):
                second_leg.get_frame().set_facecolor('None')
            else:
                leg.get_frame().set_facecolor('w')
            if(not legend_border_on):
                second_leg.get_frame().set_edgecolor('w')
                second_leg.get_frame().set_linewidth(0.0)
            else:
                second_leg.get_frame().set_edgecolor('k')
            ax.add_artist(second_leg)

    if(plot_zoomed_section):
        x_limits_zoom
        x1, x2, y1, y2 = 7.5, 10.0, 0.010, 0.0135  # subregion of the original image
        axins = ax.inset_axes(zoom_box_origin_and_extent,xticklabels=[], yticklabels=[])
        axins.set_xlim(x_limits_zoom)
        axins.set_ylim(y_limits_zoom)
        axins.set_xticks([])
        axins.set_yticks([])
        
        dummy = plot_lines(axins,xdata,ydata,ndata,lnstl,clr,mrkr,black_lines,
                        which_lines_black,which_lines_only_markers,which_lines_markers,which_lines_dashed,
                        markers,log_axes,lnstl_input,clr_input,mrkr_input,
                        error_bars_on_curve_number,yerr_below,yerr_above,
                        legend_labels_tex,leg_elements_input)

        ax.indicate_inset_zoom(axins, edgecolor="black")

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
#-----------------------------------------------------
def plotfield(xdata=[],ydata=[],udata=[],vdata=[],
            ylabel="y",xlabel="x",
            figure_filename="myfig",
            figure_filetype="pdf",
            title_label=" ",
            xlimits=[],ylimits=[],
            figure_size=(6,6),
            grid_lines_on=False,
            fig_directory = ".",
            black_lines=True,
            clr_input=[],mrkr_input=[],lnstl_input=[],
            axisTitle_FontSize=16,
            axisTickLabel_FontSize=14,
            legend_fontSize=16):
    print("---------------------------------------------")
    #-----------------------------------------------------
    # Safeguard for when empty data is passed
    #-----------------------------------------------------
    if(xdata==[] or ydata==[] or udata==[] or vdata==[]):
        print("quickplotlib error: x, y, u, or v data is empty")
        print("aborting...")
        return
    #-----------------------------------------------------
    # determine number of curves
    #-----------------------------------------------------
    if(hasattr(xdata[0],'__len__')):
        # then multiple curves
        ndata = int(len(xdata))
        print("quickplotlib error: %i sets of data were passed to quiver(), cannot pass more than 1")
        print("aborting...")
        return
    #-----------------------------------------------------
    # colors, markers, and linestyles
    #-----------------------------------------------------
    # color
    if(clr_input!=[]):
        clr = clr_input
    else:
        clr = ['tab:blue','tab:red','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    # markers
    if(mrkr_input!=[]):
        mrkr = mrkr_input
    else:
        mrkr = ['o','s','^','d','v','>','<']
    # linestyles
    if(lnstl_input!=[]):
        lnstl = lnstl_input
    else:
        lnstl = ['solid','dashed','dashdot','dotted']
    #-----------------------------------------------------
    # plotting:
    #-----------------------------------------------------
    print('Plotting: ' + fig_directory + "/" + figure_filename + "." + figure_filetype)
    fig, ax = plt.subplots(figsize=figure_size)
    if(grid_lines_on):
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
    lc = clr[0] # default line color (lc)
    if(black_lines):
        lc = 'k' # set color to black
    # quiver plot
    ax.quiver(xdata, ydata, udata, vdata, color=lc)

    plt.tight_layout()
    print('\t ... Saving figure ...')
    plt.savefig(fig_directory+"/"+figure_filename+'.'+figure_filetype,format=figure_filetype,dpi=500)
    plt.close()
    print('\t     Saved.')
    print("---------------------------------------------")
#-----------------------------------------------------
#=====================================================