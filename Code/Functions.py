##----------##
# Librairies #
##----------##
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import numpy
from scipy.interpolate import interp1d,interp2d

##---------##
# Functions #
##---------##

    # axes
rcParams['xtick.bottom'] = True
rcParams['xtick.top'] = True
rcParams['xtick.minor.visible'] = True
#rcParams['xtick.major.size'] = 10
#rcParams['xtick.minor.size'] = 7

rcParams['ytick.left'] = True
rcParams['ytick.right'] = True
rcParams['ytick.minor.visible'] = True

plt.rc('font', family='serif', size = 12)
figsize = (12, 8)

#rcParams['lines.linewidth'] = 8

def log_plot(figure_number, number_of_plot, x, y, xlabel, ylabel, symbol, linestyle, color, xmin, xmax, ymin, ymax, label_name = 'none', title = 'none', text = 'none'):
    """
    Plot a log-log graphic
    Inputs:
        figure_number   :   define the number of the figure
        number_of_plot  :   define how many plot you want on one figure (with the same axis)
        x               :   x-vector
        y               :   y-array (line = one y-array and row = each different y-plot)
        label_name      :   legend of one y (default = 'none')
        title           :   title of the plot (default = 'none')
        xlabel          :   label of the x-axis
        ylabel          :   label of the y-axis
        symbol          :   symbol of one y
        linestyle       :   style of the line (drashed, etc)
        text            :   important parameters that you will write on the figure (default = 'none')
        xmin            :   minimum of x
        xmax            :   maximum of x
        ymin            :   minimum of y
        ymax            :   maximum of y
    """
        # figure
    fig = plt.figure(figure_number, figsize = figsize)
    ax = fig.add_subplot(111)

        # limit
    plt.xlim((xmin, xmax))
    plt.ylim((ymin, ymax))

        # Plot
    if number_of_plot > 1:

        for i in range (number_of_plot):
            y_plot = y[i]

            if len(symbol) > 1:
                rcParams['lines.marker'] = symbol[i]
                rcParams['lines.linestyle'] = linestyle[i]
                col = color[i]

            else:
                rcParams['lines.marker'] = symbol
                rcParams['lines.linestyle'] = linestyle
                col = color

            if label_name == 'none':
                plt.loglog(x, y_plot, color = col)

            else:
                plt.loglog(x, y_plot, label = label_name[i], color = col)
                plt.legend(loc = 'best')

    elif label_name == 'none':

        rcParams['lines.marker'] = symbol
        rcParams['lines.linestyle'] = linestyle
        col = color

        plt.loglog(x, y, color = col)

    else:
        rcParams['lines.marker'] = symbol
        rcParams['lines.linestyle'] = linestyle
        col = color

        plt.loglog(x, y, label = label_name, col = color)
        plt.legend(loc = 'best')

    if title != 'none':

        plt.title(title)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

        # grid
    plt.grid(color = 'k', alpha = 0.15, linestyle = ':')

        # text
    if text != 'none':
        plt.text(0.5, 0.5, text, horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)

    return

def plot(figure_number, number_of_plot, x, y, xlabel, ylabel, symbol, linestyle, color, xmin, xmax, ymin, ymax, label_name = 'none', title = 'none', text = 'none'):
    """
    Function to plot a linear graphic
    Inputs:
        figure_number   :   define the number of the figure
        number_of_plot  :   define how many plot you want on one figure (with the same axis)
        x               :   x-vector
        y               :   y-array (line = one y-array and row = each different y-plot)
        label_name      :   legend of one y (default = 'none')
        title           :   title of the plot (default = 'none')
        xlabel          :   label of the x-axis
        ylabel          :   label of the y-axis
        symbol          :   symbol of one y
        linestyle       :   style of the line (drashed, etc)
        text            :   important parameters that you will write on the figure (default = 'none')
        xmin            :   minimum of x
        xmax            :   maximum of x
        ymin            :   minimum of y
        ymax            :   maximum of y
    """
        # figure
    fig = plt.figure(figure_number, figsize=figsize)
    ax = fig.add_subplot(111)

        # limit
    plt.xlim((xmin, xmax))
    plt.ylim((ymin, ymax))

        # Plot
    if number_of_plot > 1:

        for i in range (number_of_plot):
            y_plot = y[i]

            if len(symbol) > 1:
                rcParams['lines.marker'] = symbol[i]
                rcParams['lines.linestyle'] = linestyle[i]
                col = color[i]

            else:
                rcParams['lines.marker'] = symbol
                rcParams['lines.linestyle'] = linestyle
                col = color

            if label_name == 'none':
                plt.plot(x, y_plot, color = col)

            else:
                plt.plot(x, y_plot, label = label_name[i], color = col)
                plt.legend(loc = 'best')

    elif label_name == 'none':
        rcParams['lines.marker'] = symbol
        rcParams['lines.linestyle'] = linestyle
        col = color

        plt.plot(x, y, color = col)

    else:
        rcParams['lines.marker'] = symbol
        rcParams['lines.linestyle'] = linestyle
        col = color

        plt.plot(x, y, label = label_name, color = col)
        plt.legend(loc = 'best')

    if title != 'none':

        plt.title(title)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

        # grid
    plt.grid(color = 'k', alpha = 0.15, linestyle = ':')

        # text
    if text != 'none':
        plt.text(0.5, 0.5, text, horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)

    return

def semilog_plot(figure_number, number_of_plot, x, y, xlabel, ylabel, symbol, linestyle, color, xmin, xmax, ymin, ymax, label_name = 'none', title = 'none', text = 'none'):
    """
    Function to plot a linear graphic
    Inputs:
        figure_number   :   define the number of the figure
        number_of_plot  :   define how many plot you want on one figure (with the same axis)
        x               :   x-vector
        y               :   y-array (line = one y-array and row = each different y-plot)
        label_name      :   legend of one y (default = 'none')
        title           :   title of the plot (default = 'none')
        xlabel          :   label of the x-axis
        ylabel          :   label of the y-axis
        symbol          :   symbol of one y
        linestyle       :   style of the line (drashed, etc)
        text            :   important parameters that you will write on the figure (default = 'none')
        xmin            :   minimum of x
        xmax            :   maximum of x
        ymin            :   minimum of y
        ymax            :   maximum of y
    """
        # figure
    fig = plt.figure(figure_number, figsize=figsize)
    ax = fig.add_subplot(111)

        # limit
    plt.xlim((xmin, xmax))
    plt.ylim((ymin, ymax))

        # Plot
    if number_of_plot > 1:

        for i in range (number_of_plot):
            y_plot = y[i]

            if len(symbol) > 1:
                rcParams['lines.marker'] = symbol[i]
                rcParams['lines.linestyle'] = linestyle[i]
                col = color[i]

            else:
                rcParams['lines.marker'] = symbol
                rcParams['lines.linestyle'] = linestyle
                col = color

            if label_name == 'none':
                plt.plot(x, y_plot, color = col)

            else:
                plt.plot(x, y_plot, label = label_name[i], color = col)
                plt.legend(loc = 'best')

    elif label_name == 'none':
        rcParams['lines.marker'] = symbol
        rcParams['lines.linestyle'] = linestyle
        col = color

        plt.plot(x, y, color = col)

    else:
        rcParams['lines.marker'] = symbol
        rcParams['lines.linestyle'] = linestyle
        col = color

        plt.plot(x, y, label = label_name, color = col)
        plt.legend(loc = 'best')

    plt.yscale('log')

    if title != 'none':

        plt.title(title)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

        # grid
    plt.grid(color = 'k', alpha = 0.15, linestyle = ':')

        # text
    if text != 'none':
        plt.text(0.5, 0.5, text, horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)

    return

def log_plot_multi(figure_number, x, y, label_name, xlabel, ylabel, symbol, title = 'none'):
    """
    Plot a log-log graphic for two different y-axis
    Inputs:
        figure_number   :   define the number of the figure
        x               :   x-vector
        y               :   (2D) y-array
        label_name      :   legend of one y
        title           :   title of the plot (default = 'none')
        xlabel          :   label of the x-axis
        ylabel          :   labels of the y-axis
        symbol          :   symbol of one y
    """

        # figure
    fig = plt.figure(figure_number, figsize = figsize)

        # axes
    host = host_subplot(111, axes_class=AA.Axes)
    plt.subplots_adjust(right=0.75)

    par = host.twinx()

    par.axis["right"].toggle(all=True)

    host.set_xlabel(xlabel)
    host.set_ylabel(ylabel[0])
    p1, = host.loglog(x, y[0], symbol[0], label = label_name[0])

            # second axes
    par.set_ylabel(ylabel[1])
    p2, = par.loglog(x, y[1], symbol[1], label = label_name[1])

        # legend
    host.legend()
    host.axis["left"].label.set_color(p1.get_color())
    par.axis["right"].label.set_color(p2.get_color())

        # grid
    plt.grid(color = 'k', alpha = 0.15, linestyle = ':')

        # Title
    if title != 'none':

        plt.title(title)

        # Draw
    plt.draw()

    return

def plot_multi(figure_number, x, y, xlabel, ylabel, symbol, label_name, title = 'none'):
    """
    Plot a linear graphic for two different y-axis
    Inputs:
        figure_number   :   define the number of the figure
        x               :   x-vector
        y               :   (2D) y-array
        label_name      :   legend of one y
        title           :   title of the plot (default = 'none')
        xlabel          :   label of the x-axis
        ylabel          :   labels of the y-axis
        symbol          :   symbol of one y
    """

        # figure
    fig = plt.figure(figure_number, figsize = figsize)

        # axes
    host = host_subplot(111, axes_class=AA.Axes)
    plt.subplots_adjust(right=0.75)

    par = host.twinx()

    par.axis["right"].toggle(all=True)

            # host axes
    host.set_xlabel(xlabel)
    host.set_ylabel(ylabel[0])
    p1, = host.plot(x, y[0], symbol[0], label = label_name[0])

            # second axes
    par.set_ylabel(ylabel[1])
    p2, = par.plot(x, y[1], symbol[1], label = label_name[1])

        # legend
    host.legend()
    host.axis["left"].label.set_color(p1.get_color())
    par.axis["right"].label.set_color(p2.get_color())

        # grid
    plt.grid(color = 'k', alpha = 0.15, linestyle = ':')

        # Title
    if title != 'none':

        plt.title(title)

        # Draw
    plt.draw()

    return

def histogramme(figure_number, hist, xlabel, ylabel, len_bins, label_name = 'none', title = 'none'):
    """
    Return the histogramme of hist
    Inputs:
        figure_number   :   define the number of the figure
        hist            :   what you want to make the histogramme
        label_name      :   label of the data (default = 'none')
        title           :   title of the histogramme (default = 'none')
        xlabel          :   label of the x-axis
        ylabel          :   label of the y axis
    """
    plt.figure(figure_number, figsize=figsize)
    hist_max = numpy.max(hist)
    hist_min = numpy.min(hist)
    bins = int((hist_max - hist_min)/len_bins)

    if label_name == 'none':
        plt.hist(hist, histtype = 'step', bins = bins, align = 'mid')

    else:
        plt.hist(hist, histtype = 'step', bins = bins, align = 'mid', label = label_name)
        plt.legend(loc = 'best')

    if title != 'none':

        plt.title(title)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)


    return

def random_PL(xmin, xmax, alpha, size = 1):

    """
    Generate random numbers with a specific size from a pdf propto x^alpha

    Inputs:
        xmin    :   lower boundary of the range
        xmax    :   upper boundary of the range
        alpha   :   exponent of the power-law (can be negative but then [xmin, xmax] can't contain 0)
        size    :   size of the random number (default = 1)
    """

    y = numpy.random.random(size = size)
    xmin_alpha, xmax_alpha = xmin**(alpha + 1), xmax**(alpha + 1)

    return (xmin_alpha + (xmax_alpha - xmin_alpha)*y)**(1./(alpha + 1))

def interpolation1d(x, y):
    """
    Return the interpolation of a function
    Inputs:
        x       :       x-axis of the function
        y       :       y-axis of the function: y = f(x)
    """

    return interp1d(x, y, kind='linear', bounds_error = False, fill_value = 0.0)

def interpolation2d(x, y, z):
    """
    Return the interpolation of a function
    Inputs:
        x       :       x-axis of the function
        y       :       y-axis of the function
        z       :       z-axis of the function: z = f(x, y)
    """

    return interp2d(x, y, z, kind='linear', bounds_error = False, fill_value = 0.0)

def loglog_interpolation(x,y):
    """
    Return the interpolation of a function in log-log space
    Inputs:
        x       :       x-axis of the function
        y       :       y-axis of the function: y = f(x)
    """

    # Find non-zero values
    nzidx=y > 0.0

    return interp1d(numpy.log10(x[nzidx]),numpy.log10(y[nzidx]),kind='linear',bounds_error=False,fill_value=0.0)

def probability(Lum_HESS, Lum_Fermi, Lum_pwn, Lum_psr, Lum_HESS_CRb, Lum_Fermi_CRb, nit_tot, number_bin_t):

    """
    Return the probabilities

    Inputs:
        Lum_HESS        :   each iterations and time step of the gamma luminosity of CR in the HESS energy range
        Lum_Fermi       :   each iterations and time step of the gamma luminosity of CR in the Fermi energy range
        Lum_pwn         :   each iterations and time step of the gamma luminosity of PWNe in the HESS energy range
        Lum_psr         :   each iterations and time step of the gamma luminosity of PSRs in the Fermi energy range
        Lum_HESS_CRb    :   each time step of the gamma luminosity of CR background in the HESS energy range
        Lum_Fermi_CRb   :   each time step of the gamma luminosity of CR background in the Fermi energy range
        nit_tot         :   total number of iterations
        number_bin_t    :   len of the time array

    Outputs:
        Prob_HESS       :   probability to observe the SB (PWN + CR) in the HESS energy range
        Proba_HESS_CR   :   probability to observe only the CRs in the HESS energy range
        Prob_Fermi      :   probability to observe the SB (PSR + CR) in the Fermi energy range
        Proba_Fermi_CR  :   probability to observe only the CRs in the Fermi energy range
        Proba_pwn_psr   :   probability to observe no PWNe and no PSRs in the SB
    """

    nit_tot = float(nit_tot)

    Proba_HESS = numpy.zeros(number_bin_t)
    Proba_HESS_CR = numpy.zeros(number_bin_t)
    Proba_Fermi = numpy.zeros(number_bin_t)
    Proba_Fermi_CR = numpy.zeros(number_bin_t)
    Proba_pwn_psr = numpy.zeros(number_bin_t)

    for i in range (number_bin_t):

        Proba_HESS[i] = sum((Lum_HESS[:,i] + Lum_pwn[:,i]) > Lum_HESS_CRb[i])/nit_tot
        Proba_HESS_CR[i] = sum((Lum_HESS[:,i] > Lum_HESS_CRb[i]) * (Lum_pwn[:,i] <= Lum_HESS_CRb[i]))/nit_tot
        Proba_Fermi[i] = sum((Lum_Fermi[:,i] + Lum_psr[:,i]) > Lum_Fermi_CRb[i])/nit_tot
        Proba_Fermi_CR[i] = sum((Lum_Fermi[:,i] > Lum_Fermi_CRb[i]) * (Lum_psr[:,i] <= Lum_Fermi_CRb[i]))/nit_tot
        Proba_pwn_psr[i] = sum((Lum_pwn[:,i] <= 0.0) * (Lum_psr[:,i] <= 0.0))/nit_tot

    return Proba_HESS, Proba_HESS_CR, Proba_Fermi, Proba_Fermi_CR, Proba_pwn_psr
