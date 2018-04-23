"""
Here are all functions needed for the superbubble model
"""

##----------##
# Librairies #
##----------##
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy
from scipy.interpolate import interp1d

##---------##
# Functions #
##---------##

    # axes
rcParams['xtick.bottom'] = True
rcParams['xtick.top'] = True
rcParams['xtick.minor.visible'] = True

rcParams['ytick.left'] = True
rcParams['ytick.right'] = True
rcParams['ytick.minor.visible'] = True

def log_plot(figure_number, number_of_plot, x, y, label_name, title, xlabel, ylabel, symbol, linestyle, text):
    """
    Plot a log-log graphic
    Inputs:
        figure_number   :       define the number of the figure
        number_of_plot  :       define how many plot you want on one figure (with the same axis)
        x               :       x-vector
        y               :       y-array (line = one y-array and row = each different y-plot)
        label_name      :       legend of one y
        title           :       title of the plot
        xlabel          :       label of the x-axis
        ylabel          :       label of the y-axis
        symbol          :       symbol of one y
        linestyle       :       style of the line (drashed, etc)
        text            :       important parameters that you will write on the figure
    """
        # figure
    fig = plt.figure(figure_number, figsize = (8, 5))
    ax = fig.add_subplot(111)

        # Plot
    if number_of_plot > 1:

        for i in range (number_of_plot):
            y_plot = y[i]

            if len(symbol) > 1:
                rcParams['lines.marker'] = symbol[i]
                rcParams['lines.linestyle'] = linestyle[i]

            else:
                rcParams['lines.marker'] = symbol
                rcParams['lines.linestyle'] = linestyle

            if label_name == 'none':
                plt.loglog(x, y_plot)

            else:
                plt.loglog(x, y_plot, label = label_name[i])
                plt.legend(loc = 'best')

    elif label_name == 'none':

        rcParams['lines.marker'] = symbol
        rcParams['lines.linestyle'] = linestyle

        plt.loglog(x, y)

    else:
        rcParams['lines.marker'] = symbol
        rcParams['lines.linestyle'] = linestyle

        plt.loglog(x, y, label = label_name)
        plt.legend(loc = 'best')

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

        # grid
    plt.grid(color = 'k', alpha = 0.15, linestyle = ':')

        # text
    plt.text(0.5, 0.5, text, horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)

    return

def plot(figure_number, number_of_plot, x, y, label_name, title, xlabel, ylabel, symbol, linestyle, text):
    """
    Function to plot a linear graphic
    Inputs:
        figure_number   :       define the number of the figure
        number_of_plot  :       define how many plot you want on one figure (with the same axis)
        x               :       x-vector
        y               :       y-array (line = one y-array and row = each different y-plot)
        label_name      :       legend of one y
        title           :       title of the plot
        xlabel          :       label of the x-axis
        ylabel          :       label of the y-axis
        symbol          :       symbol of one y
        linestyle       :       style of the line (drashed, etc)
        text            :       important parameters that you will write on the figure
    """
        # figure
    fig = plt.figure(figure_number, figsize=(8,5))
    ax = fig.add_subplot(111)

        # Plot
    if number_of_plot > 1:

        for i in range (number_of_plot):
            y_plot = y[i]

            if len(symbol) > 1:
                rcParams['lines.marker'] = symbol[i]
                rcParams['lines.linestyle'] = linestyle[i]

            else:
                rcParams['lines.marker'] = symbol
                rcParams['lines.linestyle'] = linestyle

            if label_name == 'none':
                plt.plot(x, y_plot)

            else:
                plt.plot(x, y_plot, label = label_name[i])
                plt.legend(loc = 'best')

    elif label_name == 'none':
        rcParams['lines.marker'] = symbol
        rcParams['lines.linestyle'] = linestyle

        plt.plot(x, y)

    else:
        rcParams['lines.marker'] = symbol
        rcParams['lines.linestyle'] = linestyle

        plt.plot(x, y, label = label_name)
        plt.legend(loc = 'best')

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

        # grid
    plt.grid(color = 'k', alpha = 0.15, linestyle = ':')

        # text
    plt.text(0.5, 0.5, text, horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)

    return

def histogramme(figure_number, hist, label_name, title, xlabel, ylabel):
    """
    Return the histogramme of hist
    Inputs:
        figure_number   :       define the number of the figure
        hist            :       what you want to make the histogramme
        label_name      :       label of the data
        title           :       title of the histogramme
        xlabel          :       label of the x-axis
        ylabel          :       label of the y axis
    """
    plt.figure(figure_number, figsize=(8,5))
    plt.hist(hist, histtype = 'step', bins = 10, align = 'mid', label = label_name)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(loc = 'best')

    return

def random_PL(xmin, xmax, alpha, size = 1):

    """
    Generate random numbers with a specific size from a pdf propto x^alpha

    Inputs:
        xmin    :       lower boundary of the range
        xmax    :       upper boundary of the range
        alpha   :       exponent of the power-law (can be negative but then [xmin, xmax] can't contain 0)
        size    :       size of the random number (default = 1)
    """

    y = numpy.random.random(size = size)
    xmin_alpha, xmax_alpha = xmin**(alpha + 1), xmax**(alpha + 1)

    return (xmin_alpha + (xmax_alpha - xmin_alpha)*y)**(1./(alpha + 1))

def random_SN(xmin, xmax, size = 1):

    """
    Generate random number with a specific size from a uniform distribution on [xmin, xmax]

    Inputs:
        xmin    :       lower boundary of the range
        xmax    :       upper boundary of the range
        size    :       size of the random number (default = 1)
    """

    y = numpy.random.random(size = size)

    return (xmax-xmin) * y + xmin

def interpolation(x, y):
    """
    Return the interpolation of a function
    Inputs:
        x       :       x-axis of the function
        y       :       y-axis of the function
    """

    return interp1d(x, y, kind='linear', bounds_error = False, fill_value = 0.0)
