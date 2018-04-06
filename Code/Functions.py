"""
Here are all functions needed for the superbubble model
"""

##----------##
# Librairies #
##----------##
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoLocator, LogLocator
import numpy
import scipy.integrate as integrate
from scipy.special import erf, erfc
from scipy.interpolate import interp1d

##-----------------------------------------##
# Physical constants and Conversion factors #
##-----------------------------------------##
from Physical_constants import *
from Conversion_factors import *

##---------##
# Functions #
##---------##

def log_plot(figure_number, number_of_plot, x, y, label_name, title, xlabel, ylabel, symbol):
    """
    Function to plot a log-log graphic
    Inputs:
        figure_number   :       define the number of the figure
        number_of_plot  :       define how many plot do you want on one figure (with the same axis)
        x               :       x-vector
        y               :       y-array (line = one y-array and row = each different y-plot)
        label_name      :       legend of one plots
        title           :       title of the plot
        xlabel          :       label of the x-axis
        ylabel          :       label of the y-axis
        symbol          :       symbol
    """
        # figure
    fig, ax = plt.subplots(figsize=(10,7))
    fig = plt.figure(figure_number)

        # axes
    #xmin = min(x)
    #xmax = max(x)
    #plt.xlim(xmin, xmax)
    #plt.xticks(numpy.logspace(numpy.log10(xmin), numpy.log10(xmax), 10, endpoint = True))
    #plt.xtick.minor.visible = True
    ax.xaxis.set_minor_locator(LogLocator())

    #ymin = min(y)
    #ymax = max(y)
    #plt.ylim(ymin, ymax)
    #plt.yticks(numpy.logspace(numpy.log10(ymin), numpy.log10(ymax), 10, endpoint = True))
    ax.yaxis.set_minor_locator(LogLocator())


    if number_of_plot > 1:

        for i in range (number_of_plot):
            y_plot = y[i]

            if len(symbol) > 1:
                sym = symbol[i]

            else:
                sym = symbol

            if label_name == 'none':
                plt.loglog(x, y_plot, sym)

            else:
                plt.loglog(x, y_plot, sym, label = label_name[i])
                plt.legend(loc = 'best')

    elif label_name == 'none':
        plt.loglog(x, y, symbol)

    else:
        plt.loglog(x, y, symbol, label = label_name)
        plt.legend(loc = 'best')

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    return

def plot(figure_number, number_of_plot, x, y, label_name, title, xlabel, ylabel, symbol):
    """
    Function to plot a linear graphic
    Inputs:
        figure_number   :       define the number of the figure
        number_of_plot  :       define how many plot do you want on one figure (with the same axis)
        x               :       x-vector
        y               :       y-array (line = one y-array and row = each different y-plot)
        label_name      :       legend of one plots
        title           :       title of the plot
        xlabel          :       label of the x-axis
        ylabel          :       label of the y-axis
        symbol          :       symbol
    """
        # figure
    fig, ax = plt.subplots(figsize=(10,7))
    fig = plt.figure(figure_number)

        # axes
    #xmin = min(x)
    #xmax = max(x)
    #plt.xlim(xmin, xmax)
    #plt.xticks(numpy.linspace(xmin, xmax, 10, endpoint = True))
    ax.xaxis.set_minor_locator(AutoLocator())

    #ymin = min(y)
    #ymax = max(y)
    #plt.ylim(ymin, ymax)
    #plt.yticks(numpy.linspace(ymin, ymax, 10, endpoint = True))
    ax.yaxis.set_minor_locator(AutoLocator())

    if number_of_plot > 1:

        for i in range (number_of_plot):
            y_plot = y[i]

            if len(symbol) > 1:
                sym = symbol[i]

            if label_name == 'none':
                plt.plot(x, y_plot, symbol)

            else:
                plt.plot(x, y_plot, sym, label = label_name[i])
                plt.legend(loc = 'best')

    elif label_name == 'none':
        plt.plot(x, y, symbol)

    else:
        plt.plot(x, y, symbol, label = label_name)
        plt.legend(loc = 'best')

    #plt.tight_layout()
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    #fig.savefig(pathfigure+name_figure)
    return

def integration_log(x, y):

    """
    Return the integrale of a function in log-log scale

    Parameters :
        x           : abscisse of the function
        y           : function that we integrate along x : y = f(x)
    """

    #Looking for a and b for y = a*x^b
    def calculate_ab(xi, xf, yi, yf):
        logxi = numpy.log(xi)
        logxf = numpy.log(xf)
        logyi = numpy.log(yi)
        logyf = numpy.log(yf)
        b = (logyf - logyi)/(logxf - logxi)
        loga = logyi - b*logxi
        a = numpy.exp(loga)
        a = numpy.nan_to_num(a)
        return a, b

    #Calculate deltaS from deltaS = int from xi to xf a*x^b
    def delta_S(xi, xf, yi, yf):
        [a, b] = calculate_ab(xi, xf, yi, yf)
        return a/(b+1)*(xf**(b+1) - xi**(b+1))

    integral = 0

    # Because the function integral_log works only if there is more than two elements not zero
    idx=(y > 0.0)
    #idy=(y < 0.0) # only used if the function has negative value
    idt = idx #+ idy
    if sum(idt) > 2:

        x = x[idt]
        y = y[idt]

        #Calculate total integral from init to final a*x^b
        deltaS = 0

        for i in range (1, len(x)):
            deltaS = delta_S(x[i-1], x[i], y[i-1], y[i])
            integral = integral + deltaS

            integral = numpy.nan_to_num(integral)

    return integral

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

    return xmax * y + xmin

def interpolation(x, y):
    """
    Return the interpolation of a function
    Inputs:
        x       :       x-axis of the function
        y       :       y-axis of the function
    """

    return interp1d(x, y, kind='linear', bounds_error = False, fill_value = 0.0)
