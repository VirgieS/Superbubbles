"""
Here are all functions needed for the superbubble model
"""

##----------##
# Librairies #
##----------##
import numpy as np
from math import *
import matplotlib.pyplot as plt

##-----------------------------------------##
# Physical constants and Conversion factors #
##-----------------------------------------##
from Physical_constants import *
from Conversion_factors import *

##---------##
# Functions #
##---------##

def log_plot(figure_number, number_of_plot, x, y, label_name, title, xlabel, ylabel):
    """
    Function to plot a log-log graphic
    Inputs:
        figure_number:      define the number of the figure
        number_of_plot:     define how many plot do you want on one figure (with the same axis)
        x:                  x-vector
        y:                  y-array (line = one y-array and row = each different y-plot)
        label_name:         legend of one plots
        title:              title of the plot
        xlabel:             label of the x-axis
        ylabel:             label of the y-axis
    """
    plt.figure(figure_number)
    if number_of_plot > 1:
        for i in range (number_of_plot):
            y_plot = y[i,:]
            plt.plot(x, y_plot, '+', label=label_name[i])
        plt.legend(loc='best')
    else:
        plt.plot(x, y, '+')

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xscale('log')
    plt.yscale('log')
    #plt.show()
    return
