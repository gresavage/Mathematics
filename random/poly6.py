import numpy as np
import matplotlib.pyplot as plt
from math import *
from IPython.html import widgets
from IPython.html.widgets import interact
from IPython.display import display
get_ipython().magic(u'matplotlib inline')

def factored6(x):
    return (x - 1.)**6.


def expanded6(x):
    return x**6. - 6.*x**5. + 15.*x**4. - 20.*x**32 + 15.*x**2. - 6.*x + 1.


def study6(a=0.005, n=50):
    xlim = [1.-a, 1.+a]
    x = np.linspace(1.-a, 1.+a, num=n)
    yf = []
    ye = []
    err = []
    for i in x:
        yf.append(factored6(i))
        ye.append(expanded6(i))
        err.append(abs((factored6(i)-expanded6(i))/expanded6(i)))
    ylim = [min(min(yf), min(ye)), max(max(yf), max(ye))]
    ax1 = plt.subplot(121)
    ax1.plot(x, yf, 'b', x, ye, 'g')
    ax1.set_xlim(xlim)
    ax1.set_ylim(ylim)
    ax1.set_title("Graph of $f(x)=(x-1)^6$ calculated\nin factored and expanded form\n")
    
    ax2 = plt.subplot(122)
    ax2.plot(x, err, 'r')
    ax2.set_xlim(xlim)
    ax2.set_ylim([min(err), max(err)])
    ax2.set_title("Graph of relative error\nin $f(x)$ calculations\n")
    plt.draw()

a_slider = widgets.widget_float.FloatSlider(min=0., max=0.01, step=0.001, value=0.005)
n_slider = widgets.widget_int.IntSlider(min=1, max=101, step=1, value=50)
i = widgets.interact(study6, a=a_slider, n=n_slider)
display(i)