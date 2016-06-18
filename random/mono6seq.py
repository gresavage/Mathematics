
# coding: utf-8

# In[1]:

import numpy as np
import matplotlib.pyplot as plt
from math import *
from IPython.html import widgets
from IPython.html.widgets import interact
from IPython.display import display
get_ipython().magic(u'matplotlib inline')


# In[2]:

def mono6seq(n=20, **kwargs):
    x = [11./2., 61./11.]
    for i in range(1, n):
        x.append(111.-(1130. - 3000./x[i-1])/x[i])
    return x


# In[3]:

def plot6seq(n=20, **kwargs):
    x = mono6seq(n)
    ax = plt.subplot(111)
    ax.plot(x, 'b')
    ax.axhline(6, linestyle="--", color='violet')
    boom = plt.gca()
    ax.set_title("Plot of first %r values of monotonic sequence" %n)
    ax.annotate(s=r"$lim_{n \rightarrow \infty}=6$", fontsize=18, xy=(plt.xlim()[1]/2., 6), xytext=(0.25, 0.5), textcoords="axes fraction", arrowprops=dict(arrowstyle="->", connectionstyle="arc"))
    plt.draw()
    


# In[4]:

i = interact(plot6seq, n=(2, 100))
display(i)


# In[ ]:



