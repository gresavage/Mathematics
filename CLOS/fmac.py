__author__ = 'Tom'

import sys
import os
import math
import re
import time
from termcolor import colored
import itertools
from numpy import *
from pylab import *
from CLOSMain import *
import matplotlib.pyplot as plt

class ILC:
    # Input Line Card
    def __init__(self, threshhold, cell_time, num_inputs, num_outputs):
        self.queue = [VOQ(num_outputs) for i in range(num_inputs)]
        self.cell_time = cell_time                          # Time between connections
        self.threshhold = threshhold                        # r, number of time slots needed for packet scheduler to schedule cells/frame size
        self.frame_period = self.cell_time*self.threshhold  # Frame period (length of time for connection frame????)

class VOQ:
    # Virtual Output Queue Class
    def __init__(self, num_outputs):
        self.queue = {(j, []) for j in range(num_outputs)}

def fmac_request(network=clos_network):
    network.

