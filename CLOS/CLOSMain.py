__author__ = 'Tom'

# Fix inconsistencies between port and address conventions (port i of N vs. switch j, node k)
# Implement deque for fifo (insert(0, x), pop(0) --> appendleft(), popleft()


import sys
import os
import math
import re
import time
import profile
from termcolor import colored
import itertools
from numpy import *
from pylab import *
from collections import deque
import matplotlib.pyplot as plt
from custom_library import *

def up_to_(n1, r1):
    for i in range(n1*r1):
        # print i/r1
        # print test_network.multicast_generator(i/r1, i%r1, [i/r2], 0)
        test_network.multicast_generator(i/n1, i % n1, [i/r2], 0)

def nCr(n, r):
    f = math.factorial
    return f(n)/(f(r)*f(n-r))

def cxnplot(switch_obj, lines=None, figureNumber=None):
    # Plot the connections. inX is the x coordinates for the input nodes which are at the heights inY.
    # nin is the number of input nodes, nout is the number of output nodes. inodecolors is a list or tuple containing
    # color data for each input node. Useful if trying to view the way connections change, and cxnplot has to be called
    # many times

    xin = 0                                                         # Left endpoint of domain
    xout = switch_obj.size[1]                                       # Reft endpoint of domain

    # Creates a new figure if there is none.
    if figureNumber is None:
        figureNumber=time.time()                                    # Create a unique number for the figure
        cxnfig = plt.figure(num=figureNumber, figsize=(10, 10))     # Assign the number to the figure
        ax = cxnfig.add_subplot(111)
        # inodecolors = cxncolor()                                    # Create a class to contain the colors
        ax.scatter([0 for i in range(switch_obj.size[0])], [i for i in range(switch_obj.size[0])])                                        # Make markers for the nodes
        ax.scatter([switch_obj.size[1] for i in range(switch_obj.size[1])], [i for i in range(switch_obj.size[1])])
        lines = []

        for i in switch_obj.connections:
            if switch_obj.connections[i] is not None:
                for j in set(switch_obj.connections[i]):
                    lines.append(ax.plot([xin, xout], [i, j], color=switch_obj.colors[i]))
        plt.axis([0, switch_obj.size[0], -1, switch_obj.size[1]])
    else:
        cxnfig = plt.figure(figureNumber)
        cxnfig.clear()
        ax = cxnfig.gca()
        ax.scatter([0 for i in range(switch_obj.size[0])], [i for i in range(switch_obj.size[0])])
        ax.scatter([switch_obj.size[1] for i in range(switch_obj.size[1])], [i for i in range(switch_obj.size[1])])
        for i in switch_obj.connections:
            for j in set(switch_obj.connections[i]):
                ax.plot([0, switch_obj.size[1]], [i, j], color=switch_obj.colors[i])
        plt.axis([0, switch_obj.size[0], -1, switch_obj.size[1]])
    return figureNumber, lines

def rand_packet(in_size, out_size):
    return packet(random.randint(0, in_size), random.randint(0, out_size))

class clos_switch:
    # link refers to the input link connection state. For example, if inlink[2,2]=True, then the input from node 3
    # (since indexing starts at 0) is connected to outlet node 3
    # Each row of the link matrix denotes which outlet nodes that input is connected to. Then each column must contain
    # at most one True value.
    # Declare switch size and type upon instantiation
    def __init__(self, inputs, outputs, switchtype, switch_num):
        self.type = switchtype.lower()
        self.link = zeros((inputs, outputs), bool)
        self.shape = self.link.shape
        self.size = [self.shape[0], self.shape[1]]
        self.id = switch_num
        self.connections = dict()                                                       # Dictionary containing active connections at each input node
        # for i in range(self.size[0]):
        #     self.connections.setdefault(i, [])
        self.new_color()

    def __iter__(self):
        return iter(self.link)

    def new_color(self):
        # Randomly generates a new color
        self.colors = [(random(), random(), random()) for i in range(self.size[0])]           # List containing color information for each node

    def inherit_color(self, input, cxn_color):
        # Forces the input node to inherit the cxn_color
        self.colors[input] = cxn_color

    def clear_switch(self):
        for i, j in range(self.size[0], self.size[1]):
            self.link[i, j] = False

    def connect(self, inputs=int, outputs=(list, tuple)):
        # Method to establish connection from input(s) to output(s)
        if isinstance(inputs, (tuple, list)):
            # print "input tuple"
            for i in range(len(inputs)):
                if outputs[i] in self.taken_outs(set):
                    print "Error in %r switch %r: Output node (%r) already occupied." % (self.type, self.id, outputs[i])
                    print "Connection from input (%r) to output (%r) not established." % (inputs[i], outputs[i])
                else:
                    self.link[inputs[i], outputs[i]] = True
                    self.connections.update({inputs[i]: self.connected_to(inputs[i])})
        elif isinstance(inputs, (int, float, long)):
            # print "input int"
            if isinstance(outputs, (list, tuple)):
                # print "output tuple"
                for i in outputs:
                    if i in self.taken_outs(set):
                        print "Error in %r switch %r: Output node (%r) already occupied." % (self.type, self.id, i)
                        print "Connection from input (%r) to output (%r) not established." % (inputs, i)
                    else:
                        self.link[inputs, i] = True
                        # self.connections.update({inputs: [self.connected_to(inputs)]})
            elif isinstance(outputs, (int, float, long)):
                # print "output int"
                if outputs in self.taken_outs(set):
                    print "Error in %r switch %r: Output node (%r) already occupied." % (self.type, self.id, outputs)
                    print "Connection from input (%r) to output (%r) not established." % (inputs, outputs)
                else:
                    self.link[inputs, outputs] = True
            self.connections.update({inputs: self.connected_to(inputs)})
        # if isinstance(outputs, (list, tuple)):
        #     # print "output tuple"
        #     for i in outputs:
        #         if i in self.taken_outs(set):
        #             print "Error: Output node (%r) already occupied." % i
        #             print "Connection from input (%r) to output (%r) not established." % inputs, i
        #         else:
        #             print i
        #             self.link[inputs, i] = True
        #             self.connections.update({inputs: [self.connected_to(inputs)]})
        # elif isinstance(outputs, (int, float, long)):
        #     # print "output int"
        #     if outputs in self.taken_outs(set):
        #         print "Error: Output node (%r) already occupied.\n Connection from input (%r) to output (%r) not established."
        #     else:
        #         self.link[inputs, outputs] = True
        #         self.connections.update({inputs: self.connected_to(inputs)})
        # self.update()
        # self.update()

    def disconnect(self, inputs, outputs):
        # Method to disconnect input(s) from output(s)
        if isinstance(inputs, (tuple, list)):
            # print "input tuple"
            for i in range(len(inputs)):
                if outputs[i] not in self.taken_outs(set):
                    print "Error %r switch %r: Output node (%r) has no connection." % (self.type, self.id, outputs[i])
                    print "Connection from input (%r) to output (%r) not changed." % (inputs[i], outputs[i])
                elif inputs[i] not in self.taken_ins(set):
                    print "Error in %r switch %r: Input node (%r) has no connecton." % (self.type, self.id, outputs[i])
                    print "Connection from input (%r) to output (%r) not changed." % (inputs[i], outputs[j])
                else:
                    self.link[inputs[i], outputs[i]] = False
                    self.connections.update({inputs[i]: self.connected_to(inputs[i])})
        elif isinstance(inputs, (int, float, long)):
            # print "input int"
            if isinstance(outputs, (list, tuple)):
                # print "output tuple"
                for i in range(len(outputs)):
                    if outputs[i] not in self.taken_outs(set):
                        print "Error in %r switch %r: Output node (%r) has no connecton." % (self.type, self.id, outputs[i])
                        print "Connection from input (%r) to output (%r) not changed." % (inputs, outputs[i])
                    elif inputs not in self.taken_ins(set):
                        print "Error in %r switch %r: Input node (%r) has no connecton." % (self.type, self.id, outputs[i])
                        print "Connection from input (%r) to output (%r) not changed." % (inputs, outputs[i])
                    else:
                        self.link[inputs, outputs[i]] = False

            elif isinstance(outputs, (int, float, long)):
                # print "output int"
                if outputs not in self.taken_outs(set):
                    print "Error in %r switch %r: Output node (%r) has no connecton." % (self.type, self.id, outputs)
                    print "Connection from input (%r) to output (%r) not changed." % (inputs, outputs)
                elif inputs not in self.taken_ins(set):
                    print "Error in %r switch %r: Output node (%r) has no connecton." % (self.type, self.id, outputs)
                    print "Connection from input (%r) to output (%r) not changed." % (outputs)
                else:
                    self.link[inputs, outputs] = False
            self.connections.update({inputs: self.connected_to(inputs)})
                    # self.connections[inputs] = self.connections[inputs].remove(outputs)

    def connected_to(self, node, node_type="input", return_type=list):
        # This method returns a set, list, or tuple, depending upon returntype, denoting which of the inputs have connections
        # If the node_type is "input", the method returns the list of output nodes to which it is connected
        # If the node_type is anything else, the method returns the list of input nodes to which the output is connected
        taken = set()
        if node_type is "input":
            for i in range(self.shape[1]):
                if self.link[node, i]:
                    taken.add(i)
        else:
            for i in range(self.shape[0]):
                if self.link[i, node]:
                    taken.add(i)
        if return_type is list:
            return list(taken)
        elif return_type is set:
            return taken
        elif return_type is tuple:
            return tuple(taken)

    def taken_outs(self, return_type=list):
        # returns a list of indeces referring to the taken output nodes
        taken = set()
        for i in range(self.size[1]):
            if self.link[:, i].any():
                taken.add(i)
        if return_type is list:
            return list(taken)
        elif return_type is set:
            return taken
        elif return_type is tuple:
            return tuple(taken)

    def taken_ins(self, return_type=list):
        # returns a list of indeces referring to the taken output links
        taken = set()
        for i in range(self.size[0]):
            if self.link[i, :].any():
                taken.add(i)
        if return_type is list:
            return list(taken)
        elif return_type is set:
            return taken
        elif return_type is tuple:
            return tuple(taken)

    def free_outs(self, return_type=list):
        # returns a list of indeces referring to the available output links
        free = set()
        for i in range(self.size[1]):
            if not self.link[:, i].any():
                free.add(i)
        if return_type is list:
            return list(free)
        elif return_type is set:
            return free
        elif return_type is tuple:
            return tuple(free)

    def free_ins(self, return_type=list):
        # returns a list of indeces referring to the available output links
        free = set()
        for i in range(self.size[0]):
            if not self.link[i, :].any():
                free.add(i)
        if return_type is list:
            return list(free)
        elif return_type is set:
            return free
        elif return_type is tuple:
            return tuple(free)

    def connection_data(self):
        for i in range(self.size[0]):
            self.connections.update({i: self.connected_to(i)})      # Set the value of the node-key to a list containing all the outputs to which it is connected
        return self.connections

class packet:
    def __init__(self, port, destination):
        self.time = time.time() + datetime.time.microsecond/1000000.0
        self.port = port
        self.destination = destination
        self.data = "Hello World: %r" % random()


class clos_network:
    def __init__(self, n1, n2, r1, r2, x, m, threshhold, cell_time):
        self.input_switch_nodes = n1
        self.output_switch_nodes = n2
        self.num_input_switches = r1
        self.num_output_switches = r2
        self.mid_switch_max = x
        # self.num_mid_switches = (self.mid_switch_max-1)*self.input_switch_nodes+self.output_switch_nodes
        self.num_mid_switches = m
        self.is_non_blocking = self.num_mid_switches >= self.input_switch_nodes + self.output_switch_nodes
        # self.cell_size = cell_size                          # transferable data length
        self.cell_time = cell_time
        self.threshhold = threshhold                        # r, number of time slots needed for packet scheduler to schedule cells/frame size
        self.frame_period = self.cell_time*self.threshhold  # Frame period (length of time for connection frame???
        self.paths = []                                     # list containing tuples whose entries specify the existing connection paths (input_switch, middle_switch, output_switch)
        # self.ILC = [{(j, []) for j in range(self.output_switch_nodes*self.num_output_switches)} for i in range(self.input_switch_nodes*self.num_input_switches)]
        self.voq = {(j, k):{i: [] for i in range(self.num_output_switches)} for j in range(self.num_input_switches) for k in range(self.input_switch_nodes)}
        # self.raq = {i: [] for i in range(self.num_input_switches)}
        self.packet_queue = []

        self.input_switches = [clos_switch(self.input_switch_nodes, self.num_mid_switches, "input", i) for i in range(self.num_input_switches)]   # Input switch layer array
        self.middle_switches = [clos_switch(self.num_input_switches, self.num_output_switches, "middle", i) for i in range(self.num_mid_switches)]    # Middle switch layer array
        self.output_switches = [clos_switch(self.num_mid_switches, self.output_switch_nodes, "output", i) for i in range(self.num_output_switches)]      # Output switch layer array)
        self.outbox = {i: [] for i in range(self.num_output_switches)}
        self.terminate = False                              # Terminate network operation

        # Neural rearrangement algorithm parameters
        self.alpha = 5
        self.a = 3
        self.b = 2
        self.c = 1

    def port_connect(self, packet_=packet):
        input_switch, input_switch_port = port_to_address(packet)
        self.multicast_generator(input_switch, input_switch_port, packet.destination, packet_.data)

    def call_requests(self, input_switch, output_switch):
        # Number of ongoing calls between input_switch and output_switch
        tt = 0
        for j in range(self.num_mid_switches):
            if (input_switch, j, output_switch) in self.paths:
                tt += 1
        # for i in self.voq[(input_switch, output_switch)]:
        #     tt += self.voq[(input_switch, output_switch)][i].length()
        return tt

    def add_to_paths(self, path):
        self.paths.append(tuple(path))

    def remove_from_paths(self, path):
        self.paths.remove(tuple(path))

    def port_to_address(self, packet):
        input_switch_port = packet.port % self.num_input_switches
        input_switch = int(packet.port/self.num_input_switches)
        return input_switch, input_switch_port

    def address_to_port(self, input_switch, input_switch_port):
        return input_switch*self.num_input_switches + input_switch_port

    def transfer(self):
        # transfer data to outputs.
        do = 0

    def set_path(self, input_switch_num, middle_switch_num, output_switch_num, value):
        if value:
            self.input_switches[input_switch_num].connect(self.input_switches[input_switch_num].free_ins()[0], middle_switch_num)
            self.middle_switches[middle_switch_num].connect(input_switch_num, output_switch_num)
            self.output_switches[output_switch_num].connect(middle_switch_num, self.output_switches[output_switch_num].free_outs()[0])
            self.add_to_paths((input_switch_num, middle_switch_num, output_switch_num))
        elif not value:
            self.input_switches[input_switch_num].disconnect(self.input_switches[input_switch_num].connected_to(middle_switch_num, "output"), middle_switch_num)
            self.middle_switches[middle_switch_num].disconnect(input_switch_num, output_switch_num)
            self.output_switches[output_switch_num].disconnect(middle_switch_num, self.output_switches[output_switch_num].connected_to(middle_switch_num))
            self.remove_from_paths((input_switch_num, middle_switch_num, output_switch_num))

    def send(self, packet_):
        self.packet_queue.append(packet_)

    def operate(self):          # Normal network operation
        while True:
            if length(self.packet_queue) is not 0:
                packet_ = self.packet_queue.pop(0)
                input_switch_num, input_switch_port = self.port_to_address(packet_)
                if not self.multicast_generator(input_switch_num, input_switch_port, packet_.destinations, (packet_.time, packet_.data)):
                    NeuroSwitch(input_switch_num, packet_.destinations, self.alpha, self.a, self.b, self.c, self)
            if self.terminate:
                break
            # for i in self.voq.keys():
            #     for j in range(self.num_output_switches):
            #         if len(self.voq[i][j]) is not 0:
            #             self.multicast_generator(i[0], i[1], j, self.voq[i][j].pop(0))
            #             self.outbox[j].append()

    def multicast_generator(self, input_switch_num, input_node, output_switch_nums, packet_info):
        outs_set = set(output_switch_nums)      # Set containing all desired outputs
        outs_int = set(output_switch_nums)      # Set containing the remaining desired outputs
        union = set()                           # Set to contain the desired outputs to which the middle switches are able to connect
        mid_pairs = dict()                      # Dictionary whose keys are middle switch indeces, and values are the desired outputs to which the middle switch can connect
        mcast_est = False                       # Was the multicast established
        # combine_iter = itertools.combinations([i for i in range(self.num_mid_switches)], self.mid_switch_max)       # Iterator object containing all of the combinations of middle switch nodes up to self.mid_switch_max as tuples
        # free_outs_set = current_switch.free_outs(set)
        middle_switch_intersection = outs_set

        for g in range(1, self.mid_switch_max):
            combine_iter = itertools.combinations([i for i in range(self.num_mid_switches)], g)     # creates combinations of middle switch indeces up to length g
            while True:
                try:
                    indeces = combine_iter.next() # get the next middle switch index combination
                except StopIteration:
                    break

                # print indeces
                index_counter = 0
                for i in indeces:
                    # print "index = ", i
                    if len(outs_int & self.middle_switches[i].free_outs(set)) is not 0:             # if there are commonalities between the desired outputs and available outs from the ith middle switch add the ith middle switch and its relevant outs to mid_pairs
                        mid_pairs.update({i: list(outs_int & self.middle_switches[i].free_outs(set))})  #
                        # print colored("mid_pairs = ", "green"), mid_pairs
                    union |= self.middle_switches[i].free_outs(set)     # add the free outputs to the union
                    outs_int -= self.middle_switches[i].free_outs(set)  # remove the outputs from the list of desired outs
                    # print colored("free_outs = ", "yellow"), self.middle_switches[i].free_outs(set)
                    # print colored("union = ", "red"), union
                    # print colored("outs_int = ", "blue"), outs_int, "\n"

                if outs_set.issubset(union):        # if the set of desired outs is a subset of the union of available outs from the combination of middle switches then create the connection
                    # print "bleh"
                    # print mid_pairs
                    print input_switch_num
                    print input_node
                    print self.num_input_switches
                    cxn_color = self.input_switches[input_switch_num].colors[input_node]        # set the connection color for plotting
                    # print cxn_color
                    # print "mid_pairs.keys() = ", mid_pairs.keys()
                    self.input_switches[input_switch_num].connect(input_node, mid_pairs.keys())     # connect the input switch to the middle switches
                    # print "bauble"
                    for i in mid_pairs.keys():
                        self.middle_switches[i].connect(input_switch_num, mid_pairs[i])
                        self.middle_switches[i].inherit_color(input_switch_num, cxn_color)
                        # print self.middle_switches[i].colors[input_switch_num]
                        for j in mid_pairs[i]:
                            # print "mid_pairs[i] = ", mid_pairs[i]
                            # print "free ins %r" % j, self.output_switches[j].free_ins()
                            # print "free outs %r" % j, self.output_switches[j].free_outs()
                            # print "j= ", j
                            # print "num_output_switches ", self.output_switches.__len__()
                            # print "num free outs: ", self.output_switches[j].free_outs()
                            try:
                                self.output_switches[j].inherit_color(i, cxn_color)
                                self.output_switches[j].connect(i, self.output_switches[j].free_outs()[0])
                                self.add_to_paths((input_switch_num, i, j))
                            except IndexError:
                                break

                            # self.voq[(input_switch_num, input_node)][j].pop(0)          # Remove the request from the VOQ
                            # self.outbox[j] = packet_info                                # Collect data
                            # print "middle_switch = ", i
                            # print "output_switch = ", j, "\n"
                            # self.output_switches[j].inherit_color(i, cxn_color)
                            # print self.output_switches[j].colors[i]
                    mcast_est = True
                    break
                index_counter += 1

                if index_counter & g == 0:      # reset initial values for the next combination of middle switches
                    union = set()
                    outs_int = set(output_switch_nums)
                    mid_pairs = dict()
            if mcast_est:
                # if the connection was established, exit the loop
                break
        return mcast_est

    def stage_connection_plot(self, stage, xin, xout, yin, yout, figureNumber):
        cxnfig = plt.figure(num=figureNumber, figsize=(10, 10))     # Assign the number to the figure
        ax = cxnfig.gca()
        lines = []
        for i in range(len(stage)):
            for j in stage[i].connections.keys():
                for k in stage[i].connections[j]:
                    lines.append(ax.plot([xin, xout], [i*(yin[j+stage[i].size[0]]-yin[j])+yin[j+1], i*(yout[j+stage[i].size[1]]-yout[j])+yout[k+1]], color=stage[i].colors[j]))
        return lines

    def network_plot(self, lines=None, figureNumber=None):
    # Plot the connections. inX is the x coordinates for the input nodes which are at the heights inY.
    # nin is the number of input nodes, nout is the number of output nodes. inodecolors is a list or tuple containing
    # color data for each input node. Useful if trying to view the way connections change, and cxnplot has to be called
    # many times
        xin = 0                                                         # Left endpoint of domain
        xout = 10                                                       # Reft endpoint of domain

        yin = 0
        ymax = max([self.num_input_switches, self.num_mid_switches, self.num_output_switches])     # 8
        input_switch_height = (ymax)/(self.num_input_switches + 2)
        mid_switch_height = (ymax)/(self.num_mid_switches + 2)
        output_switch_height = (ymax)/(self.num_output_switches + 2)

        input_y = linspace(yin, ymax, self.num_input_switches + 2)
        input_single_switch_in_y = linspace(0, input_switch_height, self.input_switch_nodes + 1)
        input_in_y = linspace(yin, ymax, self.input_switch_nodes*self.num_input_switches + 2)
        input_out_y = linspace(yin, ymax, self.num_mid_switches*self.num_input_switches + 2)
        input_single_switch_out_y = linspace(0, input_switch_height, self.num_mid_switches + 1)

        middle_y = linspace(yin, ymax, self.num_mid_switches + 2)
        middle_single_switch_in_y = linspace(0, mid_switch_height, self.num_input_switches + 1)
        middle_in_y = linspace(yin, ymax, self.num_input_switches*self.num_mid_switches + 2)
        middle_out_y = linspace(yin, ymax, self.num_output_switches*self.num_mid_switches + 2)
        middle_single_switch_out_y = linspace(0, mid_switch_height, self.num_output_switches + 1)

        output_y = linspace(yin, ymax, self.num_output_switches + 2)
        output_single_switch_in_y = linspace(0, output_switch_height, self.num_mid_switches + 1)
        output_in_y = linspace(yin, ymax, self.num_mid_switches*self.num_output_switches + 2)
        output_out_y = linspace(yin, ymax, self.output_switch_nodes*self.num_output_switches + 2)
        output_single_switch_out_y = linspace(0, output_switch_height, self.num_output_switches + 1)


        input_size = (input_y[1]-input_y[0])
        middle_size = middle_y[1]-middle_y[0]
        output_size = output_y[1]-output_y[0]

        stages_x = linspace(xin, xout, 8)

        # Creates a new figure if there is none.
        if figureNumber is None:
            figureNumber=time.time()                                    # Create a unique number for the figure
            cxnfig = plt.figure(num=figureNumber, figsize=(10, 10))     # Assign the number to the figure
            ax = cxnfig.add_subplot(111)
            # for j in range(self.num_input_switches):
            #     ax.scatter([stages_x[1] for i in range(self.input_switch_nodes)], [self.num_input_switches*j + input_single_switch_in_y[k] for k in range(self.input_switch_nodes)])
            #     ax.scatter([stages_x[2] for i in range(self.num_mid_switches)], [self.num_input_switches*j + input_single_switch_out_y[k] for k in range(self.num_mid_switches)])
            # for j in range(self.num_mid_switches):
            #     ax.scatter([stages_x[3] for i in range(self.num_input_switches)], [j + middle_single_switch_in_y[k] for k in range(self.num_input_switches)])
            #     ax.scatter([stages_x[4] for i in range(self.num_output_switches)], [j + middle_single_switch_out_y[k] for k in range(self.num_output_switches)])
            # for j in range(self.num_output_switches):
            #     ax.scatter([stages_x[5] for i in range(self.num_mid_switches)], [self.num_output_switches*j + output_single_switch_in_y[k] for k in range(self.num_mid_switches)])
            #     ax.scatter([stages_x[6] for i in range(self.output_switch_nodes)], [self.num_output_switches*j + output_single_switch_out_y[k] for k in range(self.output_switch_nodes)])
            ax.scatter([stages_x[1] for i in range(self.input_switch_nodes*self.num_input_switches)], [input_in_y[i+1] for i in range(self.input_switch_nodes*self.num_input_switches)])                                       # Make markers for the nodes
            ax.scatter([stages_x[2] for i in range(self.num_mid_switches*self.num_input_switches)], [input_out_y[i+1] for i in range(self.num_mid_switches*self.num_input_switches)])
            ax.scatter([stages_x[3] for i in range(self.num_mid_switches*self.num_input_switches)], [middle_in_y[i+1] for i in range(self.num_mid_switches*self.num_input_switches)])
            ax.scatter([stages_x[4] for i in range(self.num_mid_switches*self.num_output_switches)], [middle_out_y[i+1] for i in range(self.num_mid_switches*self.num_output_switches)])
            ax.scatter([stages_x[5] for i in range(self.num_mid_switches*self.num_output_switches)], [output_in_y[i+1] for i in range(self.num_mid_switches*self.num_output_switches)])
            ax.scatter([stages_x[6] for i in range(self.output_switch_nodes*self.num_output_switches)], [output_out_y[i+1] for i in range(self.output_switch_nodes*self.num_output_switches)])
            input_lines = []
            middle_lines = []
            output_lines = []
            # self.stage_connection_plot(self.input_switches, stages_x[1], stages_x[2], input_in_y, input_out_y, figureNumber)
            # self.stage_connection_plot(self.middle_switches, stages_x[3], stages_x[4], middle_in_y, middle_out_y, figureNumber)
            # self.stage_connection_plot(self.output_switches, stages_x[5], stages_x[6], output_in_y, output_out_y, figureNumber)
            # plt.axis([xin, xout, yin - 1, yout])
        else:
            cxnfig = plt.figure(figureNumber)
            cxnfig.clear()
            # ax = cxnfig.gca()
            # ax.scatter([0 for i in range(switch_obj.size[0])], [i for i in range(switch_obj.size[0])])
            # ax.scatter([switch_obj.size[1] for i in range(switch_obj.size[1])], [i for i in range(switch_obj.size[1])])
            # for i in switch_obj.connections:
            #     for j in set(switch_obj.connections[i]):
            #         ax.plot([0, switch_obj.size[1]], [i, j], color=switch_obj.colors[i])
            # plt.axis([0, switch_obj.size[0], -1, switch_obj.size[1]])
        self.stage_connection_plot(self.input_switches, stages_x[1], stages_x[2], input_in_y, input_out_y, figureNumber)
        self.stage_connection_plot(self.middle_switches, stages_x[3], stages_x[4], middle_in_y, middle_out_y, figureNumber)
        self.stage_connection_plot(self.output_switches, stages_x[5], stages_x[6], output_in_y, output_out_y, figureNumber)
        plt.axis([xin, xout, yin - 1, ymax+1])
        plt.show()
        return figureNumber, lines

n1  = 4             # Number of internal connections at each input switch
n2  = 4             # Number of internal connections at each output switch
r1  = 6             # Number of input switches
r2  = 6             # Number of output switches
x   = 4             # Maximum number of middle switches used to establish multicast connection
m   = (x-1)*n1+n2   # Number of middle stage links

fail = 0

test_network = clos_network(n1, n2, r1, r2, x, m, 1, 1)

profile.run("up_to_(n1, r1)")
# def up_to_(n1, r1):
#     for i in range(n1*r1):
#     # print i/r1
#     # print test_network.multicast_generator(i/r1, i%r1, [i/r2], 0)
#     test_network.multicast_generator(i/r1, i%r1, [i/r2], 0)

for i in range(n1*r1):
    # print i/r1
    # print i % r1
    # print [r2 - i/r2 - 1]
    _result = test_network.multicast_generator(i/n1, i % n1, [r2 - i/r2 - 1], 0)
    if not _result:
        # print "(node, switch, output_node, output_switch): ( %r, %r, %r, %r)" %(i % r1, i/r1, i % r2, r2 - i/r2)
        fail+=1
print fail
test_network.network_plot()



# # test_network.multicast_generator(0, [1, 2, 3])
# # for k in range(r1):
# #     test_network.middle_switches[r1].connect(1, [i for i in range(n1)])
# for k in range(m):
#     test_network.middle_switches[k].connect(1, [1, 2, 3])
#     # print test_network.middle_switches[k].connections
# test_network.middle_switches[0].disconnect(1, 1)        # Example 2, 3
# # test_network.middle_switches[0].disconnect(1, (1, 2, 3)) # Example 1
# test_network.middle_switches[1].disconnect(1, 2)        # Example 3
# # test_network.middle_switches[1].disconnect(1, (2, 3))   # Example 2
# test_network.middle_switches[2].disconnect(1, 3)      # Example 3
# test_network.network_plot()
# # print test_network.middle_switches[0].connection_data()
# # print test_network.middle_switches[1].connection_data()
# # print test_network.middle_switches[2].connection_data()
#
#
# # test_network.middle_switches[0].connect(1, 2)
# # print colored("test_network.input_switches[0].link = ", "blue"), test_network.middle_switches[0].link
# # test_network.multicast(0, [1, 2, 3])
# # print "test_network.output_switches[1].free_outs(): ", test_network.output_switches[1].free_outs()
# # print "test_network.output_switches[2].free_outs(): ", test_network.output_switches[2].free_outs()
# # print "test_network.output_switches[3].free_outs(): ", test_network.output_switches[3].free_outs()
# print "test_network.multicast_generator(0, 1,  [1, 2, 3]) returns ", test_network.multicast_generator(0, 1,  [1, 2, 3])
# test_network.network_plot()
# # print test_network.multicast_generator(0, 2,  [1, 2, 3])
# # test_network.input_switches[0].connect(1, 1)
# # test_network.network_plot()
# # print "test_network.output_switches[1].connections = ", test_network.output_switches[1].connections
# # print "test_network.output_switches[2].connections = ", test_network.output_switches[2].connections
# # print "test_network.output_switches[3].connections = ", test_network.output_switches[3].connectionsn1  = 4             # Number of internal connections at each input switch
# n2  = 4             # Number of internal connections at each output switch
# r1  = 4             # Number of input switches
# r2  = 4             # Number of output switches
# x   = 3             # Maximum number of middle switches used to establish multicast connection
# m   = (x-1)*n1+n2   # Number of middle stage links
#
# test_network = clos_network(n1, n2, r1, r2, x, 1, 1)
# # test_network.multicast_generator(0, [1, 2, 3])
# # for k in range(r1):
# #     test_network.middle_switches[r1].connect(1, [i for i in range(n1)])
# for k in range(m):
#     test_network.middle_switches[k].connect(1, [1, 2, 3])
#     # print test_network.middle_switches[k].connections
# # test_network.middle_switches[0].disconnect(1, 1)        # Example 2, 3
# # test_network.middle_switches[0].disconnect(1, (1, 2, 3)) # Example 1
# # test_network.middle_switches[1].disconnect(1, 2)        # Example 3
# # test_network.middle_switches[1].disconnect(1, (2, 3))   # Example 2
# # test_network.middle_switches[2].disconnect(1, 3)      # Example 3
# test_network.network_plot()
# # print test_network.middle_switches[0].connection_data()
# # print test_network.middle_switches[1].connection_data()
# # print test_network.middle_switches[2].connection_data()
#
#
# # test_network.middle_switches[0].connect(1, 2)
# # print colored("test_network.input_switches[0].link = ", "blue"), test_network.middle_switches[0].link
# # test_network.multicast(0, [1, 2, 3])
# # print "test_network.output_switches[1].free_outs(): ", test_network.output_switches[1].free_outs()
# # print "test_network.output_switches[2].free_outs(): ", test_network.output_switches[2].free_outs()
# # print "test_network.output_switches[3].free_outs(): ", test_network.output_switches[3].free_outs()
# print "test_network.multicast_generator(0, 1,  [1, 2, 3]) returns ", test_network.multicast_generator(0, 1,  [1, 2, 3])
# test_network.network_plot()
# # print test_network.multicast_generator(0, 2,  [1, 2, 3])
# # test_network.input_switches[0].connect(1, 1)
# # test_network.network_plot()
# # print "test_network.output_switches[1].connections = ", test_network.output_switches[1].connections
# # print "test_network.output_switches[2].connections = ", test_network.output_switches[2].connections
# # print "test_network.output_switches[3].connections = ", test_network.output_switches[3].connections

