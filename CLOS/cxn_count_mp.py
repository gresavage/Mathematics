__author__ = 'Tom'

import math
import multiprocessing as mp

def nCr(n, r):
    return round(math.exp(math.lgamma(n+1)-math.lgamma(n-r+1)-math.lgamma(r+1))).__int__()

def nPr(n, r):
    return round(math.exp(math.lgamma(n+1)-math.lgamma(n-r+1))).__int__()

def accelAsc(n):
    a = [0 for i in range(n + 1)]
    k = 1
    y = n - 1
    while k != 0:
        x = a[k - 1] + 1
        k -= 1
        while 2*x <= y:
            a[k] = x
            y -= x
            k += 1
        l = k + 1
        while x <= y:
            a[k] = x
            a[l] = y
            yield a[:k + 2]
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
        yield a[:k + 1]

# a = accelAsc(3)
# while True:
#     try:
#         print next(a)
#     except StopIteration:
#         break

num_outs = 128
num_ins = 128
num_mid_switches = 15
num_input_switches = 16
num_output_switches = 16

beta_max = min(num_outs, num_mid_switches)
# beta = 0                            # Total number of connections, a network parameter

poss_cxn = 0

for beta in xrange(num_outs+1):
    beta_partition = accelAsc(beta)
    mini_sum = 0
    mini_track = 0
    while True:
        try:
            partition_list = next(beta_partition)
        except StopIteration:
            break
        num_outs_on = 0 if partition_list.__contains__(0) else len(partition_list)
        num_mids_on = num_outs_on
        # print "Number of middle switches Being Used: ", beta
        # print "Partition: ", partition_list, " middle switches active"
        # print "num_outs_on: ", num_outs_on
        # Initialize loop variables
        prod = 1                # Product to count possibilities with current configuration
        count = 0               # iteration counter
        part_track = num_outs   # Loop variable to track remaining available outs
        if num_mids_on <= num_mid_switches:
            for i in partition_list:
                # Add cxn possibilites
                # print "#MidSwitch: %r, Count: %r" % (num_mid_switches, count)
                # print "nCr(Nms - count): %r, nCr(part_track): %r" % (nCr(num_mid_switches - count, 1), nCr(part_track, i))
                # print "nCr(%r, 1)*nCr(%r, %r): " % (num_mid_switches - count, part_track, i), nCr(num_mid_switches - count, 1)*nCr(part_track, i)
                prod *= nCr(num_mid_switches - count, min(num_outs_on, 1))*nCr(part_track, i)

                # Update loop variables
                part_track -= i         # decrease the number of available outputs
                count += 1              # increment loop counter
        # else:
        #     print "-----SKIPPED-----"
        input_partition = accelAsc(num_input_switches)          # Partitions of possible input source configurations
        part_sum = 0                                            # Sum of combinations of input source configurations
        part_prod = 1                                           # Number of combinations of input source configurations for the current partition
        while True:
            try:
                input_partition_list = next(input_partition)
                # print "Input partition list: ", input_partition_list
                # print "Input partition list length: ", len(input_partition_list)
                if len(input_partition_list) <= min(beta, num_input_switches):
                    # print "boobs"
                    input_part_track = num_input_switches
                    for i in input_partition_list:
                        # Account for different configurations of input sources
                        part_prod *= nCr(input_part_track, i)
                        # print "nCr(%r, %r): " % (input_part_track, i), nCr(input_part_track, i)
                        input_part_track -= i
                    part_sum += part_prod
                #     print "Partition Product: ", part_prod
                #     print "Partition Sum: ", part_sum
                # else:
                #     print "-----SKIPPED-----"
            except StopIteration:
                break
        # print "Part Sum: ", part_sum
        # print "Prod: ", prod
        # print part_sum*prod
        poss_cxn += part_sum*prod
        # print


            # print "part track: ", part_track
            # print "number of outputs in use, num outs: ", part_track, num_outs
            # print "Partition, sum of prev parts, num_outs-sum of prev parts: ", i, part_track, num_outs-part_track
        # print beta, prod*(num_input_switches**beta)*nCr(num_mid_switches, beta)
        # poss_cxn += prod*(num_input_switches**beta)*nCr(num_mid_switches, beta)
        mini_sum += part_sum*prod
        mini_track += 1
        # if mini_track == beta or num_outs_on == 0:
            # print "Outputs in Use: ", beta
            # print "Mini Sum: ", mini_sum, "\n"
        # print prod*(num_input_switches**beta)
        # poss_cxn += prod*(num_input_switches**beta)*nCr(num_ins, num_outs_on)
        # print poss_cxn
print poss_cxn
print math.log10(poss_cxn)
print float(poss_cxn)
#
#
# n1  = 8             # Number of internal connections at each input switch
# n2  = 8             # Number of internal connections at each output switch
# r1  = 16             # Number of input switches
# r2  = 16             # Number of output switches
# x   = 3             # Maximum number of middle switches used to establish multicast connection
# # m   = (x-1)*n1+n2   # Number of middle stage links
# m = 2*n1-1
#
# print n1, n2, r1, r2, m
# # r = 16
# # m = 31
# # n = 31
# in_combin = 0
# mid_combin = 0
#
# total_paths = n1*n2*r1*r2*m
# print total_paths
# print m
# for i in range(r2+1):
#     # print i
#     # print float(nCr(total_paths, i))
#     # print nCr(r2, i)
#     # print float((n1*r1*m).__pow__(i))
#     # print float(nPr(n1*r1*m, i))
#     # print float((n1*r1*m).__pow__(i))/float(nPr(n1*r1*m, i))
#     # mid_combin += float(nCr(n1*r1*m*r2*n2, r2-i))
#     # mid_combin += float((n1*r1*m).__pow__(i))
#     # in_combin += nCr(total_paths, i)/float(mid_combin)
#     in_combin += nCr(total_paths, i)
#     # print "mid_combin: ", mid_combin
#     # print "in_combin: ", in_combin, '\n'
# print m*r2*n2
# # for i in range(32):
# #     in_combin += nCr(31, i)
# #     print i
# #
# # for i in range(17):
# #     mid_combin += nCr(16, i)
# #     print i
# #
# # combin = r*m*mid_combin*(in_combin**2.0)
#
# combin = in_combin
#
# print "Number of configurations ignoring impossibilities: ", float(combin)
# fivedays = float(30*24*60*60)
# print "Pheasible configurations/Sec: ", poss_cxn/fivedays, "tests/sec"
#
# a, b = 4, 2
# c = a
# print nCr(a, b)
# for i in range(a-b, a):
#     c *= i
# print c
# print "r!:", math.factorial(b)
# print "n!/r!: ", math.factorial(a)/math.factorial(b)
#
# print combin > poss_cxn__author__ = 'Tom'
