__author__ = 'Tom'

import math

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
def network_states(nin, nout, nis, nos, nms):
    # nin, nout, nis, nos, nms = 4, 4, 2, 2, 3
    num_outs = nout
    num_ins = nin
    num_mid_switches = nms
    num_input_switches = nis
    num_output_switches = nos

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
                # print "Mid switches in use: ", num_mids_on
                # comb_str = "["
                for i in partition_list:
                    # Add cxn possibilites
                    # print "#MidSwitch: %r, Count: %r" % (num_mid_switches, count)
                    # print "nCr(Nms - count): %r, nCr(part_track): %r" % (nCr(num_mid_switches - count, 1), nCr(part_track, i))
                    # print "nCr(%r, 1)*nCr(%r, %r): " % (num_mid_switches - count, part_track, i), nCr(num_mid_switches - count, 1)*nCr(part_track, i)
                    prod *= nCr(num_mid_switches - count, min(num_outs_on, 1))*nCr(part_track, i)
                    # comb_str += "(%r, %r)(%r, %r)" % (num_mid_switches - count, min(num_outs_on, 1), part_track, i)

                    # Update loop variables
                    part_track -= i         # decrease the number of available outputs
                    count += 1              # increment loop counter
                    # comb_str += "]"
            # else:
            #     print "-----SKIPPED-----"
            # input_partition = accelAsc(num_input_switches)          # Partitions of possible input source configurations
            in_part_sum = 0                                            # Sum of combinations of input source configurations

            for j in xrange(max(1, min(num_input_switches, num_outs_on))):
                input_partition = accelAsc(j+1)          # Partitions of possible input source configurations
                part_sum = 0                                            # Sum of combinations of input source configurations
                part_prod = 1                                           # Number of combinations of input source configurations for the current partition
                try:
                    # part_str = "["
                    input_partition_list = next(input_partition)
                    num_ins_on = 0 if input_partition_list.__contains__(0) else len(input_partition_list)
                    # print "Num of outputs in use: ", beta
                    # print "Input partition list: ", input_partition_list
                    # print "Input partition list length: ", num_ins_on
                    # choice_str = str()
                    if num_ins_on <= min(beta, num_input_switches):
                        # print "boobs"
                        input_part_track = num_input_switches
                        for i in input_partition_list:
                            # print i
                            # Account for different configurations of input sources
                            # choice_str += '(%r, %r)' % (input_part_track, i)
                            part_prod *= nCr(input_part_track, i)
                            # print "nCr(%r, %r): " % (input_part_track, i), nCr(input_part_track, i)
                            input_part_track -= i
                            # print
                            # print choice_str
                            # if num_ins_on > num_input_switches - input_part_track:
                            #     part_str += choice_str + " + "
                            # else:
                            #     part_str += choice_str
                    # print "Part product: ", part_prod
                    part_sum += part_prod
                    # print "Part sum: ", part_sum
                    # print part_str + "] * " + comb_str
                    # print part_str + "]*%r" % prod
                    # print
                    #     print "Partition Product: ", part_prod
                    #     print "Partition Sum: ", part_sum
                    # else:
                    #     print "-----SKIPPED-----"
                except StopIteration:
                    break
                in_part_sum += part_sum
                #     print "Part Sum: ", part_sum
                #     print "Input Part Sum: ", in_part_sum
                #     print "Prod: ", prod
                #     print "Iteration portion: ", part_sum*prod
                #     print "-------------------"
                # print "======================================"
                # print

                # print


                # print "part track: ", part_track
                # print "number of outputs in use, num outs: ", part_track, num_outs
                # print "Partition, sum of prev parts, num_outs-sum of prev parts: ", i, part_track, num_outs-part_track
            # print beta, prod*(num_input_switches**beta)*nCr(num_mid_switches, beta)
            # poss_cxn += prod*(num_input_switches**beta)*nCr(num_mid_switches, beta)
            mini_sum += in_part_sum*prod
            poss_cxn += in_part_sum*prod
            mini_track += 1
            # print prod*(num_input_switches**beta)
            # if mini_track == beta or num_outs_on == 0:
            #     print "Outputs in Use: ", beta
            #     print "Mini Sum: ", mini_sum, "\n"
                # poss_cxn += prod*(num_input_switches**beta)*nCr(num_ins, num_outs_on)
                # print poss_cxn
    # print poss_cxn
    # print math.log10(poss_cxn)
    # print float(poss_cxn)
    return float(poss_cxn)
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
    # print combin > poss_cxn

# timeit.timeit('network_states()', number=100)

print network_states(16, 16, 4, 4, 7)/(5*3600*24)

#
# e_pre = 1
# e_list = []
# for i in range(0, 4):
#     e_post = network_states((2**i)*8, (2**i)*8, 8, 8, 2**(i+1)-1)
#     e_list.append(math.log10(e_post/e_pre))
#     print "Rate: ", e_list[i]
#     print "Size: (%r, %r, %r)" % ((2**i)*8, 8, 2**(i+1)-1)
#     print