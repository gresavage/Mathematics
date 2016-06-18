__author__ = 'Tom'

from CLOSMain import *
from custom_library import *



def NeuroSwitch(input_switch, output_switches, alpha, a, b, c, netfabric=clos_network):
    r = max(netfabric.num_input_switches, netfabric.num_output_switches)
    u = zeros((netfabric.num_input_switches, netfabric.num_output_switches, netfabric.num_mid_switches), float)
    v = zeros((netfabric.num_input_switches, netfabric.num_output_switches, netfabric.num_mid_switches), float)

    for i in range(len(netfabric.paths)):
        v[netfabric.paths[i][0], netfabric.paths[i][1], netfabric.paths[i][2]] = 1      # Set initial value of v

    for k in range(netfabric.num_mid_switches):
        for i in range(netfabric.num_input_switches):
            for j in range(netfabric.num_output_switches):
                vpsum = 0
                vqsum = 0

                for p in range(netfabric.num_input_switches):
                    if p != i:
                        vpsum = vpsum + v[p, j, k]

                for q in range(netfabric.num_output_switches):
                    if q != j:
                        vqsum = vqsum + v[i, q, k]

                if i == input_switch and j in cast_as(output_switches, 'list'):      # Indicates there is a new call request in (i, j)
                    call_requests = netfabric.call_requests(i, j) + 1

                u[i, j, k] = u[i, j, k] + (a - (c/r)*sum(v, (1, 2))[k])*sign(call_requests - sum(v, 2)) - b*(vpsum + vqsum)

                if u[i, j, k] <= -alpha:
                    v[i, j, k] = 0
                    netfabric.set_path(i, k, j, False)
                elif u[i, j, k] >= alpha:
                    v[i, j, k] = 1
                    netfabric.set_path(i, k, j, True)


print sum([[[1, 2], [3, 4]], [[5, 6], [7, 8]]], 0)
print sum([[[1, 2], [3, 4]], [[5, 6], [7, 8]]], 1)
print sum([[[1, 2], [3, 4]], [[5, 6], [7, 8]]], 2)
print sum([[[1, 2], [3, 4]], [[5, 6], [7, 8]]], (0, 1))
print sum([[1, 2], [3, 4]], 0)
print sum([[1, 2], [3, 4]], 1)
print sum([[1, 2], [3, 4]], (0, 1))

a = [[1, 2], [3, 4]]
print a[0][0]
print a[0][1]
print a[1][0]
print a[1][1]

b = [[[1, 2], [3, 4]], [[5, 6], [7, 8]]]

print b[0][0][0]
print b[0][0][1]
print b[0][1][0]
print b[0][1][1]
print b[1][0][0]
print b[1][0][1]
print b[1][1][0]
print b[1][1][1]