#! python3
import random
import sys

# Random data to test power diagram algorithm
wps = open(sys.path[0] + "/weight_points", "w")
for i in range(0, 6):
    wps.write(
        str(random.uniform(0.1, 0.8)) + "\t" + str(random.uniform(0.1, 0.8)) +
        "\t" + str(random.uniform(0, 0.3)) + "\n")

wps.close()

# Random data for discrete distributions
dists = open(sys.path[0] + "/marginals", "w")
for j in range(1, 3):
    # j discrete marginals
    supp_size = random.randint(2, 3)
    for i in range(0, supp_size):
        # Data structure (x, y, weight)
        dists.write(
            str(random.uniform(0.1, 0.8)) + "\t" +
            str(random.uniform(0.1, 0.8)) + "\t" + str(random.uniform(0, 1)) +
            "\n")
    dists.write("\n")
dists.close()
