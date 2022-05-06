#! python3
import random
import sys
wps = open(sys.path[0] + "/weight_points", "w")
for i in range(0, 6):
    wps.write(str(random.uniform(0.1, 0.8)) + "\t" +
              str(random.uniform(0.1, 0.8)) + "\t" +
              str(random.uniform(0, 0.3)) + "\n"
              # "0.01\n"
              )

wps.close()

dists = open(sys.path[0] + "/marginals", "w")
for j in range(1, 6):
    supp_size = random.randint(2, 2)
    for i in range(0, supp_size):
        dists.write(str(random.uniform(0.1, 0.8)) + "\t" +
                    str(random.uniform(0.1, 0.8)) + "\t" +
                    str(random.uniform(0, 1)) + "\n"
                    )
    dists.write("\n")
dists.close()
