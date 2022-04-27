#! python3
import random
import sys
wps = open(sys.path[0] + "/weight_points", "w")
for i in range(0, 83):
    wps.write(str(random.uniform(0.1, 0.8)) + "\t" +
              str(random.uniform(0.1, 0.8)) + "\t" +
              # str(random.uniform(0, 0.3)) + "\n"
              "0\n"
              )

wps.close()

dists = open(sys.path[0] + "/marginals", "w")
for j in range(0, 6):
    supp_size = random.randint(1, 8)
    for i in range(0, supp_size):
        dists.write(str(random.uniform(0.1, 0.8)) + "\t" +
                    str(random.uniform(0.1, 0.8)) + "\t" +
                    str(random.uniform(0, 1)) + "\n"
                    )
    dists.write("\n")

dists.close()
