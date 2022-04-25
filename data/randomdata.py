#! python3
import random
import sys
file = open(sys.path[0] + "/weight_points", "w")
for i in range(0, 3):
    file.write(str(random.uniform(0, 0.8)) + "\t" +
               str(random.uniform(0, 0.8)) + "\t" +
               # str(random.uniform(0, 0.3)) + "\n"
               "0\n"
               )

file.close()
