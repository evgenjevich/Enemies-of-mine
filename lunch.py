
import matplotlib.pyplot as plt
import json
import hackathon_params as hp

hp.hackParams(10,200)

import itertools
import subprocess
N=(10, 40, 100, 200, 400)
steps=(2, 10, 20)
dx = (20, 5, 2, 1, .5) 
tag='benchmark1a'
bashCommand = "smt run params.json --tag=$tag"

for N, steps, dx in itertools.product(N, steps, dx):
    hp.hackParams(N, steps, dx)
    process = subprocess.Popen(bashCommand.split())
    output = process.communicate()[0]


# bashCommand = "pwd"
# process = subprocess.Popen(bashCommand.split())
# output = process.communicate()[0]
