import numpy as np
import sys

print("Hello world from topo.py", file=sys.stderr)

def topo(isd, ied, jsd, jed, max_depth):
    print("Hello world from topo.topo()", file=sys.stderr)

    D = np.empty((ied-isd+1, jed-jsd+1))
    print(f"flat bottom shape {D.shape} depth {max_depth}", file=sys.stderr)
    D[:] = max_depth

    return D
