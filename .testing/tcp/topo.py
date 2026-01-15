import numpy as np
import sys

print("Hello world from topo.py", file=sys.stderr)

def topo(isd, ied, jsd, jed, max_depth):
    print("Hello world from topo.topo()", file=sys.stderr)

    D = np.empty((ied-isd+1, jed-jsd+1))
    print(f"flat bottom shape {D.shape} depth {max_depth}", file=sys.stderr)
    D[:] = max_depth

    return D

def vel(udims, vdims, nz):
    print("Hello world from topo.vel()", file=sys.stderr)

    isd, ied, jsd, jed = udims
    u = np.empty((ied-isd+1, jed-jsd+1, nz))
    u[:,:,:] = 0.0

    isd, ied, jsd, jed = vdims
    v = np.empty((ied-isd+1, jed-jsd+1, nz))
    v[:,:,:] = 0.0

    print(f"velocity shapes: {u.shape}, {v.shape}, {u.dtype}, {v.dtype}", file=sys.stderr)

    return u, v
