#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pyopencl as cl
from timeit import default_timer

N = 50000000
N = 500
if __name__ == "__main__":
    ## Step #1. Obtain an OpenCL platform.
    platform = cl.get_platforms()[0]
     
    ## It would be necessary to add some code to check the check the support for
    ## the necessary platform extensions with platform.extensions

    ## Step #2. Obtain a device id for at least one device (accelerator).
    device = platform.get_devices()[1]
    print device

    ## It would be necessary to add some code to check the check the support for
    ## the necessary device extensions with device.extensions

    ## Step #3. Create a context for the selected device.
    ctx  = cl.Context([device])

    # why does the code breaks down when I remove astype
    a_np = np.random.rand(N).astype(np.float32)
    b_np = np.random.rand(N).astype(np.float32)

    
    start = default_timer()
    
    #ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)

    mf = cl.mem_flags
    a_g = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=a_np)
    b_g = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=b_np)

    

    prg = cl.Program(ctx, """
    __kernel void sum(__global const float *a_g, __global const float *b_g, __global float *res_g) {
      int gid = get_global_id(0);
      res_g[gid] = a_g[gid] + b_g[gid];
    }
    """).build()

    res_g = cl.Buffer(ctx, mf.WRITE_ONLY, a_np.nbytes)
    prg.sum(queue, a_np.shape, None, a_g, b_g, res_g)

    res_np = np.empty_like(a_np)
    cl.enqueue_copy(queue, res_np, res_g)


    print default_timer()-start
    print res_np[1:5]
    
    start = default_timer()
    res_npCPU = a_np + b_np
    print default_timer()-start
    print res_npCPU[1:5]
    
    # Check on CPU with Numpy:
    #print(res_np - (a_np + b_np))
    #print(np.linalg.norm(res_np - (a_np + b_np)))
    #print default_timer()-start
