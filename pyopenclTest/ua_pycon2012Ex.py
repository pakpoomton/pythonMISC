# 1st example from UA PYCON 2012, Tomasz Rybak PyOpenCL
import pyopencl
import numpy

# platform = pyopencl.get_platforms()[0]
# device = platform.get_devices(device_type=GPU)[0]
# context = pyopencl.Context(devices=[device])
# queue = pyopencl.CommandQueue(context, device=device)

#platform = pyopencl.get_platforms()[0]
#device = platform.get_devices(device_type=GPU)[0]
context = pyopencl.Context()
queue = pyopencl.CommandQueue(context)

a = numpy.zeros(1000).astype(numpy.float32)
a_gpu = pyopencl.array.Array(context, shape=a.shape)
program = pyopencl.Program(context, """
    __kernel void increase(__global float *a)
    {
        int gid = get_global_id(0);
        a[gid] = a[gid]+1.0;
    }
""").build()

pyopencl.enqueue_copy(queue, a, a_gpu)
program.increase(queue, a.shape, None, a_gpu)
a = a_gpu.get()
