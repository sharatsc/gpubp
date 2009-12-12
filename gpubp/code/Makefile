CUDA   = /usr/local/cuda
CUDASDK= /usr/local/cudasdk

CFLAGS = --maxrregcount 32 --host-compilation c++ -I $(CUDA)/include -I $(CUDASDK)/common/inc \
            -L $(CUDA)/lib
GPU_SOURCES = bp.cpp bp.cu gpu_main.cpp 
CPU_SOURCES = bp.cpp cpu_main.cpp
LIBS    = 
bp_gpu:
	nvcc $(CFLAGS) $(GPU_SOURCES) $(LIBS)
bp_cpu:
	nvcc -o bp_cpu $(CFLAGS) $(CPU_SOURCES) $(LIBS)
clean:
	rm *.o bp_gpu bp_cpu
