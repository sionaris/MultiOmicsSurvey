unset CC CXX FC F77 LD AR AS NM OBJDUMP RANLIB STRIP CFLAGS CPPFLAGS CXXFLAGS LDFLAGS

# 2) Load the same toolchain you use for jobs
module purge
module load rhel8/cclake/base
module load r/4.4.0/gcc/3z5n2tph

# 3) Compile the shared library (Linux .so)
cd Resources/MDICC
rm -f projsplx_R.o projsplx_R.so

# Use R's build flags (best)
R CMD SHLIB projsplx_R.c -o projsplx_R.so

# 4) Verify it is a Linux shared object
file projsplx_R.so
