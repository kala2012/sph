MAINTAINERCLEANFILES = Makefile.in
lib_LTLIBRARIES = libSPH.la 

#if HAVE_CUDA
#.cu.lo:
#	$(LIBTOOL) --tag=CC --mode=compile  nvcc  -v -O2 -g  -arch sm_21 -lineinfo --use_fast_math --ptxas-options=-v -ccbin $(CXX) --compiler-options `echo -fPIC -DPIC $(libMC_la_CXXFLAGS)  -DMPICH_SKIP_MPICXX -DMPICH_IGNORE_CXX_SEEK | sed -e "s/ /,/g"` -c $^ -static  
#endif

include_HEADERS = derivatives.h \
	inoutflow.h \
	inputbp.h \
	input.h \
	output.h \
	parameterfile.h \
	param.h \
	System.h	


SPH_SOURCES =  derivatives.c \
	inoutflow.c \
	inputbp.c \
	input.c \
	output.c \
	parameterfile.c \
	System.c

#	${top_srcdir}/deplib/CPU/brute.c \
#	${top_srcdir}/deplib/CPU/deprecated.c \
#	${top_srcdir}/deplib/CPU/dists.c \
#	${top_srcdir}/deplib/CPU/exactDriver.c \
#	${top_srcdir}/deplib/CPU/oneShotDriver.c \
#	${top_srcdir}/deplib/CPU/rbc.c \
#	${top_srcdir}/deplib/CPU/utils.c	

#OTHERS_SOURCES = Hamiltonian.cu \
#		 Math_cuda.cu \
#		 MonteCarlo_cuda.cu \
#	         fft.cu \
#		 Coloring.cu

libSPH_la_SOURCES = $(SPH_SOURCES)

libSPH_la_CFLAGS = -I$(top_srcdir) -I${top_srcdir}/deplib/CPU -I${top_srcdir}/deplib/GPU \
@GSL_CFLAGS@ \
@BLAS_CFLAGS@ \
@IPP_CFLAGS@ 

libSPH_la_LIBADD = @GSL_LIBS@ \
@BLAS_LIBS@ \
@IPP_LIBS@ \
-lm 

libSPH_la_LDFLAGS = -no-undefined -version-info @version_info@ 
