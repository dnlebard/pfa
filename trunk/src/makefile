G = #-g
COMPC = gcc 
CUDAC = nvcc 
OPTC =  -O3 
CUDA-OPT = ${DBG} -O3 -arch sm_13 #-deviceemu  #--ptxas-options=-v 
OBJCPU  =  gzReader.o dcdReader.o psfReader.o common.o checksz.o nxtsec.o readprm.o defines.o fastEner.o
OBJCUDA  =  dvos_gpu.o shell_calc_gpu.o pca_gpu.o ${OBJCPU}
OBJ  = ${OBJCPU} 
PGOPT = ${DBG} -O3 ${MACHINE_FLAGS} 
OPT1  =  ${PGOPT} #${DBG} 
CUDAEMURT = #-lcudartemu
CUDART = -lcudart 
CUDARTLIBLINK = ${CUDAEMURT} ${CUDART}
CUDALINK = ${DBG} ${CUDARTLIBLINK} -L${CUDALIB} -I${CUDAINC} ${CUDA_LINK_OPT} 
OPTLINK = ${DBG} 
pfa.cpu.e: $(OBJ)
	$(COMPF) ${BLASLINK} $(OPT1) ${OPTLINK} -o pfa.cpu.e $(OBJ)
pfa.cuda.e: $(OBJCUDA)
	$(COMPF) ${BLASLINK} $(OPT1) ${CUDALINK} -o pfa.cuda.e $(OBJCUDA)
pca_gpu.o: pca_gpu.cu
	$(CUDAC) $(CUDA-OPT) -o pca_gpu.o -c pca_gpu.cu
shell_calc_gpu.o: shell_calc_gpu.cu
	$(CUDAC) $(CUDA-OPT) -o shell_calc_gpu.o -c shell_calc_gpu.cu
dvos_gpu.o: dvos_gpu.cu
	$(CUDAC) $(CUDA-OPT) -o dvos_gpu.o -c dvos_gpu.cu
gzReader.o: gzReader.c
	$(COMPC) $(OPTC) -o gzReader.o -c gzReader.c
dcdReader.o: dcdReader.c
	$(COMPC) $(OPTC) -o dcdReader.o -c dcdReader.c
psfReader.o: psfReader.c
	$(COMPC) $(OPTC) -o psfReader.o -c psfReader.c
common.o: common.f90
	$(COMP) $(OPT1) -o common.o -c common.f90
defines.o: fastDefines.F 
	$(COMP) $(OPT1) $(CUDA_FLAGS) $(FORT_COMP) -o defines.o -c fastDefines.F 
fastEner.o: fastEner.f 
	$(COMP) $(OPT1) -o fastEner.o -c fastEner.f 
checksz.o: checksz.f
	$(COMP) $(OPT1) -o checksz.o -c checksz.f
readprm.o: readprm.f
	$(COMP) $(OPT1) -o readprm.o -c readprm.f
nxtsec.o: nxtsec.f
	$(COMP) $(OPT1) -o nxtsec.o -c nxtsec.f
clean:
	rm -rf *.o *.e *.mod 
cmm:
	cp ../src/*.F .
	cp ../src/*.f .
	cp ../src/*.f90 .
	cp ../src/*.c .
	cp ../src/*.cu .
	cp ../src/*.h .
	make pfa.cuda.e COMP=" gfortran -I/home/dnlebard/open-mpi/include -I/home/dnlebard/open-mpi/lib -L/home/dnlebard/open-mpi/lib -lmpi_f90 -lmpi_f77 -lmpi -lopen-rte -lopen-pal -ldl -Wl,--export-dynamic -lnsl -lutil -lm -ldl" COMPF="gfortran -I/home/dnlebard/open-mpi/include -I/home/dnlebard/open-mpi/lib -L/home/dnlebard/open-mpi/lib -lmpi_f90 -lmpi_f77 -lmpi -lopen-rte -lopen-pal -ldl -Wl,--export-dynamic -lnsl -lutil -lm -ldl" CUDA_FLAGS="-DENABLE_CUDA" FORT_COMP="-DENABLE_GFORT" CUDALIB="/usr/local/cuda/lib64" CUDAINC="/usr/local/cuda/include" MACHINE_FLAGS="-L/usr/local/cuda/lib64" BLASLINK="-llapack -lblas" 
spur:
	cp ../src/*.F .
	cp ../src/*.f .
	cp ../src/*.f90 .
	cp ../src/*.c .
	cp ../src/*.cu .
	cp ../src/*.h .
	make pfa.cuda.e COMP=mpif90 COMPF=mpif90 CUDA_FLAGS="-DENABLE_CUDA" FORT_COMP="-DENABLE_PGI" MACHINE_FLAGS="-fpic -fastsse -fast -tp barcelona-64 -Minfo" CUDALIB="/usr/local/cuda/lib" CUDAINC="/usr/local/cuda/include" BLASLINK="-Wl,-rpath,${TACC_ACML_LIB} -L${TACC_ACML_LIB} -lacml"
longhorn:
	cp ../src/*.F .
	cp ../src/*.f .
	cp ../src/*.f90 .
	cp ../src/*.c .
	cp ../src/*.cu .
	cp ../src/*.h .
	make pfa.cuda.e COMP=mpif90 COMPF=mpif90 CUDA_FLAGS="-DENABLE_CUDA" FORT_COMP="-DENABLE_GFORT" MACHINE_FLAGS="" CUDALIB=${TACC_CUDA_LIB} CUDAINC=${TACC_CUDA_INC} BLASLINK="-Wl,-rpath,${TACC_MKL_LIB} -L${TACC_MKL_LIB} -lmkl_lapack -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lguide -lpthread"
stampede:
	cp ../src/*.F .
	cp ../src/*.f .
	cp ../src/*.f90 .
	cp ../src/*.c .
	cp ../src/*.cu .
	cp ../src/*.h .
	make pfa.cuda.e COMP=mpif90 COMPF=mpif90 CUDA_FLAGS="-DENABLE_CUDA" FORT_COMP="-DENABLE_GFORT" MACHINE_FLAGS="" CUDALIB=${TACC_CUDA_LIB} CUDAINC=${TACC_CUDA_INC} BLASLINK="-Wl,-rpath,${TACC_MKL_LIB} -L${TACC_MKL_LIB} -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread"
maverick:
	cp ../src/*.F .
	cp ../src/*.f .
	cp ../src/*.f90 .
	cp ../src/*.c .
	cp ../src/*.cu .
	cp ../src/*.h .
	make pfa.cuda.e COMP=mpif90 COMPF=mpif90 CUDA_FLAGS="-DENABLE_CUDA" FORT_COMP="-DENABLE_GFORT" MACHINE_FLAGS="" CUDALIB=${TACC_CUDA_LIB} CUDAINC=${TACC_CUDA_INC} BLASLINK="-Wl,-rpath,${TACC_MKL_LIB} -L${TACC_MKL_LIB} -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread"
ncsa:
	cp ../src/*.F .
	cp ../src/*.f .
	cp ../src/*.f90 .
	cp ../src/*.c .
	cp ../src/*.cu .
	cp ../src/*.h .
	make  pfa.cuda.e COMP=mpif90 COMPF=mpif90 CUDA_FLAGS="-DENABLE_CUDA" FORT_COMP="-DENABLE_PGI" MACHINE_FLAGS="" CUDALIB="/usr/local/cuda-3.0/lib64" CUDAINC="/usr/local/cuda-3.0/include" BLASLINK="-Wl,-rpath,${MKL_HOME}/lib/em64t -L${MKL_HOME}/lib/em64t -lmkl_solver_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread"
delta:
	cp ../src/*.F .
	cp ../src/*.f .
	cp ../src/*.f90 .
	cp ../src/*.c .
	cp ../src/*.cu .
	cp ../src/*.h .
	make pfa.cuda.e COMP=mpif90 COMPF=mpif90 CUDA_FLAGS="-DENABLE_CUDA" FORT_COMP="-DENABLE_GFORT" MACHINE_FLAGS="" CUDA_LINK_OPT="-L/cmm/pkg/cuda/lib -I/cmm/pkg/cuda/include" CUDALIBi="/cmm/pkg/cuda/lib" CUDAINC="/cmm/pkg/cuda/include"
ranger:
	cp ../src/*.F .
	cp ../src/*.f .
	cp ../src/*.f90 .
	cp ../src/*.c .
	cp ../src/*.cu .
	cp ../src/*.h .
	 make pfa.cpu.e COMP=mpif90 COMPF=mpif90 FORT_COMP="-DENABLE_PGI" MACHINE_FLAGS="-fpic -fastsse -fast -tp barcelona-64 -Minfo" BLASLINK="-Wl,-rpath,${TACC_ACML_LIB} -L${TACC_ACML_LIB} -lacml"
kraken:
	cp ../src/*.F .
	cp ../src/*.f .
	cp ../src/*.f90 .
	cp ../src/*.c .
	cp ../src/*.cu .
	cp ../src/*.h .
	 make pfa.cpu.e COMP=ftn COMPF=ftn FORT_COMP="-DENABLE_PGI" MACHINE_FLAGS="-fpic -O3" BLASLINK="-lacml"
