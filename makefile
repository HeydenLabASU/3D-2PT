COMP=gcc
flags=-lm -O3 -w
FT=-lfftw3f
OMP=-fopenmp
DEPS=lib/fatal.o \
lib/geo.o \
lib/GMXtrrio.o \
lib/io.o \
lib/alloc.o \
lib/top.o \
lib/select.o \
lib/grps.o \
lib/job.o \
lib/align.o \
lib/qsort.o \
lib/mol.o

$(shell   mkdir -p bin)
$(shell   mkdir -p lib)
$(shell   chmod +x test-system/*.sh)
$(shell   chmod +x test-system/3D-2PT-files/*.sh)

all: gmxtop gmxtop-2 psf psf2 psf3 topsf dens3D water3D water3D_noRot water3D_noRot_cryst pot3D pot3D_cryst pot3D_voxAlign pot3D-wsScale pot3D-scH pot3D-OPLS pot3D-OPLS_cryst trans-vdos2entropy rot-vdos2entropy add_cubes averCubefiles weightedAverCubeFiles average process-cube

gmxtop: ${DEPS} src/core/gmxtop.c
	${COMP} ${DEPS} src/core/gmxtop.c -o bin/gmxtop.exe ${flags}

gmxtop-2: ${DEPS} src/core/gmxtop-2.c
	${COMP} ${DEPS} src/core/gmxtop-2.c -o bin/gmxtop-2.exe ${flags}

psf: ${DEPS} src/core/psf.c
	${COMP} ${DEPS} src/core/psf.c -o bin/psf.exe ${flags}

psf2: ${DEPS} src/core/psf2.c
	${COMP} ${DEPS} src/core/psf2.c -o bin/psf2.exe ${flags}

psf3: ${DEPS} src/core/psf3.c
	${COMP} ${DEPS} src/core/psf3.c -o bin/psf3.exe ${flags}

topsf: ${DEPS} src/core/topsf.c
	${COMP} ${DEPS} src/core/topsf.c -o bin/topsf.exe ${flags}

dens3D: ${DEPS} src/dens3D/dens3D.c
	${COMP} ${DEPS} src/dens3D/dens3D.c -o bin/dens3D.exe ${flags}

water3D: ${DEPS} src/water3D/water3D.c
	${COMP} ${DEPS} src/water3D/water3D.c -o bin/water3D.exe ${flags} ${OMP} ${FT}

water3D_noRot: ${DEPS} src/water3D/water3D_noRot.c
	${COMP} ${DEPS} src/water3D/water3D_noRot.c -o bin/water3D_noRot.exe ${flags} ${OMP} ${FT}

water3D_noRot_cryst: ${DEPS} src/water3D/water3D_noRot_cryst.c
	${COMP} ${DEPS} src/water3D/water3D_noRot_cryst.c -o bin/water3D_noRot_cryst.exe ${flags} ${OMP} ${FT}

pot3D: ${DEPS} src/pot3D/pot3D.c
	${COMP} ${DEPS} src/pot3D/pot3D.c -o bin/pot3D.exe ${flags} ${OMP}

pot3D_cryst: ${DEPS} src/pot3D/pot3D_cryst.c
	${COMP} ${DEPS} src/pot3D/pot3D_cryst.c -o bin/pot3D_cryst.exe ${flags} ${OMP}

pot3D_voxAlign: ${DEPS} src/pot3D/pot3D_voxAlign.c
	${COMP} ${DEPS} src/pot3D/pot3D_voxAlign.c -o bin/pot3D_voxAlign.exe ${flags} ${OMP}

pot3D-wsScale: ${DEPS} src/pot3D/pot3D-wsScale.c
	${COMP} ${DEPS} src/pot3D/pot3D-wsScale.c -o bin/pot3D-wsScale.exe ${flags} ${OMP}

pot3D-scH: ${DEPS} src/pot3D/pot3D-scH.c
	${COMP} ${DEPS} src/pot3D/pot3D-scH.c -o bin/pot3D-scH.exe ${flags} ${OMP}

pot3D-OPLS: ${DEPS} src/pot3D/pot3D-OPLS.c
	${COMP} ${DEPS} src/pot3D/pot3D-OPLS.c -o bin/pot3D-OPLS.exe ${flags} ${OMP}

pot3D-OPLS_cryst: ${DEPS} src/pot3D/pot3D-OPLS_cryst.c
	${COMP} ${DEPS} src/pot3D/pot3D-OPLS_cryst.c -o bin/pot3D-OPLS_cryst.exe ${flags} ${OMP}

trans-vdos2entropy: src/post-process/trans-vdos2entropy_lowMem.c
	${COMP} src/post-process/trans-vdos2entropy_lowMem.c -o bin/trans-vdos2entropy_lowMem.exe ${flags}

rot-vdos2entropy: src/post-process/rot-vdos2entropy_lowMem.c
	${COMP} src/post-process/rot-vdos2entropy_lowMem.c -o bin/rot-vdos2entropy_lowMem.exe ${flags}

add_cubes: src/post-process/add_cubes.c
	${COMP} src/post-process/add_cubes.c -o bin/add_cubes.exe ${flags}

averCubefiles: src/post-process/averCubefiles.c
	${COMP} src/post-process/averCubefiles.c -o bin/averCubefiles.exe ${flags}

weightedAverCubeFiles: src/post-process/weightedAverCubeFiles.c
	${COMP} src/post-process/weightedAverCubeFiles.c -o bin/weightedAverCubeFiles.exe ${flags}

average: ${DEPS} src/post-process/average.c
	${COMP} ${DEPS} src/post-process/average.c -o bin/average.exe ${flags}

process-cube: src/post-process/process-cube.c
	${COMP} src/post-process/process-cube.c -o bin/process-cube.exe ${flags}

#dependencies
lib/fatal.o: src/core/fatal.c
	${COMP} -c src/core/fatal.c -o lib/fatal.o

lib/geo.o: src/core/geo.c
	${COMP} -c src/core/geo.c -o lib/geo.o

lib/GMXtrrio.o: src/core/GMXtrrio.c
	${COMP} -c src/core/GMXtrrio.c -o lib/GMXtrrio.o

lib/io.o: src/core/io.c
	${COMP} -c src/core/io.c -o lib/io.o

lib/alloc.o: src/core/alloc.c
	${COMP} -c src/core/alloc.c -o lib/alloc.o

lib/top.o: src/core/top.c
	${COMP} -c src/core/top.c -o lib/top.o

lib/select.o: src/core/select.c
	${COMP} -c src/core/select.c -o lib/select.o

lib/grps.o: src/core/grps.c
	${COMP} -c src/core/grps.c -o lib/grps.o

lib/job.o: src/core/job.c
	${COMP} -c src/core/job.c -o lib/job.o

lib/align.o: src/core/align.c
	${COMP} -c src/core/align.c -o lib/align.o

lib/qsort.o: src/core/qsort.c
	${COMP} -c src/core/qsort.c -o lib/qsort.o

lib/mol.o: src/core/mol.c
	${COMP} -c src/core/mol.c -o lib/mol.o
