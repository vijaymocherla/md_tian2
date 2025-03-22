cmake_minimum_required(VERSION 3.10)

# Project name and language
project(md_tian2 VERSION 2.1
        DESCRIPTION "md_tian an MD code in Fortran for surface scattering"
        LANGUAGES Fortran)

enable_language(Fortran)
enable_language(CXX)

option(INTEL "compile with ifort and mkl" OFF)
option(OPENMP "compile with openmp support" OFF)
option(DEBUG "compile with debug flags" OFF)

find_package(MPI REQUIRED)
find_package(mace REQUIRED)
set(CMAKE_PREFIX_PATH ${CONDA_PREFIX})

if(INTEL)
    set(CMAKE_Fortran_COMPILER ifx)
else()
    set(CMAKE_Fortran_COMPILER gfortran)
endif()

if(INTEL)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -flto -lstdc++")
    if(DEBUG)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -qmkl -xHost -fp-model strict -check bounds -g -traceback") # -Ofast -ipo
    else()
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -flto -xHost -qmkl -Ofast")
    endif()
    if(OPENMP)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp -parallel -parallel-source-info=1 ")
    endif()
else()
    set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -std=legacy")
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    set(libraries ${libraries} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    if(DEBUG)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcheck=all -Wall -Wconversion -g -Og -fbacktrace -fbounds-check") # -ffpe-trap=zero,invalid") #-flto
    else()
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Ofast -ffast-math -march=native") #-flto
    endif()
    if(OPENMP)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
    endif()
endif()




set(SUBROUTINES
    src/md_tian2.f90
    src/RuNNer_src/abstime.f90
    src/RuNNer_src/addatoms.f90
    src/RuNNer_src/atomsymfunction1.f90
    src/RuNNer_src/atomsymfunction2.f90
    src/RuNNer_src/atomsymfunction3Andi.f90
    src/RuNNer_src/atomsymfunction3.f90
    src/RuNNer_src/atomsymfunction4.f90
    src/RuNNer_src/atomsymfunction5.f90
    src/RuNNer_src/atomsymfunction6.f90
    src/RuNNer_src/atomsymfunction8.f90
    src/RuNNer_src/atomsymfunction9.f90
    src/RuNNer_src/calconenn.f90
    src/RuNNer_src/checkelement.f90
    src/RuNNer_src/check_type.f90
    src/RuNNer_src/convert_frequencies.f90
    src/RuNNer_src/dia_hessian.f90
    src/RuNNer_src/getatomsymfunctions.f90
    src/RuNNer_src/getatomsymfunctions_sfg.f90
    src/RuNNer_src/getcutoff.f90
    src/RuNNer_src/getddnodes_values.f90
    src/RuNNer_src/getdnodes_values.f90
    src/RuNNer_src/getgroupcutoffvalues.f90
    src/RuNNer_src/getgroupneighlist.f90
    src/RuNNer_src/getlistdim.f90
    src/RuNNer_src/getneighboridxatomic_para.f90
    src/RuNNer_src/getneighborsatomic_para.f90
    src/RuNNer_src/getshortatomic.f90
    src/RuNNer_src/getsymfunctionderivatives.f90
    src/RuNNer_src/getsymfunctionsecondderivatives.f90
    src/RuNNer_src/getvibrationalfrequencies.f90
    src/RuNNer_src/getvolume.f90
    src/RuNNer_src/mpifitdistribution.f90
    src/RuNNer_src/neighbor_para_short.f90
    src/RuNNer_src/nuccharge.f90
    src/RuNNer_src/predictionshortatomic.f90
    src/RuNNer_src/sortelements.f90
    src/RuNNer_src/sortsymfunctions.f90
    src/RuNNer_src/translate.f90
    src/RuNNer_src/mpi_dummy_routines.f90
)

set(MODULES
    src/constants.f90
    src/output_mod.f90
    src/pes_nene_mod.f90
    src/fit.f90
    src/md_algo.f90
    src/pes_emt_mod.f90
    src/pes_non_mod.f90
    src/trajectory_info.f90
    src/force.f90
    src/md_init.f90
    src/pes_ho_mod.f90
    src/pes_rebo_mod.f90
    src/universe_mod.f90
    src/geometry_opt.f90
    src/pes_lj_mod.f90
    src/rpmd.f90
    src/useful_things.f90
    src/open_file.f90
    src/pes_mace_mod.f90
    src/run_config.f90
    src/RuNNer_src/fileunits.f90 
    src/RuNNer_src/predictionoptions.f90 
    src/RuNNer_src/nnflags.f90 
    src/RuNNer_src/nnconstants.f90 
    src/RuNNer_src/globaloptions.f90 
    src/RuNNer_src/mode1options.f90 
    src/RuNNer_src/nnshort_atomic.f90 
    src/RuNNer_src/nnewald.f90 
    src/RuNNer_src/saturation.f90 
    src/RuNNer_src/structures.f90 
    src/RuNNer_src/symfunctions.f90 
    src/RuNNer_src/symfunctiongroups.f90 
    src/RuNNer_src/timings.f90 
    src/RuNNer_src/fittingoptions.f90 
    src/RuNNer_src/mpi_dummy.f90

)



add_executable(md_tian2 ${SUBROUTINES} ${MODULES} ${libraries})
target_include_directories(md_tian2 PRIVATE ${mace_INCLUDE_DIRS})
target_link_libraries(md_tian2 ${libraries} PRIVATE mace)
# Add MACE include directories to md_tian2

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
