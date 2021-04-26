FROM fedora:33

RUN dnf -y update \
    && dnf -y install \
        gcc-gfortran \
        gcc-c++ \
        gcc \
        netcdf-fortran-devel \
        gsl-devel \
        metis-devel \
        lapack-devel \
        openblas-devel \
        cmake \
        make \
        git \
        gnuplot \
    && dnf clean all

# Build the SuiteSparse libraries for sparse matrix support
RUN curl -LO http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-5.1.0.tar.gz \
    && tar -zxvf SuiteSparse-5.1.0.tar.gz \
    && export CXX=/usr/bin/cc \
    && cd SuiteSparse \
    && make install INSTALL=/usr/local BLAS="-L/lib64 -lopenblas"

# install json-fortran
RUN curl -LO https://github.com/jacobwilliams/json-fortran/archive/8.2.0.tar.gz \
    && tar -zxvf 8.2.0.tar.gz \
    && cd json-fortran-8.2.0 \
    && export FC=gfortran \
    && mkdir build \
    && cd build \
    && cmake -D SKIP_DOC_GEN:BOOL=TRUE .. \
    && make install

# copy the analytical code
COPY . /trop-mam4-analytical/

# Install a modified version of CVODE
RUN tar -zxvf /trop-mam4-analytical/test/partmc/cvode-3.4-alpha.tar.gz \
    && cd cvode-3.4-alpha \
    && mkdir build \
    && cd build \
    && cmake -D CMAKE_BUILD_TYPE=release \
             -D CMAKE_C_FLAGS_DEBUG="-g -pg" \
             -D CMAKE_EXE_LINKER_FLAGS_DEBUG="-pg" \
             -D CMAKE_MODULE_LINKER_FLAGS_DEBUG="-pg" \
             -D CMAKE_SHARED_LINKER_FLAGS_DEBUG="-pg" \
             -D KLU_ENABLE:BOOL=TRUE \
             -D KLU_LIBRARY_DIR=/usr/local/lib \
             -D KLU_INCLUDE_DIR=/usr/local/include \
             .. \
    && make install

# Build PartMC
RUN mkdir pmc_build \
    && cd pmc_build \
    && export JSON_FORTRAN_HOME="/usr/local/jsonfortran-gnu-8.2.0" \
    && cmake -D CMAKE_BUILD_TYPE=release \
             -D CMAKE_C_FLAGS_DEBUG="-g -pg" \
             -D CMAKE_Fortran_FLAGS_DEBUG="-g -pg" \
             -D CMAKE_MODULE_LINKER_FLAGS="-pg" \
             -D ENABLE_SUNDIALS:BOOL=TRUE \
             -D ENABLE_GSL:BOOL=TRUE \
             -D SUNDIALS_CVODE_LIB=/usr/local/lib/libsundials_cvode.so \
             -D SUNDIALS_INCLUDE_DIR=/usr/local/include \
             /trop-mam4-analytical/test/partmc \
    && make

# build the code
RUN mkdir /build \
      && cd /build \
      && export JSON_FORTRAN_HOME="/usr/local/jsonfortran-gnu-8.2.0" \
      && cmake -D PARTMC_INCLUDE_DIR="/pmc_build" \
               -D PARTMC_LIB="/pmc_build/libpartmc.a" \
               /trop-mam4-analytical \
      && make
