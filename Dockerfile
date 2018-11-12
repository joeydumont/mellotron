FROM ubuntu:bionic

MAINTAINER joeydumont "https://github.com/joeydumont"

# Install packages for building ruby
RUN apt-get update
RUN apt-get install -y --force-yes build-essential wget git cmake clang ninja-build
RUN apt-get install -y --force-yes libhdf5-dev libboost-all-dev libgsl-dev libarmadillo-dev libgtest-dev
RUN apt-get clean

# Compile gtest as a shared library.
WORKDIR "/usr/src/googletest"
RUN mkdir build
RUN cd build
RUN cmake -G Ninja .. -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_INSTALL_LIBDIR=lib -DBUILD_SHARED_LIBS=ON
RUN cmake --build .
RUN cmake --build . --target install


# Clone the repo.
RUN git clone https://github.com/joeydumont/mellotron joeydumont/mellotron
WORKDIR "joeydumont/mellotron"
RUN git submodule update --init --recursive
