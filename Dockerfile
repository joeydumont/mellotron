FROM ubuntu:bionic

MAINTAINER joeydumont "https://github.com/joeydumont"

# Install packages for building ruby
RUN apt-get update
RUN apt-get install -y --force-yes build-essential wget git cmake clang
RUN apt-get install -y --force-yes libhdf5-dev libboost-dev libgsl-dev libarmadillo-dev
RUN apt-get clean
