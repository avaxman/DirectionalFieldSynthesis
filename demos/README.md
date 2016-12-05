[![Build Status](https://travis-ci.org/danielepanozzo/tutorial_nrosy.svg?branch=master)](https://travis-ci.org/danielepanozzo/tutorial_nrosy)
[![Build status](https://ci.appveyor.com/api/projects/status/es190vhn9ehbgvad?svg=true)](https://ci.appveyor.com/project/danielepanozzo/tutorial_nrosy)
# NRosy Demo

## Compile

Compile this project using the standard cmake routine:

    git clone --recursive https://github.com/danielepanozzo/tutorial_nrosy
    mkdir build
    cd build
    cmake ..
    make
    ./tutorial_nrosy

This should find and build the dependencies and create a `tutorial_nrosy` binary.
You must clone the repository with the --recursive option.

## Dependencies

This demo is self-contained, just clone the repo recursively:

    git clone --recursive https://github.com/danielepanozzo/tutorial_nrosy
