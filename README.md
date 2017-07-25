[![Build Status](https://travis-ci.org/avaxman/DirectionalFieldSynthesis.svg?branch=master)](https://travis-ci.org/avaxman/DirectionalFieldSynthesis)
[![Build status](https://ci.appveyor.com/api/projects/status/3h035i3nsv86lifl?svg=true)](https://ci.appveyor.com/project/danielepanozzo/directionalfieldsynthesis)

# Directional Field Synthesis, Design, and Processing

This repository contains the material for the [SIGGRAPH 2017 course](http://s2017.siggraph.org/courses/events/directional-field-synthesis-design-and-processing). The repository contains course notes, slides (marked according to order), and demos. The course is given by [Amir Vaxman](http://www.staff.science.uu.nl/~vaxma001/), [Marcel Campen](https://www.graphics.rwth-aachen.de/person/7/), and [Daniele Panozzo](http://cs.nyu.edu/~panozzo/), and was created in collaboration with [Olga Diamanti](http://web.stanford.edu/~diamanti/), [Klaus Hildebrandt](http://graphics.tudelft.nl/~klaus/), [David Bommes](https://www.aices.rwth-aachen.de/en/people/bommes), and [Mirela Ben-Chen](http://mirela.net.technion.ac.il/).

A previous version of this course appeared as a [State-of-the-art report](https://diglib.eg.org/handle/10.1111/cgf12864) in Eurographics 2016, and as a [SIGGRAPH Asia 2016 course](https://sa2016.siggraph.org/en/).

The course notes and slides can be found in the respective folders.

##Installation

To get the library, use:

```bash
git clone --recursive https://github.com/avaxman/DirectionalFieldSynthesis.git
```

The demo are given from the [libdirectional](https://github.com/avaxman/libdirectional) library that is automatically submoduled to this course (be sure to clone the library recursively as explained above). To compile the examples, go into the respective library (e.g., `libdirectional/demos/sampling`) and enter:

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make
```

##Course Schedule

The course will follow the following schedule, spanning 3:15 hours with a break in the middle:

1. Introduction - Marcel Campen.
2. Applications - Daniele Panozzo (including instant meshes demo).
3. Taxonomy - Amir Vaxman.
4. Discretization - Amir Vaxman (including sampling demo).
5. Representation - Marcel Campen.

**Break**

6. Objectives and Constraints - Daniele Panozzo.
7. Visualization - Amir Vaxman (including visualization demo).
8. Live coding demo - Daniele Panozzo.
9. Open Problems - Amir Vaxman.

**Q&A**

