# DeCa
DeCa is a high-performance, attention-based tool designed for the rapid and accurate detection of somatic variants. 

### Prerequisite
* [htslib](https://github.com/samtools/htslib) : The library to read hts files.
* [boost](https://www.boost.org/) : Boost provides free peer-reviewed portable C++ source libraries.
* [libtorch](https://pytorch.org/) : An open source machine learning framework that accelerates the path from research prototyping to production deployment.

### Building Mutect2Cpp-master

build the code from a Git repository requires following steps:

```sh
git clone https://github.com/536260668/DeCa
cd Deca
mkdir release
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_MAKE_PROGRAM=ninja -G Ninja -S /your/path/to/DeCa -B /your/path/to/DeCa/release
cmake --build /your/path/to/DeCa/release -j <threads>
```

By default, this will build against an HTSlib source tree in `../htslib`.
this program need libtorch in `/usr/libtorch` to run.
