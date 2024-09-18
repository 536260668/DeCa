# DeCa
DeCa is a high-performance, attention-based tool designed for the rapid and accurate detection of somatic variants. 

### Prerequisite
* [htslib](https://github.com/samtools/htslib) : The library to read hts files.
* [boost](https://www.boost.org/) : Boost provides free peer-reviewed portable C++ source libraries.
* [libtorch](https://pytorch.org/) : An open source machine learning framework that accelerates the path from research prototyping to production deployment.

### Building DeCa

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

### Running DeCa

#### Parameters

##### Required Arguments:

- `--input`, `-I:String`  
  BAM/SAM/CRAM file containing reads. This argument must be specified at least once.  
  **Required**

- `--output`, `-O:File`  
  The output recalibration table file to create.  
  **Required**

- `--reference`, `-R:String`  
  Reference sequence file.  
  **Required**

- `--normal:String`  
  BAM sample name of normal(s).  
  **Required**

##### Optional Arguments:

- `--callable-depth:Int`  
  Minimum depth to be considered callable for Mutect stats. Does not affect genotyping.  
  **Default:** 10

- `--max-prob-propagation-distance:Int`  
  Upper limit on how many bases away probability mass can be moved when calculating the boundaries between active and inactive assembly regions.

- `--active-probability-threshold:Float`  
  Minimum probability for a locus to be considered active.  
  **Default:** 0.002

- `--assembly-region-padding:Int`  
  Number of additional bases of context to include around each assembly region.  
  **Default:** 100

- `--max-assembly-region-size:Int`  
  Maximum size of an assembly region.  
  **Default:** 300

- `--min-assembly-region-size:Int`  
  Minimum size of an assembly region.  
  **Default:** 50

- `--bqsr-with-mutect`  
  Add `ApplyBQSR` into this tool. If set, `--tumor-table` and `--normal-table` are required.

- `--tumor-table:File`  
  Recalibration table for tumor reads, generated in the BaseRecalibrator part.

- `--normal-table:File`  
  Recalibration table for normal reads, generated in the BaseRecalibrator part.

- `-T:Int`  
  Size of thread pool.

- `-L:String`  
  Specifies the name of the chromosome to be processed.

- `-M:String`  
  Path to the machine learning model.

- `--max-reads-per-alignment-start-M:Int`  
  Maximum number of reads to retain per alignment start position. Reads above this threshold will be downsampled. Set to 0 to disable.

#### example
```sh
./Deca -R reference.fa \
-I tumor.bam \
-I normal.bam \
-O output.vcf \
-T 64 \
--normal normal \
-M trans.pt \
-L 10
```
