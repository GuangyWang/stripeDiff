# stripeDiff
stripeDiff is a software to identify differential stripes between two Hi-C samples. 

Copyright (c) 2020 Dr. Kaifu Chen lab at the Houston Methodist Hospital

Current version v1.0


## How does stripeDiff work


## Installation
The following codes download stripeDiff to the home directory. Users can install it anywhere they want by changing '~' to their desired directory.

    cd ~
    git clone https://github.com/GuangyWang/stripeDiff.git
    
To test stripeDiff, change to the "test" directory and run test.sh

    cd stripeDiff/test
    sh test.sh


## Software dependancies
Running stripeDiff requires python 3.5.1 or above, R 3.5.1 or above, and the following packages:
- [numpy](https://numpy.org)
- [bcp](https://cran.r-project.org/web/packages/bcp/index.html)
- [plotly](https://cran.r-project.org/web/packages/plotly/index.html)
- [proxy](https://cran.r-project.org/web/packages/proxy/index.html)
- [zoo](https://cran.r-project.org/web/packages/zoo/index.html)


## Execution
----------
In general, stripeCalling.sh can be executed by following command line options:

    sh /path/to/stripeCalling.sh [options]* <-a matrixA> <-b matrixB> <-n name>
    Required arguments:
        <-a matrixA> defines the input matrix A
        <-b matrixB> defines the input matrix B
        <-n name> defines two sample names and chromosome. e.g., wt,mutant,chr1
    
    Optional arguments:
        [-l length] defines row and column number of split submatrix
        [-o outDir] defines the path to output files. It must be set if the two input matrices are not under the same directory
        [-r resolution] defines the bin size. If it is not provided, the bin number of identified stripes will be output


## Contacts
If you find any bugs or have difficulties in using stripeDiff, please feel free to contact Guangyu Wang (gwang2@houstonmethodist.org).
