# stripeDiff


## Installation
The following codes install stripeDiff in home directory. Users can install it anywhere they want by changing '~' to their desired directory.

    cd ~
    git clone https://github.com/GuangyWang/stripeDiff.git
    
    


Execution
----------
In general, strapCalling.sh can be executed by following command line options:

    strapCalling.sh [options]* <-a matrixA> <-b matrixB>
    Required arguments:
        <-a matrixA> defines the input matrix A
        <-b matrixB> defines the input matrix B
        <-n name> defines two sample names and chromosome. e.g., wt,mutant,chr1
    Optional arguments:
        [-l length] defines row and column number of split submatrix
        [-o outDir] defines the path to output files. It must be set if the two input matrices are not under the same directory



