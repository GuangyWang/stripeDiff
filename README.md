# stripeDiff

If you use stripeDiff in your research, please cite: 

## Download
The following codes download stripeDiff to home directory. Users can install it anywhere they want by changing '~' to their desired directory.

    cd ~
    git clone https://github.com/GuangyWang/stripeDiff.git
    
To test stripeDiff, change to the "test" directory and run test.sh
    cd stripeDiff/test
    sh test.sh
    
    


Execution
----------
In general, stripeCalling.sh can be executed by following command line options:

    sh /path_to_stripeCalling/stripeCalling.sh [options]* <-a matrixA> <-b matrixB> <-n name>
    Required arguments:
        <-a matrixA> defines the input matrix A
        <-b matrixB> defines the input matrix B
        
        <-n name> defines two sample names and chromosome. e.g., wt,mutant,chr1
    Optional arguments:
        [-l length] defines row and column number of split submatrix
        [-o outDir] defines the path to output files. It must be set if the two input matrices are not under the same directory



