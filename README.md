# stripeDiff
stripeDiff is a software to identify differential stripes between two Hi-C samples. 

Copyright (c) 2020 Dr. Kaifu Chen lab

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

## Input
The main input are two Hi-C matrices. Each matrix could be in sparse (N * N) or verbose format. The verbose format should have three column. The first two columns are two genomic regions, and the third column is the interaction frequency. Here is an example of verbose format (the table header is not needed):
|region1|	region2|	contact frequency|
| ------| ------|----|
|3000000|	3000000|	2095.53|
|3000000	|3010000	|279.50217|
|3010000	|3010000	|2625.0593|
|3000000	|3020000	|139.41035|

## Output
stripeDiff has two output tables: stripes present in matrix A and stripes present in matrix B

| column | explaination |
| ------| ------|
| 1st | chromosome |
| 2nd | left boundary of a stripe |
| 3rd | right boundary of a stripe |
| 4th | expanded left boundary of a stripe |
| 5th | expanded right boundary of a stripe |
| 6th | left signal of a stripe in Matrix A |
| 7th | right signal of a stripe in Matrix A |
| 8th | signal fold change of a stripe in Matrix A |
| 9th | P value of a stripe in Matrix A |
| 10th | left signal of a stripe in Matrix B |
| 11th | right signal of a stripe in Matrix B |
| 12th | signal fold change of a stripe in Matrix B |
| 13th | P value of a stripe in Matrix B |
| 14th | P value for a differential stripe |
| 15th | direction of a stripe, could be left or right |


## Contacts
If you find any bugs or have difficulties in using stripeDiff, please feel free to contact Guangyu Wang (gwang2@houstonmethodist.org).
