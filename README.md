## Computes genomic relationship and NRM matrices 
### as used in Bellot et al., 2018, Can Deep Learning Improve Genomic Prediction of Complex Human Traits? GENETICS (https://doi.org/10.1534/genetics.118.301298)

### dogrm.f90
Computes genomic relationship matrices

### doA.f90
Computes numerator  relationship matrices and inverse

### To compile

```
gfortran dogrm.f90 -lblas -O4 -o dogrm
gfortran doA.f90 -O4 -o doa
```

### To run dowgrm

   `cat genotype_file | ./dogrm -nind Nind [-ploidy ploidy] [-dom] [-maf maf]`

- genotype_file contains, per row, all alleles per individual per SNP, 2 alleles for diploids, 4 for tetraploids, etc
   
- Nind: number of individuals
   
- -ploidy: ploidy [2]
   
- -dom: specifies whether dominance matrix is computed
   
- -maf: minimum allele frequency required [1e-6]
   

1. To compute genomic additive matrix

   `cat genotype_file | ./dogrm -nind Nind > add.G`

2. To compute genomic dominance matrix (Vitezica et al 2013)

   `cat genotype_file | ./dogrm -nind Nind -dom > dom.G`
   
### To run doa

1. Direct NRM
   `cat pedfile | ./doa -nind nind > add.nrm`

2. Inverse (Hendersons' rules)
`cat pedfile | ./doa -nind nind -inv > add_inverse.nrm`

where:
* `pedfile` is pedigree file with ther columns `ind id_father id_mother` with integer ids 1, ... nind, 0 for unknown parents.
* `nind` no. of inds
* `add.nrm`and `add_inverse.nrm` are output file

### dobglr.R
R script to estimate variance components using BGLR (Perez & De Los Campos et al. 2014)
