## Computes genomic relationship and NRM matrices as used in Bellot et al., 2018, Can Deep Learning Improve Genomic Prediction of Complex Human Traits? GENETICS (https://doi.org/10.1534/genetics.118.301298)

### dogrm.f90
Computes genomic relationship matrices

### doA.f90
Computes numerator  relationship matrices and inverse

### To compile

```
gfortran dogrm.f90 -lblas -O4 -o dowgrm
gfortran doA.f90 -O4 -o doa
```

### To run dowgrm

   `cat genotype_file | ./dogrm -nind Nind [-ploidy ploidy] [-dom] [-hap] [-maf maf]`

- genotype_file contains, per row, the genotypes or alleles of all individuals. genotypes are coded 0,1,2 and alleles 0,1, space separated
   
- Nind: number of individuals
   
- ploidy: ploidy [2]
   
- -dom: specifies whether dominance matrix is computed
   
- -hap specifies whether both alleles (0,1) are read (genotypes 0, 1, 2 read by default) 
   
- maf: minimum allele frequency required [1e-6]
   

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
