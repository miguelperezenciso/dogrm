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

### To run dogrm

   `cat genotype_file | ./dogrm -nind Nind [-ploidy ploidy] [-dom] [-maf maf] [-maxmiss pmiss]`

- genotype_file contains, per row, all alleles per individual per SNP, 2 alleles for diploids, 4 for tetraploids, etc
   
- Nind: number of individuals
   
- -ploidy: ploidy [2]
   
- -dom: specifies whether dominance matrix is computed
   
- -maf: minimum allele frequency required [1e-6]

- -maxmiss: max percentage of missing values [0] 

**If `maxmiss` is set to a value larger than 0, missing values are replaced by mean frequency.**

**Markers with missing values are not used in computing dominance matrix.**

1. To compute genomic additive matrix

   `cat genotype_file | ./dogrm -nind Nind > add.G`

2. To compute genomic dominance matrix (Vitezica et al 2013)

   cat genotype_file | ./dogrm -nind Nind -dom > dom.G
   
### Converting vcf into genotype files
If you have a vcf file, you can edit it as follows to feed `dogrm`

    grep -v '#' test.vcf | cut -f 10- | sed 's/\// /g'| sed 's/\./9/g' | ./dogrm -nind Nind > add.G

Here we replace missing genotypes (`.`) with ` 9`.

### Using plink ped files
Plink ped files need to be transposed to run `dogrm`. The best way is to use `plink` utility

    plink --file test --recode vcf --out test

produces a `test.vcf` file, which can be parsed as above. 

### Files
Examples to obtain add and dom matrices

```
grep -v '#' files/test.vcf | cut -f 10- | sed 's/\// /g'| sed 's/\./9/g' | ./dogrm -nind 26 -maxmiss 0 > add.G
grep -v '#' files/test.vcf | cut -f 10- | sed 's/\// /g'| sed 's/\./9/g' | ./dogrm -dom -nind 26 -maxmiss 0 > dom.G 
```
   
### To run doA
(I need to add an example pedigree file)

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
