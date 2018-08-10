## Computes genomic relationship matrices as used in Bellot et al., 2018, Genetics, in press

## dogrm.f90
Computes genomic relationship matrices

To compile

   gfortran dogrm.f90 -lblas -O4 -o dowgrm

To run

   cat genotype_file | ./dogrm -nind Nind [-ploidy ploidy] [-dom] [-hap] [-maf maf]

   genotype_file contains, per row, the genotypes or alleles of all individuals. genotypes are coded 0,1,2 and alleles 0,1, space separated
   Nind: number of individuals
   ploidy: ploidy [2]
   -dom: specifies whether dominance matrix is computed
   -hap specifies whether both alleles (0,1) are read (genotypes 0, 1, 2 read by default) 
   maf: minimum allele frequency required [1e-6]

To compute additive matrix

   cat genotype_file | ./dogrm -nind Nind > add.G

To compute dominance matrix (Vitezica et al 2013)

   cat genotype_file | ./dogrm -nind Nind -dom > dom.G

## dobglr.R
R script to estimate variance components using BGLR (Perez & De Los Campos et al. 2014)
