Detecting Causative Mutations from High-throughput Sequencing on Unordered Genomes
========================================================

Background: Mapping by HTS (Next-generation Mapping)
---------

In plants, inbred backcross lines (IBLs) refers to lines (i.e. populations) of plants derived from the repeated backcrossing of a line with artificially recombinant DNA with the wild type, operating some kind of selection that can be phenotypical or through a molecular marker (for the production of introgression lines).

Purpose
------

Modelling
-------

1. Why not use real data straight away?
 - With model data, I know the location of the causative mutation and the correct order of the fragments. This means I can tell whether or not my algorithm works.
 - I can generate datasets that model different experiments/sample data quickly: e.g. different fragment sizes, different SNP distributions, different experiment
 
A beneficial Arabidopsis phenotype is created by experimentally by inducing SNPs with EMS (e.g. a disease resistance mutation). The phenotype altering, reccesive homozygous SNP needs to be identified (cannot be hetetozygous or the selection experiment wouldn't work). The mutant individual is back-crossed with a mapping line. This mapping line is a cross between two Arabidopsis ecotypes (one of which the same ecotype as the mutant), and therefore has heterozygous SNPs across its genome.

Back-crossing like this multiple times, but always selecting for progeny with the mutant phenotype, will result in an individual that has a non-recombinant region with high homozygous to heterozygous SNP ratio near the mutant position, due to linkage. This is because with each back-cross, the EMS induced homozygous SNPs closest to mutation being selected for have a higher chance of being conserved through linkage, than those further from the location of the mutation. 

In other words, the remaining linked homozygous SNPs will be (normally) distributed around the phenotype causing SNP. Plotting the homozygous/heterozygous SNP ratio should reveal the non-recombinant region has a high ratio, because of the high density of homozygous SNPs here, when compared to the consistent distribution of heterozygous SNPs (from the mapping line); the causative SNP should be located at the peak of this ratio distribution.

Re-ordering the Genome
----------

After creating a model genome based on the expected result of the back-crossing experiment detailed above, I am designing an algorithm that will locate the causative mutation, based on the distribution of SNPs. To emulate real data, I first split the genome into fragments, which model contigs assembled from high-throughput sequencing reads. My algorithm will need to rearrange these fragments into their proper order, then locate the causative mutation.

Genetic Algorithm
--------



Permutation Fitness: Q-Q plot
------------

Finding the Peak: The Causative Mutation
--------