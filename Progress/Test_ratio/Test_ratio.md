25/04/14
========================================================

When comparing distributions in a Q-Q plot, the values of the distribution do not matter, just the distribution shape. My idea is to use typical peak shapes found in experimental data to find permutations (contig order) with a similar distribution shape.

The problem is that the shape of the distribution in experimental data is dependent on the location of the mutation. A distribution over the same range (genome size) with a peak located within that range at the equivalent point (e.g. exactly the centre) will give a good q-q plot correlation. e.g. any normal distribution created by rnorm with same standard deviation will have a good q-q plot with eachother (test this) regardless of the value centred around. 

The problem is: there may be many permuations that find the peak in different places.

So I propose an alternate use for my algorithm:

Lets say we know that a peak is in a certain place on one EMS genome, and we expect it to be there in another. We can re-order contigs from the second genome based on the distribution of the first. If contigs are correctly ordered (can be tested with alignment) then we know that the second genome has the mutation.
