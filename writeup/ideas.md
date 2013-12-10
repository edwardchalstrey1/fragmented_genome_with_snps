Rearrangement method ideas
========================================================

Machine learning
-----

Learning involves giving the machine examples, then the more examples, the better the machine gets at telling how a new thing should be catagorized. For example, a machine learning system could be trained on email messages to learn to distinguish between spam and non-spam messages.

Hill Climbing algorithms
------

An arbitrary solution to a problem is generated, then an element of the solution is changed, if this produces a better solution, it is kept. Iteratively, the solution gets better. It can however, reach a local maxima, where changing one element cannot make a better solution, whilst the global maxima (the best possible solution) is not actually reached.

What to do
--------

Looking at the [NGM website](http://bar.utoronto.ca/ngm/cgi-bin/emap.cgi), with sample data, we can see that there is a noticeable homo/hetero signal distribution around the causative mutation.

You need to modify your model genome:

- Two sets of SNPs (hetero and homozygous), that will generate a normally distributed (or similar curve to NGM sample) ratio distribution
- They should be 2 "chastity threads", i.e. values around 0.5 for hetero, and 1.0 for homozygous (and different to reference base in NGM)
- The ratio is maximized at the causative mutation, and distributed around it
- The standard deviation/ curve of the distribution is down to the kernel density. How tightly does this fit around the data? You decide: we want the distribution to look as normal as possible, as this will be what we try and regenerate when putting the fragments back together
- Use a qq plot to make sure the ratio is a normal distribution
- then use the qq plot to test the distribution of rearranged fragments' hom/het ratio (score the rearrangements based on how close to the norm dist they are)