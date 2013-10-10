Rearrangement Methods
========================================================

What I am trying to acheive:
---------------------------

- I have a list of identifiers for my fragments in their original order, and a list of those same frag ids arranged in order of SNP density (ascending)

- I need to design an algorithm that will rearrange the density ordered frags into their original order, based on SNP density/position data

- This hasn't been done before, so I will need to test a variety of rearrangement methods and I will need a way of ranking them by their effectiveness

- Once I have several methods to compare, I can plot their effectiveness scores on a histogram. The scores of the density order and random rearrangements can be used as controls

Ranking the rearrangement methods
---------------------------------

I have come up with a way of giving each rearrangement method a "Score". See the ruby method: score, in rearrangement_methods.rb. A simple way of defining what this method does is: taking away the index of each fragment in the original order, from it's new index in the rearranged list - this gives the "distance"" that the fragment has moved when re-ordered (as an absolute value) - these "distances" are then summed to give an overall score, which mathematically can be described as the ordinal similarity between the two arrangements. The higher the value of ordinal similarity, the less similar the two orders are. A perfect score would be 0, where the two orders are identical.

Mathematically the ordinal similarity can be defined as...

The highest value of ordinal similarity for my fragments is X... because...

Methods
-------

>### Control 1: Density order

> The first control was to compare the density order against the original. My hypothesis being that this should produce a high ordinal similarity score, because whilst the orginal order of the fragments is based on a normal distribution of SNPs, the density order has the fragments with the most SNPs towards the end of the list, and least near the start. Good rearrangemrent methods should have lower scores than this.

>### Control 2: Random order

> The purpose of this control is to provide a baseline for which all other rearrangement methods should perform better than. To do this, I will create a method that arranges the fragments randomly, and call the score method. This can be repeated a large number of times (e.g. 100) and an average score for randomly ordered fragments can be determined.

Results
-------

>### Ordinal Similarity Scores  (REPLACE WITH HISTOGRAM)

>1. **Density order**: 624,196

>2. **Random order**: 570,00 approx (571,040.44)

>3.

>### Discuss

>1. **Density order**: A high score as expected

>2. **Random order**: The average score of 100 random arragements generated using the shuffle method of the Array Class in ruby. Conversion of integer to float for division (to get average) is the reason for the float score. Calling this multiple times gave scores roughly from 570,000 to 574,000. A more consistent "average random score" could be ascertained with a higher repeat number than 100, but is not neccesary (would increase the running time of the method).
