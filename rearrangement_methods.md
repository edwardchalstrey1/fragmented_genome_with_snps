Rearrangement Methods
========================================================

See rearrangement_methods.rb for the code behind each method.

What I am trying to acheive:
---------------------------

- I have a list of identifiers for my fragments in their original order, and a list of those same frag ids arranged in order of SNP density (ascending)

- With real data, the frags would be assigned ids before being put into a fasta file. With my data the ids are named according to the original sequence order, and I have shuffled them before creating a fasta file, so the original order is not simply the same as the order in the fasta file.

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

>### Method 1: Even Odd Method

> Since I know the original ordered sequence has a normal distribution of SNPs, it follows that the when the fragments are in the correct order, their SNP densities follow that same normal curve. Here I assume that the fragments with the highest densities should be in the middle of the distribution, with fragments of lowest densities at either end. The method takes the fragments in the even indexes of the density ordered array and adds these to a new array, so that they are still in ascending order of density. The frags at odd indexes (of the density order) are then reversed (into descending order of density), and then added to the same new array, after (ahead of) the frags from even indexes. 

> I expect this method to perform slightly better than the controls, however, it's basic nature means there may not be much of a difference.

>### Method 2: Left Right Method

> Here I use information derived directly from the fasta and VCF files rather than the SNP densities. I use the SNP position data for each of the fragments, to work out whether the positions are "skewed" to the left or right of each fragment. My hypothesis being that fragments with more SNPs on the "left hand side", are likely to be from the "right side" of the sequence due to the normal distibution of SNPs (and vice versa). The way I work out the "skew": If the sum of the positions is < half of the value produced by multiplying half the length of a fragment by the number of positions it has, this indicates a left skew and I add this fragment to a list of the fragments that should go on the right - and vice versa. I use the Even/Odd Method to add the fragments with 0 SNPs into either the "left" or "right" list. I then add both "left" and "right" into a super-array that is flattened to give the new fragment order. 

> I don't expect this method to perform well, as I have not re-ordered the fragments within the "left" and "right" sub-arrays, which are based on the order that the frags were found in the fasta file (and so are random). It should however perform better than the controls, assuming the "skew" of SNP positions can be used as an indication to the "side" of the normal distribution that a fragment is from.

Results
-------

>### Ordinal Similarity Scores  (REPLACE WITH HISTOGRAM)

>1. **Density order**: 624,196

>2. **Random order**: 570,00 approx (571,040.44)

>3. a. **Even/Odd method**: 492,838    
>   b. **Odd/Even method**: 490,524

>4. **Left Right Method**: 570,434

>### Discuss

>1. **Density order**: A high score as expected

>2. **Random order**: The average score of 100 random arragements generated using the shuffle method of the Array Class in ruby. Conversion of integer to float for division (to get average) is the reason for the float score. Calling this multiple times gave scores roughly from 570,000 to 574,000. A more consistent "average random score" could be ascertained with a higher repeat number than 100, but is not neccesary (would increase the running time of the method).

>3. **Even Odd Method**: As expected this method performed better than the controls, but still had a very high score. I have run it twice, with even indexes of the density order first, then odd indexes reversed and vice versa. The method gives a similar score whichever way around you do it: this is an important thing to note, as I use this method within later methods i.e. there is no need to always test even/odd AND odd/even.

>4. **Left Right Method**: As expected this method performed poorly. However, with a score similar to the average random order, it is not likely that the "skew" of SNP positions is a good indicator of it's position in the original order.
