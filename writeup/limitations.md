Limitations 24/10/13
========================================================

- Because the "breaks" in the histogram method you used were a much higher number (10Kb) than the length of the fragments (50-250b), the fragments within the breaks are essentialy orderded randomly with regard to SNP density. The existing rearrangement methods could be a lot more effective if the breaks were smaller and/or fragments were bigger.

Ideas
------

- What if you could just use the Hungarian Algorithm (Munkres) to rearrange the frags into a normal order.
