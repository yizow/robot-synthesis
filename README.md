robot-synthesis
===============

07-30-14
Extract files to directory.
Change into that directory.
Run synthesis.py
Results in results.txt

Each line is a list of the xyz-coordinates of the 2nd joint.
Each line should have 15 iterations, i step through angle increments of pi/8
I give the beams lengths of 1,2,3, or 4 in all permutations, so there are 4^4 = 256 lines.

Next i'll figure out a pretty way to visualize and check this data. 