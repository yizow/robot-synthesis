robot-synthesis
===============

08-06-14
Run synthesis.py
Graphs will show up. 
Red is the crank, blue is the coupler, green is the rake(?), black is the base. 
So there's a definite issue; i thought i had eliminated all structures that don't minimize to 0, but I think i may be looking at the wrong fields of the return result. 

To exit, switch to the terminal and press ctrl+c (or whatever the escape sequence is on windows)
Or you can close 256 graphs. 
Sorry, i'll eventually figure out a better exit system.


07-30-14
Extract files to directory.
Change into that directory.
Run synthesis.py
Results in results.txt

Each line is a list of the xyz-coordinates of the 2nd joint.
Each line should have 15 iterations, i step through angle increments of pi/8
I give the beams lengths of 1,2,3, or 4 in all permutations, so there are 4^4 = 256 lines.

Next i'll figure out a pretty way to visualize and check this data. 