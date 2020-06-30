spindle is like planarg, but it tests for embeddability
in the spindle surface.  It also tests for topological
and excluded minors for the spindle surface.

gengplanar is just like geng, except it only generates
planar graphs.

gengspin is just like gengspin, except it only
generates graphs that embed in the spindle surface

The script filter_tobs.sh gives an example of
how to harness harness Gnu Parallel
to check for topological obstructions.

You can try making the spindle executable by 
itself with the enclosed Makefile.standalone, but you will
need to copy over a few files from the nauty directory.

To build the programs first extract nauty.  Then
copy makefile.in to the nauty directory (or apply
the makefile.in.diff patch) and configure.

Then copy the file files to the nauty directory:
	spindle.c
	spinutil.c
	spinutil.h
	spinprune.c

Now you should be able to make spindle and gengplanar
and gengspin with the rest of nauty.

These have only been tested with nauty27r1.

