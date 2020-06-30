spindle is like planarg, but it tests for embeddability
in the spindle surface.  it also tests for topological
and excluded minors for the spindle surface.

genplanarg is just like geng, except it only generates
planar graphs.

gengspin is just like gengspin, except it only
generates graphs that embed in the spindle surface

the script filter_tobs.sh will harness Gnu parallel
to check for topological obstructions 

You can try making the spindle executable by 
itself with the enclosed Makefile.standalone, but you will
need to copy over a few files from the nauty directory.

To build spindle and planarg, first extract nauty.  Then
copy makefile.in to the nauty directory (or apply
the makefile.in.diff patch) and configure.

Then copy the file files to the nauty directory:
	spindle.c
	spinutil.c
	spinutil.h
	spinprune.c

Now you should be able to make spindle and genplanarg
and gengspin with the rest of nauty.

This has only been tested with nauty27r1.

