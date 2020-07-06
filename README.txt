The source contained here can be used to build programs to test
for embeddability in the spindle surface.  They are based on,
and require, the source for nauty, specifically version nauty27r1,
which we do not provide.

spindle was modified from planarg.  It tests for embeddability
in the spindle surface rather than planarity.  It also tests for
being a topological obstruction or excluded minor for the spindle surface.

gengplanar is just like geng, except it only generates
planar graphs.

gengspin is just like gengspin, except it only
generates graphs that embed in the spindle surface.

splitvg generates all vertex splits (i.e. coextensions) of
a graph.  This is perhaps useful for finding topological obstructions,
particularly with the -d3 flag.

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
    splitvg.c

Now you should be able to make spindle, gengplanar, gengspin,
and splitvg with the rest of nauty.

Files containing all excluded minors and topological obstructions
known to us for the spindle surface are also provided in graph6
format.

--Thomas J. Savitsky (savitsky@gwmail.gwu.edu)
July 2, 2020
