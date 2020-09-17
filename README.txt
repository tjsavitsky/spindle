The source contained here can be used to build programs to test
for embeddability in the spindle surface.  They are based on, and
require, the source for nauty, specifically version nauty27r1, which
we do not provide.  We point to the following website for the latest
release of nauty: http://pallini.di.uniroma1.it/

spindle was modified from planarg.  It tests for embeddability in
the spindle surface rather than planarity.  It also tests for being
a topological obstruction or excluded minor for the spindle surface.

gengplanar is just like geng, except it only generates planar graphs.

gengspin is just like gengspin, except it only generates graphs that
embed in the spindle surface.

splitvg generates all vertex splits (i.e. coextensions) of a graph.
This is perhaps useful for finding topological obstructions,
particularly with the -d3 flag.

To build the programs first extract nauty.  Then copy makefile.in
to the nauty directory (or apply the makefile.in.diff patch) and
configure.

Then copy the following files to the nauty directory:
	spindle.c spinutil.c spinutil.h spinprune.c splitvg.c

Now you should be able to make spindle, gengplanar, gengspin, and
splitvg with the rest of nauty.

You can instead build the spindle executable by itself with the
enclosed Makefile.standalone, but you will need to copy over a few
files from the nauty directory.

Our typical use case is to find all topological obstructions for the
spindle surface on n vertices and m vertices by by running spindle with
the -t flag on the set of graphs outputted by nauty's geng program.
For example, the shell command
	$ geng -d3 10 20 | spindle -t
finds the 29 topological obstructions for the spindle surface with
10 vertices and 20 edges.

Since every excluded minor is also a topological obstruction, to obtain
a complete collection of excluded minors for the spindle surface on n
vertices and m edges, one needs only feed the topological obstructions
on n vertices and m edges to the spindle program with the -e flag.

The script filter_tobs.sh gives an example of how to harness GNU
Parallel to check for topological obstructions.

Files containing all excluded minors and topological obstructions
for the spindle surface known to us are provided in graph6 format.

--Thomas J. Savitsky (savitsky@gwmail.gwu.edu)
--Steven A. Schluchter (steven.schluchter@gmail.com)
August 21, 2020
