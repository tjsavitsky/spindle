/* Split all vertices of a graph in all possible ways.  Splitting a 
   vertex is the reverse operation of contracting an edge.
   By Thomas J. Savitsky.  Date started July 5, 2020. */
/* modified from addedgeg.c   nauty version 2.6; B D McKay, Jan 2013. */

#define USAGE "splitvg [-lq] [-d#] [infile [outfile]]"

#define HELPTEXT \
" For each vertex v, output all ways to split v, subject to certain\n\
  conditions.  Each subset S of N(v) correspondends to a split of v as\n\
  follows.  First, add a new vertex w=n.  Then for each x in S, remove\n\
  edge vx, and add edge wx. Finally, add the edge vw.\n\
\n\
    The output file has a header if and only if the input file does.\n\
\n\
    -l  Canonically label outputs\n\
    -d# Ensure the splitting does not introduce vertices of with degree\n\
        less than d\n\
    -q  Suppress auxiliary information\n\
\n\
  This may not work as you expect if the graph has loops.\n"

/*************************************************************************/

#include "gtools.h" 
#include "gutils.h"

/*************************************************************************/

/* modified from Bit Twiddling hacks to represent permutations in the
   opposite order */
static setword
next_perm(setword v)
{
    int a, b;
    
    a = FIRSTBITNZ( v );
    b = FIRSTBITNZ( ( ~(v << a) ) );
    ADDELEMENT( &v, ( a + b ) );
    v &= BITMASK( a + b - 1);
    v |= ALLMASK( b - 1 );

    return v;
}

/*************************************************************************/

/* Copies graph g to gcopy. If the number of vertices in gcopy
   exceeds that in g, then the excess vertices are isolated in the copy.
   Assumes space for the graphs has already been allocated. */
static void
copy_dg(graph *g, int n, int m, graph *gcopy, int ncopy, int mcopy)
{
    int i,v;
    graph *gv, *gcopyv;

    if (ncopy < n || mcopy < m)
    {
        fprintf(stderr, ">E error in copy_dg. aborting.\n");
        exit(1);
    }

    for (v = 0, gv = g, gcopyv = gcopy; v < n; ++v, gv += m, gcopyv += mcopy)
    {
        for (i = 0; i < m; ++i)
            gcopyv[i] = gv[i];
        for (i = m; i < mcopy; ++i)
            gcopyv[i] = 0;
    }

    for (v = n, gcopyv = gcopy + n*mcopy; v < ncopy; ++v, gcopyv += mcopy)
        for (i = 0; i < mcopy; i++)
            gcopyv[i] = 0;

}

/*************************************************************************/

int
main(int argc, char *argv[])
{
    char *infilename,*outfilename;
    FILE *infile,*outfile;
    boolean badargs,dolabel,quiet,dswitch;
	int i,j,m,n,v,w,nsplit,msplit,argnum;
	int codetype,outcode;
	graph *g, *gq;
	nauty_counter nin,nout;
    char *arg,sw;
	setword *gv, *gsw, *gsv, mask;
	int mindeg,degv;
	double t;
#if MAXN
	graph h[MAXN*MAXM];
	int deg[MAXN], neighbors[MAXN];
    graph gsplit[MAXN*MAXM];
#else
	DYNALLSTAT(graph,h,h_sz);
	DYNALLSTAT(int,deg,deg_sz);
	DYNALLSTAT(int,neighbors,neighbors_sz);
	DYNALLSTAT(graph,gsplit,gsplit_sz);
#endif

	HELP; PUTVERSION;

    infilename = outfilename = NULL;
	dswitch = dolabel = quiet = FALSE;

	argnum = 0;
	badargs = FALSE;
	for (j = 1; !badargs && j < argc; ++j)
	{
	    arg = argv[j];
	    if (arg[0] == '-' && arg[1] != '\0')
	    {
		++arg;
		while (*arg != '\0')
		{
		    sw = *arg++;
		         SWBOOLEAN('l',dolabel)
		    else SWBOOLEAN('q',quiet)
            else SWINT('d',dswitch,mindeg,">E splitvg -d")
		    else badargs = TRUE;
		}
	    }
	    else
	    {
		++argnum;
		if      (argnum == 1) infilename = arg;
        else if (argnum == 2) outfilename = arg;
		else                  badargs = TRUE;
	    }
	}

	if (badargs)
	{
	    fprintf(stderr,">E Usage: %s\n",USAGE);
	    GETHELP;
	    exit(1);
	}

	if (!quiet)
	{
	    fprintf(stderr,">A splitvg");
	    if (dolabel) fprintf(stderr," -l");
	    if (dswitch) fprintf(stderr," -d%d",mindeg);
	    if (argnum > 0) fprintf(stderr," %s",infilename);
	    if (argnum > 1) fprintf(stderr," %s",outfilename);
	    fprintf(stderr,"\n");
	    fflush(stderr);
	}

	if (infilename && infilename[0] == '-') infilename = NULL;
	infile = opengraphfile(infilename,&codetype,FALSE,1);
	if (!infile) exit(1);
	if (!infilename) infilename = "stdin";

	if (!outfilename || outfilename[0] == '-')
	{
	    outfilename = "stdout";
	    outfile = stdout;
	}
	else if ((outfile = fopen(outfilename,"w")) == NULL)
	{
	    fprintf(stderr,"Can't open output file %s\n",outfilename);
	    gt_abort(NULL);
	}

	if (codetype&SPARSE6) outcode = SPARSE6;
	else                  outcode = GRAPH6;

	if (codetype&HAS_HEADER)
	{
	    if (outcode == SPARSE6) writeline(outfile,SPARSE6_HEADER);
	    else    		    writeline(outfile,GRAPH6_HEADER);
	}

	if (!dswitch) mindeg = 1;
    if (mindeg < 1) mindeg = 1;

	if (dolabel) nauty_check(WORDSIZE,1,1,NAUTYVERSIONID);

	nin = nout = 0;
	t = CPUTIME;
	while (TRUE)
	{
	    if ((g = readg(infile,NULL,0,&m,&n)) == NULL) break;
	    ++nin;

#if !MAXN
	    DYNALLOC1(int,deg,deg_sz,n,"splitvg");
	    DYNALLOC1(int,neighbors,neighbors_sz,n,"splitvg");
#endif

        nsplit = n+1;
        msplit = SETWORDSNEEDED(nsplit);

#if MAXN
        if (nsplit > MAXN)
        {
            fprintf(stderr, ">E nsplit exceeds MAXN. aborting.\n");
            exit(1);
        }
#else
        DYNALLOC2(graph,gsplit,gsplit_sz,nsplit,msplit,"splitvg");
#endif

	    for (v = 0, gv = g; v < n; ++v, gv += m)
	    {
    		degv = 0;
	    	for (i = 0; i < m; ++i)
	    	    degv += POPCOUNT(gv[i]);
	    	deg[v] = degv;
	    }

	    for (v = 0, gv = g; v < n; ++v, gv += m)
	    {
            if (deg[v]/2 + 1 < mindeg)
                continue;

            if (deg[v] >= WORDSIZE)
            {
                fprintf(stderr,">E splitvg: degree too large: %d\n", deg[v]);
                exit(1);
            }

            j = 0;
            for (i = 0; i < n; ++i)
            {
                if (ISELEMENT(gv,i))
                {
                    neighbors[j] = i;
                    ++j;
                }
            }

            for (i = mindeg-1; i <= deg[v]-mindeg+1; ++i)
            {
                /* set leftmost i bits to 1's */
                mask = ALLMASK(i);
                while (!ISELEMENT(&mask,deg[v]))
                {
                    copy_dg(g,n,m,gsplit,nsplit,msplit);
                    for (j = 0; j < deg[v]; ++j)
                    {
                        if (ISELEMENT(&mask,j))
                            continue;
                        w = neighbors[j];
                        gsw = GRAPHROW(gsplit,w,msplit);
                        gsv = GRAPHROW(gsplit,v,msplit);
                        DELELEMENT(gsv,w);
                        DELELEMENT(gsw,v);
                        ADDONEEDGE(gsplit,w,nsplit-1,msplit);
                    }
                    ADDONEEDGE(gsplit,v,nsplit-1,msplit);
                    gq = gsplit;
	                if (dolabel)
	                {
#if !MAXN
	    	            DYNALLOC2(graph,h,h_sz,nsplit,msplit,"splitvg");
#endif
	     	            fcanonise(gsplit,msplit,nsplit,h,NULL,FALSE);
                            /*FIXME (loops)*/
                        gq = h;
	                }
    	            if (outcode == SPARSE6) writes6(outfile,gq,msplit,nsplit);
	                else                    writeg6(outfile,gq,msplit,nsplit);
	    	        ++nout;
                    if (i == 0) break;
                    mask = next_perm(mask);
                }
            }       
	    }
	    FREES(g);
	}
	t = CPUTIME - t;

        if (!quiet)
            fprintf(stderr,
              ">Z  " COUNTER_FMT " graphs read from %s, "
                           COUNTER_FMT " written to %s; %3.2f sec.\n",
                    nin,infilename,nout,outfilename,t);

	 exit(0);
}
