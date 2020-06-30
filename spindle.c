/* spindle : test if a graph embeds in the spindle surface 
   will not work unless the graph is simple.
   modified from the planarg.c file from nauty. */
/* author: Thomas J. Savitsky */
/* date started: June 24, 2020 */

#define USAGE "spindle [-nuV] [-t|-e] [-s] [infile [outfile]]"

#define HELPTEXT \
" For each input, write to output if it embeds spindle surface.\n\
    -n  Negate the search.  Print non-matches instead of matches.\n\
    -u  Don't write anything, just count\n\
    -V  Write report on every input\n\
    -t  Instead, write to output if the graph is a minimal topological\n\
        obstruction for spindle embeddability.\n\
    -e  Instead, write to output if the graph is an excluded minor for\n\
        spindle embeddability.\n\
    \n\
    Works for graphs with loops and multiple edges.\n\
    Will abort if given a nonplanar graph with a vertex of degree higher\n\
    than WORDSIZE. In the worst case, running time is exponential in\n\
    the number of edges.\n\
"

/*************************************************************************/

#include "nauty.h"
#include "gtools.h" 
#include "spinutil.h"

/*************************************************************************/

int
main(int argc, char *argv[])
{
    char *infilename,*outfilename;
    FILE *infile,*outfile;
    sparsegraph sg;
    boolean badargs;
    boolean verbose,nonmatch,quiet,nowrite;
    boolean obstruction,excluded;
    boolean ans;
    int i,j,k,n,argnum,ne,nloops;
    int codetype,outcode;
    nauty_counter nin,nout,nmatch;
    char *arg,sw;
    double t0,tyes,tno;

    HELP; PUTVERSION;

    infilename = outfilename = NULL;
    quiet = nowrite = verbose = nonmatch = FALSE;
    obstruction = excluded = FALSE;

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
                     SWBOOLEAN('n',nonmatch)
                else SWBOOLEAN('V',verbose)
                else SWBOOLEAN('u',nowrite)
                else SWBOOLEAN('t',obstruction)
                else SWBOOLEAN('e',excluded)
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
        exit(1);
    }

    if (obstruction && excluded)
    {
        fprintf(stderr,">E spindle: pick at most one of -t and -e\n");
        exit(1);
    }

    if (!quiet)
    {
        fprintf(stderr,">A spindle");
        if (nonmatch||nowrite||verbose)
            fprintf(stderr," -");
        if (nonmatch) fprintf(stderr,"n");
        if (nowrite) fprintf(stderr,"u");
        if (verbose) fprintf(stderr,"V");
        if (obstruction||excluded)
            fprintf(stderr," -");
        if (obstruction) fprintf(stderr,"t");
        if (excluded) fprintf(stderr,"e");

        if (argnum > 0) fprintf(stderr," %s",infilename);
        if (argnum > 1) fprintf(stderr," %s",outfilename);
        fprintf(stderr,"\n");
        fflush(stderr);
    }

    if (infilename && infilename[0] == '-') infilename = NULL;
    infile = opengraphfile(infilename,&codetype,FALSE,1);
    if (!infile) exit(1);
    if (!infilename) infilename = "stdin";

    NODIGRAPHSYET(codetype);

    if (!nowrite)
    {
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

        if (codetype&SPARSE6)   outcode = SPARSE6;
        else                       outcode = GRAPH6;

        if (codetype&HAS_HEADER)
        {
            if (outcode == SPARSE6) writeline(outfile,SPARSE6_HEADER);
            else                    writeline(outfile,GRAPH6_HEADER);
        }
    }

    nin = nout = nmatch = 0;
    SG_INIT(sg);

    tyes = tno = 0.0;
    while (TRUE)
    {
        if (read_sg_loops(infile,&sg,&nloops) == NULL) break;

        n = sg.nv;
        /* For example, the 4-cycle has sg.nde=8.
           The 4-cycle with a loop at one vertex has sg.nde=9 and loops=1.
           So sg.nde+nloops equals the degree sum of the graph. */
        ne = (sg.nde+nloops)/2;
        ++nin;
        t0 = CPUTIME;

        if (obstruction)
            ans = is_topological_obstruction_sg(&sg, nloops);
        else if (excluded)
            ans = is_excluded_minor_sg(&sg, nloops);
        else
            ans = is_spindle_sg(&sg, nloops);

        if (ans)
        {
            ++nmatch;
            tyes += CPUTIME - t0;
            if (!nowrite && !nonmatch)
            {
                writelast(outfile);
                ++nout;
            }
            if (verbose)
                fprintf(stderr,"graph " COUNTER_FMT ": n=%d ne=%d matches\n",
                        nin,n,ne);
        }
        else
        {
            tno += CPUTIME - t0;
            if (!nowrite && nonmatch)
            {
                writelast(outfile);
                ++nout;
            }
            if (verbose)
                fprintf(stderr,"graph " COUNTER_FMT ": n=%d ne=%d does not match\n",
                        nin,n,ne);
        }
    }

    if (!nowrite)
    {
        if (!quiet)
            fprintf(stderr,
            ">Z  " COUNTER_FMT " graphs read from %s, "
                  COUNTER_FMT " written to %s; %3.2f sec.\n",
                nin,infilename,nout,outfilename,tyes+tno);
    }
    else
    {
        fprintf(stderr,
                " " COUNTER_FMT " graphs input\n "
                    COUNTER_FMT " graphs matched\n",nin,nmatch);
        fprintf(stderr," cpu = %3.2f sec. ",tyes+tno);
        fprintf(stderr,"\n");
    }

    return 0;
}
