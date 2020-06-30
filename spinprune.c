/* This should be compiled with gengplanar to ensure WORDSIZE is 32. */

#include "spinutil.h"
#include "nauty.h"

extern int planar_prune_dg(graph *, int, int);
extern int spindle_prune_dg(graph *, int, int);

/*************************************************************************/

static void
va_from_dg(graph *g, t_ver_sparse_rep *V, int n, t_adjl_sparse_rep *A)
{
    int i, j, k, m, maxne;
    set *gi;    
 
    /* geng guarantees that m==1 */
    m = 1;

    j = 0;
    for (i = 0, gi = g; i < n; ++i, gi += (size_t)m)
    {
        if (POPCOUNT(*gi & ALLMASK(n)) == 0)
            V[i].first_edge = NIL;
        else
        {
            V[i].first_edge = j;
            for (k = -1; (k = nextelement(gi,m,k)) >= 0; )
            {
                A[j].end_vertex = k;
                A[j].next = j+1;
                ++j;
            }
            if (A[j-1].end_vertex == i)  /* loops go in twice */
            {
                A[j].end_vertex = i;
                A[j].next = j+1;
                ++j;
            }

        }
        A[j-1].next = NIL;
    }
}

/*************************************************************************/

/* to be used as a PRUNE function in geng
   returns 0 if the graph is planar, and 1 if it is nonplanar */

int
planar_prune(graph *g, int n, int maxn)
{
    DYNALLSTAT(t_ver_sparse_rep,V,V_sz);
    DYNALLSTAT(t_adjl_sparse_rep,A,A_sz);
    int maxne;
    boolean ans;
    
    maxne = (maxn*maxn - maxn)/2;

    DYNALLOC1(t_ver_sparse_rep,V,V_sz,maxn,"spindle");
    DYNALLOC1(t_adjl_sparse_rep,A,A_sz,2*maxne+1,"spindle");

    va_from_dg(g, V, n, A);
 
    ans =  is_planar_va(V, n, A);

    if (ans)    return 0;
    else    return 1;
}

/*************************************************************************/

/* to be used as a PRUNE function in geng
   returns 0 if the graph is spindle, and 1 if it is nonspindle */

int
spindle_prune(graph *g, int n, int maxn)
{
    DYNALLSTAT(t_ver_sparse_rep,V,V_sz);
    DYNALLSTAT(t_adjl_sparse_rep,A,A_sz);
    int maxne;
    boolean ans;
    
    maxne = (maxn*maxn - maxn)/2;

    DYNALLOC1(t_ver_sparse_rep,V,V_sz,maxn,"spindle");
    DYNALLOC1(t_adjl_sparse_rep,A,A_sz,2*maxne+1,"spindle");

    va_from_dg(g, V, n, A);
 
    ans =  is_spindle_va(V, n, A, maxne);

    if (ans)    return 0;
    else    return 1;
}

/*************************************************************************/
