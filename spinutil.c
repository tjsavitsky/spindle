/* helper functions for working with the sparse graph format in the 
   planarity.c source in nauty */
/* author: Thomas J. Savitsky */
/* date started: June 24, 2020 */

#define SPINDLE_CHECKS 1

/*************************************************************************/

#include "gtools.h" 
#include "nauty.h"
#include "nausparse.h"
#include "planarity.h"

/* Note:  You can use the following line to out a (V,A) graph to stdout:
   sparseg_adjl_print (V, n, A, FALSE);
 */

/*************************************************************************/

/* from http://rosettacode.org/wiki/Evaluate_binomial_coefficients#C */
static setword
binomial_coef(int n, int k)
{
        setword r;
        int d; 
        r = 1, 
        d = n - k;

        if (k > n) return 0;

        /* choose the smaller of k and n - k */
        if (d > k)
        {
                k = d;
                d = n - k;
        }

        while (n > k)
        {
                if (r >= UINT_MAX / n)
                        return 0; /* overflown */
                r *= n--;
                /* divide (n - k)! as soon as we can to delay overflows */
                while (d > 1 && !(r % d))
                        r /= d--;
        }
        return r;
}

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

extern int
vertex_deg_va(int v, t_ver_sparse_rep *V, t_adjl_sparse_rep *A)
{
    int i, e;
    i = 0;
    e = V[v].first_edge;
    while (e != NIL)
    {
            e = A[e].next;
            ++i;
    }
    return i;
}

/*************************************************************************/

/* Get the (V,A) structure needed for planarity testing from a nauty 
   sparse graph sg. Assumes V and A have been allocated already. */
void
va_from_sg(sparsegraph *sg_ptr, t_ver_sparse_rep *V, t_adjl_sparse_rep *A)
{
    sparsegraph     sg;
    int             n, i, j, k;

    sg = *sg_ptr;
    n = sg.nv;

    k = 0;
    for (i = 0; i < n; ++i)
        if (sg.d[i] == 0)
            V[i].first_edge = NIL;
        else
        {
            V[i].first_edge = k;
            for (j = sg.v[i]; j < sg.v[i]+sg.d[i]; ++j)
            {
                A[k].end_vertex = sg.e[j];
                A[k].next = k+1;
                ++k;
                if (A[k-1].end_vertex == i)  /* loops go in twice */
                {
                    A[k].end_vertex = i;
                    A[k].next = k+1;
                    ++k;
                }
            }
            A[k-1].next = NIL;
        }
}

/*************************************************************************/

/* If (V,A) is nonplanar and get_Ksg is TRUE, then a Kuratowski subgraph
   will be returned in (VK, AK), and its number of edges in nbr_e_obs.
   However, vertices in the subgraph that  not in the original graph will
   have degree 0 in the subgraph, rather than be truly deleted.
   
   If get_Ksg is TRUE, assumes that (VK, AK) and nbr_e_obs have been
   allocated by the caller.
*/
static boolean
boyer_myrvold(t_ver_sparse_rep *V, int n, t_adjl_sparse_rep *A,
                boolean get_Ksg, t_ver_sparse_rep **VK,
                t_adjl_sparse_rep **AK, int *nbr_e_obs)
{
    t_dlcl          **dfs_tree, **back_edges, **mult_edges;
    int             c, edge_pos, vr, wr;
    boolean         ans;
    t_ver_edge      *embed_graph;

    /*
          The input graph is given as an adjacency list:
          V: array of vertices 
          n: size of graph
          A: adjacency list
          e: number of edges
          
      The number of components is returned in c.

      If the graph is nonplanar, edge (vr,wr) is the unembedded edge.
    */
    ans = sparseg_adjl_is_planar(V, n, A, &c,
                                 &dfs_tree, &back_edges, &mult_edges,
                                 &embed_graph, &edge_pos, &vr, &wr);

    if (!ans && get_Ksg)
        embedg_obstruction(V, A, dfs_tree, back_edges,
                           embed_graph, n, &edge_pos,
                           vr, wr, VK, AK, nbr_e_obs);
        
    sparseg_dlcl_delete(dfs_tree, n);
    sparseg_dlcl_delete(back_edges, n);
    sparseg_dlcl_delete(mult_edges, n);
    embedg_VES_delete(embed_graph, n);
 
    return ans;
}

/*************************************************************************/

boolean
is_planar_va(t_ver_sparse_rep *V, int n, t_adjl_sparse_rep *A)
{
    return boyer_myrvold(V, n, A, FALSE, NULL, NULL, NULL);
}

/*************************************************************************/

void
va_copy(t_ver_sparse_rep *fromV, t_ver_sparse_rep *toV, int n,
        t_adjl_sparse_rep *fromA, t_adjl_sparse_rep *toA)
{
    int i, j, e;

    j = 0;
    for (i = 0; i < n; i++)
    {
        if (fromV[i].first_edge == NIL)
        {
            toV[i].first_edge = NIL;
            continue;
        }
        toV[i].first_edge = j;
        e = fromV[i].first_edge;
        while (e != NIL)
        {
            toA[j].end_vertex = fromA[e].end_vertex;
            toA[j].next = j+1;
            e = fromA[e].next;
            ++j;
        }
        toA[j-1].next = NIL;
    }
}

/*************************************************************************/


/* Destroy all multiple-edges incident to v, but leave loops alone. */
static void
destroy_multis(int v, t_ver_sparse_rep *V, int n, t_adjl_sparse_rep *A)
{
    int i, j, e, w;

    for (w = 0; w < n; w++)
    {
        if (w==v) continue;
        j = 0;
        e = V[w].first_edge;
        while (e != NIL)
        {
            if (A[e].end_vertex == v)
                ++j;
            e = A[e].next;
        }

        for (i = 0; i < j-1; i++)
            sparseg_adjl_remove_edge_no_red (V, A, w, v);
    }
}

/*************************************************************************/

/* Contracts the edge (u,v) in the V,A graph.
   Returns FALSE and does nothing if (u,v) is not an edge.
   If (u,v) is an edge, contracts it and returns TRUE.
   This reduces the number of edges and vertices in the graph by one,
   unless a loop was contracted.  Contracting a loop is the
   same as deleting it.

   Returns the number of vertices destroyed.
   
   Note that contracting an edge of a triangle will introduce a multiple edge.
   If nomulti is TRUE, these edges, and multiple edges incident to uv
   in the contraction will be condensed to a single edge.
*/
int
contract_edge_va(int u, int v, t_ver_sparse_rep *V, int n,
                t_adjl_sparse_rep *A, boolean nomultis)
{
    int i, j, w, e;
    boolean isedge;

    if (u > v)
    {
        /* swap u and v */
        u = u + v;
        v = u - v;
        u = u - v;
    }

    /* out of range indices */
    if (u < 0 || v >= n)
        return 0;

    isedge = FALSE;
    e = V[u].first_edge;
    while (e != NIL)
    {
        if (A[e].end_vertex == v)
        {
            isedge = TRUE;
            break;
        }
        e = A[e].next;
    }
    if (!isedge)
        return 0;

    /* contracting a loop is the same as deleting it */
    if (u==v)
    {
        sparseg_adjl_remove_edge_no_red(V, A, u, v);
        return 0;
    }

    for (i=0; i<n; i++)
    {
        e = V[i].first_edge;
        while (e != NIL)
        {
            if (A[e].end_vertex == v)
                A[e].end_vertex = u;
            e = A[e].next;
        }
    }

    /* Note that V[u].first_edge and V[v].first_edge cannot be NIL
       since (u,v) is an edge. */
    e = V[u].first_edge;  
    while (A[e].next != NIL)
        e = A[e].next;
    A[e].next = V[v].first_edge;
    V[v].first_edge = NIL;

    if(!sparseg_adjl_remove_edge_no_red (V, A, u, u))
        gt_abort(">E spindle: contraction failed.\n");

    /* Now delete the isolated vertex v.
       Then renumber every reference to a vertex w with w>v.
       Do this in place, rather than with sparseg_adjl_remove_vertex.
       The resulting (V,A) structure has n-1 vertices. */
    for (w = v; w < n-1; w++)
    {
        V[w] = V[w+1];
    }

    for (w = 0; w < n-1; w++)
    {
        e = V[w].first_edge;
        while (e != NIL)
        {
            if (A[e].end_vertex > v)
                A[e].end_vertex -= 1;
            e = A[e].next;
        }
    }

    if (nomultis)
        destroy_multis(u, V, n-1, A);
    
    return 1;
}

/*************************************************************************/

/* Deletes the edge (u,v) in the V,A graph, and suppresses any degree two
   vertices that result.  Returns the number of suppressed vertices.
   If (u,v) is not an edge, it obviously is not deleted, but u or v
   will be suppressed if they initially had degree 2.
*/
int
delete_suppress_va(int u, int v, t_ver_sparse_rep *V, int n,
                    t_adjl_sparse_rep *A)
{
    int a, b, e;

    if (u > v)
    {
        /* swap u and v */
        u = u + v;
        v = u - v;
        u = u - v;
    }

    /* out of range indices */
    if (u < 0 || v >= n)
        return 0;

    sparseg_adjl_remove_edge_no_red (V, A, u, v);
    
    if (u == v)
        return 0;

    a = b = 0;

    /* contract an edge incident to the vertex with the larger label first */
    if (vertex_deg_va(v, V, A) == 2)
    {
        e = V[v].first_edge;
        a = contract_edge_va(v, A[e].end_vertex, V, n, A, TRUE);
    }

    if (vertex_deg_va(u, V, A) == 2)
    {
        e = V[u].first_edge;
        b = contract_edge_va(u, A[e].end_vertex, V, n-a, A, TRUE);
    }

    /* this works even if (u,v) is a multiedge, and V[v].first_edge == u */

    return a+b;
}


/*************************************************************************/

/* splits (V,A) on vertex i of degree deg, given the subset of its
   neighborhood in mask.  puts the split graph into (Vsplit, Asplit).
   assumes space has already been allocated.  n is the order of (V,A).
   so the split has order n+1. */

void
split_va(t_ver_sparse_rep *V, int n, t_adjl_sparse_rep *A, int i,
            setword mask, int deg, t_ver_sparse_rep *Vsplit,
            t_adjl_sparse_rep *Asplit)
{
    int j, k, e, nprev, iprev;
    DYNALLSTAT(int,switches,switches_sz);

    DYNALLOC1(int,switches,switches_sz,n,"split_va");
    for (j = 0; j < n; j++)
        switches[j] = 0;

    va_copy(V, Vsplit, n, A, Asplit);
    Vsplit[n].first_edge = NIL;

    e = Vsplit[i].first_edge;
    nprev = -1;
    iprev = -1;
    k = 0;
    while (e != NIL)
    {
        if (ISELEMENT(&mask,k))
        {
            if (nprev < 0)
                Vsplit[n].first_edge = e;
            else
                Asplit[nprev].next = e;
            nprev = e;
            /* Later, we will change one entry in Asplit[e].end_vertex's
                neighbor list to point to n instead of i. */
            switches[Asplit[e].end_vertex] += 1;
        }
        else
        {
            if (iprev < 0)
                Vsplit[i].first_edge = e;
            else
                Asplit[iprev].next = e;
            iprev = e;
        }
        e = Asplit[e].next;
        ++k;
    }

    if (nprev >= 0)
        Asplit[nprev].next = NIL;

    if (iprev >= 0)
        Asplit[iprev].next = NIL;
    else
        Vsplit[i].first_edge = NIL;

#ifdef SPINDLE_CHECKS
    j = vertex_deg_va(i, Vsplit, Asplit);
    k = vertex_deg_va(n, Vsplit, Asplit);
    if (j + k != deg)
    {
        fprintf(stderr,">E spindle: split_va error\n");
        fprintf(stderr,"i=%d n=%d deg=%d j=%d k=%d mask=" SETWORD_FORMAT "\n",
                i,n,deg,j,k,mask);
        exit(1);
    }
#endif /* SPINDLE_CHECKS */

    for (j = 0; j < n; j++)
    {
        if (switches[j] == 0)
            continue;
        k = 0;
        e = Vsplit[j].first_edge;
        while (e != NIL)
        {
            if (Asplit[e].end_vertex == i)
            {
                Asplit[e].end_vertex = n;
                ++k;
            }
            if (k == switches[j])
                break;
            e = Asplit[e].next;
        }
    }

    return;
}

/*************************************************************************/

/* returns TRUE if (V,A) embeds in the spindle surface, FALSE otherwise */

boolean
is_spindle_va(t_ver_sparse_rep *V, int n, t_adjl_sparse_rep *A, int ne)
{
    DYNALLSTAT(t_ver_sparse_rep,Vdel,Vdel_sz);
    DYNALLSTAT(t_adjl_sparse_rep,Adel,Adel_sz);

    DYNALLSTAT(t_ver_sparse_rep,Vsplit,V_sz);
    DYNALLSTAT(t_adjl_sparse_rep,Asplit,A_sz);

    DYNALLSTAT(boolean,nosplit,nosplit_sz);

    t_ver_sparse_rep *VKsg;
    t_adjl_sparse_rep *AKsg;
    int i, j, k, e, end, ne_obs, u, d;
    setword mask;

    VKsg = NULL;
    AKsg = NULL;
    /* planar graphs are spindle */
    if (boyer_myrvold(V, n, A, TRUE, &VKsg, &AKsg, &ne_obs))
        return TRUE;

    DYNALLOC1(boolean,nosplit,nosplit_sz,n,"spindle");

    /* There's point in attempting to split a vertex that lies outside of
       some Kuratowski subgraph. */
    for (i = 0; i < n; ++i)
        if (vertex_deg_va(i, VKsg, AKsg) == 0)
            nosplit[i] = TRUE;
        else
            nosplit[i] = FALSE;
    
    /* If some deletion of the graph is planar, then the graph is spindle.
       We need only check deletions of edges in the Kuratowski subgraph
       returned above. */
    DYNALLOC1(t_ver_sparse_rep,Vdel,Vdel_sz,n,"spindle");
    DYNALLOC1(t_adjl_sparse_rep,Adel,Adel_sz,2*ne+1,"spindle");

    k = 0;
    for (i = 0; i < n; ++i)
    {
        e = VKsg[i].first_edge;
        while (e != NIL)
        {
            if (i < AKsg[e].end_vertex)
            /* Only test deletion of edges (u,v) where u<v.
               Deleting a loop can't make a nonplanar graph planar.
            */
            {
                va_copy(V, Vdel, n, A, Adel);

                end = AKsg[e].end_vertex;
                sparseg_adjl_remove_edge_no_red(Vdel, Adel, i, end);

                if (is_planar_va(Vdel, n, Adel))
                    return TRUE;

                ++k;
            }
            e = AKsg[e].next;
        }
    }

    FREES(VKsg);
    FREES(AKsg);

#ifdef SPINDLE_CHECKS
    if (k != ne_obs)
    {
        fprintf(stderr,">E spindle: is_spindle_va error\n");
        fprintf(stderr,"k=%d ne=%d\n",k,ne);
        exit(1);
    }
#endif /* SPINDLE_CHECKS */

    /* if some split of a vertex is planar, then the graph is spindle */

    DYNALLOC1(t_ver_sparse_rep,Vsplit,V_sz,n+1,"spindle");
    DYNALLOC1(t_adjl_sparse_rep,Asplit,A_sz,2*ne+1,"spindle");

    for (i = 0; i < n; ++i)
    {
        if (nosplit[i]) continue;

        /* iterate through half the subsets of N(i) 
           ignore subsets of size 1, which are equivalent to deletions
           for embedabbility.  ignore subsets of size 0. */
        
       /* first destroy any loops */
        while (sparseg_adjl_remove_edge_no_red(V, A, i, i))
            ;

       /* for vertices of degree <4, all splits are equivalent to deletions,
          for the purpose of embeddability */
        d = vertex_deg_va(i, V, A);
        if (d < 4) continue;
        if (d >= WORDSIZE)
        {
            sparseg_adjl_print (V, n, A, FALSE);
            fprintf(stderr,">E spindle: degree too large: %d\n", d);
            exit(1);
        }
    
        if (d % 2 == 0)
            u = d/2-1;
        else
            u = d/2;
        for (j = 2; j <= u; ++j)
        {
            /* set leftmost j bits to 1's */
            mask = ALLMASK(j);
            k = 0;
            while (!ISELEMENT(&mask,d))
            {
                split_va(V, n, A, i, mask, d, Vsplit, Asplit);
                if (is_planar_va(Vsplit, n+1, Asplit))
                    return TRUE;
                ++k;
                mask = next_perm(mask);
            }
#ifdef SPINDLE_CHECKS
            if (k != binomial_coef(d, j))
            {
            fprintf(stderr,">E spindle: error getting bitmasks\n");
            fprintf(stderr,"k=%d d=%d j=%d\n",k,d,j);
            exit(1);
            }
#endif /* SPINDLE_CHECKS */
        }

        /* only do half of these */
        if (d % 2 == 0)
        {
            j = d/2;
            mask = ALLMASK(j);
            for (k = 1; k <= binomial_coef(d, j) / 2; ++k)
            {
                split_va(V, n, A, i, mask, d, Vsplit, Asplit);
                if (is_planar_va(Vsplit, n+1, Asplit))
                    return TRUE;
                mask = next_perm(mask);
            }
        }
    }

    return FALSE;
}

/*************************************************************************/

boolean
is_spindle_sg(sparsegraph *sg_ptr, int nloops)
{
    DYNALLSTAT(t_ver_sparse_rep,V,V_sz);
    DYNALLSTAT(t_adjl_sparse_rep,A,A_sz);
    sparsegraph     sg;
    int             n, ne;

    sg = *sg_ptr;
    n = sg.nv;
    ne = (sg.nde+nloops)/2;

    DYNALLOC1(t_ver_sparse_rep,V,V_sz,n,"spindle");
    DYNALLOC1(t_adjl_sparse_rep,A,A_sz,2*ne+1,"spindle");

    va_from_sg(sg_ptr, V, A);

    return is_spindle_va(V, n, A, ne);
}

/*************************************************************************/

/* returns TRUE if (V,A) is a minimal topological obstruction for
   embeddability in the spindle surface, FALSE otherwise */

boolean
is_topological_obstruction_va(t_ver_sparse_rep *V, int n,
                                t_adjl_sparse_rep *A, int m)
{
    DYNALLSTAT(t_ver_sparse_rep,Vdel,Vdel_sz);
    DYNALLSTAT(t_adjl_sparse_rep,Adel,Adel_sz);
    int             i, j, k, e, end, d;

    /* if the minimum degree of some vertex is < 3, the graph cannot be 
       a minimal obstruction */
    for (i = 0; i < n; ++i)
        if (vertex_deg_va(i, V, A) < 3)
            return FALSE;

    /* the obstruction cannot itself embed in the spindle surface */
    if (is_spindle_va(V, n, A, m))
        return FALSE;

    DYNALLOC1(t_ver_sparse_rep,Vdel,Vdel_sz,n,"spindle");
    DYNALLOC1(t_adjl_sparse_rep,Adel,Adel_sz,2*m+1,"spindle");

    /* test all deletions */
    k = 0;
    for (i=0; i < n; ++i)
    {
        e = V[i].first_edge;
        while (e != NIL)
        {
            if (i < A[e].end_vertex) /* each edge appears twice in A */
            {
                va_copy(V, Vdel, n, A, Adel);

                end = A[e].end_vertex;
                d = delete_suppress_va(i, end, Vdel, n, Adel);

                /* the deletion may have even fewer than m-d-1 edges,
                   if multiedges were destroyed in a delete-suppress */
                if (!is_spindle_va(Vdel, n-d, Adel, m-d-1))
                    return FALSE;
                ++k;
            }
            e = A[e].next;
        }
    }

#ifdef SPINDLE_CHECKS
    if (k != m)
    {
        fprintf(stderr,">E spindle: is_topological_obstruction error\n");
        fprintf(stderr,"k=%d m=%d\n",k,m);
        exit(1);
    }
#endif /* SPINDLE_CHECKS */

    return TRUE;
}

/*************************************************************************/


boolean
is_topological_obstruction_sg(sparsegraph *sg_ptr, int nloops)
{
    DYNALLSTAT(t_ver_sparse_rep,V,V_sz);
    DYNALLSTAT(t_adjl_sparse_rep,A,A_sz);
    sparsegraph     sg;
    int             n, ne;

   /* embedabbility is the same after deleting a loop */
    if (nloops)
        return FALSE;

    sg = *sg_ptr;
    n = sg.nv;
    ne = (sg.nde+nloops)/2;

    DYNALLOC1(t_ver_sparse_rep,V,V_sz,n,"spindle");
    DYNALLOC1(t_adjl_sparse_rep,A,A_sz,2*ne+1,"spindle");

    va_from_sg(sg_ptr, V, A);

    return is_topological_obstruction_va(V, n, A, ne);
}

/*************************************************************************/

/* returns TRUE if (V,A) is an excluded minor
   embeddability in the spindle surface, FALSE otherwise */

boolean
is_excluded_minor_va(t_ver_sparse_rep *V, int n,
                                t_adjl_sparse_rep *A, int ne)
{
    DYNALLSTAT(t_ver_sparse_rep,Vcon,Vcon_sz);
    DYNALLSTAT(t_adjl_sparse_rep,Acon,Acon_sz);
    int i, d, k, e, end;

    if (!is_topological_obstruction_va(V, n, A, ne))
        return FALSE;

    DYNALLOC1(t_ver_sparse_rep,Vcon,Vcon_sz,n,"spindle");
    DYNALLOC1(t_adjl_sparse_rep,Acon,Acon_sz,2*ne+1,"spindle");

    /* test all contractions */
    k = 0;
    for (i=0; i < n; ++i)
    {
        e = V[i].first_edge;
        while (e != NIL)
        {
            if (i < A[e].end_vertex) /* each edge appears twice in A */
            {
                va_copy(V, Vcon, n, A, Acon);

                end = A[e].end_vertex;

                d = contract_edge_va(i, end, Vcon, n, Acon, TRUE);

                if (!is_spindle_va(Vcon, n-d, Acon, ne-1))
                    return FALSE;
                ++k;
            }
            e = A[e].next;
        }
    }

#ifdef SPINDLE_CHECKS
    if (k != ne)
    {
        fprintf(stderr,">E spindle: is_excluded_minor_va error\n");
        fprintf(stderr,"k=%d ne=%d\n",k,ne);
        exit(1);
    }
#endif /* SPINDLE_CHECKS */

    return TRUE;    
}

/*************************************************************************/

boolean
is_excluded_minor_sg(sparsegraph *sg_ptr, int nloops)
{
    DYNALLSTAT(t_ver_sparse_rep,V,V_sz);
    DYNALLSTAT(t_adjl_sparse_rep,A,A_sz);

    sparsegraph     sg;
    int             n, ne;

   /* embedabbility is the same after deleting a loop */
    if (nloops)
        return FALSE;

    sg = *sg_ptr;
    n = sg.nv;
    ne = sg.nde/2;

    DYNALLOC1(t_ver_sparse_rep,V,V_sz,n,"spindle");
    DYNALLOC1(t_adjl_sparse_rep,A,A_sz,2*ne+1,"spindle");

    va_from_sg(sg_ptr, V, A);
    
    return is_excluded_minor_va(V, n, A, ne);
}

/*************************************************************************/
