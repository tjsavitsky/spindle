/* helper functions for working with the sparse graph format in the 
   planarity.c source in nauty */
/* author: Thomas J. Savitsky */
/* date started: June 24, 2020 */

#ifndef _SPINUTIL_H_
#define _SPINUTIL_H_

#include "planarity.h"

extern int vertex_deg_va(int, t_ver_sparse_rep *, t_adjl_sparse_rep *);
extern boolean is_spindle_sg(sparsegraph *, int);
extern boolean is_spindle_va(t_ver_sparse_rep *, int, t_adjl_sparse_rep *, int);
extern boolean is_spindle_dg(graph *, int);
extern boolean is_topological_obstruction_sg(sparsegraph *, int);
extern boolean is_topological_obstruction_va(t_ver_sparse_rep *, int,
                                             t_adjl_sparse_rep *, int);
extern boolean is_excluded_minor_sg(sparsegraph *, int);
extern boolean is_excluded_minor_va(t_ver_sparse_rep *, int,
                                             t_adjl_sparse_rep *, int);
extern void va_from_sg(sparsegraph *, t_ver_sparse_rep *, t_adjl_sparse_rep *);
extern boolean is_planar_va(t_ver_sparse_rep *, int, t_adjl_sparse_rep *);
extern void va_copy(t_ver_sparse_rep *, t_ver_sparse_rep *, int,
                    t_adjl_sparse_rep *, t_adjl_sparse_rep *);
extern int contract_edge_va(int, int, t_ver_sparse_rep *, int,
                            t_adjl_sparse_rep *, boolean);
extern int delete_suppress_va(int, int, t_ver_sparse_rep *, int,
                                t_adjl_sparse_rep *);
extern void split_va(t_ver_sparse_rep *, int, t_adjl_sparse_rep *, int,
                    setword, int, t_ver_sparse_rep *, t_adjl_sparse_rep *);

#endif  /* _SPINUTIL_H_ */

