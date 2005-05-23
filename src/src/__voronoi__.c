/* Copyright (C) 2000  Kai Habel
**
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with this program; if not, write to the Free Software
** Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
*/

/*
20. Augiust 2000 - Kai Habel: first release
*/

/*
2003-12-14 Rafael Laboissiere <rafael@laboissiere.net>
Added optional second argument to pass options to the underlying
qhull command
*/

/*
extern "C" {
	#include "qhull/qhull_a.h"
}

#include <iostream>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <octave/oct.h>

#ifdef NEED_QHULL_VERSION
char qh_version[] = "__voronoi__.oct 2003-12-14";
#endif
FILE *outfile = stdout;
FILE *errfile = stderr;
char flags[250];
const char *options;
*/


#include <R.h>
#include <Rdefines.h>
#include "qhull_a.h"

/*
DEFUN_DLD (__voronoi__, args, ,
        "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {@var{tri} =} __voronoi__ (@var{point}[, @var{options}])\n\
Internal function for voronoi.\n\
@end deftypefn")
*/
SEXP __voronoi__(const SEXP p, const SEXP options)
{
	SEXP retval;
	retval = R_NilValue;

	FILE *outfile = stdout;      /* output from qh_produce_output() use NULL to skip qh_produce_output() */
	FILE *errfile = stderr;      /* error messages from qhull code */
	char flags[250];             /* option flags for qhull, see qh_opt.htm */

	if(!isString(options) || length(options) != 1){
		error("Second argument must be a single string.");
	}
	if(!isMatrix(p) || !isReal(p)){
		error("First argument should be a real matrix.");
	}

	int i=LENGTH(STRING_ELT(options,0));
	char *opts;
	opts = (char *) R_alloc( ((i>1)?i:1), sizeof(char) );
	strcpy(opts, " ");
	if(i>1) strcpy(opts, CHAR(STRING_ELT(options,0)));

	const int dim = ncols(p);
	const int np  = nrows(p);

	int j=0;
	double *pt_array;
	pt_array = (double *) R_alloc(np*dim, sizeof(double)); 
	for(i=0; i < np; i++)
		for(j=0; j < dim; j++)
			pt_array[dim*i+j] = REAL(p)[i+np*j];
	boolT ismalloc = False;   /* True if qhull should free points in qh_freeqhull() or reallocation */

	int exitcode;             /* 0 if no error from qhull */
	// hmm  lot's of options for qhull here
	sprintf(flags,"qhull v Fv T0 %s",options);
	//outfile = NULL;
	exitcode = qh_new_qhull (dim,np,pt_array,ismalloc,flags,outfile,errfile);
				/*If you want some debugging information replace the NULL
				pointer with outfile.
				*/
	if (!exitcode) {

		facetT *facet;
		vertexT *vertex;
		unsigned int i=0,n=0,k=0,ni[np],m=0,fidx=0,j=0,r=0;
		for (int i=0;i<np;i++) ni[i]=0;
		qh_setvoronoi_all();
		bool infinity_seen = false;
		facetT *neighbor,**neighborp;
		coordT *voronoi_vertex;
		FORALLfacets {
			facet->seen = False;
		}
		FORALLvertices {
			if (qh hull_dim == 3)
				qh_order_vertexneighbors(vertex);
			FOREACHneighbor_(vertex) {
				if (!neighbor->upperdelaunay) {
					if (!neighbor->seen) {
						n++;
						neighbor->seen=True;
					}
					ni[k]++;
				}
			}
			k++;
		}

		Matrix v(n,dim);
		ColumnVector AtInf(np);
		for (int i=0;i < np;i++) AtInf(i)=0;
		octave_value_list F;
		k=0;
		FORALLfacets {
			facet->seen = False;
		}
		FORALLvertices {
			if (qh hull_dim == 3)
				qh_order_vertexneighbors(vertex);
			infinity_seen = false;
			RowVector facet_list(ni[k++]);
			m = 0;
			FOREACHneighbor_(vertex) {
				if (neighbor->upperdelaunay) {
					if (!infinity_seen) {
						infinity_seen = true;
						AtInf(j) = 1;
					}
				} else {
					if (!neighbor->seen) {
						voronoi_vertex = neighbor->center;
						fidx = neighbor->id;
						for (int d=0; d<dim; d++) {
							v(i,d) = *(voronoi_vertex++);
						}
						i++;
						neighbor->seen = True;
						neighbor->visitid = i;
					}
					facet_list(m++)=neighbor->visitid;
				}
			}
			F(r++)=facet_list;
			j++;
		}

		retval(0) = v;
		retval(1) = F;
		retval(2) = AtInf;

		qh_freeqhull(!qh_ALL);
			//free long memory

		int curlong, totlong;
		qh_memfreeshort (&curlong, &totlong);
			//free short memory and memory allocator

		if (curlong || totlong) {
    		    warning("__voronoi__: did not free %d bytes of long memory (%d pieces)", totlong, curlong);
		}
	}
	return retval;
}
