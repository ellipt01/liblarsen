/*
 * activeset.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <larsen.h>

/* if item was found in cl_vector_int *v, return true  */
static bool
find_item (cl_vector_int *v, int item)
{
	int		i;
	for (i = 0; i < v->size; i++) {
		int		vi = cl_vector_int_get (v, i);
		if (vi == item) return true;
	}
	return false;
}

/* append item to cl_vector_int *v, return index of last position on A. */
static int
append_item (cl_vector_int *v, int item)
{
	int		index = v->size;
	v->size++;
	cl_vector_int_set (v, index, item);
	return index;
}

/* remove item from cl_vector_int *v, return index of removed position on A.
 * if specified item was not found in *v, return -1 */
static int
remove_item (cl_vector_int *v, int item)
{
	int		i, k;
	int		index = -1;
	k = 0;
	for (i = 0; i < v->size; i++) {
		int		vi = cl_vector_int_get (v, i);
		if (vi == item) {
			index = i;
			continue;
		}
		cl_vector_int_set (v, k++, vi);
	}
	if (index >= 0) v->size--;
	return index;
}

/* add a new item to the active set. if it was success, return true */
bool
activeset_add (larsen *l, int item)
{
	int		index;
	if (find_item (l->A, item)) return false;
	if (!find_item (l->Ac, item)) return false;

	if (l->A->size >= l->x->size2) return false;
	if (l->Ac->size <= 0) return false;

	index = append_item (l->A, item);
	remove_item (l->Ac, item);

	// store the index of item which was added to A
	l->oper.index = index;

	return true;
}

/* remove a item from the active set. if it was success, return true */
bool
activeset_remove (larsen *l, int item)
{
	int		index;
	if (!find_item (l->A, item)) return false;
	if (find_item (l->Ac, item)) return false;

	if (l->A->size <= 0) return false;
	if (l->Ac->size >= l->x->size2) return false;

	index = remove_item (l->A, item);
	append_item (l->Ac, item);

	// store the index of item which was removed from A
	l->oper.index = index;

	return true;
}
