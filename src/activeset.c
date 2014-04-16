/*
 * activeset.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <stdlib.h>
#include <larsen.h>

/* if item was found in *v, return true  */
static int
find_item (size_t size, int *v, int item)
{
	int		i;
	for (i = 0; i < size; i++) {
		if (v[i] == item) return i;
	}
	return -1;
}

/* append item to c_vector_int *v, return index of last position on A. */
static void
add_item (size_t size, int *v, int index, int item)
{
	int		i;
	if (index > size) return;
	for (i = size - 1; index <= i; i--) v[i + 1] = v[i];
	v[index] = item;
	return;
}

/* remove item from c_vector_int *v, return index of removed position on A.
 * if specified item was not found in *v, return -1 */
static void
remove_item (size_t size, int *v, int index, int item)
{
	int		i;
	if (index >= size - 1) return;
	for (i = index; i < size - 1; i++) v[i] = v[i + 1];
	return;
}

/* append a new item to the active set. if it was success, return true */
static bool
activeset_append (larsen *l, int item)
{
	int		index = l->sizeA;
	if (l->sizeA >= l->p) return false;
	if (find_item (l->sizeA, l->A, item) >= 0) return false;

	add_item (l->sizeA, l->A, index, item);
	l->sizeA++;

	// store the index of item which was added to A
	l->oper.index_of_A = index;

	return true;
}

/* remove a item from the active set. if it was success, return true */
static bool
activeset_remove (larsen *l, int item)
{
	int		index;
	if (l->sizeA <= 0) return false;
	if ((index = find_item (l->sizeA, l->A, item)) < 0) return false;

	remove_item (l->sizeA, l->A, index, item);
	l->sizeA--;

	// store the index of item which was removed from A
	l->oper.index_of_A = index;

	return true;
}

/* complement of active set A */
int *
complementA (larsen *l)
{
	int		i, k;
	int		n = l->p - l->sizeA;
	int		*Ac = (int *) malloc (n * sizeof (int));

	k = 0;
	for (i = 0; i < l->p && k < n; i++) {
		if (find_item (l->sizeA, l->A, i) < 0) Ac[k++] = i;
	}
	return Ac;
}

bool
update_activeset (larsen *l)
{
	bool	status = false;
	if (l->oper.action == ACTIVESET_ACTION_ADD) status = activeset_append (l, l->oper.column_of_X);
	else if (l->oper.action == ACTIVESET_ACTION_DROP) status = activeset_remove (l, l->oper.column_of_X);
	return status;
}
