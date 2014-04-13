/*
 * activeset.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <larsen.h>

/* if item was found in c_vector_int *v, return true  */
static bool
find_item (size_t size, int *v, int item)
{
	int		i;
	for (i = 0; i < size; i++) {
		if (v[i] == item) return true;
	}
	return false;
}

/* append item to c_vector_int *v, return index of last position on A. */
static int
append_item (size_t size, int *v, int item)
{
	int		index = size;
	v[index] = item;
	return index;
}

/* remove item from c_vector_int *v, return index of removed position on A.
 * if specified item was not found in *v, return -1 */
static int
remove_item (size_t size, int *v, int item)
{
	int		i, k;
	int		index = -1;
	k = 0;
	for (i = 0; i < size; i++) {
		if (v[i] == item) {
			index = i;
			continue;
		}
		v[k++] = v[i];
	}
	return index;
}

/* add a new item to the active set. if it was success, return true */
static bool
activeset_add (larsen *l, int item)
{
	int		index;
	if (find_item (l->sizeA, l->A, item)) return false;
	if (!find_item (l->p - l->sizeA, l->Ac, item)) return false;

	if (l->sizeA >= l->p) return false;

	index = append_item (l->sizeA, l->A, item);
	remove_item (l->p - l->sizeA, l->Ac, item);
	l->sizeA++;

	// store the index of item which was added to A
	l->oper.index = index;

	return true;
}

/* remove a item from the active set. if it was success, return true */
static bool
activeset_remove (larsen *l, int item)
{
	int		index;
	if (!find_item (l->sizeA, l->A, item)) return false;
	if (find_item (l->p - l->sizeA, l->Ac, item)) return false;

	if (l->sizeA <= 0) return false;

	index = remove_item (l->sizeA, l->A, item);
	append_item (l->p - l->sizeA, l->Ac, item);
	l->sizeA--;

	// store the index of item which was removed from A
	l->oper.index = index;

	return true;
}

bool
update_activeset (larsen *l)
{
	bool	status = false;
	if (l->oper.action == ACTIVESET_ACTION_ADD) status = activeset_add (l, l->oper.column);
	else if (l->oper.action == ACTIVESET_ACTION_DROP) status = activeset_remove (l, l->oper.column);
	return status;
}
