/*
 * activeset.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <stdlib.h>
#include <larsen.h>

/* if item was found in *v, return true  */
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
remove_item (size_t size, int *v, int index)
{
	int		i;
	if (index >= size - 1) return;
	for (i = index; i < size - 1; i++) v[i] = v[i + 1];
	return;
}

/* append a new item to the active set. if it was success, return true */
static bool
activeset_add (larsen *l, int index, int item)
{
	if (l->sizeA >= l->p) return false;
	if (find_item (l->sizeA, l->A, item)) return false;

	add_item (l->sizeA, l->A, index, item);
	l->sizeA++;

	return true;
}

/* remove a item from the active set. if it was success, return true */
static bool
activeset_remove (larsen *l, int index, int item)
{
	if (l->sizeA <= 0) return false;
	if (!find_item (l->sizeA, l->A, item)) return false;

	remove_item (l->sizeA, l->A, index);
	l->sizeA--;

	return true;
}

static int
compare (const void *_a, const void *_b)
{
	int		a = *((int *) _a);
	int		b = *((int *) _b);
	if (a == b) return 0;
	return (a > b) ? 1 : -1;
}

static void
sort_A (size_t size, int *A)
{
	qsort ((void *) A, size, sizeof (int), &compare);
	return;
}

/* complement of active set A */
int *
complementA (larsen *l)
{
	int		i, j, k;
	int		n = l->p - l->sizeA;
	int		A[l->sizeA];
	int		*Ac = (int *) malloc (n * sizeof (int));

	/* copy and sort l->A */
	for (i = 0; i < l->sizeA; i++) A[i] = l->A[i];
	sort_A (l->sizeA, A);

	/* Ac : complement of A */
	for (i = 0, j = 0, k = 0; i < l->p && k < n; i++) {
		if (find_item (l->sizeA - j, &A[j], i)) j++;
		else Ac[k++] = i;
	}
	return Ac;
}

static bool
check_action (ActiveSetAction action)
{
	return (action == ACTIVESET_ACTION_ADD || action == ACTIVESET_ACTION_DROP);
}

static bool
check_index (int index)
{
	return (0 <= index);
}

/* update active set: add / remove a variable */
bool
update_activeset (larsen *l)
{
	bool	status = false;
	if (!check_action (l->oper.action)) return false;
	if (!check_index (l->oper.column_of_X)) return false;
	if (!check_index (l->oper.index_of_A)) return false;

	if (l->oper.action == ACTIVESET_ACTION_ADD)
		status = activeset_add (l, l->oper.index_of_A, l->oper.column_of_X);
	else if (l->oper.action == ACTIVESET_ACTION_DROP)
		status = activeset_remove (l, l->oper.index_of_A, l->oper.column_of_X);

	return status;
}
