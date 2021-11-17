/*********************************************************************
 * Translated from r4_sort.f in Fortran to C++ by Wesley Miaw
 * <wesley@wesman.net>. Original code by Kevin Olson
 * <olson@jeans.gsfc.nasa.gov> is available from NASA/Goddard at
 * http://sdcd.gsfc.nasa.gov/ESS/exchange/contrib/olson/sph_1d.html.
 *
 * Modified: 05/08/2001
 ********************************************************************/
#include "r4_sort.h"
#include <stdio.h>
#include <limits.h>

static int DEBUG_SORT = 0;

/*********************************************************************
 * The subroutine r4_sort sorts a list of numbers, x. ilist is the
 * final position in the list that each number occupies after sorting
 * and can be used as a permute address for other data.
 ********************************************************************/
void r4_sort(float *x, int *ilist, int *partno, int n) {
    if (DEBUG_SORT > 0) printf("r4_sort(%i)\n", n);

    for(int i = 0; i < n; i++) {
	ilist[i] = i;
    }

    for(i = 0; i < n-1; i++) {
	for(int j = 0; j < n-1-i; j++) {
	    if (x[j+1] < x[j]) {
		float tmp = x[j];
		x[j] = x[j+1];
		x[j+1] = tmp;

		int itmp = partno[j];
		partno[j] = partno[j+1];
		partno[j+1] = itmp;

		itmp = ilist[j];
		ilist[j] = ilist[j+1];
		ilist[j+1] = itmp;
	    }
	}
    }
    return;

    // Create some local variables.
    float *xt = new float[n], *xo = new float[n];
    float size, xmax, xmin;
    int *itemp	= new int[50*n+1],	*ilistt = new int[n], 
		*ix		= new int[n],		*iflag	= new int[n+1];

    // Determine the range of the number list.
    xmax = INT_MIN;
    xmin = (float)INT_MAX;
    for(i = 0; i < n; i++) {
	if (x[i] > xmax) xmax = x[i];
	if (x[i] < xmin) xmin = x[i];
    }
    size = (xmax - xmin) / (50 * n);

    // ???
    for(i = 0; i < n; i++) {
	ix[i] = (int)((x[i] - xmin) / size) + 1;
	iflag[i+1] = 0;
    }
    for(i = 0; i < 50*n+1; i++) {
	itemp[i] = 0;
    }
    for(i = 0; i < n; i++) {
	itemp[ix[i]] = partno[i];
    }
    for(i = 0; i < 50*n+1; i++) {
	if (itemp[i] > 0) iflag[itemp[i]] = 1;
    }
    
    // ???
    int j = 1;
    while(any_eq_zero(iflag, n)) {
	for(i = 0; i < n; i++) {
	    if (iflag[i+1] == 0 &&
		ix[i] + j <= n &&
		itemp[ix[i]+j] == 0)
	    {
		itemp[ix[i]+j] = partno[i];
	    }

	    if (itemp[i] > 0) {
		iflag[itemp[i]] = 1;
	    }

	    if (iflag[i+1] == 0 &&
		ix[i] - j >= 1 &&
		itemp[ix[i]-j] == 0)
	    {
		itemp[ix[i]-j] = partno[i];
	    }

	    if (itemp[i] > 0) {
		iflag[itemp[i]] = 1;
	    }
	}
	j++;
    }

    // Pack the elements of itemp > 0 into ilist.
    pack(itemp, 50*n+1, ilist, n);

    printf("ilist: ");
    for(i = 0; i < n; i++) {
	printf("%i.", ilist[i]);
    }
    printf("\n");

    // Reposition the numbers in the list.
    float *swap_x = new float[n];
    for(i = 0; i < n; i++) {
	swap_x[i] = x[ilist[i]];
    }
    for(i = 0; i < n; i++) {
	x[i] = swap_x[i];
    }

    // Bubble sort to finish.
    for(i = 0; i < n; i++) {
	xt[i] = x[i];
	ilistt[i] = ilist[i];
	if (partno[i] > 1) xo[i] = -1.0;
	else xo[i] = 1.0;
    }
    j = 0;
    int *swap_ilist = new int[n];
    while(any_lt_zero(xo, n)) {
	for(int i = 0; i < n; i++) {
	    if (partno[i] != 0) xo[i] = x[i] - x[partno[i]-1];
	    if (xo[i] < 0 && xo[partno[i]-1] >= 0) {
		swap_x[partno[i]-1] = x[i];
		swap_ilist[partno[i]-1] = ilist[i];
	    }
	}
	for(i = 0; i < n; i++) {
	    x[i] = swap_x[i];
	    ilist[i] = swap_ilist[i];
	}
	for(i = 0; i < n; i++) {
	    if (xo[i] < 0 && xo[partno[i]-1] >= 0) {
		x[i] = xt[partno[i]-1];
		ilist[i] = ilistt[partno[i]-1];
	    }
	}
	for(i = 0; i < n; i++) {
	    xt[i] = x[i];
	    ilistt[i] = ilist[i];
	}
	j++;
    }
}

/*********************************************************************
 * The subroutine pack pushes all elements of the first array > 0 into
 * the second array.
 ********************************************************************/
void pack(int *arr, int arr_size, int *vec, int vec_size) {
    if (DEBUG_SORT > 1) printf("pack(%i,%i)\n", arr_size, vec_size);

    int vec_index = 0;
    for(int i = 0; i < arr_size && vec_index < vec_size; i++) {
	if (arr[i] > 0) {
	    vec[vec_index] = arr[i];
	    vec_index++;
	}
    }
}

/*********************************************************************
 * The subroutine any_eq_zero returns true if any of the elements in
 * the given array equal zero.
 ********************************************************************/
bool any_eq_zero(int *a, int size) {
    for(int i = 1; i < size+1; i++) {
	if (a[i] == 0) return true;
    }
    return false;
}

/*********************************************************************
 * The subroutine any_lt_zero returns true if any of the elements in
 * the given array are less than zero.
 ********************************************************************/
bool any_lt_zero(float *a, int size) {
    for(int i = 0; i < size; i++) {
	if (a[i] < 0) return true;
    }
    return false;
}
