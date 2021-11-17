/*********************************************************************
 * r4_sort.h header file for the r4_sort.cc subroutine. By Wesley Miaw
 * <wesley@wesman.net>. Original r4_sort.f code by Kevin Olson
 * <olson@jeans.gsfc.nasa.gov> available from NASA/Goddard at
 * http://sdcd.gsfc.nasa.gov/ESS/exchange/contrib/olson/sph_1d.html.
 *
 * Modified: 05/08/2001
 ********************************************************************/
#include <math.h>

void r4_sort(double *x, int *ilist, int *partno, int n);
void pack(int *arr, int arr_size, int *vec, int vec_size);
bool any_eq_zero(int *a, int size);
bool any_lt_zero(double *a, int size);
