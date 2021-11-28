/*
 *  **************************************************************************
 *  Copyright 2012, 2021 Shuai Cheng Li and Yen Kaow Ng
 *  **************************************************************************
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  **************************************************************************
 * 
 */

#include <cmath>
#include <iostream>

using namespace std;


/**
 * Computes the eigenvalues and (normalized) eigenvectors of a real symmetric
 * matrix a[3][3] in d[3] and v[3][3] respectively.
 *
 * Only the upper-right triangle of a[3][3] is referenced and over-written.
 *
 * (The method follows closely the discussions in "Numerical Recipes in C".)
 */
bool jacobi3(double a[3][3], double d[3], double v[3][3])
{
#define NUM_ITER 50
#define SMALL_NUM_ITER 4 
    // initialize v to the identify matrix
    v[0][0] = 1.0;
    v[0][1] = 0.0;
    v[0][2] = 0.0;
    v[1][0] = 0.0;
    v[1][1] = 1.0;
    v[1][2] = 0.0;
    v[2][0] = 0.0;
    v[2][1] = 0.0;
    v[2][2] = 1.0;

    // initialize d to the diagonal of a
    d[0] = a[0][0];
    d[1] = a[1][1];
    d[2] = a[2][2];

    for (int iter=0; iter < NUM_ITER; iter++)
    {
        // Find sum of non-diagonal entries in a
        double sum_a = 0.0;
        for (int i=0; i < 2; i++)
            for (int j=i+1; j < 3; j++)
                sum_a += fabs(a[i][j]);
        if (sum_a == 0.0) // Stop if sum is small
            return true;

        double threshold = (iter < SMALL_NUM_ITER)? 0.2 * sum_a / 9.: 0.;

        /**
         * i and j traverse only the upper-right part ('o') of a.
         *
         * i 0  x o o
         *   1  x x o
         *   2  x x x
         *      0 1 2
         *      j
         */
        for (int i=0; i < 2; i++)
        {
            for (int j=i+1; j < 3; j++)
            {
                double g = 100.0 * fabs(a[i][j]);
                double di = fabs(d[i]);
                double dj = fabs(d[j]);

                if (iter > SMALL_NUM_ITER && di+g == di && dj+g == dj)
                {
                    a[i][j] = 0.0; // just clean it up if it is too small
                }
                /**
                 * Use an almost diagonal matrix M to rotate a (in an axis
                 * that affects only row/column i and j) in order to eliminate
                 * a[i][j]. (Hence this elimination process is similar to
                 * Gaussian, but uses rotation instead of subtraction.)
                 *
                 * The rotation matrix M has zero for all off-diagonal
                 * elements, except for the four elements
                 *      M[i][j] = -M[j][i] = s
                 *      M[i][i] =  M[j][j] = d
                 *
                 * Under rotation M, a becomes
                 *   1. a'[k][i] = c   a[k][i] - s   a[k][j], k != i or j
                 *   2. a'[k][j] = c   a[k][j] + s   a[k][i], k != i or j
                 *   3. a'[i][i] = c^2 a[i][i] + s^2 a[j][j] - 2sc a[i][j]
                 *   4. a'[j][j] = s^2 a[i][i] + c^2 a[j][j] + 2sc a[i][j]
                 *   5. a'[i][j] = (c^2-s^2) a[i][j] + sc (a[i][i]-a[j][j])
                 *
                 * Since c and s are to be chosen so that a'[i][j] = 0,
                 *   5. 0 = (c^2-s^2) a[i][j] + sc (a[i][i]-a[j][j])
                 *   => (c^2-s^2) / sc = a[i][i]-a[j][j] / a[i][j]
                 * Let theta = cot(2 phi)
                 *   => theta = (c^2-s^2)/2sc = 0.5 * a[i][i]-a[j][j] /a[i][j]
                 * If we let t=s/c, then
                 *   =>  t^2 + 2 t theta - 1 = 0
                 *   =>  t = sign(theta) / (|theta| + sqrt(theta^2+1))
                 * When theta^2 >> 1, sqrt(theta^2+1) --> theta, then
                 *   =>  t = sign(theta) / 2 theta
                 *   =>  c = 1 / sqrt(t^2 + 1)
                 *
                 * Then, the equations become
                 *   1. a'[k][i] = a[k][i] - s(a[k][j] + tau a[k][i]),
                 *   2. a'[k][j] = a[k][j] + s(a[k][i] - tau a[k][j]),
                 *   3. a'[i][i] = a[i][i] - t a[i][j],
                 *   4. a'[j][j] = a[i][i] + t a[i][j],
                 * where tau = s/1+c.
                 */
                else if (fabs(a[i][j]) > threshold)
                {
                    double diff = d[j]-d[i]; // a[i][i]-a[j][j]
                    double t;

                    if (fabs(diff) + g == fabs(diff))
                        t = a[i][j] / diff;
                    else
                    {
                        double theta = 0.5 * diff / a[i][j];
                        t = 1.0 / ( fabs(theta) + sqrt(theta*theta + 1) );
                        if (theta < 0.0)
                            t = -t;
                    }
                    double c = 1.0 / sqrt(t*t + 1);
                    double s = t*c;
                    double tau = s / (1.0 + c);
                    double h = t*a[i][j];
                    d[i] -= h;
                    d[j] += h;
                    a[i][j] = 0.0;
#define ROTATE(a,x1,y1,x2,y2)   g = a[x1][y1],h = a[x2][y2], \
                                a[x1][y1] = g-s*(h+tau*g),   \
                                a[x2][y2] = h+s*(g-tau*h)
                    /**
                     * Since we only do the upper-right elements, we split the
                     * matrix manipulation into three parts.
                     */
                    for (int k=0; k < i; k++)
                        //         i   j
                        //         |   |
                        //         x   x
                        //     i---|---|---
                        //         |   |
                        //     j---|---|---
                        //         |   |
                        ROTATE(a, k, i, k, j);
                    for (int k=i+1; k < j; k++)
                        //         i   j                   i   j
                        //         |   |                   |   |
                        //     i---|-x-|--- (that is,) i---|-x-|---
                        //         |   x    (same as:)     |   |
                        //     j---|---|---            j---|-x-|---
                        //         |   |                   |   |
                        ROTATE(a, i, k, k, j);
                    for (int k=j+1; k < 3; k++)
                        //         i   j
                        //         |   |
                        //     i---|---|-x-
                        //         |   |
                        //     j---|---|-x-
                        //         |   |
                        ROTATE(a, i, k, j, k);
                    // Do the entire output matrix v (can we save computation
                    // here by doing this only for the final loop?)
                    for (int k=0; k < 3; k++)
                        ROTATE(v, k, i, k, j);
#undef ROTATE
                }
            }
        }
    }
    cout << "More than " << NUM_ITER << " iterations in jacobi3()" << endl;
    return false;
#undef NUM_ITER
#undef SMALL_NUM_ITER
}


#ifdef _TEST_JACOBI_
void print_matrix(double a[3][3])
{
    for (int i=0; i < 3; i++)
    {
        for (int j=0; j < 3; j++)
            cout << a[i][j] << ",";
        cout << endl;
    }
}

int main(void)
{
    double a[3][3] = {{1.1, 2.5, 5.}, {2., 3.1, 2.0}, {5., 2.5, 1.1}}; 
    double v[3][3];
    double d[3];
    jacobi3(a, d, v);
    print_matrix(a);
    print_matrix(v);
    cout << d[0] << "," << d[1] << "," << d[2] << endl;
}
#endif


