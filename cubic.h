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
#ifndef _CUBIC_H_
#define _CUBIC_H_

#define cubic_roots(a,b,c,z) cubic_roots2(a,b,c,z)

void cubic_roots1(double a0, double a1, double a2, double * z);
void cubic_roots2(double a0, double a1, double a2, double * z);

#endif
