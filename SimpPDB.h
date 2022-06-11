/*
 *  **************************************************************************
 *  Copyright 2009, 2012, 2021 Shuai Cheng Li and Yen Kaow Ng
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


#ifndef _SIMP_PDB
#define _SIMP_PDB

#define LONGEST_CHAIN 4000
#include <iostream>

#include "PreloadedPDB.h"

using namespace std;


int toInt(const string&);
float toFloat(const string&);
void center_residues(float *, int);


class PreloadedPDB;

class SimPDB
{
    public:
      /**
       * These parameters control how PDB files are to be loaded.
       * They do not apply to silent file, which are assumed to be
       * pre-processed.
       */
      static int s_residue;
      static int e_residue;
      static char * chains;
      static vector<string> atom_names;
      static vector<string> atom_matchstrs;

      /**
       * This feature allows the preloading of SimPDB objects.
       * SimPDB will then be obtained from preloadedPDB instead of from disk.
       */
      static PreloadedPDB * preloadedPDB;
      static bool preloadPDB;

    public:
      const char* mProteinFileName;
      int mNumResidue;
      //double mSquaredSum;
      float * mCAlpha;
      int read();
      static int init_atom_names(string namelist, char delimiter);

    public:
      SimPDB(char* aProteinFileName);
      SimPDB(char* aProteinFileName, int len);
      ~SimPDB();

      // Special constructors used only by PreloadedPDB. Don't touch.
      SimPDB();
      SimPDB(int len);
};

#endif
