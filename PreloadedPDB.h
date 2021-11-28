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


#ifndef _PRELOADED_PDB
#define _PRELOADED_PDB

#include <iostream>
#include <map>
#include <vector>

#include "SimpPDB.h"

using namespace std;

class SimPDB;

enum INPUT_FILE_TYPE { UNKNOWN=-1, SILENT_FILE, PDB_LIST };
INPUT_FILE_TYPE filetype(char * filename);
unsigned int num_lines_in_file(char * filename);


/**
 * This class enables two things:
 * 1. Caching of PDB file content in memory, hence reducing disk access
 * 2. Provides support for silentfile
 *
 * To use it:
 * 1. Create an instance,
 * 2. Load it up either with a silent file or with a list of PDB files,
 * 3. Set SimPDB::preloadedPDB to it, and set SimPDB::preloadPDB to true.
 *
 * The use of PreloadedPDB is compulsory for silent files.
 *
 * For a list of decoys, PreloadedPDB is used by default.
 * However, this is bad in some situations. e.g. 36,000 decoys of 100
 * residues each will require about 42M memory. In such situations,
 * using the disk might be better. Hence, PreloadedPDB is off when:
 * 1. Users override it (by switching off SimPDB::preloadPDB)
 * 2. The number of decoys used is more than ADVISED_THRESHOLD.
 */
class PreloadedPDB
{
public:
    static unsigned int ADVISED_THRESHOLD;

private:
    char * silentfilename;
    char * pdblistfilename;

public:
    int mNumResidue;
    int mNumDecoy;
    vector<char *> * mNames;
    map<char *, SimPDB *> filename2PDB; // fix this if it is deemed too slow


public:
    PreloadedPDB();
    ~PreloadedPDB();
    void loadSilentFile(char * silentfilename);
    void loadPDBFromList(char * pdblistfilename);

    SimPDB * getSimPDB(char * pdbfilename); // unused. for internal testing
};


#endif

