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


#include <iostream>
#include <fstream>

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>

#include "SimpPDB.h"
#include "PreloadedPDB.h"

using namespace std;


/**
 * Reasons for choosing 36,000 as the threshold for using disk access
 * 1. 36,000 decoys of 100 residues each is about 44MB, which is still
 *    acceptable for a workstation PC. If user has less capabled hardware,
 *    the "-d" switch may be used.
 * 2. For more than 36,000 decoys, the bottleneck is probably not in disk.
 */
unsigned int PreloadedPDB::ADVISED_THRESHOLD = 36000;


INPUT_FILE_TYPE
filetype(char * filename)
{
    char buf[400];
    ifstream input(filename);
    string line;
    if (!input)
    {
        cerr << "Can't open input file \"" << filename << "\"" << endl;
        exit(0);
    }

    input.getline(buf, 400);
    line=buf;
    if (line.substr(0, 9)=="SEQUENCE:")
    {
        input.close();
        return SILENT_FILE;
    }

    input.seekg(0);
    input.getline(buf, 400);
    char *token = strtok(buf, " ");
    if (token == NULL)
    {
        input.close();
        return UNKNOWN;
    }
    char *name = new char[strlen(token)+1];
    strcpy(name, token);

    ifstream pdbfile(name); // check if name is a file
    if (pdbfile)
    {
        pdbfile.close();
        input.close();
        return PDB_LIST;
    }
    pdbfile.close();

    input.close();

    return UNKNOWN;
}


unsigned int num_lines_in_file(char * filename)
{
    char buf[400];
    ifstream input(filename);
    if (!input)
    {
        cerr << "Can't open input file \"" << filename << "\"" << endl;
        exit(0);
    }
    unsigned int count = 0;
    while (!input.eof())
    {
        input.getline(buf, 400);
        count++;
    }
    return count-1;
}


//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

PreloadedPDB::PreloadedPDB()
{
}


/**
 * Populate the PreloadedPDB with a silentfile
 */
void
PreloadedPDB::loadSilentFile(char * filename)
{
    char buf[400];
    ifstream input(filename);
    string line;
    if (!input)
    {
        cerr << "Can't open silent file \"" << filename << "\"" << endl;
        exit(0);
    }

    silentfilename = filename;
    pdblistfilename = NULL;

    /**
     * Determine the number of residues from the silent file
     */
    input.getline(buf, 400);
    line=buf;
    if (line.substr(0, 9)!="SEQUENCE:")
    {
        cerr << "Silent file \"" << filename << "\" began with:" << endl
             << "    " << line << endl;
        exit(0);
    }

    input.getline(buf, 400);
    line=buf;
    if (line.substr(0, 6)!="SCORE:")
    {
        cerr << "Silent file \"" << filename << "\" began with:" << endl
             << "    " << line << endl;
        exit(0);
    }

    input.getline(buf, 400);
    line=buf;
    if (line.substr(0, 6)!="SCORE:")
    {
        cerr << "Silent file \"" << filename << "\" began with:" << endl
             << "    " << line << endl;
        exit(0);
    }

    int residueID;
    for (int i=1; !input.eof(); i++)
    {
        input.getline(buf, 400);
        line = buf;
        if (line.substr(0, 6) == "SCORE:")
            break;
        residueID = toInt(line.substr(1, 4));
        if (residueID != i)
        {
            cerr << "Residue id " << residueID
                 << " out of sequence in silent file" << endl; 
            exit(0);
        }
    }

    mNumResidue = residueID;
    input.seekg(0);
    input.getline(buf, 400);
    input.getline(buf, 400);
    input.getline(buf, 400);

    /**
     * Read PDBs into filename2PDB
     */
    SimPDB * pdb = new SimPDB(mNumResidue); 
    bool isNewPDB = true;
    int numResidue = 0;
    int decoyCount = 1;
    char * key;
    do
    {
        input.getline(buf, 400);
        line = buf;

        if (line.substr(0, 6) == "SCORE:") // Old PDB done
        {
            // Check if PDB has mNumResidue
            if (numResidue != mNumResidue)
            {
                cerr << "Insufficient residues in the " << decoyCount
                     << "-th decoy in silent file" << endl; 
                exit(0);
            }

            // Insert the pdb
            filename2PDB[key] = pdb;
            pdb->mProteinFileName = key;
            center_residues(pdb->mCAlpha, pdb->mNumResidue);

            // Start a new pdb
            pdb = new SimPDB(mNumResidue);
            isNewPDB = true;
            numResidue = 0;

            decoyCount++;
        }
        else if (line != "") // Sometimes an empty string is read at eof
        {
            if (isNewPDB)
            {
                // Get the filename
                string filename = line.substr(62, 400);
                string ext = filename.substr(filename.length()-4, 4);

                // If filename does not end in .pdb, generate a filename
                if (strcmp(ext.c_str(), ".pdb"))
                    key = strdup("decoy" + decoyCount);
                else
                    key = strdup(filename.c_str());

                isNewPDB = false;
            }

            residueID = toInt(line.substr( 0, 4));
            float x = toFloat(line.substr(35, 8));
            float y = toFloat(line.substr(44, 9));
            float z = toFloat(line.substr(54, 8));
            pdb->mCAlpha[numResidue*3]   = x;
            pdb->mCAlpha[numResidue*3+1] = y;
            pdb->mCAlpha[numResidue*3+2] = z;

            if (residueID != numResidue+1)
            {
                cout << residueID << "," << numResidue << endl;
                cerr << "Residue ID out of sequence in the " << decoyCount
                     << "-th decoy in silent file" << endl; 
                exit(0);
            }

            numResidue++;
        }

    }
    while (!input.eof());
    input.close();

    // Insert the final pdb
    filename2PDB[key] = pdb;
    pdb->mProteinFileName = key;

    mNumDecoy = decoyCount;

    /**
     * Generate the vector of filenames.
     * We can insert filename at the same time as inserting the SimPDB into
     * filename2PDB. However, doing this here has the advantage that if the
     * same filename is entered twice into filename2PDB, the name will not
     * be duplicated here.
     */
    mNames = new vector<char *>(0); 
    map<char *,SimPDB *>::iterator it = filename2PDB.begin();
    for (int i=0; it != filename2PDB.end(); it++, i++)
        mNames->push_back( it->first );
}


/**
 * Populate the PreloadedPDB with the PDB files specified in a list.
 *
 * Reading of the PDB files is through SimPDB.
 */
void
PreloadedPDB::loadPDBFromList(char * filename)
{
    ifstream input(filename);
    if (!input)
    {
        cerr << "Cannot find file \"" << filename << "\"" << endl;
        exit(0);
    }

    char buf[400];
    char* token;

    silentfilename = NULL;
    pdblistfilename = filename;

    /**
     * Read in names of PDB files
     */
    mNames = new vector<char *>(0); 
    while (!input.eof())
    {
        input.getline(buf, 400);
        token = strtok(buf, " ");
        if (token == NULL)
            continue;
        char * name = new char[strlen(token)+1];
        strcpy(name, token);
        mNames->push_back(name);
    }
    input.close();

    mNumDecoy = mNames->size();

    SimPDB * pdb = new SimPDB();
    pdb->mProteinFileName = (*mNames)[0];
    pdb->mNumResidue = LONGEST_CHAIN;
    pdb->mCAlpha = new float[3*LONGEST_CHAIN];
    pdb->read();

    mNumResidue = pdb->mNumResidue;

    filename2PDB[(*mNames)[0]] = pdb;

    for (int i=1; i < mNames->size(); i++)
    {
        SimPDB * pdb = new SimPDB(mNumResidue);
        pdb->mProteinFileName = (*mNames)[i];
        pdb->read();
        filename2PDB[(*mNames)[i]] = pdb;
    }
}


SimPDB *
PreloadedPDB::getSimPDB(char * filename)
{
    return filename2PDB[filename];
}


/*
int main()
{
    //char * file = "rosetta/silent_file";
    char * file = "list";

    // Preload the PDBs, either from silent file or from list of decoys
    PreloadedPDB * pdbs = new PreloadedPDB();
    switch (filetype(file))
    {
        case SILENT_FILE: pdbs->loadSilentFile(file); break;
        case PDB_LIST: pdbs->loadPDBFromList(file); break;
        default: cerr << "Unknown file type" << endl; exit(0);
    }

    // Attach the cache to SimPDB
    SimPDB::preloadedPDB = pdbs;
    SimPDB::preloadPDB = true;

    // Get the names from the loaded PDB for testing
    vector<char *> * mNames = pdbs->mNames;

    // Get the PDB for each name
    for (int i=0; i < mNames->size(); i++)
    {
        char * filename = (*mNames)[i];
        SimPDB * pdb = new SimPDB(filename);
        cout << filename << ":" << endl;
        for (int j=0; j < pdb->mNumResidue; j++)
            cout << "    " << pdb->mCAlpha[j*3] << ", "
                           << pdb->mCAlpha[j*3+1] << ", "
                           << pdb->mCAlpha[j*3+2] << endl;
    }
}
*/

