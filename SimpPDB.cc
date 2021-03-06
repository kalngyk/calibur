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


#include <iostream>
#include <fstream>
#include <sstream>

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>

using namespace std;

#include "SimpPDB.h"


//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -
// Memory and default values for the static variables

// Modifies how SimPDB reads PDB files
int SimPDB::s_residue = 1;
int SimPDB::e_residue = LONGEST_CHAIN;
char * SimPDB::chains = strdup("AC ");
vector<string> SimPDB::atom_names = {"CA"};
vector<string> SimPDB::atom_matchstrs = {" CA ", "CA  ", "  CA", "CA"};


// Set these two fields to tell SimPDB to use the PreloadedPDB mechanism
PreloadedPDB * SimPDB::preloadedPDB = NULL;
bool SimPDB::preloadPDB = true;

//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

char aa[23][4] = {"BCK", "GLY", "ALA", "SER", "CYS", "VAL", "THR", "ILE",
                  "PRO", "MET", "ASP", "ASN", "LEU",
                  "LYS", "GLU", "GLN", "ARG",
                  "HIS", "PHE", "TYR", "TRP", "CYX", "MSE"};
     
char slc[] = {'X','G','A','S','C','V','T','I',
    	      'P','M','D','N','L','K','E','Q','R',
              'H','F','Y','W','C', 'm'};


int toInt(const string& aString)
{
    static char st[20];
    int start = 0;
    for (int i=0; i < aString.size(); i++)
       if (aString[i] != ' ')
           st[start++] = aString[i];
    st[start] = '\0';
    int rev = atoi(st);
    return rev;
}


float toFloat(const string& aString)
{
    static char st[20];
    int start = 0;
    for (int i=0; i < aString.size(); i++)
        if (aString[i] != ' ')
            st[start++] = aString[i];
    st[start] = '\0';
    float rev = atof(st);
    return rev;
}


void center_residues(float * mCAlpha, int mNumResidue)
{
    float cx = 0;
    float cy = 0;
    float cz = 0;
 
    int i3 = 0;
    for (int i=0; i < mNumResidue; i++)
    {
        cx += mCAlpha[i3];
        cy += mCAlpha[i3+1];
        cz += mCAlpha[i3+2];
        i3 += 3;
    }
 
    cx /= mNumResidue;
    cy /= mNumResidue;
    cz /= mNumResidue; 
    i3 = 0;
    for (int i=0; i < mNumResidue; i++)
    {
        mCAlpha[i3]   -= cx; 
        mCAlpha[i3+1] -= cy; 
        mCAlpha[i3+2] -= cz;
        i3 += 3;
    }
}

//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

int
SimPDB::init_atom_names(string namelist, char delimiter)
{
    istringstream is(namelist);
    string name;
    atom_names.clear();
    atom_matchstrs.clear();
    while (getline(is, name, delimiter))
    {
        atom_names.push_back(name);
        switch (name.length())
        {
            case 1:
            {
                string a[4] = {"   " + name, name + "   ",
                               "  " + name + " ", " " + name + "  "};
                atom_matchstrs.insert(atom_matchstrs.end(), a, a+4);
                break;
            }
            case 2:
            {
                string a[3] = {"  " + name, name + "  ", " " + name + " "};
                atom_matchstrs.insert(atom_matchstrs.end(), a, a+3);
                break;
            }
            case 3:
            {
                string a[2] = {" " + name, name + " "};
                atom_matchstrs.insert(atom_matchstrs.end(), a, a+2);
                break;
            }
            default: // 0 or >3
            {
                cout << "Invalid atom name specification \"" << name << "\": "
                     << "atom name must be 1~3 characters long"
                     << endl << endl;
                return -1;
            }
        }
    }
    return 0;
}

//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -
// General purpose constructors which handle the preloadedPDB mechanism

SimPDB::SimPDB(char *aFileName)
{
    if (preloadPDB)
    {
        /*for(int i=0; i<len; i++)
        {
          mCAlpha[i]=new float [3]; 
        }*/
        SimPDB * pdb = preloadedPDB->filename2PDB[aFileName];
        mProteinFileName = strdup(pdb->mProteinFileName);
        mNumResidue = pdb->mNumResidue;
        mCAlpha = new float[3*mNumResidue];
        memcpy(mCAlpha, pdb->mCAlpha, 3 * mNumResidue * sizeof(float));
    }
    else
    {
        mProteinFileName = aFileName;
        mNumResidue = LONGEST_CHAIN;
        mCAlpha = new float[3*LONGEST_CHAIN];
        read();
    }
}

SimPDB::SimPDB(char *aFileName, int len)
{
    if (preloadPDB)
    {
        /*for(int i=0; i<len; i++)
        {
          mCAlpha[i]=new float [3]; 
        }*/
        SimPDB * pdb = preloadedPDB->filename2PDB[aFileName];
        mProteinFileName = strdup(pdb->mProteinFileName);
        mNumResidue = pdb->mNumResidue;
        mCAlpha = new float[3*mNumResidue];
        memcpy(mCAlpha, pdb->mCAlpha, 3 * mNumResidue * sizeof(float));
    }
    else
    {
        mProteinFileName = aFileName;
        mNumResidue = len;
        mCAlpha = new float[3*len];
        int count = read();
        if (count != mNumResidue)
        {
            cout << "Error: \"" << aFileName
                 << "\" has mismatching number of residues"
                 << " (should have " << mNumResidue
                 << " but has only " << count << ")" << endl;
            exit(0);
        }
    }
}


SimPDB::~SimPDB() 
{
    /* for(int i=0; i<mNumResidue; i++)
    {
      delete [] mCAlpha[i];
    }*/
    delete [] mCAlpha;
}


//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -
// Constructors which DO NOT handle the preloadedPDB mechanism
//
// They should be called from PreloadedPDB only, since it will need to
// bypass the mechanism

SimPDB::SimPDB() {}

SimPDB::SimPDB(int len)
{
    mNumResidue = len;
    mCAlpha = new float[3*len];
}



//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

static bool
_matches_atom_names(string line)
{
    for (int i=0; i < SimPDB::atom_matchstrs.size(); i++)
        if (line.substr(13, SimPDB::atom_matchstrs[i].size()) ==
               SimPDB::atom_matchstrs[i])
            return true;
    return false;
}

/**
 * Reads a PDB file from disk. Do not read from PDB files anywhere else.
 */
int
SimPDB::read()
{
    ifstream input(mProteinFileName);
    if (!input)
    {
        cerr << "Cannot find protein file " << mProteinFileName << endl;
        exit(0);
    }
    cout.flush();
    char buf[400];
    mNumResidue = 0;
    float x, y, z;
    char c = 'a';
    bool read = false;
    
    int prevID = -10000;
    int count = 0;
    
    int CA_number = 1;

    //float cx = 0;
    //float cy = 0;
    //float cz = 0;
    int count3 = 0; 
    //mSquaredSum=0;   
    for (int rcount=0; !input.eof(); rcount++)
    {
        if (rcount > 10000 && count == 0) // this ain't no PDB file bruh
            break;
        input.getline(buf, 400);
        string line=buf;
        if (line.substr(0, 3) == "TER" && read == true) break;
        if (line.substr(0, 6) == "ENDMDL") break;

        if (line.substr(0, 4) != "ATOM" && line.substr(0, 6) != "HETATM")
            continue;

        if (_matches_atom_names(line))
        {
            // At this point an atom in PDB::atom_names has been discovered
            // We want to further filter it based on the following two
            // criteria: chain and region.
           
            // Check if the chain which this atom belongs to is to be included
            bool include_chain = false;
            for (char * c = SimPDB::chains; *c; c++)
                if (*c == '*' || toupper(line[21]) == toupper(*c))
                {
                    include_chain = true;
                    break;
                }
            if (!include_chain)
            {
                CA_number++;
                continue;
            }

            // Check if the CA atom is within the region to analyze
            if (CA_number < SimPDB::s_residue)
            {
                CA_number++;
                continue;
            }
            else if (CA_number > SimPDB::e_residue)
                break;

            read = true;
            int residueID = toInt(line.substr(22, 6));
            if (SimPDB::atom_names.size() == 1 && residueID == prevID)
                continue;

            prevID = residueID; 
            x = toFloat(line.substr(30, 8));
            y = toFloat(line.substr(38, 8));
            z = toFloat(line.substr(46, 8));
            string AAType = line.substr();
            count3 = 3*count;   
            mCAlpha[count3]   = x; 
            mCAlpha[count3+1] = y; 
            mCAlpha[count3+2] = z;
            //mSquaredSum+=x*x+y*y+z*z;
            count++;

            CA_number++;
        }
    }//while

    mNumResidue = count;
    input.close();

    center_residues(mCAlpha, mNumResidue);
    return count;
}

/*
int main(int argc, char** argv)
{
    SimPDB::chains = strdup("ACH");
    SimPDB::s_residue = 215;
    SimPDB::e_residue = 224;
    SimPDB* sim = new SimPDB(argv[1]);
    
    cout << "numResidue=" << sim->mNumResidue << endl;
    for (int i=0; i < sim->mNumResidue; i++)
    {
        cout << sim->mCAlpha[i*3] << ","
             << sim->mCAlpha[i*3+1] << ","
             << sim->mCAlpha[i*3+2] << endl;
    }
}*/
