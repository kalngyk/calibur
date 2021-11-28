/*
 *  **************************************************************************
 *  Copyright 2009, 2021 Shuai Cheng Li and Yen Kaow Ng
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


#include <math.h>
#include <stdlib.h>
#include "InitCluster.h"
#include "SimpPDB.h"


void
usage(char * progname)
{
  cerr << "Usage: " << progname
  << " [-n] [-o] [-r #1,#2] [-c XYZ] [-t s] pdb_list [x]" << endl << endl
  << "  pdb_list is a text file which specifies the decoys. Each line in"
  << " pdb_list is" << endl
  << "    a path (relative to the working directory) to a decoy's PDB file."
  << endl << endl
  << "  -n (optional) disables the filtering of outlier decoys."
  << endl << endl
  << "  -o (optional) output all clusters instead of only the top three."
  << endl << endl
  << "  -r (optional) limits residues to only the #1-th till #2-th C-alpha atoms."
  << endl << endl
  << "  -c (optional) specifies that the chains XYZ are to be used."
  << endl
  << "                By default, XYZ=\"AC \", i.e. 'A', 'C', or unspecified"
  << endl << endl
  << "  -t (optional) specifies the threshold finding strategy." << endl
  << "    s is one of p, f, a, R, r. (default strategy: p)" << endl
  << "     p: threshold results in only x\% of \"edges\" between decoys."
  << endl
  << "         default x=100/sqrt(sqrt(#decoys))" << endl
  << "     f: threshold = min dist + x * (most frequent dist - min dist)"
  << endl
  << "         default x=0.666667 (=2/3)" << endl
  << "     a: threshold = min dist + x * min(the avarage dist of decoys from a"
  << " decoy)" << endl
  << "         default x=0.666667 (=2/3)" << endl
  << "     R: find threshold using ROSETTA's method (auto-detect parameters)."
  << endl
  << "     r: same as R, but with a sampled decoy set rather than the full set."
  << endl << endl
  << "  x (optional) specifies a floating point number" << endl
  << "    x is used according to the threshold strategy specified."
  << " (x is ignored"
  << endl
  << "         if the strategy does not use it.)" << endl
  << "    x is used as the threshold if no threshold strategy is specified."
  << endl;
}

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        usage(argv[0]);
        exit(0);
    }

    Clustering * ic = new Clustering();
    bool strategy_specified = false;
    int i;

    //SimPDB::e_residue = LONGEST_CHAIN;
    //SimPDB::s_residue = 1;
    //SimPDB::chains = strdup("AC ");

    for (i = 1; i < argc; i++)
    {
        if ('-' != *argv[i])
            break;
        switch (argv[i][1])
        {
            case 'd': // FIXME feature not revealed in Usage yet
                SimPDB::preloadPDB = false;
                break;
            case 'c':
                i++;
                if (i == argc)
                {
                    usage(argv[0]);
                    exit(0);
                }
                SimPDB::chains = strdup(argv[i]);
                break;
            case 'r':
                i++;
                if (i == argc)
                {
                    usage(argv[0]);
                    exit(0);
                }
                {
                char * segment_spec = strdup(argv[i]);
                char * p = segment_spec;
                for (; *p && *p != ','; p++)
                    ;
                if (*p == '\0') // then no ',' found. spec is "start"
                {
                    SimPDB::s_residue = atoi(segment_spec);
                    SimPDB::e_residue = LONGEST_CHAIN;
                    if (SimPDB::s_residue <= 0) // invalid spec
                    {
                        usage(argv[0]);
                        exit(0);
                    }
                }
                else if (p == segment_spec) // spec is ",end", or ","
                {
                    SimPDB::e_residue = atoi(p+1);
                    SimPDB::s_residue = 0;
                    if (SimPDB::e_residue <= 0)
                    {
                        usage(argv[0]);
                        exit(0);
                    }
                }
                else // spec must be "start," or "start,end"
                {
                    SimPDB::e_residue = (*(p+1)=='\0')?LONGEST_CHAIN:atoi(p+1);
                    *p = '\0';
                    SimPDB::s_residue = atoi(segment_spec);
                    if (SimPDB::e_residue <= 0 || SimPDB::s_residue <= 0
                        || SimPDB::e_residue < SimPDB::s_residue)
                    {
                        usage(argv[0]);
                        exit(0);
                    }
                }
                }
                cout << "Using C-alphas #" << SimPDB::s_residue << "-";
                if (SimPDB::e_residue == LONGEST_CHAIN)
                    cout << "end" << endl;
                else
                    cout << "#" << SimPDB::e_residue << endl;
                break;
            case 'n':
                Clustering::FILTER_MODE = false;
                break;
            case 'o':
                Clustering::OUTPUT_ALL = true;
                break;
            case 't':
                i++;
                if (i == argc || argv[i][1] != '\0')
                {
                    usage(argv[0]);
                    exit(0);
                }
                strategy_specified = true;
                switch (argv[i][0])
                {
                    case 'p':
                        Clustering::EST_THRESHOLD = PERCENT_EDGES;
                        break;
                    case 'f':
                        Clustering::EST_THRESHOLD = MOST_FREQ_BASED;
                        break;
                    case 'a':
                        Clustering::EST_THRESHOLD = MIN_AVG_DIST_BASED;
                        break;
                    case 'R':
                        Clustering::EST_THRESHOLD = ROSETTA;
                        break;
                    case 'r':
                        Clustering::EST_THRESHOLD = SAMPLED_ROSETTA;
                        break;
                    default:
                        usage(argv[0]);
                        exit(0);
                }
                break;
            default:
                usage(argv[0]);
                exit(0);
        }
    }

    if (i == argc)
    {
        usage(argv[0]);
        exit(0);
    }

    char * filename = strdup(argv[i]);

    float threshold = -1;
    i++;
    if (i == argc-1)
    {
        float c = atof(argv[i]);
        if (!strategy_specified)
            threshold = c;
        else
            switch (Clustering::EST_THRESHOLD)
            {
                case MOST_FREQ_BASED:
                case MIN_AVG_DIST_BASED:
                    Clustering::xFactor = c;
                    break;
                case PERCENT_EDGES:
                    Clustering::autoAdjustPercentile = false;
                    Clustering::xPercentile = c;
                    break;
                default:
                    break; // ignore
            }
    }

    /*
    cout << "filter mode=" << Clustering::FILTER_MODE << endl;
    cout << "strategy=" << Clustering::EST_THRESHOLD << endl;
    cout << "filename=" << filename << endl;
    cout << "xfactor=" << Clustering::xFactor << endl;
    cout << "xpercent=" << Clustering::xPercentile << endl;
    cout << "threshold=" << threshold << endl;
    exit(0);
    */

    ic->initialize(filename, threshold);
    ic->cluster();

    float acceptMargin = 0.15;
    if (ic->bestClusMargin < acceptMargin)
    {
        cout << "Best cluster larger than 2nd best cluster by only "
             << (ic->bestClusMargin*100) << "% (<"
             << (acceptMargin*100) << "%)" << endl
             << "Two possible clusters could be present." << endl
             << "Starting refined clustering..." << endl;

        // create new PDBs and Names out of the elements in the best two
        // clusters

        // first get the lists
        vector<AdjacentList *> * finalClusters = ic->mFinalClusters;
        vector<char *>* Names = new vector<char *>(0);
        vector<Stru *>* PDBs = new vector<Stru *>(0);

        // then add elements into them
        AdjacentList* clus;
        clus = (*finalClusters)[1];
        ic->getPDBs(Names, PDBs, clus->neigh, clus->mNumNeigh);
        clus = (*finalClusters)[0];
        ic->getPDBs(Names, PDBs, clus->neigh, clus->mNumNeigh);

        // Refined Clustering
        float minDist, maxDist, mostFreqDist, xPercentileDist;
        int numDecoys = Names->size() > 2*RANDOM_DECOY_SIZE_FOR_THRESHOLD?
                         RANDOM_DECOY_SIZE_FOR_THRESHOLD: Names->size()/2;
        ic->estimateDist(Names,
                         NUM_TRIALS_FOR_THRESHOLD,
                         numDecoys,
                         0.5,
                         &minDist,
                         &maxDist,
                         &mostFreqDist,
                         &xPercentileDist);
        ic->reinitialize(Names, PDBs, xPercentileDist);
        ic->cluster();

        if (ic->bestClusMargin < acceptMargin)
            cout << "MORE THAN ONE BEST DECOYS DETECTED!" << endl;
    }

    if (Clustering::OUTPUT_ALL)
    {
        cout << "Showing all clusters:" << endl;
        ic->showClusters(ic->mNumPDB);
    }
    else
    {
        cout << "Showing at most three largest clusters:" << endl;
        ic->showClusters(3);
    }
}

