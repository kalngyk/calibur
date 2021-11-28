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


#ifndef _INIT_CLUSTER_
#define _INIT_CLUSTER_

#include "SimpPDB.h"
//#include "sys/resource.h"
#include <vector>
#include <stdlib.h>
#include <string.h>

using namespace std;

#ifndef _LARGE_DECOY_SET_
typedef unsigned short LIST_TYPE;
#else  //number of input decoys large than 65535
typedef unsigned int LIST_TYPE;
#endif

// _MATRIX_MODE_LIMIT_ is determined by the RAM of the system
#define _MATRIX_MODE_LIMIT_ 13000

#define _OVER_RMSD_ 9999999

#define REFERENCE_SIZE 6
#define RANDOM_DECOY_SIZE_FOR_FILTERING 101
#define RANDOM_DECOY_SIZE_FOR_THRESHOLD 101
#define NUM_TRIALS_FOR_THRESHOLD 16
#define DEFAULT_PERCENTILE_FOR_THRESHOLD 10
#define MAX_PERCENTILE_FOR_THRESHOLD 50
#define MIN_PERCENTILE_FOR_THRESHOLD 3


enum ADJ_LIST_MODE { MATRIX, LIST , LITE };

enum EST_THRESHOLD_MODE {
    PERCENT_EDGES,      // % quantile pairwise distance (default)
    MIN_AVG_DIST_BASED, // t = a * (min avg dist) + b
    MOST_FREQ_BASED,    // use (2/3)*(most frequently occuring distance)
    ROSETTA,            // ROSETTA's method
    SAMPLED_ROSETTA,    // ROSETTA's method using sampled decoys
    USER_SPECIFIED,      // user supplied
};



class AdjacentList
{
public:
    static ADJ_LIST_MODE mListMode;
    int mWhich;         // index of the decoy this AdjacentList is for
    int mNumNeigh;      // synchronized with the size of neigh
    vector<int>* neigh; // keep a record of all the neighbors
    vector<float>* dist;
    LIST_TYPE* mReverseIndex; // References index of the decoy within the array
    AdjacentList();
    AdjacentList(int which, int num);
    ~AdjacentList();
    void add(int n, float d, bool ifNeigh);
    void add(int num);
    void remove(int n);
    float getD(int n);
};


class Stru
{
public:
    float* mCAlpha;
    float* mSIG; //signature
    SimPDB* mPDB;
    Stru(SimPDB* pdb, int len);
    ~Stru();
    float dist(float x, float y, float z, float *zz);
    float dist(float x, float y, float z);
};




class Clustering
{
public:
    static EST_THRESHOLD_MODE EST_THRESHOLD;
    static bool FILTER_MODE;
    static bool OUTPUT_ALL;
    static float xPercentile;
    static bool autoAdjustPercentile;
    static float xFactor;

    char* mInputFileName;   // file which contains all PDB filenames
    vector<char* >* mNames; // all decoy (file) names
    vector<Stru* >* mPDBs;  // all decoy PDBs
    int mNumPDB;            // will be set to mPDBs->size()
    int mLen;               // #residues

    float THRESHOLD;        // clustering threshold. most important parameter

    // - = - = - = - = - = - = - = - = - = - = - = - = - = -
    // for auxiliary grouping

    float CLU_RADIUS;     // cluster radius for auxClustering()
    vector<int>* mCluCen; // cluster centers found using auxClustering()
    vector<int>** mAuxCluster; // cluster elements
    int* mCen;
    float* mD2C;          // distance from decoy in auxCluster to CluCen
    int * mNumNeighbor;   // the number of neighbors of each decoy
    float bestClusMargin; // size(bestClus) -origsuze(2ndClus) /size(bestClus)
    int bestClusSize;

    // - = - = - = - = - = - = - = - = - = - = - = - = - = -
    // for clustering

    AdjacentList** mAdjacentList; // lists of all neighbors
    float* mReference;      // for {lower,upper}bounds through references

    int mFinalDecoy;
    vector<AdjacentList *> *mFinalClusters;
    int * mRemainingList;
    int * mRemainingListIndex;
    int mRemainingSize;

    //- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -
    // METHODS
    //- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

    Clustering(); // do nothing
    void initialize(char * filename, float threshold);
    void reinitialize(vector<char *>*, vector<Stru *>*, float threshold);
    void cluster();
    void showClusters(int);
    void getPDBs(vector<char *>*, vector<Stru *>*, vector<int>*, int);


    //- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

    // for reading decoys from input files
    void readDecoyNames();
    void readDecoys(vector<char *>*, vector<Stru *>*);
    vector<Stru *>* readDecoys(vector<char *>*);
    //void refilterDecoys(vector<char *>*, vector<Stru *>*);

    // - = - = - = - = - = - = - = - = - = - = - = -
    // for finding threshold

    void getThresholdAndDecoys();
    float getThreshold(float **, int, int, int, int, float, float);
    float ** getNborList(vector<Stru *> *, vector<char *> *, int);
    //float ** getNborList(vector<Stru *> *, vector<Stru *> *, int);
    float ** getNborList(vector<Stru *> *, int, float *, float *);
    void estimateDist(vector<char *>*, int, int, float,
                      float *, float *, float *, float *);
    void estimateDist(vector<Stru *>*, float,
                      float *, float *, float *, float *);

    vector<char *>* getRandomDecoyNames(vector<char *>*, int, int);
    void destroyRandomDecoys(vector<char *>*, vector<Stru *>*); 

 
    // - = - = - = - = - = - = - = - = - = - = - = -
    // for clustering

    void auxClustering();
    void buildAdjacentLists();
    void listAdjacentLists();

    void findLargestClusters();
    int findDecoyWithMostNeighbors();
    void removeDecoys(AdjacentList * adj);

    // - = - = - = - = - = - = - = - = - = - = - = -

    void initRef(int* index);
    void refBound(int i, int j, float& lower, float& upper);
    //bool find(int which, vector<int> *elements); // too slow

    // - = - = - = - = - = - = - = - = - = - = - = -

    float realignDecoys(int ref);
    float superimposeAndReplace(float* coor1, float* coor2);
    float eucD(int i, int j);

    // - = - = - = - = - = - = - = - = - = - = - = -

    // methods for rmsd computation
    float estD(int i, int j);
    float estD(Stru* a, Stru* b);
    float trueD(Stru* a, Stru* b);
    float trueD(int i, int j);
    // storage for RMSD() computation
    void allocateSpaceForRMSD(int len);
    bool spaceAllocatedForRMSD; 
    double *result_coords;
    // storage for rmsfit_() computation
    double *coord1;
    double *coord2;

    // - = - = - = - = - = - = - = - = - = - = - = -

    //double timeval_difference(struct timeval * x, struct timeval * y);
    double get_elapsed(int restart);
};
#endif

