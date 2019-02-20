#ifndef TAGE_BRANCH_PREDICTOR_H
#define TAGE_BRANCH_PREDICTOR_H

#include <cstdlib>
#include <cmath>
#include <bitset>
#include "branch_predictor.h"
#include "simulator.h"
#include "tage_base_predictor.h"

#define NUM_BITS_CTR 3
#define NUM_COMPONENTS 5
#define LOG_ENTRIES_COMP 9
#define TAGWIDTH_THRESHOLD 9
#define MAX_HIST_LEN 131 
#define MIN_HIST_LEN 5

typedef std::bitset <MAX_HIST_LEN> history_t;

/* 
 * class folded_history: Used for index and tag computation.
 * from https://github.com/masc-ucsc/esesc/blob/master/simu/libcore/IMLIBest.h#L228-L261
 */
class folded_history
{
    public:
        unsigned comp;
        int CLENGTH;
        int OLENGTH;
        int OUTPOINT;

        folded_history() {}

        void init (int original_length, int compressed_length)
        {
            comp = 0;
            OLENGTH = original_length;
            CLENGTH = compressed_length;
            //std::cout << __FUNCTION__ << "::"  << OLENGTH << ", " << CLENGTH << "\n";
            OUTPOINT = OLENGTH % CLENGTH;
        }

        void update (history_t h)
        {
            comp = (comp << 1) | h[0];
            comp ^= h[OLENGTH] << OUTPOINT;
            comp ^= (comp >> CLENGTH);
            comp &= (1 << CLENGTH) - 1;
        }
};

class TaggedComponentEntry
{
    public:
        SInt8 ctr;     // prediction depends on the sign
        UInt16 tag;   // partial tag
        SInt8 uCtr;    // useful counter
        TaggedComponentEntry();
        void ctrUpdate(bool, int);
        void uCtrUpdate(bool, bool, bool);
};

class TaggedComponent
{
    public:
        TaggedComponentEntry entry[1 << (LOG_ENTRIES_COMP)];

        TaggedComponent();

        void uCtrPeriodicReset(int);
        UInt16 matchTagAtIndex(int, int);
        bool getPredAtIndex(int index);
};

class TageBranchPredictor : public BranchPredictor
{
    public:
        int altOnNewAlloc;
        int brCtr;
        int pHist;
        history_t gHist;

        folded_history histIndex[NUM_COMPONENTS+1];
        folded_history histTag[2][NUM_COMPONENTS+1];

        TaggedComponent T[NUM_COMPONENTS+1];
        TageBasePredictor basePredictor;

        int histLength[NUM_COMPONENTS+1];
        int tagWidth[NUM_COMPONENTS+1];
        int providerComp, altComp;
        bool providerPred, altPred;
        int GI[NUM_COMPONENTS+1];
        int GTAG[NUM_COMPONENTS+1];

        TageBranchPredictor(String, core_id_t);

        int gIndex(IntPtr, int);
        UInt16 gTag(IntPtr, int);
        void uCtrPeriodicReset(void);
        void updateHistory(bool, IntPtr);
        void allocate(bool);

        bool predict(IntPtr, IntPtr);
        void update(bool, bool, IntPtr, IntPtr);
};
#endif // TAGE_BRANCH_PREDICTOR_H
