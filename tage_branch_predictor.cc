#include "tage_branch_predictor.h"

int Seed = 0;

TaggedComponentEntry::TaggedComponentEntry()
{
    ctr = 0;
    tag = 0;
    uCtr = 0;
}

/* Signed up/down counter */
    void
TaggedComponentEntry::ctrUpdate(bool actual, int nbits)
{
    if (actual) {
        if (ctr < ((1 << (nbits - 1)) - 1)) {
            ctr++;
        }
    }
    else {
        if (ctr > -(1 << (nbits - 1))) {
            ctr--;
        }
    }
}

    void
TaggedComponentEntry::uCtrUpdate(bool actual, bool providerPred, bool altPred)
{
    if ((providerPred != altPred)) {
        if (providerPred == actual) {
            if (uCtr < 3)
                uCtr++;
        } else {
            if (uCtr > 0)
                uCtr--;
        }
    }
}

TaggedComponent::TaggedComponent() {}

    UInt16
TaggedComponent::matchTagAtIndex(int index, int tag)
{
    return (entry[index].tag == tag);
}

    bool
TaggedComponent::getPredAtIndex(int index)
{
    return (entry[index].ctr >= 0);
}

    void
TaggedComponent::uCtrPeriodicReset(int brCtr)
{
    if ((brCtr & ((1 << 18) - 1)) == 0) {
        for (int i = 0; i < (1 << LOG_ENTRIES_COMP); i++)
            entry[i].uCtr = entry[i].uCtr >> 1;
    }
    return;
}

TageBranchPredictor::TageBranchPredictor(String name, core_id_t core_id):BranchPredictor(name, core_id)
{
    gHist = 0;
    pHist = 0;
    brCtr = 0;

    /* Computation of geometric history lengths */
    histLength[1] = MIN_HIST_LEN;
    histLength[NUM_COMPONENTS] = MAX_HIST_LEN;
    for (int i=2; i <= NUM_COMPONENTS; i++) {
        histLength[i] = (int) (((double) MIN_HIST_LEN *
                    std::pow ((double) (MAX_HIST_LEN) / (double) MIN_HIST_LEN,
                        (double) (i-1) / (double) ((NUM_COMPONENTS - 1)))) + 0.5);
    }

    /* Computation of partial tag widths */
    for (int i=1; i <= NUM_COMPONENTS; i++) {
        tagWidth[NUM_COMPONENTS-i+1] = TAGWIDTH_THRESHOLD - ((i + (NUM_COMPONENTS & 1)) / 2);
    }

    /*for (int i=1; i <= NUM_COMPONENTS; i++) {
      std::cout << __FUNCTION__ << " " << i << ", " << tagWidth[i] << ", " << histLength[i]  << "\n";
      }*/

    /* Initialization of index computation function */
    for (int i=1; i <= NUM_COMPONENTS; i++) {
        histIndex[i].init (histLength[i], (LOG_ENTRIES_COMP));
    }

    /* Initialization of tag computation function */
    for (int i=1; i <= NUM_COMPONENTS; i++) {
        histTag[0][i].init (histIndex[i].OLENGTH, tagWidth[i]);
        histTag[1][i].init (histIndex[i].OLENGTH, tagWidth[i]-1);
    }
}

    void
TageBranchPredictor::updateHistory(bool actual, IntPtr ip)
{
    gHist = (gHist << 1);
    if (actual)
        gHist |= (history_t) 1;

    pHist = (pHist << 1) + (ip & 1);
    pHist = (pHist & ((1 << 16) - 1));
    for (int i=1; i <= NUM_COMPONENTS; i++)
    {
        histIndex[i].update(gHist);
        histTag[0][i].update(gHist);
        histTag[1][i].update(gHist);
    }
    return;
}

/* 
 * Functions F(), gIndex() and gTag()
 * from https://github.com/masc-ucsc/esesc/blob/master/simu/libcore/IMLIBest.h#L1041-L1072
 */
    int
F(int A, int size, int bank)
{
    int A1, A2;

    A = A & ((1 << size) - 1);
    A1 = (A & ((1 << LOG_ENTRIES_COMP) - 1));
    A2 = (A >> LOG_ENTRIES_COMP);
    A2 = ((A2 << bank) & ((1 << LOG_ENTRIES_COMP) - 1)) + (A2 >> (LOG_ENTRIES_COMP - bank));
    A = A1 ^ A2;
    A = ((A << bank) & ((1 << LOG_ENTRIES_COMP) - 1)) + (A >> (LOG_ENTRIES_COMP - bank));
    return (A);
}

    int
TageBranchPredictor::gIndex(IntPtr ip, int bank)
{
    int index;
    int len = (histLength[bank] >= 16) ? 16 : histLength[bank];
    index =
        ip ^ (ip >> ((LOG_ENTRIES_COMP - (NUM_COMPONENTS - bank - 1)))) ^ histIndex[bank].
        comp ^ F(pHist, len, bank);

    return (index & ((1 << (LOG_ENTRIES_COMP)) - 1));
}

    UInt16
TageBranchPredictor::gTag(IntPtr ip, int bank)
{
    int tag = ip ^ histTag[0][bank].comp ^ (histTag[1][bank].comp << 1);
    return (tag & ((1 << tagWidth[bank]) - 1));
}

    bool
TageBranchPredictor::predict(IntPtr ip, IntPtr target)
{
    for (int i=1; i <= NUM_COMPONENTS; i++) {
        GI[i] = gIndex(ip, i);
        GTAG[i] = gTag(ip, i);
    }

    providerComp = 0;
    altComp = 0;

    /* Find hitting component with longest history */
    for (int i = NUM_COMPONENTS; i >= 1; i--) {
        if (T[i].matchTagAtIndex(GI[i], GTAG[i])) {
            providerComp = i;
            break;
        }
    }

    /* Alternate prediction component */
    for (int i = providerComp-1; i >= 1; i--) {
        if (T[i].matchTagAtIndex(GI[i], GTAG[i])) {
            if ((altOnNewAlloc < 0) || (abs (2 * T[i].entry[GI[i]].ctr + 1) > 1)) {
                altComp = i;
                break;
            }
        }
    }

    if (providerComp > 0) {
        /* If there is a hitting component */
        if (altComp > 0) {
            altPred = T[altComp].getPredAtIndex(GI[altComp]);
        } else {
            /* Use base prediction to use as alternate prediction */
            altPred = basePredictor.predict(ip, target);
        }
        if ((altOnNewAlloc < 0) || (abs (2 * T[providerComp].entry[GI[providerComp]].ctr + 1) > 1)
                || (T[providerComp].entry[GI[providerComp]].uCtr != 0)) {
            providerPred = T[providerComp].getPredAtIndex(GI[providerComp]);
        } else {
            providerPred = altPred;
        }
    } else {
        /* Only base prediction is considered */
        altPred = basePredictor.predict(ip, target);
        providerPred = altPred;
    }
    return providerPred;
}

    int
getRandom()
{
    Seed = ((1 << 2 * NUM_COMPONENTS) + 1) * Seed + 0xf3f531;
    Seed = (Seed & ((1 << (2 * (NUM_COMPONENTS))) - 1));
    return (Seed);
};

/*
 * Function allocate():
 * 1. Find alloc
 *      1.1. True if prediction was wrong
 *      1.2. If it was delivering the correct prediction, no need to allocate 
 *           a new entry even if the overall prediction was wrong, alloc = false
 *      1.3. Update altOnNewAlloc
 *
 * 2. If alloc is true, do allocation
 *  
 */
    void
TageBranchPredictor::allocate(bool actual)
{
    int NRAND = getRandom();
    bool alloc = ((providerPred != actual) & (providerComp < NUM_COMPONENTS));
    if (providerComp > 0)
    {
        bool longestPred = T[providerComp].getPredAtIndex(GI[providerComp]);
        bool pseudoNewAlloc = (abs (2 * T[providerComp].entry[GI[providerComp]].ctr + 1) <= 1);
        if (pseudoNewAlloc)
        {
            if (longestPred == actual)
                alloc = false;
            if (longestPred != altPred)
            {
                if (altPred == actual)
                {
                    if (altOnNewAlloc < 7)
                        altOnNewAlloc++;
                }
                else if (altOnNewAlloc > -8)
                    altOnNewAlloc--;
            }
        }
    }
    if (alloc)
    {
        SInt8 min = 1;
        for (int i = NUM_COMPONENTS; i > providerComp; i--) {
            if (T[i].entry[GI[i]].uCtr < min)
                min = T[i].entry[GI[i]].uCtr;
        }
        int Y = NRAND & ((1 << (NUM_COMPONENTS - providerComp - 1)) - 1);
        int X = providerComp + 1;
        if (Y & 1) {
            X++;
            if (Y & 2) {
                X++;
            }
        }

        if (min > 0) {
            T[X].entry[GI[X]].uCtr = 0;
        }
        for (int i = X; i <= NUM_COMPONENTS; i++)
        {
            if ((T[i].entry[GI[i]].uCtr == min))
            {
                T[i].entry[GI[i]].tag = GTAG[i];
                T[i].entry[GI[i]].ctr = (actual) ? 0 : -1;
                T[i].entry[GI[i]].uCtr = 0;
                break;
            }
        }
    }
}

    void
TageBranchPredictor::update(bool predicted, bool actual, IntPtr ip, IntPtr target)
{
    //std::cout << __FUNCTION__ << "::Predicted: " << predicted << ", Actual: " << actual << "\n";
    brCtr++;
    updateCounters(predicted, actual);
    allocate(actual);
    if (providerComp > 0) {
        T[providerComp].entry[GI[providerComp]].ctrUpdate(actual, NUM_BITS_CTR);
        if (T[providerComp].entry[GI[providerComp]].uCtr == 0) {
            if (altComp > 0)
                T[altComp].entry[GI[altComp]].ctrUpdate(actual, NUM_BITS_CTR);
            if (altComp == 0)
                basePredictor.update(predicted, actual, ip, target);
        }
    } else {
        basePredictor.update(predicted, actual, ip, target);
    }

    T[providerComp].entry[GI[providerComp]].uCtrUpdate(actual, providerPred, altPred);

    updateHistory(actual, ip);
    for (int i = 1; i <= NUM_COMPONENTS; i++) {
        T[i].uCtrPeriodicReset(brCtr);
    }
    return;
}
