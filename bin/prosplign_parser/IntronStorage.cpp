#include <float.h>
#include <iostream>

#include "IntronStorage.h"

using namespace std;

IntronStorage::IntronStorage() {
    clear();
}

void IntronStorage::clear() {
    minScore = DBL_MAX;
    maxScore = 0;
    introns.clear();
}

void IntronStorage::storeIntron(string protein, string gene, int start, int end,
        char strand, string spliceSites, double score, int number) {

    Intron i(protein, gene, start, end, strand, spliceSites, score, number);
    introns.push_back(i);

    if (score > maxScore) {
        maxScore = score;
    }
    if (score < minScore) {
        minScore = score;
    }

}

void IntronStorage::normalizeScores() {
    normalizeScores(minScore, maxScore);
}

void IntronStorage::normalizeScores(double min, double max) {
    double range = max - min;
    if (range == 0) {
        return;
    }
    for (unsigned int i = 0; i < introns.size(); i++) {
        if (introns[i].score > max + 0.00001) {
            cerr << "Warning, supplied maximum score is lower than a score of "
                    "intron in " << introns[i].gene << "-" << introns[i].protein << endl;
            introns[i].score = max;
        }
        if (introns[i].score < min - 0.00001) {
            cerr << "Warning, supplied minimum score is higher than a score of "
                    "intron in " << introns[i].gene << "-" << introns[i].protein << endl;
            introns[i].score = min;
        }
        introns[i].score = (introns[i].score - min) / range;
    }
}

void IntronStorage::printIntrons(string output, bool printAll) {
    ofstream ofs(output.c_str());
    for (unsigned int i = 0; i < introns.size(); i++) {
        // Discard non-canonical splice sites unless the splice site is
        // GC-AG or AT-AC with alignment score better than PRINT_THRESHOLD
        string spliceSites = introns[i].spliceSites;
        if (!printAll && spliceSites != "GT_AG") {
            if (spliceSites == "GC_AG" || spliceSites == "AT_AC") {
                if (introns[i].score < PRINT_THRESHOLD) {
                    continue;
                }
            } else {
                continue;
            }
        }

        ofs << introns[i].gene << "\tProSplign\tIntron\t";
        ofs << introns[i].start << "\t";
        ofs << introns[i].end << "\t";
        ofs << ".\t+\t.\tprot=" << introns[i].protein;
        ofs << "; intron_id=" << introns[i].number << ";";
        ofs << " splice_sites=" << spliceSites << ";";
        ofs << " al_score=" << introns[i].score << ";" << endl;
    }
    ofs.close();
}
