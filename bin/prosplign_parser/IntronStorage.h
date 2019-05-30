#ifndef INTRON_STORAGE_H
#define INTRON_STORAGE_H

#include <string>
#include <fstream>
#include <vector>

#define PRINT_THRESHOLD 0.7

using namespace std;

/**
 * Class for storing and working with the set of all introns found in multiple
 * alignments
 */

class IntronStorage {
public:
    IntronStorage();
    /**
     * Delete all introns in the storage
     */
    void clear();
    /**
     * Add an intron to the storage
     */
    void storeIntron(string protein, string gene, int start, int end,
            char strand, string spliceSites, double score, int number);
    /**
     * Normalize intron scores into 0-1 range, use maximum and minimum
     * directly from the set.
     */
    void normalizeScores();
     /**
     * Transform intron scores into to 0-1 range, use external boundaries.
     * @param min Minimum possible score, becomes 0
     * @param max Maximum possible score, becomes 1
     *
     */
    void normalizeScores(double min, double max);
    /**
     * Print all introns to the specified file
     * @param output Output file name
     * @param printAll Whether to print all splice sites.
     */
    void printIntrons(string output, bool printAll);

private:

    struct Intron {
        Intron(string protein, string gene, int start, int end,
                char strand, string spliceSites, double score, int number) :
        protein(protein), gene(gene), start(start), end(end), strand(strand),
        spliceSites(spliceSites), score(score), number(number) {
        }
        string protein;
        string gene;
        int start;
        int end;
        char strand;
        string spliceSites;
        double score;
        int number;
    };

    double minScore;
    double maxScore;
    vector<Intron> introns;
};

#endif /* INTRON_STORAGE_H */
