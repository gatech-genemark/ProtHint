#ifndef PARSER_H
#define PARSER_H

#include "Alignment.h"
#include "ScoreMatrix.h"
#include "Kernel.h"
#include <string>

#define READ_SUCCESS 0
#define OPEN_FAIL 1
#define FORMAT_FAIL 2
#define NO_MORE_ALIGNMENTS 3
#define COMPLETE_MATCH 1;
#define GOOD_MATCH 0.5;
#define BAD_MATCH 0;
#define GAP -0.5;

using namespace std;

/// Class for parsing the whole file

class Parser {
public:
    Parser();
    /**
     * Parse the alignment file, save the output into output stream.
     * @param inputFile Name of the alignment file
     * @param outputFile Name of the gff output file
     */
    int parse(string inputFile, string outputFile);

    /**
     * Set how scores from left and right intron boundary are combined
     */
    void setScoringCombination(int combination);
    /**
     * Set window length in the scoring function
     */
    void setWindowLegth(int length);

    /**
     * Print all splice sites during parsing
     */
    void printAllSites();
    /**
     * Print only canonical splice sites during parsing
     */
    void printCanonicalsites();

    /**
     * Assign scoring matrix to the parser
     */
    void setScoringMatrix(const ScoreMatrix * scoreMatrix);

    /**
     * Set which kernel will be used for scoring
     */
    void setKernel(Kernel * kernel);

    static const int BOUNDARIES_SUMMED = 1;
    static const int BOUNDARIES_MULTIPLIED = 2;

private:
    /**
     * Parse next alignment in the input file.
     * The alignment is stored in the "alignment" class variable
     */
    int parseNext();

    Alignment alignment;
    int scoreCombination;
    int windowLength;
    bool printAll;
    const ScoreMatrix * scoreMatrix;
    Kernel * kernel;

    /**
     * Return maximum possible score for an intron, depending on whether
     * a scoring matrix is used
     * @return
     */
    double maxScore();

    ifstream inputStream;
    static const char PAIR_SEPARATOR[];
    static const int NUM_SEPARATORS = 3;
};


#endif /* PARSER_H */
