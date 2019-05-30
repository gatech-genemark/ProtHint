#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <fstream>
#include <vector>
#include "ScoreMatrix.h"
#include "IntronStorage.h"
#include "Kernel.h"

using namespace std;

/// Class for parsing a single gene-protein alignment

class Alignment {
public:
    Alignment();
    /**
     * Parse a single gene-protein alignment.
     * The function checks if the general structure of the alignment
     * is OK but it does not check the validity of every single base/protein.
     *
     * @param fstream File stream starting at the position of 
     *                the alignment start.
     */
    int parse(ifstream & inputStream);
    /**
     * @return Name of the aligned gene
     */
    string getGene();
    /**
     * @return Name of the aligned protein
     */
    string getProtein();
    /**
     * @return Total alignment length in nucleotides
     */
    int getLength();
    /**
     * Print the whole alignment
     */
    void print(ostream & os);
    /**
     * Store introns found in this alignment to the IntronStorage object
     * which contains all introns
     * @param storeage The intron storage
     */
    void storeIntrons(IntronStorage & storage);
    /**
     * Score all introns in the alignment
     * @param windowWidth Number of amino acids scored in the upstream and
     *                    downstream regions
     * @param multiply    Whether to multiply or sum scores from upstream and
     *                    downstream regions
     * @param scoreMatrix Scoring matrix used for scoring amino acids. If no
     *                    matrix is given, score is determined based on the
     *                    quality indicator in the input alignment file
     */
    void scoreIntrons(int windowWidth, bool multiply,
            const ScoreMatrix * scoreMatrix, Kernel * kernel);
    /**
     * @return True if alignment contains any introns
     */
    bool hasIntrons();
private:
    /// Starting position of the alignment in the gene
    int alignmentStart;
    int i;
    string gene;
    string protein;
    /**
     * Clear the object for new alignment pair
     */
    void clear();

    /// Overall position in alignment
    int index;
    int realPositionCounter;

    /// Single nucleotide-protein pair
    struct AlignedPair {
        /**
         * Save pair and determine exon/intron
         */
        AlignedPair(char n, char tc, char q,
                char p, char type);
        /**
         * Return amino acid score. Either based on scoring matrix or alignment
         * quality from the alignment file if no matrix is given
         * @return AA score
         */
        double score(const ScoreMatrix * scoreMatrix);
        char nucleotide;
        /**
         * The proteins are saved as follows: 1P3 or ppp
         * Where P is the protein, numbers 1 and 3 fill
         * the space created by 3 nucleotide to 1 protein mapping.
         * Lowercase letters are used at splice sites.
         */
        char translatedCodon;
        /**
         * Quality of a single AA alignment.
         * '|' for complete match
         * '+' for good match
         * ' ' for bad match
         */
        double quality;
        char protein;
        /** 
         * Type of DNA base (based on alignment)
         * 'i' for intron
         * 'e' for exon
         * ' ' for unknown
         */
        char type;
        /// Position of a nucleotide in the alignment relative to a gene start
        int realPosition;
    };

    /**
     * Intron structure for the purposes of the alignment only. Final set of
     * introns from all alignments with correct start and end positions
     * relative to gene start is stored in the IntronStorage class
     */
    struct Intron {
        Intron();
        unsigned int start, end;
        double score;
        bool scoreSet;
        char donor[2];
        char acceptor[2];
        double leftScore, rightScore;
        double leftWeightSum, rightWeightSum;
    };

    /// Initial size of alignment vector
    static const int N = 3000;
    /// Array containing the actual alignment
    vector<AlignedPair> pairs;


    vector<string> header;
    static const int HEADER_SIZE = 6;
    int parseHeader(const string & headerString);

    int readAlignmentStart(const string & headerString);

    static const int BLOCK_LENGTH = 100;
    static const int BLOCK_ITEMS_CNT = 5;
    static const int BLOCK_OFFSET = 12;
    /**
     *  Parse individual block of lines containing the alignment and its properties
     */
    void parseBlock(const vector<string>& lines);

    void printLineError();

    /**
     * Fill spaces in codon to protein mapping by numbers
     * corresponding to position in the codon
     */
    void mapCodons(AlignedPair & pair);

    bool insideIntron;
    // Flag indicating that donor position of an intron is being read
    bool donorFlag;
    /**
     * Detect and save introns.
     * The function also retrieves information associated with the intron
     * such as its start/end position and donor/acceptor site.
     * @param pair
     */
    void checkForIntron(AlignedPair & pair);
    vector<Intron> introns;
    const ScoreMatrix * scoreMatrix;
    Kernel * kernel;

    /**
     * Determine score of a single intron using exon alignment in the
     * upstream and downstream region
     */
    double scoreIntron(Intron & intron, int windowWidth, bool multiply);
    /// Alignment score of amino acids in exon before intron start
    void scoreLeft(Intron & intron, int start, int windowWidth);
    /// Alignment score of amino acids in exon after intron end
    void scoreRight(Intron & intron, int start, int windowWidth);
};


#endif /* ALIGNMENT_H */


