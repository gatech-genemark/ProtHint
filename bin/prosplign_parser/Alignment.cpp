#include "Alignment.h"
#include "Parser.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <cmath>

using namespace std;

Alignment::Alignment() {
    pairs.reserve(N);
}

void Alignment::clear() {
    header.clear();
    index = 0;
    insideIntron = false;
    donorFlag = false;
    introns.clear();
    alignmentStart = 0;
}

int Alignment::parse(ifstream& inputStream) {
    clear();
    bool geneStarted = false;
    vector<string> blockLines(BLOCK_ITEMS_CNT);

    // Read header
    string line;
    int status = FORMAT_FAIL;
    if (getline(inputStream, line)) {
        status = parseHeader(line);
    }
    if (status != READ_SUCCESS) {
        cerr << "error: Invalid alignment header " << endl;
        return status;
    }

    // Read content
    while (getline(inputStream, line)) {
        if (line.empty()) {
            break;
        }

        if (!geneStarted) {
            if (line[0] != '-') {
                if (readAlignmentStart(line) != READ_SUCCESS) {
                    cerr << "warning: error in alignment " << header[1] << "-" << header[2];
                    cerr << ": Invalid alignment start" << endl;
                    return FORMAT_FAIL;
                }
                geneStarted = true;
            }
        }

        // Read block of 5 five lines containing alignment
        // The last block might be shorter than BLOCK_LENGTH
        unsigned int blockLength = line.find(" ", BLOCK_OFFSET) - BLOCK_OFFSET;
        if (line.size() <= BLOCK_OFFSET) {
            printLineError();
            return FORMAT_FAIL;
        }
        blockLines[0] = line.substr(BLOCK_OFFSET, blockLength);
        for (unsigned i = 1; i < BLOCK_ITEMS_CNT; i++) {
            if (getline(inputStream, line) && !line.empty() && line.size() > BLOCK_OFFSET) {
                blockLines[i] = line.substr(BLOCK_OFFSET, blockLength);
                if (blockLines[i].size() != blockLength) {
                    printLineError();
                    return FORMAT_FAIL;
                }
            } else {
                printLineError();
                return FORMAT_FAIL;
            }
        }
        // Parse this block
        parseBlock(blockLines);
    }

    if (!geneStarted) {
        cerr << "warning: error in alignment " << header[1] << "-" << header[2];
        cerr << ": No alignment" << endl;
        return FORMAT_FAIL;
    }
    return READ_SUCCESS;
}

int Alignment::parseHeader(const string& headerString) {
    stringstream ss(headerString);
    string column;

    while (getline(ss, column, '\t')) {
        header.push_back(column);
    }

    if (header.size() != HEADER_SIZE) {
        return FORMAT_FAIL;
    }
    return READ_SUCCESS;
}

int Alignment::readAlignmentStart(const string& headerString) {
    stringstream ss(headerString);
    ss >> alignmentStart;
    if (ss.fail()) {
        return FORMAT_FAIL;
    }
    realPositionCounter = alignmentStart;
    return READ_SUCCESS;
}

void Alignment::parseBlock(const vector<string>& lines) {
    // Parse individual pairs
    for (unsigned i = 0; i < lines[0].size(); i++) {
        AlignedPair pair(lines[0][i], lines[1][i], lines[2][i],
                lines[3][i], lines[4][i]);

        mapCodons(pair);
        checkForIntron(pair);

        if (pair.nucleotide != '-') {
            pair.realPosition = realPositionCounter++;
        }

        // Reuse space if possible
        if ((int) pairs.size() <= index) {
            pairs.push_back(pair);
        } else {
            pairs[index] = pair;
        }
        index++;
    }
}

void Alignment::printLineError() {
    cerr << "warning: error in alignment " << header[1] << "-" << header[2];
    cerr << ": corrupted alignment - wrong line length.";
    cerr << " The rest of this alignment is skipped." << endl;
}

void Alignment::mapCodons(AlignedPair& pair) {
    if (index > 0) {
        if (pair.translatedCodon >= 'A' && pair.translatedCodon <= 'Z') {
            pairs[index - 1].translatedCodon = '1';
        }
        if (pair.protein >= 'A' && pair.protein <= 'Z') {
            pairs[index - 1].protein = '1';
        }
        if (pair.translatedCodon == ' ' && pairs[index - 1].translatedCodon >= 'A'
                && pairs[index - 1].translatedCodon <= 'Z') {
            pair.translatedCodon = '3';
        }
        if (pair.protein == ' ' && pairs[index - 1].protein >= 'A'
                && pairs[index - 1].protein <= 'Z') {
            pair.protein = '3';
        }
    }
}

void Alignment::checkForIntron(AlignedPair& pair) {
    if (donorFlag) {
        introns.back().donor[1] = pair.nucleotide;
        donorFlag = false;
        // If the read is at donor position, there is nothing else to check for
        return;
    }
    if (!insideIntron && pair.type == 'i') { // intron start
        Intron i;
        i.start = index;
        i.donor[0] = pair.nucleotide;
        introns.push_back(i);
        insideIntron = true;
        donorFlag = true;
    }  else if (insideIntron && pair.type != 'i') { // intron end
        introns.back().end = index - 1;
        introns.back().acceptor[0] = pairs[index - 2].nucleotide;
        introns.back().acceptor[1] = pairs[index - 1].nucleotide;
        insideIntron = false;
    }
}

void Alignment::storeIntrons(IntronStorage& storage) {
    for (unsigned int i = 0; i < introns.size(); i++) {
        string spliceSites(introns[i].donor, 2);
        spliceSites.append("_");
        spliceSites.append(introns[i].acceptor, 2);

        if (!introns[i].scoreSet) {
            introns[i].score = 0;
        }

        storage.storeIntron(header[2], header[1],
                pairs[introns[i].start].realPosition,
                pairs[introns[i].end].realPosition,
                '+', spliceSites, introns[i].score, i + 1);
    }
}

void Alignment::print(ostream& os) {
    for (int i = 0; i < index; i += BLOCK_LENGTH) {
        for (int j = 0; (j < BLOCK_LENGTH && i + j < index); j++) {
            os << pairs[i + j].nucleotide;
        }
        os << endl;
        for (int j = 0; (j < BLOCK_LENGTH && i + j < index); j++) {
            os << pairs[i + j].translatedCodon;
        }
        os << endl;
        for (int j = 0; (j < BLOCK_LENGTH && i + j < index); j++) {
            os << pairs[i + j].protein;
        }
        os << endl;
        for (int j = 0; (j < BLOCK_LENGTH && i + j < index); j++) {
            os << pairs[i + j].type;
        }
        os << endl;
        for (int j = 0; (j < BLOCK_LENGTH && i + j < index); j++) {
            os << '*';
        }
        os << endl;
    }
}

bool Alignment::hasIntrons() {
    return introns.size() != 0;
}

string Alignment::getGene() {
    return header[1];
}

string Alignment::getProtein() {
    return header[2];
}

int Alignment::getLength() {
    return index;
}

Alignment::Intron::Intron() {
    scoreSet = false;
}

void Alignment::scoreIntrons(int windowWidth, bool multiply,
        const ScoreMatrix * scoreMatrix, Kernel * kernel) {
    this->scoreMatrix = scoreMatrix;
    this->kernel = kernel;
    this->kernel->setWidth(windowWidth);
    for (unsigned int i = 0; i < introns.size(); i++) {
        if (!introns[i].scoreSet) {
            scoreIntron(introns[i], windowWidth, multiply);
        }
    }
}

double Alignment::scoreIntron(Intron& intron, int windowWidth, bool multiply) {
    intron.leftScore = intron.rightScore = 0;
    intron.leftWeightSum = intron.rightWeightSum = 0;
    int left, right;

    // Determine if codon is split and how
    if (pairs[intron.start - 1].protein == '3' ||
            pairs[intron.start - 1].translatedCodon == '3') {
        // Codon is not split
        left = intron.start - 2;
        right = intron.end + 2;
    } else {
        if (pairs[intron.start - 2].protein == '3'
                || pairs[intron.start - 2].translatedCodon == '3') {
            // Codon is split after the first nucleotide
            left = intron.start - 3;
            right = intron.end + 4;
            // Majority of the codon belongs to the right side
            double weight = kernel->getWeight(0);
            intron.rightScore += pairs[intron.start - 1].score(scoreMatrix) * weight;
            intron.rightWeightSum += weight;
        } else {
            // Codon is split after the second nucleotide
            left = intron.start - 4;
            right = intron.end + 3;
            // Majority of the codon belongs to the left side
            double weight = kernel->getWeight(0);
            intron.leftScore += pairs[intron.start - 1].score(scoreMatrix) * weight;
            intron.leftWeightSum += weight;
        }
    }

    scoreLeft(intron, left, windowWidth);
    scoreRight(intron, right, windowWidth);

    // Normalize alignments by their length, otherwise alignments
    // in short exons between introns are penalized
    if (multiply) {
        if (intron.leftWeightSum == 0 || intron.rightWeightSum == 0 ||
                intron.leftScore < 0 || intron.rightScore < 0) {
            intron.score = 0;
        } else {
            intron.score = (intron.leftScore / (intron.leftWeightSum)) *
                    (intron.rightScore / (intron.rightWeightSum));
            intron.score = sqrt(intron.score);
        }
    } else {
        intron.score = (intron.leftScore + intron.rightScore) /
                (intron.leftWeightSum + intron.rightWeightSum);
        if (intron.score < 0) {
            intron.score = 0;
        }
    }

    intron.scoreSet = true;
    return intron.score;
}

void Alignment::scoreLeft(Intron & intron, int start, int windowWidth) {
    for (int i = start; i > (start - windowWidth * 3); i -= 3) {
        // Check for end of local alignment
        if (i < 0 || pairs[i].type != 'e') {
            return;
        }
        double weight = kernel->getWeight((i - start) / 3);
        intron.leftWeightSum += weight;
        intron.leftScore += pairs[i].score(scoreMatrix) * weight;
    }
}

void Alignment::scoreRight(Intron & intron, int start, int windowWidth) {
    for (int i = start; i < (start + windowWidth * 3); i += 3) {
        // Check for end of local alignment
        if (i >= index || pairs[i].type != 'e') {
            return;
        }
        double weight = kernel->getWeight((i - start) / 3);
        intron.rightWeightSum += weight;
        intron.rightScore += pairs[i].score(scoreMatrix) * weight;
    }
}

Alignment::AlignedPair::AlignedPair(char n, char tc, char q, char p, char type) :
nucleotide(n),
translatedCodon(tc),
protein(p) {
    if (type == '*') {
        if (protein == '.') {
            this->type = 'i';
        } else {
            this->type = 'e';
        }
    } else {
        this->type = ' ';
    }

    // Assign random amino acid to a stop codon, so it is treated as a gap
    if (translatedCodon == '*') {
        translatedCodon = 'A';
    }

    quality = BAD_MATCH;
    if (q == '|') {
        quality = COMPLETE_MATCH;
    } else if (q == '+') {
        quality = GOOD_MATCH;
    } else if (nucleotide == '-' || protein == '-') {
        quality = GAP;
    }
}

double Alignment::AlignedPair::score(const ScoreMatrix * scoreMatrix) {
    if (scoreMatrix != NULL) {
        return scoreMatrix->getScore(translatedCodon, protein);
    }
    return quality;
}
