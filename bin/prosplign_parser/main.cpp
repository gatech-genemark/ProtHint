#include "Parser.h"
#include "ScoreMatrix.h"
#include "Kernel.h"

#include <iostream>
#include <string>
#include <unistd.h>
#include <stdlib.h>

using namespace std;

void printUsage(char * name) {
    cout << "Usage: " << name << " -i input_file -o output_file "
            "[-w integer] [-am] [-s matrix_file] [-k kernel]" << endl;
    cout << "Options:" << endl;
    cout << "   -w Width of a scoring window around introns." << endl;
    cout << "   -a Print all introns. Without  this flag, only  introns\n"
            "      with canonical splite sites and perfectly scored in-\n"
            "      trons with GC-AG and AT-AC splice sites are printed." << endl;
    cout << "   -m Multiply partial scores from the left and right boun-\n"
            "      daries of the intron. Default behavior is to sum the \n"
            "      scores from both sides." << endl;
    cout << "   -s Use specified scoring matrix for intron scoring. Wi-\n"
            "      thout this option, amino acid scores are determined \n"
            "      based on the quality indicator in the input alignment \n"
            "      file." << endl;
    cout << "   -k Specify type of weighting kernel used. Available opti-\n"
            "      ons are \"triangular\", \"linear\", \"parabolic\" and \n"
            "      \"triweight\". Triangular kernel is the default option." << endl;
}

int main(int argc, char** argv) {
    int opt;
    int windowWidth = 10;
    bool printAllSites = false;
    bool multiply = true;
    string input, output;
    string matrixFile = "";
    string kernelType = "triangular";

    while ((opt = getopt(argc, argv, "i:o:w:ams:k:")) != EOF) {
        switch (opt) {
            case 'i':
                input = optarg;
                break;
            case 'o':
                output = optarg;
                break;
            case 'w':
                windowWidth = atoi(optarg);
                break;
            case 'a':
                printAllSites = true;
                break;
            case 'm':
                multiply = true;
                break;
            case 's':
                matrixFile = optarg;
                break;
            case 'k':
                kernelType = optarg;
                break;
            case '?':
                printUsage(argv[0]);
                return 1;
            default:
                return 1;
        }

    }

    if (input.size() == 0) {
        cerr << "error: Input file not specified" << endl;
        printUsage(argv[0]);
        return 1;
    }

    if (output.size() == 0) {
        cerr << "error: Output file not specified" << endl;
        printUsage(argv[0]);
        return 1;
    }

    if (kernelType != "triangular" && kernelType != "linear" &&
            kernelType != "parabolic" && kernelType != "triweight") {
        cerr << "error: Invalid kernel. Valid options are \"linear\","
                "\"triangular\", \"parabolic\" and \"triweight\" kernels." << endl;
        printUsage(argv[0]);
        return 1;
    }

    ScoreMatrix * scoreMatrix = NULL;
    if (!matrixFile.empty()) {
        scoreMatrix = new ScoreMatrix();
        if (!scoreMatrix->loadFromFile(matrixFile)) {
            cerr << "error: Could not load scoring matrix" << endl;
            printUsage(argv[0]);
            return 1;
        }
    }

    Kernel * kernel;
    if (kernelType == "triangular") {
        kernel = new TriangularKernel();
    } else if (kernelType == "linear") {
        kernel = new LinearKernel();
    } else if (kernelType == "parabolic") {
        kernel = new ParabolicKernel();
    } else if (kernelType == "triweight") {
        kernel = new TriweightKernel();
    }

    Parser fileParser;
    if (multiply) {
        fileParser.setScoringCombination(Parser::BOUNDARIES_MULTIPLIED);
    } else {
        fileParser.setScoringCombination(Parser::BOUNDARIES_SUMMED);
    }
    fileParser.setWindowLegth(windowWidth);
    printAllSites ? fileParser.printAllSites() : fileParser.printCanonicalsites();
    fileParser.setScoringMatrix(scoreMatrix);
    fileParser.setKernel(kernel);

    int result = fileParser.parse(input, output);

    delete scoreMatrix;
    delete kernel;
    return result;
}

