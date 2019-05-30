#include "common.h"
#include "catch.hpp"
#include "../ScoreMatrix.h"
#include <iostream>

using namespace std;

TEST_CASE("Test loading of same matrices from csv "
        "files with different separators") {

    string inputFile1 = PATH + "/test_files/blosum62_1.csv";
    string inputFile2 = PATH + "/test_files/blosum62_2.csv";
    string inputFile3 = PATH + "/test_files/blosum62_3.csv";

    string aminoAcids = "ARNDCQEGHILKMFPSTWYVBZX*";

    ScoreMatrix s1;
    s1.loadFromFile(inputFile1);

    ScoreMatrix s2;
    s2.loadFromFile(inputFile2);

    ScoreMatrix s3;
    s3.loadFromFile(inputFile3);

    for (unsigned int i = 0; i < aminoAcids.size(); i++) {
        for (unsigned int j = 0; j < aminoAcids.size(); j++) {
            CHECK(s1.getScore(aminoAcids[i], aminoAcids[j]) ==
                    s2.getScore(aminoAcids[i], aminoAcids[j]));
            CHECK(s1.getScore(aminoAcids[i], aminoAcids[j]) ==
                    s3.getScore(aminoAcids[i], aminoAcids[j]));
        }
    }
}

TEST_CASE("Detect corrupted file") {
    string inputFile1 = PATH + "/test_files/blosum62_corrupted.csv";
    ScoreMatrix s;
    bool result = s.loadFromFile(inputFile1);
    CHECK_FALSE(result);
}

TEST_CASE("Check random values") {
    string inputFile1 = PATH + "/test_files/blosum62_1.csv";
    ScoreMatrix s;
    s.loadFromFile(inputFile1);
    CHECK (s.getScore('A', 'A') == 4);
    CHECK (s.getScore('S', 'W') == -3);
    CHECK (s.getScore('*', '*') == 1);
    CHECK (s.getScore('X', 'T') == 0);
    CHECK (s.getScore('P', 'I') == -3);
    CHECK (s.getScore('.','A')  == 0);
    CHECK (s.getScore('A','.')  == 0);
    CHECK (s.getScore(' ','A')  == -4);
    CHECK (s.getScore('-',' ')  == 1);
}

TEST_CASE("Check max value for AA") {
    string inputFile1 = PATH + "/test_files/blosum62_1.csv";
    ScoreMatrix s;
    s.loadFromFile(inputFile1);
    CHECK (s.getMaxScore() == 11);
}
