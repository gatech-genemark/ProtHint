#include "common.h"
#include "catch.hpp"
#include "../Parser.h"
#include <stdio.h>

// system(diff) in this testing case is system dependent
// this is ok for testing purposes

int returnDiff(string expected, string result) {
    expected = PATH + "/test_files/" + expected;
    string command = "diff " + expected + " " + result + " >/dev/null";
    return system(command.c_str());
}

TEST_CASE("Test whole program with different settings") {
    Parser fileParser;
    fileParser.setWindowLegth(20);
    fileParser.printAllSites();
    string input = PATH + "/test_files/test_1.ali";
    string output = PATH + "/test_files/test_result";
    ScoreMatrix * scoreMatrix = new ScoreMatrix();
    scoreMatrix->loadFromFile(PATH + "/test_files/blosum62_1.csv");
    Kernel * linearkernel = new LinearKernel();
    Kernel * triangularkernel = new TriangularKernel();
    Kernel * parabolickernel = new ParabolicKernel();
    Kernel * triweightkernel = new TriweightKernel();

    SECTION("Multiplied score") {
        fileParser.setScoringCombination(Parser::BOUNDARIES_MULTIPLIED);
        fileParser.setKernel(linearkernel);
        fileParser.parse(input, output);
        int result = returnDiff("out_1_multiplied_linear.gff", output);
        CHECK(result == 0);

        fileParser.setKernel(triangularkernel);
        fileParser.parse(input, output);
        result = returnDiff("out_1_multiplied_triangular.gff", output);
        CHECK(result == 0);

        fileParser.setKernel(parabolickernel);
        fileParser.parse(input, output);
        result = returnDiff("out_1_multiplied_parabolic.gff", output);
        CHECK(result == 0);
    }

    SECTION("Summation score") {
        fileParser.setScoringCombination(Parser::BOUNDARIES_SUMMED);
        fileParser.setKernel(linearkernel);
        fileParser.parse(input, output);
        int result = returnDiff("out_1_summed_linear.gff", output);
        CHECK(result == 0);

        fileParser.setKernel(triangularkernel);
        fileParser.parse(input, output);
        result = returnDiff("out_1_summed_triangular.gff", output);
        CHECK(result == 0);

        fileParser.setKernel(parabolickernel);
        fileParser.parse(input, output);
        result = returnDiff("out_1_summed_parabolic.gff", output);
        CHECK(result == 0);
    }

    SECTION("BLOSUM62 Multiplication score") {
        fileParser.setScoringCombination(Parser::BOUNDARIES_MULTIPLIED);
        fileParser.setScoringMatrix(scoreMatrix);
        fileParser.setKernel(linearkernel);
        fileParser.parse(input, output);
        int result = returnDiff("out_1_blosum62_multiplied_linear.gff", output);
        CHECK(result == 0);

        fileParser.setKernel(triangularkernel);
        fileParser.parse(input, output);
        result = returnDiff("out_1_blosum62_multiplied_triangular.gff", output);
        CHECK(result == 0);

        fileParser.setKernel(parabolickernel);
        fileParser.parse(input, output);
        result = returnDiff("out_1_blosum62_multiplied_parabolic.gff", output);
        CHECK(result == 0);

        fileParser.setKernel(triweightkernel);
        fileParser.parse(input, output);
        result = returnDiff("out_1_blosum62_multiplied_triweight.gff", output);
        CHECK(result == 0);
    }

    SECTION("BLOSUM62 Summation score") {
        fileParser.setScoringCombination(Parser::BOUNDARIES_SUMMED);
        fileParser.setScoringMatrix(scoreMatrix);
        fileParser.setKernel(linearkernel);
        fileParser.parse(input, output);
        int result = returnDiff("out_1_blosum62_summed_linear.gff", output);
        CHECK(result == 0);

        fileParser.setKernel(triangularkernel);
        fileParser.parse(input, output);
        result = returnDiff("out_1_blosum62_summed_triangular.gff", output);
        CHECK(result == 0);

        fileParser.setKernel(parabolickernel);
        fileParser.parse(input, output);
        result = returnDiff("out_1_blosum62_summed_parabolic.gff", output);
        CHECK(result == 0);
    }

    delete scoreMatrix;
    delete linearkernel;
    delete triangularkernel;
    delete parabolickernel;
    delete triweightkernel;
    remove(output.c_str());
}
