#define CATCH_CONFIG_RUNNER  // This tells Catch to provide a main() - only do this in one cpp file
#include "common.h"
#include "catch.hpp"

std::string PATH;

int main(int argc, char* argv[]) {
    std::string argvStr(argv[0]);
    PATH = argvStr.substr(0, argvStr.find_last_of("/"));
    int result = Catch::Session().run(argc, argv);
    return result;
}
