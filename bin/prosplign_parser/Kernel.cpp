#include <cstdlib>
#include <cmath>
#include "Kernel.h"

using namespace std;

void Kernel::setWidth(int width) {
    this->width = width;
}

double LinearKernel::getWeight(int offset) {
    offset = abs(offset);
    if (offset < width) {
        return 1;
    }
    return 0;
}

double TriangularKernel::getWeight(int offset) {
    offset = abs(offset);
    double score = 1 - (double) offset / width;
    if (score < 0) {
        score = 0;
    }
    return score;
}


double ParabolicKernel::getWeight(int offset) {
    offset = abs(offset);
    double score = 1 - ((double) offset / width) * ((double) offset / width);
    if (score < 0) {
        score = 0;
    }
    return score;
}

double TriweightKernel::getWeight(int offset) {
    offset = abs(offset);
    double score = pow(1 - ((double) offset / width) * ((double) offset / width), 3);
    if (score < 0) {
        score = 0;
    }
    return score;
}
