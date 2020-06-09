//
// Created by Матвей on 08.06.2020.
//
#include "ksvp.h"

IMatrix generateLattice(int dimm, int diff) {
    IMatrix m(dimm, IVector(dimm, 0));
    int summ = 0;
    srand(time(0));
    for (auto &v : m) {
        for (auto &el : v) {
            float f = rand();
            el = f / RAND_MAX * (float) diff;
            summ += el;
        }
    }
    return m;
}

int main() {
    IMatrix a = generateLattice(5, 50);
    printMatrix(a);
    ksvp solver = ksvp(a);
    solver.lll();
    printMatrix(solver.lllReduced);
    return 0;
}
