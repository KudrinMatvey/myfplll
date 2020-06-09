//
// Created by Матвей on 08.06.2020.
//

#include "ksvp.h"

ksvp::ksvp() {}

#ifdef PLLL_FOUND

ksvp::ksvp(IMatrix B) {
    this->B = B;
    this->dimension = B.size();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m(this->dimension, this->dimension);
    for (size_t i = 0; i < this->dimension; i++) {
        for (size_t j = 0; j < this->dimension; j++) {
            m(i, j) = B[i][j];
        }
    }
    this->lllDelta = 0.9;
    this->norm = m.norm();
    this->sum = m.sum();
    this->determinant = m.determinant();
    this->bound = 1.01 * (sqrt(this->dimension)) * pow(abs(this->determinant), 1.0 / this->dimension);
}

#else
ksvp::ksvp(IMatrix B) {
    this->B = B;
    this->lllDelta = 0.9;
    this->dimension = B.size();
    this->bound = 1.01;
}
#endif

ksvp::~ksvp() {}


vector<FType> ksvp::proj(const vector<FType> &v1, const IVector &v2) {
    vector<FType> _proj = (dotProduct(v1, v2) / dotProduct(v1, v1)) * v1;
    return _proj;
}


void ksvp::findGSO() {
    this->GSO = findGSO(this->B);
}

Matrix<FType> ksvp::findGSO(IMatrix &B) {
    int dimension = B.size();
    Matrix<FType> uspace(dimension, vector<FType>(dimension));
    for (int i = 0; i < dimension; i++) {
        uspace[0][i] = B[0][i];
    }
    for (int i = 1; i < dimension; ++i) {
        for (int j = 0; j < i; ++j) {
            uspace.at(i) = subtract(B.at(i), proj(uspace.at(j), B.at(i)));
        }
    }
    return uspace;
}

bool ksvp::reduce(IVector &p, IMatrix &L, vector<double> &lengths) {
    double s = squaredLengthOfVector(p);
    for (size_t i = 0; i < L.size(); i++) {
        if (s >= lengths[i]) {
            auto diff = p - L[i];
            double s2 = squaredLengthOfVector(diff);
            if (s2 < s) {
                p = diff;
                return true;
            }
        }
    }
    return false;
}

IVector ksvp::gaussReduce(IVector &v_new, IMatrix &L, IMatrix &S, vector<double> &lengths) {
    bool cont = true;
    while (cont) {
        auto same = std::find_if(L.begin(), L.end(), [&v_new](IVector &v) { return v_new == v; });
        if (same != L.end()) {
            return IVector(v_new.size());
        }
        cont = reduce(v_new, L, lengths);
    }
    auto it = L.begin();
    auto it2 = lengths.begin();
    int i = 0;
    IVector a(v_new.size());
    while (it != L.end()) {
        double v_new_length = squaredLengthOfVector(v_new);
        IVector v = *it;
        if (lengths[i] > v_new_length) {
            a = v - v_new;
            if (squaredLengthOfVector(a) <= lengths[i]) {
                S.push_back(a);
                it = L.erase(it);
                it2 = lengths.erase(it2);
            } else {
                ++it;
                ++it2, ++i;
            }
        } else {
            ++it;
            ++it2, ++i;
        }
    }
    return v_new;
}

void ksvp::listReduce(IVector &p, IMatrix &L, vector<double> &lengths) {
    bool cont = true;
    while (cont) {
        auto same = std::find(L.begin(), L.end(), p);
        if (same != L.end()) {
            p = IVector(p.size());
            return;
        }
        if (!reduce(p, L, lengths))
            return;
    }
}

IVector ksvp::gauss(FType bound) {
    IMatrix L, S, B;
    B = this->B;
    IVector v_new;
    if (this->GSO.size() == 0) {
        this->findGSO();
    }
    vector<double> squarelengthsVector;
    size_t K = 0, iter = 0;

    size_t max_number_of_collistions = 200;
    do {
        iter++;
        if (S.size() > 0) {
            v_new = S.back();
            S.pop_back();
        } else {
            IVector a(B[0].size(), 0);
            v_new = kleinSample(1, a);
        }
        v_new = gaussReduce(v_new, L, S, squarelengthsVector);
        if (!v_new.size() ||
            std::count(v_new.begin(), v_new.end(), 0) == v_new.size()) {
            K++;
        } else if (sqrt(squaredLengthOfVector(v_new)) < bound) {
            return v_new;
        } else {
            L.push_back(v_new);
            squarelengthsVector.push_back(squaredLengthOfVector(v_new));
        }
        max_number_of_collistions = std::max(max_number_of_collistions, L.size());
    } while (10 * K < max_number_of_collistions + 2000 || L.size() == 0);
    IVector shortestVector = L[0];
    int count = 0;
    for (size_t i = 1; i < L.size(); i++) {
        if (sqrt(squaredLengthOfVector(L[i])) <= bound) {
            count++;
        }
        if (squaredLengthOfVector(shortestVector) > squaredLengthOfVector(L[i])) {
            shortestVector = L[i];
        }
    }
    return shortestVector;
}

void ksvp::lll() {
    IMatrix B = this->B;
    bool cont;
    FMatrix gramSmidt;
    int iter = 0;
    do {
        iter++;
        gramSmidt = findGSO(B);
        cont = false;
        for (size_t i = 1; i < B.size(); i++) {
            for (int j = 0; j < i; j++) {
                B.at(i) = lllReduce(gramSmidt.at(j), B.at(i), B.at(j));
            }
        }
        for (size_t i = 1; i < B.size() && !cont; i++) {
            float a1 = squaredLengthOfVector(gramSmidt.at(i - 1));
            float a2 = squaredLengthOfVector(gramSmidt.at(i));
            float mu1 = mu(B[i], gramSmidt[i - 1]);
            if (a2 < (this->lllDelta - mu1 * mu1) * a1) {
                B.at(i - 1).swap(B.at(i));
                cont = true;
            }
        }
    } while (cont);
    this->lllReduced = B;
}

void ksvp::list() {
    this->shortestByList = list(1);
}

void ksvp::gauss() {
    this->shortestByGauss = gauss(1);
}

#ifdef PLLL_FOUND

void ksvp::gaussWithBound() {
    this->shortestByBoundGauss = gauss(this->bound);
}

void ksvp::listWithBound() {
    this->shortestByBoundList = list(this->bound);
}


void ksvp::voronoi() {
    IVector res(this->dimension);
    plll::linalg::base_matrix<plll::arithmetic::Integer> mat(this->dimension, this->dimension);
    plll::linalg::base_matrix<plll::arithmetic::Integer>::Enumerator matrIter = mat.enumerate();

    int i = 0, x, y;
    while (matrIter.has_current()) {
        x = i / this->dimension;
        y = i % this->dimension;
        matrIter.current() = plll::arithmetic::Integer(this->B[x][y]);
        matrIter.next();
        i++;
    }
    plll::LatticeReduction lr;
    plll::linalg::math_matrix<plll::arithmetic::Integer> A(mat);
    lr.setLattice(A);
    lr.setSVPMode(plll::LatticeReduction::SVP_VoronoiCellSVP);
    lr.svp();
    mat = lr.getLattice();
    matrIter = mat.enumerate();
    for (size_t i = 0; i < this->dimension; i++) {
        res[i] = plll::arithmetic::convert<int>(matrIter.current());
        matrIter.next();
    }
    this->shortestByVoronoi = res;
}

#endif

IVector ksvp::list(FType bound) {
    if (this->GSO.size() == 0) {
        this->findGSO();
    }

    IMatrix L;
    IVector a(B[0].size(), 0);
    vector<double> lengths;
    size_t max_number_of_collistions = 200;
    int iter = 0;
    for (size_t i = 0, k = 0; 10 * i < max_number_of_collistions + 2000 || L.size() == 0; k++) {
        iter++;
        IVector v = kleinSample(1, a);
        listReduce(v, L, lengths);
        bool isZero = true;
        for (auto in : v) {
            if (in != 0) {
                isZero = false;
                break;
            }
        }
        if (isZero) {
            i++;
        } else if (sqrt(squaredLengthOfVector(v)) < bound) {
            return v;
        } else {
            L.push_back(v);
            lengths.push_back(squaredLengthOfVector(v));
        }
        max_number_of_collistions = std::max(max_number_of_collistions, L.size());
    }
    IVector shortestVector = L[0];
    int count = 0;
    for (size_t i = 1; i < L.size(); i++) {
        if (sqrt(squaredLengthOfVector(L[i])) <= bound) {
            count++;
        }
        if (squaredLengthOfVector(shortestVector) > squaredLengthOfVector(L[i])) {
            shortestVector = L[i];
        }
    }
    return shortestVector;
}


IVector ksvp::kleinSample(float deviation, IVector target) {
    IVector v(target.size(), 0);
    std::random_device mch;
    std::default_random_engine generator(mch());

    double d, z;
    for (int i = this->B.size() - 1; i >= 0; i--) {
        double length = squaredLengthOfVector(GSO.at(i));
        d = dotProduct(target, this->GSO.at(i)) / length;
        std::normal_distribution<double> distribution(d, 1);
        z = round(distribution(generator));
        target = target - z * this->B[i];
        v = v + (z * this->B[i]);
    }
    return v;
}

IVector ksvp::lllReduce(vector<FType> &gsreducing, IVector &target,
                        IVector &reducing) {
    IVector t = target - round(mu(target, gsreducing)) * reducing;
    return t;
}




