//
// Created by Матвей on 08.06.2020.
//

#ifndef PLLL_KSVP_H
#define PLLL_KSVP_H

#include <random>
#include <iostream>
#include <vector>
#include <map>
#include <numeric>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <string>
#include <ctime>
#include <cmath>
#include <regex>

#ifdef PLLL_FOUND

#include <plll.hpp>

#endif

#ifdef EIGEN_FOUND

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#endif

using std::vector;
using std::endl;
using IType = long long;
template<class T>
using Matrix = std::vector<std::vector<T>>;
using IMatrix = Matrix<IType>;
using IVector = vector<IType>;
using FType = float;
using FMatrix = Matrix<FType>;
using FVector = vector<FType>;


class ksvp {
private:
    Matrix<FType> GSO;

    IVector lllReduce(vector<FType> &gsreducing, IVector &target, IVector &reducing);

    IVector gaussReduce(IVector &v_new, IMatrix &L, IMatrix &S, vector<double> &lengths);

    void listReduce(IVector &p, IMatrix &L, vector<double> &lengths);

    void findGSO();

    Matrix<FType> findGSO(IMatrix &B);

    IVector gauss(FType bound);

    IVector list(FType bound);

    vector<FType> proj(const vector<FType> &v1, const IVector &v2);

    bool reduce(IVector &p, IMatrix &L, vector<double> &lengths);

    IVector kleinSample(float deviation, IVector target);

    template<typename T, typename T2>
    FType mu(vector<T> &a, vector<T2> &b) {
        double init = 0.0;
        double d = std::inner_product(a.begin(), a.end(), b.begin(), init);
        double dd = std::inner_product(b.begin(), b.end(), b.begin(), init);
        return d / dd;
    }


    template<typename T, typename T2>
    double dotProduct(const vector<T> &v1, const vector<T2> &v2) {
        double ret = 0;
        for (size_t i = 0; i < v1.size(); ++i)
            ret += (v1.at(i) * v2.at(i));
        return ret;
    }


    template<typename T>
    double squaredLengthOfVector(vector<T> &v) {
        return std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
    }

    template<typename T, typename T2>
    vector<T2> subtract(const vector<T> &v1, const vector<T2> &v2) {
        vector<T2> ret(v1.size());
        for (size_t i = 0; i < v1.size(); ++i)
            ret[i] = (v1[i] - v2[i]);
        return ret;
    }

public:
    ksvp();

    ksvp(IMatrix B);

    ~ksvp();

    double lllDelta;
    IMatrix B;
    size_t dimension;
#ifdef EIGEN_FOUND
    int determinant;
    FType norm;
    FType bound;
    FType sum;
#endif
    vector<IType> shortestByGauss;
    vector<IType> shortestByList;

#ifdef PLLL_FOUND
    vector<IType> shortestByBoundList;
    vector<IType> shortestByBoundGauss;
    vector<IType> shortestByVoronoi;

    void listWithBound();

    void voronoi();

    void gaussWithBound();

#endif

    IMatrix lllReduced;

    void gauss();

    void list();


    void lll();

};

template<typename T>
vector<T> operator*(const double &c, const vector<T> &v) {
    vector<T> ret = v;
    for (auto &d : ret)
        d *= c;
    return ret;
}

template<typename T, typename T2>
vector<T> operator-(const vector<T> &v1, const vector<T2> &v2) {
    vector<T> ret(v1.size());
    for (size_t i = 0; i < v1.size(); ++i)
        ret[i] = v1[i] - v2[i];
    return ret;
}


template<typename T>
vector<T> operator+(const vector<T> &v1, const vector<T> &v2) {
    vector<T> ret(v1.size());
    for (size_t i = 0; i < v1.size(); ++i)
        ret[i] = (v1.at(i) + v2.at(i));
    return ret;
}

template<typename T>
void printMatrix(const Matrix<T> &vspace, std::ostream &ost) {
    ost << '[';
    for (auto &v : vspace) {
        ost << '[';
        for (auto &d : v)
            ost << d << ' ';
        ost << "]\n";
    }
    ost << "]\n";
}

template<typename T>
void printMatrix(const Matrix<T> &vspace) {
    printMatrix(vspace, std::cout);
}

template<typename T>
void printVector(const vector<T> &v) {
    printVector(v, std::cout);
}

template<typename T>
void printVector(const vector<T> &v, std::ostream &ost) {
    ost << endl;
    ost << "[" << endl;
    for (auto &d : v)
        ost << d << ", ";
    ost << "]" << endl;
}

#endif //PLLL_KSVP_H
