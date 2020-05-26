// #include <iostream>
// #include <plll.hpp>
// // #include <matrix_enumerators.h>
// #include <random>
// #include <random>
// #include <vector>
// #include <map>
// #include <numeric>
// #include <fstream>
// #include <algorithm>
// #include <chrono>
// #include <cstdlib>
// #include <ctime>
// #include <cmath>

// using namespace std;
// //namespace plt = matplotlibcpp;
// template <class T>
// using Matrix = vector<vector<T>>;
// using _type = float;
// constexpr double _c = 0.9;

// int main()
// {
//   const int dimm = 3;
//   plll::linalg::base_matrix<plll::arithmetic::Integer> mat(dimm, dimm);
//   plll::linalg::base_matrix<plll::arithmetic::Integer>::Enumerator iterr = mat.enumerate();
//   while (iterr.has_current())
//   {
//     iterr.current() = plll::arithmetic::Integer(rand());
//     iterr.next();
//   }

//   // for (size_t i = 0; i < dimm; i++)
//   // {
//   //   for (size_t j = 0; j < dimm; j++)
//   //   {
//   //   }
//   // }

//   std::cout << mat;
//   plll::linalg::math_matrix<plll::arithmetic::Integer> A(mat);

//   plll::LatticeReduction lr;
//   lr.setLattice(mat);
//     lr.setSVPMode(plll::LatticeReduction::SVP_VoronoiCellSVP);
//   // lr.svp();

//   // std::cout << lr.getLattice() << "\n";
//   return 0;
// }
#include <plll.hpp>
#include <iostream>

int main()
{
    plll::linalg::math_matrix<plll::arithmetic::Integer> A;
    std::cin >> A;
    
    plll::LatticeReduction lr;
    lr.setLattice(A);
    // lr.lll();

    lr.setSVPMode(plll::LatticeReduction::SVP_VoronoiCellSVP);
    // cout<<"mat";

    // lr.setMinCallbackFunction(f(file, res), f_LI(file, res));
    lr.svp();
    
    std::cout << lr.getLattice() << "\n";
}
