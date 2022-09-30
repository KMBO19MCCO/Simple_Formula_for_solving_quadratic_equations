#include <iostream>
#include <vector>
#include <cmath>
#include "excerpt.h"

using namespace std;

template<typename fp_t>
void simple_formula(std::vector<fp_t> coefficients, std::vector<fp_t> &roots) {
    fp_t a, b, c, z, f_z;

    //Coefficients
    a = coefficients[0];
    b = coefficients[1];
    c = coefficients[2];
    cout << endl << "\tQuadratic equation:" << endl << "\t " << a << "*x^2 + " << b << "*x + " << c << endl;

    z = -b / (2 * a);
    f_z = a * pow(z, 2) + b * z + c;
    if (f_z <= 0) {

        roots[0] = z + sqrt(-f_z / a);
        roots[1] = z - sqrt(-f_z / a);
        cout << "\t\tRoots: " << endl;
        cout << "\t  x1 = " << roots[0] << "; x2 =  " << roots[1] << endl;
    } else
        cout << "\tComplex roots :(..." << endl;


}

int main() {
    auto p = 2;
    vector<float> roots(p);
    vector<float> roots_computed(p);
    vector<float> coefficients(p + 1);
    auto result = generate_polynomial<float>(p, 0, 2, 0, 10.0 / 5, -10, 10, roots, coefficients);

//    vector<float> koef1 = {1, 1, -2};
//    vector<float> koef2 = {1, 6, 9};
//    vector<float> koef3 = {3, -8, -35};
//    vector<float> koef4 = {1, 2, 2};
//    vector<float> koef5 = {1, 2, 1};
//    vector<float> r1(2);
//    vector<float> r2(2);
//    vector<float> r3(2);
//    vector<float> r4(2);
//    vector<float> r5(2);

    //2 different roots
//    simple_formula(koef1, r1);
//    simple_formula(koef3, r3);

    //1 root of multiplicity 2
//    simple_formula(koef2, r2);
//    simple_formula(koef5, r5);

    //complex roots
//    simple_formula(koef4, r4);
    simple_formula(coefficients, roots_computed);
    float max_deviation;
    compare_roots<float>(p, p, roots, roots, max_deviation);
    cout << max_deviation << endl;
    for (auto root: roots_computed) {
        cout << root << ' ';
    }
    cout << endl;
    for (auto root: roots) {
        cout << root << ' ';
    }
}