#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include "excerpt.h"

using namespace std;

//A simple formula for solving quadratic equations using function evaluation
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
        cout << "\t x1 = " << roots[0] << "; x2 = " << roots[1] << endl;
    }
    else{
        std::complex<fp_t> sqrt_f_z = std::sqrt(std::complex<fp_t>((-f_z / a)));
        std::complex<fp_t > root1( z, sqrt_f_z.imag());
        std::complex<fp_t> root2( z, -sqrt_f_z.imag());
        cout << "\tComplex roots: " << endl;
        cout << "x1 = " << root1 << "; x2 = " << root2 << endl;
    }
}

template<typename fp_t>
auto testPolynomial(unsigned int roots_count) {
    fp_t deviation;
    vector<fp_t> roots_computed(roots_count);
    vector<fp_t> roots(roots_count), coefficients(roots_count + 1);
    generate_polynomial<fp_t>(roots_count, 0, roots_count, 0, 1e-5, -1, 1, roots, coefficients);
    simple_formula(coefficients, roots_computed);
    auto result = compare_roots<fp_t>(roots_computed.size(), roots.size(), roots_computed, roots, deviation);
    switch (result) {
        case PR_2_INFINITE_ROOTS:
            cout << "INFINITE ROOTS";
            break;
        case PR_AT_LEAST_ONE_ROOT_IS_FAKE:
            cout << "AT LEAST ONE ROOT IS FAKE";
            break;
        case PR_AT_LEAST_ONE_ROOT_LOST:
            cout << "AT LEAST ONE ROOT LOST";
            break;
        default:
            break;
    }
    return deviation;
}

int main() {
    /*auto p = 2;
    vector<float> roots(p);
    vector<float> roots_computed(p);
    vector<float> coefficients(p + 1);
    auto result = generate_polynomial<float>(p, 0, 2, 0, 10.0 / 5, -10, 10, roots, coefficients);

    simple_formula(coefficients, roots_computed);
    float max_deviation;
    compare_roots<float>(p, p, roots, roots, max_deviation);

    cout << "max deviation: " << max_deviation << endl;

    for (auto root: roots_computed) {
        cout << root << ' ';
    }
    cout << endl;
    for (auto root: roots) {
        cout << root << ' ';
    }*/

    float deviation, max_deviation = 0;
    for (auto i = 0; i < 100; ++i) {
        deviation = testPolynomial<float>(2);
        if (deviation > max_deviation) {
            max_deviation = deviation;
        }
    }

    cout<< "max_deviation = "<< max_deviation<<endl;

   /*
    //Close roots
    vector<float> koef1 = {1, -0.00021, 0.000000011};
    vector<float> r1(2);
    simple_formula(koef1,r1);

    //Complex roots
    vector<float> koef2 = {1, 1, 1};
    vector<float> r2(2);
    simple_formula(koef2,r2);
    */
}