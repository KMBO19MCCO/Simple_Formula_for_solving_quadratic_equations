#define FLT_EPSILON 1.19209290e-07F // float
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
    a = coefficients[2];
    b = coefficients[1];
    c = coefficients[0];
    cout << endl << "\tQuadratic equation:" << endl << "\t " << a << "*x^2 + " << b << "*x + " << c << " = 0 "<< endl;

    if (a!=0) {
        z = -b / (2 * a);
        //Формула: f_z = a * pow(z, 2) + b * z + c;
        fp_t fma_bzc = fma(b, z, c);
        f_z = fma(a, pow(z, 2), fma_bzc);
        cout << "f_z = " << f_z<< endl;
        if (f_z >= 0 && f_z <= FLT_EPSILON) f_z = 0;

        //вычислим другими способами, чтоб посмотреть на точность
        /*fp_t f_z1 = a * pow(z, 2) + b * z + c;
        cout << "f_z1 = " << f_z1<< endl;

        fp_t f_z2 = a * pow(z, 2) + fma_bzc;
        cout << "f_z2 = " << f_z2<< endl;
         */

        if (f_z <= 0 ) {
            roots[0] = z + sqrt(-f_z / a);
            roots[1] = z - sqrt(-f_z / a);
            cout << "\t\tReal roots: " << endl;
            cout << "\t x1 = " << roots[0] << "; x2 = " << roots[1] << endl;
        }
        else {
            std::complex<fp_t> sqrt_f_z = std::sqrt(std::complex<fp_t>((-f_z / a)));
            std::complex<fp_t> root1(z, sqrt_f_z.imag());
            std::complex<fp_t> root2(z, -sqrt_f_z.imag());

            // выясним, действительно ли корень комплексный - попробуем отбросить комплексную часть и подставим в уравнение
            // если результат будет > FLT_EPSILON - > корень комплексный, иначе - действительный и полученная
            // комплексная часть - мусор...

            fp_t fma_b_root_c = fma(b, root1.real(), c);
            fp_t fma_res = fma(a, pow(root1.real(), 2),fma_b_root_c);
            cout << "fma_res = " << fma_res << endl;
            if (fabs(fma_res) <= FLT_EPSILON)
            {
                cout << "\t\tNot complex roots!"<<endl;
                cout<< "\t\tReal roots: "<< endl;
                cout << "x1 = " << root1.real() << "; x2 = " << root2.real() << endl;
            }
            else {
                cout << "\tComplex roots: " << endl;
                cout << "x1 = " << root1 << "; x2 = " << root2 << endl;
            }
        }
    }
    else cout << "Not Quadratic equation!"<< endl;
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

    float deviation, max_deviation = 0;
    for (auto i = 0; i < 500; ++i) {
        deviation = testPolynomial<float>(2);
        cout << "deviation = " << deviation<< endl;
        if (deviation > max_deviation) {
            max_deviation = deviation;
        }
    }
    cout<< endl<<"MAX_deviation = "<< max_deviation<<endl;



    cout << endl<< "\t\tEXAMPLES: "<< endl;

    //Close roots
    //Coef - > c, b, a
    vector<float> koef1 = {0.000000011, -0.00021, 1};
    vector<float> r1(2);
    simple_formula(koef1,r1);

    //Complex roots
    //1*x^2 + 0.163012*x + 0.00664325 = 0
    vector<float> koef2 = {0.00664325, 0.163012, 1};
    vector<float> r2(2);
    simple_formula(koef2,r2);





}