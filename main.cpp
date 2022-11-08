// КМБО-03-19 Шмелёва Полина
// почта: polya.shmeliova@gmail.com
// статья:https://www.researchgate.net/publication/337829551_A_simple_formula_for_solving_quadratic_equations_using_function_evaluation

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <limits>
#include "excerpt.h"

using namespace std;

//A simple formula for solving quadratic equations using function evaluation
template<typename fp_t>
void simple_formula(std::vector<fp_t> coefficients, std::vector<fp_t> &roots) {
    fp_t a, b, c, z, f_z;

    //Coefficients - > c, b, a
    a = coefficients[2];
    b = coefficients[1];
    c = coefficients[0];
    cout << endl << "\tQuadratic equation:" << endl << "\t " << a << "*x^2 + " << b << "*x + " << c << " = 0 "<< endl;

    if (a!=0) {
        z = -b / (2 * a);
        //Формула: f_z = a * pow(z, 2) + b * z + c;
        fp_t fma_bzc = fma(b, z, c);
        f_z = fma(a, z*z, fma_bzc);

        if (f_z <= 0 ) {
            fp_t sqrt_f = sqrt(-f_z / a);
            //z >= 0 ? roots[0] = z + sqrt_f , roots[1] =  roots[0] - 2*sqrt_f : roots[0] = z - sqrt_f, roots[1] = roots[0] + 2*sqrt_f;

            if(z >= 0) {
                roots[0] = z + sqrt_f;
                roots[1] =  roots[0] - 2*sqrt_f;
            }
            else{
                roots[0] = z - sqrt_f;
                roots[1] = roots[0] + 2*sqrt_f;
            }
            cout << "\t\tReal roots: " << endl;
            cout << "\t x1 = " << roots[0] << "; x2 = " << roots[1] << endl;
        }
        else {

            std::complex<fp_t> sqrt_f_z = std::sqrt(std::complex<fp_t>((-f_z / a)));
            std::complex<fp_t> root1(z, sqrt_f_z.imag());
            std::complex<fp_t> root2(z, -sqrt_f_z.imag());

            // выясним, действительно ли корень комплексный:
            // если fabs(root1)*std::numeric_limits<fp_t>::epsilon() < fabs(root1.imag()) - > корень комплексный, иначе - действительный и полученная
            // комплексная часть - мусор...

            if (fabs(root1)*std::numeric_limits<fp_t>::epsilon() > fabs(root1.imag()))
            {
                //значит корень действительный
                cout << "\t\tNot complex roots!"<<endl;
                cout<< "\t\tReal roots: "<< endl;
                roots[0] = roots[1] =  root1.real();

                cout << "\t x1 = " << roots[0] << "; x2 = " << roots[1] << endl;
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
    generate_polynomial<fp_t>(roots_count, 0, roots_count, 0, std::numeric_limits<fp_t>::epsilon(), -1, 1, roots, coefficients);
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
    for (auto i = 0; i < 10'000; ++i) {
        deviation = testPolynomial<float>(2);
        cout << "deviation = " << deviation<< endl;
        if (deviation > max_deviation) {
            max_deviation = deviation;
        }
    }
    cout<< endl<<"MAX_deviation = "<< max_deviation<<endl;


    //Complex roots
    //1*x^2 + 0.163012*x + 0.00664325 = 0
    /*vector<float> koef2 = {0.00664325, 0.163012, 1};
    vector<float> r2(2);
    simple_formula(koef2,r2);
    */

}