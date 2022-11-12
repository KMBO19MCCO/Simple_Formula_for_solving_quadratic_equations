// КМБО-03-19 Шмелёва Полина
// почта: polya.shmeliova@gmail.com
// статья:https://www.researchgate.net/publication/337829551_A_simple_formula_for_solving_quadratic_equations_using_function_evaluation

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <limits>
#include "excerpt.h"

//#define DEBUG
#define MAX_DISTANCE 10e-5

using namespace std;

//A simple formula for solving quadratic equations using function evaluation
template<typename fp_t>
int simple_formula(std::vector<fp_t> coefficients, std::vector<fp_t> &roots) {
    fp_t a, b, c, z, f_z;
    int cnt_real_roots = 0;

    //Coefficients - > c, b, a
    a = coefficients[2];
    b = coefficients[1];
    c = coefficients[0];
#ifdef DEBUG
    cout << endl << "\tQuadratic equation:" << endl << "\t " << a << "*x^2 + " << b << "*x + " << c << " = 0 " << endl;
#endif
    if (a != 0) {
        z = -b / (2 * a);

        if (z != std::numeric_limits<fp_t>::infinity()) {
            //Формула: f_z = a * pow(z, 2) + b * z + c;
            fp_t fma_bzc = fma(b, z, c);
            f_z = fma(a, z * z, fma_bzc);

            if (f_z <= 0) {
                fp_t sqrt_f = sqrt(-f_z / a);
                if (sqrt_f != std::numeric_limits<fp_t>::infinity()) {
                    cnt_real_roots = 2; //число действительных корней = 2
                    if (c > 0) { //уравнение имеет 2 одинаковых по знаку корня
                        if (b < 0) { // оба корня положительны (z > 0)
                            roots[0] = z + sqrt_f; // больший по модулю
                            roots[1] = c / roots[0];
                        }
                        if (b > 0) {// оба корня отрицательны (z < 0)
                            roots[0] = z - sqrt_f; // больший по модулю
                            roots[1] = c / roots[0];
                        }
                    }
                    if (c <= 0) { // 2 различных по знаку корня
                        if (b < 0) { // больший по модулю положителен (z > 0)
                            roots[0] = z + sqrt_f; // больший по модулю
                            roots[1] = c / roots[0];
                        }
                        if (b > 0) { // больший по модулю отрицателен (z < 0)
                            roots[0] = z - sqrt_f; // больший по модулю
                            roots[1] = c / roots[0];
                        }
                    }
                }
#ifdef DEBUG
                cout << "x1 = " << roots[0] << ", x2 = " << roots[1] << endl;
#endif
            } else {
                std::complex<fp_t> sqrt_f_z = std::sqrt(std::complex<fp_t>((-f_z / a)));
                std::complex<fp_t> root1(z, sqrt_f_z.imag());
                std::complex<fp_t> root2(z, -sqrt_f_z.imag());
                // выясним, действительно ли корень комплексный:
                // если std::abs(root1.imag())> std::abs(root1)*std::numeric_limits<fp_t>::epsilon() - > корень комплексный, иначе - действительный и полученная
                // комплексная часть - мусор...
                if (std::abs(root1.imag()) > std::abs(root1) * std::numeric_limits<fp_t>::epsilon()) {
                    //значит корень комплексный
                    cnt_real_roots = 0; // число действительных корней = 0
                } else { // корень действительный
#ifdef DEBUG
                    cout << "\t\tNot complex roots! Real roots: " << endl;
#endif
                    roots[0] = roots[1] = root1.real();
                    cnt_real_roots = 2; // число действительных корней = 2
                }
            }
        }
    } else {
#ifdef DEBUG
        cout << "Not Quadratic equation!" << endl;
#endif
        cnt_real_roots = 0;
    }
    return cnt_real_roots;
}

template<typename fp_t>
auto testPolynomial(unsigned int roots_count) {
    fp_t max_absolute_error, max_relative_error;
    vector<fp_t> roots_computed(roots_count);
    vector<fp_t> roots(roots_count), coefficients(roots_count + 1);
    generate_polynomial<fp_t>(roots_count, 0, roots_count, 0, MAX_DISTANCE, -1, 1, roots,
                              coefficients);
    int cnt_real_roots = simple_formula(coefficients, roots_computed);
    if (cnt_real_roots != 0) {
        auto result = compare_roots<fp_t>(roots_computed.size(), roots.size(), roots_computed, roots,
                                          max_absolute_error, max_relative_error);
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

    } else max_absolute_error = std::numeric_limits<fp_t>::infinity();
    return max_absolute_error;
}

int main() {
    float max_deviation = 0;
    for (auto i = 0; i < 100'000; ++i) {
        auto deviation = testPolynomial<float>(2);
        if (deviation != std::numeric_limits<float>::infinity()) {
//            cout << "deviation = " << deviation << endl;
            if (deviation > max_deviation) {
                max_deviation = deviation;
            }
        }
//        else cout << "\t\tComplex roots!" << endl;
    }
    cout << endl << "MAX_deviation = " << max_deviation << endl;

}