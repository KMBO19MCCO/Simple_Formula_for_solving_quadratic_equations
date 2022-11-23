// КМБО-03-19 Шмелёва Полина
// почта: polya.shmeliova@gmail.com
// статья:https://www.researchgate.net/publication/337829551_A_simple_formula_for_solving_quadratic_equations_using_function_evaluation

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <limits>
#include "excerpt.h"


//#define DEBUG // раскомментировать для вывода отладочного кода (само уравнение+корни)
#define MAX_DISTANCE 10e-5

using namespace std;

//A simple formula for solving quadratic equations using function evaluation
// МЕТОД С ВЕЩЕСТВЕННЫМИ ВЫЧИСЛЕНИЯМИ
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
                            roots[0] = z + sqrt_f; // бОльший по модулю
                            roots[1] = c / roots[0];
                        }
                        if (b > 0) {// оба корня отрицательны (z < 0)
                            roots[0] = z - sqrt_f; // бОльший по модулю
                            roots[1] = c / roots[0];
                        }
                    }
                    if (c <= 0) { // 2 различных по знаку корня
                        if (b < 0) { // больший по модулю положителен (z > 0)
                            roots[0] = z + sqrt_f; // бОльший по модулю
                            roots[1] = c / roots[0];
                        }
                        if (b > 0) { // больший по модулю отрицателен (z < 0)
                            roots[0] = z - sqrt_f; // бОльший по модулю
                            roots[1] = c / roots[0];
                        }
                    }
                }
#ifdef DEBUG
                cout << "x1 = " << roots[0] << ", x2 = " << roots[1] << endl;
#endif
            } else {
                std::complex<fp_t> sqrt_f_z = std::sqrt(std::complex<fp_t>((-f_z / a)));
                if (sqrt_f_z != std::numeric_limits<fp_t>::infinity()){
                    std::complex<fp_t> root1(z, sqrt_f_z.imag());
                    std::complex<fp_t> root2(z, -sqrt_f_z.imag());
                    // выясним, действительно ли корень комплексный:
                    // если std::abs(root1.imag())> std::abs(root1)*std::numeric_limits<fp_t>::epsilon() - > корень комплексный, иначе - действительный и полученная
                    // комплексная часть - мусор...
                    if (std::abs(root1.imag()) > std::abs(root1) * std::numeric_limits<fp_t>::epsilon()) {
                        //значит корень комплексный
                        cnt_real_roots = 0; // число действительных корней = 0
                    }else { // корень действительный
#ifdef DEBUG
                        cout << "\t\tNot complex roots! Real roots: " << endl;
#endif
                        roots[0] = roots[1] = root1.real();
                        cnt_real_roots = 2; // число действительных корней = 2
                    }
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


// МЕТОД С КОМПЛЕКСНЫМИ ВЫЧИСЛЕНИЯМИ
template<typename fp_t>
std::vector<std::complex<fp_t>>  simple_formula_complex(std::vector<fp_t> coefficients, std::vector<std::complex<fp_t>> &roots)
{
    fp_t a, b, c, z, f_z;

    //Coefficients - > c, b, a
    a = coefficients[2];
    b = coefficients[1];
    c = coefficients[0];

#ifdef DEBUG
    cout << endl << "\tQuadratic equation:" << endl << "\t " << a << "*x^2 + " << b << "*x + " << c << " = 0 " << endl;
#endif

    if (a!=0){

        z = -b / (2 * a);

        if (z != std::numeric_limits<fp_t>::infinity())
        {
            fp_t fma_bzc = fma(b, z, c);
            f_z = fma(a, z * z, fma_bzc);

            std::complex<fp_t> sqrt_f = std::sqrt(std::complex<fp_t>((-f_z / a)));
            if (sqrt_f.real() != std::numeric_limits<fp_t>::infinity() && sqrt_f.imag()!= std::numeric_limits<fp_t>::infinity())
            {
                std::complex<fp_t> z_complex = std::complex<fp_t>(z);
                std::complex<fp_t> c_complex = std::complex<fp_t>(c);

                if (c > 0) { //уравнение имеет 2 одинаковых по знаку корня
                    if (b < 0) { // оба корня положительны (z > 0)

                        roots[0] = z_complex + sqrt_f;
#ifdef DEBUG
                        cout<< "\nroots[0] = " << roots[0];
#endif
                        roots[1] = c_complex/ roots[0];
#ifdef DEBUG
                        cout<< "\nroots[1] = " << roots[1];
#endif

                    }
                    else if (b > 0) {// оба корня отрицательны (z < 0)
                        roots[0] = z_complex - sqrt_f;
#ifdef DEBUG
                        cout<< "\nroots[0] = " << roots[0];
#endif
                        roots[1] = c_complex/ roots[0];
#ifdef DEBUG
                        cout<< "\nroots[1] = " << roots[1];
#endif
                    }
                }
                else if (c <= 0) { // 2 различных по знаку корня
                    if (b < 0) { // больший по модулю положителен (z > 0)
                        roots[0] = z_complex + sqrt_f;
#ifdef DEBUG
                        cout<< "\nroots[0] = " << roots[0];
#endif
                        roots[1] = c_complex/ roots[0];
#ifdef DEBUG
                        cout<< "\nroots[1] = " << roots[1];
#endif
                    }
                    else if (b > 0) { // больший по модулю отрицателен (z < 0)
                        roots[0] = z_complex - sqrt_f;
#ifdef DEBUG
                        cout<< "\nroots[0] = " << roots[0];
#endif
                        roots[1] = c_complex/ roots[0];
#ifdef DEBUG
                        cout<< "\nroots[1] = " << roots[1];
#endif
                    }
                }

            }

//#ifdef DEBUG
           // cout << "\nx1 = " << roots[0] << ", x2 = " << roots[1] << endl;
//#endif
        }

    } return roots; // возвращаем вектор комплексных корней
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

template<typename fp_t>
fp_t testPolynomial_complex(unsigned int roots_count) {
    fp_t max_absolute_error, max_relative_error;
    std::vector<std::complex<fp_t>> roots_computed(roots_count);//вектор компл. корней
    vector<fp_t> roots(roots_count), roots_computed_real(roots_count); //вектор вещ.корней
    vector<fp_t> coefficients(roots_count + 1);
    generate_polynomial<fp_t>(roots_count, 0, roots_count, 0, numeric_limits<fp_t>::min(), -1, 1, roots, coefficients);
    roots_computed = simple_formula_complex<fp_t>(coefficients, roots_computed); //находим корни - комплексный вектор

    if (!roots_computed.empty()) // если корни найдены
    {
        roots_computed_real = return_real_roots(roots_computed); //проверяем на комплексность - > отбрасываем
        // получаем вектор вещественных корней
        if (!roots_computed_real.empty()) { //если корни действительные - сравниваем

            compare_roots<fp_t>(roots_computed_real.size(), roots.size(), roots_computed_real, roots,
                                          max_absolute_error, max_relative_error);
            return max_absolute_error;
        }
        else return std::numeric_limits<fp_t>::infinity(); //если корни все же комплексные
    }
    else return std::numeric_limits<fp_t>::infinity();
}

int main() {
    // МЕТОД С ВЕЩЕСТВЕННЫМИ ВЫЧИСЛЕНИЯМИ
    float max_deviation = 0;
    for (auto i = 0; i < 100'000; ++i) {
        auto deviation = testPolynomial<float>(2);
        if (deviation != std::numeric_limits<float>::infinity()) {
            if (deviation > max_deviation) {
                max_deviation = deviation;
            }
        }
    }
    cout << endl << "MAX_deviation = " << max_deviation << endl;



    // МЕТОД С КОМПЛЕКСНЫМИ ВЫЧИСЛЕНИЯМИ
    float max_deviation_for_complex = 0;
    for (auto i = 0; i < 100'000; ++i) {

        auto deviation_for_complex = testPolynomial_complex<float>(2);
        if (deviation_for_complex != std::numeric_limits<float>::infinity()) {
            if (deviation_for_complex > max_deviation_for_complex) {
                max_deviation_for_complex = deviation_for_complex;
            }
        }

    }
    cout << endl << "MAX_deviation_for_complex = " << max_deviation_for_complex << endl;
}