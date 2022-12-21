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

// Определение знака числа
// разумно сделать inline(встроенную функцию), т к ф-ия небольшая, вызывается часто
template <typename fp_t>
inline int sgn(fp_t val) {
    return (fp_t(0) < val) - (val < fp_t(0));
}

//A simple formula for solving quadratic equations using function evaluation
// МЕТОД С ВЕЩЕСТВЕННЫМИ ВЫЧИСЛЕНИЯМИ
template<typename fp_t>
int simple_formula(fp_t a, fp_t b, fp_t c, std::vector<fp_t> &roots) {
    fp_t  z, f_z;
    int cnt_real_roots; // число вещественных корней
#ifdef DEBUG
    cout << endl << "\tQuadratic equation:" << endl << "\t " << a << "*x^2 + " << b << "*x + " << c << " = 0 " << endl;
#endif
    // не проверяем на а!=0 , т. к. далее проверяем z на бесконечность
        z = -b / (2*a);
        if (!isinf(z)) {
            //Формула: f_z = a * pow(z, 2) + b * z + c;
            // То есть f_z  = a*b^2/4a^2 + b*z + c = (-b/2a)*(-b/2) + b*z + c = (-b/2)*z + b*z + c = (b/2)*z + c
            f_z = fma(b , z , 2*c);
            // f_z = -(b^2 - 4ac)/4a = -D/4a - - - > D>=0(f_z<=0) - real, D<0(f_z>0) - complex
            if (f_z <= 0) {
                fp_t sqrt_f = sqrt(-f_z / (2*a)); // sqrt((b^2-4ac)/4a^2) - "половина расстояния между корнями"
                if (!isinf(sqrt_f)) {
                    cnt_real_roots = 2; //число действительных корней = 2
                        // b < 0 ---> z > 0 => z + sqrt_f
                        // b > 0 ---> z < 0 => z - sqrt_f
                        roots[0] = z - sgn(b)*sqrt_f;
                        roots[1] = c/roots[0];
                    if(isinf(roots[1])) { // если 2-ой корень бесконечный - его откидываем, другой "спасаем"
                        cnt_real_roots = 1;
                        roots[1] = 0;
                    }
                }
                else return -1; // корни находятся очень далеко друг от друга - возвращаем -1(в testPolynomial идет проверка)

#ifdef DEBUG
                cout << "x1 = " << roots[0] << ", x2 = " << roots[1] << endl;
#endif
            }
            else {
                std::complex<fp_t> sqrt_f_z = std::sqrt(std::complex<fp_t>((-f_z / (2*a))));
                if (!isinf(sqrt_f_z.real())  && !isinf(sqrt_f_z.imag())){ // проверка на inf
                    std::complex<fp_t> root1(z, sqrt_f_z.imag());
                    std::complex<fp_t> root2(z, -sqrt_f_z.imag());
                    // выясним, действительно ли корень комплексный:
                    // если std::abs(root1.imag())> std::abs(root1)*std::numeric_limits<fp_t>::epsilon() - > корень комплексный, иначе - действительный и полученная
                    // комплексная часть - мусор...
                    if (std::abs(root1.imag()) > std::abs(root1) * std::numeric_limits<fp_t>::epsilon()) {
                        //значит корень комплексный
                        cnt_real_roots = 0; // число действительных корней = 0
                    }
                    else { // корень действительный
#ifdef DEBUG
                        cout << "\t\tNot complex roots! Real roots: " << endl;
#endif
                            roots[0] = roots[1] = root1.real();
                            cnt_real_roots = 2; // число действительных корней = 2
                    }
                }
                else return -1;  // корни находятся очень далеко друг от друга - возвращаем -1(в testPolynomial идет проверка)
            }
        }
        else { // b/a = inf - значит уравнение формально не квадратное, находим единственный корень
                roots[0] = -c/b; // не проверяем b!=0, т.к. делаем проверку корня на бесконечность
                if(!isinf(roots[0])) cnt_real_roots = 1;
                else cnt_real_roots = 0;
        }
    return cnt_real_roots;
}

// МЕТОД С КОМПЛЕКСНЫМИ ВЫЧИСЛЕНИЯМИ
template<typename fp_t>
std::vector<std::complex<fp_t>>  simple_formula_complex(fp_t a, fp_t b, fp_t c, std::vector<std::complex<fp_t>> &roots)
{
    fp_t z, f_z;

#ifdef DEBUG
    cout << endl << "\tQuadratic equation:" << endl << "\t " << a << "*x^2 + " << b << "*x + " << c << " = 0 " << endl;
#endif
    // тут не проверяем на а!=0 , т к далее проверка z на inf
        z = -b / (2 * a);
        if (!isinf(z)){
            f_z = fma(b , z , 2*c);
            std::complex<fp_t> sqrt_f = std::sqrt(std::complex<fp_t>((-f_z / (2*a)))); // тут впервые может появиться комплексность
            if (!isinf(sqrt_f.real()) && !isinf(sqrt_f.imag())) { // проверка на inf
                // b < 0 ---> z > 0 => z + sqrt_f
                // b > 0 ---> z < 0 => z - sqrt_f
                roots[0] = std::complex<fp_t>(z - sgn(b)*sqrt_f.real(), sqrt_f.imag());
                roots[1] = std::complex<fp_t>(c/roots[0]);
                std::complex<fp_t> tmp = roots[1];
                if(isinf(tmp.real()) || isinf(tmp.imag())) // если 2-ой корень бесконечный - его откидываем, другой "спасаем"
                    roots.pop_back();
            }
            else return roots; // если inf - просто возвращаем пустой вектор (в testPolynomial_complex проверка на пустоту)

#ifdef DEBUG
           // cout << "\nx1 = " << roots[0] << ", x2 = " << roots[1] << endl;
#endif
        }
        else { // b/a = inf -> уравнение линейное - находим корень
            roots[0] = -c/b;
            if (!isinf(roots[0].real())) return roots;
            else roots.clear(); //если корень "бесконечный", то очищаем вектор
        }
    return roots; // возвращаем вектор комплексных корней
}

// Для вещественной арифметики
template<typename fp_t>
auto testPolynomial(unsigned int roots_count) {
    fp_t max_absolute_error, max_relative_error; // макс.абсолют. и относит. погрешности
    vector<fp_t> roots_computed(roots_count); // найденные корни методом
    vector<fp_t> roots(roots_count), coefficients(roots_count + 1); // истинные корни; коэффициенты
    generate_polynomial<fp_t>(roots_count, 0, roots_count, 0, MAX_DISTANCE, -1, 1, roots,
                              coefficients);
    int cnt_real_roots = simple_formula(coefficients[2], coefficients[1],coefficients[0],  roots_computed); // находим число вещ.корней и сами корни
    if (cnt_real_roots != 0 && cnt_real_roots != -1) { // если корни найдены и мы не попали в какой-то исключ.случай,
                                                        // то сравниваем найденные корни с истинными
         compare_roots<fp_t>(roots_computed.size(), roots.size(), roots_computed, roots,
                                          max_absolute_error, max_relative_error);
         
    } else max_absolute_error = 0, max_relative_error = 0;
    return pair<fp_t, fp_t>(max_absolute_error, max_relative_error);
}

// Для комплексной арифметики
template<typename fp_t>
auto testPolynomial_complex(unsigned int roots_count) {
    fp_t max_absolute_error, max_relative_error; // максимальные абсолют. и относит. погрешности
    std::vector<std::complex<fp_t>> roots_computed(roots_count);//вектор комплексных корней
    vector<fp_t> roots(roots_count);//вектор истинных корней
    vector<fp_t> coefficients(roots_count + 1); //коэффициенты
    generate_polynomial<fp_t>(roots_count, 0, roots_count, 0, MAX_DISTANCE, -1, 1, roots, coefficients);
    simple_formula_complex<fp_t>(coefficients[2], coefficients[1] , coefficients[0],roots_computed); //находим корни - комплексный вектор
    //compare_roots<fp_t>(roots_computed.size(), roots.size(), roots_computed, roots, max_absolute_error, max_relative_error);
    // отправляем на сравнение в ф-ию compare_roots_complex, где происходит отбрасывание истинных комплексных корней
    // и вычисление абсолютной и относительной погрешности
    compare_roots_complex<fp_t>(roots_computed.size(), roots.size(), roots_computed, roots,max_absolute_error, max_relative_error);

    return pair<fp_t, fp_t>(max_absolute_error, max_relative_error);
}

int main() {
    // МЕТОД С ВЕЩЕСТВЕННЫМИ ВЫЧИСЛЕНИЯМИ
    float max_absolut_deviation = 0;
    float max_relative_deviation = 0;
    for (auto i = 0; i < 10'000'000; ++i) {
        auto deviation = testPolynomial<float>(2);

            if (deviation.first > max_absolut_deviation) {
                max_absolut_deviation = deviation.first;
            }
            if (deviation.second > max_relative_deviation) {
                max_relative_deviation = deviation.second;
            }
    }
    cout << endl << "MAX_ABSOLUT_deviation = " << max_absolut_deviation << endl;
    cout << endl << "MAX_RELATIVE_deviation = " << max_relative_deviation << endl;


    // МЕТОД С КОМПЛЕКСНЫМИ ВЫЧИСЛЕНИЯМИ
    float max_absolut_deviation_for_complex = 0;
    float max_relative_deviation_for_complex = 0;
    for (auto i = 0; i < 10'000'000; ++i) {

        auto deviation_for_complex = testPolynomial_complex<float>(2);

            if (deviation_for_complex.first > max_absolut_deviation_for_complex) {
                max_absolut_deviation_for_complex = deviation_for_complex.first;

            }
            if (deviation_for_complex.second > max_relative_deviation_for_complex) {
                max_relative_deviation_for_complex = deviation_for_complex.second;
            }
    }
    cout << endl << "MAX_ABSOLUT_deviation_for_complex = " << max_absolut_deviation_for_complex << endl;
    cout << endl << "MAX_RELATIVE_deviation_for_complex = " << max_relative_deviation_for_complex << endl;
}