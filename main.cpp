#include<vector>
#include"methods.cpp"
#include"ioData.cpp"
#include"filePath.h"
#include<iomanip>

// Процедура проверки алгоритмов
template<typename Type>
void checkTest(std::vector<std::vector<Type>> &lCoefSys, std::vector<Type> &rCoefSys,
const std::string &IN_FILE_PATH, const std::string &SIMPLE_IT_F_PATH, const std::string &JACOBI_F_PATH, const std::string &SEIDEL_F_PATH, 
Type accuracy = 1e-7, Type perturbation = 0.01){
    readData<Type>(lCoefSys, rCoefSys, IN_FILE_PATH);
    std::vector<Type> solution;
    
}

template<typename Type>
void temp_main(){

    std::vector<std::vector<Type>> lCoefSys; // Матрица левых коэффициентов
    std::vector<Type> rCoefSys; // Вектор правых коэффициентов

    // Точность 1e-4
    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_1, SIMPLE_IT_F_PATH_1_EPS4, JACOBI_F_PATH_1_EPS4, RELAXATION_F_PATH_1_EPS4, 1e-4);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_2, SIMPLE_IT_F_PATH_2_EPS4, JACOBI_F_PATH_2_EPS4, RELAXATION_F_PATH_2_EPS4, 1e-4);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_3, SIMPLE_IT_F_PATH_3_EPS4, JACOBI_F_PATH_3_EPS4, RELAXATION_F_PATH_3_EPS4, 1e-4);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_4, SIMPLE_IT_F_PATH_4_EPS4, JACOBI_F_PATH_4_EPS4, RELAXATION_F_PATH_4_EPS4, 1e-4);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_5, SIMPLE_IT_F_PATH_5_EPS4, JACOBI_F_PATH_5_EPS4, RELAXATION_F_PATH_5_EPS4, 1e-4);

    generateRandomTest<Type>(4, 1.0, 20.0, IN_FILE_PATH_6);
    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_6, SIMPLE_IT_F_PATH_6_EPS4, JACOBI_F_PATH_6_EPS4, RELAXATION_F_PATH_6_EPS4, 1e-4);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_7, SIMPLE_IT_F_PATH_7_EPS4, JACOBI_F_PATH_7_EPS4, RELAXATION_F_PATH_7_EPS4, 1e-4);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_8, SIMPLE_IT_F_PATH_8_EPS4, JACOBI_F_PATH_8_EPS4, RELAXATION_F_PATH_8_EPS4, 1e-4);

    // Точность 1e-7
    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_1, SIMPLE_IT_F_PATH_1_EPS7, JACOBI_F_PATH_1_EPS7, RELAXATION_F_PATH_1_EPS7, 1e-7);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_2, SIMPLE_IT_F_PATH_2_EPS7, JACOBI_F_PATH_2_EPS7, RELAXATION_F_PATH_2_EPS7, 1e-7);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_3, SIMPLE_IT_F_PATH_3_EPS7, JACOBI_F_PATH_3_EPS7, RELAXATION_F_PATH_3_EPS7, 1e-7);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_4, SIMPLE_IT_F_PATH_4_EPS7, JACOBI_F_PATH_4_EPS7, RELAXATION_F_PATH_4_EPS7, 1e-7);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_5, SIMPLE_IT_F_PATH_5_EPS7, JACOBI_F_PATH_5_EPS7, RELAXATION_F_PATH_5_EPS7, 1e-7);

    generateRandomTest<Type>(4, 1.0, 20.0, IN_FILE_PATH_6);
    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_6, SIMPLE_IT_F_PATH_6_EPS7, JACOBI_F_PATH_6_EPS7, RELAXATION_F_PATH_6_EPS7, 1e-7);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_7, SIMPLE_IT_F_PATH_7_EPS7, JACOBI_F_PATH_7_EPS7, RELAXATION_F_PATH_7_EPS7, 1e-7);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_8, SIMPLE_IT_F_PATH_8_EPS7, JACOBI_F_PATH_8_EPS7, RELAXATION_F_PATH_8_EPS7, 1e-7);
}

int main(){
    //temp_main<double>();
    size_t count = 5;
    std::vector<std::vector<double>> lCoefs;
    std::vector<double> rCoefs;
    readData(lCoefs, rCoefs, IN_FILE_PATH_1);
    std::size_t rows, cols = lCoefs.size();
    std::vector<double> solution;
    std::vector<double> firstVec(rCoefs.size(), 0);
    double accuracy = 1e-7;

    // Большая матрица
    std::vector<double> a, b, c, d;
    std::size_t N = 9;
    std::size_t n = 200 + N;
    for (std::size_t i = 0; i < n - 1; i++){
        a.push_back(1);
        c.push_back(1);
    }
    for (std::size_t i = 0; i < n; i++){
        b.push_back(4);
    }
    d.push_back(6);
    for (std::size_t i = 1; i < n - 1; i++){
        d.push_back(10 - 2 * (i % 2));
    }
    d.push_back(9 - 3 * ((n - 1) % 2));
    firstVec.resize(n);
    relaxationMethodFor3Diag(a, b, c, d, firstVec, solution, accuracy, 1.0, INFINITY, 1e-1);
    std::cout << std::setprecision(std::abs(std::log10(accuracy))) << solution << '\n';

    std::vector<std::vector<double>> C;
    std::vector<double> y;
    findCanonicalFormJacobi(lCoefs, rCoefs, C, y);
    std::cout << normOfMatrix(C, INFINITY) << '\n';
    simpleItMethod(lCoefs, rCoefs, firstVec, solution, 0.001, accuracy, INFINITY);
    std::cout << std::setprecision(std::abs(std::log10(accuracy))) << solution << '\n';
    JacobiMethod(lCoefs, rCoefs, firstVec, solution, accuracy, INFINITY);
    std::cout << std::setprecision(std::abs(std::log10(accuracy))) << solution << '\n';
    relaxationMethod(lCoefs, rCoefs, firstVec, solution, accuracy, 1.0, INFINITY);
    std::cout << std::setprecision(std::abs(std::log10(accuracy))) << solution << '\n';

    std::cout << C;
    
    return 0;
}