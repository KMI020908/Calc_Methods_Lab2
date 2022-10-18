#include<vector>
#include"methods.cpp"
#include"ioData.cpp"
#include"filePath.h"
#include<iomanip>

// Процедура проверки алгоритмов
template<typename Type>
void checkTest(std::vector<std::vector<Type>> &lCoefSys, std::vector<Type> &rCoefSys, const std::vector<Type> &startPoint,
const std::string &IN_FILE_PATH, const std::string &SIMPLE_IT_F_PATH, const std::string &JACOBI_F_PATH, const std::string &RELAXATION_F_PATH, 
Type accuracy = 1e-7){
    // Считываем данные
    readData<Type>(lCoefSys, rCoefSys, IN_FILE_PATH);

    // Нужные параметры и переменные
    std::vector<Type> solution;
    std::vector<std::vector<Type>> C;
    std::vector<Type> y;
    Type residual = 0.0;
    std::size_t numOfIt = 0;
    Type tao = 0.0;
    Type omega = 0.0;
    double p = INFINITY;
    Type eps0 = 1e-8;

    // Метод простой итерации
    tao = 1e-4;
    numOfIt = simpleItMethod(lCoefSys, rCoefSys, startPoint, solution, tao, accuracy, p, eps0, 10000000);
    writeData(solution, startPoint, accuracy, SIMPLE_IT_F_PATH, numOfIt, tao, 0.0);
    residual = findResidual(lCoefSys, rCoefSys, solution);
    writeResidual(residual, SIMPLE_IT_F_PATH);
    findCanonicalFormSimpleIt(lCoefSys, rCoefSys, C, y, tao);
    writeCanonicalForm(C, y, SIMPLE_IT_F_PATH);

    // Метод Якоби
    numOfIt = JacobiMethod(lCoefSys, rCoefSys, startPoint, solution, accuracy, p, eps0);
    writeData(solution, startPoint, accuracy, JACOBI_F_PATH, numOfIt);
    residual = findResidual(lCoefSys, rCoefSys, solution);
    writeResidual(residual, JACOBI_F_PATH);
    findCanonicalFormJacobi(lCoefSys, rCoefSys, C, y);
    writeCanonicalForm(C, y, JACOBI_F_PATH);

    // Метод Релаксации (Зейделя в случае омега равной 1)
    omega = 1.0;
    numOfIt = relaxationMethod(lCoefSys, rCoefSys, startPoint, solution, accuracy, omega, p, eps0);
    writeData(solution, startPoint, accuracy, RELAXATION_F_PATH, numOfIt, 0.0, omega);
    residual = findResidual(lCoefSys, rCoefSys, solution);
    writeResidual(residual, RELAXATION_F_PATH);
    omega = 1.25;
    numOfIt = relaxationMethod(lCoefSys, rCoefSys, startPoint, solution, accuracy, omega, p, eps0);
    addData(solution, startPoint, accuracy, RELAXATION_F_PATH, numOfIt, 0.0, omega);
    residual = findResidual(lCoefSys, rCoefSys, solution);
    writeResidual(residual, RELAXATION_F_PATH);
    //findCanonicalFormSimpleIt(lCoefSys, rCoefSys, C, y, tao);
    //writeCanonicalForm(C, y, SIMPLE_IT_F_PATH);

}

template<typename Type>
void temp_main(){

    std::vector<std::vector<Type>> lCoefSys; // Матрица левых коэффициентов
    std::vector<Type> rCoefSys; // Вектор правых коэффициентов
    std::vector<Type> startPoint = {0.0, 0.0, 0.0, 0.0};

    // Точность 1e-4
    checkTest(lCoefSys, rCoefSys, startPoint, IN_FILE_PATH_1, SIMPLE_IT_F_PATH_1_EPS4, JACOBI_F_PATH_1_EPS4, RELAXATION_F_PATH_1_EPS4, 1e-4);

    checkTest(lCoefSys, rCoefSys, startPoint, IN_FILE_PATH_2, SIMPLE_IT_F_PATH_2_EPS4, JACOBI_F_PATH_2_EPS4, RELAXATION_F_PATH_2_EPS4, 1e-4);

    checkTest(lCoefSys, rCoefSys, startPoint, IN_FILE_PATH_4, SIMPLE_IT_F_PATH_4_EPS4, JACOBI_F_PATH_4_EPS4, RELAXATION_F_PATH_4_EPS4, 1e-4);

    checkTest(lCoefSys, rCoefSys, startPoint, IN_FILE_PATH_5, SIMPLE_IT_F_PATH_5_EPS4, JACOBI_F_PATH_5_EPS4, RELAXATION_F_PATH_5_EPS4, 1e-4);

    // Точность 1e-7
    checkTest(lCoefSys, rCoefSys, startPoint, IN_FILE_PATH_1, SIMPLE_IT_F_PATH_1_EPS7, JACOBI_F_PATH_1_EPS7, RELAXATION_F_PATH_1_EPS7, 1e-7);

    checkTest(lCoefSys, rCoefSys, startPoint, IN_FILE_PATH_2, SIMPLE_IT_F_PATH_2_EPS7, JACOBI_F_PATH_2_EPS7, RELAXATION_F_PATH_2_EPS7, 1e-7);

    checkTest(lCoefSys, rCoefSys, startPoint, IN_FILE_PATH_4, SIMPLE_IT_F_PATH_4_EPS7, JACOBI_F_PATH_4_EPS7, RELAXATION_F_PATH_4_EPS7, 1e-7);

    checkTest(lCoefSys, rCoefSys, startPoint, IN_FILE_PATH_5, SIMPLE_IT_F_PATH_5_EPS7, JACOBI_F_PATH_5_EPS7, RELAXATION_F_PATH_5_EPS7, 1e-7);

}

int main(){
    temp_main<double>();

    std::vector<std::vector<double>> A;
    std::vector<double> f;
    std::vector<double> solution;
    readData(A, f, IN_FILE_PATH_3);
    std::vector<double> fVec = {3, 0.25};
    writePointsOfJacobiMethod(A, f, fVec, solution, JACOBI_POINTS_FILE_PATH);
    writePointsOfRelaxationMethod(A, f, fVec, solution, ZEIDEL_POINTS_FILE_PATH);
    writePointsOfRelaxationMethod(A, f, fVec, solution, RELAXATION_POINTS_FILE_PATH, 1e-7, 1.5);

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
    std::vector<double> firstVec(n, 0.0);
    double accuracy = 1e-7;
    relaxationMethodFor3Diag(a, b, c, d, firstVec, solution, accuracy, 1.0, INFINITY, 1e-4);
    std::cout << solution << '\n';
    
    return 0;
}