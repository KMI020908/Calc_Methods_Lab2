#include<vector>
#include"methods.cpp"
#include"ioData.cpp"
#include"filePath.h"
#include<iomanip>

// Процедура проверки алгоритмов
template<typename Type>
void checkTest(std::vector<std::vector<Type>> &lCoefSys, std::vector<Type> &rCoefSys, const std::vector<Type> &startPoint,
const std::string &IN_FILE_PATH, const std::string &SIMPLE_IT_F_PATH, const std::string &JACOBI_F_PATH, const std::string &RELAXATION_F_PATH, 
Type accuracy = 1e-7, std::vector<double> realSolution = {}){
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
    Type bound = 0.0;
    Type errNorm = 0.0;
    ITERATION_METHOD_FLAG method;
    std::size_t exactIt = 0;
    
    // Метод простой итерации
    /*
    method = SIMPLE_IT;
    tao = 1e-4;
    numOfIt = simpleItMethod(lCoefSys, rCoefSys, startPoint, solution, tao, accuracy, p, eps0, 100000000);
    writeData(solution, startPoint, accuracy, SIMPLE_IT_F_PATH, numOfIt, tao, 0.0);
    bound = findLowerBoundOfIterations(lCoefSys, rCoefSys, startPoint, accuracy, method, tao, omega, p);
    writeBoundOfIterations(bound, SIMPLE_IT_F_PATH);
    writeNormOfError(solution, realSolution, SIMPLE_IT_F_PATH, p);
    exactIt = findExactItersSimpleItMethod(lCoefSys, rCoefSys, startPoint, solution, realSolution, tao, accuracy, p);
    writeExactIters(exactIt, SIMPLE_IT_F_PATH);
    errNorm = findNormOfErrAfterEstIt_SIT(lCoefSys, rCoefSys, startPoint, solution, realSolution, bound, tao, accuracy, p);
    writeNormErrAfterEstIt(errNorm, SIMPLE_IT_F_PATH);
    residual = findResidual(lCoefSys, rCoefSys, solution);
    writeResidual(residual, SIMPLE_IT_F_PATH);
    findCanonicalFormSimpleIt(lCoefSys, rCoefSys, C, y, tao);
    writeCanonicalForm(C, y, SIMPLE_IT_F_PATH);
    tao = 1e-6;
    numOfIt = simpleItMethod(lCoefSys, rCoefSys, startPoint, solution, tao, accuracy, p, eps0, 100000000);
    addData(solution, startPoint, accuracy, SIMPLE_IT_F_PATH, numOfIt, tao, 0.0);
    bound = findLowerBoundOfIterations(lCoefSys, rCoefSys, startPoint, accuracy, method, tao, omega, p);
    writeBoundOfIterations(bound, SIMPLE_IT_F_PATH);
    writeNormOfError(solution, realSolution, SIMPLE_IT_F_PATH, p);
    exactIt = findExactItersSimpleItMethod(lCoefSys, rCoefSys, startPoint, solution, realSolution, tao, accuracy, p);
    writeExactIters(exactIt, SIMPLE_IT_F_PATH);
    errNorm = findNormOfErrAfterEstIt_SIT(lCoefSys, rCoefSys, startPoint, solution, realSolution, bound, tao, accuracy, p);
    writeNormErrAfterEstIt(errNorm, SIMPLE_IT_F_PATH);
    residual = findResidual(lCoefSys, rCoefSys, solution);
    writeResidual(residual, SIMPLE_IT_F_PATH);
    findCanonicalFormSimpleIt(lCoefSys, rCoefSys, C, y, tao);
    writeCanonicalForm(C, y, SIMPLE_IT_F_PATH);
    */
    // Метод Якоби
    method = JACOBI;
    numOfIt = JacobiMethod(lCoefSys, rCoefSys, startPoint, solution, accuracy, p, eps0);
    writeData(solution, startPoint, accuracy, JACOBI_F_PATH, numOfIt);
    bound = findLowerBoundOfIterations(lCoefSys, rCoefSys, startPoint, accuracy, method, tao, omega, p);
    writeBoundOfIterations(bound, JACOBI_F_PATH);
    writeNormOfError(solution, realSolution, JACOBI_F_PATH, p);
    exactIt = findExactItersJacobiMethod(lCoefSys, rCoefSys, startPoint, solution, realSolution, accuracy, p);
    writeExactIters(exactIt, JACOBI_F_PATH);
    errNorm = findNormOfErrAfterEstIt_JAC(lCoefSys, rCoefSys, startPoint, solution, realSolution, bound, accuracy, p);
    writeNormErrAfterEstIt(errNorm, JACOBI_F_PATH);
    residual = findResidual(lCoefSys, rCoefSys, solution);
    writeResidual(residual, JACOBI_F_PATH);
    findCanonicalFormJacobi(lCoefSys, rCoefSys, C, y);
    writeCanonicalForm(C, y, JACOBI_F_PATH);

    // Метод Релаксации (Зейделя в случае омега равной 1)
    method = RELAXATION;
    omega = 1.0;
    numOfIt = relaxationMethod(lCoefSys, rCoefSys, startPoint, solution, accuracy, omega, p, eps0);
    writeData(solution, startPoint, accuracy, RELAXATION_F_PATH, numOfIt, 0.0, omega);
    bound = findLowerBoundOfIterations(lCoefSys, rCoefSys, startPoint, accuracy, method, tao, omega, p);
    writeBoundOfIterations(bound, RELAXATION_F_PATH);
    writeNormOfError(solution, realSolution, RELAXATION_F_PATH, p);
    exactIt = findExactRelaxationMethod(lCoefSys, rCoefSys, startPoint, solution, realSolution, accuracy, omega, p);
    writeExactIters(exactIt, RELAXATION_F_PATH);
    errNorm = findNormOfErrAfterEstIt_REL(lCoefSys, rCoefSys, startPoint, solution, realSolution, bound, omega, accuracy, p);
    writeNormErrAfterEstIt(errNorm, RELAXATION_F_PATH);
    residual = findResidual(lCoefSys, rCoefSys, solution);
    writeResidual(residual, RELAXATION_F_PATH);
    findCanonicalFormRelaxation(lCoefSys, rCoefSys, C, y, omega);
    writeCanonicalForm(C, y, RELAXATION_F_PATH);
    omega = 0.25;
    numOfIt = relaxationMethod(lCoefSys, rCoefSys, startPoint, solution, accuracy, omega, p, eps0);
    addData(solution, startPoint, accuracy, RELAXATION_F_PATH, numOfIt, 0.0, omega);
    bound = findLowerBoundOfIterations(lCoefSys, rCoefSys, startPoint, accuracy, method, tao, omega, p);
    writeBoundOfIterations(bound, RELAXATION_F_PATH);
    writeNormOfError(solution, realSolution, RELAXATION_F_PATH, p);
    exactIt = findExactRelaxationMethod(lCoefSys, rCoefSys, startPoint, solution, realSolution, accuracy, omega, p);
    writeExactIters(exactIt, RELAXATION_F_PATH);
    errNorm = findNormOfErrAfterEstIt_REL(lCoefSys, rCoefSys, startPoint, solution, realSolution, bound, omega, accuracy, p);
    writeNormErrAfterEstIt(errNorm, RELAXATION_F_PATH);
    residual = findResidual(lCoefSys, rCoefSys, solution);
    writeResidual(residual, RELAXATION_F_PATH);
    findCanonicalFormRelaxation(lCoefSys, rCoefSys, C, y, omega);
    writeCanonicalForm(C, y, RELAXATION_F_PATH);
    omega = 0.5;
    numOfIt = relaxationMethod(lCoefSys, rCoefSys, startPoint, solution, accuracy, omega, p, eps0);
    addData(solution, startPoint, accuracy, RELAXATION_F_PATH, numOfIt, 0.0, omega);
    bound = findLowerBoundOfIterations(lCoefSys, rCoefSys, startPoint, accuracy, method, tao, omega, p);
    writeBoundOfIterations(bound, RELAXATION_F_PATH);
    writeNormOfError(solution, realSolution, RELAXATION_F_PATH, p);
    exactIt = findExactRelaxationMethod(lCoefSys, rCoefSys, startPoint, solution, realSolution, accuracy, omega, p);
    writeExactIters(exactIt, RELAXATION_F_PATH);
    errNorm = findNormOfErrAfterEstIt_REL(lCoefSys, rCoefSys, startPoint, solution, realSolution, bound, omega, accuracy, p);
    writeNormErrAfterEstIt(errNorm, RELAXATION_F_PATH);
    residual = findResidual(lCoefSys, rCoefSys, solution);
    writeResidual(residual, RELAXATION_F_PATH);
    findCanonicalFormRelaxation(lCoefSys, rCoefSys, C, y, omega);
    writeCanonicalForm(C, y, RELAXATION_F_PATH);
    omega = 1.25;
    numOfIt = relaxationMethod(lCoefSys, rCoefSys, startPoint, solution, accuracy, omega, p, eps0);
    addData(solution, startPoint, accuracy, RELAXATION_F_PATH, numOfIt, 0.0, omega);
    bound = findLowerBoundOfIterations(lCoefSys, rCoefSys, startPoint, accuracy, method, tao, omega, p);
    writeBoundOfIterations(bound, RELAXATION_F_PATH);
    writeNormOfError(solution, realSolution, RELAXATION_F_PATH, p);
    exactIt = findExactRelaxationMethod(lCoefSys, rCoefSys, startPoint, solution, realSolution, accuracy, omega, p);
    writeExactIters(exactIt, RELAXATION_F_PATH);
    errNorm = findNormOfErrAfterEstIt_REL(lCoefSys, rCoefSys, startPoint, solution, realSolution, bound, omega, accuracy, p);
    writeNormErrAfterEstIt(errNorm, RELAXATION_F_PATH);
    residual = findResidual(lCoefSys, rCoefSys, solution);
    writeResidual(residual, RELAXATION_F_PATH);
    findCanonicalFormRelaxation(lCoefSys, rCoefSys, C, y, omega);
    writeCanonicalForm(C, y, RELAXATION_F_PATH);
    omega = 1.5;
    numOfIt = relaxationMethod(lCoefSys, rCoefSys, startPoint, solution, accuracy, omega, p, eps0);
    addData(solution, startPoint, accuracy, RELAXATION_F_PATH, numOfIt, 0.0, omega);
    bound = findLowerBoundOfIterations(lCoefSys, rCoefSys, startPoint, accuracy, method, tao, omega, p);
    writeBoundOfIterations(bound, RELAXATION_F_PATH);
    writeNormOfError(solution, realSolution, RELAXATION_F_PATH, p);
    residual = findResidual(lCoefSys, rCoefSys, solution);
    writeResidual(residual, RELAXATION_F_PATH);
    findCanonicalFormRelaxation(lCoefSys, rCoefSys, C, y, omega);
    writeCanonicalForm(C, y, RELAXATION_F_PATH);
    omega = 1.75;
    numOfIt = relaxationMethod(lCoefSys, rCoefSys, startPoint, solution, accuracy, omega, p, eps0);
    addData(solution, startPoint, accuracy, RELAXATION_F_PATH, numOfIt, 0.0, omega);
    bound = findLowerBoundOfIterations(lCoefSys, rCoefSys, startPoint, accuracy, method, tao, omega, p);
    writeBoundOfIterations(bound, RELAXATION_F_PATH);
    writeNormOfError(solution, realSolution, RELAXATION_F_PATH, p);
    residual = findResidual(lCoefSys, rCoefSys, solution);
    writeResidual(residual, RELAXATION_F_PATH);
    findCanonicalFormRelaxation(lCoefSys, rCoefSys, C, y, omega);
    writeCanonicalForm(C, y, RELAXATION_F_PATH);

}

template<typename Type>
void temp_main(){

    std::vector<std::vector<Type>> lCoefSys; // Матрица левых коэффициентов
    std::vector<Type> rCoefSys; // Вектор правых коэффициентов
    std::vector<Type> startPoint = {0.0, 0.0, 0.0, 0.0};

    std::vector<double> realSolution2 = {10.0, -10.0, 12.0, 4.0};
    std::vector<double> realSolution1 = {5.0, -7.0, 12.0, 4.0};

    // Точность 1e-4
    //checkTest(lCoefSys, rCoefSys, startPoint, IN_FILE_PATH_1, SIMPLE_IT_F_PATH_1_EPS4, JACOBI_F_PATH_1_EPS4, RELAXATION_F_PATH_1_EPS4, 1e-4);

    //checkTest(lCoefSys, rCoefSys, startPoint, IN_FILE_PATH_2, SIMPLE_IT_F_PATH_2_EPS4, JACOBI_F_PATH_2_EPS4, RELAXATION_F_PATH_2_EPS4, 1e-4);

    checkTest(lCoefSys, rCoefSys, startPoint, IN_FILE_PATH_4, SIMPLE_IT_F_PATH_4_EPS4, JACOBI_F_PATH_4_EPS4, RELAXATION_F_PATH_4_EPS4, 1e-4, realSolution1);

    checkTest(lCoefSys, rCoefSys, startPoint, IN_FILE_PATH_5, SIMPLE_IT_F_PATH_5_EPS4, JACOBI_F_PATH_5_EPS4, RELAXATION_F_PATH_5_EPS4, 1e-4, realSolution2);

    // Точность 1e-7
    //checkTest(lCoefSys, rCoefSys, startPoint, IN_FILE_PATH_1, SIMPLE_IT_F_PATH_1_EPS7, JACOBI_F_PATH_1_EPS7, RELAXATION_F_PATH_1_EPS7, 1e-7);

    //checkTest(lCoefSys, rCoefSys, startPoint, IN_FILE_PATH_2, SIMPLE_IT_F_PATH_2_EPS7, JACOBI_F_PATH_2_EPS7, RELAXATION_F_PATH_2_EPS7, 1e-7);

    checkTest(lCoefSys, rCoefSys, startPoint, IN_FILE_PATH_4, SIMPLE_IT_F_PATH_4_EPS7, JACOBI_F_PATH_4_EPS7, RELAXATION_F_PATH_4_EPS7, 1e-7, realSolution1);

    checkTest(lCoefSys, rCoefSys, startPoint, IN_FILE_PATH_5, SIMPLE_IT_F_PATH_5_EPS7, JACOBI_F_PATH_5_EPS7, RELAXATION_F_PATH_5_EPS7, 1e-7, realSolution2);

}

int main(){
    //temp_main<double>();

    std::vector<std::vector<double>> A;
    std::vector<double> f;
    std::vector<double> solution;
    readData(A, f, IN_FILE_PATH_3);
    std::vector<double> fVec = {40.0, 5.0};
    writePointsOfJacobiMethod(A, f, fVec, solution, JACOBI_POINTS_FILE_PATH);
    writePointsOfRelaxationMethod(A, f, fVec, solution, ZEIDEL_POINTS_FILE_PATH);
    writePointsOfRelaxationMethod(A, f, fVec, solution, RELAXATION_POINTS_FILE_PATH, 1e-7, 1.0);
    std::vector<std::vector<double>> C;
    std::vector<double> y;
    findCanonicalFormRelaxation(A, f, C, y);
    std::cout << C;
    std::cout << '\n';
    std::cout << normOfMatrix(C, INFINITY);
    std::cout << '\n' << '\n' << '\n';

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
    //std::cout << relaxationMethodFor3Diag(a, b, c, d, firstVec, solution, accuracy, 0.25, INFINITY, 1e-4) << '\n';
    //std::cout << solution << '\n' << '\n';
    for (size_t i = 0; i < n - 1; i++){
        double temp = a[i];
        a[i] = c[i];
        c[i] = temp;
    }
    //std::cout << relaxationMethodFor3Diag(a, b, c, d, firstVec, solution, accuracy, 0.25, INFINITY, 1e-4) << '\n';
    //std::cout << solution << '\n' << '\n';

    a.clear();
    b.clear();
    c.clear();
    d.clear();
    for (size_t i = 0; i < n - 1; i++){
        a.push_back(6.0);
        c.push_back(1.0);
    }
    for (size_t i = 0; i < n; i++){
        b.push_back(8.0);
        d.push_back(1.0);
    }
    std::cout << relaxationMethodFor3Diag(a, b, c, d, firstVec, solution, accuracy, 1.0, INFINITY, 1e-7) << '\n';
    std::cout << solution << '\n' << '\n';
    std::cout << relaxationMethodFor3Diag(c, b, a, d, firstVec, solution, accuracy, 1.0, INFINITY, 1e-7) << '\n';
    std::cout << solution << '\n' << '\n';

    std::size_t dim = 200;
    std::vector<std::vector<double>> matrix(dim);
    for (size_t i = 0; i < dim; i++){
        matrix[i].resize(dim, 0.0);
    }
    matrix[0][0] = 8.0;
    matrix[0][1] = 1.0;
    matrix[dim - 1][dim - 1] = 8.0;
    matrix[dim - 1][dim - 2] = 6.0;
    for (size_t i = 1; i < dim - 1; i++){
        matrix[i][i] = 8.0;
        matrix[i][i + 1] = 1.0;
        matrix[i][i - 1] = 6.0;
    }

    
    std::vector<std::vector<double>> matrix2(dim);
    for (size_t i = 0; i < dim; i++){
        matrix2[i].resize(dim, 0.0);
    }
    matrix2[0][0] = 8.0;
    matrix2[0][1] = 6.0;
    matrix2[dim - 1][dim - 1] = 8.0;
    matrix2[dim - 1][dim - 2] = 1.0;
    for (size_t i = 1; i < dim - 1; i++){
        matrix2[i][i] = 8.0;
        matrix2[i][i + 1] = 6.0;
        matrix2[i][i - 1] = 1.0;
    }
    
    std::vector<double> l(dim);
    for (size_t i = 0; i < dim; i++){
        l[i] = 1.0;
    }

    findCanonicalFormRelaxation(matrix, l, C, firstVec, 1.0);
    double otv = normOfMatrix(C, INFINITY);
    findCanonicalFormRelaxation(matrix2, l, C, firstVec, 1.0);
    double otv2 = normOfMatrix(C, INFINITY);
    return 0;
}