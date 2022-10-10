#include<vector>
#include"methods.cpp"
#include"ioData.cpp"
#include"filePath.h"

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
    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_1, SIMPLE_IT_F_PATH_1_EPS4, JACOBI_F_PATH_1_EPS4, SEIDEL_F_PATH_1_EPS4, 1e-4);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_2, SIMPLE_IT_F_PATH_2_EPS4, JACOBI_F_PATH_2_EPS4, SEIDEL_F_PATH_2_EPS4, 1e-4);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_3, SIMPLE_IT_F_PATH_3_EPS4, JACOBI_F_PATH_3_EPS4, SEIDEL_F_PATH_3_EPS4, 1e-4);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_4, SIMPLE_IT_F_PATH_4_EPS4, JACOBI_F_PATH_4_EPS4, SEIDEL_F_PATH_4_EPS4, 1e-4);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_5, SIMPLE_IT_F_PATH_5_EPS4, JACOBI_F_PATH_5_EPS4, SEIDEL_F_PATH_5_EPS4, 1e-4);

    generateRandomTest<Type>(4, 1.0, 20.0, IN_FILE_PATH_6);
    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_6, SIMPLE_IT_F_PATH_6_EPS4, JACOBI_F_PATH_6_EPS4, SEIDEL_F_PATH_6_EPS4, 1e-4);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_7, SIMPLE_IT_F_PATH_7_EPS4, JACOBI_F_PATH_7_EPS4, SEIDEL_F_PATH_7_EPS4, 1e-4);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_8, SIMPLE_IT_F_PATH_8_EPS4, JACOBI_F_PATH_8_EPS4, SEIDEL_F_PATH_8_EPS4, 1e-4);

    // Точность 1e-7
    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_1, SIMPLE_IT_F_PATH_1_EPS7, JACOBI_F_PATH_1_EPS7, SEIDEL_F_PATH_1_EPS7, 1e-7);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_2, SIMPLE_IT_F_PATH_2_EPS7, JACOBI_F_PATH_2_EPS7, SEIDEL_F_PATH_2_EPS7, 1e-7);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_3, SIMPLE_IT_F_PATH_3_EPS7, JACOBI_F_PATH_3_EPS7, SEIDEL_F_PATH_3_EPS7, 1e-7);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_4, SIMPLE_IT_F_PATH_4_EPS7, JACOBI_F_PATH_4_EPS7, SEIDEL_F_PATH_4_EPS7, 1e-7);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_5, SIMPLE_IT_F_PATH_5_EPS7, JACOBI_F_PATH_5_EPS7, SEIDEL_F_PATH_5_EPS7, 1e-7);

    generateRandomTest<Type>(4, 1.0, 20.0, IN_FILE_PATH_6);
    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_6, SIMPLE_IT_F_PATH_6_EPS7, JACOBI_F_PATH_6_EPS7, SEIDEL_F_PATH_6_EPS7, 1e-7);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_7, SIMPLE_IT_F_PATH_7_EPS7, JACOBI_F_PATH_7_EPS7, SEIDEL_F_PATH_7_EPS7, 1e-7);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_8, SIMPLE_IT_F_PATH_8_EPS7, JACOBI_F_PATH_8_EPS7, SEIDEL_F_PATH_8_EPS7, 1e-7);
}

int main(){
    temp_main<double>();
    return 0;
}