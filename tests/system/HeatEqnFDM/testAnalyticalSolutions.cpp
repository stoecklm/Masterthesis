/* ScaFES
 * Copyright (c) 2017, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file testAnalyticalSolutions.cpp
 *  @brief Contains tests for analyticalSolutions.hpp.
 */

#include "ScaFES.hpp"
#include "analyticalSolutions.hpp"

template<typename TT, std::size_t DIM>
TT testLinFunc(ScaFES::Ntuple<TT,DIM> const& x, TT const& t) {
    if (DIM == 1) {
        double x1 = x[0];
        return (t * (1.0 + x1));
    } else if (DIM == 2) {
        double x1 = x[0];
        double x2 = x[1];
        return (t * (1.0 + x1 + x2));
    } else if (DIM == 3) {
        double x1 = x[0];
        double x2 = x[1];
        double x3 = x[2];
        return (t * (1.0 + x1 + x2 + x3));
    } else {
        std::cout << "Invalid dimension in testLinFunc()" << std::endl;
        return -1.0;
    }
}

template<typename TT, std::size_t DIM>
TT testLinFuncTimeDerivative(ScaFES::Ntuple<TT,DIM> const& x) {
    if (DIM == 1) {
        double x1 = x[0];
        return (1.0 + x1);
    } else if (DIM == 2) {
        double x1 = x[0];
        double x2 = x[1];
        return (1.0 + x1 + x2);
    } else if (DIM == 3) {
        double x1 = x[0];
        double x2 = x[1];
        double x3 = x[2];
        return (1.0 + x1 + x2 + x3);
    } else {
        std::cout << "Invalid dimension in testLinFuncTimeDerivative()" << std::endl;
        return -1.0;
    }
}

int main(void) {

    std::size_t const DIM1 = 1;
    std::size_t const DIM2 = 2;
    std::size_t const DIM3 = 3;

    double time1 = 1.25;
    double time2 = 2.25;

    ScaFES::Ntuple<double,DIM1> X_DIM1_1 = (0.15);
    ScaFES::Ntuple<double,DIM2> X_DIM2_1 = (0.25, 0.55);
    ScaFES::Ntuple<double,DIM3> X_DIM3_1 = (0.55, 1.55, 2.55);

    /* Constant analytical solutions */
    assert(consFunc(X_DIM1_1, time1) == time1);
    assert(consFunc(X_DIM2_1, time1) == time1);
    assert(consFunc(X_DIM3_1, time1) == time1);
    assert(consFunc(X_DIM1_1, time2) == time2);
    assert(consFunc(X_DIM2_1, time2) == time2);
    assert(consFunc(X_DIM3_1, time2) == time2);

    /* Time derivative for constant analytical solution */
    assert(consFuncTimeDerivative(X_DIM1_1) == 1.0);
    assert(consFuncTimeDerivative(X_DIM2_1) == 1.0);
    assert(consFuncTimeDerivative(X_DIM3_1) == 1.0);

    /* First order space derivative for constant analytical solution */
    assert(consFuncSpaceDerivative1stOrder(X_DIM1_1, time1, DIM1) == 0.0);
    assert(consFuncSpaceDerivative1stOrder(X_DIM2_1, time1, DIM1) == 0.0);
    assert(consFuncSpaceDerivative1stOrder(X_DIM3_1, time1, DIM1) == 0.0);
    assert(consFuncSpaceDerivative1stOrder(X_DIM1_1, time2, DIM1) == 0.0);
    assert(consFuncSpaceDerivative1stOrder(X_DIM2_1, time2, DIM1) == 0.0);
    assert(consFuncSpaceDerivative1stOrder(X_DIM3_1, time2, DIM1) == 0.0);

    /* Sum of second order space derivatives for constant analytical solution */
    assert(consFuncSumOfSpaceDerivatives2ndOrder(X_DIM1_1, time1) == 0.0);
    assert(consFuncSumOfSpaceDerivatives2ndOrder(X_DIM2_1, time1) == 0.0);
    assert(consFuncSumOfSpaceDerivatives2ndOrder(X_DIM3_1, time1) == 0.0);
    assert(consFuncSumOfSpaceDerivatives2ndOrder(X_DIM1_1, time2) == 0.0);
    assert(consFuncSumOfSpaceDerivatives2ndOrder(X_DIM2_1, time2) == 0.0);
    assert(consFuncSumOfSpaceDerivatives2ndOrder(X_DIM3_1, time2) == 0.0);

    /* Linear analytical solutions */
    assert(linFunc(X_DIM1_1, time1) == testLinFunc(X_DIM1_1, time1));
    assert(linFunc(X_DIM2_1, time1) == testLinFunc(X_DIM2_1, time1));
    assert(linFunc(X_DIM3_1, time1) == testLinFunc(X_DIM3_1, time1));
    assert(linFunc(X_DIM1_1, time2) == testLinFunc(X_DIM1_1, time2));
    assert(linFunc(X_DIM2_1, time2) == testLinFunc(X_DIM2_1, time2));
    assert(linFunc(X_DIM3_1, time2) == testLinFunc(X_DIM3_1, time2));

    /* Time derivative for linear analytical solution */
    assert(linFuncTimeDerivative(X_DIM1_1) == testLinFuncTimeDerivative(X_DIM1_1));
    assert(linFuncTimeDerivative(X_DIM2_1) == testLinFuncTimeDerivative(X_DIM2_1));
    assert(linFuncTimeDerivative(X_DIM3_1) == testLinFuncTimeDerivative(X_DIM3_1));

    /* First order space derivative for linear analytical solution */
    assert(linFuncSpaceDerivative1stOrder(X_DIM1_1, time1, DIM1) == time1);
    assert(linFuncSpaceDerivative1stOrder(X_DIM2_1, time1, DIM1) == time1);
    assert(linFuncSpaceDerivative1stOrder(X_DIM3_1, time1, DIM1) == time1);
    assert(linFuncSpaceDerivative1stOrder(X_DIM1_1, time2, DIM1) == time2);
    assert(linFuncSpaceDerivative1stOrder(X_DIM2_1, time2, DIM1) == time2);
    assert(linFuncSpaceDerivative1stOrder(X_DIM3_1, time2, DIM1) == time2);

    /* Sum of second order space derivatives for linear analytical solution */
    assert(linFuncSumOfSpaceDerivatives2ndOrder(X_DIM1_1, time1) == 0.0);
    assert(linFuncSumOfSpaceDerivatives2ndOrder(X_DIM2_1, time1) == 0.0);
    assert(linFuncSumOfSpaceDerivatives2ndOrder(X_DIM3_1, time1) == 0.0);
    assert(linFuncSumOfSpaceDerivatives2ndOrder(X_DIM1_1, time2) == 0.0);
    assert(linFuncSumOfSpaceDerivatives2ndOrder(X_DIM2_1, time2) == 0.0);
    assert(linFuncSumOfSpaceDerivatives2ndOrder(X_DIM3_1, time2) == 0.0);

    std::cout << "All tests passed" << std::endl;

    return 0;
}
