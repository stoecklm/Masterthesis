/* ScaFES
 * Copyright (c) 2017, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

#include "ScaFES.hpp"
#include "PennesBioheatEqnFDM.hpp"

/* Defines types of equation which can be used for validation. */
enum typesOfEqn {constant = 0, linear = 1, quadratic = 2, cubic = 3};

/* Defines types of boundary conditions. */
enum typesOfBCs {dirichlet = 1, neumann = 2, cauchy = 3};

/** Space dimension of problem. */
const int DIM = 3;

/** Main program for PennesBioheatEqnFDM. */
int main(int argc, char *argv[]) {
    std::cout << "Testing PennesBioheatEqnFDM with linear analytical solution"
              << std::endl
              << "and Dirichlet boundary conditions."
              << std::endl << std::endl;

    ScaFES::Parameters paramsCl(argc, argv);
    ScaFES::GridGlobal<DIM> gg(paramsCl);
    std::vector<std::string> nameDatafield(3);
    nameDatafield[0] = "F";
    nameDatafield[1] = "G";
    nameDatafield[2] = "U";
    std::vector<int> stencilWidth(3);
    stencilWidth[0] = 0;
    stencilWidth[1] = 0;
    stencilWidth[2] = 1;
    std::vector<bool> isKnownDf(3);
    isKnownDf[0] = true;
    isKnownDf[1] = true;
    isKnownDf[2] = false;
    std::vector<int> nLayers(3);
    nLayers[0] = 0;
    nLayers[1] = 0;
    nLayers[2] = 0;
    std::vector<double> defaultValue(3);
    defaultValue[0] = 0.0;
    defaultValue[1] = 0.0;
    defaultValue[2] = 0.0;
    std::vector<ScaFES::WriteHowOften> writeToFile(3, ScaFES::WriteHowOften::NEVER);
    std::vector<bool> computeError(3);
    computeError[0] = false;
    computeError[1] = false;
    computeError[2] = true;
    std::vector<double> geomparamsInit;

    PennesBioheatEqnFDM<double, DIM> ppp(paramsCl, gg, false, nameDatafield, stencilWidth,
                                         isKnownDf, nLayers, defaultValue, writeToFile,
                                         computeError, geomparamsInit,
                                         linear, dirichlet);

    double sumHsquared = 0.0;
    for (std::size_t pp = 0; pp < DIM; ++pp) {
        sumHsquared += 1.0/(ppp.gridsize(pp) * ppp.gridsize(pp));
    }
    /* Stability condition: 1/(\sum_p 1/h_p^2) \ge 2 \tau */
    /* nTimesteps >= 2 * dim * (nNodes - 1)^2 */
    if ((1.0/sumHsquared) < 2.0 * ppp.tau()) {
        std::cerr << "\nERROR: Stability condition is not fulfilled."
                  << std::endl;
    }

    ppp.iterateOverTime();

    return 0;
}
