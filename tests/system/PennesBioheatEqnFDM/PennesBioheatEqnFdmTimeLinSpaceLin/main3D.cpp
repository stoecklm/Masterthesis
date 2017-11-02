/* ScaFES
 * Copyright (c) 2017, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

#include "ScaFES.hpp"
#include "PennesBioheatEqnFDM.hpp"
#include "PennesBioheatEqnFdmTimeLinSpaceLin.hpp"

/** Space dimension of problem. */
const int DIM = 3;

/** Main program for PennesBioheatEqnFDM. */
int main(int argc, char *argv[]) {
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

    PennesBioheatEqnFdmTimeLinSpaceLin<double, DIM> ppp(paramsCl, gg, false, nameDatafield, stencilWidth,
                                isKnownDf, nLayers, defaultValue, writeToFile,
                                computeError, geomparamsInit);

    /* Stability condition: 1/((\sum_p 2/h_p^2) + 1) \ge \tau */
    double sumHsquared = 0.0;
    for (std::size_t pp = 0; pp < DIM; ++pp) {
        sumHsquared += 2.0/(ppp.gridsize(pp) * ppp.gridsize(pp));
    }
    sumHsquared += 1.0;
    if ((1.0/sumHsquared) < ppp.tau()) {
        std::cerr << "\nERROR: Stability condition is not fulfilled."
                  << std::endl;
    }

    ppp.iterateOverTime();

    return 0;
}
