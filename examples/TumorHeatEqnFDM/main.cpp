/* ScaFES
 * Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

#include "ScaFES.hpp"
#include "TumorHeatEqnFDM.hpp"

/** Space dimension of problem. */
const int DIM = 3;

/** Main program for TumorHeatEqnFDM. */
int main(int argc, char *argv[]) {
    ScaFES::Parameters paramsCl(argc, argv);
    ScaFES::GridGlobal<DIM> gg(paramsCl);
    std::vector<std::string> nameDatafield(6);
    nameDatafield[0] = "rho";     /* rho (density) vector.                 */
    nameDatafield[1] = "c";       /* c (specific heat capacity) vector.    */
    nameDatafield[2] = "lambda";  /* lambda (thermal conductivity) vector. */
    nameDatafield[3] = "w";       /* w (perfusion) vector.                 */
    nameDatafield[4] = "q_dot_m"; /* q_dot_m (metabolism) vector.          */
    nameDatafield[5] = "T";       /* Temperature vector.                   */
    std::vector<int> stencilWidth(6);
    stencilWidth[0] = 0;
    stencilWidth[1] = 0;
    stencilWidth[2] = 0;
    stencilWidth[3] = 0;
    stencilWidth[4] = 0;
    stencilWidth[5] = 1;
    std::vector<bool> isKnownDf(6);
    isKnownDf[0] = true;
    isKnownDf[1] = true;
    isKnownDf[2] = true;
    isKnownDf[3] = true;
    isKnownDf[4] = true;
    isKnownDf[5] = false;
    std::vector<int> nLayers(6);
    nLayers[0] = 0;
    nLayers[1] = 0;
    nLayers[2] = 0;
    nLayers[3] = 0;
    nLayers[4] = 0;
    nLayers[5] = 0;
    std::vector<double> defaultValue(6);
    defaultValue[0] = 0.0;
    defaultValue[1] = 0.0;
    defaultValue[2] = 0.0;
    defaultValue[3] = 0.0;
    defaultValue[4] = 0.0;
    defaultValue[5] = 0.0;
    std::vector<ScaFES::WriteHowOften> writeToFile(6, ScaFES::WriteHowOften::LIKE_GIVEN_AT_CL);
    std::vector<bool> computeError(6);
    computeError[0] = false;
    computeError[1] = false;
    computeError[2] = false;
    computeError[3] = false;
    computeError[4] = false;
    computeError[5] = false;
    std::vector<double> geomparamsInit;

    TumorHeatEqnFDM<double, DIM> ppp(paramsCl, gg, false, nameDatafield, stencilWidth,
                                isKnownDf, nLayers, defaultValue, writeToFile,
                                computeError, geomparamsInit);
    ppp.iterateOverTime();

    return 0;
}
