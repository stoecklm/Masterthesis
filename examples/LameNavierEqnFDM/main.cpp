/* ScaFES
 * Copyright (c) 2016, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

#include "ScaFES.hpp"
#include "LameNavierEqnFDM.hpp"

/** Space dimension of problem. */
const int DIM = 3;

/** Main program for LameNavierEqnFDM. */
int main(int argc, char *argv[]) {
    ScaFES::Parameters paramsCl(argc, argv);
    ScaFES::GridGlobal<DIM> gg(paramsCl);
    std::vector<std::string> nameDatafield(9);
    nameDatafield[0] = "F1";
    nameDatafield[1] = "F2";
    nameDatafield[2] = "F3";
    nameDatafield[3] = "G1";
    nameDatafield[4] = "G2";
    nameDatafield[5] = "G3";
    nameDatafield[6] = "Y1";
    nameDatafield[7] = "Y2";
    nameDatafield[8] = "Y3";
    std::vector<int> stencilWidth(9);
    stencilWidth[0] = 0;
    stencilWidth[1] = 0;
    stencilWidth[2] = 0;
    stencilWidth[3] = 0;
    stencilWidth[4] = 0;
    stencilWidth[5] = 0;
    stencilWidth[6] = 1;
    stencilWidth[7] = 1;
    stencilWidth[8] = 1;
    std::vector<bool> isKnownDf(9);
    isKnownDf[0] = true;
    isKnownDf[1] = true;
    isKnownDf[2] = true;
    isKnownDf[3] = true;
    isKnownDf[4] = true;
    isKnownDf[5] = true;
    isKnownDf[6] = false;
    isKnownDf[7] = false;
    isKnownDf[8] = false;
    std::vector<int> nLayers(9);
    nLayers[0] = 0;
    nLayers[1] = 0;
    nLayers[2] = 0;
    nLayers[3] = 0;
    nLayers[4] = 0;
    nLayers[5] = 0;
    nLayers[6] = 0;
    nLayers[7] = 0;
    nLayers[8] = 0;
    std::vector<double> defaultValue(9);
    defaultValue[0] = 0.0;
    defaultValue[1] = 0.0;
    defaultValue[2] = 0.0;
    defaultValue[3] = 0.0;
    defaultValue[4] = 0.0;
    defaultValue[5] = 0.0;
    defaultValue[6] = 0.0;
    defaultValue[7] = 0.0;
    defaultValue[8] = 0.0;
    std::vector<ScaFES::WriteHowOften> writeToFile(9, ScaFES::WriteHowOften::LIKE_GIVEN_AT_CL);
    std::vector<bool> computeError(9);
    computeError[0] = false;
    computeError[1] = false;
    computeError[2] = false;
    computeError[3] = false;
    computeError[4] = false;
    computeError[5] = false;
    computeError[6] = false;
    computeError[7] = false;
    computeError[8] = false;
    std::vector<double> geomparamsInit;

    LameNavierEqnFDM<double, DIM> ppp(paramsCl, gg, false, nameDatafield, stencilWidth,
                                isKnownDf, nLayers, defaultValue, writeToFile,
                                computeError, geomparamsInit);
    ppp.iterateOverTime();

    return 0;
}

