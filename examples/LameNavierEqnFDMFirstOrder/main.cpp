/* ScaFES
 * Copyright (c) 2016, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

#include "ScaFES.hpp"
#include "LameNavierEqnFDMFirstOrder.hpp"

/** Space dimension of problem. */
const int DIM = 3;

/** Main program for LameNavierEqnFDMFirstOrder. */
int main(int argc, char *argv[]) {
    ScaFES::Parameters paramsCl(argc, argv);
    ScaFES::GridGlobal<DIM> gg(paramsCl);
    std::vector<std::string> nameDatafield(12);
    nameDatafield[0] = "F1";
    nameDatafield[1] = "F2";
    nameDatafield[2] = "F3";
    nameDatafield[3] = "G1";
    nameDatafield[4] = "G2";
    nameDatafield[5] = "G3";
    nameDatafield[6] = "U1";
    nameDatafield[7] = "U2";
    nameDatafield[8] = "U3";
    nameDatafield[9] = "V1";
    nameDatafield[10] = "V2";
    nameDatafield[11] = "V3";
    std::vector<int> stencilWidth(12);
    stencilWidth[0] = 0;
    stencilWidth[1] = 0;
    stencilWidth[2] = 0;
    stencilWidth[3] = 0;
    stencilWidth[4] = 0;
    stencilWidth[5] = 0;
    stencilWidth[6] = 1;
    stencilWidth[7] = 1;
    stencilWidth[8] = 1;
    stencilWidth[9] = 1;
    stencilWidth[10] = 1;
    stencilWidth[11] = 1;
    std::vector<bool> isKnownDf(12);
    isKnownDf[0] = true;
    isKnownDf[1] = true;
    isKnownDf[2] = true;
    isKnownDf[3] = true;
    isKnownDf[4] = true;
    isKnownDf[5] = true;
    isKnownDf[6] = false;
    isKnownDf[7] = false;
    isKnownDf[8] = false;
    isKnownDf[9] = false;
    isKnownDf[10] = false;
    isKnownDf[11] = false;
    std::vector<int> nLayers(12);
    nLayers[0] = 0;
    nLayers[1] = 0;
    nLayers[2] = 0;
    nLayers[3] = 0;
    nLayers[4] = 0;
    nLayers[5] = 0;
    nLayers[6] = 0;
    nLayers[7] = 0;
    nLayers[8] = 0;
    nLayers[9] = 0;
    nLayers[10] = 0;
    nLayers[11] = 0;
    std::vector<double> defaultValue(12);
    defaultValue[0] = 0.0;
    defaultValue[1] = 0.0;
    defaultValue[2] = 0.0;
    defaultValue[3] = 0.0;
    defaultValue[4] = 0.0;
    defaultValue[5] = 0.0;
    defaultValue[6] = 0.0;
    defaultValue[7] = 0.0;
    defaultValue[8] = 0.0;
    defaultValue[9] = 0.0;
    defaultValue[10] = 0.0;
    defaultValue[11] = 0.0;
    std::vector<ScaFES::WriteHowOften> writeToFile(12, ScaFES::WriteHowOften::LIKE_GIVEN_AT_CL);
    std::vector<bool> computeError(12);
    computeError[0] = false;
    computeError[1] = false;
    computeError[2] = false;
    computeError[3] = false;
    computeError[4] = false;
    computeError[5] = false;
    computeError[6] = false;
    computeError[7] = false;
    computeError[8] = false;
    computeError[9] = false;
    computeError[10] = false;
    computeError[11] = false;
    std::vector<double> geomparamsInit;

    LameNavierEqnFDMFirstOrder<double, DIM> ppp(paramsCl, gg, false, nameDatafield, stencilWidth,
                                isKnownDf, nLayers, defaultValue, writeToFile,
                                computeError, geomparamsInit);
    ppp.iterateOverTime();

    return 0;
}

