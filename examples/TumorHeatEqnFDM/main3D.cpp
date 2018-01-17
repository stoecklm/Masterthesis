/* ScaFES
 * Copyright (c) 2017, ZIH, TU Dresden, Federal Republic of Germany.
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
    std::vector<std::string> nameDatafield(3);
    nameDatafield[0] = "dummy";  /* Will not be used.                          */
    nameDatafield[1] = "T";      /* Temperature.                               */
    nameDatafield[2] = "Region"; /* Used to specify if node is brain or tumor. */
    std::vector<int> stencilWidth(3);
    stencilWidth[0] = 0;
    stencilWidth[1] = 1;
    stencilWidth[2] = 0;
    std::vector<bool> isKnownDf(3);
    isKnownDf[0] = true;
    isKnownDf[1] = false;
    isKnownDf[2] = false;
    std::vector<int> nLayers(3);
    nLayers[0] = 0;
    nLayers[1] = 0;
    nLayers[2] = 0;
    std::vector<double> defaultValue(3);
    defaultValue[0] = 0.0;
    defaultValue[1] = 0.0;
    defaultValue[2] = 0.0;
    std::vector<ScaFES::WriteHowOften> writeToFile(3);
    writeToFile[0] = ScaFES::WriteHowOften::NEVER;
    writeToFile[1] = ScaFES::WriteHowOften::LIKE_GIVEN_AT_CL;
    writeToFile[2] = ScaFES::WriteHowOften::NEVER;
    std::vector<bool> computeError(3);
    computeError[0] = false;
    computeError[1] = false;
    computeError[2] = false;
    std::vector<double> geomparamsInit;
    std::vector<bool> checkConvergence(3);
    checkConvergence[0] = false;
    checkConvergence[1] = true;
    checkConvergence[2] = false;

    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini(paramsCl.nameConfigFile(), pt);

    TumorHeatEqnFDM<double, DIM> ppp(paramsCl, gg, false, nameDatafield,
                                     stencilWidth, isKnownDf, pt, nLayers,
                                     defaultValue, writeToFile, computeError,
                                     geomparamsInit, checkConvergence);

    ppp.iterateOverTime();

    return 0;
}
