/* ScaFES
 * Copyright (c) 2017-2018, ZIH, TU Dresden, Federal Republic of Germany.
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
    std::vector<std::string> nameDatafield(2);
    nameDatafield[0] = "dummy";  /* Will not be used.                          */
    nameDatafield[1] = "T";      /* Temperature.                               */
    std::vector<int> stencilWidth(2);
    stencilWidth[0] = 0;
    stencilWidth[1] = 1;
    std::vector<bool> isKnownDf(2);
    isKnownDf[0] = true;
    isKnownDf[1] = false;
    std::vector<int> nLayers(2);
    nLayers[0] = 0;
    nLayers[1] = 0;
    std::vector<double> defaultValue(2);
    defaultValue[0] = 0.0;
    defaultValue[1] = 0.0;
    std::vector<ScaFES::WriteHowOften> writeToFile(2);
    writeToFile[0] = ScaFES::WriteHowOften::NEVER;
    writeToFile[1] = ScaFES::WriteHowOften::AT_END;
    std::vector<bool> computeError(2);
    computeError[0] = false;
    computeError[1] = false;
    std::vector<double> geomparamsInit;
    std::vector<bool> checkConvergence(2);
    checkConvergence[0] = false;
    checkConvergence[1] = true;

    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini(paramsCl.nameConfigFile(), pt);

    TumorHeatEqnFDM<double, DIM> ppp(paramsCl, gg, false, nameDatafield,
                                     stencilWidth, isKnownDf, pt, nLayers,
                                     defaultValue, writeToFile, computeError,
                                     geomparamsInit, checkConvergence);

    ppp.iterateOverTime();

    return 0;
}
