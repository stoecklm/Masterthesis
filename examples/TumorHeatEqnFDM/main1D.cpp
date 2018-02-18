/* ScaFES
 * Copyright (c) 2017-2018, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

#include "ScaFES.hpp"
#include "TumorHeatEqnFDM.hpp"

/** Space dimension of problem. */
const int DIM = 1;

/** Main program for TumorHeatEqnFDM. */
int main(int argc, char *argv[]) {
    ScaFES::Parameters paramsCl(argc, argv);
    ScaFES::GridGlobal<DIM> gg(paramsCl);
    std::vector<std::string> nameDatafield(1);
    nameDatafield[0] = "T";      /* Temperature. */
    std::vector<int> stencilWidth(1);
    stencilWidth[0] = 1;
    std::vector<bool> isKnownDf(1);
    isKnownDf[0] = false;
    std::vector<int> nLayers(1);
    nLayers[0] = 0;
    std::vector<double> defaultValue(1);
    defaultValue[0] = 0.0;
    std::vector<ScaFES::WriteHowOften> writeToFile(1);
    writeToFile[0] = ScaFES::WriteHowOften::AT_END;
    std::vector<bool> computeError(1);
    computeError[0] = false;
    std::vector<bool> checkConvergence(1);
    checkConvergence[0] = true;
    std::vector<double> geomparamsInit;

    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini(paramsCl.nameConfigFile(), pt);

    TumorHeatEqnFDM<double, DIM> ppp(paramsCl, gg, false, nameDatafield,
                                     stencilWidth, isKnownDf, pt, nLayers,
                                     defaultValue, writeToFile, computeError,
                                     geomparamsInit, checkConvergence);

    ppp.iterateOverTime();

    return 0;
}
