/* ScaFES
 * Copyright (c) 2017, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

#include "ScaFES.hpp"
#include "InitFile.hpp"

/** Space dimension of problem. */
const int DIM = 3;

/** Main program for InitFile. */
int main(int argc, char *argv[]) {
    ScaFES::Parameters paramsCl(argc, argv);
    ScaFES::GridGlobal<DIM> gg(paramsCl);
    std::vector<std::string> nameDatafield(2);
    nameDatafield[0] = "G";
    nameDatafield[1] = "U";
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
    writeToFile[1] = ScaFES::WriteHowOften::LIKE_GIVEN_AT_CL;
    std::vector<bool> computeError(2);
    computeError[0] = false;
    computeError[1] = false;
    std::vector<double> geomparamsInit;

    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini(paramsCl.nameConfigFile(), pt);

    InitFile<double, DIM> ppp(paramsCl, gg, false, nameDatafield, stencilWidth,
                              isKnownDf, pt, nLayers, defaultValue, writeToFile,
                              computeError, geomparamsInit);

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
