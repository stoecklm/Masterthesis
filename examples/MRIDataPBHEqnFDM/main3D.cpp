/* ScaFES
 * Copyright (c) 2018, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

#include "ScaFES.hpp"
#include "MRIDataPBHEqnFDM.hpp"

/** Space dimension of problem. */
const int DIM = 3;

/** Main program for MRIDataPBHEqnFDM. */
int main(int argc, char *argv[]) {
    ScaFES::Parameters paramsCl(argc, argv);
    ScaFES::GridGlobal<DIM> gg(paramsCl);
    std::vector<std::string> nameDatafield(9);
    nameDatafield[0] = "rho";        /* density brain or tumor. */
    nameDatafield[1] = "c";          /* heat capacity. */
    nameDatafield[2] = "lambda";     /* conductivity brain or tumor. */
    nameDatafield[3] = "rho_blood";  /* density blood. */
    nameDatafield[4] = "omega";      /* blood perfusion rate brain or tumor. */
    nameDatafield[5] = "c_blood";    /* heat capacity. */
    nameDatafield[6] = "T_blood";    /* blood temperature. */
    nameDatafield[7] = "q";          /* metabolic heat rate brain or tumor. */
    nameDatafield[8] = "T";          /* Temperature. */
    std::vector<int> stencilWidth(9, 0);
    stencilWidth[8] = 1;
    std::vector<bool> isKnownDf(9, true);
    isKnownDf[8] = false;
    std::vector<int> nLayers(9, 0);
    std::vector<double> defaultValue(9, 0.0);
    std::vector<ScaFES::WriteHowOften> writeToFile(9, ScaFES::WriteHowOften::NEVER);
    writeToFile[8] = ScaFES::WriteHowOften::AT_END;
    std::vector<bool> computeError(9, false);
    std::vector<bool> checkConvergence(9, false);
    checkConvergence[8] = true;
    std::vector<double> geomparamsInit;

    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini(paramsCl.nameConfigFile(), pt);

    MRIDataPBHEqnFDM<double, DIM> ppp(paramsCl, gg, false, nameDatafield,
                                      stencilWidth, isKnownDf, pt, nLayers,
                                      defaultValue, writeToFile, computeError,
                                      geomparamsInit, checkConvergence);

    ppp.iterateOverTime();

    return 0;
}
