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
    std::vector<std::string> nameDatafield(10);
    nameDatafield[0] = "rho";        /* density. */
    nameDatafield[1] = "c";          /* specific heat capacity. */
    nameDatafield[2] = "lambda";     /* conductivity. */
    nameDatafield[3] = "rho_blood";  /* density blood. */
    nameDatafield[4] = "c_blood";    /* heat capacity. */
    nameDatafield[5] = "omega";      /* blood perfusion rate. */
    nameDatafield[6] = "T_blood";    /* blood temperature. */
    nameDatafield[7] = "q";          /* metabolic heat rate. */
    nameDatafield[8] = "surface";    /* vector to specify open skull. */
    nameDatafield[9] = "T";          /* Temperature. */
    std::vector<int> stencilWidth(10, 0);
    stencilWidth[9] = 1;
    std::vector<bool> isKnownDf(10, true);
    isKnownDf[9] = false;
    std::vector<int> nLayers(10, 0);
    std::vector<double> defaultValue(10, 0.0);
    std::vector<ScaFES::WriteHowOften> writeToFile(10, ScaFES::WriteHowOften::NEVER);
    writeToFile[9] = ScaFES::WriteHowOften::AT_END;
    std::vector<bool> computeError(10, false);
    std::vector<bool> checkConvergence(10, false);
    checkConvergence[9] = true;
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
