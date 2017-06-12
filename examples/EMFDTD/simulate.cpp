/* ScaFES
 * Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

#include <iostream>
#include <fstream>

#include <ScaFES.hpp>

#include "EMFDTD.hpp"

int main(int argc, char *argv[])
{
    ScaFES::Parameters cl(argc, argv);

    if (cl.enabledAdolc())
        std::cout
            << "WARNING: ADOL-C enabled forward computation cannot generate reference data."
            << "\nWILL NOT WRITE REFERENCE DATA." << std::endl;

    ScaFES::GridGlobal<3> grid(cl);

    // read ptree early
    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini(cl.nameConfigFile(), pt);

    EMFDTDng<double, true> EM(cl, grid, pt);

    EM.iterateOverTime();

    // write reference data only when fwd computation without AD
    if (!cl.enabledAdolc())
        EM.writeSinkData();

    return EXIT_SUCCESS;
}
