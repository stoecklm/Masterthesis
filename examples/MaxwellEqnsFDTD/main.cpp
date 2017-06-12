/** ScaFES
 *  Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
 *  For details, see the files COPYING and LICENSE in the base directory
 *  of the package.
 */
#include "ScaFES.hpp"
#include "Problem_MaxwellEqnsFDTD.hpp"

/*----------------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
    ScaFES::Parameters paramsCl(argc, argv);
    ScaFES::GridGlobal<3> gg(paramsCl);
    Problem_MaxwellEqnsFDTD<double> ppp(paramsCl, gg);
    // Ellipsoid with center at (10,10,10) and radius of 0.5 in each direction.
    std::vector<double> geomParams(6);
    geomParams[0] = 10.0;
    geomParams[1] = 10.0;
    geomParams[2] = 10.0;
    geomParams[3] =  0.5;
    geomParams[4] =  0.5;
    geomParams[5] =  0.5;
    ppp.iterate(geomParams);

    return 0;
}
