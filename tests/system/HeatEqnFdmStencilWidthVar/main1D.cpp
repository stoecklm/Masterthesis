#include <vector>
#include <string>

#include "ScaFES.hpp"
#include "HeatEqnFdmStencilWidthVar.hpp"

const int DIM = 1;

int main(int argc, char *argv[]) {
    ScaFES::Parameters paramsCl(argc, argv);
    ScaFES::GridGlobal<DIM> gg(paramsCl);

    const int N_DATAFIELDS = 3;

    char currName = 'A';
    std::vector<std::string> nameDatafield(N_DATAFIELDS);
    for (int ii = 0; ii < N_DATAFIELDS; ++ii) {
        nameDatafield[ii] = currName;
        ++currName;
    }
    std::vector<int> stencilWidth(N_DATAFIELDS, 0);
    for (int ii = 0; ii < N_DATAFIELDS; ++ii) {
        stencilWidth[ii] = ii;
    }
    std::vector<bool> isKnownDf(N_DATAFIELDS, false);
    std::vector<int> nLayers(N_DATAFIELDS, 0);
    for (int ii = 0; ii < N_DATAFIELDS; ++ii) {
        nLayers[ii] = 0;
    }
    std::vector<double> defaultValue(N_DATAFIELDS);
    for (int ii = 0; ii < N_DATAFIELDS; ++ii) {
        defaultValue[ii] = 0.0;
    }
    std::vector<ScaFES::WriteHowOften> writeToFile(N_DATAFIELDS);
    for (int ii = 0; ii < N_DATAFIELDS; ++ii) {
        writeToFile[ii] = ScaFES::WriteHowOften::NEVER;
    }
    std::vector<bool> computeError(N_DATAFIELDS, true);

    std::vector<double> geomparamsInit;

    HeatEqnFdmStencilWidthVar<double, DIM> ppp(paramsCl, gg, false,
                     nameDatafield, stencilWidth,
                     isKnownDf, nLayers, defaultValue, writeToFile,
                     computeError, geomparamsInit);
    ppp.iterateOverTime();
    return 0;
}
