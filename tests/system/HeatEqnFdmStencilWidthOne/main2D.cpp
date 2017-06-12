#include <vector>
#include <string>

#include "ScaFES.hpp"
#include "HeatEqnFdmStencilWidthOne.hpp"

const int DIM = 2;

int main(int argc, char *argv[]) {
    ScaFES::Parameters paramsCl(argc, argv);
    ScaFES::GridGlobal<DIM> gg(paramsCl);

    const int N_DATAFIELDS = 1;

    char currName = 'A';
    std::vector<std::string> nameDatafield(N_DATAFIELDS);
    for (int ii = 0; ii < N_DATAFIELDS; ++ii) {
        nameDatafield[ii] = currName;
        ++currName;
    }
    std::vector<int> stencilWidth(N_DATAFIELDS, 0);
    for (int ii = 0; ii < N_DATAFIELDS; ++ii) {
        stencilWidth[ii] = 1; // =1 for ALL fields.
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

    HeatEqnFdmStencilWidthOne<double, DIM> ppp(paramsCl, gg, false,
                     nameDatafield, stencilWidth,
                     isKnownDf, nLayers, defaultValue, writeToFile,
                     computeError, geomparamsInit);
    ppp.iterateOverTime();
    return 0;
}

