#include "gtest/gtest.h"

#include "ScaFES_GridGlobalTest.hpp"

#include "ScaFES_Ntuple.hpp"

namespace ScaFES_test
{

/*******************************************************************************
 * TESTS OF GETTER METHODS.
 ******************************************************************************/
TYPED_TEST_P(GridGlobalTest, GetNpartitionsTotalA)
{
    EXPECT_TRUE(this->mNpartitionsTotalA == this->mObjA.nPartitionsTotal());
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridGlobalTest, GetNpartitionsTotalB)
{
    EXPECT_TRUE(this->mNpartitionsTotalB == this->mObjB.nPartitionsTotal());
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridGlobalTest, GetDivideGridA)
{
    EXPECT_TRUE(this->mDivideGridA == this->mObjA.divideGrid());
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridGlobalTest, GetDivideGridB)
{
    EXPECT_TRUE(this->mDivideGridB == this->mObjB.divideGrid());
}

/*******************************************************************************
 * TESTS OF COMPARISON METHODS.
 ******************************************************************************/
TYPED_TEST_P(GridGlobalTest, EqualsOperatorComponentwise)
{
    const int DIM = 3;
    EXPECT_TRUE(this->mObjCopyA.discreteDomain().nNodes()
               == this->mObjA.discreteDomain().nNodes());
    EXPECT_TRUE(this->mObjCopyA.discreteDomain().gridsize()
               == this->mObjA.discreteDomain().gridsize());
    EXPECT_TRUE(this->mObjCopyA.discreteDomain().idxNodeFirst()
                 == this->mObjA.discreteDomain().idxNodeFirst());
    EXPECT_TRUE(this->mObjCopyA.discreteDomain().idxNodeLast()
                 == this->mObjA.discreteDomain().idxNodeLast());

    for (int ii = 0; ii < DIM; ++ii) {
        EXPECT_TRUE(this->mObjCopyA.discreteDomain().nNodes(ii)
                    == this->mObjA.discreteDomain().nNodes(ii));
    }

    for (int ii = 0; ii < DIM; ++ii) {
        EXPECT_TRUE(::fabs(this->mObjCopyA.discreteDomain().gridsize(ii)
                           - this->mObjA.discreteDomain().gridsize(ii)) < 2.2e-14);
    }
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridGlobalTest, EqualsOperator)
{
    EXPECT_TRUE(this->mObjA == this->mObjA);
    EXPECT_TRUE(this->mObjCopyA == this->mObjA);
    EXPECT_FALSE(this->mObjA == this->mObjB);
    EXPECT_FALSE(this->mObjB == this->mObjA);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridGlobalTest, UnequalsOperator)
{
    EXPECT_FALSE(this->mObjA != this->mObjA);
    EXPECT_FALSE(this->mObjCopyA != this->mObjA);
    EXPECT_TRUE(this->mObjA != this->mObjB);
    EXPECT_TRUE(this->mObjB != this->mObjA);
}

/*******************************************************************************
 * TESTS OF LIFE CYCLE METHODS.
 ******************************************************************************/
TYPED_TEST_P(GridGlobalTest, ConstructorWithParameters)
{
    const int DIM = 3;
    ScaFES::GridGlobal<DIM> objTmp(this->mNpartitionsTotalA,
                                   this->mNnodesA,
                                   this->mCoordNodeFirstA,
                                   this->mCoordNodeLastA,
                                   this->mDivideGridA);
    EXPECT_TRUE(this->mObjA == objTmp);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridGlobalTest, CopyConstructor)
{
    const int DIM = 3;
    ScaFES::GridGlobal<DIM> objTmp = this->mObjA;
    EXPECT_TRUE(this->mObjA == objTmp);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridGlobalTest, AssignmentOperator)
{
    const int DIM = 3;
    ScaFES::GridGlobal<DIM> objTmp;
    objTmp = this->mObjA;
    EXPECT_TRUE(this->mObjA == objTmp);
}

/*******************************************************************************
 * TESTS OF WORK METHODS.
 ******************************************************************************/
// TYPED_TEST_P(GridGlobalTest, ApplyRCB)
// {
//     const int DIM = 3;
//     ScaFES::GridGlobal<DIM> objTmp = this->mObjA;
//     ScaFES::Ntuple<bool, DIM> divideGrid(false, false, true);
//
//     // Check partitioning in all dimensions.
//     for (int jj = 0; jj < 1; ++jj) {
//         // TODO: Tests will fail if divideGrid will be set to "true".
//         //divideGrid[jj] = false;
//         ScaFES::Ntuple<int, DIM> nNodesPart(0);
//         std::vector< ScaFES::GridSub<DIM> > partitionsExpected;
//
//         for (int ii = 0; ii < DIM; ++ii) {
//             if (divideGrid[ii]) {
//                 nNodesPart[ii] = objTmp.nNodes(ii)
//                                  / objTmp.nPartitions();
//             }
//         }
//
//         for (int ii = 0; ii < objTmp.nPartitions(); ++ii) {
//             ScaFES::Ntuple<int, DIM> idxNodeFirstPartition
//                 = objTmp.idxNodeFirst() + ii * nNodesPart;
//             ScaFES::Ntuple<int, DIM> idxNodeLastPartition
//                 = objTmp.idxNodeLast()
//                   - (objTmp.nPartitions() - ii - 1) * nNodesPart;
//             partitionsExpected.push_back(
//                 ScaFES::GridSub<DIM>(idxNodeFirstPartition,
//                                      idxNodeLastPartition,
//                                      objTmp.idxNodeFirst(),
//                                      objTmp.idxNodeLast(),
//                                      objTmp.coordNodeFirst(),
//                                      objTmp.coordNodeLast())
//             );
//         }
//
//         objTmp.applyRCB(objTmp.nPartitions(),
//                              objTmp.idxNodeFirst(),
//                              objTmp.nNodes(),
//                              divideGrid);
//
//         for (int ii = 0; ii < objTmp.nPartitions(); ++ii) {
//             EXPECT_TRUE(partitionsExpected[ii] == objTmp.partition(ii));
//         }
//
//         // Reset cutting in jj-th dimension.
//         //divideGrid[jj] = false;
//     }
// }
// /*----------------------------------------------------------------------------*/
// TYPED_TEST_P(GridGlobalTest, SetupConnectivity)
// {
//     const int DIM = 3;
//     ScaFES::GridGlobal<DIM> objTmp = this->mObjA;
//     unsigned long int idxNodeScalar;
//     ScaFES::Ntuple<int, 2*DIM> connectExpected;
//     ScaFES::Ntuple<int, 2*DIM> connectActual;
//     objTmp.setupConnectivity();
//
//     idxNodeScalar = 0;
//     connectExpected[0] = -1;
//     connectExpected[1] =  1;
//     connectExpected[2] = -3;
//     connectExpected[3] = objTmp.nNodes(0);
//     connectExpected[4] = -5;
//     connectExpected[5] = objTmp.nNodes(0) * objTmp.nNodes(1);
//     connectActual = objTmp.connectivity(idxNodeScalar);
//     EXPECT_TRUE(connectActual == connectExpected);
//
//     idxNodeScalar = objTmp.nNodes(0) - 1;
//     connectExpected[0] =  objTmp.nNodes(0) - 2;
//     connectExpected[1] = -2;
//     connectExpected[2] = -3;
//     connectExpected[3] =  objTmp.nNodes(0) - 1 + objTmp.nNodes(0);
//     connectExpected[4] = -5;
//     connectExpected[5] = objTmp.nNodes(0) * objTmp.nNodes(1)
//                           + objTmp.nNodes(0) - 1;
//     connectActual = objTmp.connectivity(idxNodeScalar);
//     EXPECT_TRUE(connectActual == connectExpected);
// }

/*******************************************************************************
 * TESTS OF FREE METHODS.
 ******************************************************************************/
TYPED_TEST_P(GridGlobalTest, Swap)
{
    const int DIM = 3;
    ScaFES::GridGlobal<DIM> objTmpA = this->mObjA;
    ScaFES::GridGlobal<DIM> objTmpB = this->mObjB;
    ScaFES::swap(this->mObjA, this->mObjB);
    EXPECT_TRUE(objTmpA == this->mObjB);
    EXPECT_TRUE(objTmpB == this->mObjA);
}

/*******************************************************************************
 ******************************************************************************/
REGISTER_TYPED_TEST_CASE_P(GridGlobalTest,
                           GetNpartitionsTotalA,
                           GetNpartitionsTotalB,
                           GetDivideGridA,
                           GetDivideGridB,
                           EqualsOperatorComponentwise,
                           EqualsOperator,
                           UnequalsOperator,
                           //ConstructorDefault,
                           ConstructorWithParameters,
                           CopyConstructor,
                           AssignmentOperator,
//                            SetupConnectivity,
                           Swap
                          );

/*******************************************************************************
 ******************************************************************************/
typedef ::testing::Types<int> MyGridGlobalTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(My, GridGlobalTest, MyGridGlobalTestTypes);

} // End of namespace.
