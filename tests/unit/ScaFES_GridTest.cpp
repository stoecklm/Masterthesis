#include <algorithm>

#include "gtest/gtest.h"

#include "ScaFES_GridTest.hpp"

#include "ScaFES_Ntuple.hpp"
#include "ScaFES_Grid.hpp"

namespace ScaFES_test
{

const double EPS = 2.2e-14;

/*******************************************************************************
 * TESTS OF GETTER METHODS.
 ******************************************************************************/
TYPED_TEST_P(GridTest, GlobalNodeFirst)
{
    EXPECT_TRUE(this->mIdxNodeFirstA == this->mObjA.idxNodeFirst());
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, GlobalNodeLast)
{
    EXPECT_TRUE(this->mIdxNodeLastA == this->mObjA.idxNodeLast());
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, GlobalNodeFirstInOneDirection)
{
    const int DIM = 3;

    for (int ii = 0; ii < DIM; ++ii) {
        EXPECT_TRUE(this->mIdxNodeFirstA.elem(ii)
                    == this->mObjA.idxNodeFirst(ii));
    }
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, GlobalNodeLastInOneDirection)
{
    const int DIM = 3;

    for (int ii = 0; ii < DIM; ++ii) {
        EXPECT_TRUE(this->mIdxNodeLastA.elem(ii)
                    == this->mObjA.idxNodeLast(ii));
    }
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, Nnodes)
{
    EXPECT_TRUE(this->mNnodesA == this->mObjA.nNodes());
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, NnodesInOneDirection)
{
    const int DIM = 3;

    for (int ii = 0; ii < DIM; ++ii) {
        EXPECT_TRUE(this->mNnodesA.elem(ii) == this->mObjA.nNodes(ii));
    }
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, GridSize)
{
    const int DIM = 3;
    ScaFES::Ntuple<double, DIM> tmpGridSize = this->mObjA.gridsize();

    for (int ii = 0; ii < DIM; ++ii) {
        EXPECT_TRUE(::fabs(this->mGridSizeA.elem(ii)
                           - tmpGridSize.elem(ii)) < EPS);
    }
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, GridSizeInOneDirection)
{
    const int DIM = 3;

    for (int ii = 0; ii < DIM; ++ii) {
        EXPECT_TRUE(::fabs(this->mGridSizeA.elem(ii)
                           - this->mObjA.gridsize(ii)) < EPS);
    }
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, CoordNodeFirst)
{
    const int DIM = 3;
    ScaFES::Ntuple<double, DIM> tmpCoordNodeFirst
        = this->mObjA.coordNodeFirst();

    for (int ii = 0; ii < DIM; ++ii) {
        EXPECT_TRUE(this->mCoordNodeFirstA.elem(ii)
                    == tmpCoordNodeFirst.elem(ii));
    }
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, CoordNodeLast)
{
    const int DIM = 3;
    ScaFES::Ntuple<double, DIM> tmpCoordNodeLast = this->mObjA.coordNodeLast();

    for (int ii = 0; ii < DIM; ++ii) {
        EXPECT_TRUE(::fabs(this->mCoordNodeLastA.elem(ii)
                           - tmpCoordNodeLast.elem(ii)) < EPS);
    }
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, IsNumberedStyleC)
{
    EXPECT_TRUE(this->mIsNumberedStyleCA == this->mObjA.isNumberedStyleC());
}

/*******************************************************************************
 * TESTS OF operator==().
 ******************************************************************************/
TYPED_TEST_P(GridTest, OperatorEqualsComponentwise)
{
    const int DIM = 3;
    EXPECT_TRUE(this->mObjCopyA.nNodes() == this->mObjA.nNodes());
    EXPECT_TRUE(this->mObjCopyA.gridsize() == this->mObjA.gridsize());
    EXPECT_TRUE(this->mObjCopyA.idxNodeFirst() == this->mObjA.idxNodeFirst());
    EXPECT_TRUE(this->mObjCopyA.idxNodeLast() == this->mObjA.idxNodeLast());

    for (int ii = 0; ii < DIM; ++ii) {
        EXPECT_TRUE(this->mObjCopyA.nNodes(ii) == this->mObjA.nNodes(ii));
    }

    for (int ii = 0; ii < DIM; ++ii) {
        EXPECT_TRUE(::fabs(this->mObjCopyA.gridsize(ii)
                           - this->mObjA.gridsize(ii)) < EPS);
    }
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, OperatorEquals)
{
    EXPECT_TRUE(this->mObjA == this->mObjA);
    EXPECT_TRUE(this->mObjCopyA == this->mObjA);
    EXPECT_FALSE(this->mObjLess == this->mObjGreater);
    EXPECT_FALSE(this->mObjGreater == this->mObjLess);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, OperatorUnequals)
{
    EXPECT_FALSE(this->mObjA != this->mObjA);
    EXPECT_FALSE(this->mObjCopyA != this->mObjA);
    EXPECT_TRUE(this->mObjLess != this->mObjGreater);
    EXPECT_TRUE(this->mObjGreater != this->mObjLess);
}

/*******************************************************************************
 * TESTS OF LIFE CYCLE METHODS.
 ******************************************************************************/
TYPED_TEST_P(GridTest, ConstructorWithParameters)
{
    const int DIM = 3;
    ScaFES::Grid<DIM> objTmp(this->mIdxNodeFirstA,
                             this->mIdxNodeLastA,
                             this->mCoordNodeFirstA,
                             this->mCoordNodeLastA);
    EXPECT_TRUE(this->mObjA == objTmp);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, ConstructorDefault)
{
    const int DIM = 3;
    ScaFES::Grid<DIM> objTmp;
    EXPECT_TRUE((ScaFES::Ntuple<double, DIM>(1.0) == objTmp.gridsize()));
    EXPECT_TRUE((ScaFES::Ntuple<int, DIM>(0) == objTmp.idxNodeFirst()));
    EXPECT_TRUE((ScaFES::Ntuple<int, DIM>(1) == objTmp.idxNodeLast()));
    EXPECT_TRUE((ScaFES::Ntuple<int, DIM>(2) == objTmp.nNodes()));
    EXPECT_TRUE((ScaFES::Ntuple<double, DIM>(0.0) == objTmp.coordNodeFirst()));
    EXPECT_TRUE((ScaFES::Ntuple<double, DIM>(1.0) == objTmp.coordNodeLast()));
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, CopyConstructor)
{
    const int DIM = 3;
    ScaFES::Grid<DIM> objTmp = this->mObjA;
    EXPECT_TRUE(this->mObjA == objTmp);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, AssignmentOperator)
{
    const int DIM = 3;
    ScaFES::Grid<DIM> objTmp;
    objTmp = this->mObjA;
    EXPECT_TRUE(this->mObjA == objTmp);
}

/*******************************************************************************
 * TESTS OF WORK METHODS.
 ******************************************************************************/
TYPED_TEST_P(GridTest, IdxNodeTuple2Scalar)
{
    const int DIM = 3;
    ScaFES::Grid<DIM> objTmp = this->mObjA;
    unsigned long int idxNodeExpected;
    unsigned long int idxNodeActual;
    ScaFES::Ntuple<int, DIM> idxNode;

    idxNode[0] = objTmp.idxNodeFirst().elem(0);
    idxNode[1] = objTmp.idxNodeFirst().elem(1);
    idxNode[2] = objTmp.idxNodeFirst().elem(2);
    idxNodeActual = objTmp.idxNodeTuple2Scalar(idxNode);
    idxNodeExpected = 0;
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeLast().elem(0);
    idxNode[1] = objTmp.idxNodeFirst().elem(1);
    idxNode[2] = objTmp.idxNodeFirst().elem(2);
    idxNodeActual = objTmp.idxNodeTuple2Scalar(idxNode);
    idxNodeExpected = objTmp.nNodes(0) - 1;
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeLast().elem(0);
    idxNode[1] = objTmp.idxNodeLast().elem(1);
    idxNode[2] = objTmp.idxNodeFirst().elem(2);
    idxNodeActual = objTmp.idxNodeTuple2Scalar(idxNode);
    idxNodeExpected = objTmp.nNodes(0) * objTmp.nNodes(1) - 1;
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeFirst().elem(0);
    idxNode[1] = objTmp.idxNodeLast().elem(1);
    idxNode[2] = objTmp.idxNodeFirst().elem(2);
    idxNodeActual = objTmp.idxNodeTuple2Scalar(idxNode);
    idxNodeExpected = objTmp.nNodes(0) * objTmp.nNodes(1) - objTmp.nNodes(0);
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);


    idxNode[0] = objTmp.idxNodeFirst().elem(0);
    idxNode[1] = objTmp.idxNodeFirst().elem(1);
    idxNode[2] = objTmp.idxNodeLast().elem(2);
    idxNodeActual = objTmp.idxNodeTuple2Scalar(idxNode);
    idxNodeExpected = objTmp.nNodes(0) * objTmp.nNodes(1)
                      * (objTmp.nNodes(2) - 1);
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeLast().elem(0);
    idxNode[1] = objTmp.idxNodeFirst().elem(1);
    idxNode[2] = objTmp.idxNodeLast().elem(2);
    idxNodeActual = objTmp.idxNodeTuple2Scalar(idxNode);
    idxNodeExpected = objTmp.nNodes(0) * objTmp.nNodes(1)
                      * (objTmp.nNodes(2) - 1)
                      + objTmp.nNodes(0) - 1;
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeLast().elem(0);
    idxNode[1] = objTmp.idxNodeLast().elem(1);
    idxNode[2] = objTmp.idxNodeLast().elem(2);
    idxNodeActual = objTmp.idxNodeTuple2Scalar(idxNode);
    idxNodeExpected = objTmp.nNodes(0) * objTmp.nNodes(1)
                      * objTmp.nNodes(2) - 1;
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeFirst().elem(0);
    idxNode[1] = objTmp.idxNodeLast().elem(1);
    idxNode[2] = objTmp.idxNodeLast().elem(2);
    idxNodeActual = objTmp.idxNodeTuple2Scalar(idxNode);
    idxNodeExpected = objTmp.nNodes(0) * objTmp.nNodes(1) * objTmp.nNodes(2)
                      - objTmp.nNodes(0);
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    // Arbitrary node.
    idxNode[0] = objTmp.idxNodeLast().elem(0) - 1;
    idxNode[1] = objTmp.idxNodeLast().elem(1) - 1;
    idxNode[2] = objTmp.idxNodeLast().elem(2) - 1;
    idxNodeActual = objTmp.idxNodeTuple2Scalar(idxNode);
    idxNodeExpected = objTmp.nNodes(0) * objTmp.nNodes(1)
                      * (objTmp.nNodes(2) - 2)
                      +  objTmp.nNodes(0) * (objTmp.nNodes(1) - 2)
                      + objTmp.nNodes(0) - 2;
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);
}
// /*----------------------------------------------------------------------------*/
// TYPED_TEST_P(GridTest, LocNode2IdxNode)
// {
//     const int DIM = 3;
//     ScaFES::Grid<DIM> objTmp = this->mObjA;
//     unsigned long int idxNodeScalar;
//     ScaFES::Ntuple<int, DIM> idxNodeTupleActual;
//     ScaFES::Ntuple<int, DIM> idxNodeTupleExpected;
// 
//     idxNodeScalar = 0;
//     idxNodeTupleExpected[0] = objTmp.idxNodeFirst().elem(0);
//     idxNodeTupleExpected[1] = objTmp.idxNodeFirst().elem(1);
//     idxNodeTupleExpected[2] = objTmp.idxNodeFirst().elem(2);
//     idxNodeTupleActual = objTmp.idxNode2IdxNode(idxNodeScalar);
//     EXPECT_TRUE(idxNodeTupleExpected == idxNodeTupleActual);
// 
//     idxNodeScalar = objTmp.nNodes(0) - 1;
//     idxNodeTupleExpected[0] = objTmp.idxNodeLast().elem(0);
//     idxNodeTupleExpected[1] = objTmp.idxNodeFirst().elem(1);
//     idxNodeTupleExpected[2] = objTmp.idxNodeFirst().elem(2);
//     idxNodeTupleActual = objTmp.idxNode2IdxNode(idxNodeScalar);
//     EXPECT_TRUE(idxNodeTupleExpected == idxNodeTupleActual);
// 
//     idxNodeScalar = objTmp.nNodes(0) * objTmp.nNodes(1) - 1;
//     idxNodeTupleExpected[0] = objTmp.idxNodeLast().elem(0);
//     idxNodeTupleExpected[1] = objTmp.idxNodeLast().elem(1);
//     idxNodeTupleExpected[2] = objTmp.idxNodeFirst().elem(2);
//     idxNodeTupleActual = objTmp.idxNode2IdxNode(idxNodeScalar);
//     EXPECT_TRUE(idxNodeTupleExpected == idxNodeTupleActual);
// 
//     idxNodeScalar = objTmp.nNodes(0) * objTmp.nNodes(1)
//                      - objTmp.nNodes(0);
//     idxNodeTupleExpected[0] = objTmp.idxNodeFirst().elem(0);
//     idxNodeTupleExpected[1] = objTmp.idxNodeLast().elem(1);
//     idxNodeTupleExpected[2] = objTmp.idxNodeFirst().elem(2);
//     idxNodeTupleActual = objTmp.idxNode2IdxNode(idxNodeScalar);
//     EXPECT_TRUE(idxNodeTupleExpected == idxNodeTupleActual);
// 
// 
//     idxNodeScalar = objTmp.nNodes(0)
//                      * objTmp.nNodes(1) * (objTmp.nNodes(2) - 1);
//     idxNodeTupleExpected[0] = objTmp.idxNodeFirst().elem(0);
//     idxNodeTupleExpected[1] = objTmp.idxNodeFirst().elem(1);
//     idxNodeTupleExpected[2] = objTmp.idxNodeLast().elem(2);
//     idxNodeTupleActual = objTmp.idxNode2IdxNode(idxNodeScalar);
//     EXPECT_TRUE(idxNodeTupleExpected == idxNodeTupleActual);
// 
//     idxNodeScalar = objTmp.nNodes(0)
//                      * objTmp.nNodes(1) * (objTmp.nNodes(2) - 1)
//                      + objTmp.nNodes(0) - 1;
//     idxNodeTupleExpected[0] = objTmp.idxNodeLast().elem(0);
//     idxNodeTupleExpected[1] = objTmp.idxNodeFirst().elem(1);
//     idxNodeTupleExpected[2] = objTmp.idxNodeLast().elem(2);
//     idxNodeTupleActual = objTmp.idxNode2IdxNode(idxNodeScalar);
//     EXPECT_TRUE(idxNodeTupleExpected == idxNodeTupleActual);
// 
//     idxNodeScalar = objTmp.nNodes(0) * objTmp.nNodes(1)
//                      * objTmp.nNodes(2) - 1;
//     idxNodeTupleExpected[0] = objTmp.idxNodeLast().elem(0);
//     idxNodeTupleExpected[1] = objTmp.idxNodeLast().elem(1);
//     idxNodeTupleExpected[2] = objTmp.idxNodeLast().elem(2);
//     idxNodeTupleActual = objTmp.idxNode2IdxNode(idxNodeScalar);
//     EXPECT_TRUE(idxNodeTupleExpected == idxNodeTupleActual);
// 
//     idxNodeScalar = objTmp.nNodes(0) * objTmp.nNodes(1)
//                      * objTmp.nNodes(2) - objTmp.nNodes(0);
//     idxNodeTupleExpected[0] = objTmp.idxNodeFirst().elem(0);
//     idxNodeTupleExpected[1] = objTmp.idxNodeLast().elem(1);
//     idxNodeTupleExpected[2] = objTmp.idxNodeLast().elem(2);
//     idxNodeTupleActual = objTmp.idxNode2IdxNode(idxNodeScalar);
//     EXPECT_TRUE(idxNodeTupleExpected == idxNodeTupleActual);
// }
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, Coordinates)
{
    const int DIM = 3;
    ScaFES::Grid<DIM> objTmp = this->mObjA;
    ScaFES::Ntuple<double, DIM> coordNodeExpected;
    ScaFES::Ntuple<double, DIM> coordNodeActual;
    ScaFES::Ntuple<int, DIM> idxNode;

    idxNode[0] = objTmp.idxNodeFirst().elem(0);
    idxNode[1] = objTmp.idxNodeFirst().elem(1);
    idxNode[2] = objTmp.idxNodeFirst().elem(2);
    coordNodeActual = objTmp.coordinates(idxNode);
    coordNodeExpected = ScaFES::Ntuple<double, DIM>(
                            objTmp.coordNodeFirst().elem(0),
                            objTmp.coordNodeFirst().elem(1),
                            objTmp.coordNodeFirst().elem(2) );
    EXPECT_TRUE(coordNodeExpected == coordNodeActual);

    idxNode[0] = objTmp.idxNodeLast().elem(0);
    idxNode[1] = objTmp.idxNodeFirst().elem(1);
    idxNode[2] = objTmp.idxNodeFirst().elem(2);
    coordNodeActual = objTmp.coordinates(idxNode);
    coordNodeExpected = ScaFES::Ntuple<double, DIM>(
                            objTmp.coordNodeLast().elem(0),
                            objTmp.coordNodeFirst().elem(1),
                            objTmp.coordNodeFirst().elem(2) );
    EXPECT_TRUE(coordNodeExpected == coordNodeActual);

    idxNode[0] = objTmp.idxNodeLast().elem(0);
    idxNode[1] = objTmp.idxNodeLast().elem(1);
    idxNode[2] = objTmp.idxNodeFirst().elem(2);
    coordNodeActual = objTmp.coordinates(idxNode);
    coordNodeExpected = ScaFES::Ntuple<double, DIM>(
                            objTmp.coordNodeLast().elem(0),
                            objTmp.coordNodeLast().elem(1),
                            objTmp.coordNodeFirst().elem(2) );
    EXPECT_TRUE(coordNodeExpected == coordNodeActual);

    idxNode[0] = objTmp.idxNodeFirst().elem(0);
    idxNode[1] = objTmp.idxNodeLast().elem(1);
    idxNode[2] = objTmp.idxNodeFirst().elem(2);
    coordNodeActual = objTmp.coordinates(idxNode);
    coordNodeExpected = ScaFES::Ntuple<double, DIM>(
                            objTmp.coordNodeFirst().elem(0),
                            objTmp.coordNodeLast().elem(1),
                            objTmp.coordNodeFirst().elem(2) );
    EXPECT_TRUE(coordNodeExpected == coordNodeActual);


    idxNode[0] = objTmp.idxNodeFirst().elem(0);
    idxNode[1] = objTmp.idxNodeFirst().elem(1);
    idxNode[2] = objTmp.idxNodeLast().elem(2);
    coordNodeActual = objTmp.coordinates(idxNode);
    coordNodeExpected = ScaFES::Ntuple<double, DIM>(
                            objTmp.coordNodeFirst().elem(0),
                            objTmp.coordNodeFirst().elem(1),
                            objTmp.coordNodeLast().elem(2) );
    EXPECT_TRUE(coordNodeExpected == coordNodeActual);

    idxNode[0] = objTmp.idxNodeLast().elem(0);
    idxNode[1] = objTmp.idxNodeFirst().elem(1);
    idxNode[2] = objTmp.idxNodeLast().elem(2);
    coordNodeActual = objTmp.coordinates(idxNode);
    coordNodeExpected = ScaFES::Ntuple<double, DIM>(
                            objTmp.coordNodeLast().elem(0),
                            objTmp.coordNodeFirst().elem(1),
                            objTmp.coordNodeLast().elem(2) );
    EXPECT_TRUE(coordNodeExpected == coordNodeActual);

    idxNode[0] = objTmp.idxNodeLast().elem(0);
    idxNode[1] = objTmp.idxNodeLast().elem(1);
    idxNode[2] = objTmp.idxNodeLast().elem(2);
    coordNodeActual = objTmp.coordinates(idxNode);
    coordNodeExpected = ScaFES::Ntuple<double, DIM>(
                            objTmp.coordNodeLast().elem(0),
                            objTmp.coordNodeLast().elem(1),
                            objTmp.coordNodeLast().elem(2) );
    EXPECT_TRUE(coordNodeExpected == coordNodeActual);

    idxNode[0] = objTmp.idxNodeFirst().elem(0);
    idxNode[1] = objTmp.idxNodeLast().elem(1);
    idxNode[2] = objTmp.idxNodeLast().elem(2);
    coordNodeActual = objTmp.coordinates(idxNode);
    coordNodeExpected = ScaFES::Ntuple<double, DIM>(
                            objTmp.coordNodeFirst().elem(0),
                            objTmp.coordNodeLast().elem(1),
                            objTmp.coordNodeLast().elem(2) );
    EXPECT_TRUE(coordNodeExpected == coordNodeActual);

    // Arbitrary node.
    idxNode[0] = objTmp.idxNodeLast().elem(0) - 1;
    idxNode[1] = objTmp.idxNodeLast().elem(1) - 1;
    idxNode[2] = objTmp.idxNodeLast().elem(2) - 1;
    coordNodeActual = objTmp.coordinates(idxNode);
    coordNodeExpected = ScaFES::Ntuple<double, DIM>(
                          objTmp.coordNodeLast().elem(0) - objTmp.gridsize(0),
                          objTmp.coordNodeLast().elem(1) - objTmp.gridsize(1),
                          objTmp.coordNodeLast().elem(2) - objTmp.gridsize(2));
    EXPECT_TRUE(coordNodeExpected == coordNodeActual);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, MethodInside)
{
    const int DIM = 3;
    ScaFES::Grid<DIM> objTmp = this->mObjA;
    EXPECT_TRUE(objTmp.inside(objTmp.idxNodeFirst()));
    EXPECT_TRUE(objTmp.inside(objTmp.idxNodeLast()));

    ScaFES::Ntuple<int, DIM> shift(1);
    EXPECT_FALSE(objTmp.inside(objTmp.idxNodeFirst() - shift));
    EXPECT_FALSE(objTmp.inside(objTmp.idxNodeLast() + shift));
    EXPECT_TRUE(objTmp.inside(objTmp.idxNodeFirst() + shift));
    EXPECT_TRUE(objTmp.inside(objTmp.idxNodeLast() - shift));
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, UnionWith)
{
    const int DIM = 3;
    ScaFES::Ntuple<int, DIM> idxNode1(1, 2, 3);
    ScaFES::Ntuple<int, DIM> idxNode2(2, 4, 6);
    ScaFES::Ntuple<int, DIM> idxNode3(3, 6, 9);
    ScaFES::Ntuple<double, DIM> coordNode1(1.0, 2.0, 3.0);
    ScaFES::Ntuple<double, DIM> coordNode2(2.0, 4.0, 6.0);
    ScaFES::Ntuple<double, DIM> coordNode3(3.0, 6.0, 9.0);
    ScaFES::Grid<DIM> objTmpA(idxNode1,
                              idxNode3,
                              coordNode1,
                              coordNode3);
    ScaFES::Grid<DIM> objTmpB(idxNode2,
                              idxNode3,
                              coordNode2,
                              coordNode3);
    ScaFES::Grid<DIM> objExpected(idxNode1,
                                  idxNode3,
                                  coordNode1,
                                  coordNode3);
    EXPECT_TRUE(objExpected == objTmpA.unionWith(objTmpB));
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, IntersectWith)
{
    const int DIM = 3;
    ScaFES::Ntuple<int, DIM> idxNode1(1, 2, 3);
    ScaFES::Ntuple<int, DIM> idxNode2(2, 4, 6);
    ScaFES::Ntuple<int, DIM> idxNode3(3, 6, 9);
    ScaFES::Ntuple<double, DIM> coordNode1(1.0, 2.0, 3.0);
    ScaFES::Ntuple<double, DIM> coordNode2(2.0, 4.0, 6.0);
    ScaFES::Ntuple<double, DIM> coordNode3(3.0, 6.0, 9.0);
    ScaFES::Grid<DIM> objTmpA(idxNode1,
                              idxNode3,
                              coordNode1,
                              coordNode3);
    ScaFES::Grid<DIM> objTmpB(idxNode2,
                              idxNode3,
                              coordNode2,
                              coordNode3);
    ScaFES::Grid<DIM> objExpected(idxNode2,
                                  idxNode3,
                                  coordNode2,
                                  coordNode3);
    EXPECT_TRUE(objExpected == objTmpA.intersectWith(objTmpB));
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, MethodValid)
{
    const int DIM = 3;
    ScaFES::Grid<DIM> objTmp = this->mObjA;
    EXPECT_TRUE(objTmp.valid());
}
/*---------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, ExtendRight)
{
    const int DIM = 3;
    int borderwidth;
    int direction = 0;
    ScaFES::Grid<DIM> objTmp = this->mObjA;
    ScaFES::Grid<DIM> objExpected;
    ScaFES::Ntuple<int, DIM> idxNodeFirst = objTmp.idxNodeFirst();
    ScaFES::Ntuple<int, DIM> idxNodeLast = objTmp.idxNodeLast();

    for (int ii = 0; ii < 5; ++ii) {
        borderwidth = ii;
        idxNodeLast[0] += borderwidth;
        objExpected = ScaFES::Grid<DIM>(idxNodeFirst,
                                        idxNodeLast,
                                        objTmp.coordNodeFirst(),
                                        objTmp.coordNodeLast());
        EXPECT_TRUE( (objTmp.extend(direction, borderwidth) == objExpected) );
        idxNodeFirst = objTmp.idxNodeFirst();
        idxNodeLast = objTmp.idxNodeLast();
    }
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, ExtendLeft)
{
    const int DIM = 3;
    int borderwidth;
    int direction = 1;
    ScaFES::Grid<DIM> objTmp = this->mObjA;
    ScaFES::Grid<DIM> objExpected;
    ScaFES::Ntuple<int, DIM> idxNodeFirst = objTmp.idxNodeFirst();
    ScaFES::Ntuple<int, DIM> idxNodeLast = objTmp.idxNodeLast();

    for (int ii = 0; ii < 5; ++ii) {
        borderwidth = ii;
        idxNodeFirst[0] -= borderwidth;
        objExpected = ScaFES::Grid<DIM>(idxNodeFirst,
                                        idxNodeLast,
                                        objTmp.coordNodeFirst(),
                                        objTmp.coordNodeLast());
        EXPECT_TRUE( (objTmp.extend(direction, borderwidth) == objExpected) );
        idxNodeFirst = objTmp.idxNodeFirst();
        idxNodeLast = objTmp.idxNodeLast();
    }
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, ExtendTop)
{
    const int DIM = 3;
    int borderwidth;
    int direction = 2;
    ScaFES::Grid<DIM> objTmp = this->mObjA;
    ScaFES::Grid<DIM> objExpected;
    ScaFES::Ntuple<int, DIM> idxNodeFirst = objTmp.idxNodeFirst();
    ScaFES::Ntuple<int, DIM> idxNodeLast = objTmp.idxNodeLast();

    for (int ii = 0; ii < 5; ++ii) {
        borderwidth = ii;
        idxNodeLast[1] += borderwidth;
        objExpected = ScaFES::Grid<DIM>(idxNodeFirst,
                                        idxNodeLast,
                                        objTmp.coordNodeFirst(),
                                        objTmp.coordNodeLast());
        EXPECT_TRUE( (objTmp.extend(direction, borderwidth) == objExpected) );
        idxNodeFirst = objTmp.idxNodeFirst();
        idxNodeLast = objTmp.idxNodeLast();
    }
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, ExtendBottom)
{
    const int DIM = 3;
    int borderwidth;
    int direction = 3;
    ScaFES::Grid<DIM> objTmp = this->mObjA;
    ScaFES::Grid<DIM> objExpected;
    ScaFES::Ntuple<int, DIM> idxNodeFirst = objTmp.idxNodeFirst();
    ScaFES::Ntuple<int, DIM> idxNodeLast = objTmp.idxNodeLast();

    for (int ii = 0; ii < 5; ++ii) {
        borderwidth = ii;
        idxNodeFirst[1] -= borderwidth;
        objExpected = ScaFES::Grid<DIM>(idxNodeFirst,
                                        idxNodeLast,
                                        objTmp.coordNodeFirst(),
                                        objTmp.coordNodeLast());
        EXPECT_TRUE( (objTmp.extend(direction, borderwidth) == objExpected) );
        idxNodeFirst = objTmp.idxNodeFirst();
        idxNodeLast = objTmp.idxNodeLast();
    }
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, ExtendFront)
{
    const int DIM = 3;
    int borderwidth;
    int direction = 4;
    ScaFES::Grid<DIM> objTmp = this->mObjA;
    ScaFES::Grid<DIM> objExpected;
    ScaFES::Ntuple<int, DIM> idxNodeFirst = objTmp.idxNodeFirst();
    ScaFES::Ntuple<int, DIM> idxNodeLast = objTmp.idxNodeLast();

    for (int ii = 0; ii < 5; ++ii) {
        borderwidth = ii;
        idxNodeLast[2] += borderwidth;
        objExpected = ScaFES::Grid<DIM>(idxNodeFirst,
                                        idxNodeLast,
                                        objTmp.coordNodeFirst(),
                                        objTmp.coordNodeLast());
        EXPECT_TRUE( (objTmp.extend(direction, borderwidth) == objExpected) );
        idxNodeFirst = objTmp.idxNodeFirst();
        idxNodeLast = objTmp.idxNodeLast();
    }
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, ExtendBack)
{
    const int DIM = 3;
    int borderwidth;
    int direction = 5;
    ScaFES::Grid<DIM> objTmp = this->mObjA;
    ScaFES::Grid<DIM> objExpected;
    ScaFES::Ntuple<int, DIM> idxNodeFirst = objTmp.idxNodeFirst();
    ScaFES::Ntuple<int, DIM> idxNodeLast = objTmp.idxNodeLast();

    for (int ii = 0; ii < 5; ++ii) {
        borderwidth = ii;
        idxNodeFirst[2] -= borderwidth;
        objExpected = ScaFES::Grid<DIM>(idxNodeFirst,
                                        idxNodeLast,
                                        objTmp.coordNodeFirst(),
                                        objTmp.coordNodeLast());
        EXPECT_TRUE( (objTmp.extend(direction, borderwidth) == objExpected) );
        idxNodeFirst = objTmp.idxNodeFirst();
        idxNodeLast = objTmp.idxNodeLast();
    }
}
/*******************************************************************************
 * TESTS OF INTERNAL CLASS ITERATOR: GETTER METHODS.
 ******************************************************************************/
TYPED_TEST_P(GridTest, IteratorGlobalNode)
{
    const int DIM = 3;
    ScaFES::Grid<DIM> objTmp = this->mObjA;
    ScaFES::Grid<DIM>::iterator iter(&objTmp, objTmp.idxNodeFirst());
    ScaFES::Ntuple<int, DIM> idxNodeActual;
    ScaFES::Ntuple<int, DIM> idxNodeExpected;
    ScaFES::Ntuple<int, DIM> idxNode;

    for (int ii = 0; ii < 2; ++ii) {
        for (int jj = 0; jj < 2; ++jj) {
            for (int kk = 0; kk < 2; ++kk) {
                if (1 == ii) {
                    idxNode.elem(0) = objTmp.idxNodeLast().elem(0);
                } else {
                    idxNode.elem(0) = objTmp.idxNodeFirst().elem(0);
                }

                if (1 == jj) {
                    idxNode.elem(1) = objTmp.idxNodeLast().elem(1);
                } else {
                    idxNode.elem(1) = objTmp.idxNodeFirst().elem(1);
                }

                if (1 == kk) {
                    idxNode.elem(2) = objTmp.idxNodeLast().elem(2);
                } else {
                    idxNode.elem(2) = objTmp.idxNodeFirst().elem(2);
                }

                iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
                idxNodeExpected = idxNode;
                idxNodeActual = iter.idxNode();
                EXPECT_TRUE(idxNodeExpected == idxNodeActual);
            }
        }
    }

    // Arbitrary node.
    idxNode[0] = objTmp.idxNodeLast().elem(0) - 1;
    idxNode[1] = objTmp.idxNodeLast().elem(1) - 1;
    idxNode[2] = objTmp.idxNodeLast().elem(2) - 1;
    idxNodeExpected = idxNode;
    iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, IteratorLocNode)
{
    const int DIM = 3;
    ScaFES::Grid<DIM> objTmp = this->mObjA;
    ScaFES::Grid<DIM>::iterator iter(&objTmp, objTmp.idxNodeFirst());
    ScaFES::Ntuple<int, DIM> idxNode;
    unsigned long int idxNodeActual;
    unsigned long int idxNodeExpected;

    // Checking all corner points of the cuboid.
    for (int ii = 0; ii < 2; ++ii) {
        for (int jj = 0; jj < 2; ++jj) {
            for (int kk = 0; kk < 2; ++kk) {
                if (1 == ii) {
                    idxNode.elem(0) = objTmp.idxNodeLast().elem(0);
                } else {
                    idxNode.elem(0) = objTmp.idxNodeFirst().elem(0);
                }

                if (1 == jj) {
                    idxNode.elem(1) = objTmp.idxNodeLast().elem(1);
                } else {
                    idxNode.elem(1) = objTmp.idxNodeFirst().elem(1);
                }

                if (1 == kk) {
                    idxNode.elem(2) = objTmp.idxNodeLast().elem(2);
                } else {
                    idxNode.elem(2) = objTmp.idxNodeFirst().elem(2);
                }

                idxNodeExpected = objTmp.idxNodeTuple2Scalar(idxNode);
                iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
                idxNodeActual = iter.idxScalarNode();
                EXPECT_TRUE(idxNodeExpected == idxNodeActual);
            }
        }
    }

    // Arbitrary node.
    idxNode[0] = objTmp.idxNodeLast().elem(0) - 1;
    idxNode[1] = objTmp.idxNodeLast().elem(1) - 1;
    idxNode[2] = objTmp.idxNodeLast().elem(2) - 1;
    idxNodeExpected = objTmp.idxNodeTuple2Scalar(idxNode);
    iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
    idxNodeActual = iter.idxScalarNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);
}

/*******************************************************************************
 * TESTS OF INTERNAL CLASS ITERATOR: WORK METHODS.
 ******************************************************************************/
TYPED_TEST_P(GridTest, IteratorOperatorStar)
{
    const int DIM = 3;
    ScaFES::Grid<DIM> objTmp = this->mObjA;
    ScaFES::Grid<DIM>::iterator iter(&objTmp, objTmp.idxNodeFirst());
    ScaFES::Ntuple<int, DIM> idxNodeActual;
    ScaFES::Ntuple<int, DIM> idxNodeExpected;
    ScaFES::Ntuple<int, DIM> idxNode;

    for (int ii = 0; ii < 2; ++ii) {
        for (int jj = 0; jj < 2; ++jj) {
            for (int kk = 0; kk < 2; ++kk) {
                if (1 == ii) {
                    idxNode.elem(0) = objTmp.idxNodeLast().elem(0);
                } else {
                    idxNode.elem(0) = objTmp.idxNodeFirst().elem(0);
                }

                if (1 == jj) {
                    idxNode.elem(1) = objTmp.idxNodeLast().elem(1);
                } else {
                    idxNode.elem(1) = objTmp.idxNodeFirst().elem(1);
                }

                if (1 == kk) {
                    idxNode.elem(2) = objTmp.idxNodeLast().elem(2);
                } else {
                    idxNode.elem(2) = objTmp.idxNodeFirst().elem(2);
                }

                idxNodeExpected = idxNode;
                iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
                idxNodeActual = iter.operator * ();
                EXPECT_TRUE(idxNodeExpected == idxNodeActual);
            }
        }
    }

    // Arbitrary node.
    idxNode[0] = objTmp.idxNodeLast().elem(0) - 1;
    idxNode[1] = objTmp.idxNodeLast().elem(1) - 1;
    idxNode[2] = objTmp.idxNodeLast().elem(2) - 1;
    idxNodeExpected = idxNode;
    iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
    idxNodeActual = iter.operator * ();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, IteratorOperatorPlusPlusPrefix)
{
    const int DIM = 3;
    ScaFES::Grid<DIM> objTmp = this->mObjA;
    ScaFES::Grid<DIM>::iterator iter(&objTmp, objTmp.idxNodeFirst());
    ScaFES::Ntuple<int, DIM> idxNodeActual;
    ScaFES::Ntuple<int, DIM> idxNodeExpected;
    ScaFES::Ntuple<int, DIM> idxNode;

    idxNode = objTmp.idxNodeFirst();
    idxNodeExpected[0] = objTmp.idxNodeFirst().elem(0) + 1;
    idxNodeExpected[1] = objTmp.idxNodeFirst().elem(1);
    idxNodeExpected[2] = objTmp.idxNodeFirst().elem(2);
    iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++();
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeLast().elem(0);
    idxNode[1] = objTmp.idxNodeFirst().elem(1);
    idxNode[2] = objTmp.idxNodeFirst().elem(2);
    idxNodeExpected[0] = objTmp.idxNodeFirst().elem(0);
    idxNodeExpected[1] = objTmp.idxNodeFirst().elem(1) + 1;
    idxNodeExpected[2] = objTmp.idxNodeFirst().elem(2);
    iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++();
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeFirst().elem(0);
    idxNode[1] = objTmp.idxNodeLast().elem(1);
    idxNode[2] = objTmp.idxNodeFirst().elem(2);
    idxNodeExpected[0] = objTmp.idxNodeFirst().elem(0) + 1;
    idxNodeExpected[1] = objTmp.idxNodeLast().elem(1);
    idxNodeExpected[2] = objTmp.idxNodeFirst().elem(2);
    iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++();
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeLast().elem(0);
    idxNode[1] = objTmp.idxNodeLast().elem(1);
    idxNode[2] = objTmp.idxNodeFirst().elem(2);
    idxNodeExpected[0] = objTmp.idxNodeFirst().elem(0);
    idxNodeExpected[1] = objTmp.idxNodeFirst().elem(1);
    idxNodeExpected[2] = objTmp.idxNodeFirst().elem(2) + 1;
    iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++();
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeFirst().elem(0);
    idxNode[1] = objTmp.idxNodeFirst().elem(1);
    idxNode[2] = objTmp.idxNodeLast().elem(2);
    idxNodeExpected[0] = objTmp.idxNodeFirst().elem(0) + 1;
    idxNodeExpected[1] = objTmp.idxNodeFirst().elem(1);
    idxNodeExpected[2] = objTmp.idxNodeLast().elem(2);
    iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++();
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeLast().elem(0);
    idxNode[1] = objTmp.idxNodeFirst().elem(1);
    idxNode[2] = objTmp.idxNodeLast().elem(2);
    idxNodeExpected[0] = objTmp.idxNodeFirst().elem(0);
    idxNodeExpected[1] = objTmp.idxNodeFirst().elem(1) + 1;
    idxNodeExpected[2] = objTmp.idxNodeLast().elem(2);
    iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++();
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeFirst().elem(0);
    idxNode[1] = objTmp.idxNodeLast().elem(1);
    idxNode[2] = objTmp.idxNodeLast().elem(2);
    idxNodeExpected[0] = objTmp.idxNodeFirst().elem(0) + 1;
    idxNodeExpected[1] = objTmp.idxNodeLast().elem(1);
    idxNodeExpected[2] = objTmp.idxNodeLast().elem(2);
    iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++();
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeLast().elem(0);
    idxNode[1] = objTmp.idxNodeLast().elem(1);
    idxNode[2] = objTmp.idxNodeLast().elem(2);
    idxNodeExpected[0] = objTmp.idxNodeFirst().elem(0);
    idxNodeExpected[1] = objTmp.idxNodeFirst().elem(1);
    idxNodeExpected[2] = objTmp.idxNodeLast().elem(2) + 1;
    iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++();
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    // Arbitrary node.
    idxNode[0] = objTmp.idxNodeLast().elem(0) - 1;
    idxNode[1] = objTmp.idxNodeLast().elem(1) - 1;
    idxNode[2] = objTmp.idxNodeLast().elem(2) - 1;
    idxNodeExpected[0] = objTmp.idxNodeLast().elem(0);
    idxNodeExpected[1] = objTmp.idxNodeLast().elem(1) - 1;
    idxNodeExpected[2] = objTmp.idxNodeLast().elem(2) - 1;
    iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++();
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, IteratorOperatorPlusPlusPostfix)
{
    const int DIM = 3;
    ScaFES::Grid<DIM> objTmp = this->mObjA;
    ScaFES::Grid<DIM>::iterator iter(&objTmp, objTmp.idxNodeFirst());
    ScaFES::Ntuple<int, DIM> idxNodeActual;
    ScaFES::Ntuple<int, DIM> idxNodeExpected;
    ScaFES::Ntuple<int, DIM> idxNode;

    idxNode = objTmp.idxNodeFirst();
    idxNodeExpected[0] = objTmp.idxNodeFirst().elem(0) + 1;
    idxNodeExpected[1] = objTmp.idxNodeFirst().elem(1);
    idxNodeExpected[2] = objTmp.idxNodeFirst().elem(2);
    iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++(1);
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeLast().elem(0);
    idxNode[1] = objTmp.idxNodeFirst().elem(1);
    idxNode[2] = objTmp.idxNodeFirst().elem(2);
    idxNodeExpected[0] = objTmp.idxNodeFirst().elem(0);
    idxNodeExpected[1] = objTmp.idxNodeFirst().elem(1) + 1;
    idxNodeExpected[2] = objTmp.idxNodeFirst().elem(2);
    iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++(1);
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeFirst().elem(0);
    idxNode[1] = objTmp.idxNodeLast().elem(1);
    idxNode[2] = objTmp.idxNodeFirst().elem(2);
    idxNodeExpected[0] = objTmp.idxNodeFirst().elem(0) + 1;
    idxNodeExpected[1] = objTmp.idxNodeLast().elem(1);
    idxNodeExpected[2] = objTmp.idxNodeFirst().elem(2);
    iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++(1);
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeLast().elem(0);
    idxNode[1] = objTmp.idxNodeLast().elem(1);
    idxNode[2] = objTmp.idxNodeFirst().elem(2);
    idxNodeExpected[0] = objTmp.idxNodeFirst().elem(0);
    idxNodeExpected[1] = objTmp.idxNodeFirst().elem(1);
    idxNodeExpected[2] = objTmp.idxNodeFirst().elem(2) + 1;
    iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++(1);
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeFirst().elem(0);
    idxNode[1] = objTmp.idxNodeFirst().elem(1);
    idxNode[2] = objTmp.idxNodeLast().elem(2);
    idxNodeExpected[0] = objTmp.idxNodeFirst().elem(0) + 1;
    idxNodeExpected[1] = objTmp.idxNodeFirst().elem(1);
    idxNodeExpected[2] = objTmp.idxNodeLast().elem(2);
    iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++(1);
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeLast().elem(0);
    idxNode[1] = objTmp.idxNodeFirst().elem(1);
    idxNode[2] = objTmp.idxNodeLast().elem(2);
    idxNodeExpected[0] = objTmp.idxNodeFirst().elem(0);
    idxNodeExpected[1] = objTmp.idxNodeFirst().elem(1) + 1;
    idxNodeExpected[2] = objTmp.idxNodeLast().elem(2);
    iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++();
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeFirst().elem(0);
    idxNode[1] = objTmp.idxNodeLast().elem(1);
    idxNode[2] = objTmp.idxNodeLast().elem(2);
    idxNodeExpected[0] = objTmp.idxNodeFirst().elem(0) + 1;
    idxNodeExpected[1] = objTmp.idxNodeLast().elem(1);
    idxNodeExpected[2] = objTmp.idxNodeLast().elem(2);
    iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++(1);
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeLast().elem(0);
    idxNode[1] = objTmp.idxNodeLast().elem(1);
    idxNode[2] = objTmp.idxNodeLast().elem(2);
    idxNodeExpected[0] = objTmp.idxNodeFirst().elem(0);
    idxNodeExpected[1] = objTmp.idxNodeFirst().elem(1);
    idxNodeExpected[2] = objTmp.idxNodeLast().elem(2) + 1;
    iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++(1);
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    // Arbitrary node.
    idxNode[0] = objTmp.idxNodeLast().elem(0) - 1;
    idxNode[1] = objTmp.idxNodeLast().elem(1) - 1;
    idxNode[2] = objTmp.idxNodeLast().elem(2) - 1;
    idxNodeExpected[0] = objTmp.idxNodeLast().elem(0);
    idxNodeExpected[1] = objTmp.idxNodeLast().elem(1) - 1;
    idxNodeExpected[2] = objTmp.idxNodeLast().elem(2) - 1;
    iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++(1);
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, IteratorOperatorNotEquals)
{
    const int DIM = 3;
    ScaFES::Grid<DIM> objTmp = this->mObjA;
    ScaFES::Ntuple<int, DIM> idxNode1;
    ScaFES::Ntuple<int, DIM> idxNode2;
    idxNode1 = objTmp.idxNodeFirst();
    idxNode2 = objTmp.idxNodeLast();
    ScaFES::Grid<DIM>::iterator iter1(&(objTmp), idxNode1);
    ScaFES::Grid<DIM>::iterator iter2(&(objTmp), idxNode2);
    EXPECT_TRUE(iter1.operator != (iter2));
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest, IteratorOperatorSmallerThan)
{
    const int DIM = 3;
    ScaFES::Grid<DIM> objTmp = this->mObjA;
    ScaFES::Ntuple<int, DIM> idxNode1;
    ScaFES::Ntuple<int, DIM> idxNode2;
    idxNode1 = objTmp.idxNodeFirst();
    idxNode2 = objTmp.idxNodeLast();
    ScaFES::Grid<DIM>::iterator iter1(&(objTmp), idxNode1);
    ScaFES::Grid<DIM>::iterator iter2(&(objTmp), idxNode2);
    EXPECT_TRUE(iter1.operator < (iter2));
}

/*******************************************************************************
 * TESTS OF FREE METHODS.
 ******************************************************************************/
TYPED_TEST_P(GridTest, Swap)
{
    const int DIM = 3;
    ScaFES::Grid<DIM> mTmpA = this->mObjA;
    ScaFES::Grid<DIM> mTmpB = this->mObjB;
    ScaFES::swap(this->mObjA, this->mObjB);
    EXPECT_TRUE(mTmpA == this->mObjB);
    EXPECT_TRUE(mTmpB == this->mObjA);
}

/*******************************************************************************
 ******************************************************************************/
REGISTER_TYPED_TEST_CASE_P(GridTest,
                           GlobalNodeFirst,
                           GlobalNodeLast,
                           GlobalNodeFirstInOneDirection,
                           GlobalNodeLastInOneDirection,
                           Nnodes,
                           NnodesInOneDirection,
                           GridSize,
                           GridSizeInOneDirection,
                           CoordNodeFirst,
                           CoordNodeLast,
                           IsNumberedStyleC,
                           OperatorEqualsComponentwise,
                           OperatorEquals,
                           OperatorUnequals,
                           ConstructorDefault,
                           ConstructorWithParameters,
                           CopyConstructor,
                           AssignmentOperator,
                           IdxNodeTuple2Scalar,
//                            LocNode2IdxNode,
                           Coordinates,
                           MethodValid,
                           MethodInside,
                           UnionWith,
                           IntersectWith,
                           ExtendRight,
                           ExtendLeft,
                           ExtendTop,
                           ExtendBottom,
                           ExtendFront,
                           ExtendBack,
                           IteratorGlobalNode,
                           IteratorLocNode,
                           IteratorOperatorStar,
                           IteratorOperatorPlusPlusPrefix,
                           IteratorOperatorPlusPlusPostfix,
                           IteratorOperatorNotEquals,
                           IteratorOperatorSmallerThan,
                           Swap
                          );

/*******************************************************************************
 ******************************************************************************/
typedef ::testing::Types<int> MyGridTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(My, GridTest, MyGridTestTypes);

} // end namespace
