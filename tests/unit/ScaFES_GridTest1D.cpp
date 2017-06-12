#include <algorithm>

#include "gtest/gtest.h"

#include "ScaFES_GridTest1D.hpp"

#include "ScaFES_Ntuple.hpp"
#include "ScaFES_Grid.hpp"

namespace ScaFES_test
{

const double EPS = 2.2e-14;

/*******************************************************************************
 * TESTS OF GETTER METHODS.
 ******************************************************************************/
TYPED_TEST_P(GridTest1D, GlobalNodeFirst)
{
    EXPECT_TRUE(this->mIdxNodeFirstA == this->mObjA.idxNodeFirst());
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest1D, GlobalNodeLast)
{
    EXPECT_TRUE(this->mIdxNodeLastA == this->mObjA.idxNodeLast());
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest1D, GlobalNodeFirstInOneDirection)
{
    const int DIM = 1;

    for (int ii = 0; ii < DIM; ++ii) {
        EXPECT_TRUE(this->mIdxNodeFirstA.elem(ii)
                    == this->mObjA.idxNodeFirst(ii));
    }
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest1D, GlobalNodeLastInOneDirection)
{
    const int DIM = 1;

    for (int ii = 0; ii < DIM; ++ii) {
        EXPECT_TRUE(this->mIdxNodeLastA.elem(ii)
                    == this->mObjA.idxNodeLast(ii));
    }
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest1D, Nnodes)
{
    EXPECT_TRUE(this->mNnodesA == this->mObjA.nNodes());
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest1D, NnodesInOneDirection)
{
    const int DIM = 1;

    for (int ii = 0; ii < DIM; ++ii) {
        EXPECT_TRUE(this->mNnodesA.elem(ii) == this->mObjA.nNodes(ii));
    }
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest1D, GridSize)
{
    const int DIM = 1;
    ScaFES::Ntuple<double, DIM> tmpGridSize = this->mObjA.gridsize();

    for (int ii = 0; ii < DIM; ++ii) {
        EXPECT_TRUE(::fabs(this->mGridSizeA.elem(ii)
                           - tmpGridSize.elem(ii)) < EPS);
    }
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest1D, GridSizeInOneDirection)
{
    const int DIM = 1;

    for (int ii = 0; ii < DIM; ++ii) {
        EXPECT_TRUE(::fabs(this->mGridSizeA.elem(ii)
                           - this->mObjA.gridsize(ii)) < EPS);
    }
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest1D, CoordNodeFirst)
{
    const int DIM = 1;
    ScaFES::Ntuple<double, DIM> tmpCoordNodeFirst
        = this->mObjA.coordNodeFirst();

    for (int ii = 0; ii < DIM; ++ii) {
        EXPECT_TRUE(this->mCoordNodeFirstA.elem(ii)
                    == tmpCoordNodeFirst.elem(ii));
    }
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest1D, CoordNodeLast)
{
    const int DIM = 1;
    ScaFES::Ntuple<double, DIM> tmpCoordNodeLast = this->mObjA.coordNodeLast();

    for (int ii = 0; ii < DIM; ++ii) {
        EXPECT_TRUE(::fabs(this->mCoordNodeLastA.elem(ii)
                           - tmpCoordNodeLast.elem(ii)) < EPS);
    }
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest1D, IsNumberedStyleC)
{
    EXPECT_TRUE(this->mIsNumberedStyleCA == this->mObjA.isNumberedStyleC());
}

/*******************************************************************************
 * TESTS OF operator==().
 ******************************************************************************/
TYPED_TEST_P(GridTest1D, OperatorEqualsComponentwise)
{
    const int DIM = 1;
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
TYPED_TEST_P(GridTest1D, OperatorEquals)
{
    EXPECT_TRUE(this->mObjA == this->mObjA);
    EXPECT_TRUE(this->mObjCopyA == this->mObjA);
    EXPECT_FALSE(this->mObjLess == this->mObjGreater);
    EXPECT_FALSE(this->mObjGreater == this->mObjLess);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest1D, OperatorUnequals)
{
    EXPECT_FALSE(this->mObjA != this->mObjA);
    EXPECT_FALSE(this->mObjCopyA != this->mObjA);
    EXPECT_TRUE(this->mObjLess != this->mObjGreater);
    EXPECT_TRUE(this->mObjGreater != this->mObjLess);
}

/*******************************************************************************
 * TESTS OF LIFE CYCLE METHODS.
 ******************************************************************************/
TYPED_TEST_P(GridTest1D, ConstructorWithParameters)
{
    const int DIM = 1;
    ScaFES::Grid<DIM> objTmp(this->mIdxNodeFirstA,
                             this->mIdxNodeLastA,
                             this->mCoordNodeFirstA,
                             this->mCoordNodeLastA);
    EXPECT_TRUE(this->mObjA == objTmp);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest1D, ConstructorDefault)
{
    const int DIM = 1;
    ScaFES::Grid<DIM> objTmp;
    EXPECT_TRUE((ScaFES::Ntuple<double, DIM>(1.0) == objTmp.gridsize()));
    EXPECT_TRUE((ScaFES::Ntuple<int, DIM>(0) == objTmp.idxNodeFirst()));
    EXPECT_TRUE((ScaFES::Ntuple<int, DIM>(1) == objTmp.idxNodeLast()));
    EXPECT_TRUE((ScaFES::Ntuple<int, DIM>(2) == objTmp.nNodes()));
    EXPECT_TRUE((ScaFES::Ntuple<double, DIM>(0.0) == objTmp.coordNodeFirst()));
    EXPECT_TRUE((ScaFES::Ntuple<double, DIM>(1.0) == objTmp.coordNodeLast()));
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest1D, CopyConstructor)
{
    const int DIM = 1;
    ScaFES::Grid<DIM> objTmp = this->mObjA;
    EXPECT_TRUE(this->mObjA == objTmp);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest1D, AssignmentOperator)
{
    const int DIM = 1;
    ScaFES::Grid<DIM> objTmp;
    objTmp = this->mObjA;
    EXPECT_TRUE(this->mObjA == objTmp);
}

/*******************************************************************************
 * TESTS OF WORK METHODS.
 ******************************************************************************/
TYPED_TEST_P(GridTest1D, IdxNodeTuple2Scalar)
{
    const int DIM = 1;
    ScaFES::Grid<DIM> objTmp = this->mObjA;
    unsigned long int idxNodeExpected;
    unsigned long int idxNodeActual;
    ScaFES::Ntuple<int, DIM> idxNode;

    idxNode[0] = objTmp.idxNodeFirst().elem(0);
    idxNodeActual = objTmp.idxNodeTuple2Scalar(idxNode);
    idxNodeExpected = 0;
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeLast().elem(0);
    idxNodeActual = objTmp.idxNodeTuple2Scalar(idxNode);
    idxNodeExpected = objTmp.nNodes(0) - 1;
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    // Arbitrary node.
    idxNode[0] = objTmp.idxNodeLast().elem(0) - 1;
    idxNodeActual = objTmp.idxNodeTuple2Scalar(idxNode);
    idxNodeExpected = objTmp.nNodes(0) - 2;
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest1D, Coordinates)
{
    const int DIM = 1;
    ScaFES::Grid<DIM> objTmp = this->mObjA;
    ScaFES::Ntuple<double, DIM> coordNodeExpected;
    ScaFES::Ntuple<double, DIM> coordNodeActual;
    ScaFES::Ntuple<int, DIM> idxNode;

    idxNode[0] = objTmp.idxNodeFirst().elem(0);
    coordNodeActual = objTmp.coordinates(idxNode);
    coordNodeExpected = ScaFES::Ntuple<double, DIM>(
                            objTmp.coordNodeFirst().elem(0) );
    EXPECT_TRUE(coordNodeExpected == coordNodeActual);

    idxNode[0] = objTmp.idxNodeLast().elem(0);
    coordNodeActual = objTmp.coordinates(idxNode);
    coordNodeExpected = ScaFES::Ntuple<double, DIM>(
                            objTmp.coordNodeLast().elem(0) );
    EXPECT_TRUE(coordNodeExpected == coordNodeActual);

    // Arbitrary node.
    idxNode[0] = objTmp.idxNodeLast().elem(0) - 1;
    coordNodeActual = objTmp.coordinates(idxNode);
    coordNodeExpected = ScaFES::Ntuple<double, DIM>(
                           objTmp.coordNodeLast().elem(0) - objTmp.gridsize(0));
    EXPECT_TRUE(coordNodeExpected == coordNodeActual);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest1D, MethodInside)
{
    const int DIM = 1;
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
TYPED_TEST_P(GridTest1D, UnionWith)
{
    const int DIM = 1;
    ScaFES::Ntuple<int, DIM> idxNode1(1);
    ScaFES::Ntuple<int, DIM> idxNode2(2);
    ScaFES::Ntuple<int, DIM> idxNode3(3);
    ScaFES::Ntuple<double, DIM> coordNode1(1.0);
    ScaFES::Ntuple<double, DIM> coordNode2(2.0);
    ScaFES::Ntuple<double, DIM> coordNode3(3.0);
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
TYPED_TEST_P(GridTest1D, IntersectWith)
{
    const int DIM = 1;
    ScaFES::Ntuple<int, DIM> idxNode1(1);
    ScaFES::Ntuple<int, DIM> idxNode2(2);
    ScaFES::Ntuple<int, DIM> idxNode3(3);
    ScaFES::Ntuple<double, DIM> coordNode1(1.0);
    ScaFES::Ntuple<double, DIM> coordNode2(2.0);
    ScaFES::Ntuple<double, DIM> coordNode3(3.0);
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
TYPED_TEST_P(GridTest1D, MethodValid)
{
    const int DIM = 1;
    ScaFES::Grid<DIM> objTmp = this->mObjA;
    EXPECT_TRUE(objTmp.valid());
}
/*---------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest1D, ExtendRight)
{
    const int DIM = 1;
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
TYPED_TEST_P(GridTest1D, ExtendLeft)
{
    const int DIM = 1;
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

/*******************************************************************************
 * TESTS OF INTERNAL CLASS ITERATOR: GETTER METHODS.
 ******************************************************************************/
TYPED_TEST_P(GridTest1D, IteratorGlobalNode)
{
    const int DIM = 1;
    ScaFES::Grid<DIM> objTmp = this->mObjA;
    ScaFES::Grid<DIM>::iterator iter(&objTmp, objTmp.idxNodeFirst());
    ScaFES::Ntuple<int, DIM> idxNodeActual;
    ScaFES::Ntuple<int, DIM> idxNodeExpected;
    ScaFES::Ntuple<int, DIM> idxNode;

    for (int ii = 0; ii < 2; ++ii) {
                if (1 == ii) {
                    idxNode.elem(0) = objTmp.idxNodeLast().elem(0);
                } else {
                    idxNode.elem(0) = objTmp.idxNodeFirst().elem(0);
                }

                iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
                idxNodeExpected = idxNode;
                idxNodeActual = iter.idxNode();
                EXPECT_TRUE(idxNodeExpected == idxNodeActual);
    }

    // Arbitrary node.
    idxNode[0] = objTmp.idxNodeLast().elem(0) - 1;
    idxNodeExpected = idxNode;
    iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest1D, IteratorLocNode)
{
    const int DIM = 1;
    ScaFES::Grid<DIM> objTmp = this->mObjA;
    ScaFES::Grid<DIM>::iterator iter(&objTmp, objTmp.idxNodeFirst());
    ScaFES::Ntuple<int, DIM> idxNode;
    unsigned long int idxNodeActual;
    unsigned long int idxNodeExpected;

    // Checking all corner points of the cuboid.
    for (int ii = 0; ii < 2; ++ii) {
                if (1 == ii) {
                    idxNode.elem(0) = objTmp.idxNodeLast().elem(0);
                } else {
                    idxNode.elem(0) = objTmp.idxNodeFirst().elem(0);
                }

                idxNodeExpected = objTmp.idxNodeTuple2Scalar(idxNode);
                iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
                idxNodeActual = iter.idxScalarNode();
                EXPECT_TRUE(idxNodeExpected == idxNodeActual);
    }

    // Arbitrary node.
    idxNode[0] = objTmp.idxNodeLast().elem(0) - 1;
    idxNodeExpected = objTmp.idxNodeTuple2Scalar(idxNode);
    iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
    idxNodeActual = iter.idxScalarNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);
}

/*******************************************************************************
 * TESTS OF INTERNAL CLASS ITERATOR: WORK METHODS.
 ******************************************************************************/
TYPED_TEST_P(GridTest1D, IteratorOperatorStar)
{
    const int DIM = 1;
    ScaFES::Grid<DIM> objTmp = this->mObjA;
    ScaFES::Grid<DIM>::iterator iter(&objTmp, objTmp.idxNodeFirst());
    ScaFES::Ntuple<int, DIM> idxNodeActual;
    ScaFES::Ntuple<int, DIM> idxNodeExpected;
    ScaFES::Ntuple<int, DIM> idxNode;

    for (int ii = 0; ii < 2; ++ii) {
                if (1 == ii) {
                    idxNode.elem(0) = objTmp.idxNodeLast().elem(0);
                } else {
                    idxNode.elem(0) = objTmp.idxNodeFirst().elem(0);
                }

                idxNodeExpected = idxNode;
                iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
                idxNodeActual = iter.operator * ();
                EXPECT_TRUE(idxNodeExpected == idxNodeActual);
    }

    // Arbitrary node.
    idxNode[0] = objTmp.idxNodeLast().elem(0) - 1;
    idxNodeExpected = idxNode;
    iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
    idxNodeActual = iter.operator * ();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest1D, IteratorOperatorPlusPlusPrefix)
{
    const int DIM = 1;
    ScaFES::Grid<DIM> objTmp = this->mObjA;
    ScaFES::Grid<DIM>::iterator iter(&objTmp, objTmp.idxNodeFirst());
    ScaFES::Ntuple<int, DIM> idxNodeActual;
    ScaFES::Ntuple<int, DIM> idxNodeExpected;
    ScaFES::Ntuple<int, DIM> idxNode;

    idxNode = objTmp.idxNodeFirst();
    idxNodeExpected[0] = objTmp.idxNodeFirst().elem(0) + 1;
    iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++();
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeLast().elem(0);
    idxNodeExpected[0] = objTmp.idxNodeLast().elem(0) + 1;
    iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++();
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    // Arbitrary node.
    idxNode[0] = objTmp.idxNodeLast().elem(0) - 1;
    idxNodeExpected[0] = objTmp.idxNodeLast().elem(0);
    iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++();
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest1D, IteratorOperatorPlusPlusPostfix)
{
    const int DIM = 1;
    ScaFES::Grid<DIM> objTmp = this->mObjA;
    ScaFES::Grid<DIM>::iterator iter(&objTmp, objTmp.idxNodeFirst());
    ScaFES::Ntuple<int, DIM> idxNodeActual;
    ScaFES::Ntuple<int, DIM> idxNodeExpected;
    ScaFES::Ntuple<int, DIM> idxNode;

    idxNode = objTmp.idxNodeFirst();
    idxNodeExpected[0] = objTmp.idxNodeFirst().elem(0) + 1;
    iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++(1);
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeLast().elem(0);
    idxNodeExpected[0] = objTmp.idxNodeLast().elem(0) + 1;
    iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++(1);
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    // Arbitrary node.
    idxNode[0] = objTmp.idxNodeLast().elem(0) - 1;
    idxNodeExpected[0] = objTmp.idxNodeLast().elem(0);
    iter = ScaFES::Grid<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++(1);
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridTest1D, IteratorOperatorNotEquals)
{
    const int DIM = 1;
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
TYPED_TEST_P(GridTest1D, IteratorOperatorSmallerThan)
{
    const int DIM = 1;
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
TYPED_TEST_P(GridTest1D, Swap)
{
    const int DIM = 1;
    ScaFES::Grid<DIM> mTmpA = this->mObjA;
    ScaFES::Grid<DIM> mTmpB = this->mObjB;
    ScaFES::swap(this->mObjA, this->mObjB);
    EXPECT_TRUE(mTmpA == this->mObjB);
    EXPECT_TRUE(mTmpB == this->mObjA);
}

/*******************************************************************************
 ******************************************************************************/
REGISTER_TYPED_TEST_CASE_P(GridTest1D,
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
                           Coordinates,
                           MethodValid,
                           MethodInside,
                           UnionWith,
                           IntersectWith,
                           ExtendRight,
                           ExtendLeft,
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
typedef ::testing::Types<int> MyGridTest1DTypes;
INSTANTIATE_TYPED_TEST_CASE_P(My, GridTest1D, MyGridTest1DTypes);

} // end namespace
