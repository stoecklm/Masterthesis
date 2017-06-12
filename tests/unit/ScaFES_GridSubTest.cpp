#include <algorithm>

#include "gtest/gtest.h"

#include "ScaFES_GridSubTest.hpp"

#include "ScaFES_Ntuple.hpp"
#include "ScaFES_GridSub.hpp"

namespace ScaFES_test
{

/*******************************************************************************
 * TESTS OF GETTER METHODS.
 ******************************************************************************/
TYPED_TEST_P(GridSubTest, GetGlobalNodeFirstSub)
{
    EXPECT_TRUE(this->mIdxNodeFirstSubA == this->mObjA.idxNodeFirstSub());
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridSubTest, GetGlobalNodeLastSub)
{
    EXPECT_TRUE(this->mIdxNodeLastSubA == this->mObjA.idxNodeLastSub());
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridSubTest, GetGlobalNodeFirstInOneDirectionSub)
{
    const int DIM = 3;

    for (int ii = 0; ii < DIM; ++ii) {
        EXPECT_TRUE(this->mIdxNodeFirstSubA.elem(ii)
                    == this->mObjA.idxNodeFirstSub(ii));
    }
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridSubTest, GetGlobalNodeLastInOneDirectionSub)
{
    const int DIM = 3;

    for (int ii = 0; ii < DIM; ++ii) {
        EXPECT_TRUE(this->mIdxNodeLastSubA.elem(ii)
                    == this->mObjA.idxNodeLastSub(ii));
    }
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridSubTest, GetIdxNodeFirst)
{
    EXPECT_TRUE(this->mIdxNodeFirstA == this->mObjA.idxScalarNodeFirst());
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridSubTest, GetIdxNodeLast)
{
    EXPECT_TRUE(this->mIdxNodeLastA == this->mObjA.idxScalarNodeLast());
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridSubTest, GetNnodesSub)
{
    EXPECT_TRUE(this->mNnodesSubA == this->mObjA.nNodesSub());
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridSubTest, GetNnodesSubInOneDirection)
{
    const int DIM = 3;

    for (int ii = 0; ii < DIM; ++ii) {
        EXPECT_TRUE(this->mNnodesSubA.elem(ii) == this->mObjA.nNodesSub(ii));
    }
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridSubTest, GetNnodesTotal)
{
    EXPECT_TRUE(this->mNnodesTotalA == this->mObjA.nNodesTotal());
}

/*******************************************************************************
 * TESTS OF COMPARISON METHODS.
 ******************************************************************************/
TYPED_TEST_P(GridSubTest, EqualsOperatorComponentwise)
{
    EXPECT_TRUE(static_cast< ScaFES::Grid<3> >(this->mObjCopyA)
                == static_cast< ScaFES::Grid<3> >(this->mObjA));
    EXPECT_TRUE(this->mObjCopyA.nNodesSub() == this->mObjA.nNodesSub());
    EXPECT_TRUE(this->mObjCopyA.nNodesTotal() == this->mObjA.nNodesTotal());
    EXPECT_TRUE(this->mObjCopyA.idxNodeFirstSub() ==
                this->mObjA.idxNodeFirstSub());
    EXPECT_TRUE(this->mObjCopyA.idxNodeLastSub() == this->mObjA.idxNodeLastSub());
    EXPECT_TRUE(this->mObjCopyA.idxNodeFirst() == this->mObjA.idxNodeFirst());
    EXPECT_TRUE(this->mObjCopyA.idxNodeLast() == this->mObjA.idxNodeLast());
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridSubTest, EqualsOperator)
{
    EXPECT_TRUE(this->mObjA == this->mObjA);
    EXPECT_TRUE(this->mObjCopyA == this->mObjA);
    EXPECT_FALSE(this->mObjLess == this->mObjGreater);
    EXPECT_FALSE(this->mObjGreater == this->mObjLess);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridSubTest, UnequalsOperator)
{
    EXPECT_FALSE(this->mObjA != this->mObjA);
    EXPECT_FALSE(this->mObjCopyA != this->mObjA);
    EXPECT_TRUE(this->mObjLess != this->mObjGreater);
    EXPECT_TRUE(this->mObjGreater != this->mObjLess);
}

/*******************************************************************************
 * TESTS OF LIFE CYCLE METHODS.
 ******************************************************************************/
TYPED_TEST_P(GridSubTest, ConstructorWithParameters)
{
    const int DIM = 3;

    ScaFES::GridSub<DIM> objTmp(this->mIdxNodeFirstSubA,
                                this->mIdxNodeLastSubA,
                                ScaFES::Grid<3>(this->mIdxNodeFirstBaseA,
                                                this->mIdxNodeLastBaseA,
                                                this->mCoordNodeFirstA,
                                                this->mCoordNodeLastA) );
    EXPECT_TRUE(this->mObjA == objTmp);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridSubTest, ConstructorWithParameters2)
{
    const int DIM = 3;

    ScaFES::GridSub<DIM> objTmp(this->mIdxNodeFirstSubA,
                                this->mIdxNodeLastSubA,
                                this->mIdxNodeFirstBaseA,
                                this->mIdxNodeLastBaseA,
                                this->mCoordNodeFirstA,
                                this->mCoordNodeLastA);
    EXPECT_TRUE(this->mObjA == objTmp);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridSubTest, ConstructorDefault)
{
    const int DIM = 3;
    ScaFES::GridSub<DIM> objTmp;
    EXPECT_TRUE((ScaFES::Ntuple<double, DIM>(1.0) == objTmp.gridsize()));
    EXPECT_TRUE((ScaFES::Ntuple<int, DIM>(0) == objTmp.idxNodeFirst()));
    EXPECT_TRUE((ScaFES::Ntuple<int, DIM>(1) == objTmp.idxNodeLast()));
    EXPECT_TRUE((ScaFES::Ntuple<int, DIM>(2) == objTmp.nNodes()));
    EXPECT_TRUE((ScaFES::Ntuple<double, DIM>(0.0) == objTmp.coordNodeFirst()));
    EXPECT_TRUE((ScaFES::Ntuple<double, DIM>(1.0) == objTmp.coordNodeLast()));
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridSubTest, CopyConstructor)
{
    const int DIM = 3;

    ScaFES::GridSub<DIM> objTmp = this->mObjA;
    EXPECT_TRUE(this->mObjA == objTmp);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridSubTest, AssignmentOperator)
{
    const int DIM = 3;

    ScaFES::GridSub<DIM> objTmp;
    objTmp = this->mObjA;
    EXPECT_TRUE(this->mObjA == objTmp);
}

/*******************************************************************************
 * TESTS OF WORK METHODS.
 ******************************************************************************/
TYPED_TEST_P(GridSubTest, MethodValid)
{
    EXPECT_TRUE(this->mObjA.valid());
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridSubTest, MethodInside)
{
    const int DIM = 3;
    EXPECT_TRUE(this->mObjA.inside(this->mIdxNodeFirstSubA));
    EXPECT_TRUE(this->mObjA.inside(this->mIdxNodeLastSubA));

    ScaFES::Ntuple<int, DIM> shift(1);
    EXPECT_FALSE(this->mObjA.inside(this->mIdxNodeFirstSubA - shift));
    EXPECT_FALSE(this->mObjA.inside(this->mIdxNodeLastSubA + shift));
    EXPECT_TRUE(this->mObjA.inside(this->mIdxNodeFirstSubA + shift));
    EXPECT_TRUE(this->mObjA.inside(this->mIdxNodeLastSubA - shift));
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridSubTest, UnionWith)
{
    const int DIM = 3;
    ScaFES::Ntuple<int, DIM> idxNode1(0, 0, 0);
    ScaFES::Ntuple<int, DIM> idxNode2(1, 1, 1);
    ScaFES::Ntuple<int, DIM> idxNode3(2, 2, 2);
    ScaFES::Ntuple<double, DIM> coordNode1(1.0, 2.0, 3.0);
    ScaFES::Ntuple<double, DIM> coordNode2(1.0, 1.0, 1.0);
    ScaFES::Ntuple<double, DIM> coordNode3(2.0, 2.0, 2.0);
    ScaFES::Grid<DIM> objTmpA(idxNode1,
                              idxNode2,
                              coordNode1,
                              coordNode2);
    ScaFES::Grid<DIM> objTmpB(idxNode2,
                              idxNode3,
                              coordNode1,
                              coordNode2);
    ScaFES::Grid<DIM> objExpected(idxNode1,
                                  idxNode3,
                                  coordNode1,
                                  coordNode3);
    EXPECT_TRUE(objExpected == objTmpA.unionWith(objTmpB));
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridSubTest, IntersectWith)
{
    const int DIM = 3;
    ScaFES::Ntuple<int, DIM> idxNode1(2, 4, 2);
    ScaFES::Ntuple<int, DIM> idxNode2(2, 3, 1);
    ScaFES::Ntuple<int, DIM> idxNode3(2, 3, 1);
    ScaFES::Ntuple<double, DIM> coordNode1(1.0, 2.0, 3.0);
    ScaFES::Ntuple<double, DIM> coordNode2(1.0, 1.0, 1.0);
    ScaFES::Ntuple<double, DIM> coordNode3(2.0, 2.0, 2.0);
    ScaFES::Grid<DIM> objTmpA(idxNode1,
                              idxNode2,
                              coordNode1,
                              coordNode2);
    ScaFES::Grid<DIM> objTmpB(idxNode2,
                              idxNode3,
                              coordNode2,
                              coordNode3);
    ScaFES::Grid<DIM> objExpected(idxNode1,
                                  idxNode3,
                                  coordNode1,
                                  coordNode3);
    EXPECT_TRUE(objExpected == objTmpA.intersectWith(objTmpB));
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridSubTest, ExtendRight)
{
    const int DIM = 3;
    int borderwidth = 3;
    int direction = 0;
    ScaFES::Ntuple<int, DIM> idxNodeAdded(borderwidth, 0, 0);
    ScaFES::Ntuple<int, DIM> idxNodeFirstSub(this->mIdxNodeFirstSubA);
    ScaFES::Ntuple<int, DIM> idxNodeLastSub(this->mIdxNodeLastSubA +
                                             idxNodeAdded);
    ScaFES::GridSub<DIM> objTmp(idxNodeFirstSub, idxNodeLastSub,
                                this->mIdxNodeFirstBaseA,
                                this->mIdxNodeLastBaseA,
                                this->mCoordNodeFirstA, this->mCoordNodeLastA);

    EXPECT_TRUE( ((this->mObjA.extend(direction, borderwidth)) == objTmp));
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridSubTest, ExtendLeft)
{
    const int DIM = 3;
    int borderwidth = 3;
    int direction = 1;
    ScaFES::Ntuple<int, DIM> idxNodeAdded(borderwidth, 0, 0);
    ScaFES::Ntuple<int, DIM> idxNodeFirstSub(this->mIdxNodeFirstSubA -
                                              idxNodeAdded);
    ScaFES::Ntuple<int, DIM> idxNodeLastSub(this->mIdxNodeLastSubA);
    ScaFES::GridSub<DIM> objTmp(idxNodeFirstSub, idxNodeLastSub,
                                this->mIdxNodeFirstBaseA,
                                this->mIdxNodeLastBaseA,
                                this->mCoordNodeFirstA, this->mCoordNodeLastA);
    EXPECT_TRUE( (this->mObjA.extend(direction, borderwidth) == objTmp));

}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridSubTest, ExtendTop)
{
    const int DIM = 3;
    int borderwidth = 3;
    int direction = 2;
    ScaFES::Ntuple<int, DIM> idxNodeAdded(0, borderwidth, 0);
    ScaFES::Ntuple<int, DIM> idxNodeFirstSub(this->mIdxNodeFirstSubA);
    ScaFES::Ntuple<int, DIM> idxNodeLastSub(this->mIdxNodeLastSubA +
                                             idxNodeAdded);
    ScaFES::GridSub<DIM> objTmp(idxNodeFirstSub, idxNodeLastSub,
                                this->mIdxNodeFirstBaseA,
                                this->mIdxNodeLastBaseA,
                                this->mCoordNodeFirstA, this->mCoordNodeLastA);
    EXPECT_TRUE( (this->mObjA.extend(direction, borderwidth) == objTmp));
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridSubTest, ExtendBottom)
{
    const int DIM = 3;
    int borderwidth = 3;
    int direction = 3;
    ScaFES::Ntuple<int, DIM> idxNodeAdded(0, borderwidth, 0);
    ScaFES::Ntuple<int, DIM> idxNodeFirstSub(this->mIdxNodeFirstSubA -
                                              idxNodeAdded);
    ScaFES::Ntuple<int, DIM> idxNodeLastSub(this->mIdxNodeLastSubA);
    ScaFES::GridSub<DIM> objTmp(idxNodeFirstSub, idxNodeLastSub,
                                this->mIdxNodeFirstBaseA,
                                this->mIdxNodeLastBaseA,
                                this->mCoordNodeFirstA, this->mCoordNodeLastA);
    EXPECT_TRUE( (this->mObjA.extend(direction, borderwidth) == objTmp));

}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridSubTest, ExtendFront)
{
    const int DIM = 3;
    int borderwidth = 3;
    int direction = 4;
    ScaFES::Ntuple<int, DIM> idxNodeAdded(0, 0, borderwidth);
    ScaFES::Ntuple<int, DIM> idxNodeFirstSub(this->mIdxNodeFirstSubA);
    ScaFES::Ntuple<int, DIM> idxNodeLastSub(this->mIdxNodeLastSubA +
                                             idxNodeAdded);
    ScaFES::GridSub<DIM> objTmp(idxNodeFirstSub, idxNodeLastSub,
                                this->mIdxNodeFirstBaseA,
                                this->mIdxNodeLastBaseA,
                                this->mCoordNodeFirstA, this->mCoordNodeLastA);
    EXPECT_TRUE( (this->mObjA.extend(direction, borderwidth) == objTmp));

}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridSubTest, ExtendBack)
{
    const int DIM = 3;
    int borderwidth = 3;
    int direction = 5;
    ScaFES::Ntuple<int, DIM> idxNodeAdded(0, 0, borderwidth);
    ScaFES::Ntuple<int, DIM> idxNodeFirst(this->mIdxNodeFirstSubA
                                           - idxNodeAdded);
    ScaFES::Ntuple<int, DIM> idxNodeLast(this->mIdxNodeLastSubA);
    ScaFES::GridSub<DIM> objTmp(idxNodeFirst, idxNodeLast,
                                this->mIdxNodeFirstBaseA,
                                this->mIdxNodeLastBaseA,
                                this->mCoordNodeFirstA, this->mCoordNodeLastA);
    EXPECT_TRUE( (this->mObjA.extend(direction, borderwidth) == objTmp));
}

/*******************************************************************************
 * TESTS OF INTERNAL CLASS ITERATOR: GETTER METHODS.
 ******************************************************************************/
TYPED_TEST_P(GridSubTest, IteratorGlobalNode)
{
    const int DIM = 3;
    ScaFES::GridSub<DIM> objTmp = this->mObjA;
    ScaFES::GridSub<DIM>::iterator iter(&objTmp, objTmp.idxNodeFirstSub());
    ScaFES::Ntuple<int, DIM> idxNodeActual;
    ScaFES::Ntuple<int, DIM> idxNodeExpected;
    ScaFES::Ntuple<int, DIM> idxNode;

    for (int ii = 0; ii < 2; ++ii) {
        for (int jj = 0; jj < 2; ++jj) {
            for (int kk = 0; kk < 2; ++kk) {
                if (1 == ii) {
                    idxNode.elem(0) = objTmp.idxNodeLastSub().elem(0);
                } else {
                    idxNode.elem(0) = objTmp.idxNodeFirstSub().elem(0);
                }

                if (1 == jj) {
                    idxNode.elem(1) = objTmp.idxNodeLastSub().elem(1);
                } else {
                    idxNode.elem(1) = objTmp.idxNodeFirstSub().elem(1);
                }

                if (1 == kk) {
                    idxNode.elem(2) = objTmp.idxNodeLastSub().elem(2);
                } else {
                    idxNode.elem(2) = objTmp.idxNodeFirstSub().elem(2);
                }

                iter = ScaFES::GridSub<DIM>::iterator(&(objTmp), idxNode);
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
    iter = ScaFES::GridSub<DIM>::iterator(&(objTmp), idxNode);
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridSubTest, IteratorIdxNode)
{
    const int DIM = 3;
    ScaFES::GridSub<DIM> objTmp = this->mObjA;
    ScaFES::GridSub<DIM>::iterator iter(&objTmp, objTmp.idxNodeFirstSub());
    ScaFES::Ntuple<int, DIM> idxNode;
    unsigned long int idxNodeActual;
    unsigned long int idxNodeExpected;

    // Checking all corner points of the cuboid.
    for (int ii = 0; ii < 2; ++ii) {
        for (int jj = 0; jj < 2; ++jj) {
            for (int kk = 0; kk < 2; ++kk) {
                if (1 == ii) {
                    idxNode.elem(0) = objTmp.idxNodeLastSub().elem(0);
                } else {
                    idxNode.elem(0) = objTmp.idxNodeFirstSub().elem(0);
                }

                if (1 == jj) {
                    idxNode.elem(1) = objTmp.idxNodeLastSub().elem(1);
                } else {
                    idxNode.elem(1) = objTmp.idxNodeFirstSub().elem(1);
                }

                if (1 == kk) {
                    idxNode.elem(2) = objTmp.idxNodeLastSub().elem(2);
                } else {
                    idxNode.elem(2) = objTmp.idxNodeFirstSub().elem(2);
                }

                idxNodeExpected = objTmp.idxNodeTuple2Scalar(idxNode);
                iter = ScaFES::GridSub<DIM>::iterator(&(objTmp), idxNode);
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
    iter = ScaFES::GridSub<DIM>::iterator(&(objTmp), idxNode);
    idxNodeActual = iter.idxScalarNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);
}

/*******************************************************************************
 * TESTS OF INTERNAL CLASS ITERATOR: WORK METHODS.
 ******************************************************************************/
TYPED_TEST_P(GridSubTest, IteratorOperatorStar)
{
    const int DIM = 3;
    ScaFES::GridSub<DIM> objTmp = this->mObjA;
    ScaFES::GridSub<DIM>::iterator iter(&objTmp, objTmp.idxNodeFirstSub());
    ScaFES::Ntuple<int, DIM> idxNodeActual;
    ScaFES::Ntuple<int, DIM> idxNodeExpected;
    ScaFES::Ntuple<int, DIM> idxNode;

    for (int ii = 0; ii < 2; ++ii) {
        for (int jj = 0; jj < 2; ++jj) {
            for (int kk = 0; kk < 2; ++kk) {
                if (1 == ii) {
                    idxNode.elem(0) = objTmp.idxNodeLastSub().elem(0);
                } else {
                    idxNode.elem(0) = objTmp.idxNodeFirstSub().elem(0);
                }

                if (1 == jj) {
                    idxNode.elem(1) = objTmp.idxNodeLastSub().elem(1);
                } else {
                    idxNode.elem(1) = objTmp.idxNodeFirstSub().elem(1);
                }

                if (1 == kk) {
                    idxNode.elem(2) = objTmp.idxNodeLastSub().elem(2);
                } else {
                    idxNode.elem(2) = objTmp.idxNodeFirstSub().elem(2);
                }

                idxNodeExpected = idxNode;
                iter = ScaFES::GridSub<DIM>::iterator(&(objTmp), idxNode);
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
    iter = ScaFES::GridSub<DIM>::iterator(&(objTmp), idxNode);
    idxNodeActual = iter.operator * ();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridSubTest, IteratorOperatorPlusPlusPrefix)
{
    const int DIM = 3;
    ScaFES::GridSub<DIM> objTmp = this->mObjA;
    ScaFES::GridSub<DIM>::iterator iter(&objTmp, objTmp.idxNodeFirstSub());
    ScaFES::Ntuple<int, DIM> idxNodeActual;
    ScaFES::Ntuple<int, DIM> idxNodeExpected;
    ScaFES::Ntuple<int, DIM> idxNode;

    idxNode = objTmp.idxNodeFirstSub();
    idxNodeExpected[0] = objTmp.idxNodeFirstSub().elem(0) + 1;
    idxNodeExpected[1] = objTmp.idxNodeFirstSub().elem(1);
    idxNodeExpected[2] = objTmp.idxNodeFirstSub().elem(2);
    iter = ScaFES::GridSub<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++();
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeLastSub().elem(0);
    idxNode[1] = objTmp.idxNodeFirstSub().elem(1);
    idxNode[2] = objTmp.idxNodeFirstSub().elem(2);
    idxNodeExpected[0] = objTmp.idxNodeFirstSub().elem(0);
    idxNodeExpected[1] = objTmp.idxNodeFirstSub().elem(1) + 1;
    idxNodeExpected[2] = objTmp.idxNodeFirstSub().elem(2);
    iter = ScaFES::GridSub<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++();
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeFirstSub().elem(0);
    idxNode[1] = objTmp.idxNodeLastSub().elem(1);
    idxNode[2] = objTmp.idxNodeFirstSub().elem(2);
    idxNodeExpected[0] = objTmp.idxNodeFirstSub().elem(0) + 1;
    idxNodeExpected[1] = objTmp.idxNodeLastSub().elem(1);
    idxNodeExpected[2] = objTmp.idxNodeFirstSub().elem(2);
    iter = ScaFES::GridSub<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++();
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeLastSub().elem(0);
    idxNode[1] = objTmp.idxNodeLastSub().elem(1);
    idxNode[2] = objTmp.idxNodeFirstSub().elem(2);
    idxNodeExpected[0] = objTmp.idxNodeFirstSub().elem(0);
    idxNodeExpected[1] = objTmp.idxNodeFirstSub().elem(1);
    idxNodeExpected[2] = objTmp.idxNodeFirstSub().elem(2) + 1;
    iter = ScaFES::GridSub<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++();
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeFirstSub().elem(0);
    idxNode[1] = objTmp.idxNodeFirstSub().elem(1);
    idxNode[2] = objTmp.idxNodeLastSub().elem(2);
    idxNodeExpected[0] = objTmp.idxNodeFirstSub().elem(0) + 1;
    idxNodeExpected[1] = objTmp.idxNodeFirstSub().elem(1);
    idxNodeExpected[2] = objTmp.idxNodeLastSub().elem(2);
    iter = ScaFES::GridSub<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++();
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeLastSub().elem(0);
    idxNode[1] = objTmp.idxNodeFirstSub().elem(1);
    idxNode[2] = objTmp.idxNodeLastSub().elem(2);
    idxNodeExpected[0] = objTmp.idxNodeFirstSub().elem(0);
    idxNodeExpected[1] = objTmp.idxNodeFirstSub().elem(1) + 1;
    idxNodeExpected[2] = objTmp.idxNodeLastSub().elem(2);
    iter = ScaFES::GridSub<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++();
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeFirstSub().elem(0);
    idxNode[1] = objTmp.idxNodeLastSub().elem(1);
    idxNode[2] = objTmp.idxNodeLastSub().elem(2);
    idxNodeExpected[0] = objTmp.idxNodeFirstSub().elem(0) + 1;
    idxNodeExpected[1] = objTmp.idxNodeLastSub().elem(1);
    idxNodeExpected[2] = objTmp.idxNodeLastSub().elem(2);
    iter = ScaFES::GridSub<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++();
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeLastSub().elem(0);
    idxNode[1] = objTmp.idxNodeLastSub().elem(1);
    idxNode[2] = objTmp.idxNodeLastSub().elem(2);
    idxNodeExpected[0] = objTmp.idxNodeFirstSub().elem(0);
    idxNodeExpected[1] = objTmp.idxNodeFirstSub().elem(1);
    idxNodeExpected[2] = objTmp.idxNodeLastSub().elem(2) + 1;
    iter = ScaFES::GridSub<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++();
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    // Arbitrary node.
    idxNode[0] = objTmp.idxNodeLastSub().elem(0) - 1;
    idxNode[1] = objTmp.idxNodeLastSub().elem(1) - 1;
    idxNode[2] = objTmp.idxNodeLastSub().elem(2) - 1;
    idxNodeExpected[0] = objTmp.idxNodeLastSub().elem(0);
    idxNodeExpected[1] = objTmp.idxNodeLastSub().elem(1) - 1;
    idxNodeExpected[2] = objTmp.idxNodeLastSub().elem(2) - 1;
    iter = ScaFES::GridSub<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++();
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridSubTest, IteratorOperatorPlusPlusPostfix)
{
    const int DIM = 3;
    ScaFES::GridSub<DIM> objTmp = this->mObjA;
    ScaFES::GridSub<DIM>::iterator iter(&objTmp, objTmp.idxNodeFirstSub());
    ScaFES::Ntuple<int, DIM> idxNodeActual;
    ScaFES::Ntuple<int, DIM> idxNodeExpected;
    ScaFES::Ntuple<int, DIM> idxNode;

    idxNode = objTmp.idxNodeFirstSub();
    idxNodeExpected[0] = objTmp.idxNodeFirstSub().elem(0) + 1;
    idxNodeExpected[1] = objTmp.idxNodeFirstSub().elem(1);
    idxNodeExpected[2] = objTmp.idxNodeFirstSub().elem(2);
    iter = ScaFES::GridSub<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++(1);
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeLastSub().elem(0);
    idxNode[1] = objTmp.idxNodeFirstSub().elem(1);
    idxNode[2] = objTmp.idxNodeFirstSub().elem(2);
    idxNodeExpected[0] = objTmp.idxNodeFirstSub().elem(0);
    idxNodeExpected[1] = objTmp.idxNodeFirstSub().elem(1) + 1;
    idxNodeExpected[2] = objTmp.idxNodeFirstSub().elem(2);
    iter = ScaFES::GridSub<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++(1);
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeFirstSub().elem(0);
    idxNode[1] = objTmp.idxNodeLastSub().elem(1);
    idxNode[2] = objTmp.idxNodeFirstSub().elem(2);
    idxNodeExpected[0] = objTmp.idxNodeFirstSub().elem(0) + 1;
    idxNodeExpected[1] = objTmp.idxNodeLastSub().elem(1);
    idxNodeExpected[2] = objTmp.idxNodeFirstSub().elem(2);
    iter = ScaFES::GridSub<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++(1);
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeLastSub().elem(0);
    idxNode[1] = objTmp.idxNodeLastSub().elem(1);
    idxNode[2] = objTmp.idxNodeFirstSub().elem(2);
    idxNodeExpected[0] = objTmp.idxNodeFirstSub().elem(0);
    idxNodeExpected[1] = objTmp.idxNodeFirstSub().elem(1);
    idxNodeExpected[2] = objTmp.idxNodeFirstSub().elem(2) + 1;
    iter = ScaFES::GridSub<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++(1);
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeFirstSub().elem(0);
    idxNode[1] = objTmp.idxNodeFirstSub().elem(1);
    idxNode[2] = objTmp.idxNodeLastSub().elem(2);
    idxNodeExpected[0] = objTmp.idxNodeFirstSub().elem(0) + 1;
    idxNodeExpected[1] = objTmp.idxNodeFirstSub().elem(1);
    idxNodeExpected[2] = objTmp.idxNodeLastSub().elem(2);
    iter = ScaFES::GridSub<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++(1);
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeLastSub().elem(0);
    idxNode[1] = objTmp.idxNodeFirstSub().elem(1);
    idxNode[2] = objTmp.idxNodeLastSub().elem(2);
    idxNodeExpected[0] = objTmp.idxNodeFirstSub().elem(0);
    idxNodeExpected[1] = objTmp.idxNodeFirstSub().elem(1) + 1;
    idxNodeExpected[2] = objTmp.idxNodeLastSub().elem(2);
    iter = ScaFES::GridSub<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++(1);
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeFirstSub().elem(0);
    idxNode[1] = objTmp.idxNodeLastSub().elem(1);
    idxNode[2] = objTmp.idxNodeLastSub().elem(2);
    idxNodeExpected[0] = objTmp.idxNodeFirstSub().elem(0) + 1;
    idxNodeExpected[1] = objTmp.idxNodeLastSub().elem(1);
    idxNodeExpected[2] = objTmp.idxNodeLastSub().elem(2);
    iter = ScaFES::GridSub<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++(1);
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    idxNode[0] = objTmp.idxNodeLastSub().elem(0);
    idxNode[1] = objTmp.idxNodeLastSub().elem(1);
    idxNode[2] = objTmp.idxNodeLastSub().elem(2);
    idxNodeExpected[0] = objTmp.idxNodeFirstSub().elem(0);
    idxNodeExpected[1] = objTmp.idxNodeFirstSub().elem(1);
    idxNodeExpected[2] = objTmp.idxNodeLastSub().elem(2) + 1;
    iter = ScaFES::GridSub<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++(1);
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);

    // Arbitrary node.
    idxNode[0] = objTmp.idxNodeLastSub().elem(0) - 1;
    idxNode[1] = objTmp.idxNodeLastSub().elem(1) - 1;
    idxNode[2] = objTmp.idxNodeLastSub().elem(2) - 1;
    idxNodeExpected[0] = objTmp.idxNodeLastSub().elem(0);
    idxNodeExpected[1] = objTmp.idxNodeLastSub().elem(1) - 1;
    idxNodeExpected[2] = objTmp.idxNodeLastSub().elem(2) - 1;
    iter = ScaFES::GridSub<DIM>::iterator(&(objTmp), idxNode);
    iter.operator++(1);
    idxNodeActual = iter.idxNode();
    EXPECT_TRUE(idxNodeExpected == idxNodeActual);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridSubTest, IteratorOperatorNotEquals)
{
    const int DIM = 3;
    ScaFES::GridSub<DIM> objTmp = this->mObjA;
    ScaFES::Ntuple<int, DIM> idxNode1;
    ScaFES::Ntuple<int, DIM> idxNode2;
    idxNode1 = objTmp.idxNodeFirstSub();
    idxNode2 = objTmp.idxNodeLastSub();
    ScaFES::GridSub<DIM>::iterator iter1(&(objTmp), idxNode1);
    ScaFES::GridSub<DIM>::iterator iter2(&(objTmp), idxNode2);
    EXPECT_TRUE(iter1.operator != (iter2));
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(GridSubTest, IteratorOperatorSmallerThan)
{
    const int DIM = 3;
    ScaFES::GridSub<DIM> objTmp = this->mObjA;
    ScaFES::Ntuple<int, DIM> idxNode1;
    ScaFES::Ntuple<int, DIM> idxNode2;
    idxNode1 = objTmp.idxNodeFirstSub();
    idxNode2 = objTmp.idxNodeLastSub();
    ScaFES::Grid<DIM>::iterator iter1(&(objTmp), idxNode1);
    ScaFES::Grid<DIM>::iterator iter2(&(objTmp), idxNode2);
    EXPECT_TRUE(iter1.operator < (iter2));
}

/*******************************************************************************
 * TESTS OF FREE METHODS.
 ******************************************************************************/
TYPED_TEST_P(GridSubTest, Swap)
{
    const int DIM = 3;
    ScaFES::GridSub<DIM> mTmpA = this->mObjA;
    ScaFES::GridSub<DIM> mTmpB = this->mObjB;
    ScaFES::swap(this->mObjA, this->mObjB);
    EXPECT_TRUE(mTmpA == this->mObjB);
    EXPECT_TRUE(mTmpB == this->mObjA);
}

/*******************************************************************************
 ******************************************************************************/
REGISTER_TYPED_TEST_CASE_P(GridSubTest,
                           GetGlobalNodeFirstSub,
                           GetGlobalNodeLastSub,
                           GetGlobalNodeFirstInOneDirectionSub,
                           GetGlobalNodeLastInOneDirectionSub,
                           GetIdxNodeFirst,
                           GetIdxNodeLast,
                           GetNnodesSub,
                           GetNnodesSubInOneDirection,
                           GetNnodesTotal,
                           EqualsOperatorComponentwise,
                           EqualsOperator,
                           UnequalsOperator,
                           ConstructorWithParameters,
                           ConstructorWithParameters2,
                           ConstructorDefault,
                           CopyConstructor,
                           AssignmentOperator,
                           MethodValid,
                           MethodInside,
                           UnionWith,
                           IntersectWith,
                           ExtendRight,
                           ExtendLeft,
                           ExtendFront,
                           ExtendBack,
                           ExtendBottom,
                           ExtendTop,
                           IteratorGlobalNode,
                           IteratorIdxNode,
                           IteratorOperatorStar,
                           IteratorOperatorPlusPlusPrefix,
                           IteratorOperatorPlusPlusPostfix,
                           IteratorOperatorNotEquals,
                           IteratorOperatorSmallerThan,
                           Swap
                          );

/*******************************************************************************
 ******************************************************************************/
typedef ::testing::Types<int> MyGridSubTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(My, GridSubTest, MyGridSubTestTypes);

} // end namespace
