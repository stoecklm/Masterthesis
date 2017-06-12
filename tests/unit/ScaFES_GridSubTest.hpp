#ifndef SCAFES_GRIDSUBTEST_HPP_
#define SCAFES_GRIDSUBTEST_HPP_

#include "gtest/gtest.h"

#include "ScaFES_Ntuple.hpp"
#include "ScaFES_GridSub.hpp"

namespace ScaFES_test
{

/*******************************************************************************
 ******************************************************************************/
/** Test class for the class 'Grid'. */
template<typename DD>
class GridSubTest : public testing::Test
{
    public:
        /*----------------------------------------------------------------------
        | LIFE CYCLE METHODS.
        ----------------------------------------------------------------------*/
        /** Constructor. */
        GridSubTest();

        /** Copy constructor. */
        GridSubTest(GridSubTest<DD> const&) = delete;

        /** Assignment operator. */
        GridSubTest& operator= (GridSubTest<DD> const&) = delete;

        /** Destructor. */
        ~GridSubTest() = default;

    public:
        /*----------------------------------------------------------------------
        | MEMBER VARIABLES.
        ----------------------------------------------------------------------*/
        // Independent member variables(A).
        /** Start index of grid 'a'. */
        ScaFES::Ntuple<int, 3> mIdxNodeFirstSubA;

        /** End index of grid 'a'. */
        ScaFES::Ntuple<int, 3> mIdxNodeLastSubA;

        /** Start index of the base grid of 'a'. */
        ScaFES::Ntuple<int, 3> mIdxNodeFirstBaseA;

        /** End index of the base grid of 'a'. */
        ScaFES::Ntuple<int, 3> mIdxNodeLastBaseA;

        /** Coordinates of first node number of grid 'a'. */
        ScaFES::Ntuple<double, 3> mCoordNodeFirstA;

        /** Coordinates of last node number of grid 'a'. */
        ScaFES::Ntuple<double, 3> mCoordNodeLastA;

        /*--------------------------------------------------------------------*/
        // Variables depending on above variables(A).
        /** Returns the number of nodes in all dimensions. */
        ScaFES::Ntuple<int, 3> mNnodesSubA;

        /** Returns the number of nodes in all dimensions of the base grid. */
        ScaFES::Ntuple<int, 3> mNnodesBaseA;

        /** Returns the grid size in all dimensions. */
        ScaFES::Ntuple<double, 3> mGridSizeA;

        /** Returns the first local node number in the grid. */
        unsigned long int mIdxNodeFirstA;

        /** Returns the last local node number in the grid. */
        unsigned long int mIdxNodeLastA;

        /** Returns the number of all nodes of the grid. */
        int mNnodesTotalA;

        /*--------------------------------------------------------------------*/
        // Independent member variables(B)
        /** Start index of grid 'b'. */
        ScaFES::Ntuple<int, 3> mIdxNodeFirstB;

        /** End index of grid 'b'. */
        ScaFES::Ntuple<int, 3> mIdxNodeLastB;

        /*--------------------------------------------------------------------*/
        /** Index null. */
        ScaFES::Ntuple<int, 3> mIndexNull;

        /** Arbitrary element 'a'. */
        ScaFES::GridSub<3> mObjA;

        /** 1:1 copy of element 'a'. */
        ScaFES::GridSub<3> mObjCopyA;

        /** Another arbitrary element 'b'. */
        ScaFES::GridSub<3> mObjB;

        /** Null element. */
        ScaFES::GridSub<3> mObjNull;

        /** Element for lesser and greater tests. */
        ScaFES::GridSub<3> mObjLess;

        /** Element for lesser and greater tests. */
        ScaFES::GridSub<3> mObjGreater;

        /** Element for lesser and greater tests. */
        ScaFES::GridSub<3> mObjUndetermined;
};

TYPED_TEST_CASE_P(GridSubTest);

/*******************************************************************************
 * LIFE CYCLE METHODS.
 ******************************************************************************/
template<typename DD>
inline GridSubTest<DD>::GridSubTest()
    : mIdxNodeFirstSubA(1, 2, 3)
    , mIdxNodeLastSubA(6, 14, 18)
    , mIdxNodeFirstBaseA(0, 1, 2)
    , mIdxNodeLastBaseA(9, 16, 25)
    , mCoordNodeFirstA(1.1, 2.2, 3.3)
    , mCoordNodeLastA(2.3, 5.4, 7.5)
    , mNnodesSubA(6, 13, 16)
    , mNnodesBaseA(10, 16, 24)
    , mGridSizeA(1.2 / 9, 3.2 / 15, 4.2 / 23)
    , mIdxNodeFirstA(171)
    , mIdxNodeLastA(2696)
    , mNnodesTotalA(1248)
    , mIndexNull(0)
    , mObjA(mIdxNodeFirstSubA,
            mIdxNodeLastSubA,
            mIdxNodeFirstBaseA,
            mIdxNodeLastBaseA,
            mCoordNodeFirstA,
            mCoordNodeLastA)
    , mObjCopyA(mIdxNodeFirstSubA,
                mIdxNodeLastSubA,
                mIdxNodeFirstBaseA,
                mIdxNodeLastBaseA,
                mCoordNodeFirstA,
                mCoordNodeLastA)
    , mObjB(ScaFES::Ntuple<int, 3>(2, 3, 4),
            ScaFES::Ntuple<int, 3>(4, 5, 8),
            mIdxNodeFirstBaseA,
            mIdxNodeLastBaseA,
            mCoordNodeFirstA,
            mCoordNodeLastA)
    , mObjNull(mIndexNull,
               mIndexNull,
               mIdxNodeFirstBaseA,
               mIdxNodeLastBaseA,
               mCoordNodeFirstA,
               mCoordNodeLastA)
    , mObjLess(ScaFES::Ntuple<int, 3>(1, 2, 3),
               ScaFES::Ntuple<int, 3>(5, 6, 7),
               mIdxNodeFirstBaseA,
               mIdxNodeLastBaseA,
               mCoordNodeFirstA,
               mCoordNodeLastA)
    , mObjGreater(ScaFES::Ntuple<int, 3>(1, 2, 3),
                  ScaFES::Ntuple<int, 3>(8, 9, 10),
                  mIdxNodeFirstBaseA,
                  mIdxNodeLastBaseA,
                  mCoordNodeFirstA,
                  mCoordNodeLastA)
    , mObjUndetermined(ScaFES::Ntuple<int, 3>(2, 3, 4),
                       ScaFES::Ntuple<int, 3>(5, 6, 7),
                       mIdxNodeFirstBaseA,
                       mIdxNodeLastBaseA,
                       mCoordNodeFirstA,
                       mCoordNodeLastA)
{ }

} // End of namespace. //
#endif
