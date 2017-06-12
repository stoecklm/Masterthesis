#ifndef SCAFES_GRIDTEST_HPP_
#define SCAFES_GRIDTEST_HPP_

#include "gtest/gtest.h"

#include "ScaFES_Ntuple.hpp"
#include "ScaFES_Grid.hpp"

namespace ScaFES_test
{

/*******************************************************************************
 ******************************************************************************/
/**
 * Test class for the class 'Grid'.
 */
template<typename DD>
class GridTest : public testing::Test
{
    public:
        /*----------------------------------------------------------------------
        | LIFE CYCLE METHODS.
        ----------------------------------------------------------------------*/
        /** Constructor. */
        GridTest();

        /** Copy constructor. */
        GridTest(GridTest<DD> const&) = delete;

        /** Assignment operator. */
        GridTest& operator= (GridTest<DD> const&) = delete;

        /** Destructor. */
        ~GridTest() = default;

    public:
        /*----------------------------------------------------------------------
        | MEMBER VARIABLES.
        ----------------------------------------------------------------------*/
        // Independent member variables(A).
        /** Start index of grid 'a'. */
        ScaFES::Ntuple<int, 3> mIdxNodeFirstA;

        /** End index of grid 'a'. */
        ScaFES::Ntuple<int, 3> mIdxNodeLastA;

        /** Coordinates of first node number of grid 'a'. */
        ScaFES::Ntuple<double, 3> mCoordNodeFirstA;

        /** Coordinates of last node number of grid 'a'. */
        ScaFES::Ntuple<double, 3> mCoordNodeLastA;

        /*--------------------------------------------------------------------*/
        // Variables depending on above variables(A).
        /** Returns the number of nodes in all dimensions. */
        ScaFES::Ntuple<int, 3> mNnodesA;

        /** Returns the grid size in all dimensions. */
        ScaFES::Ntuple<double, 3> mGridSizeA;

        /** Returns if the grid is numbered in C style or not. */
        bool mIsNumberedStyleCA;

        /*--------------------------------------------------------------------*/
        // Independent member variables(B)
        // Remark: In order to compare the grids a and b, the coordinates
        // of the first and the last grid node have to be same.
        /** Start index of grid 'b'. */
        ScaFES::Ntuple<int, 3> mIdxNodeFirstB;

        /** End index of grid 'b'. */
        ScaFES::Ntuple<int, 3> mIdxNodeLastB;

        /*--------------------------------------------------------------------*/
        /** Index null. */
        ScaFES::Ntuple<int, 3> mIndexNull;

        /** Arbitrary element 'a'. */
        ScaFES::Grid<3> mObjA;

        /** 1:1 copy of element 'a'. */
        ScaFES::Grid<3> mObjCopyA;

        /** Another arbitrary element 'b'. */
        ScaFES::Grid<3> mObjB;

        /** Null element. */
        ScaFES::Grid<3> mObjNull;

        /** Element for lesser and greater tests. */
        ScaFES::Grid<3> mObjLess;

        /** Element for lesser and greater tests. */
        ScaFES::Grid<3> mObjGreater;

        /** Element for lesser and greater tests. */
        ScaFES::Grid<3> mObjUndetermined;
};

TYPED_TEST_CASE_P(GridTest);

/*******************************************************************************
 * LIFE CYCLE METHODS.
 ******************************************************************************/
template<typename DD>
inline GridTest<DD>::GridTest()
    : mIdxNodeFirstA(0, 1, 2)
    , mIdxNodeLastA(9, 16, 25)
    , mCoordNodeFirstA(1.1, 2.2, 3.3)
    , mCoordNodeLastA(2.3, 5.4, 7.5)
    , mNnodesA(10, 16, 24)
    , mGridSizeA(1.2 / 9, 3.2 / 15, 4.2 / 23)
    , mIsNumberedStyleCA(false)
    , mIndexNull(0)
    , mObjA(mIdxNodeFirstA,
            mIdxNodeLastA,
            mCoordNodeFirstA,
            mCoordNodeLastA)
    , mObjCopyA(mIdxNodeFirstA,
                mIdxNodeLastA,
                mCoordNodeFirstA,
                mCoordNodeLastA)
    , mObjB(ScaFES::Ntuple<int, 3>(2, 3, 4),
            ScaFES::Ntuple<int, 3>(4, 5, 8),
            mCoordNodeFirstA,
            mCoordNodeLastA)
    , mObjNull(mIndexNull,
               mIndexNull,
               mCoordNodeFirstA,
               mCoordNodeLastA)
    , mObjLess(ScaFES::Ntuple<int, 3>(1, 2, 3),
               ScaFES::Ntuple<int, 3>(5, 6, 7),
               mCoordNodeFirstA,
               mCoordNodeLastA)
    , mObjGreater(ScaFES::Ntuple<int, 3>(1, 2, 3),
                  ScaFES::Ntuple<int, 3>(8, 9, 10),
                  mCoordNodeFirstA,
                  mCoordNodeLastA)
    , mObjUndetermined(ScaFES::Ntuple<int, 3>(2, 3, 4),
                       ScaFES::Ntuple<int, 3>(5, 6, 7),
                       mCoordNodeFirstA,
                       mCoordNodeLastA)
{ }

} // End of namespace. //
#endif
