#ifndef SCAFES_GRIDTEST2D_HPP_
#define SCAFES_GRIDTEST2D_HPP_

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
class GridTest2D : public testing::Test
{
    public:
        /*----------------------------------------------------------------------
        | LIFE CYCLE METHODS.
        ----------------------------------------------------------------------*/
        /** Constructor. */
        GridTest2D();

        /** Copy constructor. */
        GridTest2D(GridTest2D<DD> const&) = delete;

        /** Assignment operator. */
        GridTest2D& operator= (GridTest2D<DD> const&) = delete;

        /** Destructor. */
        ~GridTest2D() = default;

    public:
        /*----------------------------------------------------------------------
        | MEMBER VARIABLES.
        ----------------------------------------------------------------------*/
        // Independent member variables(A).
        /** Start index of grid 'a'. */
        ScaFES::Ntuple<int, 2> mIdxNodeFirstA;

        /** End index of grid 'a'. */
        ScaFES::Ntuple<int, 2> mIdxNodeLastA;

        /** Coordinates of first node number of grid 'a'. */
        ScaFES::Ntuple<double, 2> mCoordNodeFirstA;

        /** Coordinates of last node number of grid 'a'. */
        ScaFES::Ntuple<double, 2> mCoordNodeLastA;

        /*--------------------------------------------------------------------*/
        // Variables depending on above variables(A).
        /** Returns the number of nodes in all dimensions. */
        ScaFES::Ntuple<int, 2> mNnodesA;

        /** Returns the grid size in all dimensions. */
        ScaFES::Ntuple<double, 2> mGridSizeA;

        /** Returns if the grid is numbered in C style or not. */
        bool mIsNumberedStyleCA;

        /*--------------------------------------------------------------------*/
        // Independent member variables(B)
        // Remark: In order to compare the grids a and b, the coordinates
        // of the first and the last grid node have to be same.
        /** Start index of grid 'b'. */
        ScaFES::Ntuple<int, 2> mIdxNodeFirstB;

        /** End index of grid 'b'. */
        ScaFES::Ntuple<int, 2> mIdxNodeLastB;

        /*--------------------------------------------------------------------*/
        /** Index null. */
        ScaFES::Ntuple<int, 2> mIndexNull;

        /** Arbitrary element 'a'. */
        ScaFES::Grid<2> mObjA;

        /** 1:1 copy of element 'a'. */
        ScaFES::Grid<2> mObjCopyA;

        /** Another arbitrary element 'b'. */
        ScaFES::Grid<2> mObjB;

        /** Null element. */
        ScaFES::Grid<2> mObjNull;

        /** Element for lesser and greater tests. */
        ScaFES::Grid<2> mObjLess;

        /** Element for lesser and greater tests. */
        ScaFES::Grid<2> mObjGreater;

        /** Element for lesser and greater tests. */
        ScaFES::Grid<2> mObjUndetermined;
};

TYPED_TEST_CASE_P(GridTest2D);

/*******************************************************************************
 * LIFE CYCLE METHODS.
 ******************************************************************************/
template<typename DD>
inline GridTest2D<DD>::GridTest2D()
    : mIdxNodeFirstA(0, 1)
    , mIdxNodeLastA(9, 16)
    , mCoordNodeFirstA(1.1, 2.2)
    , mCoordNodeLastA(2.3, 5.4)
    , mNnodesA(10, 16)
    , mGridSizeA(1.2 / 9, 3.2 / 15)
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
    , mObjB(ScaFES::Ntuple<int, 2>(2, 3),
            ScaFES::Ntuple<int, 2>(4, 5),
            mCoordNodeFirstA,
            mCoordNodeLastA)
    , mObjNull(mIndexNull,
               mIndexNull,
               mCoordNodeFirstA,
               mCoordNodeLastA)
    , mObjLess(ScaFES::Ntuple<int, 2>(1, 2),
               ScaFES::Ntuple<int, 2>(5, 6),
               mCoordNodeFirstA,
               mCoordNodeLastA)
    , mObjGreater(ScaFES::Ntuple<int, 2>(1, 2),
                  ScaFES::Ntuple<int, 2>(8, 9),
                  mCoordNodeFirstA,
                  mCoordNodeLastA)
    , mObjUndetermined(ScaFES::Ntuple<int, 2>(2, 3),
                       ScaFES::Ntuple<int, 2>(5, 6),
                       mCoordNodeFirstA,
                       mCoordNodeLastA)
{ }

} // End of namespace. //
#endif
