#ifndef SCAFES_GRIDTEST1D_HPP_
#define SCAFES_GRIDTEST1D_HPP_

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
class GridTest1D : public testing::Test
{
    public:
        /*----------------------------------------------------------------------
        | LIFE CYCLE METHODS.
        ----------------------------------------------------------------------*/
        /** Constructor. */
        GridTest1D();

        /** Copy constructor. */
        GridTest1D(GridTest1D<DD> const&) = delete;

        /** Assignment operator. */
        GridTest1D& operator= (GridTest1D<DD> const&) = delete;

        /** Destructor. */
        ~GridTest1D() = default;

    public:
        /*----------------------------------------------------------------------
        | MEMBER VARIABLES.
        ----------------------------------------------------------------------*/
        // Independent member variables(A).
        /** Start index of grid 'a'. */
        ScaFES::Ntuple<int, 1> mIdxNodeFirstA;

        /** End index of grid 'a'. */
        ScaFES::Ntuple<int, 1> mIdxNodeLastA;

        /** Coordinates of first node number of grid 'a'. */
        ScaFES::Ntuple<double, 1> mCoordNodeFirstA;

        /** Coordinates of last node number of grid 'a'. */
        ScaFES::Ntuple<double, 1> mCoordNodeLastA;

        /*--------------------------------------------------------------------*/
        // Variables depending on above variables(A).
        /** Returns the number of nodes in all dimensions. */
        ScaFES::Ntuple<int, 1> mNnodesA;

        /** Returns the grid size in all dimensions. */
        ScaFES::Ntuple<double, 1> mGridSizeA;

        /** Returns if the grid is numbered in C style or not. */
        bool mIsNumberedStyleCA;

        /*--------------------------------------------------------------------*/
        // Independent member variables(B)
        // Remark: In order to compare the grids a and b, the coordinates
        // of the first and the last grid node have to be same.
        /** Start index of grid 'b'. */
        ScaFES::Ntuple<int, 1> mIdxNodeFirstB;

        /** End index of grid 'b'. */
        ScaFES::Ntuple<int, 1> mIdxNodeLastB;

        /*--------------------------------------------------------------------*/
        /** Index null. */
        ScaFES::Ntuple<int, 1> mIndexNull;

        /** Arbitrary element 'a'. */
        ScaFES::Grid<1> mObjA;

        /** 1:1 copy of element 'a'. */
        ScaFES::Grid<1> mObjCopyA;

        /** Another arbitrary element 'b'. */
        ScaFES::Grid<1> mObjB;

        /** Null element. */
        ScaFES::Grid<1> mObjNull;

        /** Element for lesser and greater tests. */
        ScaFES::Grid<1> mObjLess;

        /** Element for lesser and greater tests. */
        ScaFES::Grid<1> mObjGreater;

        /** Element for lesser and greater tests. */
        ScaFES::Grid<1> mObjUndetermined;
};

TYPED_TEST_CASE_P(GridTest1D);

/*******************************************************************************
 * LIFE CYCLE METHODS.
 ******************************************************************************/
template<typename DD>
inline GridTest1D<DD>::GridTest1D()
    : mIdxNodeFirstA(0)
    , mIdxNodeLastA(9)
    , mCoordNodeFirstA(1.1)
    , mCoordNodeLastA(2.3)
    , mNnodesA(10)
    , mGridSizeA(1.2 / 9)
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
    , mObjB(ScaFES::Ntuple<int, 1>(2),
            ScaFES::Ntuple<int, 1>(4),
            mCoordNodeFirstA,
            mCoordNodeLastA)
    , mObjNull(mIndexNull,
               mIndexNull,
               mCoordNodeFirstA,
               mCoordNodeLastA)
    , mObjLess(ScaFES::Ntuple<int, 1>(1),
               ScaFES::Ntuple<int, 1>(5),
               mCoordNodeFirstA,
               mCoordNodeLastA)
    , mObjGreater(ScaFES::Ntuple<int, 1>(1),
                  ScaFES::Ntuple<int, 1>(8),
                  mCoordNodeFirstA,
                  mCoordNodeLastA)
    , mObjUndetermined(ScaFES::Ntuple<int, 1>(2),
                       ScaFES::Ntuple<int, 1>(5),
                       mCoordNodeFirstA,
                       mCoordNodeLastA)
{ }

} // End of namespace. //
#endif
