#ifndef SCAFES_GRIDGLOBALTEST_HPP_
#define SCAFES_GRIDGLOBALTEST_HPP_

#include "gtest/gtest.h"

#include "ScaFES_Ntuple.hpp"
#include "ScaFES_Grid.hpp"
#include "ScaFES_GridGlobal.hpp"

namespace ScaFES_test
{

/*******************************************************************************
 ******************************************************************************/
/**
 * Test class for the class 'GridGlobal'.
 */
template<typename DD>
class GridGlobalTest : public testing::Test
{
    public:
        /*----------------------------------------------------------------------
        | LIFE CYCLE METHODS.
        ----------------------------------------------------------------------*/
        /** Constructor. */
        GridGlobalTest();

        /** Copy constructor. */
        GridGlobalTest(GridGlobalTest<DD> const&) = delete;

        /** Assignment operator. */
        GridGlobalTest& operator= (GridGlobalTest<DD> const&) = delete;

        /** Destructor. */
        ~GridGlobalTest() = default;

    public:
        /*----------------------------------------------------------------------
        | MEMBER VARIABLES.
        ----------------------------------------------------------------------*/
        // Independent member variables(A).
        /** Start index of grid 'a'. */
        ScaFES::Ntuple<int, 3> mIdxNodeFirstA;

        /** Number of partitions. */
        ScaFES::Ntuple<int, 3> mNpartitionsTotalA;

        /** Number of Nodes. */
        ScaFES::Ntuple<int, 3> mNnodesA;

        /** Coordinates of Node. */
        ScaFES::Ntuple<double, 3> mCoordNodeFirstA;

        /** Coordinates of Node. */
        ScaFES::Ntuple<double, 3> mCoordNodeLastA;

        /** Divide grid in which directions? */
        ScaFES::Ntuple<bool, 3> mDivideGridA;

        /*--------------------------------------------------------------------*/
        // Variables depending on above variables(A).

        /*--------------------------------------------------------------------*/
        // Independent member variables(B)
        /** Number of partitions. */
        ScaFES::Ntuple<int, 3> mNpartitionsTotalB;

        /** Number of nodes. */
        ScaFES::Ntuple<int, 3> mNnodesB;

        /** Coordinates of a node. */
        ScaFES::Ntuple<double, 3> mCoordNodeFirstB;

        /** Coordinates of a node. */
        ScaFES::Ntuple<double, 3> mCoordNodeLastB;

        /** Divide grid in which directions? */
        ScaFES::Ntuple<bool, 3> mDivideGridB;

        /*--------------------------------------------------------------------*/
        /** Arbitrary element 'a'. */
        ScaFES::GridGlobal<3> mObjA;

        /** 1:1 copy of element 'a'. */
        ScaFES::GridGlobal<3> mObjCopyA;

        /** Another arbitrary element 'b'. */
        ScaFES::GridGlobal<3> mObjB;
};

TYPED_TEST_CASE_P(GridGlobalTest);

/*******************************************************************************
 * LIFE CYCLE METHODS.
 ******************************************************************************/
template<typename DD>
inline GridGlobalTest<DD>::GridGlobalTest()
    : mIdxNodeFirstA(0, 1, 2)
    , mNpartitionsTotalA(1,1,1)
    , mNnodesA(5, 10, 15)
    , mCoordNodeFirstA(1.1, 2.2, 3.3)
    , mCoordNodeLastA(2.3, 5.4, 7.5)
    , mDivideGridA(false, false, false)
    , mNpartitionsTotalB(1,1,1)
    , mNnodesB(2, 3, 4)
    , mCoordNodeFirstB(1.5, 2.4, 3.2)
    , mCoordNodeLastB(2.1, 5.0, 6.5)
    , mDivideGridB(true, true, false)
    , mObjA(mNpartitionsTotalA,
            mNnodesA,
            mCoordNodeFirstA,
            mCoordNodeLastA,
            mDivideGridA)
    , mObjCopyA(mNpartitionsTotalA,
                mNnodesA,
                mCoordNodeFirstA,
                mCoordNodeLastA,
                mDivideGridA)
    , mObjB(mNpartitionsTotalB,
            mNnodesB,
            mCoordNodeFirstB,
            mCoordNodeLastB,
            mDivideGridB)
{ }

} // End of namespace. //
#endif

