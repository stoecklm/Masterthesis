#ifndef SCAFES_NTUPLEEST_HPP_
#define SCAFES_NTUPLEEST_HPP_

#include "gtest/gtest.h"

#include "ScaFES_NTuple.hpp"
#include "ScaFES_NTuple_Operators.hpp"
#include "ScaFES_VectorSpaceTest.hpp"

namespace ScaFES_test
{

/*******************************************************************************
 ******************************************************************************/
/**
 * Test class for the class \c NTuple.
 *
 * If there exists an ordering, then the order of the objects is as follows:
 * NULL < A < scalNonZero < B.
 * Less < Greater
 * Undetermined cannot be orderd between Less and Greater.
 *
 * REMARK: Many tests of this class are based on the getter method 'elem()'
 *         'EXPECT_TRUE(this->elemV0 == TripX.elem(0));'
 * 
 * Choose characteristic values of equivalence classes
 * (so called 'representative').
 * Positive Tests (with valid inputs)    = Correctness
 * Negative Tests (with non-valid input) = Robustness 
 */
template <typename T>
class NTupleTest : public VectorSpaceTest<T>,
                 : public testing::Test
{
    public:
        /*----------------------------------------------------------------------
        | LIFECYCLE
        ----------------------------------------------------------------------*/
        /** Constructor. */
        NTupleTest();
        /** Copy constructor. */
        NTupleTest(NTupleTest<T> const&) = delete;
        /** Assignment operator. */
        NTupleTest& operator= (NTupleTest<T> const&) = delete;
        /** Destructor. */
        virtual ~NTupleTest() = default;
       
    public:
        /*----------------------------------------------------------------------
        | Initialisers.
        ----------------------------------------------------------------------*/
        /** Initialise dependent member variables. */
        int initMemberVar();

    public:
        /*----------------------------------------------------------------------
        | Check for zero.
        ----------------------------------------------------------------------*/
        /** Check for zero. */
        bool checkForZero(double a, double b);

    public:
        /*----------------------------------------------------------------------
        | Member variables.
        ----------------------------------------------------------------------*/
        /** Element for lesser and greater tests. */
        T mObjLess;
        /** Element for lesser and greater tests. */
        T mObjGreater;
        /** Element for lesser and greater tests. */
        T mObjUndetermined;
        /** Element for lesser and greater tests with a scalar. */
        T mObjScal;

        /** Components of element 'a'. */
        double* mElemA;
        /** Dimension of element 'a'. */
        double mDimA;
       
        /** Component of minimal value of element 'a'. */
        int mIdxElemAMin;
        /** Component of maximal value of element 'a'. */
        int mIdxElemAMax;
};

TYPED_TEST_CASE_P(NTupleTest);

/*******************************************************************************
 * INLINED LIFE CYCLE METHODS.
 ******************************************************************************/
template<typename T>
inline NtupleTest<T>::NtupleTest()
   : GroupTest<T>::mScalNeutrAdd(0.0)
   , GroupTest<T>::mScalNeutrMult(1.0)
   , GroupTest<T>::mTypeNull(0.0)
   , GroupTest<T>::mTypeOne(1.0)
   , mElemA0(1.0)
   , mElemA1(2.0)
   , mElemA2(3.0)
   , GroupTest<T>::mScalZero(0.0)
   , GroupTest<T>::mScalNonZero(7.0)
   , mScal1(8.0)
   , Scal2(9.0)
   , mIdxElemAMin(0)
   , mIdxElemAMax(2)
   , mObjLess(T(1.0,2.0,3.0))
   , mObjGreater(T(3.0,6.0,9.0))
   , mObjUndetermined(T(4.0,1.0,2.0))
   , mScal(mScalNonZero)
   , mObjZero(T(mTypeNull, mTypeNull, mTypeNull))
   , mObjTmp(T(mTypeNull, mTypeNull, mTypeNull))
   , mObjNeutrAdd(T(mScalNeutrAdd, mScalNeutrAdd, mScalNeutrAdd))
   , mObjNeutrMult(T(mScalNeutrMult, mElemScalMult, mScalNeutrMult))
   , mObjA(T(mElemA0, mElemA1, mElemA2))
   , mObjInvAddA(T(-mElemA0, -mElemA1, -mElemA2))
   , mObjInvMultA(T(1.0/mElemA0, 1.0/mElemA1, 1.0/mElemA2))
   , mObjB(T(4, 5, 6))
   , mObjScal(T(mScal, mScal, mScal))
{ }

} // end namespace


#endif
