#ifndef SCAFES_VECTORSPACETEST_HPP_
#define SCAFES_VECTORSPACETEST_HPP_

#include "gtest/gtest.h"

namespace ScaFES_test
{

/*******************************************************************************
 ******************************************************************************/
/**
 * Base test class for testing for group properties.
 *
 * Choose characteristic values of equivalence classes
 * (so called 'representative').
 * Positive Tests (with valid inputs)    = Correctness
 * Negative Tests (with non-valid input) = Robustness 
 */
template <typename S>
class VectorSpaceTest : public testing::Test
{
    public:
        /*----------------------------------------------------------------------
        | LIFECYCLE
        ----------------------------------------------------------------------*/
        /** Constructor. */
        VectorSpaceTest<S>();
        /** Copy constructor. */
        VectorSpaceTest<S>(VectorSpaceTest<S> const&) = delete;
        /** Assignment operator. */
        VectorSpaceTest<S>& operator= (VectorSpaceTest<S> const&) = delete;
        /** Destructor. */
        ~VectorSpaceTest<S>() = default;

    public:
        /*----------------------------------------------------------------------
        | Member variables.
        ----------------------------------------------------------------------*/
        /** Zero element. */
        S mObjZero;
        /** Neutral element w.r.t addition. */
        S mObjNeutrAdd;
        /** Neutral element w.r.t multiplication. */
        S mObjNeutrMult;
        /** Element for storing temporary results. */
        S mObjTmp;

        /** Arbitrary element 'a'. */
        S mObjA;
        /** Additive inverse of element 'a'. */
        S mObjInvAddA;
        /** Multiplicative inverse of element 'a'. */
        S mObjInvMultA;
        /** Element corresponding to a given scalar. */
        S mObjScal;

        /** Another arbitrary element 'b'. */
        S mObjB;

        /** Nonzero scalar. */
        double mScalNonZero;
        /** Zero scalar. */
        double mScalZero;
        /** Arbitrary scalar. */
        double mScal;
        /** Scalar 1 for setter method. */
        double mScal1;
        /** Scalar 2 for setter method. */
        double mScal2;
        /** Zero component. */
        double mTypeNull;
        /** One component. */
        double mTypeOne;

        /** Scalar with neutral element w.r.t addition (if the same for
         *  all components). */
        double mScalNeutrAdd;
        /** Scalar with neutral element w.r.t multiplication (if the same for
         *  all components). */
        double mScalNeutrMult;
        /** Scalar with neutral element w.r.t addition (if the same for
         *  all components). */
        double mScalInvAddA;
        /** Scalar with neutral element w.r.t multiplication (if the same for
         *  all components). */
        double mScalInvMultA;
};

/*******************************************************************************
 ******************************************************************************/
TYPED_TEST_CASE_P(VectorSpaceTest);

/*******************************************************************************
 ******************************************************************************/
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(VectorSpaceTest, ConstructorDefault)
{
    TypeParam objTmp;
    EXPECT_TRUE(this->mObjNeutrAdd == objTmp);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(VectorSpaceTest, CopyConstructor)
{
    TypeParam objTmp = this->mObjA;
    EXPECT_TRUE(this->mObjA == objTmp);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(VectorSpaceTest, AssignmentOperator)
{
    TypeParam objTmp;
    objTmp = this->mObjA;
    EXPECT_TRUE(this->mObjA == objTmp);
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(VectorSpaceTest, MultiplyAndAssignOperatorScalNeutrMult)
{
    /* tmp = A*n_M = A */
    TypeParam objTmp = this->mObjA;
    objTmp *= this->mScalNeutrMult;
    EXPECT_TRUE(this->mObjA == objTmp);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(VectorSpaceTest, MultiplyAndAssignOperatorScalTypeNull)
{
    TypeParam objTmp = this->mObjA;
    objTmp *= this->mTypeNull;
    EXPECT_TRUE(this->mObjZero == objTmp);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(VectorSpaceTest, DivideAndAssignOperatorScalNeutrMult)
{
    TypeParam objTmp = this->mObjA;
    objTmp /= this->mScalNeutrMult;
    EXPECT_TRUE(this->mObjA == objTmp);
}
// /*----------------------------------------------------------------------------*/
// Does not work for int.
// TYPED_TEST_P(VectorSpaceTest, DivideAndAssignOperatorScalInvMult)
// {
//     TypeParam objTmp = this->mObjA;
//     objTmp /= this->mScalInvMultA;
//     EXPECT_TRUE(this->mObjA*this->mObjA == objTmp);
// }
// /*----------------------------------------------------------------------------*/
// TYPED_TEST_P(VectorSpaceTest, DivideAndAssignOperatorScalTypeNull)
// {
//     TypeParam objTmp = this->mObjA;
//     this->mObjA /= this->mTypeNull;
//     EXPECT_TRUE(this->mObjA == objTmp);
// }
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(VectorSpaceTest, AddAndAssignOperatorNeutrAdd)
{
    /* A += neutrAdd. */
    TypeParam objTmp = this->mObjA;
    objTmp += this->mObjNeutrAdd;
    EXPECT_TRUE(this->mObjA == objTmp);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(VectorSpaceTest, AddAndAssignOperatorInvAdd)
{
    /* A += invAdd. */
    TypeParam objTmp = this->mObjA;
    objTmp += this->mObjInvAddA;
    EXPECT_TRUE(this->mObjNeutrAdd == objTmp);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(VectorSpaceTest, AddAndAssignOperatorCommutat)
{
    /* A+B = B+A This should be true if V is a vector space. */
    TypeParam objTmp  = this->mObjA;
    TypeParam objTmp2 = this->mObjB;
    objTmp  += this->mObjB;
    objTmp2 += this->mObjA;
    EXPECT_TRUE(objTmp2 == objTmp);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(VectorSpaceTest, SubtractAndAssignOperatorNeutrAdd)
{
    /* A -= neutrAdd. */
    TypeParam objTmp = this->mObjA;
    objTmp -= this->mObjNeutrAdd;
    EXPECT_TRUE(this->mObjA == objTmp);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(VectorSpaceTest, SubtractAndAssignOperatorInvAdd)
{
    /* A -= invAdd. */
    TypeParam objTmp = this->mObjA;
    objTmp -= this->mObjInvAddA;
    EXPECT_TRUE(2*this->mObjA == objTmp);
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(VectorSpaceTest, MultiplyAndAssignOperatorNeutrMult)
{
    /* A *= neutrMult. */
    TypeParam objTmp = this->mObjA;
    objTmp *= this->mObjNeutrMult;
    EXPECT_TRUE(this->mObjA == objTmp);
}
// /*----------------------------------------------------------------------------*/
// TYPED_TEST_P(VectorSpaceTest, MultiplyAndAssignOperatorInvMult)
// {
//     // Does not work for type 'int' since invMult(A) would have values in (0,1).
//     /* A *= invMultA. */
//     TypeParam objTmp = this->mObjA;
//     objTmp *= this->mObjInvMultA;
//     EXPECT_TRUE(this->mObjNeutrMult == objTmp);
// }
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(VectorSpaceTest, DivideAndAssignOperatorNeutrMult)
{
    /* A /= neutrMult. */
    TypeParam objTmp = this->mObjA;
    objTmp /= this->mObjNeutrMult;
    EXPECT_TRUE(this->mObjA == objTmp);
}
// /*----------------------------------------------------------------------------*/
// TYPED_TEST_P(VectorSpaceTest, DivideAndAssignOperatorInvMultA)
// {
//     // Does not work for type 'int' since invMult(A) would have values in (0,1).
//     /* A /= invMultA. */
//     TypeParam objTmp = this->mObjA;
//     objTmp /= this->mObjInvMultA;
//     EXPECT_TRUE(this->mObjA*this->mObjA == objTmp);
// }
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(VectorSpaceTest, EqualsOperator)
{
    /* A == A => TRUE */
    EXPECT_TRUE(this->mObjA == this->mObjA);
    /* A == B => FALSE*/
    EXPECT_FALSE(this->mObjA == this->mObjB);
    /* B == A => FALSE*/
    EXPECT_FALSE(this->mObjB == this->mObjA);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(VectorSpaceTest, UnequalsOperator)
{
    /* A != A => FALSE */
    EXPECT_FALSE(this->mObjA != this->mObjA);
    /* A != B => TRUE */
    EXPECT_TRUE(this->mObjA != this->mObjB);
    /* B != A => TRUE */
    EXPECT_TRUE(this->mObjB != this->mObjA);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(VectorSpaceTest, FreeFunctionOperatorAddObjObj)
{
    /* A+0 = A. */
    this->mObjTmp = this->mObjA + this->mObjNeutrAdd;
    EXPECT_TRUE(this->mObjA == this->mObjTmp);
    /* 0+A = A. */
    this->mObjTmp = this->mObjA + this->mObjNeutrAdd;
    EXPECT_TRUE(this->mObjA == this->mObjTmp);
    /* A+inv(A) = 0. */
    this->mObjTmp = this->mObjA + this->mObjInvAddA;
    EXPECT_TRUE(this->mObjNeutrAdd == this->mObjTmp);
    /* inv(A)+A = 0. */
    this->mObjTmp = this->mObjInvAddA + this->mObjA;
    EXPECT_TRUE(this->mObjNeutrAdd == this->mObjTmp);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(VectorSpaceTest, FreeFunctionOperatorSubtractObjObj)
{
    /* A-0 = A. */
    this->mObjTmp = this->mObjA - this->mObjNeutrAdd;
    EXPECT_TRUE(this->mObjA == this->mObjTmp);
    /* 0-A = A. */
    this->mObjTmp = this->mObjA - this->mObjNeutrAdd;
    EXPECT_TRUE(this->mObjA == this->mObjTmp);
    /* A-inv(A) = 2*A. */
    this->mObjTmp = this->mObjA - this->mObjInvAddA;
    EXPECT_TRUE(2*this->mObjA == this->mObjTmp);
    /* inv(A)-A = -2*A. */
    this->mObjTmp = this->mObjInvAddA - this->mObjA;
    EXPECT_TRUE(-2*this->mObjA == this->mObjTmp);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(VectorSpaceTest, FreeFunctionOperatorMultiplyObjObj)
{
//     /* A*inv(A) = 1. */
//     this->mObjTmp = this->mObjA * this->mObjInvMultA;
//     EXPECT_TRUE(this->mObjNeutrMult == this->mObjTmp);
//     /* inv(A)*A = 1. */
//     this->mObjTmp = this->mObjInvMultA * this->mObjA;
//     EXPECT_TRUE(this->mObjNeutrMult == this->mObjTmp);
    /* A*1 = A. */
    this->mObjTmp = this->mObjA * this->mObjNeutrMult;
    EXPECT_TRUE(this->mObjA == this->mObjTmp);
    /* 1*A = A. */
    this->mObjTmp = this->mObjA * this->mObjNeutrMult;
    EXPECT_TRUE(this->mObjA == this->mObjTmp);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(VectorSpaceTest, FreeFunctionOperatorMultiplyObjScalar)
{
    /* 0*A = 0. */
    this->mObjTmp = this->mObjA * this->mTypeNull;
    EXPECT_TRUE(this->mObjZero == this->mObjTmp);
    /* 1*A = A. */
    this->mObjTmp = this->mObjA * this->mScalNeutrMult;
    EXPECT_TRUE(this->mObjA == this->mObjTmp);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(VectorSpaceTest, FreeFunctionOperatorMultiplyScalarObj)
{
    /* A*0 = 0. */
    this->mObjTmp = this->mTypeNull * this->mObjA;
    EXPECT_TRUE(this->mObjZero == this->mObjTmp);
    /* A*1 = A. */
    this->mObjTmp = this->mTypeOne * this->mObjA;
    EXPECT_TRUE(this->mObjA == this->mObjTmp);
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(VectorSpaceTest, FreeFunctionOperatorDivideObjObj)
{
    /* A/A = 1. */
    this->mObjTmp = this->mObjA / this->mObjA;
    EXPECT_TRUE(this->mObjNeutrMult == this->mObjTmp);
    /* A/1 = A. */
    this->mObjTmp = this->mObjA / this->mObjNeutrMult;
    EXPECT_TRUE(this->mObjA == this->mObjTmp);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(VectorSpaceTest, FreeFunctionOperatorDivideObjScalarNonZero)
{
    /* A(alpha)/alpha = 1. */
    this->mObjTmp = this->mObjScal / this->mScal;
    EXPECT_TRUE(this->mObjNeutrMult == this->mObjTmp);
    /* A(alpha)/1 = A. */
    this->mObjTmp = this->mObjA / this->mScalNeutrMult;
    EXPECT_TRUE(this->mObjA == this->mObjTmp);
}
// /*----------------------------------------------------------------------------*/
// TYPED_TEST_P(VectorSpaceTest, FreeFunctionOperatorDivideObjScalarZero)
// {
//     this->mObjTmp = this->mObjA / this->mScalZero;
//     EXPECT_TRUE(NaN == this->mObjTmp.elem(0));
//     EXPECT_TRUE(NaN == this->mObjTmp.elem(1));
//     EXPECT_TRUE(NaN == this->mObjTmp.elem(2));
// }

/*******************************************************************************
 ******************************************************************************/
REGISTER_TYPED_TEST_CASE_P(VectorSpaceTest,
                           ConstructorDefault,
                           CopyConstructor,
                           AssignmentOperator,
                           EqualsOperator,
                           UnequalsOperator,
                           MultiplyAndAssignOperatorScalNeutrMult,
                           MultiplyAndAssignOperatorScalTypeNull,
                           DivideAndAssignOperatorScalNeutrMult,
                           AddAndAssignOperatorNeutrAdd,
                           AddAndAssignOperatorInvAdd,
                           AddAndAssignOperatorCommutat,
                           SubtractAndAssignOperatorNeutrAdd,
                           SubtractAndAssignOperatorInvAdd,
                           MultiplyAndAssignOperatorNeutrMult,
                           DivideAndAssignOperatorNeutrMult,
                           FreeFunctionOperatorAddObjObj,
                           FreeFunctionOperatorSubtractObjObj,
                           FreeFunctionOperatorMultiplyObjObj,
                           FreeFunctionOperatorMultiplyObjScalar,
                           FreeFunctionOperatorMultiplyScalarObj,
                           FreeFunctionOperatorDivideObjObj,
                           FreeFunctionOperatorDivideObjScalarNonZero
                          );
} // end namespace
#endif
