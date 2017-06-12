#include "gtest/gtest.h"

#include <algorithm>
//#include "ScaFES_GroupTest.hpp"
#include "ScaFES_NtupleTest.hpp"
#include "ScaFES_Ntuple.hpp"
#include "ScaFES_Ntuple_Operators.hpp"

namespace ScaFES_test
{

/*******************************************************************************
 ******************************************************************************/
TYPED_TEST_P(NtupleTest, GetElem)
{
    EXPECT_TRUE(this->mElemA0 == this->mObjA.elem(0));
    EXPECT_TRUE(this->mElemA1 == this->mObjA.elem(1));
    EXPECT_TRUE(this->mElemA2 == this->mObjA.elem(2));
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(NtupleTest, ConstructorWithParameters)
{
    TypeParam objTmp(this->mElemA0,this->mElemA1,this->mElemA2);
    EXPECT_TRUE(this->mObjA == objTmp);
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(NtupleTest, LessThanOperatorScalar)
{
    EXPECT_TRUE(this->mObjA < this->mScal);
    EXPECT_FALSE(this->mObjScal < this->mScal);
    EXPECT_FALSE(this->mObjB < this->mScal); // Undetermined.
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(NtupleTest, GreaterOrEqualThanOperatorScalar)
{
    EXPECT_FALSE(this->mObjA >= this->mScal);
    EXPECT_TRUE(this->mObjScal >= this->mScal);
    EXPECT_FALSE(this->mObjB >= this->mScal); // Undetermined.
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(NtupleTest, LessOrEqualThanOperatorScalar)
{
    EXPECT_TRUE(this->mObjA <= this->mScal);
    EXPECT_TRUE(this->mObjScal <= this->mScal);
    EXPECT_FALSE(this->mObjB <= this->mScal); // Undetermined.
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(NtupleTest, GreaterThanOperatorScalar)
{
    EXPECT_FALSE(this->mObjA > this->mScal);
    EXPECT_FALSE(this->mObjScal > this->mScal);
    EXPECT_FALSE(this->mObjB > this->mScal); // Undetermined.
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(NtupleTest, LessThanOperator)
{
    /* A < A => FALSE */
    EXPECT_FALSE(this->mObjA < this->mObjA);
    EXPECT_TRUE(this->mObjLess < this->mObjGreater);
    EXPECT_FALSE(this->mObjGreater < this->mObjLess);
    EXPECT_FALSE(this->mObjUndetermined < this->mObjLess);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(NtupleTest, GreaterOrEqualThanOperator)
{
    /* A >= A => TRUE */
    EXPECT_TRUE(this->mObjA >= this->mObjA);
    EXPECT_FALSE(this->mObjLess >= this->mObjGreater);
    EXPECT_TRUE(this->mObjGreater >= this->mObjLess);
    EXPECT_FALSE(this->mObjUndetermined >= this->mObjLess);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(NtupleTest, LessOrEqualThanOperator)
{
    /* A <= A => TRUE */
    EXPECT_TRUE(this->mObjA <= this->mObjA);
    EXPECT_TRUE(this->mObjLess <= this->mObjGreater);
    EXPECT_FALSE(this->mObjGreater <= this->mObjLess);
    EXPECT_FALSE(this->mObjUndetermined <= this->mObjLess);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(NtupleTest, GreaterThanOperator)
{
    /* A > A => FALSE */
    EXPECT_FALSE(this->mObjA > this->mObjA);
    EXPECT_FALSE(this->mObjLess > this->mObjGreater);
    EXPECT_TRUE(this->mObjGreater > this->mObjLess);
    EXPECT_FALSE(this->mObjUndetermined > this->mObjLess);
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(NtupleTest, addElem)
{
    this->mObjA.addElem(0,this->mScal);
    EXPECT_TRUE(this->mElemA0+this->mScal == this->mObjA.elem(0));
    this->mObjA.addElem(1,this->mScal);
    EXPECT_TRUE(this->mElemA1+this->mScal == this->mObjA.elem(1));
    this->mObjA.addElem(2,this->mScal);
    EXPECT_TRUE(this->mElemA2+this->mScal == this->mObjA.elem(2));
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(NtupleTest, set)
{
    this->mObjA.set(this->mScal,
                    this->mScal1,
                    this->mScal2);
    EXPECT_TRUE(this->mScal  == this->mObjA.elem(0));
    EXPECT_TRUE(this->mScal1 == this->mObjA.elem(1));
    EXPECT_TRUE(this->mScal2 == this->mObjA.elem(2));
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(NtupleTest, idxMaxElem)
{
    EXPECT_TRUE(this->mIdxElemAMax == this->mObjA.idxMaxElem());
    EXPECT_FALSE(this->mIdxElemAMin == this->mObjA.idxMaxElem());
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(NtupleTest, size)
{
    TypeParam sizeA = this->mElemA0 * this->mElemA1 * this->mElemA2;
    TypeParam sizeB = this->mElemB0 * this->mElemB1 * this->mElemB2;
    EXPECT_TRUE(sizeA == this->mObjA.size());
    EXPECT_FALSE(sizeB == this->mObjA.size());
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(NtupleTest, FreeFunctionSize)
{
    TypeParam sizeA = this->mElemA0 * this->mElemA1 * this->mElemA2;
    EXPECT_TRUE(sizeA == ScaFES::size(this->mObjA));
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(NtupleTest, FreeFunctionMin)
{
    this->mObjTmp = ScaFES::min(this->mObjA,this->mObjB);
    EXPECT_TRUE(std::min(this->mElemA0, this->mElemB0) == this->mObjTmp.elem(0));
    EXPECT_TRUE(std::min(this->mElemA1, this->mElemB1) == this->mObjTmp.elem(1));
    EXPECT_TRUE(std::min(this->mElemA2, this->mElemB2) == this->mObjTmp.elem(2));
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(NtupleTest, FreeFunctionMax)
{
    this->mObjTmp = ScaFES::max(this->mObjA,this->mObjB);
    EXPECT_TRUE(std::max(this->mElemA0, this->mElemB0) == this->mObjTmp.elem(0));
    EXPECT_TRUE(std::max(this->mElemA1, this->mElemB1) == this->mObjTmp.elem(1));
    EXPECT_TRUE(std::max(this->mElemA2, this->mElemB2) == this->mObjTmp.elem(2));
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(NtupleTest, FreeFunctionAbsolut)
{
    this->mObjTmp = ScaFES::absolut(this->mObjA);
    EXPECT_TRUE(std::abs(this->mElemA0) == this->mObjTmp.elem(0));
    EXPECT_TRUE(std::abs(this->mElemA1) == this->mObjTmp.elem(1));
    EXPECT_TRUE(std::abs(this->mElemA2) == this->mObjTmp.elem(2));
    /* abs(A) - abs(-A) = 0  forall A */
    this->mObjTmp = ScaFES::absolut(this->mObjA)
                                     - ScaFES::absolut(this->mObjInvAddA);
    EXPECT_TRUE(std::abs(this->mTypeNull) == this->mObjTmp.elem(0));
    EXPECT_TRUE(std::abs(this->mTypeNull) == this->mObjTmp.elem(1));
    EXPECT_TRUE(std::abs(this->mTypeNull) == this->mObjTmp.elem(2));
    /* abs(A) >= 0  forall A */
    EXPECT_TRUE(ScaFES::absolut(this->mObjA) >= this->mObjNeutrAdd);
}
/*----------------------------------------------------------------------------*/
/* Whitebox tests = tests with knowledges about internal functionality */
/* Norm: Check norm properties. */
TYPED_TEST_P(NtupleTest, FreeFunctionNorm1)
{
    /* norm(0)=0 */
    EXPECT_DOUBLE_EQ(ScaFES::norm1(this->mObjNeutrAdd),0.0);
    /* norm(1)=3 */
    EXPECT_DOUBLE_EQ(ScaFES::norm1(this->mObjNeutrMult),3.0);
    /* Positivity: norm(A) >= 0 */
    EXPECT_TRUE(fabs(ScaFES::norm1(this->mObjA) ) >= 0.0);
    /* Scalability: norm(alpha*A) = abs(alpha) * norm(A) */
    EXPECT_DOUBLE_EQ(ScaFES::norm1(this->mScal * this->mObjA),
                             fabs(this->mScal) * ScaFES::norm1(this->mObjA));
    /* Triangle inquality: norm(A+B) <= norm(A) + norm(B). */
    EXPECT_TRUE(  ScaFES::norm1(this->mObjA + this->mObjB)
               <= ScaFES::norm1(this->mObjA) + ScaFES::norm1(this->mObjB) );
    /* Reversed triangle inquality: norm(A) - norm(B) <= norm(A-B). */
    EXPECT_TRUE(  ScaFES::norm1(this->mObjA)- ScaFES::norm1(this->mObjB)
               <= ScaFES::norm1(this->mObjA - this->mObjB) );
    /* norm(A) - norm(-A) = 0 */
    EXPECT_DOUBLE_EQ(ScaFES::norm1(this->mObjA),
                                   ScaFES::norm1(this->mObjInvAddA));
}

/*******************************************************************************
 ******************************************************************************/
REGISTER_TYPED_TEST_CASE_P(NtupleTest,
                           GetElem,
                           ConstructorWithParameters,
                           LessThanOperatorScalar,
                           GreaterOrEqualThanOperatorScalar,
                           LessOrEqualThanOperatorScalar,
                           GreaterThanOperatorScalar,
                           LessThanOperator,
                           GreaterOrEqualThanOperator,
                           LessOrEqualThanOperator,
                           GreaterThanOperator,
                           addElem,
                           set,
                           idxMaxElem,
                           size,
                           FreeFunctionSize,
                           FreeFunctionMin,
                           FreeFunctionMax,
                           FreeFunctionAbsolut,
                           FreeFunctionNorm1
                          );

/*******************************************************************************
 ******************************************************************************/
/*******************************************************************************
 ******************************************************************************/
// typedef ::testing::Types<ScaFES::Ntuple<int>, ScaFES::Ntuple<double>> MyGroupTestTypes;
// INSTANTIATE_TYPED_TEST_CASE_P(My, GroupTest, MyGroupTestTypes);

typedef ::testing::Types<ScaFES::Ntuple<int>, ScaFES::Ntuple<double>> MyNtupleTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(My, NtupleTest, MyNtupleTestTypes);

} // end namespace
