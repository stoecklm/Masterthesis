#include <algorithm>

#include "gtest/gtest.h"

#include "ScaFES_Parameters.hpp"
#include "ScaFES_ParametersTest.hpp"

namespace ScaFES_test
{

/*******************************************************************************
 * TESTS OF GETTER METHODS.
 ******************************************************************************/

/*******************************************************************************
 * TEST OF operator==().
 ******************************************************************************/
TYPED_TEST_P(ParametersTest, EqualsOperator)
{
    EXPECT_TRUE(this->mObjA == this->mObjA);
    EXPECT_TRUE(this->mObjCopyA == this->mObjA);
    EXPECT_TRUE(this->mObjA == this->mObjCopyA);
    EXPECT_TRUE(this->mObjB == this->mObjB);
    EXPECT_TRUE(this->mObjCopyA == this->mObjCopyA);
}

/*******************************************************************************
 * TESTS OF LIFE CYCLE METHODS.
 ******************************************************************************/
TYPED_TEST_P(ParametersTest, CopyConstructor)
{
    ScaFES::Parameters objTmp = this->mObjA;
    EXPECT_TRUE(this->mObjA == objTmp);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(ParametersTest, AssignmentOperator)
{
    ScaFES::Parameters objTmp(0, NULL);
    objTmp = this->mObjA;
    EXPECT_TRUE(this->mObjA == objTmp);
}

/*******************************************************************************
 * TESTS OF WORK METHODS.
 ******************************************************************************/

/*******************************************************************************
 * TESTS OF FREE METHODS.
 ******************************************************************************/
TYPED_TEST_P(ParametersTest, Swap)
{
    ScaFES::Parameters objTmpA = this->mObjA;
    ScaFES::Parameters objTmpB = this->mObjB;
    ScaFES::swap(this->mObjA, this->mObjB);
    EXPECT_TRUE(objTmpA == this->mObjB);
    EXPECT_TRUE(objTmpB == this->mObjA);
}


/*******************************************************************************
 ******************************************************************************/
REGISTER_TYPED_TEST_CASE_P(ParametersTest,
                           EqualsOperator,
                           CopyConstructor,
                           AssignmentOperator,
                           Swap
                          );

/*******************************************************************************
 ******************************************************************************/
typedef ::testing::Types< int > MyParametersTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(My, ParametersTest, MyParametersTestTypes);

} // end namespace
