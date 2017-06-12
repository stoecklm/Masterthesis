#include <algorithm>

#include "gtest/gtest.h"

#ifdef SCAFES_HAVE_ADOLC
#include <adolc/adolc.h>
#include <adolc/adouble.h>
#endif

#include "ScaFES_Ntuple.hpp"
#include "ScaFES_Grid.hpp"
#include "ScaFES_DataField.hpp"

#include "ScaFES_DataFieldTest.hpp"

namespace ScaFES_test
{

const double EPS = 2.2e-12;

/*******************************************************************************
 * TESTS OF GETTER METHODS.
 ******************************************************************************/
TYPED_TEST_P(DataFieldTest, GetName)
{
    EXPECT_TRUE(this->mNameA == this->mObjA.name());
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(DataFieldTest, GetParams)
{
//    EXPECT_TRUE(this->mParamsA == this->mObjA.params());
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(DataFieldTest, GetGridGlobal)
{
    EXPECT_TRUE(this->mGgA == this->mObjA.gridGlobal());
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(DataFieldTest, GetBorderWidth)
{
    EXPECT_TRUE(this->mBorderWidthA == this->mObjA.borderWidth());
}

/*******************************************************************************
 * TESTS OF operator==().
 ******************************************************************************/
TYPED_TEST_P(DataFieldTest, EqualsOperatorComponentwise)
{
    EXPECT_TRUE(this->mObjCopyA.name() == this->mObjA.name());
//    EXPECT_TRUE(this->mObjCopyA.params() == this->mObjA.params());
    EXPECT_TRUE(this->mObjCopyA.gridGlobal() == this->mObjA.gridGlobal());
    EXPECT_TRUE(this->mObjCopyA.borderWidth() == this->mObjA.borderWidth());
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(DataFieldTest, EqualsOperator)
{
    EXPECT_TRUE(this->mObjA == this->mObjA);
    EXPECT_TRUE(this->mObjCopyA == this->mObjA);
}

/*******************************************************************************
 * TESTS OF LIFE CYCLE METHODS.
 ******************************************************************************/
TYPED_TEST_P(DataFieldTest, ConstructorWithParameters)
{
    const int DIM = 3;
    ScaFES::DataField<TypeParam,DIM> objTmp(this->mNameA,
                           this->mParamsA,
                           this->mGgA,
                           this->mBorderWidthA);
    EXPECT_TRUE(this->mObjA == objTmp);
}
/*----------------------------------------------------------------------------*/
// Check the member variables for equality, only!
TYPED_TEST_P(DataFieldTest, ConstructorDefault)
{
    const int DIM = 3;
    ScaFES::Parameters tmpParams;
    ScaFES::GridGlobal<3> tmpGg;
    ScaFES::DataField<TypeParam,DIM> objTmp;
    EXPECT_TRUE(("EMPTYNAME" == objTmp.name()));
// EXPECT_TRUE((tmpParams == objTmp.params()));
    EXPECT_TRUE((tmpGg == objTmp.gridGlobal()));
    EXPECT_TRUE((0 == objTmp.borderWidth()));
}
/*----------------------------------------------------------------------------*/
// TYPED_TEST_P(DataFieldTest, CopyConstructor)
// {
//     const int DIM = 3;
//     ScaFES::DataField<TypeParam,DIM> objTmp = this->mObjA;
//     EXPECT_TRUE(this->mObjA == objTmp);
// }
/*----------------------------------------------------------------------------*/
// TYPED_TEST_P(DataFieldTest, AssignmentOperator)
// {
//     const int DIM = 3;
//     ScaFES::DataField<TypeParam,DIM> objTmp;
//     objTmp = this->mObjA;
//     EXPECT_TRUE(this->mObjA == objTmp);
// }

/*******************************************************************************
 * TESTS OF WORK METHODS.
 ******************************************************************************/
TYPED_TEST_P(DataFieldTest, GetFuncValuesOf)
{
#ifdef SCAFES_HAVE_ADOLC
    const int DIM = 3;
    ScaFES::DataField<double, DIM> objTmp(this->mNameA,
                           funcScalA<adouble>,
                           this->mParamsA,
                           this->mGgA,
                           this->mBorderWidthA);
    ScaFES::DataField<double,DIM> objExpected(this->mNameA,
                           funcScalA<adouble>,
                           this->mParamsA,
                           this->mGgA,
                           this->mBorderWidthA);
    ScaFES::DataField<double,DIM> objActual(this->mNameA,
                           funcScalB<adouble>,
                           this->mParamsA,
                           this->mGgA,
                           this->mBorderWidthA);
    objActual.getFuncValuesOf(objTmp);
    EXPECT_TRUE(objExpected.hasSameValues(objActual, EPS));
#endif
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(DataFieldTest, GetGradOf)
{
#ifdef SCAFES_HAVE_ADOLC
    const int DIM = 3;
    ScaFES::DataField<double, DIM> objTmp(this->mNameA,
                           funcScalA<adouble>,
                           this->mParamsA,
                           this->mGgA,
                           this->mBorderWidthA);
    ScaFES::DataField< ScaFES::Ntuple<double,DIM>, DIM> objExpected(this->mNameA,
                           funcScalGradA<adouble>,
                           this->mParamsA,
                           this->mGgA,
                           this->mBorderWidthA);
    ScaFES::DataField< ScaFES::Ntuple<double,DIM>, DIM> objActual(this->mNameA,
                           funcScalGradB<adouble>,
                           this->mParamsA,
                           this->mGgA,
                           this->mBorderWidthA);
    objActual.getGradOf(objTmp);
    EXPECT_TRUE(objExpected.hasSameValues(objActual, EPS));
#endif
}

/*******************************************************************************
 ******************************************************************************/
REGISTER_TYPED_TEST_CASE_P(DataFieldTest,
                           GetName,
                           GetParams,
                           GetGridGlobal,
                           GetBorderWidth,
                           EqualsOperatorComponentwise,
                           EqualsOperator,
                           ConstructorDefault,
                           ConstructorWithParameters,
                           //CopyConstructor,
//                         AssignmentOperator,
                          GetFuncValuesOf,
                          GetGradOf
                          );

/*******************************************************************************
 ******************************************************************************/
typedef ::testing::Types< double,
                          ScaFES::Complex<double>,
                          ScaFES::Ntuple<double,3> > MyDataFieldTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(My, DataFieldTest, MyDataFieldTestTypes);

} // end namespace
