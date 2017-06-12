#include <algorithm>

#include "gtest/gtest.h"

#include "ScaFES_BufferTest.hpp"
#include "ScaFES_Buffer.hpp"
#include "ScaFES_Ntuple.hpp"
#include "ScaFES_Ntuple_FreeFunctions.hpp"
#include "ScaFES_Complex.hpp"
#include "ScaFES_Complex_FreeFunctions.hpp"

namespace ScaFES_test
{

/*******************************************************************************
 * TESTS OF GETTER METHODS.
 ******************************************************************************/
TYPED_TEST_P(BufferTest, GetIdNeighPartition)
{
    EXPECT_TRUE(this->mIdNeighPartitionA == this->mObjA.idNeighPartition());
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(BufferTest, GetSizeBuffer)
{
    EXPECT_TRUE(this->mSizeBufferA == this->mObjA.sizeBuffer());
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(BufferTest, GetUseSkeletonConcept)
{
    EXPECT_TRUE(this->mUseSkeletonConceptA == this->mObjA.useSkeletonConcept());
}

/*******************************************************************************
 * //TEST OF operator==().
 ******************************************************************************/
TYPED_TEST_P(BufferTest, EqualsOperatorComponentwise)
{
    EXPECT_TRUE(this->mObjCopyA.idNeighPartition() == this->mObjA.idNeighPartition());
    EXPECT_TRUE(this->mObjCopyA.sizeBuffer() == this->mObjA.sizeBuffer());
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(BufferTest, EqualsOperator)
{
    EXPECT_TRUE(this->mObjA == this->mObjA);
    EXPECT_TRUE(this->mObjCopyA == this->mObjA);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(BufferTest, UnequalsOperator)
{
    EXPECT_FALSE(this->mObjA != this->mObjA);
    EXPECT_FALSE(this->mObjCopyA != this->mObjA);
}

/*******************************************************************************
 * TESTS OF LIFE CYCLE METHODS.
 ******************************************************************************/
TYPED_TEST_P(BufferTest, ConstructorWithParameters)
{
    ScaFES::Buffer<TypeParam> objTmp(this->mMyWorldA,
                                     this->mIdNeighPartitionA,
                                     this->mSizeBufferA,
                                     this->mUseSkeletonConceptA
                                    );
    EXPECT_TRUE(this->mObjA == objTmp);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(BufferTest, ConstructorDefault)
{
    ScaFES::Buffer<TypeParam> objTmp;
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(BufferTest, CopyConstructor)
{
    ScaFES::Buffer<TypeParam> objTmp = this->mObjA;
    EXPECT_TRUE(this->mObjA == objTmp);
}
/*----------------------------------------------------------------------------*/
TYPED_TEST_P(BufferTest, AssignmentOperator)
{
    ScaFES::Buffer<TypeParam> objTmp;
    objTmp = this->mObjA;
    EXPECT_TRUE(this->mObjA == objTmp);
}

/*******************************************************************************
 * TESTS OF SETTER METHODS.
 ******************************************************************************/
TYPED_TEST_P(BufferTest, SetElemData)
{

    // Attention:
    // setElemData() manipulates the buffer dataSent,
    // elemData() gets the values of the buffer dataReceived.
    const TypeParam EPS = 2.2e-10;
    const TypeParam fact = 3.1415;
    for (int ii = 0; ii < this->mObjA.sizeBuffer(); ++ii) {
        TypeParam tmpObj = static_cast<double>(ii)*fact;
        this->mObjA.setElemData(ii, tmpObj);
        EXPECT_TRUE( ((this->mObjA.dataSent()[ii] - tmpObj) <  EPS)  &&
                     ((this->mObjA.dataSent()[ii] - tmpObj) > -EPS) );
    }
}

/*******************************************************************************
 * TESTS OF WORK METHODS.
 ******************************************************************************/
TYPED_TEST_P(BufferTest, ExchangeValues)
{
    // Special case:
    // Send data from rank 0 to rank 0.
    // ==> The data sent and dat received buffer should contain the same values.
    const TypeParam EPS = 2.2e-10;
    const TypeParam fact = 3.1415;
    bool upper = true;
    bool lower = true;
    ScaFES::Buffer<TypeParam> objTmp(this->mMyWorldA,
                                     0,
                                     this->mSizeBufferA,
                                     this->mUseSkeletonConceptA
                                    );
    for (int ii = 0; ii < objTmp.sizeBuffer(); ++ii) {
        objTmp.setElemData(ii, static_cast<double>(ii)*fact);
    }
    objTmp.sendValues();
    objTmp.receiveValues();
    objTmp.waitAll();
    for (int ii = 0; ii < objTmp.sizeBuffer(); ++ii) {
        upper = ((objTmp.dataSent()[ii] - objTmp.dataReceived()[ii]) <   EPS);
        lower = ((objTmp.dataSent()[ii] - objTmp.dataReceived()[ii]) >  -EPS);
        EXPECT_TRUE(upper && lower);
    }
}
/*******************************************************************************
 * TESTS OF FREE METHODS.
 ******************************************************************************/
TYPED_TEST_P(BufferTest, SwapMethod)
{
    ScaFES::Buffer<TypeParam> tmpA = this->mObjA;
    ScaFES::Buffer<TypeParam> tmpB = this->mObjB;
    ScaFES::swap(this->mObjA, this->mObjB);
    EXPECT_TRUE(tmpA == this->mObjB);
    EXPECT_TRUE(tmpB == this->mObjA);
}

/*******************************************************************************
 ******************************************************************************/
REGISTER_TYPED_TEST_CASE_P(BufferTest,
                           GetIdNeighPartition,
                           GetSizeBuffer,
                           GetUseSkeletonConcept,
                           EqualsOperatorComponentwise,
                           EqualsOperator,
                           UnequalsOperator,
                           ConstructorWithParameters,
                           ConstructorDefault,
                           CopyConstructor,
                           AssignmentOperator,
                           SetElemData,
                           ExchangeValues,
                           SwapMethod
                          );

/*******************************************************************************
 ******************************************************************************/
//typedef ::testing::Types<double> MyBufferTestTypes;
typedef ::testing::Types< float,
                          double,
                          ScaFES::Ntuple<double,1>,
                          ScaFES::Ntuple<double,2>,
                          ScaFES::Ntuple<double,3>,
                          ScaFES::Ntuple<double,4> > MyBufferTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(My, BufferTest, MyBufferTestTypes);

} // End of namespace.
