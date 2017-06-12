#include "gtest/gtest.h"

#include "ScaFES_VectorSpaceTest.hpp"

#include "ScaFES_Ntuple.hpp"
#include "ScaFES_Ntuple_FreeFunctions.hpp"

namespace ScaFES_test
{
  
/*******************************************************************************
 ******************************************************************************/
template<>
VectorSpaceTest< ScaFES::Ntuple<int,1> >::VectorSpaceTest()
 : mObjZero()
 , mObjNeutrAdd()
 , mObjNeutrMult()
 , mObjTmp()
 , mObjA()
 , mObjInvAddA()
 , mObjInvMultA()
 , mObjScal()
 , mObjB()
 , mScalNonZero(7)
 , mScalZero(0)
 , mScal(mScalNonZero)
 , mScal1(8)
 , mScal2(9)
 , mTypeNull(0)
 , mTypeOne(1)
 , mScalNeutrAdd(0)
 , mScalNeutrMult(1)
{
    int* mElemA = new int[1];
    for (int iii = 0; iii < 1; ++iii) {
        mElemA[iii] = iii+1;
    }
    mScalInvAddA  = -mElemA[0];
    mScalInvMultA = 1.0/mElemA[0];
    for (int iii = 0; iii < 1; ++iii) {
        mObjNeutrAdd[iii] = mScalNeutrAdd;
        mObjNeutrMult[iii] = mScalNeutrMult;
        mObjA[iii] = mElemA[iii];
        mObjInvAddA[iii] = -mElemA[iii];
        mObjInvMultA[iii] = 1.0/mElemA[iii];
        mObjB[iii] = mElemA[iii] + 3;
        mObjScal[iii] = mScal;
    }
}
/*******************************************************************************
 ******************************************************************************/
template<>
VectorSpaceTest<ScaFES::Ntuple<int,2>>::VectorSpaceTest()
 : mObjZero()
 , mObjNeutrAdd()
 , mObjNeutrMult()
 , mObjTmp()
 , mObjA()
 , mObjInvAddA()
 , mObjInvMultA()
 , mObjScal()
 , mObjB()
 , mScalNonZero(7)
 , mScalZero(0)
 , mScal(mScalNonZero)
 , mScal1(8)
 , mScal2(9)
 , mTypeNull(0)
 , mTypeOne(1)
 , mScalNeutrAdd(0)
 , mScalNeutrMult(1)
{
    int* mElemA = new int[2];
    for (int iii = 0; iii < 2; ++iii) {
        mElemA[iii] = iii+1;
    }
    mScalInvAddA  = -mElemA[0];
    mScalInvMultA = 1.0/mElemA[0];
    for (int iii = 0; iii < 2; ++iii) {
        mObjNeutrAdd[iii] = mScalNeutrAdd;
        mObjNeutrMult[iii] = mScalNeutrMult;
        mObjA[iii] = mElemA[iii];
        mObjInvAddA[iii] = -mElemA[iii];
        mObjInvMultA[iii] = 1.0/mElemA[iii];
        mObjB[iii] = mElemA[iii] + 3;
        mObjScal[iii] = mScal;
    }
}
/*******************************************************************************
 ******************************************************************************/
template<>
VectorSpaceTest< ScaFES::Ntuple<int,3> >::VectorSpaceTest()
 : mObjZero()
 , mObjNeutrAdd()
 , mObjNeutrMult()
 , mObjTmp()
 , mObjA()
 , mObjInvAddA()
 , mObjInvMultA()
 , mObjScal()
 , mObjB()
 , mScalNonZero(7)
 , mScalZero(0)
 , mScal(mScalNonZero)
 , mScal1(8)
 , mScal2(9)
 , mTypeNull(0)
 , mTypeOne(1)
 , mScalNeutrAdd(0)
 , mScalNeutrMult(1)
{
    int* mElemA = new int[3];
    for (int iii = 0; iii < 3; ++iii) {
        mElemA[iii] = iii+1;
    }
    mScalInvAddA  = -mElemA[0];
    mScalInvMultA = 1.0/mElemA[0];
    for (int iii = 0; iii < 3; ++iii) {
        mObjNeutrAdd[iii] = mScalNeutrAdd;
        mObjNeutrMult[iii] = mScalNeutrMult;
        mObjA[iii] = mElemA[iii];
        mObjInvAddA[iii] = -mElemA[iii];
        mObjInvMultA[iii] = 1.0/mElemA[iii];
        mObjB[iii] = mElemA[iii] + 3;
        mObjScal[iii] = mScal;
    }
}

/*******************************************************************************
 ******************************************************************************/
template<>
VectorSpaceTest< ScaFES::Ntuple<double,1> >::VectorSpaceTest()
 : mObjZero()
 , mObjNeutrAdd()
 , mObjNeutrMult()
 , mObjTmp()
 , mObjA()
 , mObjInvAddA()
 , mObjInvMultA()
 , mObjScal()
 , mObjB()
 , mScalNonZero(7.0)
 , mScalZero(0.0)
 , mScal(mScalNonZero)
 , mScal1(8.0)
 , mScal2(9.0)
 , mTypeNull(0.0)
 , mTypeOne(1.0)
 , mScalNeutrAdd(0.0)
 , mScalNeutrMult(1.0)
{
    double* mElemA = new double[1];
    for (int iii = 0; iii < 1; ++iii) {
        mElemA[iii] = iii+1;
    }
    mScalInvAddA  = -mElemA[0];
    mScalInvMultA = 1.0/mElemA[0];
    for (int iii = 0; iii < 1; ++iii) {
        mObjNeutrAdd[iii] = mScalNeutrAdd;
        mObjNeutrMult[iii] = mScalNeutrMult;
        mObjA[iii] = mElemA[iii];
        mObjInvAddA[iii] = -mElemA[iii];
        mObjInvMultA[iii] = 1.0/mElemA[iii];
        mObjB[iii] = mElemA[iii] + 3.0;
        mObjScal[iii] = mScal;
    }
}
/*******************************************************************************
 ******************************************************************************/
template<>
VectorSpaceTest< ScaFES::Ntuple<double,2> >::VectorSpaceTest()
 : mObjZero()
 , mObjNeutrAdd()
 , mObjNeutrMult()
 , mObjTmp()
 , mObjA()
 , mObjInvAddA()
 , mObjInvMultA()
 , mObjScal()
 , mObjB()
 , mScalNonZero(7.0)
 , mScalZero(0.0)
 , mScal(mScalNonZero)
 , mScal1(8.0)
 , mScal2(9.0)
 , mTypeNull(0.0)
 , mTypeOne(1.0)
 , mScalNeutrAdd(0.0)
 , mScalNeutrMult(1.0)
{
    double* mElemA = new double[2];
    for (int iii = 0; iii < 2; ++iii) {
        mElemA[iii] = iii+1;
    }
    mScalInvAddA  = -mElemA[0];
    mScalInvMultA = 1.0/mElemA[0];
    for (int iii = 0; iii < 2; ++iii) {
        mObjNeutrAdd[iii] = mScalNeutrAdd;
        mObjNeutrMult[iii] = mScalNeutrMult;
        mObjA[iii] = mElemA[iii];
        mObjInvAddA[iii] = -mElemA[iii];
        mObjInvMultA[iii] = 1.0/mElemA[iii];
        mObjB[iii] = mElemA[iii] + 3.0;
        mObjScal[iii] = mScal;
    }
}
/*******************************************************************************
 ******************************************************************************/
template<>
VectorSpaceTest <ScaFES::Ntuple<double,3> >::VectorSpaceTest()
 : mObjZero()
 , mObjNeutrAdd()
 , mObjNeutrMult()
 , mObjTmp()
 , mObjA()
 , mObjInvAddA()
 , mObjInvMultA()
 , mObjScal()
 , mObjB()
 , mScalNonZero(7.0)
 , mScalZero(0.0)
 , mScal(mScalNonZero)
 , mScal1(8.0)
 , mScal2(9.0)
 , mTypeNull(0.0)
 , mTypeOne(1.0)
 , mScalNeutrAdd(0.0)
 , mScalNeutrMult(1.0)
{
    double* mElemA = new double[3];
    for (int iii = 0; iii < 3; ++iii) {
        mElemA[iii] = iii+1;
    }
    mScalInvAddA  = -mElemA[0];
    mScalInvMultA = 1.0/mElemA[0];
    for (int iii = 0; iii < 3; ++iii) {
        mObjNeutrAdd[iii] = mScalNeutrAdd;
        mObjNeutrMult[iii] = mScalNeutrMult;
        mObjA[iii] = mElemA[iii];
        mObjInvAddA[iii] = -mElemA[iii];
        mObjInvMultA[iii] = 1.0/mElemA[iii];
        mObjB[iii] = mElemA[iii] + 3.0;
        mObjScal[iii] = mScal;
    }
}

template<class S>
VectorSpaceTest <S >::VectorSpaceTest() = default;

/*******************************************************************************
 ******************************************************************************/
typedef ::testing::Types< ScaFES::Ntuple<int,1>,
                          ScaFES::Ntuple<int,2>,
                          ScaFES::Ntuple<int,3>,
                          ScaFES::Ntuple<double,1>,
                          ScaFES::Ntuple<double,2>,
                          ScaFES::Ntuple<double,3> > MyVectorSpaceTestTypes;
//typedef ::testing::Types< ScaFES::Triple<int>, ScaFES::Triple<double> > MyVectorSpaceTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(MyScaFES_Ntuple, VectorSpaceTest, MyVectorSpaceTestTypes);
}

/*******************************************************************************
 ******************************************************************************/
int main(int argc, char** argv)
{
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}

