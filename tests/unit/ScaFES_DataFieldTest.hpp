#include "gtest/gtest.h"

#include "ScaFES_Ntuple.hpp"
#include "ScaFES_Parameters.hpp"
#include "ScaFES_GridGlobal.hpp"
#include "ScaFES_DataField.hpp"

namespace ScaFES_test
{

/*******************************************************************************
 ******************************************************************************/
/** Test class for the class \c DataField. */
template<typename CT>
class DataFieldTest : public testing::Test
{
    public:
        /*----------------------------------------------------------------------
        | LIFE CYCLE METHODS.
        ----------------------------------------------------------------------*/
        /** Constructor. */
        DataFieldTest();

        /** Copy constructor. */
        DataFieldTest(DataFieldTest<CT> const&) = delete;

        /** Assignment operator. */
        DataFieldTest& operator=(DataFieldTest<CT> const&) = delete;

        /** Destructor. */
        ~DataFieldTest() = default;

    public:
        /*----------------------------------------------------------------------
        | MEMBER VARIABLES.
        ----------------------------------------------------------------------*/
        // Independent member variables(A).
        std::string mNameA;
        ScaFES::Parameters mParamsA;
        ScaFES::GridGlobal<3> mGgA;
        int mBorderWidthA;

        /*--------------------------------------------------------------------*/
        // Independent member variables(B)
        std::string mNameB;
        ScaFES::Parameters mParamsB;
        ScaFES::GridGlobal<3> mGgB;
        int mBorderWidthB;

        /*--------------------------------------------------------------------*/
        /** Arbitrary element 'a'. */
        ScaFES::DataField<CT, 3> mObjA;

        /** 1:1 copy of element 'a'. */
        ScaFES::DataField<CT, 3> mObjCopyA;

        /** Another arbitrary element 'b'. */
        ScaFES::DataField<CT, 3> mObjB;
};

TYPED_TEST_CASE_P(DataFieldTest);

/*******************************************************************************
 * FREE METHODS.
 ******************************************************************************/
template<typename TT>
void funcScalA(TT& fx, ScaFES::Ntuple<TT, 3> const& x, TT const& t)
{
    fx = 2.0 * x[0] * x[0] + 3.0 * x[1] * x[1] + 7.0 * x[2];
}
/*----------------------------------------------------------------------------*/
template<typename TT>
inline void funcScalGradA(ScaFES::Ntuple<TT, 3>& fx,
                          ScaFES::Ntuple<TT, 3> const& x, TT const& t)
{
    fx[0] = 4.0 * x[0];
    fx[1] = 6.0 * x[1];
    fx[2] = 7.0;
}
/*----------------------------------------------------------------------------*/
template<typename TT>
inline void funcScalB(TT& fx, ScaFES::Ntuple<TT, 3> const& x, TT const& t)
{
    fx = 3.0 * x[0] * x[0] - 1.6 * x[1] * x[1] - 5.0 * x[2] * x[2] * x[2];
}
/*----------------------------------------------------------------------------*/
template<typename TT>
inline void funcScalGradB(ScaFES::Ntuple<TT, 3>& fx,
                          ScaFES::Ntuple<TT, 3> const& x, TT const& t)
{
    fx[0] =   6.0 * x[0];
    fx[1] =  -3.2 * x[1];
    fx[2] = -15.0 * x[2] * x[2];
}
/*----------------------------------------------------------------------------*/
template<typename TT>
inline void funcVectC(ScaFES::Ntuple<TT, 4>& fx,
                      ScaFES::Ntuple<TT, 3> const& x, TT const& t)
{
    fx[0] = 2.0 * x[0] * x[0] + 3.0 * x[1] * x[1] + 7.0 * x[2];
    fx[1] = 2.0 * x[1] * x[1] + 3.0 * x[2] * x[2] + 7.0 * x[0];
    fx[2] = 2.0 * x[2] * x[2] + 3.0 * x[0] * x[0] + 7.0 * x[1];
    fx[3] = 1.0 * x[2] + 2.0 * x[0] + 3.0 * x[1];
}
/*----------------------------------------------------------------------------*/
template<typename TT>
inline void funcVectGradC(ScaFES::Ntuple< ScaFES::Ntuple<TT, 3>, 4>& fx,
                          ScaFES::Ntuple<TT, 3> const& x, TT const& t)
{
    fx[0][0] = 4.0 * x[0];
    fx[0][1] = 6.0 * x[1];
    fx[0][2] = 7.0;
    fx[1][0] = 7.0;
    fx[1][1] = 4.0 * x[1];
    fx[1][2] = 6.0 * x[2];
    fx[2][0] = 6.0 * x[0];
    fx[2][1] = 7.0;
    fx[2][2] = 4.0 * x[2];
    fx[3][0] = 2.0;
    fx[3][1] = 3.0;
    fx[3][2] = 1.0;
}

/*******************************************************************************
 * LIFE CYCLE METHODS.
 ******************************************************************************/
template<typename CT>
inline DataFieldTest<CT>::DataFieldTest()
    : mNameA("mObjA")
    , mParamsA()
    , mGgA(mParamsA)
    , mBorderWidthA(1)
    , mNameB("mObjB")
    , mParamsB()
    , mGgB(mParamsB)
    , mBorderWidthB(2)
    , mObjA(mNameA, mParamsA, mGgA, mBorderWidthA)
    , mObjCopyA(mNameA, mParamsA, mGgA, mBorderWidthA)
    , mObjB(mNameB, mParamsB, mGgB, mBorderWidthB)
{ }

} // End of namespace. //
