#include "gtest/gtest.h"

#include "ScaFES_Parameters.hpp"

namespace ScaFES_test
{

/*******************************************************************************
 ******************************************************************************/
/** Test class for the class \c DataField. */
template<typename CT>
class ParametersTest : public testing::Test
{
    public:
        /*----------------------------------------------------------------------
        | LIFE CYCLE METHODS.
        ----------------------------------------------------------------------*/
        /** Constructor. */
        ParametersTest();

        /** Copy constructor. */
        ParametersTest(ParametersTest<CT> const&) = delete;

        /** Assignment operator. */
        ParametersTest& operator= (ParametersTest<CT> const&) = delete;

        /** Destructor. */
        ~ParametersTest() = default;

    public:
        /*----------------------------------------------------------------------
        | MEMBER VARIABLES.
        ----------------------------------------------------------------------*/
        // Independent member variables(A).

        /*--------------------------------------------------------------------*/
        // Independent member variables(B)

        /*--------------------------------------------------------------------*/
        /** Arbitrary element 'a'. */
        ScaFES::Parameters mObjA;

        /** 1:1 copy of element 'a'. */
        ScaFES::Parameters mObjCopyA;

        /** Another arbitrary element 'b'. */
        ScaFES::Parameters mObjB;
};

TYPED_TEST_CASE_P(ParametersTest);

/*******************************************************************************
 * LIFE CYCLE METHODS.
 ******************************************************************************/
template<typename CT>
inline ParametersTest<CT>::ParametersTest()
    : mObjA(0, NULL)
    , mObjCopyA(0, NULL)
    , mObjB(0, NULL)
{ }

} // End of namespace. //
