#include "gtest/gtest.h"

#include "ScaFES_Communicator.hpp"
#include "ScaFES_Buffer.hpp"

namespace ScaFES_test
{

/*******************************************************************************
 ******************************************************************************/
/**
 * Test class for the class 'Buffer'.
 *
 */
template<typename CT>
class BufferTest : public testing::Test
{
    public:
        /*----------------------------------------------------------------------
        | LIFECYCLE
        ----------------------------------------------------------------------*/
        /** Constructor. */
        BufferTest();

        /** Copy constructor. */
        BufferTest(BufferTest<CT> const&) = delete;

        /** Assignment operator. */
        BufferTest& operator= (BufferTest<CT> const&) = delete;

        /** Destructor. */
        ~BufferTest() = default;

    public:
        /*----------------------------------------------------------------------
        | Member variables.
        ----------------------------------------------------------------------*/
        // Independent member variables(A).
        ScaFES::Communicator mMyWorldA;

        int mIdNeighPartitionA;

        int mSizeBufferA;

        int mUseSkeletonConceptA;

        /*--------------------------------------------------------------------*/
        // Variables depending on above variables(A).

        /*--------------------------------------------------------------------*/
        // Independent member variables(B).
        ScaFES::Communicator mMyWorldB;

        int mIdNeighPartitionB;

        int mSizeBufferB;

        int mUseSkeletonConceptB;

        /*--------------------------------------------------------------------*/
        /** Arbitrary element 'a'. */
        ScaFES::Buffer<CT> mObjA;

        /** 1:1 copy of element 'a'. */
        ScaFES::Buffer<CT> mObjCopyA;

        /** Another arbitrary element 'b'. */
        ScaFES::Buffer<CT> mObjB;
};

TYPED_TEST_CASE_P(BufferTest);

/*******************************************************************************
 ******************************************************************************/
template<typename CT>
inline BufferTest<CT>::BufferTest()
    : mMyWorldA()
    , mIdNeighPartitionA(0)
    , mSizeBufferA(83)
    , mUseSkeletonConceptA(true)
    , mMyWorldB()
    , mIdNeighPartitionB(1)
    , mSizeBufferB(57)
    , mUseSkeletonConceptB(false)
    , mObjA(ScaFES::Buffer<CT>(mMyWorldA,
                               mIdNeighPartitionA,
                               mSizeBufferA,
                               mUseSkeletonConceptA))
    , mObjCopyA(ScaFES::Buffer<CT>(mMyWorldA,
                                   mIdNeighPartitionA,
                                   mSizeBufferA,
                                   mUseSkeletonConceptA))
    , mObjB(ScaFES::Buffer<CT>(mMyWorldB,
                               mIdNeighPartitionB,
                               mSizeBufferB,
                               mUseSkeletonConceptB))
{ }

} // End of namespace. //
