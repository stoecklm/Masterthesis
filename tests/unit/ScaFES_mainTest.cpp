#include "ScaFES_Config.hpp"

#include "gtest/gtest.h"

#ifdef SCAFES_HAVE_BOOST_MPI
#include "boost/mpi.hpp"
#endif

/**
 * Design of test classes:
 * If there exists an ordering, then the order of the objects is as follows:
 * - NULL < A < scalNonZero < B.
 * - Less < Greater
 * - Undetermined cannot be orderd between Less and Greater.
 *
 * \remark Many tests of this class are based on the getter method 'elem()'
 * \code
 * 'EXPECT_TRUE(this->elemV0 == TripX.elem(0));'
 * \endcode
 *
 * Choose characteristic values of equivalence classes (so called 
 * 'representatives').
 * - Positive Tests (with valid inputs)    = Correctness
 * - Negative Tests (with non-valid inputs) = Robustness
 *
 * Test order:
 * 
 * - test getter methods,
 * - test equals and unequals operator (if implemented)
 *    - Is it suffcient to check the equals operator on the basis of the getter
 *      methods of all "independent" member variables of the class
 *      because the dependent methods will be set in an own methods?!
 *    - "independent member variables"
 *        = usually, all input parameters of the class constructor.
 *    - Test "operatorEqualsComponentwise" is not necessary because it
 *      contains the (nearly) same code as the method operator==().
 * - test life cycle methods.
 *
 * - Do not use default values for the test values.
 * - Set different numbers for each element of independent members, if possible.
 * - Compute values of dependent members by hand!
 */

/*******************************************************************************
 ******************************************************************************/
int main(int argc, char** argv)
{
#ifdef SCAFES_HAVE_BOOST_MPI
   boost::mpi::environment myEnv(argc, argv);
#endif
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
