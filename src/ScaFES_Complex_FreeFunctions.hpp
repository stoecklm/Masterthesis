/* ScaFES
 * Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file ScaFES_Complex_FreeFunctions.hpp
 *  @brief Contains free functions related to the class template Complex.
 */

/** \defgroup complexFree Complex free functions
 *
 * The group consists of free functions for manipulating elements of
 * the class template \c Complex.
 */

#ifndef SCAFES_COMPLEX_FREEFUNCTIONS_HPP_
#define SCAFES_COMPLEX_FREEFUNCTIONS_HPP_

#include "ScaFES_Complex.hpp"

namespace ScaFES
{

/*******************************************************************************
 * FREE COMPARISON METHODS.
 ******************************************************************************/
/** Compares elementwise if two given complex numbers are equal.
 * \ingroup complexFree
 */
template <typename TT>
bool operator==(ScaFES::Complex<TT> const& lhs, ScaFES::Complex<TT> const& rhs);

/** Compares elementwise if two given complex numbers are not equal.
 * \ingroup complexFree
 */
template <typename TT>
bool operator!=(ScaFES::Complex<TT> const& lhs, ScaFES::Complex<TT> const& rhs);

/*******************************************************************************
 * FREE ARITHMETIC METHODS.
 ******************************************************************************/
/** Adds two given complex numbers and returns the sum of these complex numbers
 *  as a new object. (Makes uses of copy constructor and += operator).
 * \ingroup complexFree
 */
template <typename TT>
inline ScaFES::Complex<TT> operator+(ScaFES::Complex<TT> const& lhs,
                                     ScaFES::Complex<TT> const& rhs);

/** Adds a given complex number to a given real number and
    returns the sum as a new object.
 * \ingroup complexFree
 */
template <typename TT>
ScaFES::Complex<TT> operator+(ScaFES::Complex<TT> const& lhs, TT const& sca);

/** Subtracts a given complex number from another given complex number and
    returns the difference of these complex numbers as a new object.
 * \ingroup complexFree
 */
template <typename TT>
ScaFES::Complex<TT> operator-(ScaFES::Complex<TT> const& lhs,
                              ScaFES::Complex<TT> const& rhs);

/** Subtracts a given complex number from a given real number and
    returns the difference as a new object.
 * \ingroup complexFree
 */
template <typename TT>
ScaFES::Complex<TT> operator-(ScaFES::Complex<TT> const& lhs, TT const& sca);

/** Multiplies two given complex numbers and returns the product of these
 *  complex numbers as a new object.
 * \ingroup complexFree
 */
template <typename TT>
ScaFES::Complex<TT> operator*(ScaFES::Complex<TT> const& lhs,
                              ScaFES::Complex<TT> const& rhs);

/** Multiplies a given complex number by a given scalar from the right
 *  and returns the product of this arithmetic operation as a new object.
 * \ingroup complexFree
 */
template <typename TT>
ScaFES::Complex<TT> operator*(ScaFES::Complex<TT> const& lhs, TT const& sca);

/** Multiplies a given complex number by a given scalar from the left and
 * returns the product of this arithmetic operation as a new object.
 * \ingroup complexFree
 */
template <typename TT>
ScaFES::Complex<TT> operator*(TT const& sca, ScaFES::Complex<TT> const& rhs);

/** Divides a given complex number by another given complex number and
 *  returns the quotient of these complex numbers as a new object.
 * \ingroup complexFree
 */
template <typename TT>
ScaFES::Complex<TT> operator/(ScaFES::Complex<TT> const& lhs,
                              ScaFES::Complex<TT> const& rhs);

/** Divides a given complex number by a given scalar from the right
 *  and returns the quotient of this arithmetic operation as a new object.
 * \ingroup complexFree
 */
template <typename TT>
ScaFES::Complex<TT> operator/(ScaFES::Complex<TT> const& lhs, TT const& sca);

/** Sets elementwise the additive inverse a given complex number and returns the
 *  result as a new object.
 * \ingroup complexFree
 */
template <typename TT>
ScaFES::Complex<TT> operator-(ScaFES::Complex<TT> const& rhs);

/** Computes the dimension of a given complex number, i.e., all elements
 *  will be multiplied.
 * \ingroup complexFree
 */
template <typename TT> TT size(ScaFES::Complex<TT> const& tt);

/** Returns the absolute values of a given complex number.
 * \ingroup complexFree
 */
template <typename TT>
const ScaFES::Complex<TT> conj(ScaFES::Complex<TT> const& tt);

/** Returns the absolute values of a given complex number.
 * \ingroup complexFree
 */
template <typename TT> TT absolute(ScaFES::Complex<TT> const& tt);

/** Returns the absolute values of a given complex number.
 * \ingroup complexFree
 */
template <typename TT> TT fabs(ScaFES::Complex<TT> const& tt);

/** Computes the 1-norm of a given complex number.
 * \ingroup complexFree
 */
template <typename TT> inline TT norm1(ScaFES::Complex<TT> const& tt);

/*******************************************************************************
 * FREE COMPARISON METHODS.
 ******************************************************************************/
template <typename TT>
inline bool operator==(ScaFES::Complex<TT> const& lhs,
                       ScaFES::Complex<TT> const& rhs)
{
    bool res = true;
    if (!(fabs(lhs.real() - rhs.real()) < 2.2e-15))
    {
        res = false;
    }
    if (!(fabs(lhs.imag() - rhs.imag()) < 2.2e-15))
    {
        res = false;
    }
    return res;
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline bool operator!=(ScaFES::Complex<TT> const& lhs,
                       ScaFES::Complex<TT> const& rhs)
{
    return !(lhs == rhs);
}

/*******************************************************************************
 * FREE ARITHMETIC METHODS.
 ******************************************************************************/
template <typename TT>
inline ScaFES::Complex<TT> operator+(ScaFES::Complex<TT> const& lhs,
                                     ScaFES::Complex<TT> const& rhs)
{
    ScaFES::Complex<TT> res = lhs;
    res += rhs;
    return res;
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline ScaFES::Complex<TT> operator+(ScaFES::Complex<TT> const& lhs,
                                     TT const& sca)
{
    ScaFES::Complex<TT> res = lhs;
    res += sca;
    return res;
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline ScaFES::Complex<TT> operator-(ScaFES::Complex<TT> const& lhs,
                                     ScaFES::Complex<TT> const& rhs)
{
    ScaFES::Complex<TT> res = lhs;
    res -= rhs;
    return res;
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline ScaFES::Complex<TT> operator-(ScaFES::Complex<TT> const& lhs,
                                     TT const& sca)
{
    ScaFES::Complex<TT> res = lhs;
    res -= sca;
    return res;
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline ScaFES::Complex<TT> operator*(ScaFES::Complex<TT> const& lhs,
                                     ScaFES::Complex<TT> const& rhs)
{
    ScaFES::Complex<TT> res = lhs;
    res *= rhs;
    return res;
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline ScaFES::Complex<TT> operator*(ScaFES::Complex<TT> const& lhs,
                                     TT const& sca)
{
    ScaFES::Complex<TT> res = lhs;
    res *= sca;
    return res;
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline ScaFES::Complex<TT> operator*(TT const& sca,
                                     ScaFES::Complex<TT> const& rhs)
{
    ScaFES::Complex<TT> res = rhs;
    res *= sca;
    return res;
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline ScaFES::Complex<TT> operator/(ScaFES::Complex<TT> const& lhs,
                                     ScaFES::Complex<TT> const& rhs)
{
    ScaFES::Complex<TT> res = lhs;
    res /= rhs;
    return res;
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline ScaFES::Complex<TT> operator/(ScaFES::Complex<TT> const& lhs,
                                     TT const& sca)
{
    ScaFES::Complex<TT> res = lhs;
    res /= sca;
    return res;
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline ScaFES::Complex<TT> operator-(ScaFES::Complex<TT> const& rhs)
{
    return ScaFES::Complex<TT>(-rhs.real(), -rhs.imag());
}
/*----------------------------------------------------------------------------*/
template <typename TT> inline TT size(ScaFES::Complex<TT> const& tt)
{
    return (tt.size());
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline const ScaFES::Complex<TT> conj(ScaFES::Complex<TT> const& tt)
{
    return ScaFES::Complex<TT>(tt.real(), -tt.imag());
}
/*----------------------------------------------------------------------------*/
template <typename TT> inline TT absolute(ScaFES::Complex<TT> const& tt)
{
    return (tt.real() * tt.real() + tt.imag() * tt.imag());
}
/*----------------------------------------------------------------------------*/
template <typename TT> inline TT fabs(ScaFES::Complex<TT> const& tt)
{
    return absolute<TT>(tt);
}
/*----------------------------------------------------------------------------*/
template <typename TT> inline TT norm1(ScaFES::Complex<TT> const& tt)
{
    return (fabs(tt.real()) + fabs(tt.imag()));
}

} // End of namespace. //

#endif
