/* ScaFES
 * Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file ScaFES_Ntuple_FreeFunctions.hpp
 *  @brief Contains free functions related to the class template Ntuple.
 */

/** \defgroup ntupleFree Ntuple free functions
 *
 * The group consists of free functions for manipulating elements of
 * the class template \c Ntuple.
 */

#ifndef SCAFES_NTUPLE_FREEFUNCTIONS_HPP_
#define SCAFES_NTUPLE_FREEFUNCTIONS_HPP_

#include <cmath>

#include "ScaFES_Ntuple.hpp"

namespace ScaFES
{

/*******************************************************************************
 * FREE COMPARISON METHODS.
 ******************************************************************************/
/**
 * Compares elementwise if two given n-tuples are equal.
 * \ingroup ntupleFree
 */
template <typename TT, std::size_t NN>
bool operator==(ScaFES::Ntuple<TT, NN> const& lhs,
                ScaFES::Ntuple<TT, NN> const& rhs);

/**
 * Compares elementwise if two given n-tuples are not equal.
 * \ingroup ntupleFree
 */
template <typename TT, std::size_t NN>
bool operator!=(ScaFES::Ntuple<TT, NN> const& lhs,
                ScaFES::Ntuple<TT, NN> const& rhs);

/*******************************************************************************
 * FREE ARITHMETIC METHODS.
 ******************************************************************************/
/**
 *  Adds two given n-tuples and returns the sum of these n-tuples
 *  as a new object.
 * \ingroup ntupleFree
 */
template <typename TT, std::size_t NN>
ScaFES::Ntuple<TT, NN> operator+(ScaFES::Ntuple<TT, NN> const& lhs,
                                 ScaFES::Ntuple<TT, NN> const& rhs);

/**
 * Subtracts a given n-tuple from another given n-tuple and returns the
 * difference of these n-tuples as a new object.
 * \ingroup ntupleFree
 */
template <typename TT, std::size_t NN>
ScaFES::Ntuple<TT, NN> operator-(ScaFES::Ntuple<TT, NN> const& lhs,
                                 ScaFES::Ntuple<TT, NN> const& rhs);

/**
 * Multiplies two given n-tuples and returns the product of these n-tuples
 *  as a new object.
 * \ingroup ntupleFree
 */
template <typename TT, std::size_t NN>
ScaFES::Ntuple<TT, NN> operator*(ScaFES::Ntuple<TT, NN> const& lhs,
                                 ScaFES::Ntuple<TT, NN> const& rhs);

/**
 * Multiplies a given n-tuple by a given scalar from the right and returns the
 *  product of this arithmetic operation as a new object
 * \ingroup ntupleFree
 */
template <typename TT, std::size_t NN, typename CT>
ScaFES::Ntuple<TT, NN> operator*(ScaFES::Ntuple<TT, NN> const& lhs,
                                 CT const& sca);

/**
 * Multiplies a given n-tuple by a given scalar from the left and returns the
 *  product of this arithmetic operation as a new object.
 * \ingroup ntupleFree
 */
template <typename TT, std::size_t NN, typename CT>
ScaFES::Ntuple<TT, NN> operator*(CT const& sca,
                                 ScaFES::Ntuple<TT, NN> const& rhs);

/**
 *  Divides a given n-tuple by another given n-tuple and returns the
 *  quotient of these n-tuples as a new object.
 * \ingroup ntupleFree
 */
template <typename TT, std::size_t NN>
ScaFES::Ntuple<TT, NN> operator/(ScaFES::Ntuple<TT, NN> const& lhs,
                                 ScaFES::Ntuple<TT, NN> const& rhs);

/**
 * Divides a given n-tuple by a given scalar from the right and returns the
 *  quotient of this arithmetic operation as a new object.
 * \ingroup ntupleFree
 */
template <typename TT, std::size_t NN, typename CT>
ScaFES::Ntuple<TT, NN> operator/(ScaFES::Ntuple<TT, NN> const& lhs,
                                 CT const& sca);

/**
 * Sets elementwise the additive inverse of a given n-tuple and returns the
 *  result as a new object.
 * \ingroup ntupleFree
 */
template <typename TT, std::size_t NN>
inline ScaFES::Ntuple<TT, NN> operator-(ScaFES::Ntuple<TT, NN> const& rhs);

/**
 * Computes the dimension of a given n-tuple, i.e., all elements
 *  will be multiplied.
 * \ingroup ntupleFree
 */
template <typename TT, std::size_t NN> TT size(ScaFES::Ntuple<TT, NN> const& v);

/**
 * Returns the elementwise minimum of two given n-tuples.
 * \ingroup ntupleFree
 */
template <typename TT, std::size_t NN>
ScaFES::Ntuple<TT, NN> min(ScaFES::Ntuple<TT, NN> const& t,
                           ScaFES::Ntuple<TT, NN> const& u);

/**
 * Returns the elementwise maximum of two given n-tuples.
 * \ingroup ntupleFree
 */
template <typename TT, std::size_t NN>
ScaFES::Ntuple<TT, NN> max(ScaFES::Ntuple<TT, NN> const& t,
                           ScaFES::Ntuple<TT, NN> const& u);

/**
 * Returns the absolute values of a given n-tuple.
 * Compares elementwise if two given n-tuples are equal.
 * \ingroup ntupleFree
 */
template <typename TT, std::size_t NN>
ScaFES::Ntuple<TT, NN> absolut(ScaFES::Ntuple<TT, NN> const& t);

/**
 * Computes the 1-norm of a given n-tuple.
 * \ingroup ntupleFree
 */
template <typename TT, std::size_t NN>
double norm1(ScaFES::Ntuple<TT, NN> const& t);

/**
 * Computes the 1-norm of a given n-tuple.
 * \ingroup ntupleFree
 */
template <typename TT, std::size_t NN>
double fabs(ScaFES::Ntuple<TT, NN> const& t);

/*******************************************************************************
 * FREE COMPARISON METHODS.
 ******************************************************************************/
template <typename TT, std::size_t NN>
inline bool operator==(ScaFES::Ntuple<TT, NN> const& lhs,
                       ScaFES::Ntuple<TT, NN> const& rhs)
{
    bool res = true;

    for (std::size_t iii = 0; iii < NN; ++iii)
    {
        if (!(::fabs(lhs.elem(iii) - rhs.elem(iii)) < 2.2e-15))
        {
            res = false;
            break;
        }
    }

    return res;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
inline bool operator!=(ScaFES::Ntuple<TT, NN> const& lhs,
                       ScaFES::Ntuple<TT, NN> const& rhs)
{
    return !(lhs == rhs);
}

/*******************************************************************************
 * FREE ARITHMETIC METHODS.
 ******************************************************************************/
template <typename TT, std::size_t NN>
inline ScaFES::Ntuple<TT, NN> operator+(ScaFES::Ntuple<TT, NN> const& lhs,
                                        ScaFES::Ntuple<TT, NN> const& rhs)
{
    ScaFES::Ntuple<TT, NN> res = lhs;
    res += rhs;
    return res;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
inline ScaFES::Ntuple<TT, NN> operator-(ScaFES::Ntuple<TT, NN> const& lhs,
                                        ScaFES::Ntuple<TT, NN> const& rhs)
{
    ScaFES::Ntuple<TT, NN> res = lhs;
    res -= rhs;
    return res;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
inline ScaFES::Ntuple<TT, NN> operator*(ScaFES::Ntuple<TT, NN> const& lhs,
                                        ScaFES::Ntuple<TT, NN> const& rhs)
{
    ScaFES::Ntuple<TT, NN> res = lhs;
    res *= rhs;
    return res;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN, typename CT>
inline ScaFES::Ntuple<TT, NN> operator*(ScaFES::Ntuple<TT, NN> const& lhs,
                                        CT const& sca)
{
    ScaFES::Ntuple<TT, NN> res = lhs;
    res *= sca;
    return res;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN, typename CT>
inline ScaFES::Ntuple<TT, NN> operator*(CT const& sca,
                                        ScaFES::Ntuple<TT, NN> const& rhs)
{
    ScaFES::Ntuple<TT, NN> res = rhs;
    res *= sca;
    return res;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
inline ScaFES::Ntuple<TT, NN> operator/(ScaFES::Ntuple<TT, NN> const& lhs,
                                        ScaFES::Ntuple<TT, NN> const& rhs)
{
    ScaFES::Ntuple<TT, NN> res = lhs;
    res /= rhs;
    return res;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN, typename CT>
inline ScaFES::Ntuple<TT, NN> operator/(ScaFES::Ntuple<TT, NN> const& lhs,
                                        CT const& sca)
{
    ScaFES::Ntuple<TT, NN> res = lhs;
    res /= sca;
    return res;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
inline ScaFES::Ntuple<TT, NN> operator-(ScaFES::Ntuple<TT, NN> const& rhs)
{
    ScaFES::Ntuple<TT, NN> retTup(NN);

    for (std::size_t iii = 0; iii < NN; ++iii)
    {
        retTup[iii] = -rhs.elem(iii);
    }

    return retTup;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN> TT size(ScaFES::Ntuple<TT, NN> const& v)
{
    return (v.size());
}
/*----------------------------------------------------------------------------*/

template <typename TT, std::size_t NN>
inline ScaFES::Ntuple<TT, NN> min(ScaFES::Ntuple<TT, NN> const& t,
                                  ScaFES::Ntuple<TT, NN> const& u)
{
    ScaFES::Ntuple<TT, NN> retTup(NN);

    for (std::size_t iii = 0; iii < NN; ++iii)
    {
        retTup[iii] = std::min(t.elem(iii), u.elem(iii));
    }

    return retTup;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
inline ScaFES::Ntuple<TT, NN> max(ScaFES::Ntuple<TT, NN> const& t,
                                  ScaFES::Ntuple<TT, NN> const& u)
{
    ScaFES::Ntuple<TT, NN> retTup(NN);

    for (std::size_t iii = 0; iii < NN; ++iii)
    {
        retTup[iii] = std::max(t.elem(iii), u.elem(iii));
    }

    return retTup;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
inline ScaFES::Ntuple<TT, NN> absolut(ScaFES::Ntuple<TT, NN> const& t)
{
    ScaFES::Ntuple<TT, NN> retTup(NN);

    for (std::size_t iii = 0; iii < NN; ++iii)
    {
        retTup[iii] = ::fabs(t.elem(iii));
    }

    return retTup;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
inline double norm1(ScaFES::Ntuple<TT, NN> const& t)
{
    double res = ::fabs(t.elem(0));

    for (std::size_t iii = 1; iii < NN; ++iii)
    {
        res += ::fabs(t.elem(iii));
    }

    return res;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
inline double fabs(ScaFES::Ntuple<TT, NN> const& t)
{
    return norm1(t);
}

} // End of namespace. //

#endif
