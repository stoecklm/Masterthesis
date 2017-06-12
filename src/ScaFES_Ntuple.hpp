/* ScaFES
 * Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file ScaFES_Ntuple.hpp
 *  @brief Contains the class template Ntuple.
 */
#ifndef SCAFES_NTUPLE_HPP_
#define SCAFES_NTUPLE_HPP_

#include "ScaFES_Config.hpp"

//#include <cstdint>
#include <iomanip>
#include <iostream>
#include <ios>
#include <vector>
#include <cmath>
#include <type_traits>
#include <stdexcept>

#ifdef SCAFES_HAVE_BOOST
#include <boost/version.hpp>
#endif

#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
namespace boost
{
namespace serialization
{
    class access;
}
}
#include <boost/serialization/version.hpp>
#if BOOST_VERSION < 105900
   #include <boost/serialization/pfto.hpp>
#endif
#endif

namespace ScaFES
{

/*******************************************************************************
 ******************************************************************************/
/** \class Ntuple
 * @brief The class template \c Ntuple represents a
 * \c NN dimensional vector of a given type \c TT.
 *
 * All important unary and binary operators are implemented resp. overloaded
 * and behave as expected. Especially, many operators work componentwise.
 *
 * Usually, the memory for the vector elements must be allocated dynamically,
 * because the number of elements is known at run time. This is related
 * with additional effort for creating and deleting the elements on the
 * heap.
 * Due to the nontype template parameter \c NN the number of vector elements
 * is known at compile time. Thus, it is possible to allocate the memory
 * for the vector elements statically.
 *
 * Advantages:
 * - Automatic compiler SIMD vectorizations can be better applied.
 * - Using the class template \c Ntuple with three elements
 * instead of a class containing explicitly three elements does not affect
 * the performance, but will give more flexibility in solving other than
 * three-dimensional problems.
 *
 * \code
 * // Constructors.
 * ScaFES::Ntuple<double,3> a(1.1, 2.2, 3.3);  // a = (1.1, 2.2, 3.3)
 * ScaFES::Ntuple<double,3> b(6.0, 6.0, 4.0);  // b = (3.0, 3.0, 2.0)
 * ScaFES::Ntuple<double,3> c = a + 2 * b;     // c = (7.1, 8.2, 7.3)
 *
 * // Arithmetic operators.
 * a += 2.0;                                   // a = (3.1, 4.2, 5.3)
 * c = b * b;                                  // c = (9.0, 9.0, 4.0)
 *
 * // Relational operators
 * a < b;                                      // false: 3.3 > 2.0
 * b *= 2.0;                                   // b = (6.0, 6.0, 4.0)
 * b > a;                                      // true
 * \endcode
 */
template <typename TT, std::size_t NN = 3> class Ntuple
{

public:
    /*----------------------------------------------------------------------
    | TYPE DEFINITIONS.
    ----------------------------------------------------------------------*/
    /** Re-export typename TT STL-like. */
    typedef TT value_type;

    /*----------------------------------------------------------------------
    | FRIEND CLASSES.
    ----------------------------------------------------------------------*/
#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
    friend class boost::serialization::access;
#endif

    /*----------------------------------------------------------------------
    | LIFE CYCLE METHODS.
    ----------------------------------------------------------------------*/
    /** Creates the default constructor.
     *  All elements are initialized by zero.*/
    Ntuple<TT, NN>();

    /** Creates a special constructor for twodimensional vectors:
     *  The elements are initialized by given values. */
    Ntuple<TT, NN>(const TT& a, const TT& b);

    /** Creates a special constructor for threedimensional vectors:
     *  The elements are initialized by given values. */
    Ntuple<TT, NN>(const TT& a, const TT& b, const TT& c);

    /** Creates own constructor:
     *  All elements are initalized by the same value. */
    Ntuple<TT, NN>(const TT& scal);

    /** Creates own constructor: Conversion method.
     *  All elements are set by a given STL vector. */
    Ntuple<TT, NN>(const std::vector<TT>& vec);

    /** Creates copy constructor.
     * \remarks Use copy-and-swap idiom.
     */
    Ntuple<TT, NN>(const Ntuple<TT, NN>& rhs);

    /** Creates copy assignment operator.
     * All elements of the rhs. vector will be assigned. */
    Ntuple<TT, NN>& operator=(Ntuple<TT, NN> rhs);

    /** Creates the destructor. */
    ~Ntuple() {};

    /*----------------------------------------------------------------------
    | GETTER METHODS.
    ----------------------------------------------------------------------*/
    /** Returns the dimension (= number of elements) of the n-tuple. */
    size_t dim() const;

    /** Returns a reference to the element at idx \c ii in the n-tuple.
     * \remarks With range check!
     */
    const TT& elem(const std::size_t& ii) const;

    /** Returns the element at idx \c ii in the n-tuple.
     * \remarks Without range check! */
    const TT& operator[](const std::size_t& ii) const;

    /*----------------------------------------------------------------------
    | SETTER METHODS.
    ----------------------------------------------------------------------*/
    /** Returns a reference to the element at idx ii in the Ntuple.
     * \remarks With range check!
     */
    TT& elem(const std::size_t& ii);

    /** Returns a reference to the element at idx \c ii in the n-tuple.
     * \remarks Without range check! */
    TT& operator[](const std::size_t& ii);

    /*----------------------------------------------------------------------
    | ARITHMETIC METHODS.
    ----------------------------------------------------------------------*/
    /** Creates an assignment operator.
     * All elements will be assigned to a given scalar.*/
    ScaFES::Ntuple<TT, NN>& operator=(const TT& sca);

    /** Multiplies this n-tuple by a given scalar (elementwise). */
    template <typename CT> ScaFES::Ntuple<TT, NN>& operator*=(const CT& sca);

    /** Divides this n-tuple by a given scalar (elementwise). */
    template <typename CT> ScaFES::Ntuple<TT, NN>& operator/=(const CT& sca);

    /** Adds a given rhs n-tuple to this n-tuple. */
    ScaFES::Ntuple<TT, NN>& operator+=(const ScaFES::Ntuple<TT, NN>& rhs);

    /** Subtracts a given rhs n-tuple from this n-tuple. */
    ScaFES::Ntuple<TT, NN>& operator-=(const ScaFES::Ntuple<TT, NN>& rhs);

    /** Multiplies this n-tuple by a given rhs n-tuple. */
    ScaFES::Ntuple<TT, NN>& operator*=(const ScaFES::Ntuple<TT, NN>& rhs);

    /** Divides this n-tuple by a given rhs n-tuple. */
    ScaFES::Ntuple<TT, NN>& operator/=(const ScaFES::Ntuple<TT, NN>& rhs);

    /*----------------------------------------------------------------------
    | COMPARISON METHODS.
    ----------------------------------------------------------------------*/
    /** Compares elementwise if this n-tuple is smaller than a a given rhs
     *  n-tuple. */
    bool operator<(const ScaFES::Ntuple<TT, NN>& rhs) const;

    /** Compares elementwise if this n-tuple is greater or equal
     *  than a a given rhs n-tuple. */
    bool operator>=(const ScaFES::Ntuple<TT, NN>& rhs) const;

    /** Compares elementwise if this n-tuple is smaller or equal
     *  than a a given rhs n-tuple. */
    bool operator<=(const ScaFES::Ntuple<TT, NN>& rhs) const;

    /** Compares elementwise if this n-tuple is greater than a a given rhs
     *  n-tuple. */
    bool operator>(const ScaFES::Ntuple<TT, NN>& rhs) const;

    /** Compares elementwise if this n-tuple is equal to a given rhs
     *  n-tuple. */
    bool operator==(const ScaFES::Ntuple<TT, NN>& rhs) const;

    /** Compares elementwise if this n-tuple is unequal to a given rhs
     *  n-tuple. */
    bool operator!=(const ScaFES::Ntuple<TT, NN>& rhs) const;

    /*----------------------------------------------------------------------
    | WORK METHODS.
    ----------------------------------------------------------------------*/
    /** Returns the index of the element with maximal value. */
    std::size_t idxMaxElem() const;

    /** Computes the dimension of the n-tuple,
     * i.e., all elements will be multiplied. */
    TT size() const;

    /** Computes the 1-norm. */
    double norm1() const;

#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
    /** Serializes this class. */
    template <class Archive>
    void serialize(Archive& ar, unsigned int const version);
#endif

    /*----------------------------------------------------------------------
    | FREE METHODS WHICH ARE FRIENDS OF THIS CLASS.
    ----------------------------------------------------------------------*/
    /** Overloads the output operator \c operator<<.
     * An n-tuple is printed this way: '(a_0, a_1,..., a_{n-1})' */
    template <typename ST, std::size_t MM>
    friend std::ostream& operator<<(std::ostream& output,
                                    const ScaFES::Ntuple<ST, MM>& t);

    /** Method to swap members of two ScaFES::Ntuples. */
    template <typename ST, std::size_t MM>
    friend void swap(ScaFES::Ntuple<ST, MM>& first,
                     ScaFES::Ntuple<ST, MM>& second);

    /*----------------------------------------------------------------------
    | CONSTANTS.
    ----------------------------------------------------------------------*/
    /** Dimension of n-tuple. */
    static const std::size_t DIM = NN;

private:
    /*----------------------------------------------------------------------
    | MEMBER VARIABLES.
    ----------------------------------------------------------------------*/
    /** Static array for storing the elements of an n-tuple. */
    TT mElem[NN];

}; // End of class //

/*******************************************************************************
 * FREE METHODS.
 ******************************************************************************/
/** Swaps two objects of the class \c ScaFES::Ntuple. */
template <typename TT, std::size_t NN>
void swap(ScaFES::Ntuple<TT, NN>& first, ScaFES::Ntuple<TT, NN>& second);
/*----------------------------------------------------------------------------*/
/** Preprares writing an object of the class \c ScaFES::Ntuple to output. */
template <typename TT, std::size_t NN>
std::ostream& operator<<(std::ostream& output, const ScaFES::Ntuple<TT, NN>& t);

/*******************************************************************************
 * LIFE CYCLE METHODS.
 ******************************************************************************/
template <typename TT, std::size_t NN> inline Ntuple<TT, NN>::Ntuple()
{
    for (std::size_t iii = 0; iii < NN; ++iii)
    {
        this->mElem[iii] = static_cast<TT>(0);
    }
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
inline Ntuple<TT, NN>::Ntuple(const TT& scal)
{
    for (std::size_t iii = 0; iii < NN; ++iii)
    {
        this->mElem[iii] = scal;
    }
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
inline Ntuple<TT, NN>::Ntuple(const TT& a, const TT& b)
{
    static_assert((NN == 2), "Constructor is valid in case NN=2, only.");
    this->mElem[0] = a;
    this->mElem[1] = b;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
inline Ntuple<TT, NN>::Ntuple(const TT& a, const TT& b, const TT& c)
{
    static_assert((NN == 3), "Constructor is valid in case NN=3, only.");
    this->mElem[0] = a;
    this->mElem[1] = b;
    this->mElem[2] = c;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
inline Ntuple<TT, NN>::Ntuple(const std::vector<TT>& vec)
{
    for (std::size_t iii = 0; iii < NN; ++iii)
    {
        this->mElem[iii] = vec.at(iii);
    }
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
inline Ntuple<TT, NN>::Ntuple(const ScaFES::Ntuple<TT, NN>& rhs)
{
    for (std::size_t iii = 0; iii < rhs.dim(); ++iii)
    {
        this->mElem[iii] = rhs.elem(iii);
    }
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
inline ScaFES::Ntuple<TT, NN>& Ntuple<TT, NN>::
operator=(ScaFES::Ntuple<TT, NN> rhs)
{
    swap(*this, rhs);
    return *this;
}

/*******************************************************************************
 * GETTER METHODS.
 ******************************************************************************/
template <typename TT, std::size_t NN>
inline std::size_t Ntuple<TT, NN>::dim() const
{
    return NN;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
inline const TT& ScaFES::Ntuple<TT, NN>::elem(const std::size_t& idx) const
{
    return this->mElem[idx];
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
inline const TT& ScaFES::Ntuple<TT, NN>::
operator[](const std::size_t& idx) const
{
    return this->mElem[idx];
}

/*******************************************************************************
 * SETTER METHODS.
 ******************************************************************************/
template <typename TT, std::size_t NN>
inline TT& Ntuple<TT, NN>::elem(const std::size_t& idx)
{
    return this->mElem[idx];
    // const_cast<TT&>(static_cast<const
    // ScaFES::Ntuple<TT,NN>&>(*this).elem(idx));
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
inline TT& Ntuple<TT, NN>::operator[](const std::size_t& idx)
{
    return this->mElem[idx];
    // const_cast<TT&>(static_cast<const ScaFES::Ntuple<TT,NN>&>(*this)[idx]);
}

/*******************************************************************************
 * COMPARISON METHODS.
 ******************************************************************************/
template <typename TT, std::size_t NN>
inline bool Ntuple<TT, NN>::operator<(const ScaFES::Ntuple<TT, NN>& rhs) const
{
    bool res = true;

    for (std::size_t iii = 0; iii < NN; ++iii)
    {
        if (!(this->mElem[iii] < rhs.elem(iii)))
        {
            res = false;
            break;
        }
    }

    return res;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
inline bool Ntuple<TT, NN>::operator<=(const ScaFES::Ntuple<TT, NN>& rhs) const
{
    bool res = true;

    for (std::size_t iii = 0; iii < NN; ++iii)
    {
        if (this->mElem[iii] > rhs.elem(iii))
        {
            res = false;
            break;
        }
    }

    return res;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
inline bool Ntuple<TT, NN>::operator>=(const ScaFES::Ntuple<TT, NN>& rhs) const
{
    bool res = true;

    for (std::size_t iii = 0; iii < NN; ++iii)
    {
        if ((this->mElem[iii] < rhs.elem(iii)))
        {
            res = false;
            break;
        }
    }

    return res;
    // return !(*this < rhs); Does not work because we are working elementwise!
}

/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
inline bool Ntuple<TT, NN>::operator>(const ScaFES::Ntuple<TT, NN>& rhs) const
{
    bool res = true;

    for (std::size_t iii = 0; iii < NN; ++iii)
    {
        if (!(this->mElem[iii] > rhs.elem(iii)))
        {
            res = false;
            break;
        }
    }

    return res;
    // return !(*this <= rhs); Does not work because we are working elementwise!
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
inline bool Ntuple<TT, NN>::operator==(const ScaFES::Ntuple<TT, NN>& rhs) const
{
    bool res = true;

    for (std::size_t iii = 0; iii < NN; ++iii)
    {
        if (!(::fabs(this->mElem[iii] - rhs.elem(iii)) < 2.2e-15))
        {
            res = false;
            break;
        }
    }

    return res;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
inline bool Ntuple<TT, NN>::operator!=(const ScaFES::Ntuple<TT, NN>& rhs) const
{
    return !(*this == rhs);
}

/*******************************************************************************
 * WORK METHODS.
 ******************************************************************************/
template <typename TT, std::size_t NN>
inline std::size_t Ntuple<TT, NN>::idxMaxElem() const
{
    std::size_t idx = 0;
    TT maxElem = this->mElem[0];

    for (std::size_t iii = 1; iii < NN; ++iii)
    {
        if (maxElem < this->mElem[iii])
        {
            maxElem = this->mElem[iii];
            idx = iii;
        }
    }

    return idx;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN> inline TT Ntuple<TT, NN>::size() const
{
    TT prod = this->mElem[0];

    for (std::size_t iii = 1; iii < NN; ++iii)
    {
        prod *= this->mElem[iii];
    }

    return prod;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
inline double Ntuple<TT, NN>::norm1() const
{
    double res = ::fabs(this->mElem[0]);

    for (std::size_t iii = 1; iii < NN; ++iii)
    {
        res += ::fabs(this->mElem[iii]);
    }

    return res;
}
/*----------------------------------------------------------------------------*/
#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
template <typename TT, std::size_t NN>
template <class Archive>
inline void Ntuple<TT, NN>::serialize(Archive& ar, unsigned int const version)
{
    if (1 <= version)
    {
        ar&(this->mElem);
    }
}
#endif

/*******************************************************************************
 * ARITHMETIC METHODS.
 ******************************************************************************/
template <typename TT, std::size_t NN>
inline ScaFES::Ntuple<TT, NN>& Ntuple<TT, NN>::operator=(const TT& sca)
{
    for (std::size_t iii = 0; iii < NN; ++iii)
    {
        this->mElem[iii] = static_cast<TT>(sca);
    }

    return *this;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
template <typename CT>
inline ScaFES::Ntuple<TT, NN>& Ntuple<TT, NN>::operator*=(const CT& sca)
{
    for (std::size_t iii = 0; iii < NN; ++iii)
    {
        this->mElem[iii] *= static_cast<TT>(sca);
    }

    return *this;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
template <typename CT>
inline ScaFES::Ntuple<TT, NN>& Ntuple<TT, NN>::operator/=(const CT& sca)
{
    if (::fabs(sca) < 2.2e-12)
    {
        throw std::runtime_error("Given scalar is very near to zero.");
    }
    for (std::size_t iii = 0; iii < NN; ++iii)
    {
        this->mElem[iii] /= static_cast<TT>(sca);
    }

    return *this;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
inline ScaFES::Ntuple<TT, NN>& Ntuple<TT, NN>::
operator+=(const ScaFES::Ntuple<TT, NN>& rhs)
{
    for (std::size_t iii = 0; iii < NN; ++iii)
    {
        this->mElem[iii] += rhs.elem(iii);
    }

    return *this;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
inline ScaFES::Ntuple<TT, NN>& Ntuple<TT, NN>::
operator-=(const ScaFES::Ntuple<TT, NN>& rhs)
{
    for (std::size_t iii = 0; iii < NN; ++iii)
    {
        this->mElem[iii] -= rhs.elem(iii);
    }

    return *this;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
inline ScaFES::Ntuple<TT, NN>& Ntuple<TT, NN>::
operator*=(const ScaFES::Ntuple<TT, NN>& rhs)
{
    for (std::size_t iii = 0; iii < NN; ++iii)
    {
        this->mElem[iii] *= rhs.elem(iii);
    }

    return *this;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
inline ScaFES::Ntuple<TT, NN>& Ntuple<TT, NN>::
operator/=(const ScaFES::Ntuple<TT, NN>& rhs)
{
    for (std::size_t iii = 0; iii < NN; ++iii)
    {
        if (::fabs(rhs.elem(iii)) < 2.2e-12)
        {
            throw std::runtime_error("Element of rhs is very near to zero.");
        }
    }
    for (std::size_t iii = 0; iii < NN; ++iii)
    {
        this->mElem[iii] /= rhs.elem(iii);
    }

    return *this;
}

/*******************************************************************************
* FREE METHODS WHICH ARE FRIENDS OF THIS CLASS.
 ******************************************************************************/
template <typename TT, std::size_t NN>
inline void swap(ScaFES::Ntuple<TT, NN>& first, ScaFES::Ntuple<TT, NN>& second)
{
    TT tmp;

    for (std::size_t iii = 0; iii < second.dim(); ++iii)
    {
        tmp = second.mElem[iii];
        second.mElem[iii] = first.mElem[iii];
        first.mElem[iii] = tmp;
    }
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t NN>
inline std::ostream& operator<<(std::ostream& output,
                                const ScaFES::Ntuple<TT, NN>& t)
{
    output << "[";
    for (std::size_t iii = 0; iii < NN; ++iii)
    {
        output << ::std::setw(4) << ::std::right << t.elem(iii);

        if (NN - 1 > iii)
        {
            output << "; ";
        }
    }

    output << "]";
    return output;
}

} // End of namespace. //

/*******************************************************************************
 ******************************************************************************/
#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
namespace boost
{
namespace serialization
{
    /** Designed to set the boost serialization version of a class template. */
    template <typename TT, std::size_t NN>
    struct version<ScaFES::Ntuple<TT, NN>>
    {
        /** Sets the version number for serialization. */
        BOOST_STATIC_CONSTANT(unsigned long int, value = 2);
    };
}
}
#endif

#endif
