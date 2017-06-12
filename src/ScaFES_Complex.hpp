/* ScaFES
 * Copyright (c) 2011-2016, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file ScaFES_Complex.hpp
 *  @brief Contains the class template Complex.
 */

#ifndef SCAFES_COMPLEX_HPP_
#define SCAFES_COMPLEX_HPP_

#include <iostream>
#include <iomanip>
#include <ios>
#include <cmath>
#include <cstdlib>
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
//#include <boost/serialization/complex.hpp>
#endif

namespace ScaFES
{
/*******************************************************************************
 ******************************************************************************/
/** \class Complex
 * @brief The class template \c Complex represents a complex number
 * of a given type \c TT (usually of type \c double).
 *
 * All important unary and binary operators are implemented resp. overloaded
 * and behave as expected. Especially, many operators work componentwise.
 */
template <typename TT> class Complex
{
public:
    /*----------------------------------------------------------------------
    | TYPE DEFINITIONS.
    ----------------------------------------------------------------------*/
    /** Re-export typename TT STL-like. */
    typedef TT value_type;

#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
    /*----------------------------------------------------------------------
    | FRIEND CLASSES.
    ----------------------------------------------------------------------*/
    friend class boost::serialization::access;
#endif

    /*----------------------------------------------------------------------
    | LIFE CYCLE METHODS.
    ----------------------------------------------------------------------*/
    /** Creates the default constructor. */
    Complex<TT>();

    /** Creates a special constructor for threedimensional vectors:
     *  Both parts are initialized by given values. */
    Complex<TT>(TT const&, TT const&);

    /** Creates own constructor:
     *  \remark Only The REAL PART is initalized by the given value,
     *  the imaginary part remains ZERO! */
    Complex<TT>(TT const& re);

    /** Creates copy constructor. */
    Complex<TT>(ScaFES::Complex<TT> const& rhs);

    /** Creates the destructor. */
    ~Complex();

    /** Creates copy assignment operator using the copy-and-swap idiom. */
    ScaFES::Complex<TT>& operator=(ScaFES::Complex<TT> rhs);

    /*----------------------------------------------------------------------
    | GETTER METHODS.
    ----------------------------------------------------------------------*/
    /** Returns the real part of this complex number. */
    const TT& real() const;

    /** Returns the imaginary part of this complex number. */
    const TT& imag() const;

    /** Returns the element at position \c idx in the complex number. */
    const TT& operator[](std::size_t const& idx) const;

    /*----------------------------------------------------------------------
    | SETTER METHODS.
    ----------------------------------------------------------------------*/
    /** Returns a reference to the element at position \c idx in the complex
     * number.
     */
    TT& operator[](std::size_t const& idx);

    /*----------------------------------------------------------------------
    | COMPARISON METHODS.
    ----------------------------------------------------------------------*/
    /** Compares elementwise if this complex number is smaller
     *  than a given rhs complex number. */
    bool operator<(ScaFES::Complex<TT> const& rhs) const;

    /** Compares elementwise if this complex number is greater or equal
     *  than a a given rhs complex number. */
    bool operator>=(ScaFES::Complex<TT> const& rhs) const;

    /** Compares elementwise if this complex number is smaller or equal
     *  than a a given rhs complex number. */
    bool operator<=(ScaFES::Complex<TT> const& rhs) const;

    /** Compares elementwise if this complex number is greater than
     *  a given rhs complex number. */
    bool operator>(ScaFES::Complex<TT> const& rhs) const;

    /** Compares elementwise if this complex number is equal to a given rhs
     *  complex number. */
    bool operator==(ScaFES::Complex<TT> const& rhs) const;

    /** Compares elementwise if this complex number is unequal to
     * a given rhs complex number. */
    bool operator!=(ScaFES::Complex<TT> const& rhs) const;

    /*----------------------------------------------------------------------
    | WORK METHODS
    ----------------------------------------------------------------------*/
    /** Computes the dimension of the complex number,
     * i.e., all elements will be multiplied.
     */
    TT size() const;

    /** Returns the absolute value. */
    double fabs() const;

#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
    /** Serializes this class. */
    template <class Archive>
    void serialize(Archive& ar, unsigned int const version);
#endif

    /*----------------------------------------------------------------------
    | ARITHMETIC METHODS.
    ----------------------------------------------------------------------*/
    /** Creates an assignment operator.
     * Only the REAL part will be assigned to the given scalar.
     */
    ScaFES::Complex<TT>& operator=(TT const& re);

    /** Computes the additive inverse of this complex number. */
    ScaFES::Complex<TT>& operator-();

    /** Adds this complex number to a given scalar (elementwise). */
    ScaFES::Complex<TT>& operator+=(TT const& sca);

    /** Subtracts this complex number from a given scalar (elementwise). */
    ScaFES::Complex<TT>& operator-=(TT const& sca);

    /** Multiplies this complex number by a given scalar (elementwise). */
    ScaFES::Complex<TT>& operator*=(TT const& sca);

    /** Divides this complex number by a given scalar (elementwise). */
    ScaFES::Complex<TT>& operator/=(TT const& sca);

    /** Adds a given rhs complex number to this complex number. */
    ScaFES::Complex<TT>& operator+=(ScaFES::Complex<TT> const& rhs);

    /** Subtracts a given rhs complex number from this complex number. */
    ScaFES::Complex<TT>& operator-=(ScaFES::Complex<TT> const& rhs);

    /** Multiplies this complex number by a given rhs complex number. */
    ScaFES::Complex<TT>& operator*=(ScaFES::Complex<TT> const& rhs);

    /** Divides this complex number by a given rhs complex number. */
    ScaFES::Complex<TT>& operator/=(ScaFES::Complex<TT> const& rhs);

    /*----------------------------------------------------------------------
    | FREE METHODS WHICH ARE FRIENDS OF THIS CLASS.
    ----------------------------------------------------------------------*/
    /** Prints a complex number in the way: '(re, im})'. */
    template <typename RR>
    friend std::ostream& operator<<(std::ostream& output,
                                    ScaFES::Complex<RR> const& t);

    /** Method to swap two complex numbers. */
    template <typename RR>
    friend void swap(ScaFES::Complex<RR>& first, ScaFES::Complex<RR>& second);

private:
    /*----------------------------------------------------------------------
    | MEMBER VARIABLES.
    ----------------------------------------------------------------------*/
    /** Real part of complex number. */
    TT mReal;

    /** Imaginary part of complex number. */
    TT mImag;
}; // End of class //

/*******************************************************************************
 * FREE METHODS.
 ******************************************************************************/
/** Swaps two objects of the class \c Complex. */
template <typename TT>
void swap(ScaFES::Complex<TT>& first, ScaFES::Complex<TT>& second);
/*----------------------------------------------------------------------------*/
/** Preprares writing a complex number to output. */
template <typename TT>
std::ostream& operator<<(std::ostream& output, ScaFES::Complex<TT> const& rhs);

/*******************************************************************************
 * LIFE CYCLE METHODS.
 ******************************************************************************/
template <typename TT>
inline Complex<TT>::Complex()
: mReal(static_cast<TT>(0)), mImag(static_cast<TT>(0))
{
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline Complex<TT>::Complex(TT const& re, TT const& im)
: mReal(re), mImag(im)
{
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline Complex<TT>::Complex(TT const& re)
: mReal(re), mImag(static_cast<TT>(0))
{
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline ScaFES::Complex<TT>& Complex<TT>::operator=(ScaFES::Complex<TT> rhs)
{
    swap(*this, rhs);
    return *this;
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline Complex<TT>::Complex(ScaFES::Complex<TT> const& rhs)
: mReal(rhs.real()), mImag(rhs.imag())
{
}
/*----------------------------------------------------------------------------*/
template <typename TT> inline Complex<TT>::~Complex()
{
}

/*******************************************************************************
 * GETTER METHODS.
 ******************************************************************************/
template <typename TT> inline const TT& Complex<TT>::real() const
{
    return this->mReal;
}
/*----------------------------------------------------------------------------*/
template <typename TT> inline const TT& Complex<TT>::imag() const
{
    return this->mImag;
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline const TT& Complex<TT>::operator[](std::size_t const& position) const
{
    if (0 == position)
    {
        return this->mReal;
    }
    else
    {
        return this->mImag;
    }
}

/*******************************************************************************
 * SETTER METHODS.
 ******************************************************************************/
template <typename TT>
inline TT& Complex<TT>::operator[](std::size_t const& pos)
{
    return const_cast<TT&>(static_cast<const ScaFES::Complex<TT>&>(*this)[pos]);
}

/*******************************************************************************
 * WORK METHODS.
 ******************************************************************************/
template <typename TT> inline TT Complex<TT>::size() const
{
    return (this->mReal * this->mImag);
}
/*----------------------------------------------------------------------------*/
template <typename TT> inline double Complex<TT>::fabs() const
{
    return (this->mReal * this->mReal + this->mImag * this->mImag);
}
#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
/*----------------------------------------------------------------------------*/
template <typename TT>
template <class Archive>
inline void Complex<TT>::serialize(Archive& ar, unsigned int const version)
{
    if (1 <= version)
    {
        ar&(this->mReal);
        ar&(this->mImag);
    }
}
#endif
/*******************************************************************************
 * ARITHMETIC METHODS.
 ******************************************************************************/
template <typename TT> inline ScaFES::Complex<TT>& Complex<TT>::operator-()
{
    this->mReal *= static_cast<TT>(-1);
    this->mImag *= static_cast<TT>(-1);
    return *this;
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline ScaFES::Complex<TT>& Complex<TT>::operator=(TT const& sca)
{
    this->mReal = sca;
    this->mImag = sca;
    return *this;
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline ScaFES::Complex<TT>& Complex<TT>::operator+=(TT const& sca)
{
    this->mReal += sca;
    this->mImag += sca;
    return *this;
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline ScaFES::Complex<TT>& Complex<TT>::operator-=(TT const& sca)
{
    this->mReal -= sca;
    this->mImag -= sca;
    return *this;
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline ScaFES::Complex<TT>& Complex<TT>::operator*=(TT const& sca)
{
    this->mReal *= sca;
    this->mImag *= sca;
    return *this;
}
/*----------------------------------------------------------------------------*/
// TODO: Throw an exception if the given scalar is near to zero.
template <typename TT>
inline ScaFES::Complex<TT>& Complex<TT>::operator/=(TT const& sca)
{
    this->mReal /= sca;
    this->mImag /= sca;
    return *this;
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline ScaFES::Complex<TT>& Complex<TT>::
operator+=(ScaFES::Complex<TT> const& rhs)
{
    this->mReal += rhs.real();
    this->mImag += rhs.imag();
    return *this;
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline ScaFES::Complex<TT>& Complex<TT>::
operator-=(ScaFES::Complex<TT> const& rhs)
{
    this->mReal -= rhs.real();
    this->mImag -= rhs.imag();
    return *this;
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline ScaFES::Complex<TT>& Complex<TT>::
operator*=(ScaFES::Complex<TT> const& rhs)
{
    const TT tmpRe = this->mReal * rhs.real() - this->mImag * rhs.imag();
    this->mReal = tmpRe;
    const TT tmpIm = this->mReal * rhs.imag() + this->mImag * rhs.real();
    this->mImag = tmpIm;
    return *this;
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline ScaFES::Complex<TT>& Complex<TT>::
operator/=(ScaFES::Complex<TT> const& rhs)
{
    const TT denom = rhs.imag() * rhs.imag() + rhs.real() * rhs.real();
    if (::fabs(denom) < 2.2e-12)
    {
        throw std::runtime_error("Denominator is very near to zero.");
    }
    const TT tmpRe =
        (this->mReal * rhs.real() + this->mImag * rhs.imag()) / denom;
    this->mReal = tmpRe;
    const TT tmpIm =
        (this->mImag * rhs.real() - this->mReal * rhs.imag()) / denom;
    this->mImag = tmpIm;
    return *this;
}

/*******************************************************************************
 * COMPARISON METHODS.
 ******************************************************************************/
template <typename TT>
inline bool Complex<TT>::operator==(ScaFES::Complex<TT> const& rhs) const
{
    bool res = true;

    if (!(::fabs(this->mReal - rhs.real()) < 2.2e-15))
    {
        res = false;
    }

    if (!(::fabs(this->mImag - rhs.imag()) < 2.2e-15))
    {
        res = false;
    }

    return res;
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline bool Complex<TT>::operator!=(ScaFES::Complex<TT> const& rhs) const
{
    return !(*this == rhs);
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline bool Complex<TT>::operator<(ScaFES::Complex<TT> const& rhs) const
{
    bool res = true;

    if (!(this->mReal < rhs.real()))
    {
        res = false;
    }

    if (!(this->mImag < rhs.imag()))
    {
        res = false;
    }

    return res;
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline bool Complex<TT>::operator<=(ScaFES::Complex<TT> const& rhs) const
{
    bool res = true;

    if ((this->mReal > rhs.real()))
    {
        res = false;
    }

    if ((this->mImag > rhs.imag()))
    {
        res = false;
    }

    return res;
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline bool Complex<TT>::operator>=(ScaFES::Complex<TT> const& rhs) const
{
    bool res = true;

    if ((this->mReal < rhs.real()))
    {
        res = false;
    }

    if ((this->mImag < rhs.imag()))
    {
        res = false;
    }

    return res;
    // return !(*this < rhs); Does not work because we are working elementwise!
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline bool Complex<TT>::operator>(ScaFES::Complex<TT> const& rhs) const
{
    bool res = true;

    if (!(mReal > rhs.real()))
    {
        res = false;
    }

    if (!(this->mImag > rhs.imag()))
    {
        res = false;
    }

    return res;
    // return !(*this <= rhs); Does not work because we are working elementwise!
}

/*******************************************************************************
 * FREE METHODS WHICH ARE FRIENDS OF THIS CLASS.
 ******************************************************************************/
template <typename TT>
inline void swap(ScaFES::Complex<TT>& first, ScaFES::Complex<TT>& second)
{
    std::swap(first.mReal, second.mReal);
    std::swap(first.mImag, second.mImag);
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline std::ostream& operator<<(std::ostream& output,
                                ScaFES::Complex<TT> const& rhs)
{
    const TT EPSILON = static_cast<TT>(2.2e-12);

    if (fabs(rhs.real()) < EPSILON)
    {
        output << "[    0,";
    }
    else
    {
        output << "[" << ::std::setw(4) << ::std::right << rhs.real() << ";";
    }

    if (fabs(rhs.imag()) < EPSILON)
    {
        output << "    0i]";
    }
    else
    {
        output << ::std::setw(4) << ::std::right << rhs.imag() << "i]";
    }

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
    template <typename TT> struct version<ScaFES::Complex<TT>>
    {
        /** Sets the version number for serialization. */
        BOOST_STATIC_CONSTANT(unsigned long int, value = 2);
    };
}
}
#endif

#endif
