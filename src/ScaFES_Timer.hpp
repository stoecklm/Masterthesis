/* ScaFES
 * Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file ScaFES_Timer.hpp
 *  @brief Contains the class Timer.
 */

#ifndef SCAFES_TIMER_HPP_
#define SCAFES_TIMER_HPP_

#include "ScaFES_Config.hpp"

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

extern "C" {
#include <sys/time.h>
}

#ifdef _OPENMP
extern "C" {
#include <omp.h>
}
#endif

#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
#include <boost/serialization/serialization.hpp>
#endif

// #ifdef SCAFES_HAVE_BOOST_SERIALIZATION
// namespace boost
// {
// namespace serialization
// {
// class access;
// }
// }
// #include <boost/serialization/vector.hpp>
// #include <boost/serialization/version.hpp>
// #endif

#ifdef SCAFES_HAVE_BOOST_MPI
#include <boost/mpi.hpp>
/*#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/collectives.hpp>*/
#endif

#ifdef SCAFES_HAVE_BOOST_MPI
namespace boost_mpi = boost::mpi;
#endif

namespace ScaFES
{

class Timer
{
public:
    /*----------------------------------------------------------------------
    | FRIEND CLASSES.
    ----------------------------------------------------------------------*/
#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
    friend class boost::serialization::access;
#endif

    /*----------------------------------------------------------------------
    | LIFE CYCLE METHODS.
    ----------------------------------------------------------------------*/
    /** Creates own default constructor. */
    Timer();

    /** Creates own copy constructor. */
    Timer(Timer const&);

    /** Creates default destructor. */
    ~Timer() = default;

    /** Creates own assignment operator. */
    Timer& operator=(Timer rhs);

    /*----------------------------------------------------------------------
    | WORK METHODS.
     ----------------------------------------------------------------------*/
#ifdef SCAFES_HAVE_BOOST_MPI
    /** Returns Boost MPI timer. */
    boost_mpi::timer const& timer() const;
#else
    /** Returns timer. */
    int const& timer() const;
#endif

    /** Restarts the timer. */
    void restart();

    /** Returns the elapsed time. */
    double elapsed();

    /*----------------------------------------------------------------------
    | FREE METHODS.
    ----------------------------------------------------------------------*/
    /** Swaps the members of two communicators. */
    friend void swap(ScaFES::Timer& first,
                     ScaFES::Timer& second);

private:
    /*----------------------------------------------------------------------
    | WORK METHODS.
    ----------------------------------------------------------------------*/
#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
    /** Serializes class. */
    template <class Archive>
    void serialize(Archive& ar, const unsigned int version);
#endif

    /*----------------------------------------------------------------------
    | MEMBER VARIABLES.
    ----------------------------------------------------------------------*/
    /** Elapsed wall clock time. */
    double mWtime;

#ifdef SCAFES_HAVE_BOOST_MPI
    /** Boost MPI timer. */
    boost_mpi::timer mTimer;
#else
    int mTimer;
#endif
};

/*******************************************************************************
 * FREE METHODS.
 ******************************************************************************/
/** Swaps two objects of the class \c Timer. */
void swap(ScaFES::Timer& first, ScaFES::Timer& second);

/*******************************************************************************
 * LIFE CYCLE METHODS.
 ******************************************************************************/
inline Timer::Timer()
: mWtime(0.0)
, mTimer()
{ }
/*----------------------------------------------------------------------------*/
inline Timer::Timer(Timer const& rhs)
: mWtime(rhs.mWtime)
, mTimer(rhs.timer())
{
}
/*----------------------------------------------------------------------------*/
inline Timer& Timer::operator=(Timer rhs)
{
    swap(*this, rhs);
    return *this;
}
/*----------------------------------------------------------------------------*/
#ifdef SCAFES_HAVE_BOOST_MPI
inline boost_mpi::timer const& Timer::timer() const
#else
inline int const& Timer::timer() const
#endif
{
    return this->mTimer;
}
/*----------------------------------------------------------------------------*/
inline void Timer::restart()
{
#ifdef SCAFES_HAVE_MPI
#ifdef SCAFES_HAVE_BOOST_MPI
    this->mTimer.restart();
    this->mWtime = 0.0;
#else
    this->mWtime = MPI_Wtime();
#endif
#else
#ifdef _OPENMP
    this->mWtime = omp_get_wtime();
#else
    struct timeval ct;
    gettimeofday(&ct, NULL);
    this->mWtime = static_cast<double>(ct.tv_sec)
                   + (static_cast<double>(ct.tv_usec) / 1000000.0);
#endif
#endif
}
/*----------------------------------------------------------------------------*/
inline double Timer::elapsed()
{
#ifdef SCAFES_HAVE_MPI
#ifdef SCAFES_HAVE_BOOST_MPI
    this->mWtime = this->mTimer.elapsed();
#else
    this->mWtime = (MPI_Wtime() - this->mWtime);
#endif
#else
#ifdef _OPENMP
    this->mWtime = (omp_get_wtime() - this->mWtime);
#else
    struct timeval ct;
    gettimeofday(&ct, NULL);
    this->mWtime = static_cast<double>(ct.tv_sec)
                   + (static_cast<double>(ct.tv_usec) / 1000000.0)
                   - this->mWtime;
#endif
#endif
    return this->mWtime;
}
/*----------------------------------------------------------------------------*/
#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
template <class Archive>
void Timer::serialize(Archive& ar, const unsigned int version)
{
    if (1 <= version)
    {
        ar&(this->mWtime);
        ar&(this->mTimer);
    }
}
#endif

/*******************************************************************************
 * INLINED OPERATORS (FREE FUNCTIONS)
 ******************************************************************************/

/*******************************************************************************
 * FREE METHODS WHICH ARE FRIENDS OF THIS CLASS.
 ******************************************************************************/
inline void swap(Timer& first, Timer& second)
{
    std::swap(first.mWtime, second.mWtime);
#ifdef SCAFES_HAVE_BOOST_MPI
    boost_mpi::timer tmp;
    tmp = first.mTimer;
    first.mTimer = second.mTimer;
    second.mTimer = tmp;
#else
    std::swap(first.mTimer, second.mTimer);
#endif
}

} // End of namespace. //


/*******************************************************************************
 ******************************************************************************/
#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
/** Boost serialization version of class \c Timer. */
BOOST_CLASS_VERSION(ScaFES::Timer, 2);
#endif

#endif
