/* ScaFES
 * Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file ScaFES_Environment.hpp
 *  @brief Contains the class Environment.
 */

#ifndef SCAFES_ENVIRONMENT_HPP_
#define SCAFES_ENVIRONMENT_HPP_

#include "ScaFES_Config.hpp"
// TODO
#ifdef SCAFES_HAVE_BOOST_MPI
#define SCAFES_HAVE_MPI 1
#endif

#include <cstdlib>
#include <iostream>

#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
#include <boost/serialization/serialization.hpp>
#endif

#ifdef SCAFES_HAVE_MPI
#include <mpi.h>
#endif

#ifdef SCAFES_HAVE_BOOST_MPI
#include <boost/mpi.hpp>
#endif

#ifdef SCAFES_HAVE_BOOST_MPI
namespace boost_mpi = boost::mpi;
#endif

#include "ScaFES_Communicator.hpp"

namespace ScaFES
{

class Environment
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
#ifdef SCAFES_HAVE_BOOST_MPI
     /** Creates own default constructor. */
    Environment();
#endif

    /** Creates own constructor from command line parameters. */
    Environment(int argc, char* argv[]);

    /** Creates own copy constructor. */
    /*Environment(Environment const&);*/

    /** Creates default destructor. */
    ~Environment();

    /** Creates own assignment operator. */
    Environment& operator=(const Environment& rhs) = default;

    /** Aborts all MPI processes. */
    void abort(ScaFES::Communicator const& myWorld, int errcode);

    /*----------------------------------------------------------------------
    | FREE METHODS.
    ----------------------------------------------------------------------*/

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
#ifdef SCAFES_HAVE_MPI
#ifdef SCAFES_HAVE_BOOST_MPI
    /** Boost MPI communicator. */
    boost_mpi::environment mMyEnv;
#endif
#endif
};

/*******************************************************************************
 * FREE METHODS.
 ******************************************************************************/

/*******************************************************************************
 * LIFE CYCLE METHODS.
 ******************************************************************************/
#ifdef SCAFES_HAVE_BOOST_MPI
inline Environment::Environment()
: mMyEnv()
{ }
#endif
/*----------------------------------------------------------------------------*/
inline Environment::Environment(int argc, char* argv[])
#ifdef SCAFES_HAVE_MPI
#ifdef SCAFES_HAVE_BOOST_MPI
:
  mMyEnv(argc, argv)
#endif
#endif
{
#ifdef SCAFES_HAVE_MPI
#ifdef SCAFES_HAVE_BOOST_MPI
#else
    int isInitialized = 0;

    MPI_Initialized(&isInitialized);
    if (!isInitialized) {
        MPI_Init(&argc, &argv);
    }
#endif
#else
    std::ignore = argc;
    std::ignore = argv;
#endif
}
/*----------------------------------------------------------------------------*/
inline Environment::~Environment()
{
#ifdef SCAFES_HAVE_MPI
#ifdef SCAFES_HAVE_BOOST_MPI
#else
    int isFinalized = 0;

    MPI_Finalized(&isFinalized);
    if (!isFinalized) {
        MPI_Finalize();
    }
#endif
#endif
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
template <class Archive>
void Environment::serialize(Archive& ar, const unsigned int version)
{
    if (1 <= version)
    {
    }
}
#endif
/*----------------------------------------------------------------------------*/
inline void Environment::abort(ScaFES::Communicator const& myComm, int errcode)
{
#ifdef SCAFES_HAVE_MPI
#ifdef SCAFES_HAVE_BOOST_MPI
    std::ignore = myComm;
    this->mMyEnv.abort(errcode);
#else
    MPI_Abort(myComm.myWorld(), errcode);
#endif
#else
    std::ignore = myComm;
    exit(errcode);
#endif
}

/*******************************************************************************
 * FREE METHODS WHICH ARE FRIENDS OF THIS CLASS.
 ******************************************************************************/
} // End of namespace. //


/*******************************************************************************
 ******************************************************************************/
#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
/** Boost serialization version of class \c Environment. */
BOOST_CLASS_VERSION(ScaFES::Environment, 2);
#endif

#endif
