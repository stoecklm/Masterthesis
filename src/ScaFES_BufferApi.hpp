/* ScaFES
 * Copyright (c) 2011-2016, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file ScaFES_BufferApi.hpp
 *  @brief Contains the class template Buffer.
 */

#ifndef SCAFES_BUFFERAPI_HPP_
#define SCAFES_BUFFERAPI_HPP_

#include "ScaFES_Config.hpp"
#include "ScaFES_Communicator.hpp"

#include <vector>

#ifdef SCAFES_HAVE_BOOST
#include <boost/version.hpp>
#endif

#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
#include <boost/version.hpp>
#if BOOST_VERSION < 105900
    #include <boost/serialization/pfto.hpp>
#endif
#include <boost/serialization/vector.hpp>
#endif

#ifdef VTRACE
#include "vt_user.h"
#endif

namespace ScaFES
{

/*******************************************************************************
 ******************************************************************************/
/** \class Buffer
 * @brief The class \c Buffer provides buffers for sending and receiving values
 * from one grid partition to a neighboured one.
 *
 * The Boost.MPI skeleton concept can be used in order to improve the
 * performance of the data exchange.
 *
 * \remark
 * File \c ScaFES_BufferDummy.hpp will be created by shorten this file.
 * Type of member variable related to \c boost::mpi will be changed to \c int
 * in order to shorten as much as necessary, only.
 *
 * \remark
 * Boost.MPI cannot send and receive data of type
 * \code std::complex<double> \endcode.
 * Thus, one should use the ScaFES class \c Complex.
 */
template <typename TT> class Buffer
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
    /** Creates the default constructor. */
    Buffer();

    /** Creates own constructor. */
    Buffer(ScaFES::Communicator const& myWorld, int const& idNeighPartition,
           int const& sizeBuffer, bool const& useSkeletonConcept);

    /** Creates own copy constructor. */
    Buffer(Buffer<TT> const&);

    /** Creates own assignment operator. */
    Buffer& operator=(Buffer<TT> rhs);

    /** Creates own destructor. */
    ~Buffer();

    /*----------------------------------------------------------------------
    | GETTER METHODS.
    ----------------------------------------------------------------------*/
    /** Returns the MPI communicator. */
    ScaFES::Communicator const& myWorld() const;

    /** Returns the identifier of the neighboured grid partition. */
    int const& idNeighPartition() const;

    /** Returns the number of the elements within the buffers. */
    int const& sizeBuffer() const;

    /** Returns if the Boost.MPI skeleton concept will be used or not. */
    bool const& useSkeletonConcept() const;

    /** Returns the buffer for storing the send values. */
    std::vector<TT> const& dataReceived() const;

    /** Returns the buffer for storing the receive values. */
    std::vector<TT> const& dataSent() const;

    /** Returns the value of the buffer at a given index. */
    TT const& elemData(int const& idx) const;

    /*----------------------------------------------------------------------
    | SETTER METHODS.
    ----------------------------------------------------------------------*/
    /** Sets a given value to a given position in the send buffer.
     * \param memoryPos Given position.
     * \param elemData Given value.
     */
    void setElemData(unsigned long int const& memoryPos, TT const& elemData);

    /*----------------------------------------------------------------------
    | COMPARISON METHODS.
    ----------------------------------------------------------------------*/
    /** Tests if two buffers are equal. */
    bool operator==(Buffer<TT> const& rhs) const;

    /** Tests if two buffers are not equal. */
    bool operator!=(Buffer<TT> const& rhs) const;

    /*----------------------------------------------------------------------
    | WORK METHODS.
    ----------------------------------------------------------------------*/
    /** Sends the values stored in the send buffer to the given neighbour
     * partition. */
    void sendValues();

    /** Receives the values of send buffer from the given neighbour */
    void receiveValues();

    /** Waits until all communication calls have bee finished. */
    void waitAll();

#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
    /** Serializes this class. */
    template <class Archive>
    void serialize(Archive& ar, unsigned int const version);
#endif

    /*----------------------------------------------------------------------
    | FREE METHODS WHICH ARE FRIENDS OF THIS CLASS.
    ----------------------------------------------------------------------*/
    /** Swaps two objecty by applying the copy-and-swap idiom. */
    template <typename SS>
    friend void swap(Buffer<SS>& first, Buffer<SS>& second);

private:
    /*----------------------------------------------------------------------
    | MEMBER VARIABLES.
    ----------------------------------------------------------------------*/
    /** Boost.MPI communicator. */
    ScaFES::Communicator mMyWorld;

    /** Identifiers of the (neighboured) grid partitions
    * to which data are sent and received. */
    int mIdNeighPartition;

    /** Size of underlying grid of memory = Number of elements for buffer.*/
    int mSizeBuffer;

    /** Uses the Boost.MPI skeleton concept? */
    bool mUseSkeletonConcept;

    /** Buffer for storing the send values. */
    std::vector<TT> mDataSent;

    /** Buffer for storing the receive values. */
    std::vector<TT> mDataReceived;

/** List of interfaces to the neighbours. */
// std::vector<Communicator> vCommunicator;

#ifdef SCAFES_HAVE_BOOST_MPI
    /** Content of the send buffer that is separated from the structure. */
    boost_mpi::content mvContSend;

    /** Content of receive buffer that is separated from the structure. */
    boost_mpi::content mvContRecv;
#endif

#ifdef SCAFES_HAVE_BOOST_MPI
    /** List of request states of sending and receiving. */
    std::vector<boost_mpi::request> mReqs;
#else
    /** List of request states of sending and receiving. */
    std::vector<MPI_Request> mReqs;
#endif
}; // End of class. //

/*******************************************************************************
 * FREE METHODS.
 ******************************************************************************/
/** Swaps two objects of the class \c Buffer. */
template <typename TT> void swap(Buffer<TT>& first, Buffer<TT>& second);

/*******************************************************************************
 * LIFE CYCLE METHODS.
 ******************************************************************************/
template <typename TT>
inline Buffer<TT>::Buffer()
: mMyWorld(), mIdNeighPartition(-1), mSizeBuffer(0), mUseSkeletonConcept(false),
  mDataSent(), mDataReceived(), mvContSend(), mvContRecv(), mReqs()
{
    mReqs.resize(2);
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline Buffer<TT>::Buffer(ScaFES::Communicator const& myWorld,
                          int const& idNeighPartition, int const& sizeBuffer,
                          bool const& useSkeletonConcept)
: mMyWorld(myWorld), mIdNeighPartition(idNeighPartition),
  mSizeBuffer(sizeBuffer), mUseSkeletonConcept(useSkeletonConcept), mDataSent(),
  mDataReceived(), mvContSend(), mvContRecv(), mReqs()
{
    // At the beginning, the content of the buffers should be empty.
    this->mDataSent.resize(sizeBuffer);
    this->mDataReceived.resize(sizeBuffer);
    this->mReqs.resize(2);

#ifdef VTRACE
    VT_TRACER("Buffer()");
#endif
    if (this->useSkeletonConcept())
    {
#ifdef SCAFES_HAVE_BOOST_MPI
        boost_mpi::request reqs[2];
        reqs[0] = this->myWorld().myWorld().isend(this->idNeighPartition(), 0,
                                        boost_mpi::skeleton(this->mDataSent));
        this->mvContSend = boost_mpi::get_content(mDataSent);
        reqs[1] = this->myWorld().myWorld().irecv(this->idNeighPartition(), 0,
                                        boost_mpi::skeleton(this->mDataReceived));
        this->mvContRecv = boost_mpi::get_content(this->mDataReceived);
        boost_mpi::wait_all(reqs, reqs + 2);
#else
#error "Boost.MPI is required for using the skeleton concept."
#endif
    }
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline Buffer<TT>::Buffer(Buffer<TT> const& rhs)
: mMyWorld(rhs.myWorld()), mIdNeighPartition(rhs.mIdNeighPartition),
  mSizeBuffer(rhs.sizeBuffer()), mUseSkeletonConcept(rhs.useSkeletonConcept()),
  mDataSent(rhs.dataSent()), mDataReceived(rhs.dataReceived()),
  mvContSend(rhs.mvContSend), mvContRecv(rhs.mvContRecv), mReqs(rhs.mReqs)
{
}
/*----------------------------------------------------------------------------*/
template <typename TT> inline Buffer<TT>& Buffer<TT>::operator=(Buffer<TT> rhs)
{
    swap(*this, rhs);
    return *this;
}
/*----------------------------------------------------------------------------*/
template <typename TT> inline Buffer<TT>::~Buffer()
{
}

/*******************************************************************************
 * GETTER METHODS.
 ******************************************************************************/
template <typename TT>
inline ScaFES::Communicator const& Buffer<TT>::myWorld() const
{
    return this->mMyWorld;
}
/*----------------------------------------------------------------------------*/
template <typename TT> inline int const& Buffer<TT>::idNeighPartition() const
{
    return this->mIdNeighPartition;
}
/*----------------------------------------------------------------------------*/
template <typename TT> inline int const& Buffer<TT>::sizeBuffer() const
{
    return this->mSizeBuffer;
}
/*----------------------------------------------------------------------------*/
template <typename TT> inline bool const& Buffer<TT>::useSkeletonConcept() const
{
    return this->mUseSkeletonConcept;
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline std::vector<TT> const& Buffer<TT>::dataSent() const
{
    return this->mDataSent;
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline std::vector<TT> const& Buffer<TT>::dataReceived() const
{
    return this->mDataReceived;
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline TT const& Buffer<TT>::elemData(int const& idx) const
{
    return this->mDataReceived[idx];
}

/*******************************************************************************
 * SETTER METHODS.
 ******************************************************************************/
template <typename TT>
inline void Buffer<TT>::setElemData(unsigned long int const& memoryPos,
                                    TT const& elemData)
{
    this->mDataSent[memoryPos] = elemData;
}

/*******************************************************************************
 * COMPARISON METHODS.
 ******************************************************************************/
template <typename TT>
inline bool Buffer<TT>::operator==(Buffer<TT> const& rhs) const
{
    bool isEqual = true;

    if (this->idNeighPartition() != rhs.idNeighPartition())
    {
        isEqual = false;
    }
    if (this->sizeBuffer() != rhs.sizeBuffer())
    {
        isEqual = false;
    }
    if (this->dataSent() != rhs.dataSent())
    {
        isEqual = false;
    }
    if (this->dataReceived() != rhs.dataReceived())
    {
        isEqual = false;
    }
    return isEqual;
}
/*----------------------------------------------------------------------------*/
template <typename TT>
inline bool Buffer<TT>::operator!=(Buffer<TT> const& rhs) const
{
    return !(*this == rhs);
}

/*******************************************************************************
 * WORK METHODS.
 ******************************************************************************/
template <typename TT> void Buffer<TT>::sendValues()
{
#ifdef VTRACE
    VT_TRACER("Buffer::send()");
#endif
    if (this->useSkeletonConcept())
    {
#ifdef SCAFES_HAVE_BOOST_MPI
        this->mvContSend = boost_mpi::get_content(mDataSent);
        this->mReqs[0] = this->myWorld().myWorld().isend(this->mIdNeighPartition, 0, this->mvContSend);
#else
#error "Boost.MPI is required for using the skeleton concept."
#endif
    }
    else
    {
        this->mReqs[0] = this->myWorld().myWorld().isend(this->mIdNeighPartition, 0, this->mDataSent);
    }
}
/*----------------------------------------------------------------------------*/
template <typename TT> void Buffer<TT>::receiveValues()
{
#ifdef VTRACE
    VT_TRACER("Buffer::receive()");
#endif
    if (this->useSkeletonConcept())
    {
#ifdef SCAFES_HAVE_BOOST_MPI
        this->mvContRecv = boost_mpi::get_content(this->mDataReceived);
        this->mReqs[1] = this->myWorld().myWorld().irecv(this->mIdNeighPartition, 0, this->mvContRecv);
#else
#error "Boost.MPI is required for using the skeleton concept."
#endif
    }
    else
    {
        this->mReqs[1] = this->myWorld().myWorld().irecv(this->mIdNeighPartition, 0, this->mDataReceived);
    }
}
/*----------------------------------------------------------------------------*/
template <typename TT> void Buffer<TT>::waitAll()
{
#ifdef VTRACE
    VT_TRACER("Buffer::waitAll()");
#endif
#ifdef SCAFES_HAVE_BOOST_MPI
    boost_mpi::wait_all(this->mReqs.begin(), this->mReqs.end());
#else
    MPI_Waitall(2, this->mReqs, MPI_STATUSES_IGNORE);
#endif
}
/*----------------------------------------------------------------------------*/
#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
template <typename TT>
template <class Archive>
inline void Buffer<TT>::serialize(Archive& ar, unsigned int const version)
{
    if (1 <= version)
    {
        ar& (this->mMyWorld);
        ar& (this->mIdNeighPartition);
        ar& (this->mSizeBuffer);
        ar& (this->mUseSkeletonConcept);
        ar& (this->mDataSent);
        ar& (this->mDataReceived);
        ar& (this->mvContSend);
        ar& (this->mvContRecv);
        ar& (this->mReqs);
    }
}
#endif

/*******************************************************************************
 * FREE METHODS.
 ******************************************************************************/
template <typename TT> inline void swap(Buffer<TT>& first, Buffer<TT>& second)
{
    std::swap(first.mMyWorld, second.mMyWorld);
    std::swap(first.mIdNeighPartition, second.mIdNeighPartition);
    std::swap(first.mSizeBuffer, second.mSizeBuffer);
    std::swap(first.mUseSkeletonConcept, second.mUseSkeletonConcept);
    std::swap(first.mDataSent, second.mDataSent);
    std::swap(first.mDataReceived, second.mDataReceived);
    std::swap(first.mvContSend, second.mvContSend);
    std::swap(first.mvContRecv, second.mvContRecv);
    std::swap(first.mReqs, second.mReqs);
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
    template <typename TT> struct version<ScaFES::Buffer<TT>>
    {
#if BOOST_VERSION < 105900
        /** Sets the version number for serialization. */
        BOOST_STATIC_CONSTANT(unsigned BOOST_PFTO int, value = 2);
#endif
    };
} // namespace serialization
} // namespace boost
#endif

#endif
