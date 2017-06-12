/* ScaFES
 * Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file ScaFES_Communicator.hpp
 *  @brief Contains the class Communicator.
 */

#ifndef SCAFES_COMMUNICATOR_HPP_
#define SCAFES_COMMUNICATOR_HPP_

#include "ScaFES_Config.hpp"
// TODO
#ifdef SCAFES_HAVE_BOOST_MPI
#define SCAFES_HAVE_MPI 1
#endif

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>
#include <tuple>
#include <string>

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

namespace ScaFES
{

class Communicator
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
    Communicator();

    /** Creates own constructor from command line parameters. */
    Communicator(int argc, char* argv[]);

    /** Creates own copy constructor. */
    Communicator(Communicator const&);

    /** Creates default destructor. */
    ~Communicator() = default;

    /** Creates own assignment operator. */
    Communicator& operator=(Communicator rhs);

    /*----------------------------------------------------------------------
    | GETTER METHODS.
     ----------------------------------------------------------------------*/
#ifdef SCAFES_HAVE_MPI
#ifdef SCAFES_HAVE_BOOST_MPI
    /** Returns the Boost MPI communicator. */
    boost_mpi::communicator const& myWorld() const;
#else
    /** Returns the MPI communicator. */
    MPI_Comm const& myWorld() const;
#endif
#else
    /** Returns the communicator. */
    int const& myWorld() const;
#endif

    /** Returns the MPI rank. */
    int rank() const;

    /** Returns the number of MPI processes. */
    int size() const;

    /** Barrier for all MPI processes. */
    void barrier() const;

    /** Collects a scalar from each grid partition and sums it up
     * at the given MPI process. */
    template<typename CT>
    void collectScalar(CT& scalarGlobal, CT const& scalarLocal,
                       int const& root) const;

    /** Collects a vector from each grid partition and sums it up
     * at the given MPI process. */
    template<typename CT>
    void collectVector(std::vector<CT>& vectGlobal,
                       std::vector<CT> const& vectLocal,
                       int const& nElems,
                       int const& root) const;

    /** Collects an one-dimensional MPI local array from all grid
     * partitions in an one-dimensional MPI global array. */
    template<typename CT>
    void assembleSubvectors(std::vector<CT>& arrayGlobal,
                        std::vector<CT>& arrayLocal, int** const& rcounts,
                        int** const& displs, std::vector<int> const& nVectLocal,
                        int const& nTotalLocal, int const& nDataFields,
                        int const& root) const;

    /** Distributes an one-dimensional MPI-global arrary
    *  to one-dimensional MPI local arrays
    *  defined on the grid partitions. */
    template<typename CT>
    void partitionGlobalvector(std::vector<CT>& arrayLocal,
                           std::vector<CT>& arrayGlobal, int** const& rcounts,
                           int** const& displs,
                           std::vector<int> const& nVectLocal,
                           int const& nTotalLocal, int const& nDataFields,
                           int const& root) const;

    /** Computes a matrix containing the rcounts from each data field
     * which should be optimized. */
    void compMatRcounts(
        int**& rcounts,
        std::vector<int> const& nParamsToOptimizeLocalOfDf,
        int const& nDataFields
    ) const;

    /** Computes a matrix containing the displacements from each data field
     * which should be optimized. */
    void compMatDisplacements(
        int**& displs,
        int** const& rcounts,
        int const& nDataFields
    ) const;

    /*----------------------------------------------------------------------
    | FREE METHODS.
    ----------------------------------------------------------------------*/
    /** Swaps the members of two communicators. */
    friend void swap(ScaFES::Communicator& first,
                     ScaFES::Communicator& second);

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
    boost_mpi::communicator mMyWorld;
#else
    /** MPI communicator. */
    MPI_Comm mMyWorld;
#endif
#else
    int mMyWorld;
#endif
};

/*******************************************************************************
 * FREE METHODS.
 ******************************************************************************/
/** Swaps two objects of the class \c Communicator. */
void swap(ScaFES::Communicator& first, ScaFES::Communicator& second);
/*----------------------------------------------------------------------------*/
/** Broadcast values of type \c TT from a given rank. */
template<typename TT>
inline void broadcast(ScaFES::Communicator const& comm, TT& value,
                      int const& rank);
/*----------------------------------------------------------------------------*/
/** Broadcast vector of type \c TT from a given rank. */
template<typename TT>
inline void broadcastVector(ScaFES::Communicator const& comm,
                            std::vector<TT>& value,
                            int const& rank);
/*----------------------------------------------------------------------------*/
/** Broadcast array of characters. from a given rank. */
inline void broadcastChars(ScaFES::Communicator const& comm, char* value,
                           int const& nChars, int const& rank);

/*******************************************************************************
 * LIFE CYCLE METHODS.
 ******************************************************************************/
inline Communicator::Communicator()
: 
  mMyWorld()
{    
 }
/*----------------------------------------------------------------------------*/
inline Communicator::Communicator(int argc, char* argv[])
:
  mMyWorld()
{
#ifdef SCAFES_HAVE_MPI
#ifdef SCAFES_HAVE_BOOST_MPI
#else
    MPI_Comm_dup(MPI_COMM_WORLD, &(this->mMyWorld) );
#endif
#endif
    std::ignore = argc;
    std::ignore = argv;
}
/*----------------------------------------------------------------------------*/
inline Communicator::Communicator(Communicator const& rhs)
: mMyWorld(rhs.myWorld())
{
#ifdef SCAFES_HAVE_MPI
#ifdef SCAFES_HAVE_BOOST_MPI
#else
    MPI_Comm_dup(rhs.myWorld(), &(this->mMyWorld) );
    
#endif
#endif
}
/*----------------------------------------------------------------------------*/
inline Communicator& Communicator::operator=(Communicator rhs)
{
    swap(*this, rhs);
    return *this;
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
#ifdef SCAFES_HAVE_MPI
#ifdef SCAFES_HAVE_BOOST_MPI
inline boost_mpi::communicator const& Communicator::myWorld() const
{
    return this->mMyWorld;
}
#else
/*----------------------------------------------------------------------------*/
inline MPI_Comm const& Communicator::myWorld() const
{
    //std::cout << "comm myworld " << std::endl;
    return this->mMyWorld;
}
#endif
#else
/*----------------------------------------------------------------------------*/
inline int const& Communicator::myWorld() const
{
    return this->mMyWorld;
}
#endif
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
inline int Communicator::rank() const
{
#ifdef SCAFES_HAVE_MPI
#ifdef SCAFES_HAVE_BOOST_MPI
    return this->myWorld().rank();
#else
    int rk = 0;
    MPI_Comm_rank(this->myWorld(), &rk);
    return rk;
#endif
#else
    return 0;
#endif
}
/*----------------------------------------------------------------------------*/
inline int Communicator::size() const
{
#ifdef SCAFES_HAVE_MPI
#ifdef SCAFES_HAVE_BOOST_MPI
    return this->myWorld().size();
#else
    int sz = 0;
    MPI_Comm_size(this->myWorld(), &sz);
    return sz;
#endif
#else
    return 1;
#endif
}
/*----------------------------------------------------------------------------*/
inline void Communicator::barrier() const
{
#ifdef SCAFES_HAVE_MPI
#ifdef SCAFES_HAVE_BOOST_MPI
    this->myWorld().barrier();
#else
    MPI_Barrier(this->myWorld());
#endif
#endif
}
#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
/*----------------------------------------------------------------------------*/
template <class Archive>
void Communicator::serialize(Archive& ar, const unsigned int version)
{
    if (1 <= version)
    {
        ar&(this->mMyWorld);
    }
}
#endif

/*******************************************************************************
 * INLINED OPERATORS (FREE FUNCTIONS)
 ******************************************************************************/
template<typename TT>
inline void broadcast(ScaFES::Communicator const& comm, TT& value, int const& rank)
{
    if (1 < comm.size())
    {
#ifdef SCAFES_HAVE_MPI
#ifdef SCAFES_HAVE_BOOST_MPI
        boost::mpi::broadcast(comm.myWorld(), value, rank);
#else
        // TODO: Has to be implemented for broadcasting all ScaFES parameters.
        // MPI_Bcast(void *buffer, int count, MPI_Datatype datatype,
        // rank, comm.myWorld());
        std::ignore = value;
        std::ignore = rank;
#endif
#else
    std::ignore = value;
    std::ignore = rank;
#endif
    }
    // Else: Do nothing.
}
/*----------------------------------------------------------------------------*/
template<typename TT>
inline void broadcastVector(ScaFES::Communicator const& comm,
                            std::vector<TT>& value, int const& rank)
{
    if (1 < comm.size())
    {
#ifdef SCAFES_HAVE_MPI
#ifdef SCAFES_HAVE_BOOST_MPI
        boost::mpi::broadcast(comm.myWorld(), value, rank);
#else
        // TODO: Has to be implemented for broadcasting all ScaFES parameters.
        // MPI_Bcast(void *buffer, int count, MPI_Datatype datatype,
        // rank, comm.myWorld());
        std::ignore = value;
        std::ignore = rank;
#endif
#else
    std::ignore = value;
    std::ignore = rank;
#endif
    }
    // Else: Do nothing.
}
/*----------------------------------------------------------------------------*/
inline void broadcastChars(ScaFES::Communicator const& comm, char* value,
                           int const& nChars, int const& rank)
{
    if (1 < comm.size())
    {
#ifdef SCAFES_HAVE_MPI
#ifdef SCAFES_HAVE_BOOST_MPI
        MPI_Bcast(value, nChars, MPI_CHAR, rank, comm.myWorld());
#endif
#else
    std::ignore = value;
    std::ignore = rank;
    std::ignore = nChars;
#endif
    }
    // Else: Do nothing.
}
/*----------------------------------------------------------------------------*/
template <typename CT>
inline void Communicator::collectVector(std::vector<CT>& vectGlobal,
                                        std::vector<CT> const& vectLocal,
                                        int const & nElems,
                                        int const& root) const
{
    if (1 < this->size())
    {
#ifdef SCAFES_HAVE_MPI
#ifdef SCAFES_HAVE_BOOST_MPI
        const int N_COMPS = nElems;
        double* sendbuff = new double[nElems];
        double* recvbuff = new double[nElems];
        for (int ii = 0; ii < nElems; ++ii)
        {
            sendbuff[ii] = vectGlobal[ii];
            recvbuff[ii] = 0.0;
        }
        MPI_Reduce(sendbuff, recvbuff, N_COMPS, MPI_DOUBLE, MPI_SUM, root,
                   this->myWorld());
        for (int ii = 0; ii < nElems; ++ii)
        {
            vectGlobal[ii] = recvbuff[ii];
        }
        delete[] sendbuff;
        delete[] recvbuff;
#else
        // TODO: Check this implementation.
//         for (int ii = 0; ii < nElems, ++ii) {
//             boost::mpi::reduce(this->myWorld(), vectLocal, vectGlobal,
//                                std::plus<CT>(), root);
//         }

        const int N_COMPS = nElems;
        double* sendbuff = new double[nElems];
        double* recvbuff = new double[nElems];
        for (int ii = 0; ii < nElems; ++ii)
        {
            sendbuff[ii] = vectGlobal[ii];
            recvbuff[ii] = 0.0;
        }
        MPI_Reduce(sendbuff, recvbuff, N_COMPS, MPI_DOUBLE, MPI_SUM, root,
                   this->myWorld());
        for (int ii = 0; ii < nElems; ++ii)
        {
            vectGlobal[ii] = recvbuff[ii];
        }
        delete[] sendbuff;
        delete[] recvbuff;
#endif
#else
    std::ignore = root;
#endif
    }
    else
    {
        for (int ii = 0; ii < nElems; ++ii)
        {
            vectGlobal[ii] = vectLocal[ii];
        }
    }
}

/*----------------------------------------------------------------------------*/
template <typename CT>
inline void Communicator::collectScalar(CT& scalarGlobal,
                                        CT const& scalarLocal,
                                        int const& root) const
{
    scalarGlobal = 0.0; // default value.

    if (1 < this->size())
    {
#ifdef SCAFES_HAVE_MPI
#ifdef SCAFES_HAVE_BOOST_MPI
        boost::mpi::reduce(
            this->myWorld(),
            scalarLocal,
            scalarGlobal,
            std::plus<CT>(),
            root
        );
#else
        const int N_COMPS = 1;
        double sendbuff[N_COMPS];
        double recvbuff[N_COMPS];
        sendbuff[0] = scalarLocal;
        recvbuff[0] = scalarGlobal;
        MPI_Reduce(sendbuff, recvbuff, N_COMPS, MPI_DOUBLE, MPI_SUM, root,
                   this->myWorld());
        scalarGlobal = recvbuff[0];
#endif
#else
    std::ignore = root;
#endif
    }
    else
    {
        scalarGlobal = scalarLocal;
    }
}
/*----------------------------------------------------------------------------*/
template <typename CT>
inline void Communicator::assembleSubvectors(
    std::vector<CT>& arrayGlobal, std::vector<CT>& arrayLocal,
    int** const& rcounts, int** const& displs,
    std::vector<int> const& nVectLocal,
    int const& nTotalLocal, int const& nDataFields, int const& root) const
{
    if (1 < this->size())
    {
#ifdef SCAFES_HAVE_MPI
        CT* sendbufferStart = arrayLocal.data();
        const int* ptr2nVectLocal = nVectLocal.data();
        for (int jj = 0; jj < nDataFields; ++jj)
        {
#ifdef SCAFES_HAVE_BOOST_MPI
            //TODO: Check for correctness.
            //boost_mpi::gather(this->myWorld, arrayLocal,
            //                  arrayGlobal, root);
            MPI_Gatherv(sendbufferStart, ptr2nVectLocal[jj], MPI_DOUBLE,
                        arrayGlobal.data(), rcounts[jj], displs[jj],
                        MPI_DOUBLE,
                        root, this->myWorld());
#else
            MPI_Gatherv(sendbufferStart, ptr2nVectLocal[jj], MPI_DOUBLE,
                        arrayGlobal.data(), rcounts[jj], displs[jj],
                        MPI_DOUBLE,
                        root, this->myWorld());
#endif
            this->barrier();
            sendbufferStart += nVectLocal[jj];
        }
#else
    std::ignore = rcounts;
    std::ignore = displs;
    std::ignore = nVectLocal;
    std::ignore = nTotalLocal;
    std::ignore = nDataFields;
    std::ignore = root;
#endif
    }
    else
    {
        for (int ii = 0; ii < nTotalLocal; ++ii)
        {
            arrayGlobal[ii] = arrayLocal[ii];
        }
    }
}
/*----------------------------------------------------------------------------*/
template <typename CT>
inline void Communicator::partitionGlobalvector(
    std::vector<CT>& arrayLocal, std::vector<CT>& arrayGlobal,
    int** const& rcounts, int** const& displs,
    std::vector<int> const& nVectLocal,
    int const& nTotalLocal, int const& nDataFields, int const& root) const
{
    if (1 < this->size())
    {
#ifdef SCAFES_HAVE_MPI
        // Distribute arrays from process for optimization to
        // all MPI processes.
        // Case: Current rank == root
        // ---> Send local array portion to all other processes.
        // Case: Current rank is NOT root.
        // ---> Receive local array from root.
        CT* ptrToArrayGlobal = arrayGlobal.data();
        CT* recvBufferStart = arrayLocal.data();
        const int* ptr2nVectLocal = nVectLocal.data();
        for (int jj = 0; jj < nDataFields; ++jj)
        {
#ifdef SCAFES_HAVE_BOOST_MPI
            //TODO: Check for correctness.
            //boost_mpi::scatter(this->myWorld, arrayGlobal,
            //                   arrayLocal, root);
            MPI_Scatterv(ptrToArrayGlobal, rcounts[jj], displs[jj], MPI_DOUBLE,
                         recvBufferStart, ptr2nVectLocal[jj], MPI_DOUBLE,
                         root,
                         this->myWorld());
#else
            MPI_Scatterv(ptrToArrayGlobal, rcounts[jj], displs[jj], MPI_DOUBLE,
                         recvBufferStart, ptr2nVectLocal[jj], MPI_DOUBLE,
                         root,
                         this->myWorld());
#endif
            this->barrier();
            recvBufferStart += nVectLocal[jj];
        }
#else
    std::ignore = rcounts;
    std::ignore = displs;
    std::ignore = nVectLocal;
    std::ignore = nTotalLocal;
    std::ignore = nDataFields;
    std::ignore = root;
#endif
    }
    else
    {
        for (int ii = 0; ii < nTotalLocal; ++ii)
        {
            arrayLocal[ii] = arrayGlobal[ii];
        }
    }
}
/*----------------------------------------------------------------------------*/
inline void Communicator::compMatRcounts(
    int**& rcounts,
    std::vector<int> const& nParamsToOptimizeLocalOfDf,
    int const& nDataFields
) const
{
    int nProcesses = this->size();
    rcounts = new int* [nDataFields];
    for (int jj = 0; jj < nDataFields; ++jj)
    {
        rcounts[jj] = new int[nProcesses];
    }

    for (int jj = 0; jj < nDataFields; ++jj)
    {
        // Collect portions of local arrays and send it to the
        // MPI process for optimization "root".
        // Exchange number of nodes of each process (=grid partition).
        // Case: Current rank != root
        // ---> Send local portion of array to root.
        // if (root != this->myRank()) {
        // this->params().myWorld().send(root, 0,
        //                               nParamsToOptimizeLocalOfDf[jj]);
        //}
        int nUnknownsLocalCurr;
        // Case: Current rank == root
        // ---> Receive local portion of array from all other processes.
        // if (root == this->myRank()) {
        for (int ii = 0; ii < nProcesses; ++ii)
        {
            // if (root != ii) {
            nUnknownsLocalCurr = nParamsToOptimizeLocalOfDf[jj];
#ifdef SCAFES_HAVE_MPI
#ifdef SCAFES_HAVE_BOOST_MPI
            boost::mpi::broadcast(this->myWorld(), nUnknownsLocalCurr,
                                  ii);
#endif
#endif
            // this->params().myWorld().recv(ii, 0, nUnknownsLocalCurr);
            rcounts[jj][ii] = nUnknownsLocalCurr;
            //  }
            //  else {
            //      rcounts[jj][ii] = nParamsToOptimizeLocalOfDf[jj];
            //  }
        }
        //}
        this->barrier();
    }
}
/*----------------------------------------------------------------------------*/
inline void Communicator::compMatDisplacements(
    int**& displs,
    int** const& rcounts,
    int const& nDataFields
) const
{
    int nProcesses = this->size();
    displs = new int*[nDataFields];
    for (int jj = 0; jj < nDataFields; ++jj)
    {
        displs[jj] = new int[nProcesses];
    }

    int nCurrDf = 0;
    int nPreviousDfs = 0;
    for (int jj = 0; jj < nDataFields; ++jj)
    {
        for (int ii = 0; ii < nProcesses; ++ii)
        {
            nPreviousDfs = 0;
            for (int ll = 0; ll < jj; ++ll)
            {
                for (int kk = 0; kk < nProcesses; ++kk)
                {
                    nPreviousDfs += rcounts[ll][kk];
                }
            }
            nCurrDf = 0;
            for (int kk = 0; kk < ii; ++kk)
            {
                nCurrDf += rcounts[jj][kk];
            }
            displs[jj][ii] = nPreviousDfs + nCurrDf;
        }
    }
}


/*******************************************************************************
 * FREE METHODS WHICH ARE FRIENDS OF THIS CLASS.
 ******************************************************************************/
inline void swap(Communicator& first, Communicator& second)
{
#ifdef SCAFES_HAVE_MPI
#ifdef SCAFES_HAVE_BOOST_MPI
    boost_mpi::communicator tmp;
#else
    MPI_Comm tmp;
#endif
#endif
#ifdef SCAFES_HAVE_MPI
    tmp = first.mMyWorld;
    first.mMyWorld = second.mMyWorld;
    second.mMyWorld = tmp;
#else
    std::ignore = first;
    std::ignore = second;
#endif
}
} // End of namespace. //


/*******************************************************************************
 ******************************************************************************/
#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
/** Boost serialization version of class \c Communicator. */
BOOST_CLASS_VERSION(ScaFES::Communicator, 2);
#endif

#endif
