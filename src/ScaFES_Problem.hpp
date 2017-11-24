/* ScaFES
 * Copyright (c) 2011-2017, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file ScaFES_Problem.hpp
 *  @brief Contains the class template Problem.
 */

#ifndef SCAFES_PROBLEM_HPP_
#define SCAFES_PROBLEM_HPP_

#include "ScaFES_Config.hpp"

#include <tuple>
#include <iostream>
#include <iomanip>
#include <ios>
#include <fstream>
#include <string>
#include <cassert>
#include <vector>
#include <stdexcept>
#include <unordered_map>

extern "C" {
#ifdef _OPENMP
#include <omp.h>
#endif
}

#ifdef SCAFES_HAVE_ADOLC
#include <adolc/adolc.h>
#include <adolc/adouble.h>
#ifdef _OPENMP
#include <adolc/adolc_openmp.h>
#endif // _OPENMP
#endif // SCAFES_HAVE_ADOLC

#ifdef SCAFES_HAVE_BOOST
#include <boost/version.hpp>
#include <boost/type_traits.hpp>
#include <boost/assign/list_of.hpp>
#endif

#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
namespace boost
{
namespace serialization
{
    class access;
}
}
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/version.hpp>
#if BOOST_VERSION < 105900
   #include <boost/serialization/pfto.hpp>
#endif
#include <boost/serialization/vector.hpp>
#endif

#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
#ifdef SCAFES_HAVE_ADOLC
namespace boost
{
namespace serialization
{
    template <class Archive>
    /** Serialization for type adouble. */
    void serialize(Archive& ar, adouble& g, const unsigned int /*version*/)
    {
        double tmp = g.getValue();
        ar&(tmp);
#if defined(TAPELESS)
        adtl::ADVAL_TYPE tmp2 = g.getADValue();
        ar&(tmp2);
#endif
    }
}
}
#endif // SCAFES_HAVE_ADOLC
#endif // SCAFES_HAVE_BOOST_SERIALIZATION

#ifdef VTRACE
#include "vt_user.h"
#endif

#include "ScaFES_Ntuple.hpp"
#include "ScaFES_Parameters.hpp"
#include "ScaFES_GridSub.hpp"
#include "ScaFES_GridGlobal.hpp"
#include "ScaFES_Timer.hpp"
#include "ScaFES_DataFile.hpp"
#include "ScaFES_DataField.hpp"
#ifdef SCAFES_HAVE_CUDA
#include "ScaFES_ProblemCuda.hpp"
#endif

namespace ScaFES
{

/*******************************************************************************
 ******************************************************************************/
/**
 * \class Problem
 *  @brief Base class for initial boundary value problems (IBVP).
 *
 * The class \c Problem is designed to describe an initial boundary value
 * problem (IBVP) on a given MPI partition. An IBVP usually consists of
 * several data fields (some of them are unknown, some of them are given).
 *
 * In order to solve a particular IBVP,
 * the user has to implement an own class that is inherited by this class.
 * The derived class has to contain the following methods:
 * - \c evalInner
 * - \c evalBorder
 * - \c initInner
 * - \c initBorder
 * - \c updateInner
 * - \c updateBorder
 *
 * If the leap frog scheme is used then the user has to provide additionally
 * the following methods:
 * - \c updateInner2
 * - \c updateBorder2
 *
 * These methods must be declared public in the sub class
 * as they will be accessed within in this class via the CRTP.
 *
 * The method \c iterateOverTime computes the new iterates for all time
 * steps.
 */

/**
 * \internal
 * Case: ADOL-C enabled, synchronous MPI communication
 * One tape is created in each computation phase such that the MPI
 * communication can be synchronously.
 * 1 tape for the computation phase at all grid partition nodes.
 *
 * \n Case: TypeOne method:
 * Tape 1:
 * new1(idxNode) = defaultValue  if idxNode in SetGridAll and not explicitly
 *                               set in user-written method funcTypeOneBorder
 *                               or funcTypeOneInner.
 * new1(idxNode) = funcTypeTwoInner(idxNode)
 *                               if idxNode in SetGridInner and explicitly
 *                               set in user-written method funcTypeOneInner.
 * new1(idxNode) = funcTypeTwoBorder(idxNode)
 *                               if idxNode in SetGridBorder and explicitly
 *                               set in user-written method funcTypeOneBorder.
 * ==> new(All) = new1(All)
 *              = funcTypeOne(Border) + funcTypeOne(Inner)
 *              = funcTypeOne(All)
 *
 * \n Case: TypeTwo method:
 * Tape 1:
 * new1(idxNode) = old(idxNode)  if idxNode in SetGridAll and not explicitly
 *                               set in user-written method funcTypeTwoBorder
 *                               or funcTypeTwoInner.
 * new1(idxNode) = funcTypeTwoInner(idxNode)
 *                               if idxNode in SetGridInner and explicitly
 *                               set in user-written method funcTypeTwoInner.
 * new1(idxNode) = funcTypeTwoBorder(idxNode)
 *                               if idxNode in SetGridBorder and explicitly
 *                               set in user-written method funcTypeTwoBorder.
 * ==> new(All) = new1(All)
 *              = funcTypeOne(Border) + funcTypeOne(Inner)
 *              = funcTypeOne(All)
 *
 * \n Case: TypeThree method:
 * No MPI communication is necessary because each MPI process computes one
 * local value from values which are defined at the local grid partition nodes.
 * Tape 1:
 * new1 = 0
 * new1 = sum_(idxNode) funcTypeThreeBorder(idxNode)
 *                    with idxNode in SetGridBorder and
 *                    explicitly set in user-written method funcTypeThreeBorder.
 *        + sum_(idxNode) funcTypeThreeInner(idxNode)
 *                    with idxNode in SetGridInner and
 *                    explicitly set in user-written method funcTypeThreeInner.
 *
 * ==> new(All) = new1(All)
 *              = sum(i in border) old(i) + sum(i in inner) old(i)
 *              = sum(i in all) old(i)
 */

/**
 * \internal
 * Case: ADOL-C enabled, asynchronous MPI communication
 * Two tapes are created in each computation phase such that the MPI
 * communication can be asynchronously.
 * Two tapes will be created:
 * 1 tape for the computation phase at all border grid partition nodes.
 * 1 tape for the computation phase at all inner grid partition nodes.
 *
 * \n Case: TypeOne method:
 * Tape 1:
 * new1(idxNode) = defaultValue  if idxNode in SetGridBorder and not explicitly
 *                               set in user-written method funcTypeOneBorder.
 * new1(idxNode) = funcTypeTwoBorder(idxNode)
 *                               if idxNode in SetGridBorder and explicitly
 *                               set in user-written method funcTypeOneBorder.
 * new2(idxNode) = 0             if idxNode in SetGridInner.
 *
 * Tape 2:
 * new2(idxNode) = 0             if idxNode in SetGridBorder
 * new2(idxNode) = defaultValue  if idxNode in SetGridInner and not explicitly
 *                               set in user-written method funcTypeOneInner.
 * new2(idxNode) = funcTypeTwoInner(idxNode)
 *                               if idxNode in SetGridInner and explicitly
 *                               set in user-written method funcTypeOneInner.
 * ==> new(All) = new1(All) + new2(All) - defaultValue(All)
 *              = defaultValue(Inner) + funcTypeOne(Border)
 *                 + defaultValue(Border) + funcTypeOne(Inner)
 *                 - defaultValue(All)
 *              = funcTypeOne(All)
 *
 * \n Case: TypeTwo method:
 * Tape 1:
 * new1(idxNode) = old(idxNode)  if idxNode in SetGridBorder and not explicitly
 *                               set in user-written method funcTypeTwoBorder.
 * new1(idxNode) = funcTypeTwoBorder(idxNode)
 *                               if idxNode in SetGridBorder and explicitly
 *                               set in user-written method funcTypeTwoBorder.
 * new1(idxNode) = 0             if idxNode in SetGridInner
 *                               set in user-written method funcTypeTwoBorder.
 * Tape 2:
 * new2(idxNode) = 0             if idxNode in SetGridBorder
 * new2(idxNode) = old(idxNode)  if idxNode in SetGridInner and not explicitly
 *                               set in user-written method funcTypeTwoInner.
 * new2(idxNode) = funcTypeTwoInner(idxNode)
 *                               if idxNode in SetGridInner and explicitly
 *                               set in user-written method funcTypeTwoInner.
 * ==> new(All) = new1(All) + new2(All) - old(All)
 *              = old(Inner) + funcTypeTwo(Border)
 *                + old(Border) + funcTypeTwo(Inner) - old(All)
 *              = funcTypeTwo(All)
 *
 * \n Case: TypeThree method:
 * No MPI communication is necessary because each MPI process computes one
 * local value from values which are defined at the local grid partition nodes.
 * Tape 1:
 * new1 = 0
 * new1 = sum_(idxNode) funcTypeThreeBorder(idxNode)
 *                    with idxNode in SetGridBorder and
 *                    explicitly set in user-written method funcTypeThreeBorder.
 *
 * Tape 2:
 * new2 = 0
 * new2 = sum_(idxNode) funcTypeThreeInner(idxNode)
 *                    with idxNode in SetGridInner and
 *                    explicitly set in user-written method funcTypeThreeInner.
 *
 * ==> new(All) = new1(All) + new2(All)
 *              = sum(i in border) old(i) + sum(i in inner) old(i)
 *              = sum(i in all) old(i)
 */

template <class OWNPRBLM, typename CT = double, std::size_t DIM = 3>
class Problem
{
public:
    /*--------------------------------------------------------------------------
    | TYPE DEFINITIONS.
    --------------------------------------------------------------------------*/
    /** Re-export typename CT STL-like. */
    typedef CT value_type;

    /** Alias declaration for vector of data fields. */
    template <typename TT>
    using ScaFES_VectDf = std::vector<ScaFES::DataField<TT, DIM>>;

    /** Type definition for integer tuple. */
    typedef ScaFES::Ntuple<int, DIM> ScaFES_IntNtuple;

    /*--------------------------------------------------------------------------
    | FRIEND CLASSES.
    --------------------------------------------------------------------------*/
#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
    friend class boost::serialization::access;
#endif

    /*--------------------------------------------------------------------------
    | LIFE CYCLE METHODS.
    --------------------------------------------------------------------------*/
    /** No default constructor. */
    // Problem() = delete;

    /** Creates own constructor.
     * All data fields which are related to the underlying problem
     * have to be added in terms of an entry of the parameters of
     * type \c std::vector.
     * @param params Set of ScaFES parameters.
     * @param gg Global grid.
     * @param useLeapfrog Should the leap frog scheme be used?
     * @param nameDatafield Name of the data fields.
     * @param stencilWidth Stencil width of the data fields.
     * @param isKnownDf Is the data field are known or unknown one?
     * @param nLayers Number of layers at the global boundary.
     * @param defaultValue Default value of data fields.
     * @param writeToFile How often should the data field be written
     *                    to file.
     * @param computeError Should the Linf error between the numerical
     *                     and exact solution be computed?
     * @param geomparamsInit Initial guess of geometrical parameters.
     */
    Problem(const ScaFES::Parameters & params, const ScaFES::GridGlobal<DIM>& gg,
            const bool& useLeapfrog,
            const std::vector<std::string>& nameDatafield,
            const std::vector<int>& stencilWidth,
            const std::vector<bool>& isKnownDf,
            const std::vector<int>& nLayers = std::vector<int>(),
            const std::vector<CT>& defaultValue = std::vector<CT>(),
            const std::vector<ScaFES::WriteHowOften>& writeToFile =
                std::vector<ScaFES::WriteHowOften>(),
            const std::vector<bool>& computeError = std::vector<bool>(),
            const std::vector<CT>& geomparamsInit = std::vector<CT>());

    /** Creates compiler provided copy constructor. */
    Problem(const Problem<OWNPRBLM, CT, DIM>& /*ppp*/) = default;

    /** Creates compiler provided assignment operator. */
    Problem& operator=(const Problem<OWNPRBLM, CT, DIM>& /*ppp*/) = default;

    /** Creates own destructor. */
    virtual ~Problem() = default;

    /*--------------------------------------------------------------------------
    | GETTER METHODS.
    --------------------------------------------------------------------------*/
    /** Returns the ScaFES parameters.
     * @return ScaFES parameters.
     */
    const ScaFES::Parameters& params() const;

    /** Returns the grid size in a given direction. */
    const double& gridsize(const int& idx) const;

    /** Returns the time step size. */
    const double& tau() const;

    /** Returns the number of time steps. */
    const int& nTimesteps() const;

    /** Returns the pointer to global grid. */
    const ScaFES::GridGlobal<DIM>& globalGrid() const;

    /** Returns the node type of a given global node. */
    const short int& kind(const ScaFES_IntNtuple& idxNode) const;

    /** Returns the node type of a given global node. */
    const ScaFES::DataField<short int, DIM>& kind() const;

    /** Returns the number of nodes in each direction. */
    const ScaFES_IntNtuple& nNodes() const;

    /** Returns the number of nodes in a given direction. */
    const int& nNodes(const int& dim) const;

    /** Returns the coordinates of a given global node. */
    ScaFES::Ntuple<double, DIM>
    coordinates(const ScaFES_IntNtuple& idxNode) const;

    /** Returns the current time for a given time step. */
    double time(const int& timestep) const;

    /** Returns the global local node of a given global node which
     * resides in a given direction. */
    ScaFES_IntNtuple connect(const ScaFES_IntNtuple& idxNode, const int& dir);

    /** Returns the numerical solution of parameters. */
    const std::vector<CT>& vectParamsCurr() const;

    /** Returns the analytical solution vector of all data fields. */
    const ScaFES_VectDf<CT>& vectSol() const;

    /** Returns the value of analytical solution data field given by \c idx
     * at the grid node with node number \c idxNode.
     */
    const CT& vectSol(const int& idx, const ScaFES_IntNtuple& idxNode) const;

    /** Returns the value of the known data field given by \c idx
     * at the grid node with node number \c idxNode.
     */
    const CT& knownDf(const int& idx, const ScaFES_IntNtuple& idxNode) const;

    /*--------------------------------------------------------------------------
    | WORK METHODS.
    --------------------------------------------------------------------------*/
    /** Iterates over all time steps.
     * - Computes the new values of all data fields at all grid nodes,
     * - synchronises values with neighboured grid partitions,
     * - redefines new situation (new data fields will be the old ones
     *   in the next time step).
     * \remark It is sufficient to swap the pointers to the memories of
     * the old and the new unknown data fields.
     * */
    void iterateOverTime();

protected:
    /*--------------------------------------------------------------------------
    | TYPE DEFINITIONS + ALIAS DECLARATIONS (C++11).
    --------------------------------------------------------------------------*/
    /** Alias declaration for type1 methods. These methods must be
     *  implemented in the derived problem class \c OWNPRBLM. */
    template <typename TT>
    using ScaFES_FuncTypeOneToImplement =
        void (OWNPRBLM::*)(ScaFES_VectDf<TT>&, const std::vector<TT>&,
                           const ScaFES_IntNtuple&, const int&);

    /** Alias declaration for type2 methods. These methods must be
     *  implemented in the derived problem class \c OWNPRBLM. */
    template <typename TT>
    using ScaFES_FuncTypeTwoToImplement =
        void (OWNPRBLM::*)(ScaFES_VectDf<TT>&, const ScaFES_VectDf<TT>&,
                           const ScaFES_IntNtuple&, const int&);

    /** Alias declaration for type3 methods. These methods must be
     *  implemented in the derived problem class \c OWNPRBLM. */
    template <typename TT>
    using ScaFES_FuncTypeThreeToImplement =
        void (OWNPRBLM::*)(TT&, const ScaFES_VectDf<TT>&,
                           const ScaFES_IntNtuple&, const int&);

    /*--------------------------------------------------------------------------
    | SETTER METHODS.
    --------------------------------------------------------------------------*/
    /** Increases output level by one. */
    void increaseLevel();

    /** Decreases output level by one. */
    void decreaseLevel();

    /*--------------------------------------------------------------------------
    | GETTER METHODS.
    --------------------------------------------------------------------------*/
    /** Returns if the initialization phase was already traced. */
    const bool& isTracedInit() const;

    /** Returns the tape identifier corresponding to the tape
     * containing the initialization at all inner partition grid nodes. */
    short int tapeIdInitPartitionInner();

    /** Returns the tape identifier corresponding to the tape
     * containing the initialization at all local partition grid nodes. */
    short int tapeIdInitPartitionBorder();

    /** Returns the tape identifier corresponding to the tape of the
     * containing the initialization at all partition grid nodes. */
    short int tapeIdInitPartitionAll();

    /** Returns if the update phase was already traced. */
    const bool& isTracedUpdatePhase() const;

    /** Returns the tape identifier corresponding to the tape
     * containing the update step at all inner partition grid nodes. */
    const short int& tapeIdUpdatePartitionInner() const;

    /** Returns the tape identifier corresponding to the tape
     * containing the update step at all border partition grid nodes. */
    const short int& tapeIdUpdatePartitionBorder() const;

    /** Returns the tape identifier corresponding to the tape
     * containing the update step at all all partition grid nodes. */
    const short int& tapeIdUpdatePartitionAll() const;

    /** Returns if the update2 phase was already traced. */
    const bool& isTracedUpdate2Phase() const;

    /** Returns the tape identifier corresponding to the tape
     * containing the update2 step at all inner partition grid nodes. */
    const short int& tapeIdUpdate2PartitionInner() const;

    /** Returns the tape identifier corresponding to the tape
     * containing the update2 step at all border partition grid nodes. */
    const short int& tapeIdUpdate2PartitionBorder() const;

    /** Returns the tape identifier corresponding to the tape
     * containing the update2 step at all all partition grid nodes. */
    const short int& tapeIdUpdate2PartitionAll() const;

    /** Returns the rank of underlying MPI process. */
    int myRank() const;

    /** Returns the total number of nodes of the grid partition. */
    const int& nNodesTotalPartition() const;

    /** Returns the flags if a data field is a known one or not. */
    const std::vector<bool>& isKnownDf() const;

    /** Returns the flags where a data field is defined. */
    const std::vector<int>& nLayers() const;

    /** Returns the default values of all data fields. */
    const std::vector<CT>& defaultValue() const;

    /** Returns the number of elements of all unknown data fields. */
    const int& nElemsUnknownDfs() const;

    /** Returns the number of elements of all known data fields. */
    const int& nElemsKnownDfs() const;

    /** Returns the number of those data fields which are unknown and
     *  defined at the boundary. */
    const int& nUnknownBdryDfs() const;

    /** Returns the vector of all unknown data fields defined at domain. */
    const ScaFES_VectDf<CT>& vectUnknownDfsDomNew() const;

    /** Returns the old vector of all data fields defined at domain. */
    const ScaFES_VectDf<CT>& vectUnknownDfsDomOld() const;

    /** Returns the new vector of all data fields defined at boundary. */
    const ScaFES_VectDf<CT>& vectUnknownDfsBdryNew() const;

    /** Returns the old vector of all data fields defined at boundary. */
    const ScaFES_VectDf<CT>& vectUnknownDfsBdryOld() const;

    /** Returns the number of directions during tracing. */
    const int& nColumns() const;

    /** Returns if gradients should be computed during the
     * initialization and the update phase. */
    const bool& computeGradients() const;

    /** Returns if asynchron MPI communication will be used. */
    const bool& useAsynchronMode() const;

    /** Returns if leap frog scheme will be used. */
    const bool& useLeapfrogIntegration() const;

    /** Returns if the parameter based ansatz will be used. */
    const bool& isParamsBased() const;

    /*--------------------------------------------------------------------------
    | WORK METHODS.
    --------------------------------------------------------------------------*/
    /** Creates new vector of unknown data fields.
     * Splits up allocated memory and distributes it to all data fields. */
    void createVectorOfUnknownDataFields(
        ScaFES_VectDf<CT>& dfInitDomain, ScaFES_VectDf<CT>& dfInitBdry,
        const std::vector<CT>& memoryOfDfs,
        const std::vector<std::string>& nameDatafield,
        const std::string& nameSuffix, const std::vector<bool>& isKnownDf,
        const std::vector<int>& nLayers, const std::vector<int>& stencilWidth,
        const std::vector<ScaFES::Grid<DIM>>& memAll,
        const std::vector<ScaFES::GridSub<DIM>>& memHolePart,
        const std::vector<int>& nRowsDf,
        const std::vector<ScaFES::WriteHowOften>& writeToFile,
        const std::vector<CT>& defaultValue, const int& nColumns);

    /** Copies structure of data field vector. */
    template <typename TT>
    void copyDataFieldStructure(ScaFES_VectDf<TT>& aDfDep,
                                const ScaFES_VectDf<CT>& dfDep,
                                const std::vector<TT>& memoryAdfDep);

    /*------------------------------------------------------------------------*/
    /** Initialise all unknown data fields. */
    void evalInitPbPhase();

    /** Update values at all grid partition nodes of all unknown
     * data fields (1st cycle). */
    void evalUpdatePhase(const int& timeIter);

    /** Update values at all grid partition nodes of all data fields
     *  (2nd cycle) using the leapfrog scheme. */
    void evalUpdate2Phase(const int& timeIter);

    /** Reset values at all grid partition nodes of all unknown
     * data fields to its default values given as constructor parameter. */
    void resetValsAtNodesOfUnknownDfs(const int& timeIter);

    /** Reset times of all unknown data fields. */
    void resetTimesOfUnknownDfs();

    /** Update times of all unknown data fields. */
    void updateTimesOfUnknownDfs();

    /** Swap pointers of all unknown data fields. */
    void swapPointersOfUnknownDfs();

    /** Swap times of all unknown data fields. */
    void swapTimesOfUnknownDfs();

    /** Evaluate known data fields. */
    void evalKnownDfs(const int& timeIter);

    /** Compute errors. */
    void compErrOfUnknownDfs();

    /** Write data fields to file. */
    void writeAllDfs(const int& timeIter);

    /** Determines which MPI communication mode will be used. */
    void setCommunicationMode();

    /** Determines if gradients should be computed during
     * initialization and update phase. */
    void setCompGradients();

    /*------------------------------------------------------------------------*/
    /** Calls barrier for MPI-communicator.
     * Returns, if all processes have entered it. */
    void worldBarrier();

    /** Sets up the kind field w.r.t. the global grid
     * (global inner node, global border node). */
    void setGridGlobalFlags();

    /** Sets up the kind field w.r.t. the underlying grid partition:
     * (local inner node, local border node). */
    void setGridPartitionFlags(const int& maxStencilWidth);

    /** Determines the number of grid partition nodes which are
     * inner local AND inner global. */
    int getNnodesLocInnerGlobInner() const;

    /** Determines the number of grid partition nodes which are
     * inner local AND border global. */
    int getNnodesLocInnerGlobBorder() const;

    /** Determines the number of grid partition nodes which are
     * border local AND inner global. */
    int getNnodesLocBorderGlobInner() const;

    /** Determines the number of grid partition nodes which are
     * border local AND border global. */
    int getNnodesLocBorderGlobBorder() const;

    /** Determines the number of grid partition nodes which are
     * inner local. */
    int getNnodesLocInner() const;

    /** Determines the number of grid partition nodes which are
     * border local. */
    int getNnodesLocBorder() const;

    /*--------------------------------------------------------------------*/
    /** Evaluates a given function at all global inner nodes of the
     * grid partition. */
    template <typename TT>
    void evalGlobInner(ScaFES_VectDf<TT>& dfDep, const int& timeIter);

    /** Evaluates a given function at all global border nodes of the
     * grid partition. */
    template <typename TT>
    void evalGlobBorder(ScaFES_VectDf<TT>& dfDep, const int& timeIter);

    /*------------------------------------------------------------------------*/
    /** Writes all data fields to one NetCDF file depending
     * on the parameter of the constructor. */
    void writeDfsToFile(const int& timeIter);

    /*------------------------------------------------------------------------*/
    /*------------------------------------------------------------------------*/
    /** Evaluates phase type 1 using asynchronous MPI communication. */
    void evalTypeOneAsynchronously(
        const std::string& whichPhase, ScaFES_VectDf<CT>& dfDep,
        ScaFES_VectDf<CT>& dfGradDep, const std::vector<CT>& dfIndep,
        const std::vector<CT>& dfGradIndep, const int& timeIter,
        const int& tapeIdLocalInner, const int& tapeIdLocalBorder,
        ScaFES_FuncTypeOneToImplement<CT> funcPhaseInner,
        ScaFES_FuncTypeOneToImplement<CT> funcPhaseBorder);

    /** Evaluates phase type 2 using asynchronous MPI communication. */
    void evalTypeTwoAsynchronously(
        const std::string& whichPhase, ScaFES_VectDf<CT>& dfDep,
        ScaFES_VectDf<CT>& dfGradDep, const ScaFES_VectDf<CT>& dfIndep,
        const ScaFES_VectDf<CT>& dfGradIndep, const int& timeIter,
        const int& tapeIdLocalInner, const int& tapeIdLocalBorder,
        ScaFES_FuncTypeTwoToImplement<CT> funcPhaseInner,
        ScaFES_FuncTypeTwoToImplement<CT> funcPhaseBorder);

    /** Computes target functional using asynchronous MPI communication. */
    void evalTypeThreeAsynchronously(
        const std::string& whichPhase, CT& dfDep, std::vector<CT>& dfGradDep,
        const ScaFES_VectDf<CT>& dfIndep, const ScaFES_VectDf<CT>& dfGradIndep,
        const int& timeIter, const int& tapeIdLocalInner,
        const int& tapeIdLocalBorder,
        ScaFES_FuncTypeThreeToImplement<CT> funcPhaseInner,
        ScaFES_FuncTypeThreeToImplement<CT> funcPhaseBorder);

    /** Evaluates type 1 using synchronous MPI communication. */
    void evalTypeOneSynchronously(
        const std::string& whichPhase, ScaFES_VectDf<CT>& dfDepDomain,
        ScaFES_VectDf<CT>& dfDepGradDomain, ScaFES_VectDf<CT>& dfDepBdry,
        ScaFES_VectDf<CT>& dfDepGradBdry, const std::vector<CT>& dfIndep,
        const std::vector<CT>& dfGradIndep, const int& timeIter,
        const int& tapeIdLocalAll,
        ScaFES_FuncTypeOneToImplement<CT> funcPhaseInner,
        ScaFES_FuncTypeOneToImplement<CT> funcPhaseBorder);

    /** Evaluates type 2 using synchronous MPI communication. */
    void evalTypeTwoSynchronously(
        const std::string& whichPhase, ScaFES_VectDf<CT>& dfDepDomain,
        ScaFES_VectDf<CT>& dfDepGradDomain, ScaFES_VectDf<CT>& dfDepBdry,
        ScaFES_VectDf<CT>& dfDepGradBdry,
        const ScaFES_VectDf<CT>& dfIndepDomain,
        const ScaFES_VectDf<CT>& dfGradIndepDomain,
        const ScaFES_VectDf<CT>& dfIndepBdry,
        const ScaFES_VectDf<CT>& dfGradIndepBdry, const int& timeIter,
        const int& tapeIdLocalAll,
        ScaFES_FuncTypeTwoToImplement<CT> funcPhaseInner,
        ScaFES_FuncTypeTwoToImplement<CT> funcPhaseBorder);

    /** Computes type 3 using synchronous MPI communication. */
    void evalTypeThreeSynchronously(
        const std::string& whichPhase, CT& dfDep, std::vector<CT>& dfGradDep,
        const ScaFES_VectDf<CT>& dfIndepDomain,
        const ScaFES_VectDf<CT>& dfGradIndepDomain,
        const ScaFES_VectDf<CT>& dfIndepBdry,
        const ScaFES_VectDf<CT>& dfGradIndepBdry, const int& timeIter,
        const int& tapeIdLocalAll,
        ScaFES_FuncTypeThreeToImplement<CT> funcPhaseInner,
        ScaFES_FuncTypeThreeToImplement<CT> funcPhaseBorder);

    /** Evaluates the methods of the given phase
     * at all inner grid partition nodes at a given time step. */
    template <typename TT>
    void evalTypeOneFuncAtPartInner(
        const std::string& whichPhase, ScaFES_VectDf<TT>& dfDep,
        const std::vector<TT>& dfIndep, const int& timeIter,
        ScaFES_FuncTypeOneToImplement<TT> funcPhaseInner,
        ScaFES_FuncTypeOneToImplement<TT> funcPhaseBorder);

    /** Evaluates the methods of the given phase
     * at all inner grid partition nodes at a given time step. */
    template <typename TT>
    void evalTypeTwoFuncAtPartInner(
        const std::string& whichPhase, ScaFES_VectDf<TT>& dfDep,
        const ScaFES_VectDf<TT>& dfIndep, const int& timeIter,
        ScaFES_FuncTypeTwoToImplement<TT> funcPhaseInner,
        ScaFES_FuncTypeTwoToImplement<TT> funcPhaseBorder);

    /** Evaluates the methods of the given phase
     * at all inner grid partition nodes at a given time step. */
    template <typename TT>
    void evalTypeThreeFuncAtPartInner(
        const std::string& whichPhase, TT& dfDep,
        const ScaFES_VectDf<TT>& dfIndep, const int& timeIter,
        ScaFES_FuncTypeThreeToImplement<TT> funcPhaseInner,
        ScaFES_FuncTypeThreeToImplement<TT> funcPhaseBorder);

    /** Evaluates the methods of the given phase
     * at all border grid partition nodes at a given time step. */
    template <typename TT>
    void evalTypeOneFuncAtPartBorder(
        const std::string& whichPhase, ScaFES_VectDf<TT>& dfDep,
        const std::vector<TT>& dfIndep, const int& timeIter,
        ScaFES_FuncTypeOneToImplement<TT> funcPhaseInner,
        ScaFES_FuncTypeOneToImplement<TT> funcPhaseBorder);

    /** Evaluates the methods of the given phase
     * at all border grid partition nodes at a given time step. */
    template <typename TT>
    void evalTypeTwoFuncAtPartBorder(
        const std::string& whichPhase, ScaFES_VectDf<TT>& dfDep,
        const ScaFES_VectDf<TT>& dfIndep, const int& timeIter,
        ScaFES_FuncTypeTwoToImplement<TT> funcPhaseInner,
        ScaFES_FuncTypeTwoToImplement<TT> funcPhaseBorder);

    /** Evaluates the methods of the given phase
     * at all border grid partition nodes at a given time step. */
    template <typename TT>
    void evalTypeThreeFuncAtPartBorder(
        const std::string& whichPhase, TT& dfDep,
        const ScaFES_VectDf<TT>& dfIndep, const int& timeIter,
        ScaFES_FuncTypeThreeToImplement<TT> funcPhaseInner,
        ScaFES_FuncTypeThreeToImplement<TT> funcPhaseBorder);

    /** Evaluates the methods of the given phase
     * at all inner grid global nodes at a given time step. */
    template <typename TT>
    void evalTypeOneFuncAtPartAll(
        const std::string& whichPhase, ScaFES_VectDf<TT>& dfDep,
        const std::vector<TT>& dfIndep, const int& timeIter,
        ScaFES_FuncTypeOneToImplement<TT> funcPhaseInner,
        ScaFES_FuncTypeOneToImplement<TT> funcPhaseBorder);

    /** Evaluates the methods of the given phase
     * at all inner grid global nodes at a given time step. */
    template <typename TT>
    void evalTypeTwoFuncAtPartAll(
        const std::string& whichPhase, ScaFES_VectDf<TT>& dfDep,
        const ScaFES_VectDf<TT>& dfIndep, const int& timeIter,
        ScaFES_FuncTypeTwoToImplement<TT> funcPhaseInner,
        ScaFES_FuncTypeTwoToImplement<TT> funcPhaseBorder);

    /** Evaluates the methods of the given phase
     * at all inner grid global nodes at a given time step. */
    template <typename TT>
    void evalTypeThreeFuncAtPartAll(
        const std::string& whichPhase, TT& dfDep,
        const ScaFES_VectDf<TT>& dfIndep, const int& timeIter,
        ScaFES_FuncTypeThreeToImplement<TT> funcPhaseInner,
        ScaFES_FuncTypeThreeToImplement<TT> funcPhaseBorder);

    /** Evaluates the data fields at all grid partition nodes
     *  at a given time. */
    template <typename TT>
    void evalTypeOneAt(const std::string& whichPhase,
                       const std::string& whereToEval,
                       ScaFES_FuncTypeOneToImplement<TT>,
                       ScaFES_VectDf<TT>& dfDep, const std::vector<TT>& dfIndep,
                       std::vector<int> setGridNodes, const int& timeIter);

    /** Evaluates the data fields at all grid partition nodes
     *  at a given time. */
    template <typename TT>
    void
    evalTypeTwoAt(const std::string& whichPhase, const std::string& whereToEval,
                  ScaFES_FuncTypeTwoToImplement<TT>,
                  ScaFES_FuncTypeTwoToImplement<TT>,
                  ScaFES_VectDf<TT>& dfDep,
                  const ScaFES_VectDf<TT>& dfIndep,
                  std::vector<int> setGridNodes, const int& timeIter);

    /** Evaluates the data fields at all grid partition nodes
     *  at a given time. */
    template <typename TT>
    void evalTypeThreeAt(const std::string& whichPhase,
                         const std::string& whereToEval,
                         ScaFES_FuncTypeThreeToImplement<TT>, TT& dfDep,
                         const ScaFES_VectDf<TT>& dfIndep,
                         std::vector<int> setGridNodes, const int& timeIter);

    /*------------------------------------------------------------------------*/
    /*------------------------------------------------------------------------*/
    /** Traces the init(pb) (parameter based) phase. */
    void traceInitPbPhase(const int& timeIter);

    /** Traces the init(Gb) (grid based) phase. */
    void traceInitGbPhase(const int& timeIter);

    /** Traces the update phase. */
    void traceUpdatePhase(const int& timeIter);

    /** Traces the update2 phase. */
    void traceUpdate2Phase(const int& timeIter);

    /** Traces functions of a given phase containing TypeOne methods. */
    template <typename TT>
    void
    traceTypeOne(const std::string& whichPhase, ScaFES_VectDf<CT>& dfDepDomain,
                 ScaFES_VectDf<CT>& dfDepBdry, const std::vector<CT>& dfIndep,
                 const int& timeIter, const int& tapeIdLocalInner,
                 const int& tapeIdLocalBorder, const int& tapeIdLocalAll,
                 ScaFES_FuncTypeOneToImplement<TT> funcPhaseInner,
                 ScaFES_FuncTypeOneToImplement<TT> funcPhaseBorder);

    /** Traces functions of a given phase containing TypeTwo methods. */
    template <typename TT>
    void traceTypeTwo(const std::string& whichPhase,
                      ScaFES_VectDf<CT>& dfDepDomain,
                      ScaFES_VectDf<CT>& dfDepBdry,
                      const ScaFES_VectDf<CT>& dfIndepDomain,
                      const ScaFES_VectDf<CT>& dfIndepBdry, const int& timeIter,
                      const int& tapeIdLocalInner, const int& tapeIdLocalBorder,
                      const int& tapeIdLocalAll,
                      ScaFES_FuncTypeTwoToImplement<TT> funcPhaseInner,
                      ScaFES_FuncTypeTwoToImplement<TT> funcPhaseBorder);

    /** Traces functions of a given phase containing TypeThree methods. */
    template <typename TT>
    void traceTypeThree(const std::string& whichPhase, CT& dfDep,
                        const ScaFES_VectDf<CT>& dfIndepDomain,
                        const ScaFES_VectDf<CT>& dfIndepBdry,
                        const int& timeIter, const int& tapeIdLocalInner,
                        const int& tapeIdLocalBorder, const int& tapeIdLocalAll,
                        ScaFES_FuncTypeThreeToImplement<TT> funcPhaseInner,
                        ScaFES_FuncTypeThreeToImplement<TT> funcPhaseBorder);

    /*------------------------------------------------------------------------*/
    /** Sets independent variables. */
    template <typename TT> void setIndependentScalar(TT& aDf, const CT& df);

    /** Sets independent variables. */
    template <typename TT>
    void setIndependentDf(ScaFES::DataField<TT, DIM>& aDf,
                          const ScaFES::DataField<CT, DIM>& df);

    /** Sets independent variables. */
    template <typename TT>
    void setIndependentVectDf(ScaFES_VectDf<TT>& aDf,
                              const ScaFES_VectDf<CT>& df);

    /** Sets independent variables. */
    template <typename TT>
    void setIndependentVector(std::vector<TT>& aDf, const std::vector<CT>& df);

    /** Sets dependent variables. */
    template <typename TT> void setDependentScalar(TT& aDf, CT& df);

    /** Sets dependent variables. */
    template <typename TT>
    void setDependentDf(ScaFES::DataField<TT, DIM>& aDf,
                        ScaFES::DataField<CT, DIM>& df);

    /** Sets dependent variables. */
    template <typename TT>
    void setDependentVectDf(ScaFES_VectDf<TT>& aDf, ScaFES_VectDf<CT>& df);

    /** Sets dependent variables. */
    template <typename TT>
    void setDependentVector(std::vector<TT>& aDf, std::vector<CT>& df);

    /** Sets a dependent vector of data fields to an independent vector
     * of data fields at all grid partition nodes at a given time. */
    template <typename TT>
    void setVectDfDepToDfIndep(ScaFES_VectDf<TT>& dfDep,
                               const ScaFES_VectDf<TT>& dfIndep,
                               const int& timeIter);

    /** Sets a dependent vector of data fields to a given default value
     * at all grid partition nodes at a given time. */
    template <typename TT>
    void setVectDfDepToDefValue(ScaFES_VectDf<TT>& dfDep,
                                const std::vector<CT>& defVal,
                                const int& timeIter);

    /** Computes new iterate using ADOL-C tape
     * and the AD forward mode. */
    void evalZosForwardTypeOne(const std::string& whichPhase,
                               ScaFES_VectDf<CT>& dfDep,
                               const short int& tapeId,
                               const std::vector<CT>& dfIndep,
                               const int& timeIter);

    /** Computes new iterate using ADOL-C tape
     * and the AD forward mode. */
    void evalZosForwardTypeTwo(const std::string& whichPhase,
                               ScaFES_VectDf<CT>& dfDep,
                               const short int& tapeId,
                               const ScaFES_VectDf<CT>& dfIndep,
                               const int& currTimeIter);

    /** Evaluates the target function
     *  at all grid partition nodes using the AD forward mode. */
    void evalZosForwardTypeThree(const std::string& whichPhase, CT& dfDep,
                                 const short int& tapeId,
                                 const ScaFES_VectDf<CT>& dfIndep,
                                 const int& timeIter);

    /** Computes new iterate and its gradient using ADOL-C tape
     * and the AD forward mode. */
    void evalFovForwardTypeOne(const std::string& whichPhase,
                               ScaFES_VectDf<CT>& dfDep,
                               ScaFES_VectDf<CT>& dfGradDep,
                               const short int& tapeId,
                               const std::vector<CT>& dfIndep,
                               const std::vector<CT>& dfGradIndep,
                               const int& timeIter);

    /** Computes new iterate and its gradient using ADOL-C tape
     * and the AD forward mode. */
    void evalFovForwardTypeTwo(const std::string& whichPhase,
                               ScaFES_VectDf<CT>& dfDep,
                               ScaFES_VectDf<CT>& dfGradDep,
                               const short int& tapeId,
                               const ScaFES_VectDf<CT>& dfIndep,
                               const ScaFES_VectDf<CT>& dfGradIndep,
                               const int& currTimeIter);

    /** Evaluates the target function and the its gradient
    *  at all grid partition nodes using the AD forward mode. */
    void evalFovForwardTypeThree(const std::string& whichPhase, CT& dfDep,
                                 std::vector<CT>& dfGradDep,
                                 const short int& tapeId,
                                 const ScaFES_VectDf<CT>& dfIndep,
                                 const ScaFES_VectDf<CT>& dfGradIndep,
                                 const int& timeIter);

    /*--------------------------------------------------------------------------
    | MEMBER VARIABLES.
    --------------------------------------------------------------------------*/
    /** Parameters of the problem. */
    ScaFES::Parameters mParams;

    /** Use leapfrog integration scheme. */
    const bool mUseLeapfrogIntegration;

    /** List of interfaces to the neighbours. */
    /** Pointer to related global grid to access partitions and more. */
    ScaFES::GridGlobal<DIM> mGG;

    /** Names of all data fields. */
    std::vector<std::string> mNameDataField;

    /** Stencil widths of all data fields. */
    std::vector<int> mStencilWidth;

    /** Flags for all data fields if a data field is a known one or not. */
    std::vector<bool> mIsKnownDf;

    /** Flags for all dfs where a data field is defined. */
    std::vector<int> mNlayers;

    /** Default values for all dfs. */
    std::vector<CT> mDefaultValue;

    /** Flags for all dfs if a data field should be written or not. */
    std::vector<WriteHowOften> mWriteToFile;

    /** Flags for all data fields if error should be computed or not. */
    std::vector<bool> mComputeError;

    /** Number of given parameters. */
    int mNparams;

    /** Compute gradients? */
    bool mComputeGradients;

    /** Use asynchronous MPI communication? */
    bool mUseAsynchronMode;

    /** (Geometrical) parameters based ansatz? */
    bool mIsParamsBased;

    /** Tape identifier for initialization step. */
    short int mTapeIdInitPartitionInner;

    /** Tape identifier for initialization step. */
    short int mTapeIdInitPartitionBorder;

    /** Tape identifier for initialization step. */
    short int mTapeIdInitPartitionAll;

    /** Are methods \c initInner and \c initBorder traced? */
    bool mIsTracedInit;

    /** Tape identifier for update step at all local inner nodes. */
    short int mTapeIdUpdatePartitionInner;

    /** Tape identifier for update step at all local border nodes. */
    short int mTapeIdUpdatePartitionBorder;

    /** Tape identifier for update step at all global inner nodes. */
    short int mTapeIdUpdatePartitionAll;

    /** Tape identifier for update step 2 at all local inner nodes. */
    short int mTapeIdUpdate2PartitionInner;

    /** Tape identifier for update step 2 at all local border nodes. */
    short int mTapeIdUpdate2PartitionBorder;

    /** Tape identifier for update step 2 at all global inner nodes. */
    short int mTapeIdUpdate2PartitionAll;

    /** Are methods \c updateInner and \c updateBorder traced? */
    bool mIsTracedUpdatePhase;

    /** Are methods \c updateInner2 and \c updateBorder2 traced? */
    bool mIsTracedUpdate2Phase;

    /** Number of total grid nodes of the grid partition. */
    int mNnodesTotalPartition;

    /** Number of elements of all KNOWN data fields. */
    int mNelemsKnownDfs;

    /** Number of elements of all UNKNOWN data fields. */
    int mNelemsUnknownDfs;

    /** Number of UNKNOWN AND BOUNDARY data fields. */
    int mNunknownBdryDfs;

    /** Wall clock time for different sections. */
    std::unordered_map<std::string, double> mWallClockTimes;

    /** Global output file. */
    ScaFES::DataFile<CT, DIM> mCommonFile;

    /** Type of each grid node. */
    ScaFES::DataField<short int, DIM> mKind;

    /** Continuous memory lump for kind field. */
    std::vector<short int> mMemoryKind;

    /** List of all grid partition nodes which are inner grid partition
     * nodes AND which are inner grid global nodes. */
    std::vector<int> mNodesLocInnerGlobInner;

    /** List of all grid partition nodes which are inner grid partition
    * nodes AND which are border grid global nodes. */
    std::vector<int> mNodesLocInnerGlobBorder;

    /** List of all grid partition nodes which are border grid partition
    * nodes AND which are inner grid global nodes. */
    std::vector<int> mNodesLocBorderGlobInner;

    /** List of all grid partition nodes which are border grid partition
     * nodes AND which are border grid global nodes. */
    std::vector<int> mNodesLocBorderGlobBorder;

    /** List of all grid partition nodes which are inner grid partition
     * nodes. */
    std::vector<int> mNodesLocInner;

    /** List of all grid partition nodes which are border grid partition
     * nodes. */
    std::vector<int> mNodesLocBorder;

    /** List of all grid partition nodes. */
    std::vector<int> mNodes;

    /** Vector of all known data fields defined at the whole domain. */
    ScaFES_VectDf<CT> mVectKnownDfDom;

    /** Vector of all known data fields defined at the boundary, only. */
    ScaFES_VectDf<CT> mVectKnownDfBdry;

    /** Continuous memory lump for all known data fields. */
    std::vector<CT> mMemoryVectKnownDf;

    /** Numerical solution at last time step defined at the domain. */
    ScaFES_VectDf<CT> mVectUnknownDfsDomOld;

    /** Numerical solution at last time step defined at the boundary. */
    ScaFES_VectDf<CT> mVectUnknownDfsBdryOld;

    /** Numerical solution at current time step defined at the domain. */
    ScaFES_VectDf<CT> mVectUnknownDfsDomNew;

    /** Numerical solution at current time step defined at the boundary. */
    ScaFES_VectDf<CT> mVectUnknownDfsBdryNew;

    /** Continuous memory lump #1 for storing values of all unknown
     *  data fields at all grid partition nodes (incl. halos)
     *  (old / new scheme). */
    std::vector<CT> mMemoryVectUnknownDfsOne;

    /** Continuous memory lump #2 for storing values of all unknown
     *  data fields at all grid partition nodes (incl. halos)
     *  (old / new scheme). */
    std::vector<CT> mMemoryVectUnknownDfsTwo;

    /** Continuous memory lump #3 for storing values of all unknown
     *  data fields at all grid partition nodes (incl. halos)
     *  This memory lump will be allocated only if a temp. vector
     *  of data fields is needed (case: async. MPI comm. + ADOL-C). */
    std::vector<CT> mMemoryVectUnknownDfsThree;

    /** Gradient of numerical solution at last time step
     *  defined at the domain. */
    ScaFES_VectDf<CT> mVectGradUnknownDfsDomOld;

    /** Gradient of numerical solution at last time step. */
    ScaFES_VectDf<CT> mVectGradUnknownDfsBdryOld;

    /** Gradient of numerical solution at current time step. */
    ScaFES_VectDf<CT> mVectGradUnknownDfsDomNew;

    /** Gradient of numerical solution at current time step. */
    ScaFES_VectDf<CT> mVectGradUnknownDfsBdryNew;

    /** Continuous memory lump #1 for storing values of gradients
     *  of all unknown data fields at all grid partition nodes (incl. halos)
     *  (old / new scheme). */
    std::vector<CT> mMemoryVectGradUnknownDfsOne;

    /** Continuous memory lump #2 for storing values of gradients
     *  of all unknown data fields at all grid partition nodes (incl. halos)
     *  (old / new scheme). */
    std::vector<CT> mMemoryVectGradUnknownDfsTwo;

    /** Continuous memory lump #2 for storing values of gradients
     *  of all unknown data fields at all grid partition nodes (incl. halos)
     *  \remark This memory lump will be allocated only if a temp. vector
     *  of data fields is needed (case: async. MPI comm. + ADOL-C). */
    std::vector<CT> mMemoryVectGradUnknownDfsThree;

    /** Continuous memory lump #1 for storing pointers to rows of gradients
     *  of all unknown data fields at all grid partition nodes (incl. halos)
     *  \remark This memory lump will be allocated only if gradients should
     *  be computed. ADOL-C needs variables of type double** for independent
     * and dependent gradients. */
    std::vector<CT*> mMemoryVectGradUnknownDfsRowsOne;

    /** Continuous memory lump #2 for storing pointers to rows of gradients
     *  of all unknown data fields at all grid partition nodes (incl. halos)
     *   \remark This memory lump will be allocated only if gradients should
     *  be computed. ADOL-C needs variables of type double** for independent
     * and dependent gradients. */
    std::vector<CT*> mMemoryVectGradUnknownDfsRowsTwo;

    /** Temporary vector defined at the domain
     * in order to implement asynchronous communication with ADOL-C. */
    ScaFES_VectDf<CT> mVectUnknownDfsDomTmp;

    /** Temporary vector defined at the boundary
     * in order to implement asynchronous communication with ADOL-C. */
    ScaFES_VectDf<CT> mVectUnknownDfsBdryTmp;

    /** Temporary gradient in order to implement asynchronous communication
     * when using ADOL-C tapes. */
    ScaFES_VectDf<CT> mVectGradUnknownDfsDomTmp;

    /** Temporary gradient in order to implement asynchronous communication
     * when using ADOL-C tapes. */
    ScaFES_VectDf<CT> mVectGradUnknownDfsBdryTmp;

    /** Vector containing current set of (geom.) parameters. */
    std::vector<CT> mVectParamsCurr;

    /** Vector containing gradient of current set of (geom.) parameters. */
    std::vector<CT> mVectGradParamsCurr;

    /** Pointer to memory for storing values of corresponding gradient. */
    std::vector<CT*> mMemoryVectGradParamsAsMatrix;

}; // End of class. //

/*******************************************************************************
 * LIFE CYCLE METHODS.
 ******************************************************************************/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline Problem<OWNPRBLM, CT, DIM>::Problem(
    const ScaFES::Parameters & params, const ScaFES::GridGlobal<DIM>& gg,
    const bool& useLeapfrog, const std::vector<std::string>& nameDatafield,
    const std::vector<int>& stencilWidth, const std::vector<bool>& isKnownDf,
    const std::vector<int>& nLayers, const std::vector<CT>& defaultValue,
    const std::vector<ScaFES::WriteHowOften>& writeToFile,
    const std::vector<bool>& computeError,
    const std::vector<CT>& geomparamsInit)
: mParams(params), mUseLeapfrogIntegration(useLeapfrog), mGG(gg),
  mNameDataField(nameDatafield), mStencilWidth(stencilWidth),
  mIsKnownDf(isKnownDf), mNlayers(nLayers), mDefaultValue(defaultValue),
  mWriteToFile(writeToFile), mComputeError(computeError),
  mNparams(1) // = 1 because gradients are of type R^{m,nParams}.
  ,
  mComputeGradients(false), mUseAsynchronMode(false),
  mIsParamsBased(true) // Grid based ansatz currently not supported!
  ,
  mTapeIdInitPartitionInner(0 + params.myWorld().rank()),
  mTapeIdInitPartitionBorder(1592 + params.myWorld().rank()),
  mTapeIdInitPartitionAll(3184 + params.myWorld().rank()), mIsTracedInit(false),
  mTapeIdUpdatePartitionInner(10922 + params.myWorld().rank()),
  mTapeIdUpdatePartitionBorder(12922 + params.myWorld().rank()),
  mTapeIdUpdatePartitionAll(14922 + params.myWorld().rank()),
  mTapeIdUpdate2PartitionInner(16922 + params.myWorld().rank()),
  mTapeIdUpdate2PartitionBorder(18922 + params.myWorld().rank()),
  mTapeIdUpdate2PartitionAll(20922 + params.myWorld().rank()),
  mIsTracedUpdatePhase(false), mIsTracedUpdate2Phase(false),
  mNnodesTotalPartition(0), mNelemsKnownDfs(0), mNelemsUnknownDfs(0),
  mNunknownBdryDfs(0), mWallClockTimes(), mCommonFile(), mKind(&mParams), mMemoryKind(),
  mNodesLocInnerGlobInner(), mNodesLocInnerGlobBorder(),
  mNodesLocBorderGlobInner(), mNodesLocBorderGlobBorder(),
  mNodesLocInner(), mNodesLocBorder(), mNodes(),
  mVectKnownDfDom(),
  mVectKnownDfBdry(), mMemoryVectKnownDf(), mVectUnknownDfsDomOld(),
  mVectUnknownDfsBdryOld(), mVectUnknownDfsDomNew(), mVectUnknownDfsBdryNew(),
  mMemoryVectUnknownDfsOne(), mMemoryVectUnknownDfsTwo(),
  mMemoryVectUnknownDfsThree(), mVectGradUnknownDfsDomOld(),
  mVectGradUnknownDfsBdryOld(), mVectGradUnknownDfsDomNew(),
  mVectGradUnknownDfsBdryNew(), mMemoryVectGradUnknownDfsOne(),
  mMemoryVectGradUnknownDfsTwo(), mMemoryVectGradUnknownDfsThree(),
  mMemoryVectGradUnknownDfsRowsOne(), mMemoryVectGradUnknownDfsRowsTwo(),
  mVectUnknownDfsDomTmp(), mVectUnknownDfsBdryTmp(),
  mVectGradUnknownDfsDomTmp(), mVectGradUnknownDfsBdryTmp(),
  mVectParamsCurr(geomparamsInit), mVectGradParamsCurr(),
  mMemoryVectGradParamsAsMatrix()
{
#ifdef VTRACE
    VT_TRACER("ScaFES::Problem Ctor");
#endif
    this->mWallClockTimes["trace"] = 0.0;
    this->mWallClockTimes["write"] = 0.0;
    this->mWallClockTimes["evalKnownDf"] = 0.0;
    this->mWallClockTimes["sync"] = 0.0;
    this->mWallClockTimes["update"] = 0.0;
    this->mWallClockTimes["init"] = 0.0;

    /*------------------------------------------------------------------------*/
    this->setCommunicationMode();
    this->setCompGradients();

    if (this->params().asynchronMode() && this->params().enabledAdolc() &&
        this->computeGradients())
    { // AFTER calling setCompGradients().
        throw std::runtime_error(
            "Async. MPI + ADOL-C + compGradients currently not supported!");
    }

    /*------------------------------------------------------------------------*/
    // Check sizes of parameters (basis: size of nameDatafield).
    // It suffices to check if one parameter has another size than the others.
    if (nameDatafield.size() != stencilWidth.size())
    {
        throw std::runtime_error("Mismatch in sizes of parameters!");
    }
    if (nameDatafield.size() != isKnownDf.size())
    {
        throw std::runtime_error("Mismatch in sizes of parameters!");
    }
    if (nameDatafield.size() != nLayers.size())
    {
        throw std::runtime_error("Mismatch in sizes of parameters!");
    }
    if (nameDatafield.size() != defaultValue.size())
    {
        throw std::runtime_error("Mismatch in sizes of parameters!");
    }
    if (nameDatafield.size() != writeToFile.size())
    {
        throw std::runtime_error("Mismatch in sizes of parameters!");
    }
    if (nameDatafield.size() != computeError.size())
    {
        throw std::runtime_error("Mismatch in sizes of parameters!");
    }

    /*------------------------------------------------------------------------*/
    // Determine number of all data fields.
    int nDataFields = stencilWidth.size();

    /*------------------------------------------------------------------------*/
    for (int ii = 0; ii < nDataFields; ++ii)
    {
        for (std::size_t jj = 0; jj < DIM; ++jj)
        {
            if (nLayers.at(ii) > gg.discreteDomain().nNodes()[jj] / 2)
            {
                throw std::runtime_error(
                    "Number of layers does not fit to global grid!");
            }
        }
    }

    /*------------------------------------------------------------------------*/
    if (0 < geomparamsInit.size())
    {
        this->mIsParamsBased = true;
        this->mVectGradParamsCurr.resize(
            geomparamsInit.size() * geomparamsInit.size(), static_cast<CT>(0));
        this->mMemoryVectGradParamsAsMatrix.resize(
            geomparamsInit.size() * geomparamsInit.size(), static_cast<CT*>(0));
    }

    /*------------------------------------------------------------------------*/
    if (this->isParamsBased())
    {
        this->mNparams = 1; // Default = 1.
        if (0 < this->mVectParamsCurr.size())
        {
            this->mNparams = this->mVectParamsCurr.size();
        }
    }

    /*------------------------------------------------------------------------*/
    // Determine number of total nodes of the grid partition.
    this->mNnodesTotalPartition =
        gg.partition(this->myRank()).nNodesSub().size();

    /*------------------------------------------------------------------------*/
    // Determine number of unknown data fields.
    int nUnknownDfs = 0;
    for (int ii = 0; ii < nDataFields; ++ii)
    {
        if (!this->isKnownDf().at(ii))
        {
            ++nUnknownDfs;
        }
    }

    /*------------------------------------------------------------------------*/
    // Determine number of unknown data fields which are defined at boundary.
    this->mNunknownBdryDfs = 0;
    for (int ii = 0; ii < nDataFields; ++ii)
    {
        if (!(this->isKnownDf().at(ii)))
        {
            if (0 < nLayers.at(ii))
            {
                ++this->mNunknownBdryDfs;
            }
        }
    }

    /*------------------------------------------------------------------------*/
    // Determine memAll and memNormal wrt. to the grid partition
    // for all data fields.
    std::vector<ScaFES::GridSub<DIM>> memNormal(nDataFields);
    std::vector<ScaFES::Grid<DIM>> memAll(nDataFields);
    for (int ii = 0; ii < nDataFields; ++ii)
    {
        ScaFES_IntNtuple idxNodeFirst =
            gg.partition(this->myRank()).idxNodeFirstSub() -
            ScaFES_IntNtuple(stencilWidth.at(ii));
        ScaFES_IntNtuple idxNodeLast =
            gg.partition(this->myRank()).idxNodeLastSub() +
            ScaFES_IntNtuple(stencilWidth.at(ii));
        ScaFES_IntNtuple idxNodeFirstNormal =
            gg.partition(this->myRank()).idxNodeFirstSub();
        ScaFES_IntNtuple idxNodeLastNormal =
            gg.partition(this->myRank()).idxNodeLastSub();
        memAll[ii] = ScaFES::Grid<DIM>(idxNodeFirst, idxNodeLast,
                                       gg.discreteDomain().coordinates(idxNodeFirst),
                                       gg.discreteDomain().coordinates(idxNodeLast));
        memNormal[ii] = ScaFES::GridSub<DIM>(idxNodeFirstNormal,
                                             idxNodeLastNormal, memAll[ii]);
    }

    /*------------------------------------------------------------------------*/
    // Determine memHole wrt. to global grid for all data fields.
    std::vector<ScaFES::GridSub<DIM>> memHoleGlobal(nDataFields);
    for (int ii = 0; ii < nDataFields; ++ii)
    {
        ScaFES_IntNtuple idxNodeFirst =
            gg.partition(this->myRank()).idxNodeFirstSub() -
            ScaFES_IntNtuple(stencilWidth.at(ii));
        ScaFES_IntNtuple idxNodeLast =
            gg.partition(this->myRank()).idxNodeLastSub() +
            ScaFES_IntNtuple(stencilWidth.at(ii));
        memHoleGlobal[ii] = ScaFES::GridSub<DIM>(
            gg.discreteDomain().idxNodeFirst() + ScaFES_IntNtuple(nLayers.at(ii)),
            gg.discreteDomain().idxNodeLast() - ScaFES_IntNtuple(nLayers.at(ii)),
            gg.discreteDomain());
    }

    /*------------------------------------------------------------------------*/
    // Determine number of elements of all data fields.
    std::vector<ScaFES::GridSub<DIM>> memHolePart(nDataFields);
    for (int ii = 0; ii < nDataFields; ++ii)
    {
        memHolePart[ii] = ScaFES::GridSub<DIM>(
            ScaFES::max(memNormal.at(ii).idxNodeFirstSub(),
                        memHoleGlobal.at(ii).idxNodeFirstSub()),
            ScaFES::min(memNormal.at(ii).idxNodeLastSub(),
                        memHoleGlobal.at(ii).idxNodeLastSub()),
            memAll.at(ii));
    }
    /*------------------------------------------------------------------------*/
    // Determine number of elements of all data fields (grid nodes incl. time).
    std::vector<int> nRowsDf(nDataFields);
    for (int ii = 0; ii < nDataFields; ++ii)
    {
        if (0 < nLayers.at(ii))
        {
            nRowsDf[ii] = memAll[ii].nNodes().size() -
                          memHolePart[ii].nNodesSub().size() + 1;
        }
        else
        {
            nRowsDf[ii] = memAll[ii].nNodes().size() + 1;
        }
    }

    /*------------------------------------------------------------------------*/
    // Determine number of elements of all unknown data fields.
    this->mNelemsUnknownDfs = 0;
    for (int ii = 0; ii < nDataFields; ++ii)
    {
        if (!(this->isKnownDf().at(ii)))
        {
            this->mNelemsUnknownDfs += nRowsDf[ii];
        }
    }
    if (0 == this->nElemsUnknownDfs())
    {
        throw std::runtime_error("#(unknown data fields) = 0.");
    }
    /*------------------------------------------------------------------------*/
    // Determine number of elements of all known data fields.
    this->mNelemsKnownDfs = 0;
    for (int ii = 0; ii < nDataFields; ++ii)
    {
        if (this->isKnownDf().at(ii))
        {
            this->mNelemsKnownDfs += nRowsDf[ii];
        }
    }

    /*------------------------------------------------------------------------*/
    // Check if reference vector for all UNKNOWN data fields should be created.
    bool createVectKnownDf = false;
    for (int ii = 0; ii < nDataFields; ++ii)
    {
        if (isKnownDf.at(ii) || computeError.at(ii))
        {
            createVectKnownDf = true;
            break;
        }
    }

    /*------------------------------------------------------------------------*/
    // Create known data fields only if solution is known and some error
    // should be computed.
    // Process all KNOWN data fields (including exact solutions of known dfs).
    // Create one large memory block for all KNOWN data fields.
    const int N_COLUMNS_KNOWN_DF = 1;
    int nElemsAllDfs = this->nElemsKnownDfs();
    if (createVectKnownDf)
    {
        nElemsAllDfs += this->nElemsUnknownDfs();
    }
    if (0 < nElemsAllDfs)
    {
        this->mMemoryVectKnownDf.resize(N_COLUMNS_KNOWN_DF * nElemsAllDfs,
                                        static_cast<CT>(0));
        this->mVectKnownDfDom.reserve(nDataFields);
        this->mVectKnownDfBdry.reserve(nDataFields);
    }
    CT* elemDataKnownDfCurr = this->mMemoryVectKnownDf.data();
    /*------------------------------------------------------------------------*/
    // Firstly, all data fields which are defined at the whole domain.
    // Secondly, all data fields which are defined at the boundary, only.
    // This reordering of the data fields is important for the ADOL-C tracing
    // and evaluation as the fields must be available as a continuous memory
    // lump.
    for (int ii = 0; ii < nDataFields; ++ii)
    {
        if (0 == nLayers.at(ii))
        {
            if (this->isKnownDf().at(ii))
            {
                this->mVectKnownDfDom.push_back(ScaFES::DataField<CT, DIM>(
                    nameDatafield.at(ii), &(this->mParams), this->globalGrid(),
                    memAll.at(ii), stencilWidth.at(ii), elemDataKnownDfCurr,
                    N_COLUMNS_KNOWN_DF, memHolePart.at(ii), nLayers.at(ii),
                    writeToFile.at(ii)));
                elemDataKnownDfCurr += nRowsDf.at(ii);
            }
            else
            {
                // Create reference data fields regardless if error should be
                // computed or not. This important for the correct order of
                // the data fields.
                // BUT: Memory for reference data fields will be created only
                // if error should be computed.
                std::string nameDf = nameDatafield.at(ii) + "SolDom";
                this->mVectKnownDfDom.push_back(ScaFES::DataField<CT, DIM>(
                    nameDf, &(this->mParams), this->globalGrid(), memAll.at(ii),
                    0, elemDataKnownDfCurr, N_COLUMNS_KNOWN_DF,
                    memHolePart.at(ii), nLayers.at(ii), writeToFile.at(ii)));
                if (computeError.at(ii))
                {
                    elemDataKnownDfCurr += nRowsDf.at(ii);
                }
            }
        }
    }
    /*------------------------------------------------------------------------*/
    for (int ii = 0; ii < nDataFields; ++ii)
    {
        if (0 < nLayers.at(ii))
        {
            if (this->isKnownDf().at(ii))
            {
                this->mVectKnownDfBdry.push_back(ScaFES::DataField<CT, DIM>(
                    nameDatafield.at(ii), &(this->mParams), this->globalGrid(),
                    memAll.at(ii), stencilWidth.at(ii), elemDataKnownDfCurr,
                    N_COLUMNS_KNOWN_DF, memHolePart.at(ii), nLayers.at(ii),
                    writeToFile.at(ii)));
                elemDataKnownDfCurr += nRowsDf.at(ii);
            }
            else
            {
                // Create reference data field regardless if error should be
                // computed or not. This important for the correct order of
                // the data fields.
                // BUT: Memory for reference data field will be created only in
                // case OF COMPUTING THE ERROR.
                std::string nameDf = nameDatafield.at(ii) + "SolBdry";
                this->mVectKnownDfBdry.push_back(ScaFES::DataField<CT, DIM>(
                    nameDf, &(this->mParams), this->globalGrid(), memAll.at(ii),
                    0, elemDataKnownDfCurr, N_COLUMNS_KNOWN_DF,
                    memHolePart.at(ii), nLayers.at(ii), writeToFile.at(ii)));
                if (computeError.at(ii))
                {
                    elemDataKnownDfCurr += nRowsDf.at(ii);
                }
            }
        }
    }
    /*------------------------------------------------------------------------*/
    // Set default value for each known data field.
    // //     int idDf = 0;
    // //     for (int ii = 0; ii < nDataFields; ++ii) {
    // //         if (0 == nLayers.at(ii)) {
    // //             if (this->isKnownDf().at(ii)) {
    // //                 this->mVectKnownDfDom.at(idDf) = defaultValue.at(ii);
    // //                 ++idDf;
    // //             }
    // //         }
    // //     }
    // //     idDf = 0;
    // //     for (int ii = 0; ii < nDataFields; ++ii) {
    // //         if (0 < nLayers.at(ii)) {
    // //             if (this->isKnownDf().at(ii)) {
    // //                 this->mVectKnownDfBdry.at(idDf) = defaultValue.at(ii);
    // //                 ++idDf;
    // //             }
    // //         }
    // //     }

    /*------------------------------------------------------------------------*/
    /*------------------------------------------------------------------------*/
    // Process all UNKNOWN data fields.
    // Create one large memory block for all UNKNOWN data fields.
    // Use old / new scheme.
    // Firstly, allocate memory for all old data fields (independent ones),
    // then, allocate memory for all new data fields (dependent ones).
    const int N_COLUMNS_VECTUNKNOWN = 1;
    if (0 < this->nElemsUnknownDfs())
    {
        this->mMemoryVectUnknownDfsOne.resize(N_COLUMNS_VECTUNKNOWN *
                                                  this->nElemsUnknownDfs(),
                                              static_cast<CT>(0));
    }
    if (0 < this->nElemsUnknownDfs())
    {
        this->mMemoryVectUnknownDfsTwo.resize(N_COLUMNS_VECTUNKNOWN *
                                                  this->nElemsUnknownDfs(),
                                              static_cast<CT>(0));
    }
    /*------------------------------------------------------------------------*/
    std::string nameSuffixOld = "Old";
    this->createVectorOfUnknownDataFields(
        this->mVectUnknownDfsDomOld, this->mVectUnknownDfsBdryOld,
        this->mMemoryVectUnknownDfsOne, nameDatafield, nameSuffixOld, isKnownDf,
        nLayers, stencilWidth, memAll, memHolePart, nRowsDf, writeToFile,
        defaultValue, N_COLUMNS_VECTUNKNOWN);
    /*------------------------------------------------------------------------*/
    std::string nameSuffixNew = "New";
    this->createVectorOfUnknownDataFields(
        this->mVectUnknownDfsDomNew, this->mVectUnknownDfsBdryNew,
        this->mMemoryVectUnknownDfsTwo, nameDatafield, nameSuffixNew, isKnownDf,
        nLayers, stencilWidth, memAll, memHolePart, nRowsDf, writeToFile,
        defaultValue, N_COLUMNS_VECTUNKNOWN);
    /*------------------------------------------------------------------------*/
    // Create temporary vector only in case of
    // ADOL-C + asynchronously MPI communication.
    if (params.enabledAdolc() && this->useAsynchronMode())
    {
        if (0 < this->nElemsUnknownDfs())
        {
            this->mMemoryVectUnknownDfsThree.resize(
                N_COLUMNS_VECTUNKNOWN * this->nElemsUnknownDfs(),
                static_cast<CT>(0));
        }
        std::string nameSuffixTmp = "Tmp";
        this->createVectorOfUnknownDataFields(
            this->mVectUnknownDfsDomTmp, this->mVectUnknownDfsBdryTmp,
            this->mMemoryVectUnknownDfsThree, nameDatafield, nameSuffixTmp,
            isKnownDf, nLayers, stencilWidth, memAll, memHolePart, nRowsDf,
            writeToFile, defaultValue, N_COLUMNS_VECTUNKNOWN);
    }

    /*------------------------------------------------------------------------*/
    /*------------------------------------------------------------------------*/
    if (this->computeGradients())
    {
        // Process GRADIENTS of all UNKNOWN data fields.
        // Use old / new scheme.
        // Gradients are \in R^{nxm}:
        // Nevertheless, create an one-dimensional array:
        // +++ one continuous memory lump
        // +++ Methods of class DataField can be used,
        //     if pointer to one-dimensional
        //     array will be passed to the constructor of the class DataField.
        const int N_COLUMNS_VECTGRADUNKNOWN = this->nColumns();
        if (0 < this->nElemsUnknownDfs())
        {
            this->mMemoryVectGradUnknownDfsOne.resize(
                N_COLUMNS_VECTGRADUNKNOWN * this->nElemsUnknownDfs(),
                static_cast<CT>(0));
        }
        if (0 < this->nElemsUnknownDfs())
        {
            this->mMemoryVectGradUnknownDfsTwo.resize(
                (N_COLUMNS_VECTGRADUNKNOWN * (this->nElemsUnknownDfs())),
                static_cast<CT>(0));
        }
        /*--------------------------------------------------------------------*/
        std::string nameSuffixGradOld = "GradOld";
        this->createVectorOfUnknownDataFields(
            this->mVectGradUnknownDfsDomOld, this->mVectGradUnknownDfsBdryOld,
            this->mMemoryVectGradUnknownDfsOne, nameDatafield,
            nameSuffixGradOld, isKnownDf, nLayers, stencilWidth, memAll,
            memHolePart, nRowsDf, writeToFile, defaultValue,
            N_COLUMNS_VECTGRADUNKNOWN);
        /*--------------------------------------------------------------------*/
        std::string nameSuffixGradNew = "GradNew";
        this->createVectorOfUnknownDataFields(
            this->mVectGradUnknownDfsDomNew, this->mVectGradUnknownDfsBdryNew,
            this->mMemoryVectGradUnknownDfsTwo, nameDatafield,
            nameSuffixGradNew, isKnownDf, nLayers, stencilWidth, memAll,
            memHolePart, nRowsDf, writeToFile, defaultValue,
            N_COLUMNS_VECTGRADUNKNOWN);
        /*--------------------------------------------------------------------*/
        // Create temporary vector of data fields only in case of
        // ADOL-C + asynchronously MPI communication.
        if (params.enabledAdolc() && this->useAsynchronMode())
        {
            if (0 < this->nElemsUnknownDfs())
            {
                this->mMemoryVectGradUnknownDfsThree.resize(
                    (N_COLUMNS_VECTGRADUNKNOWN * (this->nElemsUnknownDfs())),
                    static_cast<CT>(0));
            }
            std::string nameSuffixGradTmp = "GradTmp";
            this->createVectorOfUnknownDataFields(
                this->mVectGradUnknownDfsDomTmp,
                this->mVectGradUnknownDfsBdryTmp,
                this->mMemoryVectGradUnknownDfsThree, nameDatafield,
                nameSuffixGradTmp, isKnownDf, nLayers, stencilWidth, memAll,
                memHolePart, nRowsDf, writeToFile, defaultValue,
                N_COLUMNS_VECTGRADUNKNOWN);
        }
    }

    /*------------------------------------------------------------------------*/
    /*------------------------------------------------------------------------*/
    // Create memory for pointers only if gradients should be computed
    // using ADOL-C.
    if (params.enabledAdolc() && this->computeGradients())
    {
        if (0 < this->nElemsUnknownDfs())
        {
            this->mMemoryVectGradUnknownDfsRowsOne.resize(
                1 * this->nElemsUnknownDfs(), static_cast<CT*>(0));
            this->mMemoryVectGradUnknownDfsRowsTwo.resize(
                1 * this->nElemsUnknownDfs(), static_cast<CT*>(0));
        }
    }

    /*------------------------------------------------------------------------*/
    /*------------------------------------------------------------------------*/
    int nDataFieldsToWrite = 0;
    for (std::size_t ii = 0; ii < this->mWriteToFile.size(); ++ii)
    {
        if (ScaFES::WriteHowOften::NEVER != this->mWriteToFile[ii])
        {
            ++nDataFieldsToWrite;
        }
    }

    /*------------------------------------------------------------------------*/
    std::vector<std::string> tmpNames;
    tmpNames.reserve(nDataFieldsToWrite);
    int idxKnownDf = 0;
    int idxUnknownDf = 0;
    for (std::size_t ii = 0; ii < this->mWriteToFile.size(); ++ii)
    {
        if (0 == nLayers.at(ii))
        {
            if (this->isKnownDf().at(ii))
            {
                if (ScaFES::WriteHowOften::NEVER != this->mWriteToFile[ii])
                {
                    tmpNames.push_back(
                        this->mVectKnownDfDom[idxKnownDf].name());
                }
                ++idxKnownDf;
            }
            else
            {
                if (ScaFES::WriteHowOften::NEVER != this->mWriteToFile[ii])
                {
                    tmpNames.push_back(
                        this->vectUnknownDfsDomNew()[idxUnknownDf].name());
                }
                ++idxUnknownDf;
            }
        }
    }
    idxKnownDf = 0;
    idxUnknownDf = 0;
    for (std::size_t ii = 0; ii < this->mWriteToFile.size(); ++ii)
    {
        if (0 < nLayers.at(ii))
        {
            if (this->isKnownDf().at(ii))
            {
                if (ScaFES::WriteHowOften::NEVER != this->mWriteToFile[ii])
                {
                    tmpNames.push_back(
                        this->mVectKnownDfBdry[idxKnownDf].name());
                }
                ++idxKnownDf;
            }
            else
            {
                if (ScaFES::WriteHowOften::NEVER != this->mWriteToFile[ii])
                {
                    tmpNames.push_back(
                        this->vectUnknownDfsBdryNew()[idxUnknownDf].name());
                }
                ++idxUnknownDf;
            }
        }
    }

    /*------------------------------------------------------------------------*/
    std::vector<ScaFES::GridSub<DIM>> tmpMemNormal;
    tmpMemNormal.reserve(nDataFieldsToWrite);
    idxKnownDf = 0;
    idxUnknownDf = 0;
    for (std::size_t ii = 0; ii < this->mWriteToFile.size(); ++ii)
    {
        if (0 == nLayers.at(ii))
        {
            if (this->isKnownDf().at(ii))
            {
                if (ScaFES::WriteHowOften::NEVER != this->mWriteToFile[ii])
                {
                    tmpMemNormal.push_back(
                        this->mVectKnownDfDom[idxKnownDf].memNormal());
                }
                ++idxKnownDf;
            }
            else
            {
                if (ScaFES::WriteHowOften::NEVER != this->mWriteToFile[ii])
                {
                    tmpMemNormal.push_back(
                        this->vectUnknownDfsDomNew()[idxUnknownDf].memNormal());
                }
                ++idxUnknownDf;
            }
        }
    }
    idxKnownDf = 0;
    idxUnknownDf = 0;
    for (std::size_t ii = 0; ii < this->mWriteToFile.size(); ++ii)
    {
        if (0 < nLayers.at(ii))
        {
            if (this->isKnownDf().at(ii))
            {
                if (ScaFES::WriteHowOften::NEVER != this->mWriteToFile[ii])
                {
                    tmpMemNormal.push_back(
                        this->mVectKnownDfBdry[idxKnownDf].memNormal());
                }
                ++idxKnownDf;
            }
            else
            {
                if (ScaFES::WriteHowOften::NEVER != this->mWriteToFile[ii])
                {
                    tmpMemNormal.push_back(
                        this->vectUnknownDfsBdryNew()[idxUnknownDf]
                            .memNormal());
                }
                ++idxUnknownDf;
            }
        }
    }

    std::ostringstream tmpStringstream;
    tmpStringstream << this->params().nameDataFile();

    this->mCommonFile = ScaFES::DataFile<CT, DIM>(
        tmpStringstream.str(), tmpNames, this->params().myWorld(),
        this->globalGrid().discreteDomain().nNodes(), tmpMemNormal);

    /*------------------------------------------------------------------------*/
    // Read in partition file.
    if (this->params().readPartitionFile())
    {
        try
        {
#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
            std::ifstream ifs(this->params().nameReadPartitionFile().c_str());
            boost::archive::binary_iarchive ia(ifs);
            ia >> (this->mGG);
#else
            throw std::runtime_error("Boost.serialization not enabled.");
#endif

            // The partition file is okay.
            // Check if number of partitions == number of MPI procs.
            if (this->params().myWorld().size() !=
                this->globalGrid().nPartitionsTotal())
            {
                throw std::runtime_error("Number of partitions != mpiSize().");
            }
        }
        catch (std::exception& e)
        {
            throw std::runtime_error("Could not open partition file.");
        }
    }

    /*------------------------------------------------------------------------*/
    // Serialize the global grid into an archive and write file.
    if (this->params().writePartitionFile())
    {
        std::ostringstream tmpStringstream;
        tmpStringstream << this->params().nameWritePartitionFile()
                        << ".partfile";
        std::string tmpStr = tmpStringstream.str();
        std::ofstream ofs(tmpStr.c_str());
#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
        boost::archive::binary_oarchive oa(ofs);
        try
        {
            oa << (this->mGG);
        }
        catch (...)
        {
            throw std::runtime_error("Global grid could be not be serialized.");
        }
#else
        throw std::runtime_error("Boost.serialization not enabled.");
#endif


        ofs.close();
    }

    /*------------------------------------------------------------------------*/
    ScaFES_IntNtuple idxNodeFirstKindAll =
        gg.partition(this->myRank()).idxNodeFirstSub();
    ScaFES_IntNtuple idxNodeLastKindAll =
        gg.partition(this->myRank()).idxNodeLastSub();
    ScaFES_IntNtuple idxNodeFirstKindNormal =
        gg.partition(this->myRank()).idxNodeFirstSub();
    ScaFES_IntNtuple idxNodeLastKindNormal =
        gg.partition(this->myRank()).idxNodeLastSub();
    ScaFES::Grid<DIM> memAllKind(idxNodeFirstKindAll, idxNodeLastKindAll,
                                 gg.discreteDomain().coordinates(idxNodeFirstKindAll),
                                 gg.discreteDomain().coordinates(idxNodeLastKindAll));
    ScaFES::GridSub<DIM> memNormalKind(idxNodeFirstKindNormal,
                                       idxNodeLastKindNormal, memAllKind);
    ScaFES::GridSub<DIM> memHolePartKind(memNormalKind.idxNodeFirstSub(),
                                         memNormalKind.idxNodeLastSub(),
                                         memAllKind);

    if (0 < this->nNodesTotalPartition())
    {
        this->mMemoryKind.resize(this->nNodesTotalPartition() + 1, // incl. time
                                 static_cast<short int>(0));
    }
    short int* elemDataKind = this->mMemoryKind.data();
    ScaFES::WriteHowOften writeKind = ScaFES::WriteHowOften::NEVER;
    if (this->params().writeKindFile())
    {
        writeKind = ScaFES::WriteHowOften::LIKE_GIVEN_AT_CL;
    }
    this->mKind = ScaFES::DataField<short int, DIM>(
        "kind", &(this->mParams), this->globalGrid(), memNormalKind, 0,
        elemDataKind, 1, memHolePartKind, 0, writeKind);

    /*------------------------------------------------------------------------*/
    // Read in a given kind file?
    if (this->params().readKindFile())
    {
        this->mKind.read(0);
    }
    else
    {
        this->setGridGlobalFlags();
    }
    int maxStencilWidth = 0;
    for (std::size_t ii = 0; ii < this->vectUnknownDfsDomNew().size(); ++ii)
    {
        if (maxStencilWidth < this->vectUnknownDfsDomNew()[ii].stencilWidth())
        {
            maxStencilWidth = this->vectUnknownDfsDomNew()[ii].stencilWidth();
        }
    }
    for (std::size_t ii = 0; ii < this->vectUnknownDfsBdryNew().size(); ++ii)
    {
        if (maxStencilWidth < this->vectUnknownDfsBdryNew()[ii].stencilWidth())
        {
            maxStencilWidth = this->vectUnknownDfsBdryNew()[ii].stencilWidth();
        }
    }
    this->setGridPartitionFlags(maxStencilWidth);

    /*------------------------------------------------------------------------*/
    // Determine node list of all global inner and all global border nodes.
    this->mNodesLocInnerGlobInner.resize(this->getNnodesLocInnerGlobInner());
    int nNodes = 0;
    for (int ii = 0; ii < this->nNodesTotalPartition(); ++ii)
    {
        if (this->mKind(ii) & this->globalGrid().SCAFES_PARTITION_INNER_NODE)
        {
            if (this->mKind(ii) & this->globalGrid().SCAFES_GLOBAL_INNER_NODE)
            {
                this->mNodesLocInnerGlobInner[nNodes] = ii;
                ++nNodes;
            }
        }
    }
    this->mNodesLocInnerGlobBorder.resize(this->getNnodesLocInnerGlobBorder());
    nNodes = 0;
    for (int ii = 0; ii < this->nNodesTotalPartition(); ++ii)
    {
        if (this->mKind(ii) & this->globalGrid().SCAFES_PARTITION_INNER_NODE)
        {
            if (this->mKind(ii) & this->globalGrid().SCAFES_GLOBAL_BORDER_NODE)
            {
                this->mNodesLocInnerGlobBorder[nNodes] = ii;
                ++nNodes;
            }
        }
    }
    this->mNodesLocBorderGlobInner.resize(this->getNnodesLocBorderGlobInner());
    nNodes = 0;
    for (int ii = 0; ii < this->nNodesTotalPartition(); ++ii)
    {
        if (this->mKind(ii) & this->globalGrid().SCAFES_PARTITION_BORDER_NODE)
        {
            if (this->mKind(ii) & this->globalGrid().SCAFES_GLOBAL_INNER_NODE)
            {
                this->mNodesLocBorderGlobInner[nNodes] = ii;
                ++nNodes;
            }
        }
    }
    this->mNodesLocBorderGlobBorder.resize(
        this->getNnodesLocBorderGlobBorder());
    nNodes = 0;
    for (int ii = 0; ii < this->nNodesTotalPartition(); ++ii)
    {
        if (this->mKind(ii) & this->globalGrid().SCAFES_PARTITION_BORDER_NODE)
        {
            if (this->mKind(ii) & this->globalGrid().SCAFES_GLOBAL_BORDER_NODE)
            {
                this->mNodesLocBorderGlobBorder[nNodes] = ii;
                ++nNodes;
            }
        }
    }
    /*------------------------------------------------------------------------*/
    this->mNodesLocInner.resize(this->getNnodesLocInner());
    nNodes = 0;
    for (int ii = 0; ii < this->nNodesTotalPartition(); ++ii)
    {
        if (this->mKind(ii) & this->globalGrid().SCAFES_PARTITION_INNER_NODE)
        {
            this->mNodesLocInner[nNodes] = ii;
            ++nNodes;
        }
    }
    this->mNodesLocBorder.resize(this->getNnodesLocBorder());
    nNodes = 0;
    for (int ii = 0; ii < this->nNodesTotalPartition(); ++ii)
    {
        if (this->mKind(ii) & this->globalGrid().SCAFES_PARTITION_BORDER_NODE)
        {
            this->mNodesLocBorder[nNodes] = ii;
            ++nNodes;
        }
    }
    this->mNodes.resize(this->nNodesTotalPartition());
    for (int ii = 0; ii < this->nNodesTotalPartition(); ++ii)
    {
        this->mNodes[ii] = ii;
    }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::createVectorOfUnknownDataFields(
    ScaFES_VectDf<CT>& dfDomain, ScaFES_VectDf<CT>& dfBdry,
    const std::vector<CT>& memoryOfDfs,
    const std::vector<std::string>& nameDatafield,
    const std::string& nameSuffix, const std::vector<bool>& isKnownDf,
    const std::vector<int>& nLayers, const std::vector<int>& stencilWidth,
    const std::vector<ScaFES::Grid<DIM>>& memAll,
    const std::vector<ScaFES::GridSub<DIM>>& memHolePart,
    const std::vector<int>& nRowsDf,
    const std::vector<ScaFES::WriteHowOften>& writeToFile,
    const std::vector<CT>& defaultValue, const int& nColumns)
{
    CT* elemDataCurr = const_cast<CT*>(memoryOfDfs.data());
    int idDf = 0;

    for (std::size_t ii = 0; ii < nLayers.size(); ++ii)
    {
        if (!(isKnownDf.at(ii)))
        {
            if (0 == nLayers.at(ii))
            {
                std::string nameDf = nameDatafield.at(ii) + nameSuffix + "Dom";
                dfDomain.push_back(ScaFES::DataField<CT, DIM>(
                    nameDf, &(this->mParams), this->globalGrid(), memAll.at(ii),
                    stencilWidth.at(ii), elemDataCurr, nColumns,
                    memHolePart.at(ii), nLayers.at(ii), writeToFile.at(ii)));
                elemDataCurr += (nRowsDf[ii] * nColumns);
                ++idDf;
            }
        }
    }

    idDf = 0;
    for (std::size_t ii = 0; ii < nLayers.size(); ++ii)
    {
        if (!isKnownDf.at(ii))
        {
            if (0 < nLayers.at(ii))
            {
                std::string nameDf = nameDatafield.at(ii) + nameSuffix + "Bdry";
                dfBdry.push_back(ScaFES::DataField<CT, DIM>(
                    nameDf, &(this->mParams), this->globalGrid(), memAll.at(ii),
                    stencilWidth.at(ii), elemDataCurr, nColumns,
                    memHolePart.at(ii), nLayers.at(ii), writeToFile.at(ii)));
                elemDataCurr += (nRowsDf[ii] * nColumns);
                ++idDf;
            }
        }
    }

    /*------------------------------------------------------------------------*/
    // Set default value for each unknown data field.
    idDf = 0;
    for (std::size_t ii = 0; ii < nLayers.size(); ++ii)
    {
        if (0 == nLayers.at(ii))
        {
            if (!(isKnownDf.at(ii)))
            {
                dfDomain.at(idDf) = defaultValue.at(ii);
                ++idDf;
            }
        }
    }
    idDf = 0;
    for (std::size_t ii = 0; ii < nLayers.size(); ++ii)
    {
        if (0 < nLayers.at(ii))
        {
            if (!(isKnownDf.at(ii)))
            {
                dfBdry.at(idDf) = defaultValue.at(ii);
                ++idDf;
            }
        }
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
template <typename TT>
inline void Problem<OWNPRBLM, CT, DIM>::copyDataFieldStructure(
    ScaFES_VectDf<TT>& aDfDep, const ScaFES_VectDf<CT>& dfDep,
    const std::vector<TT>& memoryAdfDep)
{
    // Copy constructor would not work because the types of the data fields
    // are different (TT vs. CT).
    aDfDep.reserve(dfDep.size());
    TT* elemDataDepCurr = const_cast<TT*>(memoryAdfDep.data());
    for (std::size_t ii = 0; ii < dfDep.size(); ++ii)
    {
        aDfDep.push_back(ScaFES::DataField<TT, DIM>(
            dfDep[ii].name(), dfDep[ii].params(), dfDep[ii].gridGlobal(),
            dfDep[ii].memAll(), dfDep[ii].stencilWidth(), elemDataDepCurr, 1,
            dfDep[ii].memHole(), dfDep[ii].nLayers()));
        elemDataDepCurr += dfDep[ii].nElemsAllocated();
    }
}

/*******************************************************************************
 * GETTER METHODS.
 ******************************************************************************/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const ScaFES::Parameters& Problem<OWNPRBLM, CT, DIM>::params() const
{
    return this->mParams;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const int& Problem<OWNPRBLM, CT, DIM>::nTimesteps() const
{
    return this->params().nTimesteps();
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const double& Problem<OWNPRBLM, CT, DIM>::gridsize(const int& idx) const
{
    return (this->mGG.discreteDomain().gridsize(idx));
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const double& Problem<OWNPRBLM, CT, DIM>::tau() const
{
    return this->params().tau();
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const ScaFES::GridGlobal<DIM>&
Problem<OWNPRBLM, CT, DIM>::globalGrid() const
{
    return this->mGG;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const short int&
Problem<OWNPRBLM, CT, DIM>::kind(const ScaFES_IntNtuple& idxNode) const
{
    return this->mKind[idxNode];
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const ScaFES::DataField<short int, DIM>&
Problem<OWNPRBLM, CT, DIM>::kind() const
{
    return this->mKind;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const typename ScaFES::Problem<OWNPRBLM, CT, DIM>::ScaFES_IntNtuple&
Problem<OWNPRBLM, CT, DIM>::nNodes() const
{
    return this->globalGrid().discreteDomain().nNodes();
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const int& Problem<OWNPRBLM, CT, DIM>::nNodes(const int& dim) const
{
    return this->nNodes().elem(dim);
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline ScaFES::Ntuple<double, DIM>
Problem<OWNPRBLM, CT, DIM>::coordinates(const ScaFES_IntNtuple& idxNode) const
{
    return this->globalGrid().discreteDomain().coordinates(idxNode);
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline double Problem<OWNPRBLM, CT, DIM>::time(const int& timestep) const
{
    return this->params().timeIntervalStart() + timestep * this->params().tau();
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline typename ScaFES::Problem<OWNPRBLM, CT, DIM>::ScaFES_IntNtuple
Problem<OWNPRBLM, CT, DIM>::connect(const ScaFES_IntNtuple& idxNode,
                                    const int& dir)
{
    return this->mGG.discreteDomain().neigh(idxNode, dir);
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const std::vector<CT>& Problem<OWNPRBLM, CT, DIM>::vectParamsCurr() const
{
    return this->mVectParamsCurr;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const CT&
Problem<OWNPRBLM, CT, DIM>::vectSol(const int& idx,
                                    const ScaFES_IntNtuple& idxNode) const
{
    return this->mVectKnownDfDom[idx](idxNode);
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const CT&
Problem<OWNPRBLM, CT, DIM>::knownDf(const int& idx,
                                    const ScaFES_IntNtuple& idxNode) const
{
    return this->mVectKnownDfDom[idx](idxNode);
}

/*******************************************************************************
 * WORK METHODS.
 ******************************************************************************/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::iterateOverTime()
{
    ScaFES::Timer timerIterate;
    int currTimeIter = 0;
    this->worldBarrier();
    timerIterate.restart();

    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Iterate over time..." << std::endl;
    }
    this->mParams.decreaseLevel();

    /*------------------------------------------------------------------------*/
#ifdef VTRACE
    unsigned int mid;
    mid = VT_MARKER_DEF("Prblm_iterate", VT_MARKER_TYPE_HINT);
    VT_MARKER(mid, "Prblm_iterate_start");
#endif

    for (currTimeIter = 0; currTimeIter <= this->nTimesteps(); ++currTimeIter)
    {

        /*--------------------------------------------------------------------*/
        if ((this->params().rankOutput() == this->myRank()) &&
            (0 < this->params().indentDepth()))
        {
            std::cout << this->mParams.getPrefix()
                      << " * Process TIME ITERATE = " << currTimeIter << "..."
                      << std::endl;
        }
        this->mParams.decreaseLevel();

        /*--------------------------------------------------------------------*/
        this->evalKnownDfs(currTimeIter);
        if (0 < currTimeIter)
        {
            if (this->params().enabledAdolc())
            {
                this->traceUpdatePhase(currTimeIter);
            }
            this->evalUpdatePhase(currTimeIter);
            if (this->useLeapfrogIntegration())
            {
                // Swap old / new data fields, BUT DO NOT SWAP times.
                // ==> Swap times 2x.
                this->swapPointersOfUnknownDfs(); // 1st swap of times.
                this->swapTimesOfUnknownDfs();    // 2nd swap of times.
                if (this->params().enabledAdolc())
                {
                    this->traceUpdate2Phase(currTimeIter);
                }
                this->evalUpdate2Phase(currTimeIter);
            }
        }
        else
        {
            this->resetValsAtNodesOfUnknownDfs(currTimeIter);
            this->resetTimesOfUnknownDfs();
            if (this->params().enabledAdolc())
            {
                this->traceInitPbPhase(currTimeIter);
            }
            this->evalInitPbPhase();
        }
        this->writeAllDfs(currTimeIter);
        this->compErrOfUnknownDfs();
        this->swapPointersOfUnknownDfs();
        // Update times of old / new data fields AFTER swapping.
        this->updateTimesOfUnknownDfs();

        /*--------------------------------------------------------------------*/
        this->mParams.increaseLevel();
        if ((this->params().rankOutput() == this->myRank()) &&
            (0 < this->params().indentDepth()))
        {
            std::cout << this->mParams.getPrefix() << "   Processed."
                      << std::endl;
        }
    }

    /*------------------------------------------------------------------------*/
    this->worldBarrier();
    double tIterate = timerIterate.elapsed();
    // Print some statistics.
    double tSum = 0.0;
    for (auto it = this->mWallClockTimes.begin();
         it != this->mWallClockTimes.end(); ++it)
    {
        tSum += it->second;
    }
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix() << std::scientific
                  << std::fixed
                  << " * Summary of Wall clock times(iterateOverTime):"
                  << std::endl;
        for (auto it = this->mWallClockTimes.begin();
             it != this->mWallClockTimes.end(); ++it)
        {
            std::cout << this->mParams.getPrefix() << "       t("
                      << std::setw(11) << it->first << ") = " << std::setw(8)
                      << std::setprecision(5) << it->second << " s ("
                      << std::setw(5) << std::setprecision(2)
                      << it->second / tSum * 100.0 << " %%)" << std::endl;
        }
        std::cout << this->mParams.getPrefix()
                  << "       --------------------------------------"
                  << std::endl
                  << this->mParams.getPrefix()
                  << "       t(    iterate) = " << std::setw(8)
                  << std::setprecision(5) << tSum << " s."
                  << std::resetiosflags(std::ios::fixed | std::ios::showpoint)
                  << std::endl
                  << this->mParams.getPrefix()
                  << " t(iterate) w.barriers= " << std::setw(8)
                  << std::setprecision(5) << tIterate << " s."
                  << std::resetiosflags(std::ios::fixed | std::ios::showpoint)
                  << std::endl;
    }
    this->mParams.increaseLevel();
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix() << "   Iterated."
                  << std::endl;
    }
#ifdef VTRACE
    mid = VT_MARKER_DEF("EndLoop", VT_MARKER_TYPE_HINT);
    VT_MARKER(mid, "Prblm_iterate_end");
#endif
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::setCommunicationMode()
{
    this->mUseAsynchronMode = this->params().asynchronMode();
    if (0 < this->nUnknownBdryDfs() && this->mUseAsynchronMode)
    {
        this->mUseAsynchronMode = false; // Force synchronous MPI communication.
        if (this->params().rankOutput() == this->myRank())
        {
            std::cout
                << "\n *** WARNING: There exists at least one border df!"
                << "===> Synchronous MPI communication must be used! ***\n"
                << std::endl;
        }
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::setCompGradients()
{
    this->mComputeGradients = this->params().computeGradients();
    if (!(this->params().enabledAdolc()) && this->mComputeGradients)
    {
        this->mComputeGradients = false;
        if (this->params().rankOutput() == this->myRank())
        {
            std::cout << "\n *** WARNING: ADOL-C is not enabled!"
                      << "===> Gradients can not be computed! ***\n"
                      << std::endl;
        }
    }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::evalInitPbPhase()
{
    std::string whichPhase = "init";

    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix() << " * Evaluate phase "
                  << whichPhase << "..." << std::endl;
    }
    this->mParams.decreaseLevel();

    int timeIter = 0;

    if (mParams.readInitFile() == false)
    {
        if (this->useAsynchronMode())
        {
            this->evalTypeOneAsynchronously(
                whichPhase, this->mVectUnknownDfsDomNew,
                this->mVectGradUnknownDfsDomNew, this->vectParamsCurr(),
                this->mVectGradParamsCurr, timeIter,
                this->tapeIdInitPartitionInner(), this->tapeIdInitPartitionBorder(),
                &OWNPRBLM::template initInner<CT>,
                &OWNPRBLM::template initBorder<CT>);
        }
        else
        {
            this->evalTypeOneSynchronously(
                whichPhase, this->mVectUnknownDfsDomNew,
                this->mVectGradUnknownDfsDomNew, this->mVectUnknownDfsBdryNew,
                this->mVectGradUnknownDfsBdryNew, this->vectParamsCurr(),
                this->mVectGradParamsCurr, timeIter, this->tapeIdInitPartitionAll(),
                &OWNPRBLM::template initInner<CT>,
                &OWNPRBLM::template initBorder<CT>);
        }
    }
    else
    {
        /* Code for reading init file. */
    }

    for (std::size_t ii = 0; ii < this->mVectUnknownDfsDomNew.size(); ++ii)
    {
        this->mVectUnknownDfsDomOld[ii]
            .setValuesAtMemAll(this->mVectUnknownDfsDomNew[ii]);
    }
    for (std::size_t ii = 0; ii < this->mVectUnknownDfsBdryNew.size(); ++ii)
    {
        this->mVectUnknownDfsBdryOld[ii]
            .setValuesAtMemAll(this->mVectUnknownDfsBdryNew[ii]);
    }

    this->mParams.increaseLevel();
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix() << "   Evaluated."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::evalUpdatePhase(const int& timeIter)
{
    std::string whichPhase = "update";

    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix() << " * Evaluate phase "
                  << whichPhase << "..." << std::endl;
    }
    this->mParams.decreaseLevel();

    if (this->useAsynchronMode())
    {
        this->evalTypeTwoAsynchronously(
            whichPhase, this->mVectUnknownDfsDomNew,
            this->mVectGradUnknownDfsDomNew, this->vectUnknownDfsDomOld(),
            this->mVectGradUnknownDfsDomOld, timeIter,
            this->tapeIdUpdatePartitionInner(),
            this->tapeIdUpdatePartitionBorder(),
            &OWNPRBLM::template updateInner<CT>,
            &OWNPRBLM::template updateBorder<CT>);
    }
    else
    {
        this->evalTypeTwoSynchronously(
            whichPhase, this->mVectUnknownDfsDomNew,
            this->mVectGradUnknownDfsDomNew, this->mVectUnknownDfsBdryNew,
            this->mVectGradUnknownDfsBdryNew, this->vectUnknownDfsDomOld(),
            this->mVectGradUnknownDfsDomOld, this->vectUnknownDfsBdryOld(),
            this->mVectGradUnknownDfsBdryOld, timeIter,
            this->tapeIdUpdatePartitionAll(),
            &OWNPRBLM::template updateInner<CT>,
            &OWNPRBLM::template updateBorder<CT>);
    }

    this->mParams.increaseLevel();
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix() << "   Evaluated."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::evalUpdate2Phase(const int& timeIter)
{
    std::string whichPhase = "update2";

    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix() << " * Evaluate phase "
                  << whichPhase << "..." << std::endl;
    }
    this->mParams.decreaseLevel();

    if (this->useAsynchronMode())
    {
        this->evalTypeTwoAsynchronously(
            whichPhase, this->mVectUnknownDfsDomNew,
            this->mVectGradUnknownDfsDomNew, this->vectUnknownDfsDomOld(),
            this->mVectGradUnknownDfsDomOld, timeIter,
            this->tapeIdUpdate2PartitionInner(),
            this->tapeIdUpdate2PartitionBorder(),
            &OWNPRBLM::template updateInner2<CT>,
            &OWNPRBLM::template updateBorder2<CT>);
    }
    else
    {
        this->evalTypeTwoSynchronously(
            whichPhase, this->mVectUnknownDfsDomNew,
            this->mVectGradUnknownDfsDomNew, this->mVectUnknownDfsBdryNew,
            this->mVectGradUnknownDfsBdryNew, this->vectUnknownDfsDomOld(),
            this->mVectGradUnknownDfsDomOld, this->vectUnknownDfsBdryOld(),
            this->mVectGradUnknownDfsBdryOld, timeIter,
            this->tapeIdUpdate2PartitionAll(),
            &OWNPRBLM::template updateInner2<CT>,
            &OWNPRBLM::template updateBorder2<CT>);
    }

    this->mParams.increaseLevel();
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix() << "   Evaluated."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void
Problem<OWNPRBLM, CT, DIM>::resetValsAtNodesOfUnknownDfs(const int& /*timeIter*/
                                                         )
{
    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Reset all unknown data fields..." << std::endl;
    }
    this->mParams.decreaseLevel();

    /*------------------------------------------------------------------------*/
    // Set default value for each unknown data field.
    int idDf = 0;
    for (std::size_t ii = 0; ii < this->nLayers().size(); ++ii)
    {
        if (0 == this->nLayers().at(ii))
        {
            if (!(this->isKnownDf().at(ii)))
            {
                this->mVectUnknownDfsDomOld.at(idDf) =
                    this->defaultValue().at(ii);
                ++idDf;
            }
        }
    }
    idDf = 0;
    for (std::size_t ii = 0; ii < this->nLayers().size(); ++ii)
    {
        if (0 < this->nLayers().at(ii))
        {
            if (!(this->isKnownDf().at(ii)))
            {
                this->mVectUnknownDfsBdryOld.at(idDf) =
                    this->defaultValue().at(ii);
                ++idDf;
            }
        }
    }
    idDf = 0;
    for (std::size_t ii = 0; ii < this->nLayers().size(); ++ii)
    {
        if (0 == this->nLayers().at(ii))
        {
            if (!(this->isKnownDf().at(ii)))
            {
                this->mVectUnknownDfsDomNew.at(idDf) =
                    this->defaultValue().at(ii);
                ++idDf;
            }
        }
    }
    idDf = 0;
    for (std::size_t ii = 0; ii < this->nLayers().size(); ++ii)
    {
        if (0 < this->nLayers().at(ii))
        {
            if (!(this->isKnownDf().at(ii)))
            {
                this->mVectUnknownDfsBdryNew.at(idDf) =
                    this->defaultValue().at(ii);
                ++idDf;
            }
        }
    }

    this->mParams.increaseLevel();
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix() << "   Reset."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::resetTimesOfUnknownDfs()
{
    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Reset times of all unknown data fields..."
                  << std::endl;
    }
    this->mParams.decreaseLevel();

    /*------------------------------------------------------------------------*/
    // Update current time of "new" data fields AFTER swapping.
    for (std::size_t ii = 0; ii < this->mVectUnknownDfsDomNew.size(); ++ii)
    {
        this->mVectUnknownDfsDomNew[ii].time() =
            this->params().timeIntervalStart();
    }
    for (std::size_t ii = 0; ii < this->mVectUnknownDfsBdryNew.size(); ++ii)
    {
        this->mVectUnknownDfsBdryNew[ii].time() =
            this->params().timeIntervalStart();
    }
    for (std::size_t ii = 0; ii < this->mVectGradUnknownDfsDomNew.size(); ++ii)
    {
        this->mVectGradUnknownDfsDomNew[ii] = 1.0;
    }
    for (std::size_t ii = 0; ii < this->mVectGradUnknownDfsBdryNew.size(); ++ii)
    {
        this->mVectGradUnknownDfsBdryNew[ii] = 1.0;
    }

    this->mParams.increaseLevel();
    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix() << "   Reset."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::updateTimesOfUnknownDfs()
{
    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Update times of all unknown data fields..."
                  << std::endl;
    }
    this->mParams.decreaseLevel();

    /*------------------------------------------------------------------------*/
    // Update current time of "new" data fields AFTER swapping.
    for (std::size_t ii = 0; ii < this->mVectUnknownDfsDomNew.size(); ++ii)
    {
        this->mVectUnknownDfsDomNew[ii].time() =
            this->mVectUnknownDfsDomOld[ii].time() + this->params().tau();
    }
    for (std::size_t ii = 0; ii < this->mVectUnknownDfsBdryNew.size(); ++ii)
    {
        this->mVectUnknownDfsBdryNew[ii].time() =
            this->mVectUnknownDfsBdryOld[ii].time() + this->params().tau();
    }
    for (std::size_t ii = 0; ii < this->mVectGradUnknownDfsDomNew.size(); ++ii)
    {
        this->mVectGradUnknownDfsDomNew[ii] = 1.0;
    }
    for (std::size_t ii = 0; ii < this->mVectGradUnknownDfsBdryNew.size(); ++ii)
    {
        this->mVectGradUnknownDfsBdryNew[ii] = 1.0;
    }

    this->mParams.increaseLevel();
    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix() << "   Updated."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::swapPointersOfUnknownDfs()
{
    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Swap pointers of all unknown data fields..."
                  << std::endl;
    }
    this->mParams.decreaseLevel();

    /*------------------------------------------------------------------------*/
    // Swap old and new iterate.
    for (std::size_t ii = 0; ii < this->vectUnknownDfsDomNew().size(); ++ii)
    {
        this->mVectUnknownDfsDomNew[ii]
            .swapPointerToMemoryWith(this->mVectUnknownDfsDomOld[ii]);
    }
    for (std::size_t ii = 0; ii < this->vectUnknownDfsBdryNew().size(); ++ii)
    {
        this->mVectUnknownDfsBdryNew[ii]
            .swapPointerToMemoryWith(this->mVectUnknownDfsBdryOld[ii]);
    }
    /*------------------------------------------------------------------------*/
    // Swap gradients of old and new iterate.
    for (std::size_t ii = 0; ii < this->mVectGradUnknownDfsDomNew.size(); ++ii)
    {
        this->mVectGradUnknownDfsDomNew[ii]
            .swapPointerToMemoryWith(this->mVectGradUnknownDfsDomOld[ii]);
    }
    for (std::size_t ii = 0; ii < this->mVectGradUnknownDfsBdryNew.size(); ++ii)
    {
        this->mVectGradUnknownDfsBdryNew[ii]
            .swapPointerToMemoryWith(this->mVectGradUnknownDfsBdryOld[ii]);
    }

    this->mParams.increaseLevel();
    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix() << "   Swapped."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::swapTimesOfUnknownDfs()
{
    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Swap times of all unknown data fields..." << std::endl;
    }
    this->mParams.decreaseLevel();

    double tmpTime = 0.0;
    for (std::size_t ii = 0; ii < this->mVectUnknownDfsDomNew.size(); ++ii)
    {
        tmpTime = this->mVectUnknownDfsDomOld[ii].time();
        this->mVectUnknownDfsDomOld[ii].time() =
            this->mVectUnknownDfsDomNew[ii].time();
        this->mVectUnknownDfsDomNew[ii].time() = tmpTime;
    }
    for (std::size_t ii = 0; ii < this->mVectUnknownDfsBdryNew.size(); ++ii)
    {
        tmpTime = this->mVectUnknownDfsBdryOld[ii].time();
        this->mVectUnknownDfsBdryOld[ii].time() =
            this->mVectUnknownDfsBdryNew[ii].time();
        this->mVectUnknownDfsBdryNew[ii].time() = tmpTime;
    }
    for (std::size_t ii = 0; ii < this->mVectGradUnknownDfsDomNew.size(); ++ii)
    {
        tmpTime = this->mVectGradUnknownDfsDomOld[ii].time();
        this->mVectGradUnknownDfsDomOld[ii].time() =
            this->mVectGradUnknownDfsDomNew[ii].time();
        this->mVectGradUnknownDfsDomNew[ii].time() = tmpTime;
    }
    for (std::size_t ii = 0; ii < this->mVectGradUnknownDfsBdryNew.size(); ++ii)
    {
        tmpTime = this->mVectGradUnknownDfsBdryOld[ii].time();
        this->mVectGradUnknownDfsBdryOld[ii].time() =
            this->mVectGradUnknownDfsBdryNew[ii].time();
        this->mVectGradUnknownDfsBdryNew[ii].time() = tmpTime;
    }

    this->mParams.increaseLevel();
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << "   Swapped."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::writeAllDfs(const int& timeIter)
{
    ScaFES::Timer timerWrite;
    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Write all data fields to file..." << std::endl;
    }
    this->mParams.decreaseLevel();

    timerWrite.restart();
    this->writeDfsToFile(timeIter);
    this->mKind.write(timeIter);
    this->mWallClockTimes["write"] += timerWrite.elapsed();
    this->mParams.increaseLevel();
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                                       << "   Written."
                  << std::endl;
    }
}

/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::evalKnownDfs(const int& timeIter)
{
    bool evalDf = false;
    for (std::size_t ii = 0; ii < this->mComputeError.size(); ++ii)
    {
        if (this->mComputeError.at(ii) || this->isKnownDf().at(ii))
        {
            evalDf = true;
            break;
        }
    }
    if (!evalDf)
    {
        return;
    }

    ScaFES::Timer timerEval;
    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Evaluate all known data fields..." << std::endl;
    }
    this->mParams.decreaseLevel();

    timerEval.restart();

    this->evalGlobInner(this->mVectKnownDfDom, timeIter);
    this->evalGlobBorder(this->mVectKnownDfDom, timeIter);
    for (std::size_t ii = 0; ii < this->mVectKnownDfDom.size(); ++ii)
    {
        this->mVectKnownDfDom[ii].sync(timeIter);
    }
    if (0 < this->nUnknownBdryDfs())
    {
        this->evalGlobBorder(this->mVectKnownDfBdry, timeIter);
        for (std::size_t ii = 0; ii < this->mVectKnownDfBdry.size(); ++ii)
        {
            this->mVectKnownDfBdry[ii].sync(timeIter);
        }
    }
    this->mWallClockTimes["evalKnownDf"] += timerEval.elapsed();

    this->mParams.increaseLevel();
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                                       << "   Evaluated."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::compErrOfUnknownDfs()
{
    ScaFES::Timer timerCompErr;
    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Compute error for all unknown data fields..."
                  << std::endl;
    }
    this->mParams.decreaseLevel();

    timerCompErr.restart();
    // Compute nodal error.
    int idDf = 0;
    for (std::size_t ii = 0; ii < this->mComputeError.size(); ++ii)
    {
        if (0 == this->nLayers().at(ii))
        {
            if (!(this->isKnownDf().at(ii)))
            {
                if (this->mComputeError.at(ii))
                {
                    this->mVectUnknownDfsDomNew[idDf]
                        .compErrLinf(this->mVectKnownDfDom[ii]);
                }
                ++idDf;
            }
        }
    }
    if (0 < this->nUnknownBdryDfs())
    {
        idDf = 0;
        for (std::size_t ii = 0; ii < this->mComputeError.size(); ++ii)
        {
            if (0 < this->nLayers().at(ii))
            {
                if (!(this->isKnownDf().at(ii)))
                {
                    if (this->mComputeError.at(ii))
                    {
                        this->mVectUnknownDfsBdryNew[idDf]
                            .compErrLinf(this->mVectKnownDfBdry[ii]);
                    }
                    ++idDf;
                }
            }
        }
    }
    this->mWallClockTimes["compErr"] += timerCompErr.elapsed();

    this->mParams.increaseLevel();
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                                       << "   Computed."
                  << std::endl;
    }
}

/*******************************************************************************
 * SETTER METHODS.
 ******************************************************************************/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::increaseLevel()
{
    this->mParams.increaseLevel();
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::decreaseLevel()
{
    this->mParams.decreaseLevel();
}

/*******************************************************************************
 * GETTER METHODS.
 ******************************************************************************/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const bool& Problem<OWNPRBLM, CT, DIM>::isTracedInit() const
{
    return this->mIsTracedInit;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline short int Problem<OWNPRBLM, CT, DIM>::tapeIdInitPartitionInner()
{
    return this->mTapeIdInitPartitionInner;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline short int Problem<OWNPRBLM, CT, DIM>::tapeIdInitPartitionBorder()
{
    return this->mTapeIdInitPartitionBorder;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline short int Problem<OWNPRBLM, CT, DIM>::tapeIdInitPartitionAll()
{
    return this->mTapeIdInitPartitionAll;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const bool& Problem<OWNPRBLM, CT, DIM>::isTracedUpdatePhase() const
{
    return this->mIsTracedUpdatePhase;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const short int&
Problem<OWNPRBLM, CT, DIM>::tapeIdUpdatePartitionInner() const
{
    return this->mTapeIdUpdatePartitionInner;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const short int&
Problem<OWNPRBLM, CT, DIM>::tapeIdUpdatePartitionBorder() const
{
    return this->mTapeIdUpdatePartitionBorder;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const short int&
Problem<OWNPRBLM, CT, DIM>::tapeIdUpdatePartitionAll() const
{
    return this->mTapeIdUpdatePartitionAll;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const bool& Problem<OWNPRBLM, CT, DIM>::isTracedUpdate2Phase() const
{
    return this->mIsTracedUpdate2Phase;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const short int&
Problem<OWNPRBLM, CT, DIM>::tapeIdUpdate2PartitionInner() const
{
    return this->mTapeIdUpdate2PartitionInner;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const short int&
Problem<OWNPRBLM, CT, DIM>::tapeIdUpdate2PartitionBorder() const
{
    return this->mTapeIdUpdate2PartitionBorder;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const short int&
Problem<OWNPRBLM, CT, DIM>::tapeIdUpdate2PartitionAll() const
{
    return this->mTapeIdUpdate2PartitionAll;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline int Problem<OWNPRBLM, CT, DIM>::myRank() const
{
    return this->params().myWorld().rank();
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const std::vector<bool>& Problem<OWNPRBLM, CT, DIM>::isKnownDf() const
{
    return this->mIsKnownDf;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const std::vector<int>& Problem<OWNPRBLM, CT, DIM>::nLayers() const
{
    return this->mNlayers;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const std::vector<CT>& Problem<OWNPRBLM, CT, DIM>::defaultValue() const
{
    return this->mDefaultValue;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const int& Problem<OWNPRBLM, CT, DIM>::nElemsUnknownDfs() const
{
    return this->mNelemsUnknownDfs;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const int& Problem<OWNPRBLM, CT, DIM>::nElemsKnownDfs() const
{
    return this->mNelemsKnownDfs;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const int& Problem<OWNPRBLM, CT, DIM>::nUnknownBdryDfs() const
{
    return this->mNunknownBdryDfs;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const int& Problem<OWNPRBLM, CT, DIM>::nNodesTotalPartition() const
{
    return this->mNnodesTotalPartition;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const std::vector<DataField<CT, DIM>>&
Problem<OWNPRBLM, CT, DIM>::vectUnknownDfsDomNew() const
{
    return this->mVectUnknownDfsDomNew;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const std::vector<DataField<CT, DIM>>&
Problem<OWNPRBLM, CT, DIM>::vectUnknownDfsDomOld() const
{
    return this->mVectUnknownDfsDomOld;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const std::vector<DataField<CT, DIM>>&
Problem<OWNPRBLM, CT, DIM>::vectUnknownDfsBdryNew() const
{
    return this->mVectUnknownDfsBdryNew;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const std::vector<DataField<CT, DIM>>&
Problem<OWNPRBLM, CT, DIM>::vectUnknownDfsBdryOld() const
{
    return this->mVectUnknownDfsBdryOld;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const int& Problem<OWNPRBLM, CT, DIM>::nColumns() const
{
    return this->mNparams;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const bool& Problem<OWNPRBLM, CT, DIM>::computeGradients() const
{
    return this->mComputeGradients;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const bool& Problem<OWNPRBLM, CT, DIM>::useAsynchronMode() const
{
    return this->mUseAsynchronMode;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const bool& Problem<OWNPRBLM, CT, DIM>::isParamsBased() const
{
    return this->mIsParamsBased;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline const bool& Problem<OWNPRBLM, CT, DIM>::useLeapfrogIntegration() const
{
    return this->mUseLeapfrogIntegration;
}

/*******************************************************************************
 * WORK METHODS.
 ******************************************************************************/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::worldBarrier()
{
    mParams.myWorld().barrier();
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::setGridGlobalFlags()
{
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Set grid global related flags..." << std::endl;
    }
    this->mParams.decreaseLevel();

    int nLayersAtBorder = this->params().nLayersAtBorder();
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * nLayersAtGlobalBorder = " << nLayersAtBorder
                  << std::endl;
    }

    int comp = 0;
    int ii;
// * ii Iteration variable ---> private
// * idxNode WRITE access Corresponding global node number ---> private
#ifdef _OPENMP
#pragma omp parallel for schedule(static) private(ii) shared(comp)
#endif
    // Default: Set all nodes of the grid partition as MPI global inner nodes.
    for (ii = 0; ii < this->nNodesTotalPartition(); ++ii)
    {
        ScaFES_IntNtuple idxNode;
        this->mKind.memoryPos2IdxNode(idxNode, comp, ii);
        this->mKind(idxNode) = this->globalGrid().SCAFES_GLOBAL_INNER_NODE;
    }

    ScaFES_IntNtuple idxNodeFirstGlobal = this->globalGrid().discreteDomain().idxNodeFirst();
    ScaFES_IntNtuple idxNodeLastGlobal = this->globalGrid().discreteDomain().idxNodeLast();
    ScaFES_IntNtuple idxNodeFirstPart =
        this->globalGrid().partition(this->myRank()).idxNodeFirstSub();
    ScaFES_IntNtuple idxNodeLastPart =
        this->globalGrid().partition(this->myRank()).idxNodeLastSub();

// * ii Iteration variable ---> private
// * idxNode WRITE access Corresponding global node number ---> private
// * nLayersAtBorder READ access --> shared
// * idxNodeFirst... READ access --> shared
#ifdef _OPENMP
#pragma omp parallel for schedule(static) private(ii)            \
    shared(comp, nLayersAtBorder, idxNodeFirstGlobal, idxNodeLastGlobal)
#endif
    for (ii = 0; ii < this->nNodesTotalPartition(); ++ii)
    {
        ScaFES_IntNtuple idxNode;
        this->mKind.memoryPos2IdxNode(idxNode, comp, ii);

        for (std::size_t jj = 0; jj < DIM; ++jj)
        {
            if (((idxNodeFirstGlobal.elem(jj) + nLayersAtBorder) >
                 idxNode.elem(jj)) ||
                ((idxNodeLastGlobal.elem(jj) - nLayersAtBorder) <
                 idxNode.elem(jj)))
            {
                this->mKind(idxNode) +=
                    this->globalGrid().SCAFES_GLOBAL_BORDER_NODE;
                this->mKind(idxNode) -=
                    this->globalGrid().SCAFES_GLOBAL_INNER_NODE;
                break;
            }
        }
    }
    this->mParams.increaseLevel();
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                                       << "   Set."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void
Problem<OWNPRBLM, CT, DIM>::setGridPartitionFlags(const int& maxStencilWidth)
{
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Set grid partition related flags..." << std::endl;
    }
    this->mParams.decreaseLevel();

    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * maxStencilWidth = " << maxStencilWidth << std::endl;
    }

    int comp = 0;
    int ii;
// * ii Iteration variable ---> private
// * idxNode WRITE access Corresponding global node number ---> private
#ifdef _OPENMP
#pragma omp parallel for schedule(static) private(ii) shared(comp)
#endif
    // Default: Set all nodes of the grid partition as MPI local inner nodes.
    for (ii = 0; ii < this->nNodesTotalPartition(); ++ii)
    {
        ScaFES_IntNtuple idxNode;
        this->mKind.memoryPos2IdxNode(idxNode, comp, ii);
        this->mKind(idxNode) += this->globalGrid().SCAFES_PARTITION_INNER_NODE;
    }
    ScaFES_IntNtuple idxNodeFirstGlobal = this->globalGrid().discreteDomain().idxNodeFirst();
    ScaFES_IntNtuple idxNodeLastGlobal = this->globalGrid().discreteDomain().idxNodeLast();
    ScaFES_IntNtuple idxNodeFirstPart =
        this->globalGrid().partition(this->myRank()).idxNodeFirstSub();
    ScaFES_IntNtuple idxNodeLastPart =
        this->globalGrid().partition(this->myRank()).idxNodeLastSub();

// * ii Iteration variable ---> private
// * idxNode WRITE access Corresponding global node number ---> private
// * maxStencilWidth READ access --> shared
// * idxNodeFirst... READ access --> shared
#ifdef _OPENMP
#pragma omp parallel for schedule(static) private(ii)            \
    shared(comp, maxStencilWidth, idxNodeFirstPart, idxNodeLastPart)
#endif
    for (ii = 0; ii < this->nNodesTotalPartition(); ++ii)
    {
        ScaFES_IntNtuple idxNode;
        this->mKind.memoryPos2IdxNode(idxNode, comp, ii);

        for (std::size_t jj = 0; jj < DIM; ++jj)
        {
            if (((idxNodeFirstPart.elem(jj) + maxStencilWidth) >
                 idxNode.elem(jj)) ||
                ((idxNodeLastPart.elem(jj) - maxStencilWidth) <
                 idxNode.elem(jj)))
            {
                this->mKind(idxNode) +=
                    this->globalGrid().SCAFES_PARTITION_BORDER_NODE;
                this->mKind(idxNode) -=
                    this->globalGrid().SCAFES_PARTITION_INNER_NODE;
                break;
            }
        }
    }
    this->mParams.increaseLevel();
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                                       << "   Set."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
template <typename TT>
inline void Problem<OWNPRBLM, CT, DIM>::evalGlobInner(ScaFES_VectDf<TT>& dfDep,
                                                      const int& timeIter)
{
    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Evaluate data fields " << dfDep.at(0).name()
                  << "... at global inner at timeIter " << timeIter << "..."
                  << std::endl;
    }
    this->mParams.decreaseLevel();

    int comp = 0;
    std::size_t ii;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) shared(dfDep, timeIter) private(ii, comp)
#endif
    for (ii = 0; ii < this->mNodesLocBorderGlobInner.size(); ++ii)
    {
        ScaFES_IntNtuple idxNode;
        int idxNodeScalar = this->mNodesLocBorderGlobInner[ii];
        this->mKind.memoryPos2IdxNode(idxNode, comp, idxNodeScalar);
        static_cast<OWNPRBLM*>(this)->evalInner(dfDep, idxNode, timeIter);
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) shared(dfDep, timeIter) private(ii, comp)
#endif
    for (ii = 0; ii < this->mNodesLocInnerGlobInner.size(); ++ii)
    {
        ScaFES_IntNtuple idxNode;
        int idxNodeScalar = this->mNodesLocInnerGlobInner[ii];
        this->mKind.memoryPos2IdxNode(idxNode, comp, idxNodeScalar);
        static_cast<OWNPRBLM*>(this)->evalInner(dfDep, idxNode, timeIter);
    }

    this->mParams.increaseLevel();
    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << "   Evaluated."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
template <typename TT>
inline void Problem<OWNPRBLM, CT, DIM>::evalGlobBorder(ScaFES_VectDf<TT>& dfDep,
                                                       const int& timeIter)
{
    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Evaluate data fields " << dfDep.at(0).name()
                  << "... at global border at timeIter " << timeIter << "..."
                  << std::endl;
    }
    this->mParams.decreaseLevel();

    int comp = 0;
    std::size_t ii;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) private(ii, comp)      \
    shared(dfDep, timeIter)
#endif
    for (ii = 0; ii < this->mNodesLocBorderGlobBorder.size(); ++ii)
    {
        ScaFES_IntNtuple idxNode;
        int idxNodeScalar = this->mNodesLocBorderGlobBorder[ii];
        this->mKind.memoryPos2IdxNode(idxNode, comp, idxNodeScalar);
        static_cast<OWNPRBLM*>(this)->evalBorder(dfDep, idxNode, timeIter);
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) private(ii, comp)      \
    shared(dfDep, timeIter)
#endif
    for (ii = 0; ii < this->mNodesLocInnerGlobBorder.size(); ++ii)
    {
        ScaFES_IntNtuple idxNode;
        int idxNodeScalar = this->mNodesLocInnerGlobBorder[ii];
        this->mKind.memoryPos2IdxNode(idxNode, comp, idxNodeScalar);
        static_cast<OWNPRBLM*>(this)->evalBorder(dfDep, idxNode, timeIter);
    }

    this->mParams.increaseLevel();
    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << "   Evaluated."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline int Problem<OWNPRBLM, CT, DIM>::getNnodesLocInnerGlobInner() const
{
    int nNodes = 0;
    for (int ii = 0; ii < this->nNodesTotalPartition(); ++ii)
    {
        if (this->mKind(ii) & this->globalGrid().SCAFES_PARTITION_INNER_NODE)
        {
            if (this->mKind(ii) & this->globalGrid().SCAFES_GLOBAL_INNER_NODE)
            {
                ++nNodes;
            }
        }
    }
    return nNodes;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline int Problem<OWNPRBLM, CT, DIM>::getNnodesLocInnerGlobBorder() const
{
    int nNodes = 0;
    for (int ii = 0; ii < this->nNodesTotalPartition(); ++ii)
    {
        if (this->mKind(ii) & this->globalGrid().SCAFES_PARTITION_INNER_NODE)
        {
            if (this->mKind(ii) & this->globalGrid().SCAFES_GLOBAL_BORDER_NODE)
            {
                ++nNodes;
            }
        }
    }
    return nNodes;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline int Problem<OWNPRBLM, CT, DIM>::getNnodesLocBorderGlobInner() const
{
    int nNodes = 0;
    for (int ii = 0; ii < this->nNodesTotalPartition(); ++ii)
    {
        if (this->mKind(ii) & this->globalGrid().SCAFES_PARTITION_BORDER_NODE)
        {
            if (this->mKind(ii) & this->globalGrid().SCAFES_GLOBAL_INNER_NODE)
            {
                ++nNodes;
            }
        }
    }
    return nNodes;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline int Problem<OWNPRBLM, CT, DIM>::getNnodesLocBorderGlobBorder() const
{
    int nNodes = 0;
    for (int ii = 0; ii < this->nNodesTotalPartition(); ++ii)
    {
        if (this->mKind(ii) & this->globalGrid().SCAFES_PARTITION_BORDER_NODE)
        {
            if (this->mKind(ii) & this->globalGrid().SCAFES_GLOBAL_BORDER_NODE)
            {
                ++nNodes;
            }
        }
    }
    return nNodes;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline int Problem<OWNPRBLM, CT, DIM>::getNnodesLocInner() const
{
    int nNodes = 0;
    for (int ii = 0; ii < this->nNodesTotalPartition(); ++ii)
    {
        if (this->mKind(ii) & this->globalGrid().SCAFES_PARTITION_INNER_NODE)
        {
            ++nNodes;
        }
    }
    return nNodes;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline int Problem<OWNPRBLM, CT, DIM>::getNnodesLocBorder() const
{
    int nNodes = 0;
    for (int ii = 0; ii < this->nNodesTotalPartition(); ++ii)
    {
        if (this->mKind(ii) & this->globalGrid().SCAFES_PARTITION_BORDER_NODE)
        {
            ++nNodes;
        }
    }
    return nNodes;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::writeDfsToFile(const int& timeIter)
{
    // Check if file should be written at current time step.
    int intervalOfWrites = this->params().nTimesteps();

    if (0 < this->params().nSnapshots())
    {
        if (this->params().nSnapshots() <= this->params().nTimesteps())
        {
            intervalOfWrites =
                this->params().nTimesteps() / this->params().nSnapshots();
        }
    }

    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Write data fields to file "
                  << mCommonFile.nameDataFile() << " ..." << std::endl;
    }
    this->mParams.decreaseLevel();

    // Do not consider those data fields which should be never written.
    int nDataFieldsToWrite = 0;
    for (std::size_t ii = 0; ii < this->mWriteToFile.size(); ++ii)
    {
        if (ScaFES::WriteHowOften::NEVER != this->mWriteToFile[ii])
        {
            ++nDataFieldsToWrite;
        }
    }

    if (0 < nDataFieldsToWrite)
    {
        int kk = 0;
        std::vector<bool> writeData(nDataFieldsToWrite, false);
        for (std::size_t ii = 0; ii < this->mWriteToFile.size(); ++ii)
        {
            if (ScaFES::WriteHowOften::NEVER != this->mWriteToFile[ii])
            {
                switch (this->mWriteToFile.at(ii))
                {
                case ScaFES::WriteHowOften::AT_START:
                    if (0 == timeIter)
                    {
                        writeData.at(kk) = true;
                    }
                    break;
                case ScaFES::WriteHowOften::AT_START_AND_END:
                    if (0 == timeIter ||
                        this->params().nTimesteps() == timeIter)
                    {
                        writeData.at(kk) = true;
                    }
                    break;
                case ScaFES::WriteHowOften::LIKE_GIVEN_AT_CL:
                    if (0 == timeIter % intervalOfWrites)
                    {
                        writeData.at(kk) = true;
                    }
                    break;
                case ScaFES::WriteHowOften::ALWAYS:
                default:
                    writeData.at(kk) = true;
                    break;
                }
                if (this->params().rankOutput() == this->myRank() &&
                    (0 < this->params().indentDepth()))
                {
                    std::cout << this->mParams.getPrefix()
                              << " * Status of data field #" << ii << ": "
                              << ((writeData.at(kk)) ? "write" : "do not write")
                              << std::endl;
                }
                ++kk;
            }
        }

        std::vector<CT*> tmpElemData;
        tmpElemData.reserve(nDataFieldsToWrite);
        int idxKnownDf = 0;
        int idxUnknownDf = 0;
        for (std::size_t ii = 0; ii < this->mWriteToFile.size(); ++ii)
        {
            if (0 == this->nLayers().at(ii))
            {
                if (this->isKnownDf().at(ii))
                {
                    if (ScaFES::WriteHowOften::NEVER != this->mWriteToFile[ii])
                    {
                        tmpElemData.push_back(
                            this->mVectKnownDfDom[idxKnownDf].elemData());
                        ++idxKnownDf;
                    }
                }
                else
                {
                    if (ScaFES::WriteHowOften::NEVER != this->mWriteToFile[ii])
                    {
                        tmpElemData.push_back(
                            this->vectUnknownDfsDomNew()[idxUnknownDf]
                                .elemData());
                    }
                    ++idxUnknownDf;
                }
            }
        }
        idxKnownDf = 0;
        idxUnknownDf = 0;
        for (std::size_t ii = 0; ii < this->mWriteToFile.size(); ++ii)
        {
            if (0 < this->nLayers().at(ii))
            {
                if (this->isKnownDf().at(ii))
                {
                    if (ScaFES::WriteHowOften::NEVER != this->mWriteToFile[ii])
                    {
                        tmpElemData.push_back(
                            this->mVectKnownDfBdry[idxKnownDf].elemData());
                        ++idxKnownDf;
                    }
                }
                else
                {
                    if (ScaFES::WriteHowOften::NEVER != this->mWriteToFile[ii])
                    {
                        tmpElemData.push_back(
                            this->vectUnknownDfsBdryNew()[idxUnknownDf]
                                .elemData());
                    }
                    ++idxUnknownDf;
                }
            }
        }

        this->mCommonFile.write(tmpElemData, writeData, timeIter);
    }

    this->mParams.increaseLevel();
    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << "   Written."
                  << std::endl;
    }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::evalTypeOneAsynchronously(
    const std::string& whichPhase, ScaFES_VectDf<CT>& dfDep,
    ScaFES_VectDf<CT>& dfGradDep, const std::vector<CT>& dfIndep,
    const std::vector<CT>& dfGradIndep, const int& timeIter,
    const int& tapeIdLocalInner, const int& tapeIdLocalBorder,
    ScaFES_FuncTypeOneToImplement<CT> funcPhaseInner,
    ScaFES_FuncTypeOneToImplement<CT> funcPhaseBorder)
{
    ScaFES::Timer timerSync;
    ScaFES::Timer timerPhase;

    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Evaluate phase "
                  << whichPhase << " asynchronously: " << dfDep.at(0).name()
                  << "..." << std::endl;
    }
    this->mParams.decreaseLevel();
    /*------------------------------------------------------------------------*/
    timerPhase.restart();
    if (this->params().enabledAdolc())
    {
        if (this->computeGradients())
        {
            this->evalFovForwardTypeOne(whichPhase, dfDep, dfGradDep,
                                        tapeIdLocalBorder, dfIndep, dfGradIndep,
                                        timeIter);
        }
        else
        {
            this->evalZosForwardTypeOne(whichPhase, dfDep, tapeIdLocalBorder,
                                        dfIndep, timeIter);
        }
    }
    else
    {
        this->evalTypeOneFuncAtPartBorder(whichPhase, dfDep, dfIndep, timeIter,
                                          funcPhaseInner, funcPhaseBorder);
    }
    this->mWallClockTimes[whichPhase] += timerPhase.elapsed();
    /*------------------------------------------------------------------------*/
    timerSync.restart();
    for (std::size_t ii = 0; ii < dfDep.size(); ++ii)
    {
        dfDep[ii].copyValuesFromMemCommToSendBuffer(timeIter);
    }
    for (std::size_t ii = 0; ii < dfGradDep.size(); ++ii)
    {
        dfGradDep[ii].copyValuesFromMemCommToSendBuffer(timeIter);
    }
    for (std::size_t ii = 0; ii < dfDep.size(); ++ii)
    {
        dfDep[ii].exchangeValuesInBuffers(timeIter);
    }
    for (std::size_t ii = 0; ii < dfGradDep.size(); ++ii)
    {
        dfGradDep[ii].exchangeValuesInBuffers(timeIter);
    }
    this->mWallClockTimes["sync"] += timerSync.elapsed();
    /*------------------------------------------------------------------------*/
    timerPhase.restart();
    if (this->params().enabledAdolc())
    {
        if (this->computeGradients())
        {
            this->evalFovForwardTypeOne(whichPhase, this->mVectUnknownDfsDomTmp,
                                        this->mVectGradUnknownDfsDomTmp,
                                        tapeIdLocalInner, dfIndep, dfGradIndep,
                                        timeIter);
        }
        else
        {
            this->evalZosForwardTypeOne(whichPhase, this->mVectUnknownDfsDomTmp,
                                        tapeIdLocalInner, dfIndep, timeIter);
        }
        for (std::size_t ii = 0; ii < dfDep.size(); ++ii)
        {
            dfDep[ii] += this->mVectUnknownDfsDomTmp[ii];
        }
        int idx = 0;
        for (std::size_t ii = 0; ii < this->isKnownDf().size(); ++ii)
        {
            if (!(this->isKnownDf().at(ii)))
            {
                dfDep[idx] -= this->defaultValue()[ii];
                ++idx;
            }
        }
        for (std::size_t ii = 0; ii < dfGradDep.size(); ++ii)
        {
            dfGradDep[ii] += this->mVectGradUnknownDfsDomTmp[ii]; // - 0;
        }
    }
    else
    {
        this->evalTypeOneFuncAtPartInner(whichPhase, dfDep, dfIndep, timeIter,
                                         funcPhaseInner, funcPhaseBorder);
    }
    this->mWallClockTimes[whichPhase] += timerPhase.elapsed();
    /*------------------------------------------------------------------------*/
    timerSync.restart();
    for (std::size_t ii = 0; ii < dfDep.size(); ++ii)
    {
        dfDep[ii].waitAll();
    }
    for (std::size_t ii = 0; ii < dfGradDep.size(); ++ii)
    {
        dfGradDep[ii].waitAll();
    }
    for (std::size_t ii = 0; ii < dfDep.size(); ++ii)
    {
        dfDep[ii].copyValuesFromReceiveBufferToMemGhost(timeIter);
    }
    for (std::size_t ii = 0; ii < dfGradDep.size(); ++ii)
    {
        dfGradDep[ii].copyValuesFromReceiveBufferToMemGhost(timeIter);
    }
    this->mWallClockTimes["sync"] += timerSync.elapsed();

    /*------------------------------------------------------------------------*/
    this->mParams.increaseLevel();
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << "   Evaluated."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::evalTypeTwoAsynchronously(
    const std::string& whichPhase, ScaFES_VectDf<CT>& dfDep,
    ScaFES_VectDf<CT>& dfGradDep, const ScaFES_VectDf<CT>& dfIndep,
    const ScaFES_VectDf<CT>& dfGradIndep, const int& timeIter,
    const int& tapeIdLocalInner, const int& tapeIdLocalBorder,
    ScaFES_FuncTypeTwoToImplement<CT> funcPhaseInner,
    ScaFES_FuncTypeTwoToImplement<CT> funcPhaseBorder)
{
    ScaFES::Timer timerSync;
    ScaFES::Timer timerPhase;

    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Evaluate phase "
                  << whichPhase << " asynchronously: " << dfDep.at(0).name()
                  << "..." << std::endl;
    }
    this->mParams.decreaseLevel();
    /*------------------------------------------------------------------------*/
    timerPhase.restart();
    if (this->params().enabledAdolc())
    {
        if (this->computeGradients())
        {
            this->evalFovForwardTypeTwo(whichPhase, dfDep, dfGradDep,
                                        tapeIdLocalBorder, dfIndep, dfGradIndep,
                                        timeIter);
        }
        else
        {
            this->evalZosForwardTypeTwo(whichPhase, dfDep, tapeIdLocalBorder,
                                        dfIndep, timeIter);
        }
    }
    else
    {
        this->template setVectDfDepToDfIndep<CT>(dfDep, dfIndep, timeIter);
        this->evalTypeTwoFuncAtPartBorder(whichPhase, dfDep, dfIndep, timeIter,
                                          funcPhaseInner, funcPhaseBorder);
    }
    this->mWallClockTimes[whichPhase] += timerPhase.elapsed();
    /*------------------------------------------------------------------------*/
    timerSync.restart();
    for (std::size_t ii = 0; ii < dfDep.size(); ++ii)
    {
        dfDep[ii].copyValuesFromMemCommToSendBuffer(timeIter);
    }
    for (std::size_t ii = 0; ii < dfGradDep.size(); ++ii)
    {
        dfGradDep[ii].copyValuesFromMemCommToSendBuffer(timeIter);
    }
    for (std::size_t ii = 0; ii < dfDep.size(); ++ii)
    {
        dfDep[ii].exchangeValuesInBuffers(timeIter);
    }
    for (std::size_t ii = 0; ii < dfGradDep.size(); ++ii)
    {
        dfGradDep[ii].exchangeValuesInBuffers(timeIter);
    }
    this->mWallClockTimes["sync"] += timerSync.elapsed();
    /*------------------------------------------------------------------------*/
    timerPhase.restart();
    if (this->params().enabledAdolc())
    {
        if (this->computeGradients())
        {
            this->evalFovForwardTypeTwo(whichPhase, this->mVectUnknownDfsDomTmp,
                                        this->mVectGradUnknownDfsDomTmp,
                                        tapeIdLocalInner, dfIndep, dfGradIndep,
                                        timeIter);
        }
        else
        {
            this->evalZosForwardTypeTwo(whichPhase, this->mVectUnknownDfsDomTmp,
                                        tapeIdLocalInner, dfIndep, timeIter);
        }
        for (std::size_t ii = 0; ii < this->mVectUnknownDfsDomTmp.size(); ++ii)
        {
            dfDep[ii] += (this->mVectUnknownDfsDomTmp[ii] - dfIndep[ii]);
        }
        for (std::size_t ii = 0; ii < this->mVectGradUnknownDfsDomTmp.size();
             ++ii)
        {
            dfGradDep[ii] += this->mVectGradUnknownDfsDomTmp[ii];
            // Subtract identity matrix from dfGradDep [KF, 15.05.2014].
            //             // dfGradDep[ii] -= IdMat;
            //
        }
    }
    else
    {
        this->evalTypeTwoFuncAtPartInner(whichPhase, dfDep, dfIndep, timeIter,
                                         funcPhaseInner, funcPhaseBorder);
    }
    this->mWallClockTimes[whichPhase] += timerPhase.elapsed();
    /*------------------------------------------------------------------------*/
    timerSync.restart();
    for (std::size_t ii = 0; ii < dfDep.size(); ++ii)
    {
        dfDep[ii].waitAll();
    }
    for (std::size_t ii = 0; ii < dfGradDep.size(); ++ii)
    {
        dfGradDep[ii].waitAll();
    }
    for (std::size_t ii = 0; ii < dfDep.size(); ++ii)
    {
        dfDep[ii].copyValuesFromReceiveBufferToMemGhost(timeIter);
    }
    for (std::size_t ii = 0; ii < dfGradDep.size(); ++ii)
    {
        dfGradDep[ii].copyValuesFromReceiveBufferToMemGhost(timeIter);
    }
    this->mWallClockTimes["sync"] += timerSync.elapsed();

    /*------------------------------------------------------------------------*/
    this->mParams.increaseLevel();
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << "   Evaluated."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::evalTypeThreeAsynchronously(
    const std::string& whichPhase, CT& dfDep, std::vector<CT>& dfGradDep,
    const ScaFES_VectDf<CT>& dfIndep, const ScaFES_VectDf<CT>& dfGradIndep,
    const int& timeIter, const int& tapeIdLocalInner,
    const int& tapeIdLocalBorder,
    ScaFES_FuncTypeThreeToImplement<CT> funcPhaseInner,
    ScaFES_FuncTypeThreeToImplement<CT> funcPhaseBorder)
{
    ScaFES::Timer timerSync;
    ScaFES::Timer timerPhase;

    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Evaluate phase "
                  << whichPhase << " asynchronously: "
                  << " * " << whichPhase
                  << " asynchronously: " << dfIndep.at(0).name() << "..."
                  << std::endl;
    }
    this->mParams.decreaseLevel();
    /*------------------------------------------------------------------------*/
    dfDep = static_cast<CT>(0);
    /*------------------------------------------------------------------------*/
    timerPhase.restart();
    if (this->params().enabledAdolc())
    {
        if (this->computeGradients())
        {
            this->evalFovForwardTypeThree(whichPhase, dfDep, dfGradDep,
                                          tapeIdLocalBorder, dfIndep,
                                          dfGradIndep, timeIter);
        }
        else
        {
            this->evalZosForwardTypeThree(whichPhase, dfDep, tapeIdLocalBorder,
                                          dfIndep, timeIter);
        }
    }
    else
    {
        this->evalTypeThreeFuncAtPartBorder(whichPhase, dfDep, dfIndep,
                                            timeIter, funcPhaseInner,
                                            funcPhaseBorder);
    }
    this->mWallClockTimes[whichPhase] += timerPhase.elapsed();
    CT dfTmp = static_cast<CT>(0);
    std::vector<CT> dfGradTmp(dfGradDep.size());
    timerPhase.restart();
    if (this->params().enabledAdolc())
    {
        if (this->computeGradients())
        {
            this->evalFovForwardTypeThree(whichPhase, dfTmp, dfGradTmp,
                                          tapeIdLocalInner, dfIndep,
                                          dfGradIndep, timeIter);
        }
        else
        {
            this->evalZosForwardTypeThree(whichPhase, dfTmp, tapeIdLocalInner,
                                          dfIndep, timeIter);
        }
    }
    else
    {
        this->evalTypeThreeFuncAtPartInner(whichPhase, dfTmp, dfIndep, timeIter,
                                           funcPhaseInner, funcPhaseBorder);
    }
    this->mWallClockTimes[whichPhase] += timerPhase.elapsed();
    dfDep += dfTmp;
    for (std::size_t ii = 0; ii < dfGradDep.size(); ++ii)
    {
        dfGradDep[ii] += dfGradTmp[ii];
    }

    /*------------------------------------------------------------------------*/
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * (Local) value = " << dfDep << std::endl;
    }
    /*------------------------------------------------------------------------*/
    this->mParams.increaseLevel();
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << "   Done."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::evalTypeOneSynchronously(
    const std::string& whichPhase, ScaFES_VectDf<CT>& dfDepDomain,
    ScaFES_VectDf<CT>& dfDepGradDomain, ScaFES_VectDf<CT>& dfDepBdry,
    ScaFES_VectDf<CT>& dfDepGradBdry, const std::vector<CT>& dfIndep,
    const std::vector<CT>& dfGradIndep, const int& timeIter,
    const int& tapeIdLocalAll, ScaFES_FuncTypeOneToImplement<CT> funcPhaseInner,
    ScaFES_FuncTypeOneToImplement<CT> funcPhaseBorder)
{
    ScaFES::Timer timerSync;
    ScaFES::Timer timerPhase;

    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Evaluate phase "
                  << whichPhase
                  << " synchronously: " << dfDepDomain.at(0).name() << "..."
                  << std::endl;
    }
    this->mParams.decreaseLevel();

    /*------------------------------------------------------------------------*/
    ScaFES_VectDf<CT> vectAllTmp(dfDepDomain);
    for (std::size_t ii = 0; ii < dfDepBdry.size(); ++ii)
    {
        vectAllTmp.push_back(dfDepBdry[ii]);
    }
    ScaFES_VectDf<CT> vectGradAllTmp(dfDepGradDomain);
    for (std::size_t ii = 0; ii < dfDepGradBdry.size(); ++ii)
    {
        vectGradAllTmp.push_back(dfDepGradBdry[ii]);
    }
    timerPhase.restart();
    if (this->params().enabledAdolc())
    {
        if (this->computeGradients())
        {
            this->evalFovForwardTypeOne(whichPhase, vectAllTmp, vectGradAllTmp,
                                        tapeIdLocalAll, dfIndep, dfGradIndep,
                                        timeIter);
        }
        else
        {
            this->evalZosForwardTypeOne(whichPhase, vectAllTmp, tapeIdLocalAll,
                                        dfIndep, timeIter);
        }
    }
    else
    {
        this->evalTypeOneFuncAtPartAll(whichPhase, vectAllTmp, dfIndep,
                                       timeIter, funcPhaseInner,
                                       funcPhaseBorder);
    }
    this->mWallClockTimes[whichPhase] += timerPhase.elapsed();
    /*------------------------------------------------------------------------*/
    timerSync.restart();
    for (std::size_t ii = 0; ii < dfDepDomain.size(); ++ii)
    {
        dfDepDomain[ii].sync(timeIter);
    }
    for (std::size_t ii = 0; ii < dfDepBdry.size(); ++ii)
    {
        dfDepBdry[ii].sync(timeIter);
    }
    this->mWallClockTimes["sync"] += timerSync.elapsed();
    /*------------------------------------------------------------------------*/
    this->mParams.increaseLevel();
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << "   Evaluated."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::evalTypeTwoSynchronously(
    const std::string& whichPhase, ScaFES_VectDf<CT>& dfDepDomain,
    ScaFES_VectDf<CT>& dfDepGradDomain, ScaFES_VectDf<CT>& dfDepBdry,
    ScaFES_VectDf<CT>& dfDepGradBdry, const ScaFES_VectDf<CT>& dfIndepDomain,
    const ScaFES_VectDf<CT>& dfIndepGradDomain,
    const ScaFES_VectDf<CT>& dfIndepBdry,
    const ScaFES_VectDf<CT>& dfIndepGradBdry, const int& timeIter,
    const int& tapeIdLocalAll, ScaFES_FuncTypeTwoToImplement<CT> funcPhaseInner,
    ScaFES_FuncTypeTwoToImplement<CT> funcPhaseBorder)
{
    ScaFES::Timer timerSync;
    ScaFES::Timer timerPhase;

    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Evaluate phase "
                  << whichPhase
                  << " synchronously: " << dfDepDomain.at(0).name() << "..."
                  << std::endl;
    }
    this->mParams.decreaseLevel();

    /*------------------------------------------------------------------------*/
    // 1) Update DOMAIN data fields at all GLOBAL INNER nodes.
    // Remark: Only those data fields are involved which are defined at
    // the global inner nodes.
    // In detail: the dependent data fields are only dependent on data fields
    // defined on the global inner nodes, i.e. "domain data fields".
    // 2) Update ALL (=DOMAIN+BOUNDARY) data fields at all GLOBAL BORDER nodes.
    ScaFES_VectDf<CT> vectIndepAllTmp(dfIndepDomain);
    for (std::size_t ii = 0; ii < dfIndepBdry.size(); ++ii)
    {
        vectIndepAllTmp.push_back(dfIndepBdry[ii]);
    }
    ScaFES_VectDf<CT> vectIndepGradAllTmp(dfIndepGradDomain);
    for (std::size_t ii = 0; ii < dfIndepGradBdry.size(); ++ii)
    {
        vectIndepGradAllTmp.push_back(dfIndepGradBdry[ii]);
    }
    ScaFES_VectDf<CT> vectDepAllTmp(dfDepDomain);
    for (std::size_t ii = 0; ii < dfDepBdry.size(); ++ii)
    {
        vectDepAllTmp.push_back(dfDepBdry[ii]);
    }
    ScaFES_VectDf<CT> vectDepGradAllTmp(dfDepGradDomain);
    for (std::size_t ii = 0; ii < dfDepGradBdry.size(); ++ii)
    {
        vectDepGradAllTmp.push_back(dfDepGradBdry[ii]);
    }
    timerPhase.restart();
    if (this->params().enabledAdolc())
    {
        if (this->computeGradients())
        {
            this->evalFovForwardTypeTwo(
                whichPhase, vectDepAllTmp, vectDepGradAllTmp, tapeIdLocalAll,
                vectIndepAllTmp, vectIndepGradAllTmp, timeIter);
        }
        else
        {
            this->evalZosForwardTypeTwo(whichPhase, vectDepAllTmp,
                                        tapeIdLocalAll, vectIndepAllTmp,
                                        timeIter);
        }
    }
    else
    {
        this->template setVectDfDepToDfIndep<CT>(vectDepAllTmp, vectIndepAllTmp,
                                                 timeIter);
        this->evalTypeTwoFuncAtPartAll(whichPhase, vectDepAllTmp,
                                       vectIndepAllTmp, timeIter,
                                       funcPhaseInner, funcPhaseBorder);
    }
    this->mWallClockTimes[whichPhase] += timerPhase.elapsed();
    /*------------------------------------------------------------------------*/
    timerSync.restart();
    for (std::size_t ii = 0; ii < dfDepDomain.size(); ++ii)
    {
        dfDepDomain[ii].sync(timeIter);
    }
    for (std::size_t ii = 0; ii < dfDepGradDomain.size(); ++ii)
    {
        dfDepGradDomain[ii].sync(timeIter);
    }
    for (std::size_t ii = 0; ii < dfDepBdry.size(); ++ii)
    {
        dfDepBdry[ii].sync(timeIter);
    }
    for (std::size_t ii = 0; ii < dfDepGradBdry.size(); ++ii)
    {
        dfDepGradBdry[ii].sync(timeIter);
    }
    this->mWallClockTimes["sync"] += timerSync.elapsed();
    /*------------------------------------------------------------------------*/
    this->mParams.increaseLevel();
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << "   Done."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::evalTypeThreeSynchronously(
    const std::string& whichPhase, CT& dfDep, std::vector<CT>& dfGradDep,
    const ScaFES_VectDf<CT>& dfIndepDomain,
    const ScaFES_VectDf<CT>& dfGradIndepDomain,
    const ScaFES_VectDf<CT>& dfIndepBdry,
    const ScaFES_VectDf<CT>& dfGradIndepBdry, const int& timeIter,
    const int& tapeIdLocalAll,
    ScaFES_FuncTypeThreeToImplement<CT> funcPhaseInner,
    ScaFES_FuncTypeThreeToImplement<CT> funcPhaseBorder)
{
    ScaFES::Timer timerPhase;

    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Evaluate phase " << whichPhase
                  << " synchronously: " << dfIndepDomain.at(0).name() << "..."
                  << std::endl;
    }
    this->mParams.decreaseLevel();
    /*------------------------------------------------------------------------*/
    // Compute scalar dfDep for ALL (=DOMAIN+BOUNDARY) data fields
    // at all GLOBAL BORDER nodes.
    dfDep = static_cast<CT>(0);
    CT dfTmp = static_cast<CT>(0);
    std::vector<CT> dfGradTmp(dfGradDep.size(), static_cast<CT>(0));
    ScaFES_VectDf<CT> vectAllTmp(dfIndepDomain);
    for (std::size_t ii = 0; ii < dfIndepBdry.size(); ++ii)
    {
        vectAllTmp.push_back(dfIndepBdry[ii]);
    }
    ScaFES_VectDf<CT> vectGradAllTmp(dfGradIndepDomain);
    for (std::size_t ii = 0; ii < dfGradIndepBdry.size(); ++ii)
    {
        vectGradAllTmp.push_back(dfGradIndepBdry[ii]);
    }
    timerPhase.restart();
    if (this->params().enabledAdolc())
    {
        if (this->computeGradients())
        {
            this->evalFovForwardTypeThree(whichPhase, dfTmp, dfGradTmp,
                                          tapeIdLocalAll, vectAllTmp,
                                          vectGradAllTmp, timeIter);
        }
        else
        {
            this->evalZosForwardTypeThree(whichPhase, dfTmp, tapeIdLocalAll,
                                          vectAllTmp, timeIter);
        }
    }
    else
    {
        this->evalTypeThreeFuncAtPartAll(whichPhase, dfTmp, vectAllTmp,
                                         timeIter, funcPhaseInner,
                                         funcPhaseBorder);
    }
    this->mWallClockTimes[whichPhase] += timerPhase.elapsed();
    dfDep += dfTmp;
    for (std::size_t ii = 0; ii < dfGradTmp.size(); ++ii)
    {
        dfGradDep[ii] += dfGradTmp[ii];
    }
    /*------------------------------------------------------------------------*/
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Local value = " << dfDep << std::endl;
    }
    /*------------------------------------------------------------------------*/
    this->mParams.increaseLevel();
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << "   Done."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
template <typename TT>
inline void Problem<OWNPRBLM, CT, DIM>::evalTypeOneFuncAtPartInner(
    const std::string& whichPhase, ScaFES_VectDf<TT>& dfDep,
    const std::vector<TT>& dfIndep, const int& timeIter,
    ScaFES_FuncTypeOneToImplement<TT> funcPhaseInner,
    ScaFES_FuncTypeOneToImplement<TT> funcPhaseBorder)
{
    this->template evalTypeOneAt<TT>(whichPhase, "loc inner AND glob inner",
                                     funcPhaseInner, dfDep, dfIndep,
                                     this->mNodesLocInnerGlobInner, timeIter);
    this->template evalTypeOneAt<TT>(whichPhase, "loc inner AND glob border",
                                     funcPhaseBorder, dfDep, dfIndep,
                                     this->mNodesLocInnerGlobBorder, timeIter);
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
template <typename TT>
inline void Problem<OWNPRBLM, CT, DIM>::evalTypeTwoFuncAtPartInner(
    const std::string& whichPhase, ScaFES_VectDf<TT>& dfDep,
    const ScaFES_VectDf<TT>& dfIndep, const int& timeIter,
    ScaFES_FuncTypeTwoToImplement<TT> funcPhaseInner,
    ScaFES_FuncTypeTwoToImplement<TT> funcPhaseBorder)
{
    this->template evalTypeTwoAt<TT>(whichPhase, "loc inner",
                                     funcPhaseInner, funcPhaseBorder, dfDep, dfIndep,
                                     this->mNodesLocInner, timeIter);
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
template <typename TT>
inline void Problem<OWNPRBLM, CT, DIM>::evalTypeThreeFuncAtPartInner(
    const std::string& whichPhase, TT& dfDep, const ScaFES_VectDf<TT>& dfIndep,
    const int& timeIter, ScaFES_FuncTypeThreeToImplement<TT> funcPhaseInner,
    ScaFES_FuncTypeThreeToImplement<TT> funcPhaseBorder)
{
    this->template evalTypeThreeAt<TT>(whichPhase, "loc inner AND glob inner",
                                       funcPhaseInner, dfDep, dfIndep,
                                       this->mNodesLocInnerGlobInner, timeIter);
    TT dfDepTmp = static_cast<TT>(0);
    this->template evalTypeThreeAt<TT>(
        whichPhase, "loc inner AND glob border", funcPhaseBorder, dfDepTmp,
        dfIndep, this->mNodesLocInnerGlobBorder, timeIter);
    dfDep = dfDep + dfDepTmp;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
template <typename TT>
inline void Problem<OWNPRBLM, CT, DIM>::evalTypeOneFuncAtPartBorder(
    const std::string& whichPhase, ScaFES_VectDf<TT>& dfDep,
    const std::vector<TT>& dfIndep, const int& timeIter,
    ScaFES_FuncTypeOneToImplement<TT> funcPhaseInner,
    ScaFES_FuncTypeOneToImplement<TT> funcPhaseBorder)
{
    this->template evalTypeOneAt<TT>(whichPhase, "loc border AND glob inner",
                                     funcPhaseInner, dfDep, dfIndep,
                                     this->mNodesLocBorderGlobInner, timeIter);
    this->template evalTypeOneAt<TT>(whichPhase, "loc border AND glob border",
                                     funcPhaseBorder, dfDep, dfIndep,
                                     this->mNodesLocBorderGlobBorder, timeIter);
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
template <typename TT>
inline void Problem<OWNPRBLM, CT, DIM>::evalTypeTwoFuncAtPartBorder(
    const std::string& whichPhase, ScaFES_VectDf<TT>& dfDep,
    const ScaFES_VectDf<TT>& dfIndep, const int& timeIter,
    ScaFES_FuncTypeTwoToImplement<TT> funcPhaseInner,
    ScaFES_FuncTypeTwoToImplement<TT> funcPhaseBorder)
{
    this->template evalTypeTwoAt<TT>(whichPhase, "loc border",
                                     funcPhaseInner, funcPhaseBorder,
                                     dfDep, dfIndep,
                                     this->mNodesLocBorder, timeIter);
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
template <typename TT>
inline void Problem<OWNPRBLM, CT, DIM>::evalTypeThreeFuncAtPartBorder(
    const std::string& whichPhase, TT& dfDep, const ScaFES_VectDf<TT>& dfIndep,
    const int& timeIter, ScaFES_FuncTypeThreeToImplement<TT> funcPhaseInner,
    ScaFES_FuncTypeThreeToImplement<TT> funcPhaseBorder)
{
    this->template evalTypeThreeAt<TT>(
        whichPhase, "loc border AND glob inner", funcPhaseInner, dfDep, dfIndep,
        this->mNodesLocBorderGlobInner, timeIter);
    TT dfDepTmp = static_cast<TT>(0);
    this->template evalTypeThreeAt<TT>(
        whichPhase, "loc border AND glob border", funcPhaseBorder, dfDepTmp,
        dfIndep, this->mNodesLocBorderGlobBorder, timeIter);
    dfDep = dfDep + dfDepTmp;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
template <typename TT>
inline void Problem<OWNPRBLM, CT, DIM>::evalTypeOneFuncAtPartAll(
    const std::string& whichPhase, ScaFES_VectDf<TT>& dfDep,
    const std::vector<TT>& dfIndep, const int& timeIter,
    ScaFES_FuncTypeOneToImplement<TT> funcPhaseInner,
    ScaFES_FuncTypeOneToImplement<TT> funcPhaseBorder)
{
    this->template evalTypeOneFuncAtPartInner<TT>(
        whichPhase, dfDep, dfIndep, timeIter, funcPhaseInner, funcPhaseBorder);
    this->template evalTypeOneFuncAtPartBorder<TT>(
        whichPhase, dfDep, dfIndep, timeIter, funcPhaseInner, funcPhaseBorder);
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
template <typename TT>
inline void Problem<OWNPRBLM, CT, DIM>::evalTypeTwoFuncAtPartAll(
    const std::string& whichPhase, ScaFES_VectDf<TT>& dfDep,
    const ScaFES_VectDf<TT>& dfIndep, const int& timeIter,
    ScaFES_FuncTypeTwoToImplement<TT> funcPhaseInner,
    ScaFES_FuncTypeTwoToImplement<TT> funcPhaseBorder)
{
    this->template evalTypeTwoAt<TT>(whichPhase, "local",
                                     funcPhaseInner, funcPhaseBorder,
                                     dfDep, dfIndep,
                                     this->mNodes, timeIter);
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
template <typename TT>
inline void Problem<OWNPRBLM, CT, DIM>::evalTypeThreeFuncAtPartAll(
    const std::string& whichPhase, TT& dfDep, const ScaFES_VectDf<TT>& dfIndep,
    const int& timeIter, ScaFES_FuncTypeThreeToImplement<TT> funcPhaseInner,
    ScaFES_FuncTypeThreeToImplement<TT> funcPhaseBorder)
{
    this->template evalTypeThreeFuncAtPartInner<TT>(
        whichPhase, dfDep, dfIndep, timeIter, funcPhaseInner, funcPhaseBorder);
    TT dfDepTmp = static_cast<TT>(0);
    this->template evalTypeThreeFuncAtPartBorder<TT>(
        whichPhase, dfDepTmp, dfIndep, timeIter, funcPhaseInner,
        funcPhaseBorder);
    dfDep = dfDep + dfDepTmp;
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
template <typename TT>
inline void Problem<OWNPRBLM, CT, DIM>::evalTypeOneAt(
    const std::string& whichPhase, const std::string& whereToEval,
    ScaFES_FuncTypeOneToImplement<TT> funcPhase, ScaFES_VectDf<TT>& dfDep,
    const std::vector<TT>& dfIndep, std::vector<int> setGridNodes,
    const int& timeIter)
{
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * " << whichPhase
                  << " at all " << whereToEval
                  << " nodes: " << dfDep.at(0).name() << ",...." << std::endl;
    }
    this->mParams.decreaseLevel();

    // Set t = t_S.
    // The rhs need not to be an independent variable the time at the
    // initialization will always be t=t_S.
    // ==> Variable of type double for rhs is sufficient.
    for (std::size_t ii = 0; ii < dfDep.size(); ++ii)
    {
        dfDep[ii].time() = this->params().timeIntervalStart();
    }

// * ii Iteration variable ---> private
// * idxNode WRITE access Corresponding global node number ---> private
// * comp WRITE access Corresponding component at idxNode ---> private
// * dfIndep READ access --> shared
// * dfDep WRITE access, Sum of all grid node portions. ---> shared
//   but must be protected using "critical".
//   "reduction" is not possible as type "adouble" is not suported.
// * timeIter READ access --> shared
#ifdef _OPENMP
#pragma omp parallel for schedule(static) \
    shared(dfDep, dfIndep, timeIter, funcPhase, setGridNodes)
#endif
    for (std::size_t ii = 0; ii < setGridNodes.size(); ++ii)
    {
        int comp = 0;
        ScaFES_IntNtuple idxNode;
        int idxNodeScalar = setGridNodes[ii]; // = memoryPos for kind field.
        this->mKind.memoryPos2IdxNode(idxNode, comp, idxNodeScalar);
        //        for (std::size_t jj = 0; jj < dfDep.size(); ++jj) {
        //            dfDep[jj](idxNode) = this->defaultValue().at(jj);
        //        }
        (static_cast<OWNPRBLM*>(this)->*funcPhase)(dfDep, dfIndep, idxNode,
                                                   timeIter);
    }

    this->mParams.increaseLevel();
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << "   Done."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
template <typename TT>
inline void Problem<OWNPRBLM, CT, DIM>::evalTypeTwoAt(
    const std::string& whichPhase, const std::string& whereToEval,
    ScaFES_FuncTypeTwoToImplement<TT> funcPhaseInner,
    ScaFES_FuncTypeTwoToImplement<TT> funcPhaseBorder,
    ScaFES_VectDf<TT>& dfDep,
    const ScaFES_VectDf<TT>& dfIndep, std::vector<int> setGridNodes,
    const int& timeIter)
{
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                                       << " * " << whichPhase
                  << " at all " << whereToEval
                  << " nodes: " << dfIndep.at(0).name() << ",...." << std::endl;
    }
    this->mParams.decreaseLevel();

    // Set t = t_S + timeIter * \tau = t_l + \tau.
    // Attention: Update time of dep. var. BEFORE calls of update methods.
    for (std::size_t ii = 0; ii < dfIndep.size(); ++ii)
    {
        // \tau is a constant. ==> \tau need not to be of type adouble. ==> OK.
        dfDep[ii].time() = dfIndep[ii].time() + this->params().tau();
    }

// for (std::size_t jj = 0; jj < dfDep.size(); ++jj) {
//    for (std::size_t ii = 0; ii < setGridNodes.size(); ++ii) {
//        int comp = 0;
//        ScaFES_IntNtuple idxNode;
//        int memPos = setGridNodes[ii];
//        this->mKind.memoryPos2IdxNode(idxNode, comp, memPos);
//        dfDep[jj](idxNode) = dfIndep[jj](idxNode);
//    }
//}

// * ii Iteration variable ---> private
// * idxNode WRITE access Corresponding global node number ---> private
// * comp WRITE access Corresponding component at idxNode ---> private
// * dfIndep READ access --> shared
// * dfDep WRITE access, Sum of all grid node portions. ---> shared
//   but must be protected using "critical".
//   "reduction" is not possible as type "adouble" is not suported.
// * timeIter READ access --> shared
#ifdef _OPENMP
#pragma omp parallel for schedule(static) \
    shared(dfDep, dfIndep, timeIter, funcPhaseInner, funcPhaseBorder, setGridNodes)
#endif
    for (std::size_t ii = 0; ii < setGridNodes.size(); ++ii)
    {
        int comp = 0;
        ScaFES_IntNtuple idxNode;
        int idxNodeScalar = setGridNodes[ii]; // = memoryPos for kind field.
        this->mKind.memoryPos2IdxNode(idxNode, comp, idxNodeScalar);
        // for (std::size_t jj = 0; jj < dfDep.size(); ++jj) {
        //    dfDep[jj](idxNode) = dfIndep[jj](idxNode);
        //}
        if (this->mKind(idxNodeScalar) & this->globalGrid().SCAFES_GLOBAL_INNER_NODE)
        {
            (static_cast<OWNPRBLM*>(this)->*funcPhaseInner)(dfDep, dfIndep, idxNode,
                                                            timeIter);
        }
        else if (this->mKind(idxNodeScalar) & this->globalGrid().SCAFES_GLOBAL_BORDER_NODE)
        {
            (static_cast<OWNPRBLM*>(this)->*funcPhaseBorder)(dfDep, dfIndep, idxNode,
                                                            timeIter);
        }
    }

    this->mParams.increaseLevel();
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << "   Done."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
template <typename TT>
inline void Problem<OWNPRBLM, CT, DIM>::evalTypeThreeAt(
    const std::string& whichPhase, const std::string& whereToEval,
    ScaFES_FuncTypeThreeToImplement<TT> funcPhase, TT& dfDep,
    const ScaFES_VectDf<TT>& dfIndep, std::vector<int> setGridNodes,
    const int& timeIter)
{
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * " << whichPhase
                  << " at all " << whereToEval
                  << " nodes: " << dfIndep.at(0).name() << ",...." << std::endl;
    }
    this->mParams.decreaseLevel();

    // Reset value of dependent variable.
    dfDep = static_cast<TT>(0);

// * ii Iteration variable ---> private
// * idxNode WRITE access Corresponding global node number ---> private
// * comp WRITE access Corresponding component at idxNode ---> private
// * dfIndep READ access --> shared
// * dfDep WRITE access, Sum of all grid node portions. ---> shared
//   but must be protected using "critical".
//   "reduction" is not possible as type "adouble" is not suported.
// * timeIter READ access --> shared

#ifdef _OPENMP
#pragma omp parallel for schedule(static) \
    shared(dfDep, dfIndep, timeIter, funcPhase, setGridNodes)
#endif
    for (std::size_t ii = 0; ii < setGridNodes.size(); ++ii)
    {
        int comp = 0;
        ScaFES_IntNtuple idxNode;
        int idxNodeScalar = setGridNodes[ii]; // = memoryPos for kind field.
        this->mKind.memoryPos2IdxNode(idxNode, comp, idxNodeScalar);
        TT dfDepTmp = static_cast<TT>(0);
        (static_cast<OWNPRBLM*>(this)->*funcPhase)(dfDepTmp, dfIndep, idxNode,
                                                   timeIter);
#ifdef _OPENMP
#pragma omp critical(sumOfdfDep)
        {
#endif
            dfDep = dfDep + dfDepTmp;
#ifdef _OPENMP
        }
#endif
    }

    this->mParams.increaseLevel();
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << "   Done."
                  << std::endl;
    }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::traceInitPbPhase(const int& timeIter)
{
    if (this->isTracedInit())
    {
        return;
    }

    std::string whichPhase = "init";

    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Trace "
                  << whichPhase << " phase..." << std::endl;
    }
    this->mParams.decreaseLevel();

#ifdef SCAFES_HAVE_ADOLC
    this->template traceTypeOne<adouble>(
        whichPhase, this->mVectUnknownDfsDomNew, this->mVectUnknownDfsBdryNew,
        this->vectParamsCurr(), timeIter, this->tapeIdInitPartitionInner(),
        this->tapeIdInitPartitionBorder(), this->tapeIdInitPartitionAll(),
        &OWNPRBLM::template initInner<adouble>,
        &OWNPRBLM::template initBorder<adouble>);
    this->mIsTracedInit = true;
#else
    std::ignore = timeIter;
#endif

    this->mParams.increaseLevel();
    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << "   Traced."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::traceInitGbPhase(const int& timeIter)
{
    if (this->isTracedInit())
    {
        return;
    }

    std::string whichPhase = "init";

    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Trace "
                  << whichPhase << " phase..." << std::endl;
    }
    this->mParams.decreaseLevel();

#ifdef SCAFES_HAVE_ADOLC
    this->template traceTypeTwo<adouble>(
        whichPhase, this->mVectUnknownDfsDomNew, this->mVectUnknownDfsBdryNew,
        this->vectParamsCurr(), timeIter, this->tapeIdInitPartitionInner(),
        this->tapeIdInitPartitionBorder(), this->tapeIdInitPartitionAll(),
        &OWNPRBLM::template initInner<adouble>,
        &OWNPRBLM::template initBorder<adouble>);
    this->mIsTracedInit = true;
#else
    std::ignore = timeIter;
#endif

    this->mParams.increaseLevel();
    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << "   Traced."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::traceUpdatePhase(const int& timeIter)
{
    if (this->isTracedUpdatePhase())
    {
        return;
    }

    std::string whichPhase = "update";

    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Trace "
                  << whichPhase << " phase..." << std::endl;
    }
    this->mParams.decreaseLevel();

#ifdef SCAFES_HAVE_ADOLC
    this->template traceTypeTwo<adouble>(
        whichPhase, this->mVectUnknownDfsDomNew, this->mVectUnknownDfsBdryNew,
        this->vectUnknownDfsDomOld(), this->mVectUnknownDfsBdryOld, timeIter,
        this->tapeIdUpdatePartitionInner(), this->tapeIdUpdatePartitionBorder(),
        this->tapeIdUpdatePartitionAll(),
        &OWNPRBLM::template updateInner<adouble>,
        &OWNPRBLM::template updateBorder<adouble>);
    this->mIsTracedUpdatePhase = true;
#else
    std::ignore = timeIter;
#endif

    this->mParams.increaseLevel();
    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                                        << "   Traced."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::traceUpdate2Phase(const int& timeIter)
{
    if (this->isTracedUpdate2Phase())
    {
        return;
    }

    std::string whichPhase = "update2";

    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Trace "
                  << whichPhase << " phase..." << std::endl;
    }
    this->mParams.decreaseLevel();

#ifdef SCAFES_HAVE_ADOLC
    this->template traceTypeTwo<adouble>(
        whichPhase, this->mVectUnknownDfsDomNew, this->mVectUnknownDfsBdryNew,
        this->vectUnknownDfsDomOld(), this->mVectUnknownDfsBdryOld, timeIter,
        this->tapeIdUpdate2PartitionInner(),
        this->tapeIdUpdate2PartitionBorder(), this->tapeIdUpdate2PartitionAll(),
        &OWNPRBLM::template updateInner2<adouble>,
        &OWNPRBLM::template updateBorder2<adouble>);
    this->mIsTracedUpdate2Phase = true;
#else
    std::ignore = timeIter;
#endif

    this->mParams.increaseLevel();
    if (this->params().rankOutput() == this->myRank() &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << "   Traced."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
template <typename TT>
inline void Problem<OWNPRBLM, CT, DIM>::traceTypeOne(
    const std::string& whichPhase, ScaFES_VectDf<CT>& dfDepDomain,
    ScaFES_VectDf<CT>& /*dfDepBdry*/, const std::vector<CT>& dfIndep,
    const int& timeIter, const int& tapeIdLocalInner,
    const int& tapeIdLocalBorder, const int& tapeIdLocalAll,
    ScaFES_FuncTypeOneToImplement<TT> funcPhaseInner,
    ScaFES_FuncTypeOneToImplement<TT> funcPhaseBorder)
{
#ifdef SCAFES_HAVE_ADOLC
    TT tmp;
    static_assert(std::is_same<decltype(tmp), adouble>::value,
                  "Only for TT=adouble!");

    if (!(this->params().enabledAdolc()))
    {
        throw std::runtime_error("Tracing not possible. ADOL-C not enabled.");
    }

    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Trace "
                  << whichPhase << ": " << dfDepDomain.at(0).name() << ",...)."
                  << std::endl;
    }
    this->mParams.decreaseLevel();

    ScaFES::Timer timerTrace;
    timerTrace.restart();

// REMARK: Currently, tracing can only be done using 1 OpenMP thread.
#ifdef _OPENMP
    const int NTHREADS_ENV = omp_get_num_threads();
    omp_set_num_threads(1);
#endif

    /*-------------------------------------------------&-----------------------*/
    std::vector<TT> aDfIndep;
    aDfIndep.reserve(dfIndep.size());
    for (std::size_t ii = 0; ii < dfIndep.size(); ++ii)
    {
        aDfIndep.push_back(dfIndep.at(ii));
    }

    /*------------------------------------------------------------------------*/
    unsigned long int nElemsDfIndepTotal = 1;
    for (std::size_t ii = 0; ii < dfDepDomain.size(); ++ii)
    {
        nElemsDfIndepTotal += dfDepDomain.at(ii).nElemsAllocated();
    }
    std::vector<TT> memoryAdfDep;
    if (0 < nElemsDfIndepTotal)
    {
        memoryAdfDep.resize(nElemsDfIndepTotal);
    }
    ScaFES_VectDf<TT> aDfDep;
    this->template copyDataFieldStructure<TT>(aDfDep, dfDepDomain,
                                              memoryAdfDep);

/*------------------------------------------------------------------------*/
#ifdef _OPENMP
    int idxThread = omp_get_thread_num();
#else
    int idxThread = 0;
#endif
    if (this->useAsynchronMode())
    {
        trace_on(tapeIdLocalBorder + idxThread, 1);
        this->template setIndependentVector<TT>(aDfIndep, dfIndep);
        this->template setVectDfDepToDefValue<TT>(aDfDep, this->defaultValue(),
                                                  timeIter);
        this->template evalTypeOneFuncAtPartBorder<TT>(
            whichPhase, aDfDep, aDfIndep, timeIter, funcPhaseInner,
            funcPhaseBorder);
        this->template setDependentVectDf<TT>(aDfDep, dfDepDomain);
        trace_off();

        trace_on(tapeIdLocalInner + idxThread, 1);
        this->template setIndependentVector<TT>(aDfIndep, dfIndep);
        this->template setVectDfDepToDefValue<TT>(aDfDep, this->defaultValue(),
                                                  timeIter);
        this->template evalTypeOneFuncAtPartInner<TT>(
            whichPhase, aDfDep, aDfIndep, timeIter, funcPhaseInner,
            funcPhaseBorder);
        this->template setDependentVectDf<TT>(aDfDep, dfDepDomain);
        trace_off();
    }
    else
    {
        trace_on(tapeIdLocalAll + idxThread, 1);
        this->template setIndependentVector<TT>(aDfIndep, dfIndep);
        this->template setVectDfDepToDefValue<TT>(aDfDep, this->defaultValue(),
                                                  timeIter);
        this->template evalTypeOneFuncAtPartAll<TT>(
            whichPhase, aDfDep, aDfIndep, timeIter, funcPhaseInner,
            funcPhaseBorder);
        this->template setDependentVectDf<TT>(aDfDep, dfDepDomain);
        trace_off();
    }
#ifdef _OPENMP
    omp_set_num_threads(NTHREADS_ENV);
#endif

    this->mWallClockTimes["trace"] += timerTrace.elapsed();
#else
    std::ignore = whichPhase;
    std::ignore = dfIndep;
    std::ignore = dfDepDomain;
    std::ignore = timeIter;
    std::ignore = tapeIdLocalInner;
    std::ignore = tapeIdLocalBorder;
    std::ignore = tapeIdLocalAll;
    std::ignore = funcPhaseInner;
    std::ignore = funcPhaseBorder;
#endif

    // REMARK: No synchronization is necessary as one just wants to trace
    // the mathematical operations.
    this->mParams.increaseLevel();
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << "   Traced."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
template <typename TT>
inline void Problem<OWNPRBLM, CT, DIM>::traceTypeTwo(
    const std::string& whichPhase, ScaFES_VectDf<CT>& dfDepDomain,
    ScaFES_VectDf<CT>& /*dfDepBdry*/, const ScaFES_VectDf<CT>& dfIndepDomain,
    const ScaFES_VectDf<CT>& /*dfIndepBdry*/, const int& timeIter,
    const int& tapeIdLocalInner, const int& tapeIdLocalBorder,
    const int& tapeIdLocalAll, ScaFES_FuncTypeTwoToImplement<TT> funcPhaseInner,
    ScaFES_FuncTypeTwoToImplement<TT> funcPhaseBorder)
{
#ifdef SCAFES_HAVE_ADOLC
    TT tmp;
    static_assert(std::is_same<decltype(tmp), adouble>::value,
                  "Only for TT=adouble!");

    if (!(this->params().enabledAdolc()))
    {
        throw std::runtime_error("Tracing not possible. ADOL-C not enabled.");
    }

    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Trace "
                  << whichPhase << ": " << dfIndepDomain.at(0).name()
                  << ",...)." << std::endl;
    }
    this->mParams.decreaseLevel();

    ScaFES::Timer timerTrace;
    timerTrace.restart();

// REMARK: Currently, tracing can only be done using 1 OpenMP thread.
#ifdef _OPENMP
    const int NTHREADS_ENV = omp_get_num_threads();
    omp_set_num_threads(1);
#endif

    /*------------------------------------------------------------------------*/
    unsigned long int nElemsDfIndepTotal = 1;
    for (std::size_t ii = 0; ii < dfIndepDomain.size(); ++ii)
    {
        nElemsDfIndepTotal += dfIndepDomain.at(ii).nElemsAllocated();
    }
    std::vector<TT> memoryAdfIndep;
    if (0 < nElemsDfIndepTotal)
    {
        memoryAdfIndep.resize(nElemsDfIndepTotal);
    }
    ScaFES_VectDf<TT> aDfIndep;
    this->template copyDataFieldStructure<TT>(aDfIndep, dfIndepDomain,
                                              memoryAdfIndep);
    std::vector<TT> memoryAdfDep;
    if (0 < nElemsDfIndepTotal)
    {
        memoryAdfDep.resize(nElemsDfIndepTotal);
    }
    ScaFES_VectDf<TT> aDfDep;
    this->template copyDataFieldStructure<TT>(aDfDep, dfDepDomain,
                                              memoryAdfDep);

/*------------------------------------------------------------------------*/
#ifdef _OPENMP
    int idxThread = omp_get_thread_num();
#else
    int idxThread = 0;
#endif
    if (this->useAsynchronMode())
    {
        trace_on(tapeIdLocalBorder + idxThread, 1);
        this->template setIndependentVectDf<TT>(aDfIndep, dfIndepDomain);
        this->template setVectDfDepToDfIndep<TT>(aDfDep, aDfIndep, timeIter);
        this->template evalTypeTwoFuncAtPartBorder<TT>(
            whichPhase, aDfDep, aDfIndep, timeIter, funcPhaseInner,
            funcPhaseBorder);
        this->template setDependentVectDf<TT>(aDfDep, dfDepDomain);
        trace_off();

        trace_on(tapeIdLocalInner + idxThread, 1);
        this->template setIndependentVectDf<TT>(aDfIndep, dfIndepDomain);
        this->template setVectDfDepToDfIndep<TT>(aDfDep, aDfIndep, timeIter);
        this->template evalTypeTwoFuncAtPartInner<TT>(
            whichPhase, aDfDep, aDfIndep, timeIter, funcPhaseInner,
            funcPhaseBorder);
        this->template setDependentVectDf<TT>(aDfDep, dfDepDomain);
        trace_off();
    }
    else
    {
        trace_on(tapeIdLocalAll + idxThread, 1);
        this->template setIndependentVectDf<TT>(aDfIndep, dfIndepDomain);
        this->template setVectDfDepToDfIndep<TT>(aDfDep, aDfIndep, timeIter);
        this->template evalTypeTwoFuncAtPartAll<TT>(
            whichPhase, aDfDep, aDfIndep, timeIter, funcPhaseInner,
            funcPhaseBorder);
        this->template setDependentVectDf<TT>(aDfDep, dfDepDomain);
        trace_off();
    }
#ifdef _OPENMP
    omp_set_num_threads(NTHREADS_ENV);
#endif

    this->mWallClockTimes["trace"] += timerTrace.elapsed();
#else
    std::ignore = whichPhase;
    std::ignore = dfDepDomain;
    std::ignore = dfIndepDomain;
    std::ignore = tapeIdLocalInner;
    std::ignore = tapeIdLocalBorder;
    std::ignore = tapeIdLocalAll;
    std::ignore = timeIter;
    std::ignore = funcPhaseInner;
    std::ignore = funcPhaseBorder;
#endif

    // REMARK: No synchronization is necessary as one just wants to trace
    // the mathematical operations.
    this->mParams.increaseLevel();
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << "   Traced."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
template <typename TT>
inline void Problem<OWNPRBLM, CT, DIM>::traceTypeThree(
    const std::string& whichPhase, CT& dfDep,
    const ScaFES_VectDf<CT>& dfIndepDomain,
    const ScaFES_VectDf<CT>& /*dfIndepBdry*/, const int& timeIter,
    const int& tapeIdLocalInner, const int& tapeIdLocalBorder,
    const int& tapeIdLocalAll,
    ScaFES_FuncTypeThreeToImplement<TT> funcPhaseInner,
    ScaFES_FuncTypeThreeToImplement<TT> funcPhaseBorder)
{
#ifdef SCAFES_HAVE_ADOLC
    TT tmp;
    static_assert(std::is_same<decltype(tmp), adouble>::value,
                  "Only for TT=adouble!");

    if (!(this->params().enabledAdolc()))
    {
        throw std::runtime_error("Tracing not possible. ADOL-C not enabled.");
    }

    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Trace "
                  << whichPhase << ": " << dfIndepDomain.at(0).name()
                  << ",...)." << std::endl;
    }
    this->mParams.decreaseLevel();

    ScaFES::Timer timerTrace;
    timerTrace.restart();

// REMARK: Currently, tracing can only done using 1 OpenMP thread.
#ifdef _OPENMP
    const int NTHREADS_ENV = omp_get_num_threads();
    omp_set_num_threads(1);
#endif

    /*------------------------------------------------------------------------*/
    unsigned long int nElemsDfIndepTotal = 1;
    for (std::size_t ii = 0; ii < dfIndepDomain.size(); ++ii)
    {
        nElemsDfIndepTotal += dfIndepDomain.at(ii).nElemsAllocated();
    }
    std::vector<TT> memoryAdfIndep;
    if (0 < nElemsDfIndepTotal)
    {
        memoryAdfIndep.resize(nElemsDfIndepTotal);
    }
    ScaFES_VectDf<TT> aDfIndep;
    this->template copyDataFieldStructure<TT>(aDfIndep, dfIndepDomain,
                                              memoryAdfIndep);

    /*------------------------------------------------------------------------*/
    dfDep = static_cast<CT>(0);
    TT aDfDep = static_cast<TT>(dfDep);

/*------------------------------------------------------------------------*/
#ifdef _OPENMP
    int idxThread = omp_get_thread_num();
#else
    int idxThread = 0;
#endif

    if (this->useAsynchronMode())
    {
        trace_on(tapeIdLocalBorder + idxThread, 1);
        this->template setIndependentVectDf<TT>(aDfIndep, dfIndepDomain);
        this->template evalTypeThreeFuncAtPartBorder<adouble>(
            whichPhase, aDfDep, aDfIndep, timeIter, funcPhaseInner,
            funcPhaseBorder);
        this->template setDependentScalar<TT>(aDfDep, dfDep);
        trace_off();

        trace_on(tapeIdLocalInner + idxThread, 1);
        this->template setIndependentVectDf<TT>(aDfIndep, dfIndepDomain);
        this->template evalTypeThreeFuncAtPartInner<adouble>(
            whichPhase, aDfDep, aDfIndep, timeIter, funcPhaseInner,
            funcPhaseBorder);
        this->template setDependentScalar<TT>(aDfDep, dfDep);
        trace_off();
    }
    else
    {
        trace_on(tapeIdLocalAll + idxThread, 1);
        this->template setIndependentVectDf<TT>(aDfIndep, dfIndepDomain);
        this->template evalTypeThreeFuncAtPartAll<adouble>(
            whichPhase, aDfDep, aDfIndep, timeIter, funcPhaseInner,
            funcPhaseBorder);
        this->template setDependentScalar<TT>(aDfDep, dfDep);
        trace_off();
    }

#ifdef _OPENMP
    omp_set_num_threads(NTHREADS_ENV);
#endif

    this->mWallClockTimes["trace"] += timerTrace.elapsed();
#else
    std::ignore = whichPhase;
    std::ignore = dfDep;
    std::ignore = dfIndepDomain;
    std::ignore = tapeIdLocalInner;
    std::ignore = tapeIdLocalBorder;
    std::ignore = tapeIdLocalAll;
    std::ignore = timeIter;
    std::ignore = funcPhaseInner;
    std::ignore = funcPhaseBorder;
#endif

    // REMARK: No synchronization is necessary as one just wants to trace
    // the mathematical operations.
    this->mParams.increaseLevel();
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << "   Traced."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
template <typename TT>
inline void Problem<OWNPRBLM, CT, DIM>::setIndependentScalar(TT& aDf,
                                                             const CT& df)
{
#ifdef SCAFES_HAVE_ADOLC
    TT tmp;
    static_assert(std::is_same<decltype(tmp), adouble>::value,
                  "Only for TT=adouble!");
    aDf <<= df;
#else
    std::ignore = aDf;
    std::ignore = df;
#endif
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
template <typename TT>
inline void Problem<OWNPRBLM, CT, DIM>::setDependentScalar(TT& aDf, CT& df)
{
#ifdef SCAFES_HAVE_ADOLC
    TT tmp;
    static_assert(std::is_same<decltype(tmp), adouble>::value,
                  "Only for TT=adouble!");
    aDf >>= df;
#else
    std::ignore = aDf;
    std::ignore = df;
#endif
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
template <typename TT>
inline void Problem<OWNPRBLM, CT, DIM>::setIndependentDf(
    ScaFES::DataField<TT, DIM>& aDf, const ScaFES::DataField<CT, DIM>& df)
{
    int comp = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) shared(df, aDf)
#endif
    for (int jj = 0; jj < df.memAll().nNodesTotal(); ++jj)
    {
        ScaFES_IntNtuple idxNode;
        df.memoryPos2IdxNode(idxNode, comp, jj);
        this->template setIndependentScalar<TT>(aDf(idxNode), df(idxNode));
    }
    this->template setIndependentScalar<TT>(aDf.time(), df.time());
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
template <typename TT>
inline void
Problem<OWNPRBLM, CT, DIM>::setDependentDf(ScaFES::DataField<TT, DIM>& aDf,
                                           ScaFES::DataField<CT, DIM>& df)
{
    int comp = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) shared(df, aDf)
#endif
    for (int jj = 0; jj < df.memAll().nNodesTotal(); ++jj)
    {
        ScaFES_IntNtuple idxNode;
        df.memoryPos2IdxNode(idxNode, comp, jj);
        this->template setDependentScalar<TT>(aDf(idxNode), df(idxNode));
    }
    this->template setDependentScalar<TT>(aDf.time(), df.time());
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
template <typename TT>
inline void
Problem<OWNPRBLM, CT, DIM>::setIndependentVectDf(ScaFES_VectDf<TT>& aDf,
                                                 const ScaFES_VectDf<CT>& df)
{
    for (std::size_t ii = 0; ii < df.size(); ++ii)
    {
        this->template setIndependentDf<TT>(aDf[ii], df[ii]);
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
template <typename TT>
inline void
Problem<OWNPRBLM, CT, DIM>::setDependentVectDf(ScaFES_VectDf<TT>& aDf,
                                               ScaFES_VectDf<CT>& df)
{
    for (std::size_t ii = 0; ii < df.size(); ++ii)
    {
        this->template setDependentDf<TT>(aDf[ii], df[ii]);
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
template <typename TT>
inline void
Problem<OWNPRBLM, CT, DIM>::setIndependentVector(std::vector<TT>& aDf,
                                                 const std::vector<CT>& df)
{
    for (std::size_t ii = 0; ii < df.size(); ++ii)
    {
        this->template setIndependentScalar<TT>(aDf[ii], df[ii]);
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
template <typename TT>
inline void Problem<OWNPRBLM, CT, DIM>::setDependentVector(std::vector<TT>& aDf,
                                                           std::vector<CT>& df)
{
    for (std::size_t ii = 0; ii < df.size(); ++ii)
    {
        this->template setDependentScalar(aDf[ii], df[ii]);
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
template <typename TT>
inline void Problem<OWNPRBLM, CT, DIM>::setVectDfDepToDefValue(
    ScaFES_VectDf<TT>& dfDep, const std::vector<CT>& defVal,
    const int& /*timeIter*/
    )
{
    int idx = 0;
    for (std::size_t ii = 0; ii < this->isKnownDf().size(); ++ii)
    {
        if (!(this->isKnownDf().at(ii)))
        {
            dfDep[idx] = defVal[ii];
            ++idx;
        }
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
template <typename TT>
inline void Problem<OWNPRBLM, CT, DIM>::setVectDfDepToDfIndep(
    ScaFES_VectDf<TT>& dfDep, const ScaFES_VectDf<TT>& dfIndep,
    const int& /*timeIter*/
    )
{
    for (std::size_t jj = 0; jj < dfDep.size(); ++jj)
    {
#ifdef _OPENMP
#pragma omp parallel shared(dfDep, dfIndep)
        {
        #pragma omp for schedule(static)
#endif
        for (int ii = 0; ii < dfDep[jj].memAll().nNodesTotal(); ++ii)
        {
            int comp = 0;
            ScaFES_IntNtuple idxNode;
            dfDep[jj].memoryPos2IdxNode(idxNode, comp, ii);
            dfDep[jj](idxNode) = dfIndep[jj](idxNode);
        }
#ifdef _OPENMP
        }
#endif
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::evalZosForwardTypeOne(
    const std::string& whichPhase, ScaFES_VectDf<CT>& dfDep,
    const short int& tapeId, const std::vector<CT>& dfIndep,
    const int& /*timeIter*/
    )
{
    CT tmp;
    static_assert(std::is_same<decltype(tmp), double>::value,
                  "Only for CT=double!");

    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Evaluate phase "
                  << whichPhase << " using ADOL-C zos_forward() and tape #"
                  << tapeId << "..." << std::endl;
    }
    this->mParams.decreaseLevel();

#ifdef SCAFES_HAVE_ADOLC
    ScaFES::Timer timerPhase;
    timerPhase.restart();

#ifdef _OPENMP
    const int NTHREADS_ENV = omp_get_max_threads();
    omp_set_num_threads(1);
#endif
    // // #ifdef _OPENMP
    // // #pragma omp parallel // firstprivate(ADOLC_OpenMP_Handler)
    // //     {
    // // #endif
    const int nIndependents = dfIndep.size();
    int nDependents = 0;
    for (std::size_t ii = 0; ii < dfDep.size(); ++ii)
    {
        nDependents += dfDep[ii].nElemsAllocated();
    }

    /*------------------------------------------------------------------------*/
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * nIndependents = " << nIndependents << std::endl
                  << this->mParams.getPrefix()
                  << " * nDependents   = " << nDependents << std::endl;
    }

    /*------------------------------------------------------------------------*/
    // Points to first element of first vector of old iterate.
    //    double* vectIndep = dfIndep[0].elemData();
    const CT* vectIndep = static_cast<const CT*>(dfIndep.data());
    // Points to first element of first vector of new iterate.
    CT* vectDep = static_cast<CT*>(dfDep[0].elemData());

    /*------------------------------------------------------------------------*/
    ::zos_forward(tapeId, nDependents, nIndependents, 0, // Do not keep.
                  vectIndep, vectDep);
// // #ifdef _OPENMP
// //     }
// // #endif
#ifdef _OPENMP
    omp_set_num_threads(NTHREADS_ENV);
#endif

    this->mWallClockTimes[whichPhase] += timerPhase.elapsed();
#else
    std::ignore = dfDep;
    std::ignore = dfIndep;
#endif // HAVE_ADOLC

    this->mParams.increaseLevel();
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << "   Evaluated."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::evalZosForwardTypeTwo(
    const std::string& whichPhase, ScaFES_VectDf<CT>& dfDep,
    const short int& tapeId, const ScaFES_VectDf<CT>& dfIndep,
    const int& /*timeIter*/
    )
{
    CT tmp;
    static_assert(std::is_same<decltype(tmp), double>::value,
                  "Only for CT=double!");

    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Evaluate phase "
                  << whichPhase << " using ADOL-C zos_forward() and tape #"
                  << tapeId << "..." << std::endl;
    }
    this->mParams.decreaseLevel();

#ifdef SCAFES_HAVE_ADOLC
    ScaFES::Timer timerPhase;
    timerPhase.restart();

#ifdef _OPENMP
    const int NTHREADS_ENV = omp_get_max_threads();
    omp_set_num_threads(1);
#endif
    // // #ifdef _OPENMP
    // // #pragma omp parallel // firstprivate(ADOLC_OpenMP_Handler)
    // //     {
    // // #endif
    int nDependents = 0;
    int nIndependents = 0;
    for (std::size_t ii = 0; ii < dfIndep.size(); ++ii)
    {
        nIndependents += dfIndep[ii].nElemsAllocated();
    }
    for (std::size_t ii = 0; ii < dfDep.size(); ++ii)
    {
        nDependents += dfDep[ii].nElemsAllocated();
    }

    /*------------------------------------------------------------------------*/
    // Rearrange one-dimensional arrays to two-dimensional ones
    // such that these arrays can be passed to ADOL-C functions.
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * nIndependents = " << nIndependents << std::endl
                  << this->mParams.getPrefix()
                  << " * nDependents   = " << nDependents << std::endl;
    }

    /*------------------------------------------------------------------------*/
    // Points to first element of first vector of old iterate.
    const CT* vectIndep = static_cast<const CT*>(dfIndep[0].elemData());
    // Points to first element of first vector of new iterate.
    CT* vectDep = static_cast<CT*>(dfDep[0].elemData());

    /*------------------------------------------------------------------------*/
    ::zos_forward(tapeId, nDependents, nIndependents, 0, // Do not keep.
                  vectIndep, vectDep);

// // #ifdef _OPENMP
// //     }
// // #endif
#ifdef _OPENMP
    omp_set_num_threads(NTHREADS_ENV);
#endif

    this->mWallClockTimes[whichPhase] += timerPhase.elapsed();
#else
    std::ignore = dfDep;
    std::ignore = dfIndep;
#endif // HAVE_ADOLC

    this->mParams.increaseLevel();
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << "   Done."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::evalZosForwardTypeThree(
    const std::string& whichPhase, CT& dfDep, const short int& tapeId,
    const ScaFES_VectDf<CT>& dfIndep, const int& /*timeIter*/
    )
{
    CT tmp;
    static_assert(std::is_same<decltype(tmp), double>::value,
                  "Only for CT=double!");

    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Evaluate phase "
                  << whichPhase
                  << " using ADOL-C zos_forward() and tape #" << tapeId << "..."
                  << std::endl;
    }
    this->mParams.decreaseLevel();

#ifdef SCAFES_HAVE_ADOLC
    ScaFES::Timer timerPhase;
    timerPhase.restart();

    dfDep = static_cast<CT>(0);
    const int nDirections = this->nColumns();
    const int nDependents = 1;
    int nIndependents = 0;
    for (std::size_t ii = 0; ii < dfIndep.size(); ++ii)
    {
        nIndependents += dfIndep[ii].nElemsAllocated();
    }

    // Rearrange one-dimensional arrays to two-dimensional ones
    // such that these arrays can be passed to ADOL-C functions.
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * nIndependents = " << nIndependents << std::endl
                  << this->mParams.getPrefix()
                  << " * nDependents   = " << nDependents << std::endl
                  << this->mParams.getPrefix()
                  << " * nDirections   = " << nDirections << std::endl;
    }

    // Points to first element of first vector of old iterate.
    CT* vectIndep = static_cast<CT*>(dfIndep[0].elemData());
    // Points to first element of a new arrary.
    std::vector<CT> tmpVectDep(1);
    CT* vectDep = tmpVectDep.data();

    /*------------------------------------------------------------------------*/
    ::zos_forward(tapeId, nDependents, nIndependents, 0, // Do not keep.
                  vectIndep, vectDep);
    // Assign value of local array.
    dfDep = vectDep[0];

    this->mWallClockTimes[whichPhase] += timerPhase.elapsed();
#else
    std::ignore = dfDep;
    std::ignore = dfIndep;
#endif

    this->mParams.increaseLevel();
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << "   Evaluated."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::evalFovForwardTypeOne(
    const std::string& whichPhase, ScaFES_VectDf<CT>& dfDep,
    ScaFES_VectDf<CT>& dfGradDep, const short int& tapeId,
    const std::vector<CT>& dfIndep, const std::vector<CT>& dfGradIndep,
    const int& /*timeIter*/
    )
{
    CT tmp;
    static_assert(std::is_same<decltype(tmp), double>::value,
                  "Only for CT=double!");

    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Evaluate phase "
                  << whichPhase << " using ADOL-C fov_forward() and tape #"
                  << tapeId << "..." << std::endl;
    }
    this->mParams.decreaseLevel();

#ifdef SCAFES_HAVE_ADOLC
    ScaFES::Timer timerPhase;
    timerPhase.restart();

#ifdef _OPENMP
    const int NTHREADS_ENV = omp_get_max_threads();
    omp_set_num_threads(1);
#endif
    const int nDirections = this->nColumns();
    const int nIndependents = dfIndep.size();
    int nDependents = 0;
    for (std::size_t ii = 0; ii < dfDep.size(); ++ii)
    {
        nDependents += dfDep[ii].nElemsAllocated();
    }

    /*------------------------------------------------------------------------*/
    // Rearrange one-dimensional arrays to two-dimensional ones
    // such that these arrays can be passed to ADOL-C functions.
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout
            << this->mParams.getPrefix()
            << " * nIndependents = " << nIndependents << std::endl
            << this->mParams.getPrefix()
            << " * nDependents   = " << nDependents << std::endl
            << this->mParams.getPrefix()
            << " * nDirections   = " << nDirections << std::endl;
    }

    /*------------------------------------------------------------------------*/
    // Points to first element of gradient of parameter vector.
    CT* elemCurrIndep = const_cast<CT*>(dfGradIndep.data());
    for (int ii = 0; ii < nIndependents; ++ii)
    {
        this->mMemoryVectGradParamsAsMatrix[ii] = elemCurrIndep;
        elemCurrIndep += nDirections;
    }
    CT* elemCurrDep = dfGradDep[0].elemData();
    // Just take memory lump #1 OR memory lump #2.
    for (int ii = 0; ii < nDependents; ++ii)
    {
        this->mMemoryVectGradUnknownDfsRowsTwo[ii] = elemCurrDep;
        elemCurrDep += nDirections;
    }

    /*------------------------------------------------------------------------*/
    // Set gradient of independent variables to I.
    for (int ii = 0; ii < nIndependents; ++ii)
    {
        for (int jj = 0; jj < nDirections; ++jj)
        {
            this->mMemoryVectGradParamsAsMatrix[ii][jj] = static_cast<CT>(0);
            if (jj == ii)
            {
                this->mMemoryVectGradParamsAsMatrix[ii][jj] =
                    static_cast<CT>(1);
            }
        }
    }

    /*------------------------------------------------------------------------*/
    // Points to first element of first vector.
    const CT* vectIndep = dfIndep.data();
    CT* vectDep = dfDep[0].elemData();
    CT** vectIndepGrad = this->mMemoryVectGradParamsAsMatrix.data();
    CT** vectDepGrad = this->mMemoryVectGradUnknownDfsRowsTwo.data();

    /*------------------------------------------------------------------------*/
    ::fov_forward(tapeId, nDependents, nIndependents, nDirections, vectIndep,
                  vectIndepGrad, vectDep, vectDepGrad);

// // #ifdef _OPENMP
// //     }
// // #endif
#ifdef _OPENMP
    omp_set_num_threads(NTHREADS_ENV);
#endif

    this->mWallClockTimes[whichPhase] += timerPhase.elapsed();
#else
    std::ignore = dfDep;
    std::ignore = dfGradDep;
    std::ignore = dfIndep;
    std::ignore = dfGradIndep;
#endif // HAVE_ADOLC

    /*------------------------------------------------------------------------*/
    // Rearrange back two-dimensional arrays (gradients)
    // to one-dimensional arrays.
    // NOT necessary: Values are already at correct position in memory.
    // Local variables vectDep, vectGradDep etc. can be thrown away.
    // These variables are just pointers to the memory of the data fields.

    this->mParams.increaseLevel();
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << "   Evaluated."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::evalFovForwardTypeTwo(
    const std::string& whichPhase, ScaFES_VectDf<CT>& dfDep,
    ScaFES_VectDf<CT>& dfGradDep, const short int& tapeId,
    const ScaFES_VectDf<CT>& dfIndep, const ScaFES_VectDf<CT>& dfGradIndep,
    const int& /*timeIter*/
    )
{
    CT tmp;
    static_assert(std::is_same<decltype(tmp), double>::value,
                  "Only for CT=double!");

    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Evaluate phase "
                  << whichPhase << " using ADOL-C fov_forward() and tape #"
                  << tapeId << "..." << std::endl;
    }
    this->mParams.decreaseLevel();

#ifdef SCAFES_HAVE_ADOLC
    ScaFES::Timer timerPhase;
    timerPhase.restart();

#ifdef _OPENMP
    const int NTHREADS_ENV = omp_get_max_threads();
    omp_set_num_threads(1);
#endif
    // // #ifdef _OPENMP
    // // #pragma omp parallel // firstprivate(ADOLC_OpenMP_Handler)
    // //     {
    // // #endif
    //    int kk = 0;
    const int nDirections = this->nColumns();
    int nDependents = 0;
    int nIndependents = 0;
    for (std::size_t ii = 0; ii < dfIndep.size(); ++ii)
    {
        nIndependents += dfIndep[ii].nElemsAllocated();
    }
    for (std::size_t ii = 0; ii < dfDep.size(); ++ii)
    {
        nDependents += dfDep[ii].nElemsAllocated();
    }

    // Rearrange one-dimensional arrays to two-dimensional ones
    // such that these arrays can be passed to ADOL-C functions.
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout
            << this->mParams.getPrefix()
            << " * nIndependents = " << nIndependents << std::endl
            << this->mParams.getPrefix()
            << " * nDependents   = " << nDependents << std::endl
            << this->mParams.getPrefix()
            << " * nDirections   = " << nDirections << std::endl;
    }

    CT* elemCurr = const_cast<CT*>(dfGradIndep[0].elemData());
    for (int ii = 0; ii < nIndependents; ++ii)
    {
        this->mMemoryVectGradUnknownDfsRowsOne[ii] = elemCurr;
        elemCurr += nDirections;
    }
    elemCurr = dfGradDep[0].elemData();
    for (int ii = 0; ii < nDependents; ++ii)
    {
        this->mMemoryVectGradUnknownDfsRowsTwo[ii] = elemCurr;
        elemCurr += nDirections;
    }

    // Points to first element of first vector of old iterate.
    CT* vectIndep = dfIndep[0].elemData();
    // Points to first element of first vector of new iterate.
    CT* vectDep = dfDep[0].elemData();
    CT** vectIndepGrad = this->mMemoryVectGradUnknownDfsRowsOne.data();
    CT** vectDepGrad = this->mMemoryVectGradUnknownDfsRowsTwo.data();

    /*------------------------------------------------------------------------*/
    ::fov_forward(tapeId, nDependents, nIndependents, nDirections, vectIndep,
                  vectIndepGrad, vectDep, vectDepGrad);

/*------------------------------------------------------------------------*/
// Rearrange back two-dimensional arrays (gradients)
// to one-dimensional arrays.
// NOT necessary: Values are already at correct position in memory.
// Local variables vectDep, vectGradDep etc. can be thrown away.
// These variables are just pointers to the memory of the data fields.
// // #ifdef _OPENMP
// //     }
// // #endif
#ifdef _OPENMP
    omp_set_num_threads(NTHREADS_ENV);
#endif

    this->mWallClockTimes[whichPhase] += timerPhase.elapsed();
#else
    std::ignore = dfDep;
    std::ignore = dfGradDep;
    std::ignore = dfIndep;
    std::ignore = dfGradIndep;
#endif // HAVE_ADOLC

    this->mParams.increaseLevel();
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << "   Done."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <class OWNPRBLM, typename CT, std::size_t DIM>
inline void Problem<OWNPRBLM, CT, DIM>::evalFovForwardTypeThree(
    const std::string& whichPhase, CT& dfDep, std::vector<CT>& dfGradDep,
    const short int& tapeId, const ScaFES_VectDf<CT>& dfIndep,
    const ScaFES_VectDf<CT>& dfGradIndep, const int& /*timeIter*/
    )
{
    CT tmp;
    static_assert(std::is_same<decltype(tmp), double>::value,
                  "Only for CT=double!");

    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * Evaluate phase "
                  << whichPhase
                  << " using ADOL-C fov_forward() and tape #" << tapeId << "..."
                  << std::endl;
    }
    this->mParams.decreaseLevel();

#ifdef SCAFES_HAVE_ADOLC
    ScaFES::Timer timerPhase;
    timerPhase.restart();

    dfDep = static_cast<CT>(0);
    const int nDirections = this->nColumns();
    const int nDependents = 1;
    int nIndependents = 0;
    for (std::size_t ii = 0; ii < dfIndep.size(); ++ii)
    {
        nIndependents += dfIndep[ii].nElemsAllocated();
    }

    // Rearrange one-dimensional arrays to two-dimensional ones
    // such that these arrays can be passed to ADOL-C functions.
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << " * nIndependents = " << nIndependents << std::endl
                  << this->mParams.getPrefix()
                  << " * nDependents   = " << nDependents << std::endl
                  << this->mParams.getPrefix()
                  << " * nDirections   = " << nDirections << std::endl;
    }

    // Use "this->mMemoryVectGradUnknownDfsRowsTwo"
    // because data fields will not be swapped between updateUnknownDfs()
    // and evalTfPhase().
    CT* elemCurr = const_cast<CT*>(dfGradIndep[0].elemData());
    for (int ii = 0; ii < nIndependents; ++ii)
    {
        this->mMemoryVectGradUnknownDfsRowsTwo[ii] = elemCurr;
        elemCurr += nDirections;
    }

    // Points to first element of first vector of old iterate.
    CT* vectIndep = dfIndep[0].elemData();
    CT** vectIndepGrad = this->mMemoryVectGradUnknownDfsRowsTwo.data();
    // Points to first element of first vector of new iterate.
    std::vector<CT> tmpVectDep(1);
    CT* vectDep = tmpVectDep.data();
    CT** vectDepGrad = new CT* [1];
    vectDepGrad[0] = dfGradDep.data();

    ::fov_forward(tapeId, nDependents, nIndependents, nDirections, vectIndep,
                  vectIndepGrad, vectDep, vectDepGrad);

    dfDep = vectDep[0];
    delete[] vectDepGrad;
    /*------------------------------------------------------------------------*/
    // Rearrange back two-dimensional arrays (gradients)
    // to one-dimensional arrays.
    // NOT necessary: Values are already at correct position in memory.
    // Local variables vectDep, vectGradDep etc. can be thrown away.
    // These variables are just pointers to the memory of the data fields.
    // Not necessary here:
    // vectDepGrad[0] points already to this->mVectGradTf.

    this->mWallClockTimes[whichPhase] += timerPhase.elapsed();
#else
    std::ignore = dfDep;
    std::ignore = dfGradDep;
    std::ignore = dfIndep;
    std::ignore = dfGradIndep;
#endif

    this->mParams.increaseLevel();
    if ((this->params().rankOutput() == this->myRank()) &&
        (0 < this->params().indentDepth()))
    {
        std::cout << this->mParams.getPrefix()
                  << "   Evaluated."
                  << std::endl;
    }
}
} // End of namespace.
#endif
