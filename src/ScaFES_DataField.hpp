/* ScaFES
 * Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file ScaFES_DataField.hpp
 *  @brief Contains the class template DataField.
 */

#ifndef SCAFES_DATAFIELD_HPP_
#define SCAFES_DATAFIELD_HPP_

#include "ScaFES_Config.hpp"

#include <iostream>
#include <algorithm>
#include <utility>
#include <vector>
#include <cstring>
#include <string>
#include <type_traits>
#include <stdexcept>

#ifdef VTRACE
#include "vt_user.h"
#endif

#ifdef _OPENMP
extern "C" {
#include <omp.h>
}
#endif

#include "ScaFES_Ntuple.hpp"
#include "ScaFES_Parameters.hpp"
#include "ScaFES_Grid.hpp"
#include "ScaFES_GridSub.hpp"
#include "ScaFES_GridGlobal.hpp"
#include "ScaFES_Buffer.hpp"
#include "ScaFES_DataFile.hpp"

namespace ScaFES
{

/*******************************************************************************
 ******************************************************************************/
/**
 * \class DataField
 * @brief The class template \c DataField represents a physical field on a
 * grid partition.
*
* The global grid is subdivided into so-called grid partitions.
* Physical fields which are defined on the global grid are subdivided
* into subfields and mapped to the corresponding grid partitions.
* An object of this class stores the function values at all grid partition
* nodes. The function values at an n-dimensional grid are stored in
* an one-dimensional array. This ensures that one memory lump is allocated
* for the data field, only.
*
* An object of type \c DataField can have an extended grid (so called ghost
* cells with ghost nodes) where neighbouring processes store their
* function values. This means, all function values at these ghost nodes
* are stored at multiple locations to enable and speed up the computation.
*
* The node numbering starts with zero at the first ghost node of the first
* ghost cell. If there are no ghost cells, then the first local node number
* will be zero.
*
* In order to access the elements of a three-dimensional data field
* named \c E, the following possibilities exist:
* <ul>
* <li> \code int i, j, k;                 E(i,j,k) = ...\endcode </li>
* <li> \code unsigned long int memoryPos; E(memoryPos) = ...\endcode </li>
* <li> \code Ntuple<int,DIM> idxNode;     E(idxNode) = ...\endcode </li>
* </ul>
*
* \remarks \c memoryPos is the MPI process local memory position.
* This memory position is only the same for those data fields which have
* the same ghost cell width.
* \c idxNode returns the MPI global node number.
*  ==> Using \c memoryPos is better for performance reasons:
*  It avoids unnecessary recalculations (using the method
*  \c idxNode2MemoryPos()) for accessing the correct memory locations,
*  but: if one is unsure about indexing, one should use \c idxNode:
*  It always gives the global node number and thus, the correct memory location
*  can be accessed.
*
* \remarks There are two variants for each access method:
* According to the class \c std::vector<CT>, the \c operator(node) is overloaded
* and returns the addressed component. One the other hand, the method
* \code at(node) \endcode performs a range check and returns the addressed
* component.
*
* \remarks DataField itself has no communicator attribute. For the actual communication,
* the communicator of the \c Parameter attribute is used.
*/
template <typename CT, std::size_t DIM> class DataField
{
public:
    /*--------------------------------------------------------------------------
    | TYPE DEFINITONS.
    --------------------------------------------------------------------------*/
    /** Type definition for integer d-dimensional tuple. */
    typedef ScaFES::Ntuple<int, DIM> ScaFES_IntNtuple;

    /** Re-export typename CT STL-like. */
    typedef CT value_type;

    /*--------------------------------------------------------------------------
    | LIFE CYCLE METHODS.
    ----------------------------------------------------------------------*/
    /** Creates default constructor. */
    DataField(ScaFES::Parameters* const params);

    /** Creates another own constructor. */
    DataField(const std::string& nameDataField,
              ScaFES::Parameters* const params,
              const ScaFES::GridGlobal<DIM>& gg,
              const ScaFES::Grid<DIM>& memAll, const int& stencilWidth,
              CT* const elemData, const int& nColumns = 1,
              const ScaFES::GridSub<DIM>& memHole = ScaFES::GridSub<DIM>(),
              const int& nLayers = 0,
              const ScaFES::WriteHowOften& writeToFile =
                  ScaFES::WriteHowOften::LIKE_GIVEN_AT_CL);

    /** Creates compiler provided copy constructor. */
    DataField(const ScaFES::DataField<CT, DIM>& /*from*/) = default;

    /** Creates own destructor. */
    ~DataField() = default;

    /** Creates compiler provided assignment operator. */
    ScaFES::DataField<CT, DIM>&
    operator=(const ScaFES::DataField<CT, DIM>& /*rhs*/) = default;

    /*--------------------------------------------------------------------------
    | GETTER METHODS.
    --------------------------------------------------------------------------*/
    /** Returns the parameter set. */
    const std::string& name() const;

    /** Returns pointer to parameters. */
    ScaFES::Parameters* params() const;

    /** Returns the global grid which includes all partitions. */
    const ScaFES::GridGlobal<DIM>& gridGlobal() const;

    /** Returns the number of ghost cells. */
    const int& stencilWidth() const;

    /** Returns the grid incl. ghost nodes of the data field. */
    const ScaFES::Grid<DIM>& memAll() const;

    /** Returns the normal grid of the data field. */
    const ScaFES::GridSub<DIM>& memNormal() const;

    /** Returns the number of layers at the boundary. */
    const int& nLayers() const;

    /** Returns the grid where the data field is NOT defined. */
    const ScaFES::GridSub<DIM>& memHole() const;

    /** Returns the interior grid of the data field. */
    const ScaFES::GridSub<DIM>& memInner() const;

    /** Returns the element size of the underlying type. */
    std::size_t elementSize() const;

    /** Returns the memory. */
    CT* elemData() const;

    /** Returns the element \c idx of the memory. */
    const CT& elemData(const int& ii) const;

    /** Returns the border grid of the data field. */
    const std::vector<ScaFES::GridSub<DIM>>& memBorder() const;

    /** Returns the communication grid of the data field. */
    const std::vector<ScaFES::GridSub<DIM>>& memComm() const;

    /** Returns the ghost grid of the data field. */
    const std::vector<ScaFES::GridSub<DIM>>& memGhost() const;

    /** Returns the number of columns (in case of a 2D array). */
    const int& nColumns() const;

    /** Returns the number of rows. */
    const int& nRows() const;

    /** Returns the number of allocated elements. */
    int nElemsAllocated() const;

    /*--------------------------------------------------------------------------
    | SETTER METHODS.
    --------------------------------------------------------------------------*/
    /** Set the current time. */
    CT& time() const;

    /** Address elements via local node numbers.
     * Returns a reference, what makes an assignment like \e A(0) = 2
     * possible.
     * \param memoryPos Local node number.
     */
    CT& operator()(const unsigned long int& memoryPos) const;

    /** Address elements via global node numbers.
     * \param idxNode Global node number.
     * \param comp Given component at node number
     */
    CT& operator()(const ScaFES_IntNtuple& idxNode, const int& comp = 0) const;

    /** Address elements via global node number given explicitly
     *  as (ii,jj,kj).
     * \param ii First entry of global node number.
     * \param jj Second entry of global node number.
     * \param kk Third entry of global node number.
     * \param comp Given component at node number
     * \remark Only valid in 3D.
     */
    CT& operator()(const int& ii, const int& jj, const int& kk,
                   const int& comp = 0) const;

    /** Addresses elements via global node number given explicitly
     *  as (ii,jj).
     * \param ii First entry of global node number.
     * \param jj Second entry of global node number.
     * \param comp Given component at node number
     * \remark Only valid in 2D.
     */
    CT& operator()(const int& ii, const int& jj, const int& comp = 0) const;

    /** Addresses elements via local (=grid partition) scalar node numbers.
     * This method includes a range check.
     * \param memoryPos Local node number.
     */
    CT& at(const unsigned long int& memoryPos) const;

    /** Addresses elements via global node numbers.
     * This method includes a range check.
     * \param idxNode Given node number.
     * \param comp Given component at node number
     */
    CT& at(const ScaFES_IntNtuple& idxNode, const int& comp = 0) const;

    /** Overloads ()-operator to address elements via global node number
     * (ii,jj,kk). This method includes a range check.
     * \param ii First entry of global node number.
     * \param jj Second entry of global node number.
     * \param kk Third entry of global node number.
     * \param comp Given component at node number
     * \remark Only valid in 3D.
     */
    CT& at(const int& ii, const int& jj, const int& kk,
           const int& comp = 0) const;

    /*--------------------------------------------------------------------------
    | COMPARISON METHODS.
     -------------------------------------------------------------------------*/
    /** Compares two data fields for equality.
     * First, geometrical sizes are compared.
     * In a second step, the method haveSameValues(DF1, DF2) compares
     * node values.
     */
    bool operator==(const ScaFES::DataField<CT, DIM>& second);

    /** Checks if this data field and the rhs data field have the same
     * geometrical size.
     */
    bool hasSameSize(const ScaFES::DataField<CT, DIM>& rhs) const;

    /** Checks if this data field and the rhs data field have the same
     *  values.
     */
    bool hasSameValues(const ScaFES::DataField<CT, DIM>& rhs,
                       const double& eps) const;

    /*--------------------------------------------------------------------------
    | WORK METHODS.
    --------------------------------------------------------------------------*/
    /** Synchronizes the data exchange of the ghost cell values with all
     *  neighbours by executing the following steps:
     *  1. Copy values to send buffers.
     *  2. Communicate values.
     *  3. Wait to finish communication.
     *  4. Copy values from receive buffers to ghost cells.
     */
    void sync(const int& timeIter);

    /** Sums up and synchronizes the data exchange of the ghost cell values
     *  with all neighbours.
     */
    void collectValuesAtMemComm(const int& timeIter);

    /** Copies the values of the memory of the communication nodes
     *  of this data field to the send buffer.
     */
    void copyValuesFromMemCommToSendBuffer(const int& timeIter);

    /** Copies the values of the memory of the ghost nodes
     *  of this data field to the send buffer.
     */
    void copyValuesFromMemGhostToSendBuffer(const int& timeIter);

    /** Communicates asynchronously with direct neighbours.
     * Manages the sending and receiving values of the ghost
     * nodes of a data field to the corresponding neighboured data field.
     */
    void exchangeValuesInBuffers(const int& timeIter);

    /** Waits until all communication calls have been finished. */
    void waitAll();

    /** Copies the values of the receive buffer to the memory
     *  of the ghost nodes of this data field.
     */
    void copyValuesFromReceiveBufferToMemGhost(const int& timeIter);

    /** Copies the values of the receive buffer to the memory
     *  of the communication nodes of this data field.
     */
    void copyValuesFromReceiveBufferToMemComm(const int& timeIter);

    /** Copies and sums the values of the receive buffer to the memory
     *  of the communication nodes of this data field.
     */
    void copyAndSumValuesFromReceiveBufferToMemComm(const int& timeIter);

    /** Writes the data field to a NetCDF file. */
    void write(const int& timeIter);

    /** Reads the data field from a NetCDF file. */
    void read(const int& timeIter);

    /** Adds the values at all grid nodes of a given data field
     *  to the existing one. */
    ScaFES::DataField<CT, DIM>&
    operator+=(const ScaFES::DataField<CT, DIM>& rhs);

    /** Subtracts the values at all grid nodes of a given data field
     *  to the existing one. */
    ScaFES::DataField<CT, DIM>&
    operator-=(const ScaFES::DataField<CT, DIM>& rhs);

    /** Adds a given value at all grid nodes from the underlying
     * one. */
    ScaFES::DataField<CT, DIM>& operator+=(const CT& rhs);

    /** Subtracts a given value at all grid nodes from the underlying
     * one. */
    ScaFES::DataField<CT, DIM>& operator-=(const CT& rhs);

    /** Sets the elements of a data field at all nodes to a given
     * value \c val.
     * \param val Given value.
     */
    ScaFES::DataField<CT, DIM>& operator=(const CT& val);

    /** Computes the corresponding memory position from a given
     * node number.
     * \param memoryPos Memory position (return value).
     * \param idxNode Given node number.
     * \param comp Given component at node number.
     */
    void idxNode2MemoryPos(unsigned long int& memoryPos,
                           const ScaFES_IntNtuple& idxNode,
                           const int& comp = 0) const;

    /** Computes the corresponding memory position from a given
     *  scalar node number.
     * \param memoryPos Memory position (return value).
     * \param idxNodeScalar Given scalar node number.
     * \param comp Given component at node number.
     */
    void idxNodeScalar2MemoryPos(unsigned long int& memoryPos,
                                 const unsigned long int& idxNodeScalar,
                                 const int& comp = 0) const;

    /** Gets the memory position with considered hole.
     * \param memoryPosWithHole Memory position (return value).
     * \param memoryPos Given memory position (hole not considered).
     */
    void getMemoryPosWithHole(unsigned long int& memoryPosWithHole,
                              const unsigned long int& memoryPos) const;

    /** Computes the corresponding global node number from a given
     * local memory position.
     * \param idxNode Node number (return value).
     * \param comp Component at node number (return value).
     * \param memoryPos Given memory position.
     */
    void memoryPos2IdxNode(ScaFES_IntNtuple& idxNode, int& comp,
                           const unsigned long int& memoryPos) const;

    /** Computes the corresponding global node number from a given
     * local memory position.
     * \param idxNodeScalar Scalar node number (return value).
     * \param comp Component at node number (return value).
     * \param memoryPos Given memory position.
     */
    void memoryPos2IdxNodeScalar(unsigned long int& idxNodeScalar, int& comp,
                                 const unsigned long int& memoryPos) const;

    /** Gets the memory position without considered hole.
     * \param memoryPosWithoutHole Memory position (return value).
     * \param memoryPos Given memory position (hole considered).
     */
    void getMemoryPosWithoutHole(unsigned long int& memoryPosWithoutHole,
                                 const unsigned long int& memoryPos) const;

    /** Swaps the values of two data fields.
     * \remarks For time dependent equations, it is necessary to change the
     * meaning of the new and the old data field.
     * The old data field will be reused as new data field in the next
     * time step.
     */
    void swapPointerToMemoryWith(ScaFES::DataField<CT, DIM>& df);

    /** Assigns node values from one data field to another data field.
     * \remarks Copies values of all (normal + ghost) nodes!
     */
    void assignValues(const ScaFES::DataField<CT, DIM>& from);

    /** Checks if this data field and another given data field have the
     * same values. */
    void check(const ScaFES::DataField<CT, DIM>& df2);

    /** Set values at all normal nodes from a given data field. */
    void setValuesAtMemNormal(const ScaFES::DataField<CT, DIM>& df);

    /** Set values at all nodes from a given data field. */
    void setValuesAtMemAll(const ScaFES::DataField<CT, DIM>& df);

    /** Computes the L^{inf} norm (maximum norm) of this data field. */
    double normLinf() const;

    /** Computes the error of this and another given data field
     * using the maximum norm. */
    double compErrLinf(const ScaFES::DataField<CT, DIM>& dfAna);

    /** Overloads the output operator \c operator<<. */
    template <typename ST, std::size_t DD>
    friend std::ostream& operator<<(std::ostream& output,
                                    const ScaFES::DataField<ST, DD>& df);

protected:
    /*--------------------------------------------------------------------------
    | MEMBER VARIABLES.
    --------------------------------------------------------------------------*/
    /** Parameter set. */
    ScaFES::Parameters* mParams;

    /** Global grid including all grid partitions. */
    ScaFES::GridGlobal<DIM> mGG;

    /** Underlying grid where the data field lives:
     * Roughly spoken: grid partition + grid(stencilWidth).
     * A data field can be restricted to a partition of the global grid.
     * Thus, \c mMemAll.idxNodeFirst() indicates the first global node
     * number of the data field (numbering incl. nodes of ghost cells).
     */
    ScaFES::Grid<DIM> mMemAll;

    /** The normal grid: All grid nodes without the nodes at the ghost
     * cells.
     */
    ScaFES::GridSub<DIM> mMemNormal;

    /** (Global) grid where the data field is NOT defined
     * (derived from the global grid).
     */
    ScaFES::GridSub<DIM> mMemHole;

    /** Grid of all border nodes. */
    std::vector<ScaFES::GridSub<DIM>> mMemBorder;

    /** Memory of all communication nodes. */
    std::vector<ScaFES::GridSub<DIM>> mMemComm;

    /** Memory of all ghost nodes. */
    std::vector<ScaFES::GridSub<DIM>> mMemGhost;

    /** Number of ghost cells (uniform in all directions). */
    int mStencilWidth;

    /** Name of data field. */
    std::string mName;

    /** Used in assertions to make sure that only elements in grid are
     * addressed. */
    std::size_t mMemoryPosGridAllFirst;

    /** Used in assertions to make sure that only elements in grid are
     * addressed. */
    std::size_t mMemoryPosGridAllLast;

    /** Memory position of first node of hole. */
    std::size_t mMemoryPosGridHoleFirst;

    /** Memory position of last node of hole. */
    std::size_t mMemoryPosGridHoleLast;

    /** Pointer to linear continuous memory lump
     * for the elements of a data field.
     * A data field can only store values of one type.
     * Moreover, the type remains constant during the life time
     * of the data field.
     */
    CT* mElemData;

    /** Number of allocated elements. */
    int mNrows;

    /** Number of columns (important for matrices). */
    int mNcolumns;

    /** Buffer for all elements of the data field which should be sent
     * to and received from the neighbour.
     */
    std::vector<ScaFES::Buffer<CT>> mValuesToExchange;

    /** Data file. */
    ScaFES::DataFile<CT, DIM> mOutput;

    /** Number layers at the boundary. */
    int mNlayers;

    /** Write status. */
    ScaFES::WriteHowOften mWriteToFile;

    /** List of interfaces to the neighbours. */
    // std::vector<Communicator> vCommunicator;

private:
    /*--------------------------------------------------------------------------
    | WORK METHODS.
    --------------------------------------------------------------------------*/
    /** Initializes dependent member variables. */
    void initDependentMembers();

    /** Performs range check. */
    void checkRange(const unsigned long int& memoryPos) const;

    /** Computes the corresponding memory position from a given
     * global node number using the C style numbering.
     * \param idxNodeScalar Scalar node number (return value).
     * \param idxNode Given global node number.
     * \param comp Given comp of at node number.
     */
    void idxNode2IdxNodeScalarStyleC(unsigned long int& idxNodeScalar,
                                     const ScaFES_IntNtuple& idxNode,
                                     const int& comp = 0) const;

    /** Computes the corresponding memory position from a given
     * global node number using the FORTRAN style numbering.
     * \param idxNodeScalar Scalar node number (return value).
     * \param idxNode Given node number.
     * \param comp Given comp of at node number.
     */
    void idxNode2IdxNodeScalarStyleFortran(unsigned long int& idxNodeScalar,
                                           const ScaFES_IntNtuple& idxNode,
                                           const int& comp = 0) const;

    /** Computes the corresponding global node number from a given
     * local memory position using the C style numbering.
     * \param idxNode Global node number (return value).
     * \param comp Comp of at node number (return value).
     * \param idxNodeScalar Given scalar node number.
     */
    void idxNodeScalar2IdxNodeStyleC(ScaFES_IntNtuple& idxNode, int& comp,
                                     const unsigned long int& idxNodeScalar) const;

    /** Computes the corresponding global node number from a given
     * local memory position using the FORTRAN style numbering.
     * \param idxNode Global node number (return value).
     * \param comp Comp of at node number (return value).
     * \param idxNodeScalar Given scalar node number.
     */
    void idxNodeScalar2IdxNodeStyleFortran(ScaFES_IntNtuple& idxNode, int& comp,
                                           const unsigned long int& idxNodeScalar) const;

}; // End of class. //

/*******************************************************************************
 * FREE METHODS.
 ******************************************************************************/
/** Preprares writing a data field to output. */
template <typename CT, std::size_t DIM>
std::ostream& operator<<(std::ostream& output,
                         const ScaFES::DataField<CT, DIM>& df);
/*----------------------------------------------------------------------------*/
/** Adds the elements of two data fields. */
template <typename CT, std::size_t DIM>
const ScaFES::DataField<CT, DIM>
operator+(const ScaFES::DataField<CT, DIM>& lhs,
          const ScaFES::DataField<CT, DIM>& rhs);
/*----------------------------------------------------------------------------*/
/** Subtracts the elements of two data fields. */
template <typename CT, std::size_t DIM>
const ScaFES::DataField<CT, DIM>
operator-(const ScaFES::DataField<CT, DIM>& lhs,
          const ScaFES::DataField<CT, DIM>& rhs);

/*******************************************************************************
 * LIFE CYCLE METHODS.
 ******************************************************************************/
template <typename CT, std::size_t DIM>
inline DataField<CT, DIM>::DataField(ScaFES::Parameters* const params)
: mParams(params), mGG(), mMemAll(), mMemNormal(), mMemHole(), mMemBorder(),
  mMemComm(), mMemGhost(), mStencilWidth(0), mName("EMPTYNAME"),
  mMemoryPosGridAllFirst(0L), mMemoryPosGridAllLast(0L),
  mMemoryPosGridHoleFirst(0L), mMemoryPosGridHoleLast(0L), mElemData(0),
  mNrows(1), mNcolumns(1), mValuesToExchange(), mOutput(), mNlayers(0),
  mWriteToFile(ScaFES::WriteHowOften::LIKE_GIVEN_AT_CL)
{
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline DataField<CT, DIM>::DataField(
    const std::string& nameDataField, ScaFES::Parameters* const params,
    const ScaFES::GridGlobal<DIM>& gg, const ScaFES::Grid<DIM>& memAll,
    const int& stencilWidth, CT* const elemData, const int& nColumns,
    const ScaFES::GridSub<DIM>& memHole, const int& nLayers,
    const ScaFES::WriteHowOften& writeToFile)
: mParams(params), mGG(gg), mMemAll(memAll),
  mMemNormal(memAll.idxNodeFirst() + ScaFES_IntNtuple(stencilWidth),
             memAll.idxNodeLast() - ScaFES_IntNtuple(stencilWidth), memAll),
  mMemHole(memHole), mMemBorder(), mMemComm(), mMemGhost(),
  mStencilWidth(stencilWidth), mName(nameDataField), mMemoryPosGridAllFirst(0L),
  mMemoryPosGridAllLast(0L), mMemoryPosGridHoleFirst(0L),
  mMemoryPosGridHoleLast(0L), mElemData(elemData),
  mNrows(memAll.nNodes().size() + 1) // Default: Adapted below! (incl. t)
  ,
  mNcolumns(nColumns), mValuesToExchange(), mOutput(), mNlayers(nLayers),
  mWriteToFile(writeToFile)
{
    initDependentMembers();
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline void DataField<CT, DIM>::initDependentMembers()
{
    assert(0 <= this->stencilWidth());
    assert(
        this->gridGlobal().discreteDomain().idxNodeFirst() <=
        this->gridGlobal().partition(this->params()->rank()).idxNodeFirstSub());
    assert(
        this->gridGlobal().partition(this->params()->rank()).idxNodeLastSub() <=
        this->gridGlobal().discreteDomain().idxNodeLast());

    assert(0 == this->mMemoryPosGridAllFirst);

    if (0 < nLayers())
    {
        this->mNrows = this->memAll().nNodes().size() -
                       this->memHole().nNodesSub().size() + 1; // incl. time
    }

    std::ostringstream tmpStringstream;
    tmpStringstream << this->params()->nameDataFile();
    tmpStringstream << "_" << this->name();
    std::vector<std::string> tmpName;
    tmpName.push_back(this->name());
    std::vector<GridSub<DIM>> tmpMemNormal;
    tmpMemNormal.push_back(this->memNormal());
    this->mOutput = ScaFES::DataFile<CT, DIM>(
        tmpStringstream.str(), tmpName, this->params()->myWorld(),
        this->gridGlobal().discreteDomain().nNodes(), tmpMemNormal);

    // Check if the value of border width is too large.
    // Is the grid too small to set all faces?
    // There should be at least 3*stencilWidth()() nodes in every direction.
    ScaFES_IntNtuple itStart = this->memNormal().idxNodeFirstSub();
    ScaFES_IntNtuple itEnd = this->memNormal().idxNodeLastSub();
    ScaFES_IntNtuple diff = this->memNormal().nNodesSub();
    int bw = std::max(this->stencilWidth(), 1);

    // Create iterators for all boundary faces of the current underlying
    // grid partition.
    for (std::size_t ii = 0; ii < DIM; ++ii)
    {
        // Does a face fit?
        if (diff.elem(ii) >= bw)
        {
            // Left, front, bottom, ...
            ScaFES_IntNtuple tmpIdxNodeFirst(itStart);
            ScaFES_IntNtuple tmpIdxNodeLast(itEnd);
            tmpIdxNodeLast[ii] = itStart.elem(ii) + bw;

            for (std::size_t jj = 0; jj < ii; ++jj)
            {
                tmpIdxNodeFirst[jj] += bw;
                tmpIdxNodeLast[jj] -= bw;
            }

            ScaFES::GridSub<DIM> gridToAdd(tmpIdxNodeFirst, tmpIdxNodeLast,
                                           this->memAll());
            this->mMemBorder.push_back(gridToAdd);

            // Does an additional face fit?
            if (diff.elem(ii) > 2 * bw)
            {
                // Right, back, top,...
                ScaFES_IntNtuple tmpIdxNodeFirst2(itStart);
                ScaFES_IntNtuple tmpIdxNodeLast2(itEnd);
                tmpIdxNodeFirst2[ii] = itEnd.elem(ii) - bw;

                for (std::size_t jj = 0; jj < ii; ++jj)
                {
                    tmpIdxNodeFirst2[jj] += bw;
                    tmpIdxNodeLast2[jj] -= bw;
                }

                ScaFES::GridSub<DIM> gridToAddFit(
                    tmpIdxNodeFirst2, tmpIdxNodeLast2, this->memAll());
                this->mMemBorder.push_back(gridToAddFit);
            }
        }
    }

    // Add buffers for the communication.
    if (0 < this->stencilWidth())
    {
        int dirNeighbour;
        int idNeighbour;

        ScaFES::GridSub<DIM> memFaceCurr;
        ScaFES::GridSub<DIM> memCommCurr;
        ScaFES::GridSub<DIM> memExtendedNeighbour;

        // Loop over all neighbours of the current grid partition.
        int nNeighbours =
            this->gridGlobal().neighbourId(this->params()->rank()).size();

        for (int ii = 0; ii < nNeighbours; ++ii)
        {
            // Create grid of neighbour.
            idNeighbour =
                this->gridGlobal().neighbourId(this->params()->rank(), ii);
            dirNeighbour =
                this->gridGlobal().neighbourDir(this->params()->rank(), ii);
            ScaFES::GridSub<DIM> memNormalNeighbour(
                this->gridGlobal().partition(idNeighbour).idxNodeFirstSub(),
                this->gridGlobal().partition(idNeighbour).idxNodeLastSub(),
                this->memAll());

            // The grid of the data field can be larger / smaller than the one
            // of the domain.
            // ==> The communication grid must be larger / smaller, too.
            memExtendedNeighbour =
                memNormalNeighbour.extend(dirNeighbour, this->stencilWidth());
            memCommCurr = memExtendedNeighbour.intersectWith(this->memNormal());

            for (std::size_t idxFace = 0; idxFace < 2 * DIM; ++idxFace)
            {
                memFaceCurr = memExtendedNeighbour.intersectWith(
                    this->memBorder().at(idxFace));

                if (memFaceCurr.valid())
                {
                    // Check if memFaceCurr is equal to the global
                    // dimension of the data field in at least one direction.
                    // If the data field ends within a gird partition,
                    // then a new grid will be created, BUT:
                    // The corresponding neighbour grid partition
                    // would NOT create a grid!
                    int idx = dirNeighbour / 2;
                    int rem = dirNeighbour % 2;

                    // LEFT, BACK, BOTTOM...
                    if ((0 == rem) &&
                        (this->gridGlobal().discreteDomain().idxNodeFirst().elem(idx) ==
                         memCommCurr.idxNodeFirstSub().elem(idx)))
                    {
                        std::cout << __FILE__ << " " << __LINE__
                                  << " myRank = " << this->params()->rank()
                                  << "\t end(df) reached: "
                                  << " this->gridGlobal().idxNodeFirst()= "
                                  << this->gridGlobal().discreteDomain().idxNodeFirst()
                                  << " memCommCurr.idxNodeFirstSub= "
                                  << memCommCurr.idxNodeFirstSub() << "\n";
                        break;
                    }
                    // RIGHT, FRONT, TOP...
                    else if ((1 == rem) &&
                             (this->gridGlobal().discreteDomain().idxNodeLast().elem(idx) + 1 ==
                              memCommCurr.idxNodeLastSub().elem(idx)))
                    {
                        std::cout << __FILE__ << " " << __LINE__
                                  << "\t end(df) reached: "
                                  << " this->gridGlobal().idxNodeLast()= "
                                  << this->gridGlobal().discreteDomain().idxNodeLast()
                                  << " memCommCurr.idxNodeLastSub= "
                                  << memCommCurr.idxNodeLastSub() << "\n";
                        break;
                    }
                    else
                    {
                        this->mMemComm.push_back(memCommCurr);
                        this->mValuesToExchange.push_back(ScaFES::Buffer<CT>(
                            this->params()->myWorld(), idNeighbour,
                            nColumns() * memCommCurr.nNodesSub().size(),
                            this->params()->useBoostMpiSkeletonConcept() ));
                        // Create ghost grids.
                        int idx;
                        ScaFES_IntNtuple idxNodeFirstGridGhost;
                        ScaFES_IntNtuple idxNodeLastGridGhost;
                        idxNodeFirstGridGhost = memCommCurr.idxNodeFirstSub();
                        idxNodeLastGridGhost = memCommCurr.idxNodeLastSub();

                        idx = dirNeighbour / 2;

                        if (0 == dirNeighbour % 2)
                        {
                            // Left, back, bottom.
                            idxNodeFirstGridGhost[idx] -= this->stencilWidth();
                            idxNodeLastGridGhost[idx] -= this->stencilWidth();
                        }
                        else
                        {
                            // Right, front, top.
                            idxNodeFirstGridGhost[idx] += this->stencilWidth();
                            idxNodeLastGridGhost[idx] += this->stencilWidth();
                        }

                        ScaFES::GridSub<DIM> memCurr(idxNodeFirstGridGhost,
                                                     idxNodeLastGridGhost,
                                                     this->memAll());
                        this->mMemGhost.push_back(memCurr);
                        break;
                    }
                }
            }
        }
    }
}

/*******************************************************************************
 * GETTER METHODS.
 ******************************************************************************/
template <typename CT, std::size_t DIM>
inline const std::string& DataField<CT, DIM>::name() const
{
    return this->mName;
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline ScaFES::Parameters* DataField<CT, DIM>::params() const
{
    return this->mParams;
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline const ScaFES::GridGlobal<DIM>& DataField<CT, DIM>::gridGlobal() const
{
    return this->mGG;
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline const int& DataField<CT, DIM>::stencilWidth() const
{
    return this->mStencilWidth;
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline const ScaFES::Grid<DIM>& DataField<CT, DIM>::memAll() const
{
    return this->mMemAll;
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline const ScaFES::GridSub<DIM>& DataField<CT, DIM>::memNormal() const
{
    return this->mMemNormal;
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline const int& DataField<CT, DIM>::nLayers() const
{
    return this->mNlayers;
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline const ScaFES::GridSub<DIM>& DataField<CT, DIM>::memHole() const
{
    return this->mMemHole;
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline const std::vector<ScaFES::GridSub<DIM>>&
DataField<CT, DIM>::memBorder() const
{
    return this->mMemBorder;
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline const std::vector<ScaFES::GridSub<DIM>>&
DataField<CT, DIM>::memComm() const
{
    return this->mMemComm;
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline const std::vector<ScaFES::GridSub<DIM>>&
DataField<CT, DIM>::memGhost() const
{
    return this->mMemGhost;
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline const CT& DataField<CT, DIM>::elemData(const int& idx) const
{
    return this->mElemData[idx];
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline CT* DataField<CT, DIM>::elemData() const
{
    return this->mElemData;
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline const int& DataField<CT, DIM>::nColumns() const
{
    return this->mNcolumns;
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline const int& DataField<CT, DIM>::nRows() const
{
    return this->mNrows;
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline int DataField<CT, DIM>::nElemsAllocated() const
{
    return (this->nColumns() * this->nRows());
}

/*******************************************************************************
 * SETTER METHODS.
 ******************************************************************************/
template <typename CT, std::size_t DIM>
inline CT& DataField<CT, DIM>::time() const
{
    // Time is stored at last component of vector.
    // [ x0 x1 x2 | x3 x4 x5 | x6 x7 x8 | t | 0 | 0 ]
    return (this->mElemData[this->nElemsAllocated() - this->nColumns()]);
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline CT& DataField<CT, DIM>::
operator()(const unsigned long int& memoryPos) const
{
    return this->mElemData[memoryPos];
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline CT& DataField<CT, DIM>::operator()(const ScaFES_IntNtuple& idxNode,
                                          const int& comp) const
{
    unsigned long int memoryPos = 0L;
    this->idxNode2MemoryPos(memoryPos, idxNode, comp);
    return this->mElemData[memoryPos];
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline CT& DataField<CT, DIM>::operator()(const int& ii, const int& jj,
                                          const int& kk, const int& comp) const
{
    static_assert((DIM == 3), "Method is valid in case DIM=3, only.");
    ScaFES_IntNtuple idxNode(ii, jj, kk);
    unsigned long int memoryPos = 0;
    this->idxNode2MemoryPos(memoryPos, idxNode, comp);
    return this->mElemData[memoryPos];
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline CT& DataField<CT, DIM>::operator()(const int& ii, const int& jj,
                                          const int& comp) const
{
    static_assert((DIM == 2), "Method is valid in case DIM=2, only.");
    ScaFES_IntNtuple idxNode(ii, jj);
    unsigned long int memoryPos = 0;
    this->idxNode2MemoryPos(memoryPos, idxNode, comp);
    return this->mElemData[memoryPos];
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline CT& DataField<CT, DIM>::at(const unsigned long int& memoryPos) const
{
    assert(memoryPos <= this->mMemoryPosGridAllLast);
    assert(this->mMemoryPosGridAllFirst <= memoryPos);
    // throw(out-of-range)
    this->checkRange(memoryPos);
    return this->mElemData[memoryPos];
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline CT& DataField<CT, DIM>::at(const ScaFES_IntNtuple& idxNode,
                                  const int& comp) const
{
    assert(idxNode <= this->memAll().idxNodeLast());
    assert(this->memAll().idxNodeFirst() <= idxNode);
    // throw(out-of-range)
    unsigned long int memoryPos = 0;
    this->idxNode2MemoryPos(memoryPos, idxNode, comp);
    this->checkRange(memoryPos);
    return this->mElemData[memoryPos];
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline CT& DataField<CT, DIM>::at(const int& ii, const int& jj, const int& kk,
                                  const int& comp) const
{
#ifdef SCAFES_HAVE_BOOST
    BOOST_STATIC_ASSERT_MSG(
        (3 == DIM), "Please use ScaFES::Ntuple objects for data access.");
#endif
    ScaFES_IntNtuple idxNode(ii, jj, kk);
    assert(this->memAll().idxNodeFirst() <= idxNode);
    assert(idxNode <= this->memAll().idxNodeLast());
    unsigned long int memoryPos = 0;
    this->idxNode2MemoryPos(memoryPos, idxNode, comp);
    return this->mElemData[memoryPos];
}

/*******************************************************************************
 * COMPARISON METHODS.
 ******************************************************************************/
template <typename CT, std::size_t DIM>
inline bool DataField<CT, DIM>::
operator==(const ScaFES::DataField<CT, DIM>& second)
{
    double eps = 1E-10;
    return this->hasSameValues(second, eps);
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline bool
DataField<CT, DIM>::hasSameSize(const ScaFES::DataField<CT, DIM>& rhs) const
{
    if (this->stencilWidth() != rhs.stencilWidth())
    {
        return false;
    }

    if (this->memAll() != rhs.memAll())
    {
        return false;
    }

    return true;
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline bool
DataField<CT, DIM>::hasSameValues(const ScaFES::DataField<CT, DIM>& second,
                                  const double& eps) const
{
    if (this->hasSameSize(second))
    {
        bool isSame = true;
        ScaFES_IntNtuple idxNode;
        double tmpDiff = 0.0;

        //#pragma omp parallel for
        for (typename ScaFES::GridSub<DIM>::iterator
                 it = second.memNormal().begin(),
                 et = second.memNormal().end();
             it < et; ++it)
        {
            idxNode = it.idxNode();
            tmpDiff = ::fabs(df1(idxNode) - second(idxNode));

            if (tmpDiff > eps)
            {
                isSame = false;
                break;
            }
        }

        if (!isSame)
        {
            return false;
        }

        return true;
    }

    return false;
}

/*******************************************************************************
 * WORK METHODS.
 ******************************************************************************/
template <typename CT, std::size_t DIM>
inline void
DataField<CT, DIM>::checkRange(const unsigned long int& memoryPos) const
{
    assert(memoryPos <= this->mMemoryPosGridAllLast);
    assert(this->mMemoryPosGridAllFirst <= memoryPos);
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline void
DataField<CT, DIM>::swapPointerToMemoryWith(ScaFES::DataField<CT, DIM>& df)
{
    if (this->params()->rankOutput() == this->params()->rank() &&
        (0 < this->params()->indentDepth()))
    {
        std::cout << mParams->getPrefix()
                  << " * Swap pointers to memory of data fields "
                  << this->name() << " and " << df.name() << "..." << std::endl;
    }
    this->mParams->decreaseLevel();

    std::swap(this->mElemData, df.mElemData);

    this->mParams->increaseLevel();
    if (this->params()->rankOutput() == this->params()->rank() &&
        (0 < this->params()->indentDepth()))
    {
        std::cout << mParams->getPrefix()
                  << "   Swapped."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline ScaFES::DataField<CT, DIM>& DataField<CT, DIM>::operator=(const CT& val)
{
    for (int ii = 0; ii < this->nElemsAllocated(); ++ii)
    {
        this->mElemData[ii] = val;
    }
    return *this;
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline void
DataField<CT, DIM>::assignValues(const ScaFES::DataField<CT, DIM>& from)
{
    if (this->hasSameSize(from))
    {
        for (int ii = 0; ii < this->nElemsAllocated(); ++ii)
        {
            this->mElemData[ii] = from.elemData(ii);
        }
    }
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline void DataField<CT, DIM>::idxNodeScalar2MemoryPos(
    unsigned long int& memoryPos, const unsigned long int& idxNodeScalar,
    const int& comp) const
{
    memoryPos = (this->nColumns() * idxNodeScalar + comp);
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline void DataField<CT, DIM>::getMemoryPosWithHole(
    unsigned long int& memoryPosWithHole,
    const unsigned long int& memoryPos) const
{
    if (memoryPos > this->mMemoryPosGridHoleLast)
    {
        memoryPosWithHole = memoryPos - this->memHole().nNodesSub().size();
    }
    else if (memoryPos < this->mMemoryPosGridHoleLast)
    {
        memoryPosWithHole = memoryPos;
    }
    else
    {
        throw std::runtime_error("DataField has hole at memoryPos.");
    }
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline void
DataField<CT, DIM>::idxNode2MemoryPos(unsigned long int& memoryPos,
                                      const ScaFES_IntNtuple& idxNode,
                                      const int& comp) const
{
    unsigned long int memoryPosWithoutHole = 0L;
    unsigned long int idxNodeScalar = 0L;
    this->memAll().idxNodeTuple2Scalar(idxNodeScalar, idxNode);
    this->idxNodeScalar2MemoryPos(memoryPosWithoutHole, idxNodeScalar, comp);
    //     this->getMemoryPosWithHole(memoryPos, memoryPosWithoutHole);
    memoryPos = memoryPosWithoutHole;
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline void DataField<CT, DIM>::getMemoryPosWithoutHole(
    unsigned long int& memoryPos,
    const unsigned long int& memoryPosWithHole) const
{
    if (memoryPosWithHole > this->mMemoryPosGridHoleLast)
    {
        memoryPos = memoryPosWithHole + this->mMemHole.nNodesSub().size();
    }
    else if (memoryPosWithHole < this->mMemoryPosGridHoleLast)
    {
        memoryPos = memoryPosWithHole;
    }
    else
    {
        throw std::runtime_error("DataField has hole at memoryPos.");
    }
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline void DataField<CT, DIM>::memoryPos2IdxNodeScalar(
    unsigned long int& idxNodeScalar, int& comp,
    const unsigned long int& memoryPos) const
{
    comp = memoryPos % this->nColumns();
    idxNodeScalar = memoryPos / this->nColumns();
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline void
DataField<CT, DIM>::memoryPos2IdxNode(ScaFES_IntNtuple& idxNode, int& comp,
                                      const unsigned long int& memoryPos) const
{
    unsigned long int memoryPosWithoutHole = 0L;
    // this->getMemoryPosWithoutHole(memoryPosWithoutHole, memoryPos);
    memoryPosWithoutHole = memoryPos;
    unsigned long int idxNodeScalar = 0L;
    this->memoryPos2IdxNodeScalar(idxNodeScalar, comp, memoryPosWithoutHole);
    this->memAll().idxNodeScalar2Tuple(idxNode, idxNodeScalar);
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline std::size_t DataField<CT, DIM>::elementSize() const
{
    // Cast from size_t to int.
    return sizeof(CT);
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline ScaFES::DataField<CT, DIM>& DataField<CT, DIM>::
operator+=(const ScaFES::DataField<CT, DIM>& rhs)
{
    unsigned long int idxScalarNode;

    for (typename ScaFES::GridSub<DIM>::iterator it = this->mMemNormal.begin(),
                                                 et = this->mMemNormal.end();
         it < et; ++it)
    {
        idxScalarNode = it.idxScalarNode();
        for (int kk = 0; kk < this->nColumns(); ++kk)
        {
            this->mElemData[idxScalarNode + kk] +=
                rhs.elemData(idxScalarNode + kk);
        }
    }

    return *this;
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline ScaFES::DataField<CT, DIM>& DataField<CT, DIM>::
operator-=(const ScaFES::DataField<CT, DIM>& rhs)
{
    unsigned long int idxScalarNode;

    for (typename ScaFES::GridSub<DIM>::iterator it = this->mMemNormal.begin(),
                                                 et = this->mMemNormal.end();
         it < et; ++it)
    {
        idxScalarNode = it.idxScalarNode();
        for (int kk = 0; kk < this->nColumns(); ++kk)
        {
            this->mElemData[idxScalarNode + kk] -=
                rhs.elemData(idxScalarNode + kk);
        }
    }

    return *this;
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline ScaFES::DataField<CT, DIM>& DataField<CT, DIM>::operator+=(const CT& rhs)
{
    unsigned long int idxScalarNode;

    for (typename ScaFES::GridSub<DIM>::iterator it = this->mMemNormal.begin(),
                                                 et = this->mMemNormal.end();
         it < et; ++it)
    {
        idxScalarNode = it.idxScalarNode();
        for (int kk = 0; kk < this->nColumns(); ++kk)
        {
            this->mElemData[idxScalarNode + kk] += rhs;
        }
    }

    return *this;
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline ScaFES::DataField<CT, DIM>& DataField<CT, DIM>::operator-=(const CT& rhs)
{
    unsigned long int idxScalarNode;

    for (typename ScaFES::GridSub<DIM>::iterator it = this->mMemNormal.begin(),
                                                 et = this->mMemNormal.end();
         it < et; ++it)
    {
        idxScalarNode = it.idxScalarNode();
        for (int kk = 0; kk < this->nColumns(); ++kk)
        {
            this->mElemData[idxScalarNode + kk] -= rhs;
        }
    }

    return *this;
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline double DataField<CT, DIM>::normLinf() const
{
    double res = 0.0;

    for (typename ScaFES::GridSub<DIM>::iterator it = this->mMemNormal.begin(),
                                                 et = this->mMemNormal.end();
         it < et; ++it)
    {
        for (int kk = 0; kk < this->nColumns(); ++kk)
        {
            if (res < fabs(this->elemData(it.idxScalarNode() + kk)))
            {
                res = this->elemData(it.idxScalarNode() + kk);
            }
        }
    }

    return res;
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
void
DataField<CT, DIM>::setValuesAtMemNormal(const ScaFES::DataField<CT, DIM>& df)
{
    if ((this->params()->rankOutput() == this->params()->rank()) &&
        (0 < this->params()->indentDepth()))
    {
        std::cout << mParams->getPrefix()
                  << " * Set values of " << df.name() << " to " << this->name()
                  << " at memNormal... " << std::endl;
    }
    this->mParams->decreaseLevel();

    unsigned long int pos = 0L;
    for (typename ScaFES::GridSub<DIM>::iterator it = df.memNormal().begin(),
                                                 et = df.memNormal().end();
         it < et; ++it)
    {
        for (int kk = 0; kk < this->nColumns(); ++kk)
        {
            // Compute corresponding memory position of current node number
            // of normal grid in this data field,
            this->idxNode2MemoryPos(pos, it.idxNode(), kk);
            this->mElemData[pos] = df(it.idxNode());
        }
    }

    this->mParams->increaseLevel();
    if ((this->params()->rankOutput() == this->params()->rank()) &&
        (0 < this->params()->indentDepth()))
    {
        std::cout << mParams->getPrefix()
                  << "   Set."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
void DataField<CT, DIM>::setValuesAtMemAll(const ScaFES::DataField<CT, DIM>& df)
{
    if ((this->params()->rankOutput() == this->params()->rank()) &&
        (0 < this->params()->indentDepth()))
    {
        std::cout << mParams->getPrefix()
                  << " * Set values of " << df.name() << " to " << this->name()
                  << " at memAll... " << std::endl;
    }
    this->mParams->decreaseLevel();

    memcpy(this->mElemData, df.elemData(),
           (this->nElemsAllocated() * sizeof(this->elemData())));

    this->mParams->increaseLevel();
    if ((this->params()->rankOutput() == this->params()->rank()) &&
        (0 < this->params()->indentDepth()))
    {
        std::cout << mParams->getPrefix()
                  << "   Set."
                  << std::endl;
    }
}

/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
void DataField<CT, DIM>::check(const DataField<CT, DIM>& df)
{
    double eps = 1E-10;
    if (this->hasSameValues(df, eps))
    {
        if (this->params()->rankOutput() == this->params()->rank() &&
            (0 < this->params()->indentDepth()))
        {
            std::cout << mParams->getPrefix()
                      << ": Check of (" << this->name() << ", " << df.name()
                      << ")  SUCCESSFUL." << std::endl;
        }
    }
    else
    {
        if (this->params()->rankOutput() == this->params()->rank() &&
            (0 < this->params()->indentDepth()))
        {
            std::cout << mParams->getPrefix()
                      << ": Check of (" << this->name() << ", " << df.name()
                      << ")  FAILED." << std::endl;
        }
    }
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
double DataField<CT, DIM>::compErrLinf(const ScaFES::DataField<CT, DIM>& df)
{
    if ((this->params()->rankOutput() == this->params()->rank()) &&
        (0 < this->params()->indentDepth()))
    {
        std::cout << mParams->getPrefix()
                  << " * Compute Linf error(" << name() << ", " << df.name()
                  << ")... " << std::endl;
    }
    this->mParams->decreaseLevel();
    ScaFES_IntNtuple idxNode;
    ScaFES_IntNtuple idxNodeFirstSub =
        this->gridGlobal().partition(df.params()->rank()).idxNodeFirstSub();
    ScaFES_IntNtuple idxNodeLastSub =
        this->gridGlobal().partition(df.params()->rank()).idxNodeLastSub();

    int comp = 0;
    double errLinf = 0.0;
    bool isIdxNodeAtMemNormal = true;

    int nNodesTotal = this->memAll().nNodes().size();
    int ii;
    // * ii Iteration variable ---> private
    // * idxNode WRITE access Corresponding global node number ---> private
    // * df1 READ access --> shared
    // * df2 READ access --> shared
    // Parallelize the region, assign portions of the loop to each thread
    // and compute for each thread a local maximum.
    // At the end, compute the global maximum of all threads.

    // #ifdef _OPENMP
    // #pragma omp parallel for
    //     schedule(static)
    //     private(ii, idxNode, isIdxNodeAtMemNormal)
    //     shared(df, nNodesTotal, idxNodeFirstSub, idxNodeLastSub)
    //     reduction(max: errLinf)
    // #endif
    for (ii = 0; ii < nNodesTotal; ++ii)
    {
        this->memoryPos2IdxNode(idxNode, comp, ii);
        isIdxNodeAtMemNormal = true;
        for (std::size_t jj = 0; jj < DIM; ++jj)
        {
            if ((idxNodeFirstSub.elem(jj) > idxNode.elem(jj)) ||
                (idxNodeLastSub.elem(jj) < idxNode.elem(jj)))
            {
                isIdxNodeAtMemNormal = false;
            }
        }

        if (isIdxNodeAtMemNormal)
        {
            for (int kk = 0; kk < this->nColumns(); ++kk)
            {
                /* For debugging purposes, only. */
                /*
                std::cout << mParams->getPrefix()
                          << " * compare fields ( "
                          << name() << " " << df.name() << "): "
                          << std::scientific << this->elemData(ii + kk)
                          << "  " << this->elemData(ii + kk)  << "  "
                          << ::fabs(this->elemData(ii + kk) - df(idxNode, kk))
                          << std::fixed
                          << std::endl;
                */
                if (::fabs(this->elemData(ii + kk) - df(idxNode, kk)) > errLinf)
                {
                    errLinf = ::fabs(this->elemData(ii + kk) - df(idxNode, kk));
                }
            }
        }
    }

    if ((this->params()->rankOutput() == this->params()->rank()) &&
        (0 < this->params()->indentDepth()))
    {
        std::cout << mParams->getPrefix()
                  << " * Local Linf error(" << name() << ", " << df.name()
                  << "): " << std::scientific << errLinf << std::fixed
                  << std::endl;
    }
    this->mParams->increaseLevel();
    if ((this->params()->rankOutput() == this->params()->rank()) &&
        (0 < this->params()->indentDepth()))
    {
        std::cout << mParams->getPrefix()
                  << "   Computed."
                  << std::endl;
    }
    return errLinf;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM>
void DataField<TT, DIM>::exchangeValuesInBuffers(const int& /*timeIter*/)
{
    if (this->params()->rankOutput() == this->params()->rank() &&
        (0 < this->params()->indentDepth()))
    {
        std::cout << mParams->getPrefix()
                  << " * Exchange buffer values of data field " << this->name()
                  << "..." << std::endl;
    }
    this->mParams->decreaseLevel();
    for (std::size_t ii = 0; ii < this->mValuesToExchange.size(); ++ii)
    {
        this->mValuesToExchange[ii].sendValues();
        this->mValuesToExchange[ii].receiveValues();
    }
    this->mParams->increaseLevel();
    if (this->params()->rankOutput() == this->params()->rank() &&
        (0 < this->params()->indentDepth()))
    {
        std::cout << mParams->getPrefix()
                  << "   Exchanged."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM> void DataField<TT, DIM>::waitAll()
{
    for (std::size_t ii = 0; ii < this->mValuesToExchange.size(); ++ii)
    {
        this->mValuesToExchange[ii].waitAll();
    }
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM>
void
DataField<TT, DIM>::copyValuesFromMemCommToSendBuffer(const int& /*timeIter*/
                                                      )
{
#ifdef VTRACE
    VT_TRACER("Df::copyValuesFromMemCommToSendBuffer()");
#endif
    unsigned long int posBuffer;

    if (this->params()->rankOutput() == this->params()->rank() &&
        (0 < this->params()->indentDepth()))
    {
        std::cout << mParams->getPrefix()
                  << " * Copy values from mem comm. to send buffer of df "
                  << this->name() << "..." << std::endl;
    }
    this->mParams->decreaseLevel();
    for (std::size_t ii = 0; ii < memComm().size(); ++ii)
    {
        posBuffer = 0;

        for (typename ScaFES::GridSub<DIM>::iterator
                 it = memComm().at(ii).begin(),
                 et = memComm().at(ii).end();
             it < et; ++it)
        {
            for (int kk = 0; kk < nColumns(); ++kk)
            {
                unsigned long int memPos = 0L;
                this->idxNode2MemoryPos(memPos, it.idxNode(), kk);
                this->mValuesToExchange.at(ii)
                    .setElemData(posBuffer, this->elemData(memPos));
                ++posBuffer;
            }
        }
    }
    this->mParams->increaseLevel();
    if (this->params()->rankOutput() == this->params()->rank() &&
        (0 < this->params()->indentDepth()))
    {
        std::cout << mParams->getPrefix()
                  << "   Copied."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM>
void
DataField<TT, DIM>::copyValuesFromMemGhostToSendBuffer(const int& /*timeIter*/
                                                       )
{
#ifdef VTRACE
    VT_TRACER("Df::copyValuesFromMemGhostToSendBuffer()");
#endif
    unsigned long int posBuffer;

    if (this->params()->rankOutput() == this->params()->rank() &&
        (0 < this->params()->indentDepth()))
    {
        std::cout << mParams->getPrefix()
                  << " * Copy values from mem ghost to send buffer of df "
                  << this->name() << "..." << std::endl;
    }
    this->mParams->decreaseLevel();
    //#pragma omp parallel for
    for (std::size_t ii = 0; ii < memGhost().size(); ++ii)
    {
        posBuffer = 0;

        for (typename ScaFES::GridSub<DIM>::iterator
                 it = memGhost().at(ii).begin(),
                 et = memGhost().at(ii).end();
             it < et; ++it)
        {
            for (int kk = 0; kk < this->nColumns(); ++kk)
            {
                unsigned long int memPos = 0L;
                this->idxNode2MemoryPos(memPos, it.idxNode(), kk);
                this->mValuesToExchange.at(ii)
                    .setElemData(posBuffer, this->elemData(memPos));
                ++posBuffer;
            }
        }
    }
    this->mParams->increaseLevel();
    if (this->params()->rankOutput() == this->params()->rank() &&
        (0 < this->params()->indentDepth()))
    {
        std::cout << mParams->getPrefix()
                  << "   Copied."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM>
void DataField<TT, DIM>::copyValuesFromReceiveBufferToMemGhost(
    const int& /*timeIter*/
    )
{
#ifdef VTRACE
    VT_TRACER("Df::copyValuesFromReceiveBufferToMemGhost()");
#endif
    unsigned long int posBuffer;

    if (this->params()->rankOutput() == this->params()->rank() &&
        (0 < this->params()->indentDepth()))
    {
        std::cout << mParams->getPrefix()
                  << " * Copy values from receive buffer to mem ghost of df "
                  << this->name() << "..." << std::endl;
    }
    this->mParams->decreaseLevel();

    //#pragma omp parallel for
    for (std::size_t ii = 0; ii < memGhost().size(); ++ii)
    {
        posBuffer = 0;

        for (typename ScaFES::GridSub<DIM>::iterator
                 it = memGhost().at(ii).begin(),
                 et = memGhost().at(ii).end();
             it < et; ++it)
        {
            for (int kk = 0; kk < nColumns(); ++kk)
            {
                unsigned long int memPos = 0L;
                this->idxNode2MemoryPos(memPos, it.idxNode(), kk);
                this->mElemData[memPos] =
                    this->mValuesToExchange.at(ii).elemData(posBuffer);
                ++posBuffer;
            }
        }
    }
    this->mParams->increaseLevel();
    if (this->params()->rankOutput() == this->params()->rank() &&
        (0 < this->params()->indentDepth()))
    {
        std::cout << mParams->getPrefix()
                  << "   Copied."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM>
void
DataField<TT, DIM>::copyValuesFromReceiveBufferToMemComm(const int& /*timeIter*/
                                                         )
{
#ifdef VTRACE
    VT_TRACER("Df::copyValuesFromReceiveBufferToMemComm()");
#endif
    unsigned long int posBuffer;

    if (this->params()->rankOutput() == this->params()->rank() &&
        (0 < this->params()->indentDepth()))
    {
        std::cout << mParams->getPrefix()
                  << " * Copy values from receive buffer to mem comm. of df "
                  << this->name() << "..." << std::endl;
    }
    this->mParams->decreaseLevel();

    //#pragma omp parallel for
    for (std::size_t ii = 0; ii < memComm().size(); ++ii)
    {
        posBuffer = 0;

        for (typename ScaFES::GridSub<DIM>::iterator
                 it = memComm().at(ii).begin(),
                 et = memComm().at(ii).end();
             it < et; ++it)
        {
            for (int kk = 0; kk < nColumns(); ++kk)
            {
                this->mElemData[it.idxScalarNode() + kk] =
                    this->mValuesToExchange.at(ii).elemData(posBuffer);
                ++posBuffer;
            }
        }
    }
    this->mParams->increaseLevel();
    if (this->params()->rankOutput() == this->params()->rank() &&
        (0 < this->params()->indentDepth()))
    {
        std::cout << mParams->getPrefix()
                  << "   Copied."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM>
void DataField<TT, DIM>::copyAndSumValuesFromReceiveBufferToMemComm(
    const int& /*timeIter*/
    )
{
#ifdef VTRACE
    VT_TRACER("Df::copyAndSumValuesFromReceiveBufferToMemComm()");
#endif
    unsigned long int posBuffer;

    if (this->params()->rankOutput() == this->params()->rank() &&
        (0 < this->params()->indentDepth()))
    {
        std::cout
            << mParams->getPrefix()
            << " * Copy and sum values from receive buffer to mem comm. of df "
            << this->name() << "..." << std::endl;
    }
    this->mParams->decreaseLevel();

    //#pragma omp parallel for
    for (std::size_t ii = 0; ii < memComm().size(); ++ii)
    {
        posBuffer = 0;

        for (typename ScaFES::GridSub<DIM>::iterator
                 it = memComm().at(ii).begin(),
                 et = memComm().at(ii).end();
             it < et; ++it)
        {
            for (int kk = 0; kk < nColumns(); ++kk)
            {
                this->mElemData[it.idxScalarNode() + kk] +=
                    this->mValuesToExchange.at(ii).elemData(posBuffer);
                ++posBuffer;
            }
        }
    }
    this->mParams->increaseLevel();
    if (this->params()->rankOutput() == this->params()->rank() &&
        (0 < this->params()->indentDepth()))
    {
        std::cout << mParams->getPrefix()
                  << "   Copied."
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM>
void DataField<TT, DIM>::sync(const int& timeIter)
{
#ifdef VTRACE
    VT_TRACER("Df::sync()");
#endif

    // For the following special cases, no communication is necessary.
    // Non-parallel execution of ScaFES application.
    if (1 == this->mParams->myWorld().size())
    {
        return;
    }
    // This DataField has no ghost cells.
    if (0 == this->stencilWidth())
    {
        return;
    }

    copyValuesFromMemCommToSendBuffer(timeIter);
    exchangeValuesInBuffers(timeIter);
    waitAll();
    copyValuesFromReceiveBufferToMemGhost(timeIter);
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM>
void DataField<TT, DIM>::collectValuesAtMemComm(const int& timeIter)
{
#ifdef VTRACE
    VT_TRACER("Df::collectValuesAtMemComm()");
#endif

    // For the following special cases, no communication is necessary.
    if (1 == this->mParams->myWorld().size())
    {
        return;
    }

    if (0 == this->stencilWidth())
    {
        return;
    }

    copyValuesFromMemGhostToSendBuffer(timeIter);
    exchangeValuesInBuffers(timeIter);
    waitAll();
    copyAndSumValuesFromReceiveBufferToMemComm(timeIter);
}

/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM>
void DataField<TT, DIM>::write(const int& timeIter)
{
    // Check if file should be written at current time step.
    int intervalOfWrites = this->params()->nTimesteps();

    if (0 < this->params()->nSnapshots())
    {
        if (this->params()->nSnapshots() <= this->params()->nTimesteps())
        {
            intervalOfWrites =
                this->params()->nTimesteps() / this->params()->nSnapshots();
        }
    }
    if (this->params()->rankOutput() == this->params()->rank() &&
        (0 < this->params()->indentDepth()))
    {
        std::cout << mParams->getPrefix()
                  << " * Write data field " << name() << " to NetCDF file..."
                  << std::endl;
    }
    this->mParams->decreaseLevel();

    std::vector<bool> writeData(1, false);
    switch (mWriteToFile)
    {
    case WriteHowOften::NEVER:
        break;
    case WriteHowOften::AT_START:
        if (0 == timeIter)
        {
            writeData.at(0) = true;
        }
        break;
    case WriteHowOften::AT_START_AND_END:
        if (0 == timeIter || this->params()->nTimesteps() == timeIter)
        {
            writeData.at(0) = true;
        }
        break;
    case WriteHowOften::LIKE_GIVEN_AT_CL:
        if (0 == timeIter % intervalOfWrites)
        {
            writeData.at(0) = true;
        }
        break;
    case WriteHowOften::ALWAYS:
    default:
        writeData.at(0) = true;
        break;
    }
    if (this->params()->rankOutput() == this->params()->rank() &&
        (0 < this->params()->indentDepth()))
    {
        std::cout << mParams->getPrefix()
                  << " * Write status of " << this->name() << ": "
                  << writeData.at(0) << std::endl;
    }

    std::vector<TT*> tmpElemData(1, static_cast<TT*>(0));
    tmpElemData[0] = this->elemData();
    if (writeData.at(0))
    {
        std::vector<bool> writeToFile(1, true);
        mOutput.write(tmpElemData, writeToFile, timeIter);
    }

    this->mParams->increaseLevel();
    if (this->params()->rankOutput() == this->params()->rank() &&
        (0 < this->params()->indentDepth()))
    {
        std::cout << mParams->getPrefix()
                  << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM>
void DataField<TT, DIM>::read(const int& timeIter)
{
    // Check if file should be written at current time step.
    int intervalOfWrites = this->params()->nTimesteps();

    if (0 < this->params()->nSnapshots())
    {
        if (this->params()->nSnapshots() <= this->params()->nTimesteps())
        {
            intervalOfWrites =
                this->params()->nTimesteps() / this->params()->nSnapshots();
        }
    }

    std::vector<TT*> tmpElemData;
    tmpElemData.push_back(this->elemData());
    std::vector<ScaFES::GridSub<DIM>> tmpMemNormal;
    tmpMemNormal.push_back(this->memNormal());

    if (0 == timeIter % intervalOfWrites)
    {
        if (this->params()->rankOutput() == this->params()->rank() &&
            (0 < this->params()->indentDepth()))
        {
            std::cout << mParams->getPrefix()
                      << ": * Read data field " << name()
                      << "from NetCDF file..." << std::endl;
        }
        this->mParams->decreaseLevel();

        mOutput.read(tmpElemData, timeIter);

        this->mParams->increaseLevel();
        if (this->params()->rankOutput() == this->params()->rank() &&
            (0 < this->params()->indentDepth()))
        {
            std::cout << mParams->getPrefix() << ":   Read."
                      << std::endl;
        }
    }
}

/*******************************************************************************
 * FREE METHODS.
 ******************************************************************************/
template <typename CT, std::size_t DIM>
inline std::ostream& operator<<(std::ostream& output,
                                const ScaFES::DataField<CT, DIM>& df)
{
    //#pragma omp parallel for
    for (typename ScaFES::Grid<DIM>::iterator it = df.memAll().begin(),
                                              et = df.memAll().end();
         it < et; ++it)
    {
        for (int kk = 0; kk < df.nColumns(); ++kk)
        {
            output << it.idxScalarNode() + kk << "   " << it.idxNode() << "   "
                   << df(it.idxScalarNode() + kk) << std::endl;
        }
    }

    return output;
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline const ScaFES::DataField<CT, DIM>
operator+(const ScaFES::DataField<CT, DIM>& lhs,
          const ScaFES::DataField<CT, DIM>& rhs)
{
    ScaFES::DataField<CT, DIM> tmp(lhs); // Creates a copy of the left operand.
    tmp += rhs;                          // Implementation via the += operator.
    return tmp;
}
/*----------------------------------------------------------------------------*/
template <typename CT, std::size_t DIM>
inline const ScaFES::DataField<CT, DIM>
operator-(const ScaFES::DataField<CT, DIM>& lhs,
          const ScaFES::DataField<CT, DIM>& rhs)
{
    ScaFES::DataField<CT, DIM> tmp(lhs); // Creates a copy of the left operand.
    tmp -= rhs;                          // Implementation via the -= operator.
    return tmp;
}

} // End of namespace. //

#endif
