/* ScaFES
 * Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file ScaFES_GridGlobal.hpp
 *  @brief Contains the class template GridGlobal.
 */

#ifndef SCAFES_GRIDGLOBAL_HPP_
#define SCAFES_GRIDGLOBAL_HPP_

#include "ScaFES_Config.hpp"

//#include <iomanip>
#include <stdexcept>
#include <vector>

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
// #include <boost/serialization/version.hpp>
// #include <boost/serialization/pfto.hpp>
// #include <boost/serialization/vector.hpp>
// #include <boost/serialization/base_object.hpp>
// #endif

#include "ScaFES_Parameters.hpp"
#include "ScaFES_Timer.hpp"
#include "ScaFES_Ntuple.hpp"
#include "ScaFES_Grid.hpp"
#include "ScaFES_GridSub.hpp"

namespace ScaFES
{

/*******************************************************************************
 ******************************************************************************/
/**
* \class GridGlobal
* @brief The class \c GridGlobal represents a global grid
* of dimension \c DIM, where \c DIM is the template parameter of the class
* template.
*
* \n Each grid partition gets a unique identifier.
* Currently, this identifier will be equal to the MPI process identifier,
* i.e. the MPI process with id (=rank) k will work on grid partition k.
*
* \n The global grid will be divided into a given number of grid partitions.
* This can be done using
* - a RCB algorithm. This version of the algorithm delivers the minimal
*   number of cutting faces.
* - a uniform decomposition of the computational domain.
*
* \n The class contains two grids:
* - the discretized computational domain and
* - the partitioned discretized computational domain.
*
* \n The order of the direct neighbours within in the connectivity matrix
* is as follows:
* LEFT   0
* RIGHT  1
* BACK   2
* FRONT  3
* BOTTOM 4
* TOP    5
*
* \n The global grid contains information about the global computational
* domain and the number of partitions, into which the global grid is divided.
* The computational domain will be divided into sub domains. This process
* will be usually called domain decomposition.
* \n The partitioning can be controlled via the command line option
* \code divideGrid(a,b,c) \endcode with a,b,c in {0,1}.
* If a component of \c divideGrid is set to 1, cutting the computational domain
* into this direction is allowed. Analogously, if a component of \c divideGrid
* is set to 0, then cutting into this direction is forbidden.
* \remarks
* Currently, the domain decomposition works not in parallel.
* <ol>
* <li>
*  Grid partitions are created and stored in the file indicated.
*  The decomposition of the domain uses an algorithm that delivers minimal
*  cutting faces. If one do not want to cut in a certain dimension, one has
*  to disable the cutting via the variable
*  \code Index divideGrid(a,b,c) \endcode
* with \code a,b,c \in {0,1}'. \endcode
*  The value 0 disables cutting in this direction.
* The numerical grid can be given via the command line parameter
*  '--nNodes=width x depth x height.
* This will create a grid with w*d*h nodes.
* </li>
* </ol>
* Case: Domain decomposition will be done uniformly in each direction:
* Grid partitions p_k are numbered lexicographically:
*    *-------*--------*
*    |  p_2  |  p_3   |
*    |       |        |
*    *-------*--------*
*    |  p_0  |  p_1   |
*    |       |        |
*    *-------*--------*
*
* \remarks
* - As the grid partitions will be sorted by an RCB algorithm,
*   the order of the grid partitions (i.e. their identifiers) is not known,
*   and therefore, the grid partitions have to be stored in a member variable
*   (vector). Thus, the neighbour ids and neighbour directions have to be
*   stored in member variables, too.
* - If the order of the grid partitions is known (see uniform decomposition
*   of the global grid, e.g.), then the partitions, the neighbours ids, and the
*   neighbour directions could be computed ONLINE via calling appropiate
*   methods in order to save memory.
*
* \n Assumption:
* The vector containing the number of grid nodes in each direction
* has to be a multiple (elementwise) of the vector containing the number of
* partitions:
* nNodes(i) = alpha_i * nPartitions(i) with alpha_i != 0 integers for all i.
*/
template <std::size_t DIM = 3>
class GridGlobal
{
public:
    /*--------------------------------------------------------------------------
    | ENUMERATIONS.
    --------------------------------------------------------------------------*/
    /** Possible types of a grid node. */
    enum scafes_gridnode_type
    {
        SCAFES_GLOBAL_ALL_NODE = 1,
        SCAFES_GLOBAL_NORMAL_NODE = 2,
        SCAFES_GLOBAL_INNER_NODE = 4,
        SCAFES_GLOBAL_BORDER_NODE = 8,
        SCAFES_GLOBAL_GHOST_NODE = 16,
        SCAFES_GLOBAL_COMM_NODE = 32,
        SCAFES_PARTITITON_ALL_NODE = 64,
        SCAFES_PARTITION_INNER_NODE = 128,
        SCAFES_PARTITION_BORDER_NODE = 256,
        SCAFES_PARTITION_GHOST_NODE = 512,
        SCAFES_PARTITION_COMM_NODE = 1024
    };

    /*--------------------------------------------------------------------------
    | LIFE CYCLE METHODS.
    --------------------------------------------------------------------------*/
    /** Creates own default constructor. */
    GridGlobal();

    /** Creates own constructor.
     * \param params All parameters of ScaFES.
     */
    GridGlobal(ScaFES::Parameters const& params);

    /** Creates another constructor.
     * \param nPartitions Number of grid partitions in each dimension.
     * \param nNodes Number of grid nodes in each dimension.
     * \param coordNodeFirst Coordinates of the first node.
     * \param coordNodeLast Coordinates of the last node.
     * \param divideGridIntoDir Responsible for cutting the grid.
     */
    GridGlobal(ScaFES::Ntuple<int, DIM> const& nPartitions,
               ScaFES::Ntuple<int, DIM> const& nNodes,
               ScaFES::Ntuple<double, DIM> const& coordNodeFirst,
               ScaFES::Ntuple<double, DIM> const& coordNodeLast,
               ScaFES::Ntuple<bool, DIM> const& divideGridIntoDir);

    /** Creates the compiler provided copy constructor. */
    GridGlobal(ScaFES::GridGlobal<DIM> const& rhs);

    /** Creates the compiler provided assignment operator. */
    ScaFES::GridGlobal<DIM>& operator=(ScaFES::GridGlobal<DIM> from);

    /** Tests if two global grids are not equal. */
    bool operator!=(GridGlobal<DIM> const& rhs) const;

    /** Creates own destructor. */
    ~GridGlobal();

    /*--------------------------------------------------------------------------
    | GETTER METHODS.
    --------------------------------------------------------------------------*/
    /** Returns the number of partitions. */
    ScaFES::Ntuple<int, DIM> const& nPartitions() const;

    /** Returns the discretized computational (global) domain. */
    ScaFES::Grid<DIM> const& discreteDomain() const;

    /** Returns the decomposition of the discretized computational domain. */
    ScaFES::Grid<DIM> const& decompDomain() const;

    /** Returns the type of the domain decomposition. */
    std::string const& typeDomainDecomp() const;

    /** Returns the decision in which direction the grid
     * should be divided. */
    ScaFES::Ntuple<bool, DIM> const& divideGrid() const;

    /** Returns the vector of all partitions. */
    std::vector<ScaFES::GridSub<DIM>> const& partition() const;

    /** Returns the partition \c idx. */
    ScaFES::GridSub<DIM> const& partition(int const& ii) const;

    /** Returns the identifiers of all neighbours of all partitions. */
    std::vector<std::vector<int>> const& neighbourId() const;

    /** Returns the identifiers of all neighbours of partition \c ii. */
    std::vector<int> const& neighbourId(int const& ii) const;

    /** Returns the identifier of neighbour \c idx of partition \c ii. */
    int const& neighbourId(int const& ii, int const& idx) const;

    /** Returns the directions of all neighbours of all partitions. */
    std::vector<std::vector<int>> const& neighbourDir() const;

    /** Returns the directions of all neighbours of partition \c ii. */
    std::vector<int> const& neighbourDir(int const& ii) const;

    /** Returns the direction of neighbour \c idx of partition \c ii. */
    int const& neighbourDir(int const& ii, int const& idx) const;

    /*----------------------------------------------------------------------
    | COMPARISON METHODS.
    ----------------------------------------------------------------------*/
    /** Tests if two grids are equal. */
    bool operator==(GridGlobal<DIM> const& rhs) const;

    /*--------------------------------------------------------------------------
    | WORK METHODS.
    --------------------------------------------------------------------------*/
    /** Returns the number of total partitions. */
    int nPartitionsTotal() const;

    /*--------------------------------------------------------------------------
    | FREE FRIENDS METHODS.
    --------------------------------------------------------------------------*/
    /** Swaps two objects of this class. */
    template <std::size_t DD>
    friend void swap(ScaFES::GridGlobal<DD>& first,
                     ScaFES::GridGlobal<DD>& second);

    /** Prints all grid partitions of the global grid together with its
     * the identifiers and the directions of the direct neighboured grid
     * partitions.
     * \remark This method should be called by one process, only! */
    template <std::size_t DD>
    friend std::ostream& operator<<(std::ostream& output,
                                    ScaFES::GridGlobal<DD> const& gg);

private:
    /*--------------------------------------------------------------------------
    | WORK METHODS.
    --------------------------------------------------------------------------*/
    /** Divides the computational domain into a given number of partitions.
     * \param nPartitions Number of desired partitions in each direction.
     * \param divideGridIntoDir If a component is set to \c false,
     * the domain will not be partitioned into the
     * corresponding direction.
     */
    void partitionDomainUniformly();

    /** Divides the computational domain into partitions.
     * \param nPartitions Number of partitions.
     * \param idxNodeFirst Node number of first grid point.
     * \param nNodes Number of nodes of the grid.
     * \param divideGridIntoDir If a component is set to \c false,
     * the domain will not be partitioned into the
     * corresponding direction.
     * \remarks The algorithm works recursively.
     * It stops if the number of partitions is equal to one.
     */
    void
    partitionDomainUsingRCB(int const& nPartitions,
                            ScaFES::Ntuple<int, DIM> const& idxNodeFirst,
                            ScaFES::Ntuple<int, DIM> const& nNodes,
                            ScaFES::Ntuple<bool, DIM> const& divideGridIntoDir =
                                ScaFES::Ntuple<bool, DIM>(true));

    /** Searches all neighbours in a given direction.
     * \param direction Direction.
     */
    void searchNeighbours(int const& direction);

    /** Searches neighbours in all grid dimensions. */
    void searchNeighbours();

    /** Determines the identifiers and directions of all neighbours of all
     * partitions in case of uniform domain decomposition.
     */
    void searchNeighboursUniform();

    /** Returns the sum over all CPU loads.
     * \param sum Sum of all CPU loads within a grid along a coordinate.
     * \param direction Direction.
     * \param idxNodeFirst First node number.
     * \param nNodes Number of grid nodes in each dimension.
     */
    void sumCpuLoads(std::vector<float>& sum, int const& direction,
                     ScaFES::Ntuple<int, DIM> const& idxNodeFirst,
                     ScaFES::Ntuple<int, DIM> const& nNodes) const;

    /** Initialises dependent member variables. */
    void initDependentMembers();

    /*--------------------------------------------------------------------------
    | FRIEND CLASSES.
    --------------------------------------------------------------------------*/
#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
    /** Be friend for serializing this class. */
    friend class boost::serialization::access;

    /** Serializes the class. */
    template <class Archive>
    void serialize(Archive& ar, const unsigned int version);
#endif

    /*--------------------------------------------------------------------------
    | MEMBER VARIABLES.
    --------------------------------------------------------------------------*/
    /** Number of partitions in each dimension. */
    ScaFES::Ntuple<int, DIM> mNpartitions;

    /** Discretized computational domain = global grid. */
    ScaFES::Grid<DIM> mDiscreteDomain;

    /** Decomposition of discretized computational domain. */
    ScaFES::Grid<DIM> mDecompDomain;

    /** Type of domain decomposition. */
    std::string mTypeDomainDecomp;

    /** Divide grid in which directions? */
    ScaFES::Ntuple<bool, DIM> mDivideGridIntoDirs;

    /** Vector of all partitions. */
    std::vector<ScaFES::GridSub<DIM>> mPartitions;

    /** Identifiers of all neighbours of all partitions. */
    std::vector<std::vector<int>> mNeighbourIds;

    /** Directions of all neighbours of all partitions. */
    std::vector<std::vector<int>> mNeighbourDirs;
}; // End of class. //

/*******************************************************************************
 * FREE METHODS.
 ******************************************************************************/
/** Swaps two objects of the class \c GridGlobal. */
template <std::size_t DIM>
void swap(ScaFES::GridGlobal<DIM>& first, ScaFES::GridGlobal<DIM>& second);
/*----------------------------------------------------------------------------*/
/** Prints information about all grid partitions. */
template <std::size_t DIM>
std::ostream& operator<<(std::ostream& output,
                         ScaFES::GridGlobal<DIM> const& gg);

/*******************************************************************************
 * LIFE CYCLE METHODS.
 ******************************************************************************/
template <std::size_t DIM>
inline ScaFES::GridGlobal<DIM>::GridGlobal()
: mNpartitions()
, mDiscreteDomain()
, mDecompDomain()
, mTypeDomainDecomp("RCB")
, mDivideGridIntoDirs(ScaFES::Ntuple<bool, DIM>(false))
, mPartitions()
, mNeighbourIds()
, mNeighbourDirs()
{
    this->initDependentMembers();
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline GridGlobal<DIM>::GridGlobal(ScaFES::Parameters const& params)
: mNpartitions(params.nPartitions())
, mDiscreteDomain(ScaFES::Ntuple<int, DIM>(0),
                  ScaFES::Ntuple<int, DIM>(params.nNodes()) -
                      ScaFES::Ntuple<int, DIM>(1),
                  ScaFES::Ntuple<double, DIM>(params.coordNodeFirst()),
                  ScaFES::Ntuple<double, DIM>(params.coordNodeLast()))
, mDecompDomain(ScaFES::Ntuple<int, DIM>(0),
                ScaFES::Ntuple<int, DIM>(params.nPartitions()),
                ScaFES::Ntuple<double, DIM>(params.coordNodeFirst()),
                ScaFES::Ntuple<double, DIM>(params.coordNodeLast()))
, mTypeDomainDecomp(params.typeDomainDecomp())
, mDivideGridIntoDirs(ScaFES::Ntuple<bool, DIM>(params.divideGrid()))
, mPartitions()
, mNeighbourIds()
, mNeighbourDirs()
{
    if (params.rankOutput() == params.rank() && (0 < params.indentDepth()))
    {
        std::cout << "ScaFES : Rank=" << params.rankOutput() << "  * Partition the global grid..." << std::endl;
    }
    ScaFES::Timer timerTotal;
    params.myWorld().barrier();
    timerTotal.restart();
    this->initDependentMembers();
    params.myWorld().barrier();
    double tTotal = timerTotal.elapsed();
    if (params.rankOutput() == params.rank() && (0 < params.indentDepth()))
    {
        std::cout << std::scientific << std::fixed
                  << "ScaFES : Rank=" << params.rankOutput() << "  * time(gridGlobal):" << std::endl
                  << "ScaFES : Rank=" << params.rankOutput() << "    ----> t(gridGlobal) with barriers = " << std::setw(8)
                  << std::setprecision(5) << tTotal << " s. <-------"
                  << std::resetiosflags(std::ios::fixed | std::ios::showpoint)
                  << std::endl
                  << "ScaFES : Rank=" << params.rankOutput() << "     Partitioned." << std::endl;
    }
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline GridGlobal<DIM>::GridGlobal(
    ScaFES::Ntuple<int, DIM> const& nPartitions,
    ScaFES::Ntuple<int, DIM> const& nNodes,
    ScaFES::Ntuple<double, DIM> const& coordNodeFirst,
    ScaFES::Ntuple<double, DIM> const& coordNodeLast,
    ScaFES::Ntuple<bool, DIM> const& divideGridIntoDir)
: mNpartitions(nPartitions)
, mDiscreteDomain(ScaFES::Ntuple<int, DIM>(0),
                  nNodes - ScaFES::Ntuple<int, DIM>(1), coordNodeFirst,
                  coordNodeLast)
, mDecompDomain(ScaFES::Ntuple<int, DIM>(0), nPartitions, coordNodeFirst,
                coordNodeLast)
, mTypeDomainDecomp("RCB")
, mDivideGridIntoDirs(divideGridIntoDir)
, mPartitions()
, mNeighbourIds()
, mNeighbourDirs()
{
    initDependentMembers();
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline GridGlobal<DIM>::GridGlobal(ScaFES::GridGlobal<DIM> const& from)
: mNpartitions(from.nPartitions())
, mDiscreteDomain(from.discreteDomain())
, mDecompDomain(from.decompDomain())
, mTypeDomainDecomp(from.typeDomainDecomp())
, mDivideGridIntoDirs(from.divideGrid())
, mPartitions(from.partition())
, mNeighbourIds(from.neighbourId())
, mNeighbourDirs(from.neighbourDir())
{
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline ScaFES::GridGlobal<DIM>& GridGlobal<DIM>::
operator=(ScaFES::GridGlobal<DIM> rhs)
{
    swap(*this, rhs);
    return *this;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline GridGlobal<DIM>::~GridGlobal()
{
}

/*******************************************************************************
 * GETTER METHODS.
 ******************************************************************************/
template <std::size_t DIM>
inline ScaFES::Ntuple<int, DIM> const& GridGlobal<DIM>::nPartitions() const
{
    return this->mNpartitions;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline ScaFES::Grid<DIM> const& GridGlobal<DIM>::discreteDomain() const
{
    return this->mDiscreteDomain;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline ScaFES::Grid<DIM> const& GridGlobal<DIM>::decompDomain() const
{
    return this->mDecompDomain;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline std::string const& GridGlobal<DIM>::typeDomainDecomp() const
{
    return this->mTypeDomainDecomp;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline ScaFES::Ntuple<bool, DIM> const& GridGlobal<DIM>::divideGrid() const
{
    return this->mDivideGridIntoDirs;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline std::vector<GridSub<DIM>> const& GridGlobal<DIM>::partition() const
{
    return this->mPartitions;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline GridSub<DIM> const& GridGlobal<DIM>::partition(int const& ii) const
{
    return this->mPartitions.at(ii);
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline int const& GridGlobal<DIM>::neighbourId(int const& ii,
                                               int const& idx) const
{
    return this->mNeighbourIds.at(ii).at(idx);
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline std::vector<int> const& GridGlobal<DIM>::neighbourId(int const& ii) const
{
    return this->mNeighbourIds.at(ii);
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline std::vector<std::vector<int>> const& GridGlobal<DIM>::neighbourId() const
{
    return this->mNeighbourIds;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline int const& GridGlobal<DIM>::neighbourDir(int const& ii,
                                                int const& idx) const
{
    return this->mNeighbourDirs.at(ii).at(idx);
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline std::vector<int> const&
GridGlobal<DIM>::neighbourDir(int const& ii) const
{
    return this->mNeighbourDirs.at(ii);
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline std::vector<std::vector<int>> const&
GridGlobal<DIM>::neighbourDir() const
{
    return this->mNeighbourDirs;
}

/*******************************************************************************
 * WORK METHODS.
 ******************************************************************************/
template <std::size_t DIM>
inline int GridGlobal<DIM>::nPartitionsTotal() const
{
    int nPartitionsTotal = 1;
    for (std::size_t ii = 0; ii < this->nPartitions().dim(); ++ii)
    {
        nPartitionsTotal *= ((this->nPartitions())[ii]);
    }
    return nPartitionsTotal;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline bool GridGlobal<DIM>::operator==(GridGlobal<DIM> const& rhs) const
{
    // TODO: Some more comparisons are still missing. [KF, 23.10.2014]
    return (this->discreteDomain() == rhs.discreteDomain());
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline bool GridGlobal<DIM>::operator!=(GridGlobal<DIM> const& rhs) const
{
    return !(*this == rhs);
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline void GridGlobal<DIM>::partitionDomainUniformly()
{
    ScaFES::Ntuple<int, DIM> idxNodeFirstPart(
        this->discreteDomain().idxNodeFirst());
    ScaFES::Ntuple<int, DIM> idxNodeLastPart(
        this->discreteDomain().idxNodeFirst());
    this->mPartitions.resize(this->nPartitionsTotal());

    for (int kk = 0; kk < this->nPartitionsTotal(); ++kk)
    {
        ScaFES::Ntuple<int, DIM> vectIndices(this->nPartitions().dim());
        // Determine index tuple [i_0,i_1,...,i_{d-1}] from scalar index kk.
        int tmp = kk;
        for (std::size_t ii = this->nPartitions().dim() - 1; ii > 0; --ii)
        {
            vectIndices[ii] = (tmp % this->nPartitions()[ii]);
            tmp /= (this->nPartitions()[ii]);
        }
        vectIndices[0] = tmp;

        for (std::size_t pp = 0; pp < this->nPartitions().dim(); ++pp)
        {
            idxNodeFirstPart[pp] = this->discreteDomain().idxNodeFirst()[pp] +
                                   vectIndices[pp] *
                                       this->discreteDomain().nNodes()[pp] /
                                       this->nPartitions()[pp];
        }
        for (std::size_t pp = 0; pp < this->nPartitions().dim(); ++pp)
        {
            idxNodeLastPart[pp] = this->discreteDomain().idxNodeFirst()[pp] +
                                  (vectIndices[pp] + 1) *
                                      this->discreteDomain().nNodes()[pp] /
                                      this->nPartitions()[pp] -
                                  1;
        }
        this->mPartitions[kk] = GridSub<DIM>(idxNodeFirstPart, idxNodeLastPart,
                                             this->discreteDomain());
    }
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline void GridGlobal<DIM>::searchNeighboursUniform()
{
    unsigned long int gap = 1;

    for (int kk = 0; kk < this->nPartitionsTotal(); ++kk)
    {
        gap = 1;
        for (std::size_t pp = 0; pp < DIM; ++pp)
        {
            if (this->partition(kk).idxNodeFirstSub().elem(pp)
                != this->discreteDomain().idxNodeFirst(pp))
            {
                this->mNeighbourIds.at(kk).push_back(kk - gap);
                this->mNeighbourDirs.at(kk).push_back(2 * pp);

            }
            if (this->partition(kk).idxNodeLastSub().elem(pp)
                != this->discreteDomain().idxNodeLast(pp))
            {
                this->mNeighbourIds.at(kk).push_back(kk + gap);
                this->mNeighbourDirs.at(kk).push_back(2 * pp + 1);
            }
            // Compare numbering of partitions in partitionDomainUniformly().
            gap *= this->nPartitions()[DIM - 1 + pp];
        }
    }
}

/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline void GridGlobal<DIM>::partitionDomainUsingRCB(
    int const& nPartitionsTotal, ScaFES::Ntuple<int, DIM> const& idxNodeFirst,
    ScaFES::Ntuple<int, DIM> const& nNodes,
    ScaFES::Ntuple<bool, DIM> const& divideGrid)
{
    // Recursion tail: One partition.
    if (1 == nPartitionsTotal)
    {
        ScaFES::Ntuple<int, DIM> idxNodeLast =
            idxNodeFirst + nNodes - ScaFES::Ntuple<int, DIM>(1);
        assert(idxNodeFirst < (idxNodeFirst + nNodes));

        if (idxNodeLast <= idxNodeFirst)
        {
            throw std::runtime_error("idxNodeLast <= idxNodeFirst");
        }

        this->mPartitions.push_back(GridSub<DIM>(
            idxNodeFirst, idxNodeLast, this->discreteDomain().idxNodeFirst(),
            this->discreteDomain().idxNodeLast(),
            this->discreteDomain().coordNodeFirst(),
            this->discreteDomain().coordNodeLast()));
        return;
    }

    // Problem is decomposed into more than one partition.
    ScaFES::Ntuple<int, DIM> enDir(0);

    for (std::size_t kk = 0; kk < DIM; ++kk)
    {
        enDir.elem(kk) = static_cast<int>(divideGrid.elem(kk));
    }

    int maxnNodes = (enDir * nNodes).idxMaxElem();
    //     int maxnNodes = (divideGrid * nNodes).idxMaxElem();
    int max = nNodes.elem(maxnNodes);
    int p2 = nPartitionsTotal / 2;
    int c;
    std::vector<float> sum(nNodes.elem(maxnNodes), 0.0f);
    this->sumCpuLoads(sum, maxnNodes, idxNodeFirst, nNodes);

    float border = p2 * sum[max - 1] / nPartitionsTotal;

    for (c = 0; (c < max) && (sum[c] < border); ++c)
        ;
    {
        ++c;
    }
    ScaFES::Ntuple<int, DIM> s2 = idxNodeFirst;
    s2.elem(maxnNodes) += c;

    ScaFES::Ntuple<int, DIM> d1 = nNodes;
    ScaFES::Ntuple<int, DIM> d2 = nNodes;
    d1.elem(maxnNodes) += c - max;
    d2.elem(maxnNodes) -= c;

    // Recursive call.
    this->partitionDomainUsingRCB(p2, idxNodeFirst, d1, divideGrid);
    this->partitionDomainUsingRCB(nPartitionsTotal - p2, s2, d2, divideGrid);
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline void GridGlobal<DIM>::searchNeighbours(int const& direction)
{
    // Direction (0|1|...|DIM-1) for grid dimensions.
    ScaFES::Ntuple<int, DIM - 1> dirs;

    const int DIM_MINUS_ONE = DIM - 1;
    for (int kk = 0; kk < DIM_MINUS_ONE; ++kk)
    {
        dirs[kk] = (direction + 1) % DIM;
    }

    bool addGrid = false;
    ScaFES::Ntuple<int, DIM> tempSt;
    int idJJ;

    // Loop over all partitions.
    for (std::size_t ii = 0; ii < this->partition().size(); ++ii)
    {
        // Loop over all other partitions.
        for (std::size_t jj = 0; jj < this->partition().size(); ++jj)
        {
            idJJ = jj;

            // Ignore reference cuboid.
            if (jj != ii)
            {
                ScaFES::Ntuple<int, DIM> startII =
                    this->partition(ii).idxNodeFirstSub();
                ScaFES::Ntuple<int, DIM> endII =
                    this->partition(ii).idxNodeLastSub();
                ScaFES::Ntuple<int, DIM> startJJ =
                    this->partition(jj).idxNodeFirstSub();
                ScaFES::Ntuple<int, DIM> endJJ =
                    this->partition(jj).idxNodeLastSub();

                // Do the cuboids have a common face?
                // in der gleichen Ebene?
                if (startII.elem(direction) == (endJJ.elem(direction) + 1) ||
                    (endII.elem(direction) + 1) == startJJ.elem(direction))
                {

                    // Case DIM=1: There DOES exist a neighboured partition.
                    // Overlap the cuboids with this layer?
                    // Quick hack: Partitions are in x-y-layer.
                    // TODO: How does it look in general??
                    addGrid = true;

                    for (int ll = 0; ll < DIM_MINUS_ONE; ++ll)
                    {
                        if (std::min(endII.elem(dirs[ll]),
                                     endJJ.elem(dirs[ll])) <
                            std::max(startII.elem(dirs[ll]),
                                     startJJ.elem(dirs[ll])))
                        {
                            addGrid = false;
                            break;
                        }
                    }

                    if (addGrid)
                    {
                        tempSt = ScaFES::Ntuple<int, DIM>(0);
                        tempSt[direction] = std::max(startII.elem(direction),
                                                     startJJ.elem(direction));

                        for (int kk = 0; kk < DIM_MINUS_ONE; ++kk)
                        {
                            tempSt[dirs[kk]] = std::max(startII.elem(dirs[kk]),
                                                        startJJ.elem(dirs[kk]));
                        }

                        // TODO: Why?
                        int dn = 2 * direction;

                        if (!(tempSt < (endII + ScaFES::Ntuple<int, DIM>(1))))
                        {
                            ++dn;
                        }

                        this->mNeighbourIds.at(ii).push_back(idJJ);
                        this->mNeighbourDirs.at(ii).push_back(dn);
                    }
                }
            }
        }
    }
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline void GridGlobal<DIM>::searchNeighbours()
{
    for (std::size_t ii = 0; ii < DIM; ++ii)
    {
        this->searchNeighbours(ii);
    }
}

/*******************************************************************************
 * PRIVATE WORK METHODS.
 ******************************************************************************/
template <std::size_t DIM>
inline void
GridGlobal<DIM>::sumCpuLoads(std::vector<float>& sum, int const& direction,
                             ScaFES::Ntuple<int, DIM> const& idxNodeFirst,
                             ScaFES::Ntuple<int, DIM> const& nNodes) const
{
    assert(0 <= direction && direction <= static_cast<int>(DIM));

    int d1 = (direction + 1) % DIM;
    int d2 = (direction + 2) % DIM;

    ScaFES::Ntuple<int, DIM> idxNode;
    double tmpSum = 0.0;

    for (int ii = 0; ii < nNodes.elem(direction); ++ii)
    {
        idxNode[direction] = ii + idxNodeFirst.elem(direction);

        for (int jj = idxNodeFirst.elem(d1);
             jj < idxNodeFirst.elem(d1) + nNodes.elem(d1); ++jj)
        {
            idxNode[d1] = jj;

            for (int kk = idxNodeFirst.elem(d2);
                 kk < idxNodeFirst.elem(d2) + nNodes.elem(d2); ++kk)
            {
                idxNode[d2] = kk;
                // tmpSum += static_cast<double>(loadComp(index));
                tmpSum += 1.0;
            }
        }
        sum[ii] = static_cast<float>(tmpSum);
    }
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline void GridGlobal<DIM>::initDependentMembers()
{
    if (0 < this->nPartitionsTotal())
    {
        this->mNeighbourIds.resize(this->nPartitionsTotal());
        this->mNeighbourDirs.resize(this->nPartitionsTotal());
    }

    if (ScaFES::Ntuple<int, DIM>(0) < this->discreteDomain().nNodes())
    {
        if (0 < this->nPartitionsTotal())
        {
            if (0 == this->typeDomainDecomp().compare("UNI"))
            {
                this->partitionDomainUniformly();
                this->searchNeighboursUniform();
            }
            else
            {
                this->partitionDomainUsingRCB(
                    this->nPartitionsTotal(), ScaFES::Ntuple<int, DIM>(0),
                    this->discreteDomain().nNodes(), this->divideGrid());
                this->searchNeighbours();
            }
        }
    }
}
/*----------------------------------------------------------------------------*/
#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
template <std::size_t DIM>
template <class Archive>
inline void GridGlobal<DIM>::serialize(Archive& ar, const unsigned int version)
{
    if (1 <= version)
    {
        ar&(this->mDiscreteDomain);
        ar&(this->mDecompDomain);
        ar&(this->mTypeDomainDecomp);
        ar&(this->mNpartitions);
        ar&(this->mDivideGridIntoDirs);
        ar&(this->mPartitions);
        ar&(this->mNeighbourIds);
        ar&(this->mNeighbourDirs);
    }
}
#endif

/*******************************************************************************
 * FREE METHODS.
 ******************************************************************************/
template <std::size_t DIM>
inline void swap(ScaFES::GridGlobal<DIM>& first,
                 ScaFES::GridGlobal<DIM>& second)
{
    ScaFES::swap(first.mDiscreteDomain, second.mDiscreteDomain);
    ScaFES::swap(first.mDecompDomain, second.mDecompDomain);
    std::swap(first.mTypeDomainDecomp, second.mTypeDomainDecomp);
    ScaFES::swap(first.mNpartitions, second.mNpartitions);
    ScaFES::swap(first.mDivideGridIntoDirs, second.mDivideGridIntoDirs);
    std::swap(first.mPartitions, second.mPartitions);
    std::swap(first.mNeighbourIds, second.mNeighbourIds);
    std::swap(first.mNeighbourDirs, second.mNeighbourDirs);
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline std::ostream& operator<<(std::ostream& output,
                                ScaFES::GridGlobal<DIM> const& gg)
{
    for (int ii = 0; ii < gg.nPartitionsTotal(); ++ii)
    {
        output << std::fixed;
        output << ii << ":   ";
        for (std::size_t jj = 0; jj < gg.neighbourId(ii).size(); ++jj)
        {
             output << "(" << gg.neighbourId(ii, jj) << ";"
                    << gg.neighbourDir(ii, jj) << "), ";
        }
        output << std::endl;
    }

    return output << std::endl << std::endl;
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
    template <std::size_t DIM>
    struct version<ScaFES::GridGlobal<DIM>>
    {
        /** Sets the version number for serialization. */
        BOOST_STATIC_CONSTANT(unsigned long int, value = 2);
    };
} // namespace serialization
} // namespace boost
#endif

#endif
