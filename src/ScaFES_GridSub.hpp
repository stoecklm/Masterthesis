/* ScaFES
 * Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file ScaFES_GridSub.hpp
 *  @brief Contains the class template GridSub.
 */

#ifndef SCAFES_GRIDSUB_HPP_
#define SCAFES_GRIDSUB_HPP_

#include "ScaFES_Config.hpp"

#include <iterator>
#include <vector>

#ifdef SCAFES_HAVE_BOOST
#include <boost/version.hpp>
#endif

#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
// Forward declaration.
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
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>

#include <boost/config/suffix.hpp>
#endif


#include "ScaFES_Ntuple.hpp"
#include "ScaFES_Ntuple_FreeFunctions.hpp"
#include "ScaFES_Grid.hpp"

namespace ScaFES
{

/*******************************************************************************
 ******************************************************************************/
/** \class GridSub
 * @brief The class template \c GridSub represents a sub grid
 * of dimension \c DIM, where \c DIM is the template parameter of the class
 * template.
 *
 * The grid nodes of the sub grid are lexicographically numbered, too.
 * Therefore, only the first and the last node number of the sub grid together
 * with its corresponding coordinates have to be stored.
 *
 * Furthermore, the class contains an iterator class. The iterator can be
 * used like an iterator of the STL.
 *
 * It is assumed that the sub grid and the base grid have the same
 * grid size.
 */
template <std::size_t DIM = 3> class GridSub : public Grid<DIM>
{
public:
    /*----------------------------------------------------------------------
    | INTERNAL CLASSES.
    ----------------------------------------------------------------------*/
    /** Class for iterating through the grid.*/
    class Iterator : public std::iterator<std::random_access_iterator_tag,
                                          ScaFES::Ntuple<int, DIM>>
    {

    public:
        /*--------------------------------------------------------------
        | LIFE CYCLE METHODS.
        --------------------------------------------------------------*/
        /** Creates own constructor. */
        Iterator(const ScaFES::GridSub<DIM>* gd,
                 const ScaFES::Ntuple<int, DIM>& currIdxNode);

        Iterator();

        //                 /** Creates own copy constructor. */
        //                 Iterator(Iterator const&); // = default;
        //
        //                 /** Creates own assignment operator. */
        //                 Iterator& operator=(Iterator const&); // = default;

        /*--------------------------------------------------------------
        | GETTER METHODS.
        --------------------------------------------------------------*/
        /** Returns the current global node number. */
        const ScaFES::Ntuple<int, DIM>& idxNode() const;

        /** Returns the current local node number. */
        const unsigned long int& idxScalarNode() const;

        /*--------------------------------------------------------------
        | COMPARISON METHODS.
        --------------------------------------------------------------*/
        /** Checks if the local node number of this iterator is unequal
         * to the one of the given iterator. */
        bool operator!=(const Iterator& other) const;

        /** Checks if the local node number of this iterator is smaller
         * than the one of the given iterator. */
        bool operator<(const Iterator& other) const;

        /*--------------------------------------------------------------
        | WORK METHODS.
        --------------------------------------------------------------*/
        /** Returns the current global node number. */
        const ScaFES::Ntuple<int, DIM>& operator*();

        /** Updates the member variables of the iterator class:
         * Increases the current local and the current global
         * node number.
         * \remarks Prefix variant.
         */
        Iterator& operator++();

        /** Updates the member variables of the iterator class:
         * Increases the current local and the current global node
         * number.
         * \remarks Postfix variant.
         */
        Iterator operator++(int /*incr*/);

    private:
        /*--------------------------------------------------------------
        | FRIEND CLASSES.
        --------------------------------------------------------------*/
        friend class GridSub<DIM>;

        /*--------------------------------------------------------------
        | MEMBER VARIABLES.
        --------------------------------------------------------------*/
        /** Grid over which should be iterated. */
        const ScaFES::GridSub<DIM>* mGridIter;

        /** Current global node number. */
        ScaFES::Ntuple<int, DIM> mIdxNode;

        /** Current local node number. */
        unsigned long int mIdxScalarNode;

        /** Helper variable in order to improve the performance.
         * Jumps in numbering while iterating rhs one grid dimension
         * to the next one. These jumps occur due to the fact that
         * the base grid is normally not the sub grid. */
        ScaFES::Ntuple<int, DIM> mJumpInNumbering;

        /*--------------------------------------------------------------
        | INTERNAL MEMBER METHODS.
        --------------------------------------------------------------*/
        /** Updates the member variables of the iterator class:
         * Increases the current local and the current global
         * node number iterating "columnwise" through the grid.
         */
        Iterator& nextStyleC();

        /** Updates the member variables of the iterator class:
         * Increases the current local and the current global
         * node number iterating "rowwise" through the grid.
         */
        Iterator& nextStyleFortran();
    };

    /*----------------------------------------------------------------------
    | TYPE DEFINITIONS.
    ----------------------------------------------------------------------*/
    /** Type definition for iterator class analogously to STL class. */
    typedef Iterator iterator;

    /** Type definition for iterator class analogously to STL class. */
    typedef const Iterator const_iterator;

    /** Type definition for iterator class analogously to STL class. */
    typedef ScaFES::Ntuple<int, DIM>& reference_type;

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
    GridSub();

    /** Creates a new sub grid with given first global node number and
     * given last global node related to a given base grid.
     * \param idxNodeFirstSub First node number of the sub grid.
     * \param idxNodeLastSub Last node number of the sub grid.
     * \param idxNodeFirst First node number of the BASE GRID.
     * \param idxNodeLast Last node number of the BASE GRID.
     * \param coordNodeFirst Coordinates of first node number.
     * \param coordNodeLast Coordinates of last node number.
     */
    GridSub(const ScaFES::Ntuple<int, DIM>& idxNodeFirstSub,
            const ScaFES::Ntuple<int, DIM>& idxNodeLastSub,
            const ScaFES::Ntuple<int, DIM>& idxNodeFirst,
            const ScaFES::Ntuple<int, DIM>& idxNodeLast,
            ScaFES::Ntuple<double, DIM> coordNodeFirst,
            ScaFES::Ntuple<double, DIM> coordNodeLast);

    /** Creates own constructor.
     * The first and the last node number are explicitly set.
     * \param idxNodeFirstSub First node number of the sub grid.
     * \param idxNodeLastSub Last node number of the sub grid.
     * \param ggBase BASE GRID.
     */
    GridSub(const ScaFES::Ntuple<int, DIM>& idxNodeFirstSub,
            const ScaFES::Ntuple<int, DIM>& idxNodeLastSub,
            const Grid<DIM>& ggBase);

    /** Creates own copy constructor. */
    GridSub(const GridSub<DIM>& rhs);

    /** Creates own assignment operator using the copy-and-swap idiom.
     */
    GridSub& operator=(GridSub<DIM> rhs);

    /** Creates own destructor.  */
    ~GridSub();

    /*----------------------------------------------------------------------
    | GETTER METHODS.
    ----------------------------------------------------------------------*/
    /** Returns the first node number. */
    const ScaFES::Ntuple<int, DIM>& idxNodeFirstSub() const;

    /** Returns the component \c idx of the first node number. */
    const int& idxNodeFirstSub(const int& idx) const;

    /** Returns the last node number. */
    const ScaFES::Ntuple<int, DIM>& idxNodeLastSub() const;

    /** Returns the component \c idx of the last node number. */
    const int& idxNodeLastSub(const int& idx) const;

    /** Returns the first local node number in the grid. */
    const unsigned long int& idxScalarNodeFirst() const;

    /** Returns the last local node number in the grid. */
    const unsigned long int& idxScalarNodeLast() const;

    /** Returns the number of all nodes of the grid. */
    const int& nNodesTotal() const;

    /** Returns the number of nodes in each direction. */
    const ScaFES::Ntuple<int, DIM>& nNodesSub() const;

    /** Returns the number of nodes in one given dimension. */
    const int& nNodesSub(const int& idx) const;

    /*----------------------------------------------------------------------
    | COMPARISON METHODS.
    ----------------------------------------------------------------------*/
    /** Tests if two grids are equal. */
    bool operator==(const ScaFES::GridSub<DIM>& rhs) const;

    /** Tests if two grids are not equal. */
    bool operator!=(const ScaFES::GridSub<DIM>& rhs) const;

    /*----------------------------------------------------------------------
    | WORK METHODS.
    ----------------------------------------------------------------------*/
    /** Returns the iterator corresponding to the underlying grid. */
    iterator begin() const;

    /** Returns the default iterator. */
    iterator end() const;

    /** Returns true if the first global node number is equal or smaller the
     * the last global node number, i.e. if the grid consists of at
     * least one node. */
    bool valid() const;

    /** Returns true if a given global node number is inside the grid.
     */
    bool inside(const ScaFES::Ntuple<int, DIM>& idxNode) const;

    /** Extends the grid by a given border width into a given direction.
     * \param dir Direction to expand.
     * \param nLayers Number of layers to expand.
     */
    GridSub<DIM> extend(const int& dir, const int& nLayers);

    /** Returns the union of this grid and another given grid.
     * Assumption: Both grids must have the same base grid.
     * \return Grid gg(min(mIdxNodeFirst, rhs.idxNodeFirstSub),
     *                 max(mIdxNodeLast, rhs.idxNodeLastSub))
     */
    GridSub<DIM> unionWith(const GridSub<DIM>& rhs) const;

    /** Returns the intersection of this grid and another given grid.
     * Assumption: Both grids must have the same base grid.
     * \return Grid gg(max(mIdxNodeFirst, rhs.idxNodeFirstSub),
     *                 min(mIdxNodeLast, rhs.idxNodeLastSub))
     */
    GridSub<DIM> intersectWith(const GridSub<DIM>& rhs) const;

    /** Serializes the class. */
    template <class Archive>
    void serialize(Archive& ar, const unsigned int version);

    /*----------------------------------------------------------------------
    | FREE METHODS.
    ----------------------------------------------------------------------*/
    /** Swaps the members of two grids. */
    template <std::size_t SS>
    friend void swap(GridSub<SS>& first, GridSub<SS>& second);

    /** Overloads the output operator '<<'. */
    template <std::size_t SS>
    friend std::ostream& operator<<(std::ostream& output,
                                    const GridSub<SS>& rhs);

private:
    /*----------------------------------------------------------------------
    | MEMBER VARIABLES.
    ----------------------------------------------------------------------*/
    /** First global node number of the grid. */
    ScaFES::Ntuple<int, DIM> mIdxNodeFirstSub;

    /** Last global node number of the grid.
     * \remarks The number of nodes in each direction is given by the
     difference of the last and first node number + 1. */
    ScaFES::Ntuple<int, DIM> mIdxNodeLastSub;

    /** Number of nodes of the grid in all dimensions. */
    ScaFES::Ntuple<int, DIM> mNnodesSub;

    /** Total number of grid nodes. */
    int mNnodesTotal;

    /** First MPI local node number of the grid. */
    unsigned long int mIdxScalarNodeFirst;

    /** Last MPI local node number of the grid. */
    unsigned long int mIdxScalarNodeLast;

    /*----------------------------------------------------------------------
    | INTERNAL METHODS.
    ----------------------------------------------------------------------*/
    /** Computes the corresponding local node number of a given global
     * node number of the grid using the C style.
     * \param idxNode Given global node number
     */
    unsigned long int
    idxNodeTuple2ScalarStyleC(const ScaFES::Ntuple<int, DIM>& idxNode) const;

    /** Computes the corresponding local node numbervof a given global
     * node number of the grid using the FORTRAN style.
     * \param idxNode Given global node number
     */
    unsigned long int idxNodeTuple2ScalarStyleFortran(
        const ScaFES::Ntuple<int, DIM>& idxNode) const;
}; // End of class. //

/*******************************************************************************
 * FREE METHODS.
 ******************************************************************************/
/** Swaps two objects of the class \c GridSub. */
template <std::size_t DIM>
void swap(ScaFES::GridSub<DIM>& first, ScaFES::GridSub<DIM>& second);
/*----------------------------------------------------------------------------*/
/** Preprares writing an object of the class \c GridSub to output. */
template <std::size_t DIM>
std::ostream& operator<<(std::ostream& output, const ScaFES::GridSub<DIM>& gg);

/*******************************************************************************
 * LIFE CYCLE METHODS.
 ******************************************************************************/
template <std::size_t DIM>
inline GridSub<DIM>::GridSub()
: Grid<DIM>(), mIdxNodeFirstSub(ScaFES::Ntuple<int, DIM>(0)),
  mIdxNodeLastSub(ScaFES::Ntuple<int, DIM>(-1)),
  mNnodesSub(ScaFES::Ntuple<int, DIM>(0)), mNnodesTotal(0),
  mIdxScalarNodeFirst(0), mIdxScalarNodeLast(0)
{
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline GridSub<DIM>::GridSub(const ScaFES::Ntuple<int, DIM>& idxNodeFirstSub,
                             const ScaFES::Ntuple<int, DIM>& idxNodeLastSub,
                             const ScaFES::Ntuple<int, DIM>& idxNodeFirst,
                             const ScaFES::Ntuple<int, DIM>& idxNodeLast,
                             ScaFES::Ntuple<double, DIM> coordNodeFirst,
                             ScaFES::Ntuple<double, DIM> coordNodeLast)
: ScaFES::Grid<DIM>(idxNodeFirst, idxNodeLast, coordNodeFirst, coordNodeLast),
  mIdxNodeFirstSub(idxNodeFirstSub), mIdxNodeLastSub(idxNodeLastSub),
  mNnodesSub(mIdxNodeLastSub - mIdxNodeFirstSub + ScaFES::Ntuple<int, DIM>(1)),
  mNnodesTotal(mNnodesSub.size()),
  mIdxScalarNodeFirst(this->idxNodeTuple2Scalar(mIdxNodeFirstSub)),
  mIdxScalarNodeLast(this->idxNodeTuple2Scalar(mIdxNodeLastSub))
{
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline GridSub<DIM>::GridSub(const ScaFES::Ntuple<int, DIM>& idxNodeFirstSub,
                             const ScaFES::Ntuple<int, DIM>& idxNodeLastSub,
                             const ScaFES::Grid<DIM>& ggBase)
: ScaFES::Grid<DIM>(ggBase), mIdxNodeFirstSub(idxNodeFirstSub),
  mIdxNodeLastSub(idxNodeLastSub),
  mNnodesSub(mIdxNodeLastSub - mIdxNodeFirstSub + ScaFES::Ntuple<int, DIM>(1)),
  mNnodesTotal(mNnodesSub.size()),
  mIdxScalarNodeFirst(this->idxNodeTuple2Scalar(mIdxNodeFirstSub)),
  mIdxScalarNodeLast(this->idxNodeTuple2Scalar(mIdxNodeLastSub))
{
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline GridSub<DIM>::GridSub(const GridSub<DIM>& rhs)
: Grid<DIM>(rhs), mIdxNodeFirstSub(rhs.idxNodeFirstSub()),
  mIdxNodeLastSub(rhs.idxNodeLastSub()), mNnodesSub(rhs.nNodesSub()),
  mNnodesTotal(rhs.nNodesTotal()),
  mIdxScalarNodeFirst(rhs.idxScalarNodeFirst()),
  mIdxScalarNodeLast(rhs.idxScalarNodeLast())
{
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline GridSub<DIM>& GridSub<DIM>::operator=(GridSub<DIM> rhs)
{
    swap(*this, rhs);
    return *this;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM> inline GridSub<DIM>::~GridSub()
{
}

/*******************************************************************************
 * GETTER METHODS.
 ******************************************************************************/
template <std::size_t DIM>
inline const ScaFES::Ntuple<int, DIM>& GridSub<DIM>::idxNodeFirstSub() const
{
    return this->mIdxNodeFirstSub;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline const int& GridSub<DIM>::idxNodeFirstSub(const int& idx) const
{
    return this->mIdxNodeFirstSub.elem(idx);
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline const ScaFES::Ntuple<int, DIM>& GridSub<DIM>::idxNodeLastSub() const
{
    return this->mIdxNodeLastSub;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline const int& GridSub<DIM>::idxNodeLastSub(const int& idx) const
{
    return this->mIdxNodeLastSub.elem(idx);
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline const unsigned long int& GridSub<DIM>::idxScalarNodeFirst() const
{
    return this->mIdxScalarNodeFirst;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline const unsigned long int& GridSub<DIM>::idxScalarNodeLast() const
{
    return this->mIdxScalarNodeLast;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline const ScaFES::Ntuple<int, DIM>& GridSub<DIM>::nNodesSub() const
{
    return this->mNnodesSub;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline const int& GridSub<DIM>::nNodesSub(const int& idx) const
{
    return this->mNnodesSub.elem(idx);
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM> inline const int& GridSub<DIM>::nNodesTotal() const
{
    return this->mNnodesTotal;
}

/*******************************************************************************
 * COMPARISON METHODS.
 ******************************************************************************/
template <std::size_t DIM>
inline bool GridSub<DIM>::operator==(const GridSub<DIM>& rhs) const
{
    return (this->idxNodeFirstSub() == rhs.idxNodeFirstSub() &&
            this->idxNodeLastSub() == rhs.idxNodeLastSub());
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline bool GridSub<DIM>::operator!=(const GridSub<DIM>& rhs) const
{
    return !(*this == rhs);
}

/*******************************************************************************
 * WORK METHODS.
 ******************************************************************************/
template <std::size_t DIM>
inline typename GridSub<DIM>::iterator GridSub<DIM>::begin() const
{
    return iterator(this, this->idxNodeFirstSub());
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline typename GridSub<DIM>::iterator GridSub<DIM>::end() const
{
    return iterator(this, this->idxNodeLastSub());
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM> inline bool GridSub<DIM>::valid() const
{
    return (this->idxNodeFirstSub() <= this->idxNodeLastSub());
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline bool GridSub<DIM>::inside(const ScaFES::Ntuple<int, DIM>& idxNode) const
{
    return ((this->idxNodeFirstSub() <= idxNode) &&
            (idxNode <= this->idxNodeLastSub()));
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline GridSub<DIM> GridSub<DIM>::extend(const int& dir, const int& nLayers)
{
    ScaFES::Ntuple<int, DIM> idxNodeFirstSub(this->mIdxNodeFirstSub);
    ScaFES::Ntuple<int, DIM> idxNodeLastSub(this->mIdxNodeLastSub);
    int idx = dir / 2;

    if (0 == (dir % 2))
    {
        idxNodeLastSub[idx] += nLayers;
    }
    else
    {
        idxNodeFirstSub[idx] -= nLayers;
    }

    return ScaFES::GridSub<DIM>(idxNodeFirstSub, idxNodeLastSub,
                                static_cast<Grid<DIM>>(*this));
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline GridSub<DIM> GridSub<DIM>::unionWith(const GridSub<DIM>& rhs) const
{
    return ScaFES::GridSub<DIM>(
        min(this->idxNodeFirstSub(), rhs.idxNodeFirstSub()),
        max(this->idxNodeLastSub(), rhs.idxNodeLastSub()), this->idxNodeFirst(),
        this->idxNodeLast(), this->coordNodeFirst(), this->coordNodeLast());
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline GridSub<DIM> GridSub<DIM>::intersectWith(const GridSub<DIM>& rhs) const
{
    return GridSub<DIM>(max(this->idxNodeFirstSub(), rhs.idxNodeFirstSub()),
                        min(this->idxNodeLastSub(), rhs.idxNodeLastSub()),
                        this->idxNodeFirst(), this->idxNodeLast(),
                        this->coordNodeFirst(), this->coordNodeLast());
}
/*----------------------------------------------------------------------------*/
#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
template <std::size_t DIM>
template <class Archive>
inline void GridSub<DIM>::serialize(Archive& ar, const unsigned int version)
{
    if (1 <= version)
    {
        ar& boost::serialization::base_object<ScaFES::Grid<DIM>>(*this);
        ar&(this->mIdxNodeFirstSub);
        ar&(this->mIdxNodeLastSub);
        ar&(this->mNnodesSub);
        ar&(this->mNnodesTotal);
        ar&(this->mIdxScalarNodeFirst);
        ar&(this->mIdxScalarNodeLast);
    }
}
#endif
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline unsigned long int GridSub<DIM>::idxNodeTuple2ScalarStyleC(
    const ScaFES::Ntuple<int, DIM>& idxNode) const
{
    unsigned long int idxScalarNode =
        idxNode.elem(0) - this->idxNodeFirstSub(0);

    for (std::size_t ii = 1; ii < DIM; ++ii)
    {
        idxScalarNode = idxScalarNode * this->nNodesSub(ii) + idxNode.elem(ii) -
                        this->idxNodeFirst(ii);
    }

    return idxScalarNode;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline unsigned long int GridSub<DIM>::idxNodeTuple2ScalarStyleFortran(
    const ScaFES::Ntuple<int, DIM>& idxNode) const
{
    unsigned long int idxScalarNode =
        idxNode.elem(DIM - 1) - this->idxNodeFirstSub(DIM - 1);
    const int dimLocal = DIM;

    for (int ii = dimLocal - 2; ii >= 0; --ii)
    {
        idxScalarNode = idxScalarNode * this->nNodesSub(ii) + idxNode.elem(ii) -
                        this->idxNodeFirst(ii);
    }

    return idxScalarNode;
}

/*******************************************************************************
 * FREE METHODS.
 ******************************************************************************/
template <std::size_t DIM>
inline void swap(GridSub<DIM>& first, GridSub<DIM>& second)
{
    swap(static_cast<Grid<DIM>&>(first), static_cast<Grid<DIM>&>(second));
    ScaFES::swap(first.mIdxNodeFirstSub, second.mIdxNodeFirstSub);
    ScaFES::swap(first.mIdxNodeLastSub, second.mIdxNodeLastSub);
    ScaFES::swap(first.mNnodesSub, second.mNnodesSub);
    int tmp = first.mNnodesTotal;
    first.mNnodesTotal = second.mNnodesTotal;
    second.mNnodesTotal = tmp;
    std::swap(first.mIdxScalarNodeFirst, second.mIdxScalarNodeFirst);
    std::swap(first.mIdxScalarNodeLast, second.mIdxScalarNodeLast);
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline std::ostream& operator<<(std::ostream& output,
                                const ScaFES::GridSub<DIM>& gg)
{
    output << "Sub Grid   <-" << gg.idxNodeFirstSub() << "-"
           << gg.idxNodeLastSub() << "->\n"
           << "Base Grid: <-" << gg.idxNodeFirst() << "-" << gg.idxNodeLast()
           << "->";
    return output;
}

/*******************************************************************************
 * INTERNAL CLASS ITERATOR: LIFE CYCLE METHODS.
 ******************************************************************************/
template <std::size_t DIM>
inline GridSub<DIM>::Iterator::Iterator()
: mGridIter(NULL), mIdxNode(ScaFES::Ntuple<int, DIM>(0)), mIdxScalarNode(0L)
{
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline GridSub<DIM>::Iterator::Iterator(
    const ScaFES::GridSub<DIM>* gd, const ScaFES::Ntuple<int, DIM>& currIdxNode)
: mGridIter(gd), mIdxNode(currIdxNode),
  mIdxScalarNode(gd->idxNodeTuple2Scalar(currIdxNode)),
  mJumpInNumbering()
{
    // TODO: Check if \c currIdxNode is inside underlying grid.
    if (gd->isNumberedStyleC())
    {
        this->mJumpInNumbering[DIM - 1] = 1;
        this->mJumpInNumbering[DIM - 2] =
            (this->mGridIter->idxNodeLast(DIM - 1) -
             this->mGridIter->idxNodeLastSub(DIM - 1)) -
            (this->mGridIter->idxNodeFirst(DIM - 1) -
             this->mGridIter->idxNodeFirstSub(DIM - 1));
        const int dimLocal = DIM;
        for (int ii = dimLocal - 2; ii >= 1; --ii)
        {
            this->mJumpInNumbering[ii - 1] =
                this->mJumpInNumbering.elem(ii) +
                this->mGridIter->nNodes(ii + 1) *
                    ((this->mGridIter->idxNodeLast(ii) -
                      this->mGridIter->idxNodeLastSub(ii)) -
                     (this->mGridIter->idxNodeFirst(ii) -
                      this->mGridIter->idxNodeFirstSub(ii)));
        }
    }
    else
    {
        this->mJumpInNumbering[0] = 1;
        this->mJumpInNumbering[1] = (this->mGridIter->idxNodeLast(0) -
                                     this->mGridIter->idxNodeLastSub(0)) -
                                    (this->mGridIter->idxNodeFirst(0) -
                                     this->mGridIter->idxNodeFirstSub(0));

        const int DIM_MINUS_ONE = DIM - 1;
        for (int ii = 1; ii < DIM_MINUS_ONE; ++ii)
        {
            this->mJumpInNumbering[ii + 1] =
                this->mJumpInNumbering.elem(ii) +
                this->mGridIter->nNodes(ii - 1) *
                    ((this->mGridIter->idxNodeLast(ii) -
                      this->mGridIter->idxNodeLastSub(ii)) -
                     (this->mGridIter->idxNodeFirst(ii) -
                      this->mGridIter->idxNodeFirstSub(ii)));
        }
    }
}

/*******************************************************************************
 * INTERNAL CLASS ITERATOR: GETTER METHODS.
 ******************************************************************************/
template <std::size_t DIM>
inline const ScaFES::Ntuple<int, DIM>& GridSub<DIM>::Iterator::idxNode() const
{
    return this->mIdxNode;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline const unsigned long int& GridSub<DIM>::Iterator::idxScalarNode() const
{
    return this->mIdxScalarNode;
}

/*******************************************************************************
 * INTERNAL CLASS ITERATOR: COMPARISON METHODS.
 ******************************************************************************/
template <std::size_t DIM>
inline bool GridSub<DIM>::Iterator::operator!=(const Iterator& other) const
{
    return (*this < other) || (other < *this);
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline bool GridSub<DIM>::Iterator::operator<(const Iterator& it) const
{
    return (this->idxScalarNode() <= it.idxScalarNode());
}

/*******************************************************************************
 * INTERNAL CLASS ITERATOR: WORK METHODS.
 ******************************************************************************/
template <std::size_t DIM>
inline const ScaFES::Ntuple<int, DIM>& GridSub<DIM>::Iterator::operator*()
{
    return this->idxNode();
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline typename GridSub<DIM>::Iterator& GridSub<DIM>::Iterator::operator++()
{
    if (this->mGridIter->isNumberedStyleC())
    {
        return this->nextStyleC();
    }

    return this->nextStyleFortran();
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline typename GridSub<DIM>::Iterator GridSub<DIM>::Iterator::operator++(int /*incr*/)
{
    Iterator tmp = *this;
    ++(*this);
    return tmp;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline typename GridSub<DIM>::Iterator& GridSub<DIM>::Iterator::nextStyleC()
{
    ++(this->mIdxScalarNode);
    // Case: Next node in first dimension?
    ++(this->mIdxNode[DIM - 1]);

    if (this->mIdxNode.elem(DIM - 1) > this->mGridIter->idxNodeLastSub(DIM - 1))
    {
        int dimLocal = DIM;

        // Case: Next node in (ii-1)-th dimension?
        for (int ii = dimLocal - 2; ii >= 0; --ii)
        {
            this->mIdxNode[ii + 1] = this->mGridIter->idxNodeFirstSub(ii + 1);
            ++(this->mIdxNode[ii]);

            if (this->mIdxNode.elem(ii) <= this->mGridIter->idxNodeLastSub(ii))
            {
                this->mIdxScalarNode += this->mJumpInNumbering.elem(ii);
                break;
            }
        }
    }

    return *this;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline typename GridSub<DIM>::Iterator&
GridSub<DIM>::Iterator::nextStyleFortran()
{
    ++(this->mIdxScalarNode);
    // Case: Next node in first dimension?
    ++(this->mIdxNode[0]);

    if (this->mIdxNode.elem(0) > this->mGridIter->idxNodeLastSub(0))
    {
        // Case: Next node in (ii+1)-th dimension?
        for (std::size_t ii = 1; ii < DIM; ++ii)
        {
            this->mIdxNode[ii - 1] = this->mGridIter->idxNodeFirstSub(ii - 1);
            ++(this->mIdxNode[ii]);

            if (this->mIdxNode.elem(ii) <= this->mGridIter->idxNodeLastSub(ii))
            {
                this->mIdxScalarNode += this->mJumpInNumbering.elem(ii);
                break;
            }
        }
    }

    return *this;
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
    template <std::size_t DIM> struct version<ScaFES::GridSub<DIM>>
    {
        /** Sets the version number for serialization. */
        BOOST_STATIC_CONSTANT(unsigned long int, value = 2);
    };
} // namespace serialization
} // namespace boost
#endif

#endif
