/* ScaFES
 * Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file ScaFES_Grid.hpp
 *  @brief Contains the class template Grid.
 */

#ifndef SCAFES_GRID_HPP_
#define SCAFES_GRID_HPP_

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

#include <boost/config/suffix.hpp>
#endif

#include "ScaFES_Ntuple.hpp"
#include "ScaFES_Ntuple_FreeFunctions.hpp"

namespace ScaFES
{

/*******************************************************************************
 ******************************************************************************/
/**
* \class Grid
* @brief The class template \c Grid represents a structured grid
* of dimension \c DIM, where \c DIM is the template parameter of the class
* template. A grid is a \c DIM - dimensional hypercuboid.
*
* The grid nodes are numbered in two different ways:
* <ul>
* <li> globally, i.e. all grid nodes of the WHOLE MPI domain are numbered
*      lexicographically, </li>
* <li> locally, i.e. all grid nodes of ONE MPI partition are numbered
*      lexicographically. </li>
* </ul>
*
* As the nodes of the grid are lexicographically numbered, the first and
* the last node number of the grid together with its corresponding coordinates
* have to be stored, only.
*
* Due to the lexicographical node numbering, the coordinates of all nodes
* can be easily computed. For this purpose, the method \c coordinates() was
* provided.
*
* Furthermore, the class contains an iterator class. As the iterator is
* designed like an STL iterator, it can be applied completely analogously.

* The most important components of the inner class \c Iterator are the member
* methods \c idxNode() and \c idxScalarNode(). The first one returns
* the current node number in tuple notation, the second on returns
* the current node number as a scalar number.
*
* Numbering of grid nodes for DIM = 3 of a grid with n^3 nodes:
* Gap in y-direction: n     (1, n+1, 2*n+1, 3*n+1,...)
* Gap in z-direction: n^2   (n^2, 2*n^2, 3*n^2,....)
*        *-----------------* n^3
*       /                 /|
*      /                 / |
*     /                 /  |
*    /     (n-1)*n^2+n /   |
*    *----------------*    |
*    |                |    * n^2
*    |                |   /
*    |                |  /
*    |                | /
*    |                |/
*  1 *----------------* n
*
* \n Numbering of grid faces for DIM = 3:
*                            ^          ^
*                            |         / Back = 3
*                            | Top = 4/
*                                    /
*                     *-----------------*
*                    /                 /|
*                   /                 / |
*                  /                 /  |
*                 /                 /   |
*                *-----------------*    |
*                |                 |    *
*       <------  |                 |   /   ------->
*        Left=1  |                 |  /    Right = 0
*                |                 | /
*                |                 |/
*                *-----------------*
* Right:  direction = 0
* Left:   direction = 1      | Bottom = 5
* Top:    direction = 2      |
* Bottom: direction = 3      |
* Front:  direction = 4      v
* Back:   direction = 5
*/
template <std::size_t DIM = 3> class Grid
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
        /** Creates own constructor.
         * \param gridIter Given grid over which will be iterated.
         * \param idxNode Iteration will start from this node.
         */
        Iterator(const Grid<DIM>* gridIter,
                 const ScaFES::Ntuple<int, DIM>& idxNode);

        /*--------------------------------------------------------------
        | GETTER METHODS.
        --------------------------------------------------------------*/
        /** Returns the index of the current node. */
        const ScaFES::Ntuple<int, DIM>& idxNode() const;

        /** Returns the scalar index of the current node. */
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
        /** Returns the index of the current node. */
        const ScaFES::Ntuple<int, DIM>& operator*();

        /** Increases the current local and the current global
         * node number.
         * \remarks Prefix variant.
         */
        Iterator& operator++();

        /** Increases the current local and the current global
         * node number.
         * \remarks Postfix variant.
         */
        Iterator operator++(int /*incr*/);

    private:
        /*--------------------------------------------------------------
        | MEMBER VARIABLES.
        --------------------------------------------------------------*/
        /** Grid over which should be iterated. */
        const Grid<DIM>* mGridIter;

        /** Index of current node. */
        ScaFES::Ntuple<int, DIM> mIdxNode;

        /** Scalar index of current node. */
        unsigned long int mIdxScalarNode;

        /*--------------------------------------------------------------
        | INTERNAL WORK METHODS.
        --------------------------------------------------------------*/
        /** Increases the current local and the current global node
         * number iterating "columnwise" through the grid.
         */
        Iterator& nextStyleC();

        /** Increases the current local and the current global node
         * number iterating "rowwise" through the grid.
         */
        Iterator& nextStyleFortran();

        /*--------------------------------------------------------------
        | FRIEND CLASSES.
        --------------------------------------------------------------*/
        friend class Grid<DIM>;
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
    Grid();

    /** Creates a new grid with given first global node number and
     * given last global node together with its corresponding coordinates.
     *  TODO: Reason for pass-by-value of coordinates?! [KF], 2015-11-17
     * \param idxNodeFirst Index of first node.
     * \param idxNodeLast Index of last node.
     * \param coordNodeFirst Coordinates of first node.
     * \param coordNodeLast Coordinates of last node.
     */
    Grid(const ScaFES::Ntuple<int, DIM>& idxNodeFirst,
         const ScaFES::Ntuple<int, DIM>& idxNodeLast,
         ScaFES::Ntuple<double, DIM> coordNodeFirst,
         ScaFES::Ntuple<double, DIM> coordNodeLast);

    /** Creates own copy constructor. */
    Grid(const Grid<DIM>& rhs);

    /** Creates own assignment operator using the copy-and-swap idiom. */
    Grid& operator=(Grid<DIM> rhs);

    /** Creates own destructor. */
    ~Grid();

    /*----------------------------------------------------------------------
    | GETTER METHODS.
    ----------------------------------------------------------------------*/
    /** Returns the dimension of the grid. */
    std::size_t dimGrid() const;

    /** Returns the index of the first node. */
    const ScaFES::Ntuple<int, DIM>& idxNodeFirst() const;

    /** Returns the component \c idx of the index of the first node. */
    const int& idxNodeFirst(const int& idx) const;

    /** Returns the index of the last node. */
    const ScaFES::Ntuple<int, DIM>& idxNodeLast() const;

    /** Returns the component \c idx of the index of the last node. */
    const int& idxNodeLast(const int& idx) const;

    /** Returns the number of nodes of the grid in each direction. */
    const ScaFES::Ntuple<int, DIM>& nNodes() const;

    /** Returns the number of nodes of the grid in one given direction. */
    const int& nNodes(const int& idx) const;

    /** Returns the total number of nodes of the grid. */
    int nNodesTotal() const;

    /** Returns the coordinates of the first node. */
    const ScaFES::Ntuple<double, DIM>& coordNodeFirst() const;

    /** Returns the coordinates of the last node. */
    const ScaFES::Ntuple<double, DIM>& coordNodeLast() const;

    /** Returns the grid size in all directions. */
    const ScaFES::Ntuple<double, DIM>& gridsize() const;

    /** Returns the grid size in one given direction. */
    const double& gridsize(const int& idx) const;

    /** Returns if the grid is numbered in C style notation. */
    const bool& isNumberedStyleC() const;

    /*----------------------------------------------------------------------
    | COMPARISON METHODS.
    ----------------------------------------------------------------------*/
    /** Tests if two grids are equal. */
    bool operator==(const Grid<DIM>& rhs) const;

    /** Tests if two grids are not equal. */
    bool operator!=(const Grid<DIM>& rhs) const;

    /*----------------------------------------------------------------------
    | WORK METHODS.
    ----------------------------------------------------------------------*/
    /** Returns an iterator which will start with the first node of the
     * underlying grid. */
    iterator begin() const;

    /** Returns an iterator which will start with the last node of the
     * underlying grid. */
    iterator end() const;

    /** Computes the corresponding scalar node number of a given node number
     * of the grid.
     * \param idxNode ScaFES::Ntuple<int,DIM> of a given node.
     */
    unsigned long int
    idxNodeTuple2Scalar(const ScaFES::Ntuple<int, DIM>& idxNode) const;

    /** Computes the corresponding scalar node number of a given node number
     * of the grid.
     * \param idxNodeScalar Scalar index of a given node.
     * \param idxNode ScaFES::Ntuple<int,DIM> of a given node.
     */
    void idxNodeTuple2Scalar(unsigned long int& idxNodeScalar,
                             const ScaFES::Ntuple<int, DIM>& idxNode) const;

    /** Computes the corresponding node number from a given scalar
     * node number.
     * \param idxScalarNode Scalar index of a given node number.
     */
    ScaFES::Ntuple<int, DIM>
    idxNodeScalar2Tuple(const unsigned long int& idxScalarNode) const;

    /** Computes the corresponding node number from a given scalar
     * node number.
     * \param idxNode ScaFES::Ntuple<int,DIM> tuple of a given node.
     * \param idxScalarNode Scalar index of a given node number.
     */
    void idxNodeScalar2Tuple(ScaFES::Ntuple<int, DIM>& idxNode,
                             const unsigned long int& idxScalarNode) const;

    /** Computes the coordinates of a given node number of the grid.
     * \param idxNode Given node number
     */
    ScaFES::Ntuple<double, DIM>
    coordinates(const ScaFES::Ntuple<int, DIM>& idxNode) const;

    /** Returns true if the first node number is equal or smaller the
     * the last node number, i.e. if the grid consists of at
     * least one node. */
    bool valid() const;

    /** Returns true if a given node number is inside the grid.
     */
    bool inside(const ScaFES::Ntuple<int, DIM>& idxNode) const;

    /** Extends the underlying grid by a given border width into a
     * given direction and returns the new modified grid.
     * The coordinates of the extended grid will remain the same as before,
     * whereas the number of nodes within the same physical domain will
     * increases.
     * \param dir Direction to expand.
     * \param nLayers Number of layers to expand into the given dir.
     * \remarks The ordering of the direction is as given in the
     * description of this class.
     */
    Grid<DIM> extend(const int& dir, const int& nLayers);

    /** Returns the union of this grid and another given grid.
     * It will be assumed that both grids have the same grid sizes!
     * \return Grid gg(min(mIdxNodeFirst, rhs.idxNodeFirst),
     *                 max(mIdxNodeLast, rhs.idxNodeLast))
     */
    Grid<DIM> unionWith(const Grid<DIM>& rhs) const;

    /** Returns the intersection of this grid and another given grid.
     * It will be assumed that both grids have the same grid sizes!
     * \return Grid gg(max(mIdxNodeFirst, rhs.idxNodeFirst),
     *                 min(mIdxNodeLast, rhs.idxNodeLast))
     */
    Grid<DIM> intersectWith(const Grid<DIM>& rhs) const;

#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
    /** Serializes the class. */
    template <class Archive>
    void serialize(Archive& ar, const unsigned int version);
#endif

    /** Returns all direct neighboured LOCAL node numbers
     * of a given node number. */
    Ntuple<unsigned long int, 2 * DIM>
    connect(const ScaFES::Ntuple<int, DIM>& idxNode);

    /** Determines the neighbour global node into a given direction
     * of a global node number. */
    ScaFES::Ntuple<int, DIM> neigh(const ScaFES::Ntuple<int, DIM>& idxNode,
                                   const int& dir) const;

    /*----------------------------------------------------------------------
    | FREE METHODS.
    ----------------------------------------------------------------------*/
    /** Swaps the members of two grids. */
    template <std::size_t SS>
    friend void swap(ScaFES::Grid<SS>& first, ScaFES::Grid<SS>& second);

    /** Overloads the output operator. */
    template <std::size_t SS>
    friend std::ostream& operator<<(std::ostream& output,
                                    const ScaFES::Grid<SS>& rhs);

protected:
    /*----------------------------------------------------------------------
    | MEMBER VARIABLES.
    ----------------------------------------------------------------------*/
    /** First node number of the grid. */
    ScaFES::Ntuple<int, DIM> mIdxNodeFirst;

    /** Last node number of the grid. */
    ScaFES::Ntuple<int, DIM> mIdxNodeLast;

    /** Number of nodes of the grid in all directions. */
    ScaFES::Ntuple<int, DIM> mNnodes;

    /** Coordinate of the first node number. */
    ScaFES::Ntuple<double, DIM> mCoordNodeFirst;

    /** Coordinate of the last node number. */
    ScaFES::Ntuple<double, DIM> mCoordNodeLast;

    /** Grid size in all directions. */
    ScaFES::Ntuple<double, DIM> mGridsize;

    /** Is grid numbered in C style notation (alternative: FORTRAN style).*/
    bool mIsNumberedStyleC;

private:
    /*----------------------------------------------------------------------
    | INTERNAL METHODS.
    ----------------------------------------------------------------------*/
    /** Computes the corresponding scalar node number of a given
     * node number of the grid using the C style.
     * \param idxNode Given node number
     */
    unsigned long int
    idxNodeTuple2ScalarStyleC(const ScaFES::Ntuple<int, DIM>& idxNode) const;

    /** Computes the corresponding scalar node number of a given
     * node number of the grid using the C style.
     * \param idxScalarNode Corresponding scalar node number (return value).
     * \param idxNode Given node number.
     */
    void
    idxNodeTuple2ScalarStyleC(unsigned long int& idxScalarNode,
                              const ScaFES::Ntuple<int, DIM>& idxNode) const;

    /** Computes the corresponding scalar node number of a given
     * node number of the grid using the FORTRAN style.
     * \param idxNode Given node number.
     * \return Corresponding scalar node number.
     */
    unsigned long int idxNodeTuple2ScalarStyleFortran(
        const ScaFES::Ntuple<int, DIM>& idxNode) const;

    /** Computes the corresponding scalar node number of a given
     * node number of the grid using the FORTRAN style.
     * \param idxScalarNode Corresponding scalar node number (return value).
     * \param idxNode Given node number.
     */
    void idxNodeTuple2ScalarStyleFortran(
        unsigned long int& idxScalarNode,
        const ScaFES::Ntuple<int, DIM>& idxNode) const;

    /** Computes the corresponding node number of a given scalar
     * node number of the grid using the C style.
     * \param idxScalarNode Given scalar node number
     * \return Corresponding node number.
     */
    ScaFES::Ntuple<int, DIM>
    idxNodeScalar2TupleStyleC(const unsigned long int& idxScalarNode) const;

    /** Computes the corresponding node number of a given scalar
     * node number of the grid using the C style.
     * \param idxNode Corresponding node number (return value).
     * \param idxScalarNode Given scalar node number
     */
    void
    idxNodeScalar2TupleStyleC(ScaFES::Ntuple<int, DIM>& idxNode,
                              const unsigned long int& idxScalarNode) const;

    /** Computes the corresponding node number of a given scalar
     * node number of the grid using the Fortran style.
     * \param idxScalarNode Given scalar node number.
     * \return Corresponding node number.
     */
    ScaFES::Ntuple<int, DIM> idxNodeScalar2TupleStyleFortran(
        const unsigned long int& idxScalarNode) const;

    /** Computes the corresponding node number of a given scalar
     * node number of the grid using the Fortran style.
     * \param idxNode Corresponding node number (return value).
     * \param idxScalarNode Given scalar node number.
     */
    void idxNodeScalar2TupleStyleFortran(ScaFES::Ntuple<int, DIM>& idxNode,
                                         const unsigned long int& idxScalarNode)
        const;
}; // End of class. //

/*******************************************************************************
 * FREE METHODS.
 ******************************************************************************/
/** Swaps two objects of the class \c Grid. */
template <std::size_t DIM>
void swap(ScaFES::Grid<DIM>& first, ScaFES::Grid<DIM>& second);
/*----------------------------------------------------------------------------*/
/** Preprares writing an object of the class \c Grid to output. */
template <std::size_t DIM>
std::ostream& operator<<(std::ostream& output, const ScaFES::Grid<DIM>& gg);

/*******************************************************************************
 * LIFE CYCLE METHODS.
 ******************************************************************************/
template <std::size_t DIM>
inline Grid<DIM>::Grid()
: mIdxNodeFirst(0), mIdxNodeLast(1), mNnodes(2), mCoordNodeFirst(0.0),
  mCoordNodeLast(1.0), mGridsize(1.0), mIsNumberedStyleC(false)
{
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline Grid<DIM>::Grid(const ScaFES::Ntuple<int, DIM>& idxNodeFirst,
                       const ScaFES::Ntuple<int, DIM>& idxNodeLast,
                       ScaFES::Ntuple<double, DIM> coordNodeFirst,
                       ScaFES::Ntuple<double, DIM> coordNodeLast)
: mIdxNodeFirst(idxNodeFirst), mIdxNodeLast(idxNodeLast),
  mNnodes(idxNodeLast - idxNodeFirst + ScaFES::Ntuple<int, DIM>(1)),
  mCoordNodeFirst(coordNodeFirst), mCoordNodeLast(coordNodeLast), mGridsize(),
  mIsNumberedStyleC(false)
{
    if (0 < this->nNodes().size())
    {
        for (std::size_t ii = 0; ii < DIM; ++ii)
        {
            this->mGridsize[ii] =
                (coordNodeLast.elem(ii) - coordNodeFirst.elem(ii)) /
                (this->nNodes().elem(ii) - 1);
        }
    }
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline Grid<DIM>::Grid(const Grid<DIM>& from)
: mIdxNodeFirst(from.idxNodeFirst()), mIdxNodeLast(from.idxNodeLast()),
  mNnodes(from.nNodes()), mCoordNodeFirst(from.coordNodeFirst()),
  mCoordNodeLast(from.coordNodeLast()), mGridsize(from.gridsize()),
  mIsNumberedStyleC(from.isNumberedStyleC())
{
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM> inline Grid<DIM>& Grid<DIM>::operator=(Grid<DIM> rhs)
{
    ScaFES::swap(*this, rhs);
    return *this;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM> inline Grid<DIM>::~Grid()
{
}

/*******************************************************************************
 * GETTER METHODS.
 ******************************************************************************/
template <std::size_t DIM> inline std::size_t Grid<DIM>::dimGrid() const
{
    return DIM;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline const ScaFES::Ntuple<int, DIM>& Grid<DIM>::idxNodeFirst() const
{
    return this->mIdxNodeFirst;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline const int& Grid<DIM>::idxNodeFirst(const int& idx) const
{
    return this->mIdxNodeFirst.elem(idx);
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline const ScaFES::Ntuple<int, DIM>& Grid<DIM>::idxNodeLast() const
{
    return this->mIdxNodeLast;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline const int& Grid<DIM>::idxNodeLast(const int& idx) const
{
    return this->mIdxNodeLast.elem(idx);
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline const ScaFES::Ntuple<int, DIM>& Grid<DIM>::nNodes() const
{
    return this->mNnodes;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline const int& Grid<DIM>::nNodes(const int& idx) const
{
    return this->mNnodes.elem(idx);
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM> inline int Grid<DIM>::nNodesTotal() const
{
    return this->mNnodes.size();
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline const ScaFES::Ntuple<double, DIM>& Grid<DIM>::coordNodeFirst() const
{
    return this->mCoordNodeFirst;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline const ScaFES::Ntuple<double, DIM>& Grid<DIM>::coordNodeLast() const
{
    return this->mCoordNodeLast;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline const ScaFES::Ntuple<double, DIM>& Grid<DIM>::gridsize() const
{
    return this->mGridsize;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline const double& Grid<DIM>::gridsize(const int& ii) const
{
    return this->mGridsize.elem(ii);
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline const bool& Grid<DIM>::isNumberedStyleC() const
{
    return this->mIsNumberedStyleC;
}

/*******************************************************************************
 * WORK METHODS.
 ******************************************************************************/
template <std::size_t DIM>
inline typename Grid<DIM>::iterator Grid<DIM>::begin() const
{
    return iterator(this, this->mIdxNodeFirst);
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline typename Grid<DIM>::iterator Grid<DIM>::end() const
{
    return iterator(this, this->mIdxNodeLast);
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline unsigned long int
Grid<DIM>::idxNodeTuple2Scalar(const ScaFES::Ntuple<int, DIM>& idxNode) const
{
    unsigned long int idxScalarNode = 0L;
    this->idxNodeTuple2Scalar(idxScalarNode, idxNode);
    return idxScalarNode;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline void
Grid<DIM>::idxNodeTuple2Scalar(unsigned long int& idxScalarNode,
                               const ScaFES::Ntuple<int, DIM>& idxNode) const
{
    if (this->mIsNumberedStyleC)
    {
        this->idxNodeTuple2ScalarStyleC(idxScalarNode, idxNode);
    }
    else
    {
        this->idxNodeTuple2ScalarStyleFortran(idxScalarNode, idxNode);
    }
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline ScaFES::Ntuple<int, DIM>
Grid<DIM>::idxNodeScalar2Tuple(const unsigned long int& idxScalarNode) const
{
    ScaFES::Ntuple<int, DIM> idxNode;
    this->idxNodeScalar2Tuple(idxNode, idxScalarNode);
    return idxNode;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline void
Grid<DIM>::idxNodeScalar2Tuple(ScaFES::Ntuple<int, DIM>& idxNode,
                               const unsigned long int& idxScalarNode) const
{
    if (this->mIsNumberedStyleC)
    {
        this->idxNodeScalar2TupleStyleC(idxNode, idxScalarNode);
    }
    else
    {
        this->idxNodeScalar2TupleStyleFortran(idxNode, idxScalarNode);
    }
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline ScaFES::Ntuple<double, DIM>
Grid<DIM>::coordinates(const ScaFES::Ntuple<int, DIM>& idxNode) const
{
    ScaFES::Ntuple<double, DIM> coord(0.0);

    for (std::size_t ii = 0; ii < DIM; ++ii)
    {
        coord[ii] = this->coordNodeFirst().elem(ii) +
                    this->gridsize(ii) *
                        (idxNode.elem(ii) - this->idxNodeFirst().elem(ii));
    }

    return coord;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM> inline bool Grid<DIM>::valid() const
{
    return (this->mIdxNodeFirst <= this->mIdxNodeLast);
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline bool Grid<DIM>::inside(const ScaFES::Ntuple<int, DIM>& idxNode) const
{
    return ((this->mIdxNodeFirst <= idxNode) &&
            (idxNode <= this->mIdxNodeLast));
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline Grid<DIM> Grid<DIM>::extend(const int& dir, const int& nLayers)
{
    ScaFES::Ntuple<int, DIM> idxNodeFirst(this->mIdxNodeFirst);
    ScaFES::Ntuple<int, DIM> idxNodeLast(this->mIdxNodeLast);
    int idx = dir / 2;

    if (0 == (dir % 2))
    {
        idxNodeLast[idx] += nLayers;
    }
    else
    {
        idxNodeFirst[idx] -= nLayers;
    }

    return Grid<DIM>(idxNodeFirst, idxNodeLast, this->coordNodeFirst(),
                     this->coordNodeLast());
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline Grid<DIM> Grid<DIM>::unionWith(const Grid<DIM>& rhs) const
{
    ScaFES::Ntuple<double, DIM> coordNodeFirstNew;
    ScaFES::Ntuple<double, DIM> coordNodeLastNew;
    ScaFES::Ntuple<int, DIM> idxNodeFirstNew =
        ScaFES::min(idxNodeFirst(), rhs.idxNodeFirst());
    ScaFES::Ntuple<int, DIM> idxNodeLastNew =
        ScaFES::max(idxNodeLast(), rhs.idxNodeLast());
    if (idxNodeFirstNew == this->idxNodeFirst())
    {
        coordNodeFirstNew = this->coordNodeFirst();
    }
    else
    {
        coordNodeFirstNew = rhs.coordNodeFirst();
    }
    if (idxNodeLastNew == this->idxNodeLast())
    {
        coordNodeLastNew = this->coordNodeLast();
    }
    else
    {
        coordNodeLastNew = rhs.coordNodeLast();
    }
    return ScaFES::Grid<DIM>(idxNodeFirstNew, idxNodeLastNew, coordNodeFirstNew,
                             coordNodeLastNew);
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline Grid<DIM> Grid<DIM>::intersectWith(const Grid<DIM>& rhs) const
{
    ScaFES::Ntuple<double, DIM> coordNodeFirstNew;
    ScaFES::Ntuple<double, DIM> coordNodeLastNew;
    ScaFES::Ntuple<int, DIM> idxNodeFirstNew =
        ScaFES::max(this->idxNodeFirst(), rhs.idxNodeFirst());
    ScaFES::Ntuple<int, DIM> idxNodeLastNew =
        ScaFES::min(this->idxNodeLast(), rhs.idxNodeLast());
    if (idxNodeFirstNew == this->idxNodeFirst())
    {
        coordNodeFirstNew = this->coordNodeFirst();
    }
    else
    {
        coordNodeFirstNew = rhs.coordNodeFirst();
    }
    if (idxNodeLastNew == this->idxNodeLast())
    {
        coordNodeLastNew = this->coordNodeLast();
    }
    else
    {
        coordNodeLastNew = rhs.coordNodeLast();
    }
    return ScaFES::Grid<DIM>(idxNodeFirstNew, idxNodeLastNew, coordNodeFirstNew,
                             coordNodeLastNew);
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline bool Grid<DIM>::operator==(const Grid<DIM>& rhs) const
{
    return (this->idxNodeFirst() == rhs.idxNodeFirst() &&
            this->idxNodeLast() == rhs.idxNodeLast());
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline bool Grid<DIM>::operator!=(const Grid<DIM>& rhs) const
{
    return !(*this == rhs);
}
/*----------------------------------------------------------------------------*/
#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
template <std::size_t DIM>
template <class Archive>
inline void Grid<DIM>::serialize(Archive& ar, const unsigned int version)
{
    if (1 <= version)
    {
        ar&(this->mIdxNodeFirst);
        ar&(this->mIdxNodeLast);
        ar&(this->mNnodes);
        ar&(this->mCoordNodeFirst);
        ar&(this->mCoordNodeLast);
        ar&(this->mGridsize);
        ar&(this->mIsNumberedStyleC);
    }
}
#endif
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline Ntuple<unsigned long int, 2 * DIM>
Grid<DIM>::connect(const ScaFES::Ntuple<int, DIM>& idxNode)
{
    // First, set up connectivity for ALL nodes (also for the boundary nodes).
    // Then, correct the connectivity entries at the boundary nodes.
    unsigned long int gap;
    unsigned long int idxScalarNode;
    ScaFES::Ntuple<unsigned long int, 2 * DIM> neighIdxNode;
    idxScalarNode = this->idxNodeTuple2Scalar(idxNode);
    gap = 1;
    for (std::size_t kk = 0; kk < DIM; ++kk)
    {
        neighIdxNode[2 * kk] = idxScalarNode - gap;
        neighIdxNode[2 * kk + 1] = idxScalarNode + gap;
        // Adapt boundary of hypercube.
        if (idxNode.elem(kk) == this->idxNodeFirst(kk))
        {
            neighIdxNode[2 * kk] = -(2 * kk + 1);
        }
        else if (idxNode.elem(kk) == this->idxNodeLast(kk))
        {
            neighIdxNode[2 * kk + 1] = -(2 * kk + 2);
        }
        gap *= this->nNodes(kk);
    }
    return neighIdxNode;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline ScaFES::Ntuple<int, DIM>
Grid<DIM>::neigh(const ScaFES::Ntuple<int, DIM>& idxNode, const int& dir) const
{
    // Adapt boundary of hypercube.
    ScaFES::Ntuple<int, DIM> neighIdxNode(idxNode);
    int idx = dir / 2;
    if (0 == (dir % 2))
    {
        if (idxNode.elem(idx) == this->idxNodeFirst(idx))
        {
            for (std::size_t ii = 0; ii < DIM; ++ii)
            {
                neighIdxNode[ii] = -(dir + 1);
            }
        }
        else
        {
            --(neighIdxNode[idx]);
        }
    }
    else
    {
        if (idxNode.elem(idx) == this->idxNodeLast(idx))
        {
            for (std::size_t ii = 0; ii < DIM; ++ii)
            {
                neighIdxNode[ii] = -(dir + 1);
            }
        }
        else
        {
            ++(neighIdxNode[idx]);
        }
    }
    return neighIdxNode;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline unsigned long int Grid<DIM>::idxNodeTuple2ScalarStyleC(
    const ScaFES::Ntuple<int, DIM>& idxNode) const
{
    unsigned long int idxScalarNode = 0L;
    this->idxNodeTuple2ScalarStyleC(idxScalarNode, idxNode);
    return idxScalarNode;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline void Grid<DIM>::idxNodeTuple2ScalarStyleC(
    unsigned long int& idxScalarNode,
    const ScaFES::Ntuple<int, DIM>& idxNode) const
{
    idxScalarNode = idxNode.elem(0) - this->idxNodeFirst(0);

    for (std::size_t ii = 1; ii < DIM; ++ii)
    {
        idxScalarNode = idxScalarNode * this->nNodes(ii) + idxNode.elem(ii) -
                        this->idxNodeFirst(ii);
    }
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline unsigned long int Grid<DIM>::idxNodeTuple2ScalarStyleFortran(
    const ScaFES::Ntuple<int, DIM>& idxNode) const
{
    unsigned long int idxScalarNode = 0L;
    this->idxNodeTuple2ScalarStyleFortran(idxScalarNode, idxNode);
    return idxScalarNode;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline void Grid<DIM>::idxNodeTuple2ScalarStyleFortran(
    unsigned long int& idxScalarNode,
    const ScaFES::Ntuple<int, DIM>& idxNode) const
{
    idxScalarNode = idxNode.elem(DIM - 1) - this->idxNodeFirst(DIM - 1);
    const int dimLocal = DIM;

    for (int ii = dimLocal - 2; ii >= 0; --ii)
    {
        idxScalarNode = idxScalarNode * this->nNodes(ii) + idxNode.elem(ii) -
                        this->idxNodeFirst(ii);
    }
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline ScaFES::Ntuple<int, DIM> Grid<DIM>::idxNodeScalar2TupleStyleC(
    const unsigned long int& idxScalarNode) const
{
    ScaFES::Ntuple<int, DIM> idxNode;
    this->idxNodeScalar2TupleStyleC(idxNode, idxScalarNode);
    return idxNode;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline void Grid<DIM>::idxNodeScalar2TupleStyleC(
    ScaFES::Ntuple<int, DIM>& idxNode,
    const unsigned long int& idxScalarNode) const
{
    int tmp = idxScalarNode;
    for (std::size_t ii = DIM - 1; ii > 0; --ii)
    {
        idxNode[ii] = (tmp % this->nNodes(ii)) + this->idxNodeFirst(ii);
        tmp /= (this->nNodes(ii));
    }
    idxNode[0] = tmp;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline ScaFES::Ntuple<int, DIM> Grid<DIM>::idxNodeScalar2TupleStyleFortran(
    const unsigned long int& idxScalarNode) const
{
    ScaFES::Ntuple<int, DIM> idxNode;
    this->idxNodeScalar2TupleStyleFortran(idxNode, idxScalarNode);
    return idxNode;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline void Grid<DIM>::idxNodeScalar2TupleStyleFortran(
    ScaFES::Ntuple<int, DIM>& idxNode,
    const unsigned long int& idxScalarNode) const
{
    int tmp = idxScalarNode;
    for (std::size_t ii = 0; ii < DIM; ++ii)
    {
        idxNode[ii] = (tmp % this->nNodes(ii)) + this->idxNodeFirst(ii);
        tmp /= (this->nNodes(ii));
    }
}

/*******************************************************************************
 * INLINED OPERATORS (FREE FUNCTIONS)
 ******************************************************************************/
template <std::size_t DIM> inline void swap(Grid<DIM>& first, Grid<DIM>& second)
{
    ScaFES::swap(first.mIdxNodeFirst, second.mIdxNodeFirst);
    ScaFES::swap(first.mIdxNodeLast, second.mIdxNodeLast);
    ScaFES::swap(first.mNnodes, second.mNnodes);
    ScaFES::swap(first.mCoordNodeFirst, second.mCoordNodeFirst);
    ScaFES::swap(first.mCoordNodeLast, second.mCoordNodeLast);
    ScaFES::swap(first.mGridsize, second.mGridsize);
    std::swap(first.mIsNumberedStyleC, second.mIsNumberedStyleC);
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline std::ostream& operator<<(std::ostream& output, const Grid<DIM>& gg)
{
    output << "* GRID:" << std::endl << "  -Node First:" << gg.idxNodeFirst()
           << std::endl << "  -Node Last: " << gg.idxNodeLast() << std::endl
           << "  -Coordinates(Node First): " << gg.coordNodeFirst() << std::endl
           << "  -Coordinates(Node Last):  " << gg.coordNodeLast() << std::endl;
    return output;
}

/*******************************************************************************
 * INTERNAL CLASS ITERATOR: LIFE CYCLE METHODS.
 ******************************************************************************/
template <std::size_t DIM>
inline Grid<DIM>::Iterator::Iterator(const Grid<DIM>* gridIter,
                                     const ScaFES::Ntuple<int, DIM>& idxNode)
: mGridIter(gridIter), mIdxNode(idxNode),
  mIdxScalarNode(gridIter->idxNodeTuple2Scalar(idxNode))
{
    // TODO: Check if \c idxNode is inside underlying grid.
}

/*******************************************************************************
 * INTERNAL CLASS ITERATOR: GETTER METHODS.
 ******************************************************************************/
template <std::size_t DIM>
inline const ScaFES::Ntuple<int, DIM>& Grid<DIM>::Iterator::idxNode() const
{
    return this->mIdxNode;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline const unsigned long int& Grid<DIM>::Iterator::idxScalarNode() const
{
    return this->mIdxScalarNode;
}

/*******************************************************************************
 * INTERNAL CLASS ITERATOR: COMPARISON METHODS.
 ******************************************************************************/
template <std::size_t DIM>
inline bool Grid<DIM>::Iterator::operator!=(const Iterator& other) const
{
    return (*this < other) || (other < *this);
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline bool Grid<DIM>::Iterator::operator<(const Iterator& it) const
{
    return (this->idxScalarNode() <= it.idxScalarNode());
}

/*******************************************************************************
 * INTERNAL CLASS ITERATOR: WORK METHODS.
 ******************************************************************************/
template <std::size_t DIM>
inline const ScaFES::Ntuple<int, DIM>& Grid<DIM>::Iterator::operator*()
{
    return this->idxNode();
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline typename Grid<DIM>::Iterator& Grid<DIM>::Iterator::operator++()
{
    if (this->mGridIter->isNumberedStyleC())
    {
        return this->nextStyleC();
    }

    return this->nextStyleFortran();
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline typename Grid<DIM>::Iterator Grid<DIM>::Iterator::operator++(int /*incr*/)
{
    Iterator tmp = *this;
    ++(*this);
    return tmp;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline typename Grid<DIM>::Iterator& Grid<DIM>::Iterator::nextStyleC()
{
    ++(this->mIdxScalarNode);
    ++(this->mIdxNode[DIM - 1]);
    const int DIM_MINUS_ONE = DIM - 1;

    for (int ii = DIM_MINUS_ONE; ii > 0; --ii)
    {
        if (this->mIdxNode.elem(ii) > this->mGridIter->idxNodeLast(ii))
        {
            this->mIdxNode[ii] = this->mGridIter->idxNodeFirst(ii);
            ++(this->mIdxNode[ii - 1]);
        }
    }
    return *this;
}
/*----------------------------------------------------------------------------*/
template <std::size_t DIM>
inline typename Grid<DIM>::Iterator& Grid<DIM>::Iterator::nextStyleFortran()
{
    ++(this->mIdxScalarNode);
    ++(this->mIdxNode[0]);
    const int DIM_MINUS_ONE = DIM - 1;

    for (int ii = 0; ii < DIM_MINUS_ONE; ++ii)
    {
        if (this->mIdxNode.elem(ii) > this->mGridIter->idxNodeLast(ii))
        {
            this->mIdxNode[ii] = this->mGridIter->idxNodeFirst(ii);
            ++(this->mIdxNode[ii + 1]);
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
    template <std::size_t DIM> struct version<ScaFES::Grid<DIM>>
    {
        /** Sets the version number for serialization. */
        BOOST_STATIC_CONSTANT(unsigned long int, value = 2);
    };
} // namespace serialization
} // namespace boost
#endif

#endif
