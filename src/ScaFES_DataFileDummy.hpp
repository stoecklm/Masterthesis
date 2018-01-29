/* ScaFES
 * Copyright (c) 2011-2017, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file ScaFES_DataFileDummy.hpp
 *  @brief Contains the class template Datafile with a dummy implementation.
 */

#ifndef SCAFES_DATAFILEDUMMY_HPP_
#define SCAFES_DATAFILEDUMMY_HPP_

#include "ScaFES_Config.hpp"

#ifdef VTRACE
#include "vt_user.h"
#endif

#include <iostream>
#include <vector>
#include <tuple>

#include "ScaFES_Communicator.hpp"
#include "ScaFES_Ntuple.hpp"
#include "ScaFES_GridSub.hpp"

namespace ScaFES
{

/** Version number of data file. */
const int DATA_FILE_VERSION = 3;

/*******************************************************************************
 ******************************************************************************/
/**
* \class DataFile
* @brief The class template \c DataFile is responsible for writing
* (per default: in parallel) a data field to a given file in the NetCDf
* output format.
*
* This is the dummy version of the class template.
*/
template <typename TT, std::size_t DIM> class DataFile
{
public:
    /*----------------------------------------------------------------------
    | TYPE DEFINITIONS.
    ----------------------------------------------------------------------*/
    /** Re-export typename TT STL-like. */
    typedef TT value_type;

    /** Type definition for integer d-dimensional tuple. */
    typedef Ntuple<int, DIM> IntNtuple;

    /*----------------------------------------------------------------------
    | LIFE CYCLE METHODS.
    ----------------------------------------------------------------------*/
    /** Creates the default constructor. */
    DataFile();

    /** Creates own constructor.
     *  \remarks This constructor only sets the names of the file and the
     *  variable which should be written, the file itself WILL NOT be
     *  created here. */
    DataFile(std::string const& nameDataFile,
             std::vector<std::string> const& nameVariables,
             ScaFES::Communicator const& myWorld, IntNtuple const& nNodes,
             std::vector<GridSub<DIM>> const& memNormal);

    /** Creates default copy constructor. */
    DataFile(DataFile<TT, DIM> const&) = default;

    /** Creates default assignment operator. */
    DataFile& operator=(DataFile<TT, DIM> const& /*rhs*/) = default;

    /** Creates own destructor: Closes the file. */
    ~DataFile();

    /*----------------------------------------------------------------------
    | GETTER METHODS.
    ----------------------------------------------------------------------*/
    /** Returns the name of the data file. */
    std::string const& nameDataFile() const;

    /** Returns the name of the variable which should be written. */
    std::vector<std::string> const& nameVariables() const;

    /** Returns the identifier of the file. */
    int const& idFile() const;

    /** Returns the identifiers of the components of the variable. */
    std::vector<std::vector<int>> const& idVariable() const;

    /** Returns if the file is already open or not. */
    bool const& isOpen() const;

    /** Returns the number the data field was written to the file. */
    int const& nWritesCurr() const;

    /** Returns the memory grid. */
    std::vector<GridSub<DIM>> const& memNormal() const;

    /** Returns a pointer to the data which should be written. */
    std::vector<TT*> const& elemData() const;

    /** Returns the data for a given memory position. */
    TT const& elemData(int const& df, int const& idx) const;

    /** Returns the number of elements stored in the memory \c mNcMemory. */
    int const& nElemsNcMem() const;

    /*----------------------------------------------------------------------
    | WORK METHODS.
    ----------------------------------------------------------------------*/
    /** Creates the data file and sets the file layout of the
     *  underlying data field (depending on the template type). */
    void create();

    /** Writes the elements in the memory to which the pointer points to
     *  the file at the given time step. */
    void write(std::vector<TT*> const& elemData,
               std::vector<bool> const& writeToFile, int const& timestep);

    /** Reads the elements in the memory to which the pointer points to
     *  from the file at the given time step. */
    void read(std::vector<TT*>& elemData, int const& timestep);

    /** Initializes the elements in the memory to which the pointer points to
     *  from the file. */
    void init(std::vector<TT*>& elemData);

private:
    /*----------------------------------------------------------------------
    | TYPE DEFINITIONS.
    ----------------------------------------------------------------------*/
    /** Type of error variables of NetCDF lib. */
    typedef int NCErrorType;

    /*----------------------------------------------------------------------
    | STATIC VARIABLES AND METHODS.
    ----------------------------------------------------------------------*/
    /** A map containing all error messages. */
    static std::unordered_map<NCErrorType, std::string> mNCErrorStrings;

    /** Handles all errors in conjunction with NetCDF.
    * \param status Return value of NetCDF method.
    * \param file Name of file where error occurs.
    * \param line Number of line where error occurs.
    */
    static void wrapNCCall(NCErrorType const& status, std::string const& file,
                           int const& line);

    /*----------------------------------------------------------------------
    | MEMBER VARIABLES.
    ----------------------------------------------------------------------*/
    /** Name of NetCDF data file. */
    std::string mNameDataFile;

    /** Name of the data field which should be written to the file. */
    std::vector<std::string> mNameVariables;

    /** MPI communicator. */
    ScaFES::Communicator mMyWorld;

    /** Number of nodes in each dimension. */
    IntNtuple mNnodes;

    /** Identifier of a NetCDF dataset.
     * "A netCDF ID that can subsequently be used to refer to the netCDF
     * dataset in other netCDF function calls."
     * [http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-c.html#nc_005fcreate]
     */
    int mIdFile;

    /** Every variable has a unique identifier. */
    std::vector<std::vector<int>> mIdVariable;

    /** Identifiers of space and time dimensions
     * = idNetCDF + nameDimension + lengthDimension. */
    // int mIdDims[DIM + 1];
    std::vector<int> mIdDims;

    /** Is file open? */
    bool mIsOpen;

    /** Number of writes. */
    int mNwritesCurr;

    /** Memory grid = Normal grid partition of each data field. */
    std::vector<GridSub<DIM>> mMemNormal;

    /** Pointers to the memory of all data fields. */
    std::vector<TT*> mElemData;

    /** Number of elements which will be allocated. */
    int mNelemsNcMem;

    /** Pointer to linear continuous memory for the elements of
     * a data field.
     * \remarks Memory for all normal grid nodes will be allocated.
     * NetCDF needs the memory at all normal grid nodes.
     * Thus, the memory for all normal grid nodes will be allocated.
     * In contrast to this, elemData is a pointer to the memory at
     * really all nodes including the halos.
     */
    std::vector<double> mNcMemory;

}; // End of class. //

/*******************************************************************************
 * DECLARATION OF FREE METHODS.
 ******************************************************************************/
/** Free function that handles all errors in conjunction with NetCDF.
 * \param status Return value of NetCDF method.
 * \param file Name of file where error occurs.
 * \param line Number of line where error occurs.
 */
void handle_netCDF_error(int const& status, std::string const& file,
                         int const& line);

/*******************************************************************************
 * LIFE CYCLE METHODS.
 ******************************************************************************/
template <typename TT, std::size_t DIM>
inline DataFile<TT, DIM>::DataFile()
: mNameDataFile("EMPTYFILE.nc"), mNameVariables(), mMyWorld(), mNnodes(0),
  mIdFile(0), mIdVariable(), mIdDims(), mIsOpen(false), mNwritesCurr(0),
  mMemNormal(), mElemData(), mNelemsNcMem(0), mNcMemory()
{
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM>
inline DataFile<TT, DIM>::DataFile(
    std::string const& nameDataFile,
    std::vector<std::string> const& nameVariables,
    ScaFES::Communicator const& myWorld, IntNtuple const& nNodes,
    std::vector<GridSub<DIM>> const& memNormal)
: mNameDataFile(), mNameVariables(nameVariables), mMyWorld(myWorld),
  mNnodes(nNodes), mIdFile(0), mIdVariable(), mIdDims(DIM + 1, 0),
  mIsOpen(false), mNwritesCurr(0), mMemNormal(memNormal),
  mElemData(nameVariables.size(), static_cast<TT*>(0)), mNelemsNcMem(0),
  mNcMemory()
{
    std::ostringstream tmpStringstream;
    tmpStringstream << nameDataFile << ".nc";
    this->mNameDataFile = tmpStringstream.str();
    this->mIdVariable.resize(nameVariables.size());
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM> inline DataFile<TT, DIM>::~DataFile()
{
}

/*******************************************************************************
 * GETTER METHODS.
 ******************************************************************************/
template <typename TT, std::size_t DIM>
inline std::string const& DataFile<TT, DIM>::nameDataFile() const
{
    return mNameDataFile;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM>
inline std::vector<std::string> const& DataFile<TT, DIM>::nameVariables() const
{
    return mNameVariables;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM>
inline int const& DataFile<TT, DIM>::idFile() const
{
    return mIdFile;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM>
inline std::vector<std::vector<int>> const&
DataFile<TT, DIM>::idVariable() const
{
    return mIdVariable;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM>
inline bool const& DataFile<TT, DIM>::isOpen() const
{
    return mIsOpen;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM>
inline int const& DataFile<TT, DIM>::nWritesCurr() const
{
    return mNwritesCurr;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM>
inline std::vector<GridSub<DIM>> const& DataFile<TT, DIM>::memNormal() const
{
    return mMemNormal;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM>
inline std::vector<TT*> const& DataFile<TT, DIM>::elemData() const
{
    return mElemData;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM>
inline TT const& DataFile<TT, DIM>::elemData(int const& jj,
                                             int const& idx) const
{
    return mElemData[jj][idx];
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM>
inline int const& DataFile<TT, DIM>::nElemsNcMem() const
{
    return mNelemsNcMem;
}

/*******************************************************************************
 * WORK METHODS.
 ******************************************************************************/
template <typename TT, std::size_t DIM> void DataFile<TT, DIM>::create()
{
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM>
void DataFile<TT, DIM>::write(std::vector<TT*> const& /*elemData*/,
                              std::vector<bool> const& /*writeToFile*/,
                              int const& /*timestep*/
                              )
{
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM>
void DataFile<TT, DIM>::read(std::vector<TT*>& /*elemData*/,
                             int const& /*timestep*/
                             )
{
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM>
void DataFile<TT, DIM>::init(std::vector<TT*>& /*elemData*/
                            )
{
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM>
void DataFile<TT, DIM>::wrapNCCall(int const& status, std::string const& file,
                                   int const& line)
{
    std::ignore = status;
    std::ignore = file;
    std::ignore = line;
}

} // End of namespace. //

#endif
