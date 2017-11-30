/* ScaFES
 * Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file ScaFES_DataFileApi.hpp
 *  @brief Contains the class template Datafile.
 */

#ifndef SCAFES_DATAFILEAPI_HPP_
#define SCAFES_DATAFILEAPI_HPP_

#include "ScaFES_Config.hpp"


// MPI must be available for parallel NetCDF.
// Boost.MPI is not necessary for parallel NetCDF.
#ifdef SCAFES_HAVE_MPI
    #ifdef SCAFES_HAVE_NETCDF_PAR
        #include <netcdf_par.h>
    #endif
#endif
#include <netcdf.h>

#ifdef VTRACE
#include "vt_user.h"
#endif

#include <iostream>
#include <vector>
#include <unordered_map>

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
* The call to the NetCDF creation of the data file is not set
* within the constructor.
* Instead, it will be created if the user calls the method \c write() the
* first time. Thus, in case that no data will be written to the output
* (using the method \c write()), then no data file will be created, too.
*
* \n The command \code ncdump -h <netcdffile>\endcode provides some information
* about the stored data fields and parameters.
* \code login@host:~/> ncdump -h out.nc
  dimensions:
       Nx = 130 ;
       Ny = 400 ;
       Nz = 200 ;
       time = UNLIMITED ; // (52 currently)
  variables:
       int version ;
       double D_0(time, Nz, Ny, Nx) ;
       double D_1(time, Nz, Ny, Nx) ;
       double D_2(time, Nz, Ny, Nx) ;
       int writeSteps ;
               writeSteps:_FillValue = 10 ;
       int timesteps ;
               timesteps:_FillValue = 500 ;
       double tau ;
               tau:_FillValue = 1.e-07 ;
 \endcode
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

    /** Inits the elements in the memory to which the pointer points to
     *  from the file. */
    void init(std::vector<TT*>& elemData);

private:
    /*----------------------------------------------------------------------
    | TYPE DEFINITIONS.
    ----------------------------------------------------------------------*/
    /** Type of error variables of NetCDF lib. */
    typedef decltype(NC_NOERR) NCErrorType;

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

mNCErrorStrings[NC_EBADID] = "Not a NetCDF id.";
mNCErrorStrings[NC_ENFILE] = "Too many NetCDF files open";
mNCErrorStrings[NC_EEXIST] = "NetCDF file exists, and NC_NOCLOBBER is set.";
mNCErrorStrings[NC_EINVAL] = "Invalid Argument." ;
mNCErrorStrings[NC_EPERM] = "Trying to write to read-only file.";
mNCErrorStrings[NC_ENOTINDEFINE] = "Operation not allowed in data mode.";
mNCErrorStrings[NC_EINDEFINE] = "Operation not allowed in define mode.";
mNCErrorStrings[NC_EINVALCOORDS] = "Index exceeds dimension bound.";
mNCErrorStrings[NC_EMAXDIMS] = "NC_MAX_DIMS exceeded.";
mNCErrorStrings[NC_ENAMEINUSE] = "String matches name in use.";
mNCErrorStrings[NC_ENOTATT] = "Attribute not found.";
mNCErrorStrings[NC_EMAXATTS] = "NC_MAX_ATTRS exceeded.";
mNCErrorStrings[NC_EBADTYPE] = "Not a NetCDF data type.";
mNCErrorStrings[NC_EBADDIM] = "Invalid dimension id or name.";
mNCErrorStrings[NC_EUNLIMPOS] = "NC_UNLIMITED in wrong index.";
mNCErrorStrings[NC_EMAXVARS] = "NC_MAX_VARS exceeded.";
mNCErrorStrings[NC_ENOTVAR] = "Variable not found";
mNCErrorStrings[NC_EGLOBAL] = "Action prohibited on NC_GLOBAL variable id.";
mNCErrorStrings[NC_ENOTNC] = "Not a NetCDF file.";
mNCErrorStrings[NC_EMAXNAME] = "NC_MAX_NAME exceeded.";
mNCErrorStrings[NC_EUNLIMIT] = "NC_UNLIMITED size already in use.";
mNCErrorStrings[NC_ENORECVARS] = "nc_rec op when there are no record vars.";
mNCErrorStrings[NC_ECHAR] = "Attempt to convert between text and numbers.";
mNCErrorStrings[NC_EEDGE] = "Edge+start exceeds dimension bound.";
mNCErrorStrings[NC_ESTRIDE] = "Illegal stride.";
mNCErrorStrings[NC_EBADNAME] = "Attribute or variable name contains illegal characters.";
mNCErrorStrings[NC_ERANGE] = "Math result not representable.";
mNCErrorStrings[NC_ENOMEM] = "Memory allocation failure.";
mNCErrorStrings[NC_EVARSIZE] = "One or more variable sizes violate format constraints.";
mNCErrorStrings[NC_EDIMSIZE] = "Invalid dimension size";
mNCErrorStrings[NC_ETRUNC] = "File likely truncated or possibly corrupted.";
mNCErrorStrings[NC_EAXISTYPE] = "Unknown axis type.";


}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM> inline DataFile<TT, DIM>::~DataFile()
{
    if (this->mIsOpen)
    {
        if (0 < this->mIdFile)
        {
            wrapNCCall(nc_close(this->mIdFile), __FILE__, __LINE__);
        }
    }
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
#ifdef VTRACE
    VT_TRACER("DataFile::create()");
#endif
    /*------------------------------------------------------------------------*/
    int version = DATA_FILE_VERSION;
    int versionID;
    int status; // Return value of netCDF methods.
                // Create netCDF dataset.

#ifdef SCAFES_HAVE_MPI
    #ifdef SCAFES_HAVE_NETCDF_PAR
    // Flag NC_NETCDF4 specifies that HDF5 library is used for parallel I/O //
    MPI_Info info = MPI_INFO_NULL;
    status = nc_create_par(mNameDataFile.c_str(), NC_NETCDF4 | NC_MPIIO,
                           mMyWorld.myWorld(), info, &mIdFile);

    wrapNCCall(status, __FILE__, __LINE__);
    #else
    status = nc_create(nameDataFile().c_str(), NC_FORMAT_NETCDF4, &mIdFile);

    wrapNCCall(status, __FILE__, __LINE__);
    #endif
#else // Case MPI not available.
    status = nc_create(nameDataFile().c_str(), NC_FORMAT_NETCDF4, &mIdFile);

    wrapNCCall(status, __FILE__, __LINE__);
#endif

    this->mIsOpen = true;

    // Add the new variable "version" to an open netCDF dataset in define mode.
    // Variable is of type NC_INT and a scalar.
    wrapNCCall(nc_def_var(mIdFile, "version", NC_INT, 0, 0, &versionID),
               __FILE__, __LINE__);

    /*------------------------------------------------------------------------*/
    // Create the local memory mNcMemory for writing the data to the file once:
    // Here: During the creation process.
    mNelemsNcMem = memNormal().at(0).nNodesTotal();

    if (0 < mNelemsNcMem)
    {
        // Initialise memory at all nodes in order to be sure that
        // there is no dump in the memory.
        mNcMemory = std::vector<double>(mNelemsNcMem, 0.0);
    }

    /*------------------------------------------------------------------------*/
    // Add a new dimension to an open netCDF dataset.
    std::string tmpName = "nNodes";

    for (std::size_t ii = 0; ii < DIM; ++ii)
    {
        std::ostringstream oss;
        oss << tmpName << "_" << ii;
        std::string nameOfDim(oss.str());
        status = nc_def_dim(mIdFile, nameOfDim.c_str(), mNnodes.elem(ii),
                            &mIdDims[DIM - ii]);
        wrapNCCall(status, __FILE__, __LINE__);
    }

    status = nc_def_dim(mIdFile, "time", NC_UNLIMITED, &mIdDims[0]);
    wrapNCCall(status, __FILE__, __LINE__);

    int vID;

    int tmpIdDims[DIM + 1];
    for (std::size_t ii = 0; ii < DIM + 1; ++ii)
    {
        tmpIdDims[ii] = mIdDims[ii];
    }
    for (std::size_t ii = 0; ii < mNameVariables.size(); ++ii)
    {
        status = nc_def_var(mIdFile, (mNameVariables[ii]).c_str(), NC_DOUBLE,
                            DIM + 1, tmpIdDims, &vID);

        wrapNCCall(status, __FILE__, __LINE__);
        mIdVariable[ii].push_back(vID);
    }

    wrapNCCall(nc_enddef(mIdFile), __FILE__, __LINE__);

    /*------------------------------------------------------------------------*/
    // Write the version number to file.
    wrapNCCall(nc_put_var_int(mIdFile, versionID, &version), __FILE__,
               __LINE__);
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM>
void DataFile<TT, DIM>::write(std::vector<TT*> const& elemData,
                              std::vector<bool> const& writeToFile,
                              int const& /*timestep*/
                              )
{
#ifdef VTRACE
    VT_TRACER("DataFile::write()");
#endif

    int status = 0;

    // Update pointers to memory of data fiels which should be written to file.
    // Remark: Pointers must be updated each time step because the old
    // and the new iterate are swapped in each time step.
    for (std::size_t ii = 0; ii < elemData.size(); ++ii)
    {
        this->mElemData[ii] = elemData[ii];
    }

    // Create the file if it was not done before.
    if (!(this->mIsOpen))
    {
        this->create();
#ifdef SCAFES_HAVE_MPI
    #ifdef SCAFES_HAVE_NETCDF_PAR
        for (std::size_t ii = 0; ii < this->nameVariables().size(); ++ii)
        {
            status = nc_var_par_access(this->mIdFile, this->mIdVariable[ii][0],
                                       NC_COLLECTIVE);
            wrapNCCall(status, __FILE__, __LINE__);
        }
    #endif
#endif
    }

    const int nDimensions = 1 + DIM;
    IntNtuple idxNodeFirst = memNormal().at(0).idxNodeFirstSub();
    IntNtuple nNodes = memNormal().at(0).nNodesSub();

    // Set dimensions and transform member variables to NC input format.
    std::size_t elemFirst[nDimensions];
    std::size_t nElems[nDimensions];
    elemFirst[0] = mNwritesCurr; // Has to start with 0.
    ++(this->mNwritesCurr);

    for (int ii = 1; ii < nDimensions; ++ii)
    {
        elemFirst[ii] = (unsigned long int)idxNodeFirst.elem(DIM - ii);
    }

    nElems[0] = 1;
    for (int ii = 1; ii < nDimensions; ++ii)
    {
        nElems[ii] = (unsigned long int)nNodes.elem(DIM - ii);
    }

    double* ptrToMem = 0;
    TT value = 0;

    for (std::size_t ii = 0; ii < this->nameVariables().size(); ++ii)
    {
        // Reset pointer to memory for each data field.
        ptrToMem = this->mNcMemory.data();
        int kk = 0;
        if (writeToFile.at(ii))
        {
            //#pragma omp parallel for
            for (typename GridSub<DIM>::iterator
                     it = this->memNormal().at(ii).begin(),
                     et = this->memNormal().at(ii).end();
                 it < et; ++it)
            {
                // Get real number and store it in
                // mNcMemory.data() = ptrToMem[kk].
                value = this->elemData(ii, it.idxScalarNode());
                ptrToMem[kk] = static_cast<double>(value);
                ++kk;
            }
            // mIdVariable[0] means real_part!
            status =
                nc_put_vara_double(this->mIdFile, this->mIdVariable[ii][0],
                                   elemFirst, nElems, this->mNcMemory.data());
            wrapNCCall(status, __FILE__, __LINE__);
        }
    }
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM>
void DataFile<TT, DIM>::read(std::vector<TT*>& elemData, int const& /*timestep*/
                             )
{
    int idVersion;
    int version = 0;
    int status; // Return value of netCDF methods.

    const int nDimensions = 1 + DIM;
    std::size_t nElems[nDimensions];

    for (std::size_t ii = 0; ii < elemData.size(); ++ii)
    {
        this->mElemData[ii] = elemData[ii];
    }

#ifdef VTRACE
    VT_TRACER("DataFile::read()");
#endif

    /*------------------------------------------------------------------------*/
    // Open the file with read-only access if it was not done before.
    if (!(this->mIsOpen))
    {
        status = nc_open(mNameDataFile.c_str(), NC_NOWRITE, &mIdFile);
        wrapNCCall(status, __FILE__, __LINE__);
        mIsOpen = true;
    }

    /*------------------------------------------------------------------------*/
    // Read in version number.
    status = nc_get_att_int(mIdFile, NC_GLOBAL, "version", &version);
    if (status != NC_NOERR)
    {
        status = nc_inq_varid(mIdFile, "version", &idVersion);
        if (status != NC_NOERR)
        {
            version = 0;
        }
        else
        {
            nc_get_var_int(mIdFile, idVersion, &version);
        }
    }

    /*------------------------------------------------------------------------*/
    // Get the varid of the data variable, based on its name.
    for (std::size_t ii = 0; ii < elemData.size(); ++ii)
    {
        status = nc_inq_varid(mIdFile, mNameVariables[ii].c_str(),
                              &(mIdVariable[ii][0]));
        wrapNCCall(status, __FILE__, __LINE__);
    }

    /*------------------------------------------------------------------------*/
    // Get the var id and the lengths of the grid and time dimensions.
    std::string tmpName = "nNodes";
    for (std::size_t ii = 0; ii < DIM; ++ii)
    {
        std::ostringstream oss;
        oss << tmpName << "_" << ii;
        std::string nameOfDim(oss.str());
        status = nc_inq_dimid(mIdFile, nameOfDim.c_str(), &mIdDims[DIM - ii]);

        wrapNCCall(status, __FILE__, __LINE__);
        status = nc_inq_dimlen(mIdFile, mIdDims[DIM - ii], &(nElems[ii]));

        wrapNCCall(status, __FILE__, __LINE__);

        if (memNormal().at(0).nNodes().elem(ii) != static_cast<int>(nElems[ii]))
        {
            std::cerr << " Grid dimensions does not fit to given grid. "
                      << std::endl;
            break;
        }
    }

    wrapNCCall(nc_inq_dimid(mIdFile, "time", &mIdDims[0]), __FILE__, __LINE__);

    status = nc_inq_dimlen(mIdFile, mIdDims[0], &(nElems[nDimensions - 1]));

    wrapNCCall(status, __FILE__, __LINE__);

    /*------------------------------------------------------------------------*/
    /* Read data. */
    double* ptrToMem = 0;
    double value = 0.0;
    // Get real number and store it in mNcMemory.
    double tmpVal = 0.0;

    for (std::size_t ii = 0; ii < this->nameVariables().size(); ++ii)
    {
        // mIdVariable[0] means real_part!
        status = nc_get_var_double(this->mIdFile, this->mIdVariable[ii][0],
                                   this->mNcMemory.data());
        wrapNCCall(status, __FILE__, __LINE__);

        ptrToMem = this->mNcMemory.data();
        //#pragma omp parallel for
        for (typename GridSub<DIM>::iterator
                 it = this->memNormal().at(ii).begin(),
                 et = this->memNormal().at(ii).end();
             it < et; ++it)
        {
            tmpVal = *ptrToMem;
            value = tmpVal;
            // Will be written to elemData.
            this->mElemData[ii][it.idxScalarNode()] = static_cast<TT>(value);
            ++ptrToMem;
        }
    }
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM>
void DataFile<TT, DIM>::init(std::vector<TT*>& elemData)
{
    int idVersion;
    int version = 0;
    int status; // Return value of netCDF methods.

    const int nDimensions = 1 + DIM;
    std::size_t nElems[nDimensions];
    std::size_t elemFirst[nDimensions];
    elemFirst[0] = 0; // elemFirst[0] = timestep, in this case time = 0

    for (std::size_t ii = 0; ii < elemData.size(); ++ii)
    {
        this->mElemData[ii] = elemData[ii];
    }

#ifdef VTRACE
    VT_TRACER("DataFile::init()");
#endif

    /*------------------------------------------------------------------------*/
    // Open the file with read-only access if it was not done before.
    if (!(this->mIsOpen))
    {
#ifdef SCAFES_HAVE_MPI
    #ifdef SCAFES_HAVE_NETCDF_PAR
        MPI_Info info = MPI_INFO_NULL;
        status = nc_open_par(mNameDataFile.c_str(), NC_NOWRITE | NC_MPIIO,
                             mMyWorld.myWorld(), info, &mIdFile);
    #else
        std::cerr << "\nERROR: Cannot read file in parallel without "
                  << "parallel netCDF (PnetCDF)." << std::endl;
        status = nc_open(mNameDataFile.c_str(), NC_NOWRITE, &mIdFile);
    #endif
#else
        status = nc_open(mNameDataFile.c_str(), NC_NOWRITE, &mIdFile);
#endif
        wrapNCCall(status, __FILE__, __LINE__);
        mIsOpen = true;
    }

    /*------------------------------------------------------------------------*/
    // Read in version number.
    status = nc_get_att_int(mIdFile, NC_GLOBAL, "version", &version);
    if (status != NC_NOERR)
    {
        status = nc_inq_varid(mIdFile, "version", &idVersion);
        if (status != NC_NOERR)
        {
            version = 0;
        }
        else
        {
            nc_get_var_int(mIdFile, idVersion, &version);
        }
    }

    /*------------------------------------------------------------------------*/
    // Initialise memory at all nodes in order to be sure that
    // there is no dump in the memory.
    mNelemsNcMem = memNormal().at(0).nNodesTotal();

    if (0 < mNelemsNcMem)
    {
        mNcMemory = std::vector<double>(mNelemsNcMem, 0.0);
    }

    /*------------------------------------------------------------------------*/
    // Get the varid of the data variable, based on its name.
    for (std::size_t ii = 0; ii < elemData.size(); ++ii)
    {
        mIdVariable[ii].reserve(1);
        status = nc_inq_varid(mIdFile, mNameVariables[ii].c_str(),
                              &(mIdVariable[ii][0]));
        wrapNCCall(status, __FILE__, __LINE__);
    }

    /*------------------------------------------------------------------------*/
    // Get the var id and the lengths of the grid and time dimensions.
    std::string tmpName = "nNodes";
    for (std::size_t ii = 0; ii < DIM; ++ii)
    {
        std::ostringstream oss;
        oss << tmpName << "_" << ii;
        std::string nameOfDim(oss.str());
        status = nc_inq_dimid(mIdFile, nameOfDim.c_str(), &mIdDims[DIM - ii]);

        wrapNCCall(status, __FILE__, __LINE__);
        status = nc_inq_dimlen(mIdFile, mIdDims[DIM - ii], &(nElems[ii]));

        wrapNCCall(status, __FILE__, __LINE__);

        if (memNormal().at(0).nNodes().elem(ii) != static_cast<int>(nElems[ii]))
        {
            std::cerr << " Grid dimensions does not fit to given grid. "
                      << std::endl;
            break;
        }
    }
    status = nc_inq_dimid(mIdFile, "time", &mIdDims[0]);
    wrapNCCall(status, __FILE__, __LINE__);

    status = nc_inq_dimlen(mIdFile, mIdDims[0], &(nElems[nDimensions - 1]));
    wrapNCCall(status, __FILE__, __LINE__);

    /*------------------------------------------------------------------------*/
    // Set dimensions and transform member variables to NC input format.
    IntNtuple idxNodeFirst = memNormal().at(0).idxNodeFirstSub();
    IntNtuple nNodes = memNormal().at(0).nNodesSub();

    for (int ii = 1; ii < nDimensions; ++ii)
    {
        elemFirst[ii] = (unsigned long int)idxNodeFirst.elem(DIM - ii);
    }

    nElems[0] = 1;
    for (int ii = 1; ii < nDimensions; ++ii)
    {
        nElems[ii] = (unsigned long int)nNodes.elem(DIM - ii);
    }

    /*------------------------------------------------------------------------*/
    /* Read data. */
    double* ptrToMem = 0;
    double value = 0.0;
    // Get real number and store it in mNcMemory.
    double tmpVal = 0.0;

    for (std::size_t ii = 0; ii < this->nameVariables().size(); ++ii)
    {
        // mIdVariable[0] means real_part!
#ifdef SCAFES_HAVE_MPI
    #ifdef SCAFES_HAVE_NETCDF_PAR
        status = nc_var_par_access(this->mIdFile, this->mIdVariable[ii][0],
                                   NC_COLLECTIVE);
        wrapNCCall(status, __FILE__, __LINE__);
    #endif
#endif
        status = nc_get_vara_double(this->mIdFile, this->mIdVariable[ii][0],
                                    elemFirst, nElems, this->mNcMemory.data());
        wrapNCCall(status, __FILE__, __LINE__);

        ptrToMem = this->mNcMemory.data();
        //#pragma omp parallel for
        for (typename GridSub<DIM>::iterator
                 it = this->memNormal().at(ii).begin(),
                 et = this->memNormal().at(ii).end();
             it < et; ++it)
        {
            tmpVal = *ptrToMem;
            value = tmpVal;
            // Will be written to elemData.
            this->mElemData[ii][it.idxScalarNode()] = static_cast<TT>(value);
            ++ptrToMem;
        }
    }
    status = nc_close(this->mIdFile);
    wrapNCCall(status, __FILE__, __LINE__);
    this->mIsOpen = false;
}
/*----------------------------------------------------------------------------*/
template <typename TT, std::size_t DIM>
void DataFile<TT, DIM>::wrapNCCall(int const& status, std::string const& file,
                                   int const& line)
{
    // nothing to do: return
    if (status == NC_NOERR)
    {
        return;
    }
    // print useful error message
    std::cout << file << ":" << line << "\t NetCDF status code = " << status
              << ", \"";
    if (mNCErrorStrings.count(status) == 0)
    {
        std::cout << "Unknown error.";
    }
    else
    {
        std::cout << mNCErrorStrings[status];
    }
    std::cout << "\". Exiting." << std::endl;
    // can throw integers, not sure if useful though
    throw status;
}
/*----------------------------------------------------------------------------*/
// fill mNCErrorStrings
// TODO: push this into libscafes.so
template <typename TT, std::size_t DIM>
std::unordered_map<decltype(NC_NOERR), std::string>
DataFile<TT, DIM>::mNCErrorStrings;

} // End of namespace. //

#endif
