/* ScaFES
 * Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file ScaFES_Parameters.hpp
 *  @brief Contains the class Parameters.
 */

#ifndef SCAFES_PARAMETERS_HPP_
#define SCAFES_PARAMETERS_HPP_

#include "ScaFES_Config.hpp"
#include "ScaFES_Communicator.hpp"
#include "ScaFES_Environment.hpp"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>

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
// #include <boost/serialization/vector.hpp>
// #include <boost/serialization/version.hpp>
// #endif


#ifdef SCAFES_HAVE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
#endif

#ifdef SCAFES_HAVE_BOOST
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#endif

#ifdef SCAFES_HAVE_BOOST_REGEX
#include <boost/regex.hpp>
#endif

#include "ScaFES_Ntuple.hpp"

namespace ScaFES
{

#ifdef SCAFES_HAVE_BOOST_PROGRAM_OPTIONS
namespace po = boost::program_options;
#endif

/*******************************************************************************
 ******************************************************************************/
/**
 * \class Parameters
 * @brief The class \c Parameters contains all parameters related to ScaFES.
 *
 * Parameters can be given via the command line.
 * If a command line parameter is invalid, a help message will be displayed.
 *
 * The MPI environment lives as long as an instance of this class is living.
 * In particular, if the destructor of this class is called,
 * the MPI environment will be finalized. Thus, if the descructor is called
 * too early, the instance of the calls goes out of scope. So, please be careful
 * with calls to the destructor of this class.
 *
 * \remarks The master MPI process reads in the command line arguments
 * and distributed the parameters to all other MPI processes.
 */
class Parameters
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
// //     /** Creates own default constructor. */
// //     Parameters();

    /** Creates own constructor from command line parameters.
     * Parses the command line of a program call and provides a help
     * message.
     */
    Parameters(int argc, char* argv[]);

    /** Creates own copy constructor. */
    Parameters(const Parameters& rhs);

    /** Creates own destructor. */
    ~Parameters();

    /** Creates own assignment operator. */
    Parameters& operator=(Parameters rhs);

    /*----------------------------------------------------------------------
    | GETTER METHODS.
     ----------------------------------------------------------------------*/
    /** Returns the space dimension. */
    const int& dim() const;

    /** Returns the number of grid nodes in each dimension. */
    const std::vector<int>& nNodes() const;

    /** Returns the coordinates of the first grid point. */
    const std::vector<double>& coordNodeFirst() const;

    /** Returns the coordinates of the last grid point. */
    const std::vector<double>& coordNodeLast() const;

    /** Returns the number of partitions in each direction. */
    const std::vector<int>& nPartitions() const;

    /** Returns the directions in which the grid should be partitioned. */
    const std::vector<bool>& divideGrid() const;

    /** Returns the decision if the kind field should be read in. */
    const bool& readKindFile() const;

    /** Returns if the kind field should be written to file. */
    const bool& writeKindFile() const;

    /** Returns the name of the kind file. */
    const std::string& kindFile() const;

    /** Returns the name of the domain decomposition approach. */
    const std::string& typeDomainDecomp() const;

    /** Returns the name of the partition file for writing. */
    const std::string& nameWritePartitionFile() const;

    /** Returns the decision if a partition file should be written. */
    const bool& writePartitionFile() const;

    /** Returns the name of the partition file for reading. */
    const std::string& nameReadPartitionFile() const;

    /** Returns the decision if a partition file should be written. */
    const bool& readPartitionFile() const;

    /** Returns the start of the time interval. */
    const double& timeIntervalStart() const;

    /** Returns the end of the time interval. */
    const double& timeIntervalEnd() const;

    /** Returns the number of time steps. */
    const int& nTimesteps() const;

    /** Returns time step width. */
    const double& tau() const;

    /** Returns the number of snapshots. */
    const int& nSnapshots() const;

    /** Returns the number of layers at the border. */
    const int& nLayersAtBorder() const;

    /** Returns the decision if an output file should be written. */
    const bool& writeDataFile() const;

    /** Returns the name of the output file. */
    const std::string& nameDataFile() const;

    /** Returns the name of the config file. */
    const std::string& nameConfigFile() const;

    /** Returns if ADOL-C is enabled. */
    const bool& enabledAdolc() const;

    /** Returns if asynchronous MPI communication is enabled. */
    const bool& asynchronMode() const;

    /** Returns if the Boost.MPI skeleton concept for the MPI communication
     * will be used. */
    const bool& useBoostMpiSkeletonConcept() const;

    /** Returns if gradients should be computed. */
    const bool& computeGradients() const;

    /** Returns if gradients should be checked. */
    const bool& checkGradients() const;

    /** Returns the MPI rank for output. */
    int rankOutput() const;

    /** Returns the MPI communicator. */
    const ScaFES::Communicator& myWorld() const;

    /** Returns the environment. */
    const ScaFES::Environment& myEnv() const;

    /** Returns the MPI rank. */
    int rank() const;

    /** Returns the maximal output level at the console. */
    int outputLevelMax() const;

    /** Returns the current output level at the console. */
    int outputLevelCurr() const;

    /** Returns the maximum level of console output. */
    const int& indentDepth() const;

    /*----------------------------------------------------------------------
    | SETTER METHODS.
     ----------------------------------------------------------------------*/
    /** Increase output level by one. */
    void increaseLevel();

    /** Decrease output level by one. */
    void decreaseLevel();

    /*----------------------------------------------------------------------
    | COMPARISON METHODS.
    ----------------------------------------------------------------------*/
    /** Tests if two grids are equal. */
    bool operator==(const Parameters& rhs) const;

    /** Tests if two grids are not equal. */
    bool operator!=(const Parameters& rhs) const;

    /*----------------------------------------------------------------------
    | WORK METHODS.
    ----------------------------------------------------------------------*/
    /** Sets a prefix to the output of process with MPI rank. */
    std::string getPrefix(const std::string& nameFramework =
                                 std::string("ScaFES ")) const;

    /*----------------------------------------------------------------------
    | FREE METHODS WHICH ARE FRIENDS OF THIS CLASS.
    ----------------------------------------------------------------------*/
    /** Swaps to elements of this class. */
    friend void swap(Parameters& first, Parameters& second);

private:
    /*----------------------------------------------------------------------
    | WORK METHODS.
    ----------------------------------------------------------------------*/
#ifdef SCAFES_HAVE_BOOST
    /** Checks a parameter given by its name. */
    template <typename TT>
    void checkParameter(const po::variables_map& vm,
                        const std::string& nameParam, std::vector<TT>& param,
                        bool isNecessary);

    /** Checks a parameter with floating point values given by its name. */
    template <typename TT>
    void checkParameterDouble(const po::variables_map& vm,
                        const std::string& nameParam, std::vector<TT>& param,
                        bool isNecessary);
#endif

#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
    /** Serializes class. */
    template <class Archive>
    void serialize(Archive& ar, const unsigned int version);
#endif

    /*----------------------------------------------------------------------
    | MEMBER VARIABLES.
    ----------------------------------------------------------------------*/
    /** Space dimension. */
    int mDim;

    /** Number of grid nodes in each dimension. */
    std::vector<int> mNnodes;

    /** Coordinate of first grid point (left lower point of grid). */
    std::vector<double> mCoordNodeFirst;

    /** Coordinate of last grid point (right upper point of grid). */
    std::vector<double> mCoordNodeLast;

    /** Number of partitions in each direction. */
    std::vector<int> mNpartitions;

    /** Divide grid into partitions: Which directions are enabled? */
    std::vector<bool> mDivideGrid;

    /** Read kind field from a given kind file? */
    bool mReadKindFile;

    /** Write kind field? */
    bool mWriteKindFile;

    /** Name of the kind file. */
    std::string mKindFile;

    /** Type of domain decomposition. */
    std::string mTypeDomainDecomp;

    /** Write partition table to a partition file? */
    bool mWritePartitionFile;

    /** Name of partition file which should be written. */
    std::string mNameWritePartitionFile;

    /** Read partition table from a partition file? */
    bool mReadPartitionFile;

    /** Name of the partition file which should be read in. */
    std::string mNameReadPartitionFile;

    /** Start of time interval [tstart;tend]. */
    double mTimeIntervalStart;

    /** End of time interval [tstart;t_end]. */
    double mTimeIntervalEnd;

    /** Number of timesteps. */
    int mNtimesteps;

    /** Time step width. */
    double mTau;

    /** Number of snapshots: Write data each m timesteps. */
    int mNsnapshots;

    /** Number of layers at the border. */
    int mNlayersAtBorder;

    /** Write output file? */
    bool mWriteDataFile;

    /** Name of the output file. */
    std::string mNameOfDataFile;

    /** Name of the configuration file. */
    std::string mNameOfConfigFile;

    /** Is ADOL-C enabled? */
    bool mEnabledAdolc;

    /** Is asynchronous MPI communication enabled? */
    bool mAsynchronMode;

    /** Use Boost.MPI skeleton concept for the MPI communication? */
    bool mUseBoostMpiSkeletonConcept;

    /** Should gradients be computed? */
    bool mComputeGradients;

    /** Should gradients be checked? */
    bool mCheckGradients;

    /** MPI rank for output. */
    int mRankOutput;

    /** MPI environment. */
    ScaFES::Environment mMyEnv;

    /** MPI communicator. */
    ScaFES::Communicator mMyWorld;

    /** Maximal indentation depth level of console output. */
    int gOutputLevelMax;

    /** Current indentation depth level of console output. */
    int gOutputLevelCurr;

}; // End of class. //

/*******************************************************************************
 * LIFE CYCLE METHODS.
 ******************************************************************************/
inline Parameters::Parameters(int argc, char* argv[])
: mDim(3)
, mNnodes(3)
, mCoordNodeFirst(3)
, mCoordNodeLast(3)
, mNpartitions(3)
, mDivideGrid(3)
, mReadKindFile(false)
, mWriteKindFile(false)
, mKindFile("EMPTY")
, mTypeDomainDecomp("RCB")
, mWritePartitionFile(false)
, mNameWritePartitionFile("EMPTY")
, mReadPartitionFile(false)
, mNameReadPartitionFile("EMPTY")
, mTimeIntervalStart(0.0)
, mTimeIntervalEnd(1.0)
, mNtimesteps(1)
, mTau((mTimeIntervalEnd-mTimeIntervalStart)/mNtimesteps)
, mNsnapshots(1)
, mNlayersAtBorder(1)
, mWriteDataFile(false)
, mNameOfDataFile("EMPTYDATAFILE")
, mNameOfConfigFile("EMPTYCONFIGFILE")
, mEnabledAdolc(false)
, mAsynchronMode(true)
, mUseBoostMpiSkeletonConcept(true)
, mComputeGradients(false)
, mCheckGradients(false)
, mRankOutput(0)
, mMyEnv(argc, argv)
, mMyWorld(argc, argv)
, gOutputLevelMax(2)
, gOutputLevelCurr(gOutputLevelMax)
{
//     if (0 == this->myWorld().rank())
//     {
#ifdef SCAFES_HAVE_BOOST
        po::variables_map vm;

        /*----------------------------------------------------------------------
        | Describe the options of the command line.
        ----------------------------------------------------------------------*/
        po::options_description desc("Available ScaFES parameters");
        desc.add_options()("help", "Displays this help message.")(
            "dim", po::value<int>(), "Sets the space dimension.")(
            "nPartitions", po::value<std::string>(),
            "Sets the number of partitions in each dimension (Np1x..xNpd).")(
            "typeDomainDecomposition", po::value<std::string>(),
            "Sets the type of domain decomposition: RCB (RCB algorithm), UNI "
            "(Uniform decomposition in each direction).")(
            "nNodes", po::value<std::string>(),
            "Sets the number of grid nodes in each dimension (N1x..xNd).")(
            "coordNodeFirst", po::value<std::string>(),
            "Sets the coordinates of the first grid node (c1x..xcd).")(
            "coordNodeLast", po::value<std::string>(),
            "Sets the coordinates of the last grid node (c1x..xcd).")(
            "divideGrid", po::value<std::string>(),
            "Divides the grid in which dimension? (g1x..xgd).")(
            "readKindFile", po::value<std::string>(),
            "Reads the kind field from a given kind file.")(
            "writeKindFile", po::value<int>(),
            "Sets if the kind field should be written to a NetCDF data file, yes=1, no=0 (default=0).")(
            "readPartitionfile", po::value<std::string>(),
            "Reads the partition table from a given partition file.")(
            "writePartitionfile", po::value<std::string>(),
            "Writes the partition table to a given partition file.")(
            "starttime", po::value<double>(), "Sets the start time.")(
            "endtime", po::value<double>(), "Sets the end time.")(
            "nTimesteps", po::value<int>(),
            "Sets the number of time steps (=#{time intervals} (default=0).")(
            "tau", po::value<double>(), "Sets the time step width.")(
            "nSnapshots", po::value<int>(), "Sets how many times the data "
                                            "fields are written to output "
                                            "(default=1).")(
            "nLayersAtBorder", po::value<int>(),
            "Sets of how many layers the border should consist (default=1).")(
            "outputfile", po::value<std::string>(),
            "Sets the name of the output data file.")(
            "configfile", po::value<std::string>(),
            "Sets the name of the config file.")(
            "enabledAdolc", po::value<int>(),
            "Sets if ADOL-C should be enabled, 1=yes, 0=no (default=0).")(
            "asynchronMode", po::value<int>(), "Sets if asynchronous MPI "
                                               "communication should be used, 1=yes, 0=no "
                                               "(default=1).")(
            "useBoostMpiSkeletonConcept", po::value<int>(), "Sets if the Boost.MPI skeleton "
                                               "concept should be used for the MPI communication, "
                                               "1=yes, 0=no "
                                               "(default=1).")(
            "computeGradients", po::value<int>(),
            "Sets if gradients should be computed, 1=yes, 0=no (default=0).")(
            "checkGradients", po::value<int>(), "Sets if gradients of forward "
                                                "computation should be checked, "
                                                "1=yes, 0=no "
                                                "(default=0).")(
            "rankOutput", po::value<int>(),
            "Sets the rank of the MPI process whose output should be displayed "
            "(default=0).")("outputlevel", po::value<int>(),
                            "Controls the amount of console output.");

        /*----------------------------------------------------------------------
        | Try to parse the command line.
        ----------------------------------------------------------------------*/
        try
        {
            po::parsed_options parsed =
                po::command_line_parser(argc, argv).options(desc).allow_unregistered().run();
            po::store(parsed, vm);
            po::notify(vm);
        }
        catch (std::exception& e)
        {
            throw std::runtime_error("Wrong command line options.");
        }

        /*----------------------------------------------------------------------
        | Actions of command line options.
        ----------------------------------------------------------------------*/
        if (vm.count("help"))
        {
            std::cout << desc << "\n";
            this->mMyEnv.abort(this->mMyWorld,0);
        }
        else {

        /*----------------------------------------------------------------------
        | NECESSARY PARAMETERS.
        ----------------------------------------------------------------------*/
        if (vm.count("dim"))
        {
            this->mDim = vm["dim"].as<int>();

            if (0 >= this->dim())
            {
                std::cerr << "\nERROR: Space dimension is less or equal"
                             " to zero. Exit...\n\n";
                this->mMyEnv.abort(this->mMyWorld,-1);
            }
        }
        else
        {
            std::cerr << "\nERROR: Use --dim=<integral value>.\n\n";
            this->mMyEnv.abort(this->mMyWorld,-1);
        }

        /*--------------------------------------------------------------------*/
        this->mNnodes.resize(this->dim());
        this->checkParameter(vm, "nNodes", mNnodes, true);

        /*--------------------------------------------------------------------*/
        if (vm.count("nTimesteps"))
        {
            this->mNtimesteps = vm["nTimesteps"].as<int>();

            if (0 > this->nTimesteps())
            {
                std::cerr << "\nERROR: Number of time steps is negative."
                             " Exit...\n\n";
                this->mMyEnv.abort(this->mMyWorld,-1);
            }
        }
        else
        {
            std::cerr << "\nREMARK: Use --nTimesteps=<integral value>.\n\n";
        }

        /*--------------------------------------------------------------------*/
        this->mCoordNodeLast.resize(this->dim());
        this->checkParameterDouble(vm, "coordNodeLast", this->mCoordNodeLast, true);

        /*----------------------------------------------------------------------
        | OPTIONAL PARAMETERS.
        ----------------------------------------------------------------------*/
        this->mCoordNodeFirst.resize(this->dim());
        this->checkParameterDouble(vm, "coordNodeFirst", this->mCoordNodeFirst,
                             false);

        /*--------------------------------------------------------------------*/
        if (vm.count("typeDomainDecomposition"))
        {
            this->mTypeDomainDecomp =
                vm["typeDomainDecomposition"].as<std::string>();
        }
        else
        {
            std::cerr
                << "\nREMARK: Set default --typeDomainDecomposition=RCB.\n\n";
        }

        /*--------------------------------------------------------------------*/
        this->mNpartitions.resize(this->dim());
        if (vm.count("nPartitions"))
        {
            this->checkParameter(vm, "nPartitions", this->mNpartitions, false);
            int nPartitionsTotal = 1;
            for (std::size_t ii = 0; ii < this->nPartitions().size(); ++ii)
            {
                nPartitionsTotal *= (this->mNpartitions.at(ii));
            }
            if (nPartitionsTotal != this->mMyWorld.size())
            {
                std::cerr << "\nERROR: #{Partitions} != #{MPI processes}."
                             " Exit...\n\n";
                this->mMyEnv.abort(this->mMyWorld,-1);
            }
        }
        else
        {
            std::cerr << "\nREMARK: Set nPartitions = [nMPIprocs 1,...,1].\n\n";
            /* If the number of partitions is not given at the command line,
            * set the number from the number of generated MPI processes:
            * nPartitions = [MPI_size(), 1, ..., 1]. */
            this->mNpartitions.at(0) = this->mMyWorld.size();
            for (std::size_t ii = 1; ii < this->nPartitions().size(); ++ii)
            {
                this->mNpartitions.at(ii) = 1;
            }
        }

        /*--------------------------------------------------------------------*/
        this->mDivideGrid.resize(this->dim());
        this->checkParameter(vm, "divideGrid", this->mDivideGrid, false);

        /*--------------------------------------------------------------------*/
        if (vm.count("readKindFile"))
        {
            this->mKindFile = vm["readKindFile"].as<std::string>();
            this->mReadKindFile = true;
        }
        else
        {
            std::cerr << "\nREMARK: Use --readKindFile=<name>.\n\n";
        }
        /*--------------------------------------------------------------------*/
        if (vm.count("writeKindFile"))
        {
            int tmpWriteKindFile = vm["writeKindFile"].as<int>();
            if (0 == tmpWriteKindFile)
            {
                this->mWriteKindFile = false;
            }
            else
            {
                this->mWriteKindFile = true;
            }
        }
        else
        {
            std::cerr << "\nREMARK: Use --writeKindFile=<0/1>.\n\n";
        }
        /*--------------------------------------------------------------------*/
        if (vm.count("useBoostMpiSkeletonConcept"))
        {
            int tmpUseBoostMpiSkeletonConcept = vm["useBoostMpiSkeletonConcept"].as<int>();
            if (0 == tmpUseBoostMpiSkeletonConcept)
            {
                this->mUseBoostMpiSkeletonConcept = false;
            }
            else
            {
                this->mUseBoostMpiSkeletonConcept = true;
            }
        }
        else
        {
            std::cerr << "\nREMARK: Use --useBoostMpiSkeletonConcept=<0/1>.\n\n";
        }

        /*--------------------------------------------------------------------*/
        if (vm.count("readPartitionFile"))
        {
            this->mNameReadPartitionFile =
                vm["readPartitionfile"].as<std::string>();
            this->mReadPartitionFile = true;
        }
        else
        {
            std::cerr << "\nREMARK: Use --readPartitionfile=<name>.\n\n";
        }

        /*--------------------------------------------------------------------*/
        if (vm.count("writePartitionfile"))
        {
            this->mNameWritePartitionFile =
                vm["writePartitionfile"].as<std::string>();
            this->mWritePartitionFile = true;
        }
        else
        {
            std::cerr << "\nREMARK: Use --writePartitionfile=<name>.\n\n";
        }

        /*--------------------------------------------------------------------*/
        if (vm.count("starttime"))
        {
            this->mTimeIntervalStart = vm["starttime"].as<double>();

            if (0.0 > this->timeIntervalStart())
            {
                std::cerr << "\nERROR: Start time < 0.0."
                             " Exit...\n\n";
                this->mMyEnv.abort(this->mMyWorld,-1);
            }
        }
        else
        {
            std::cerr << "\nREMARK: Use --starttime=<real value>.\n\n";
        }

        /*--------------------------------------------------------------------*/
        if (vm.count("endtime"))
        {
            this->mTimeIntervalEnd = vm["endtime"].as<double>();

            if (0.0 >= this->timeIntervalEnd())
            {
                std::cerr << "\nERROR: End time <= 0.0."
                             " Exit...\n\n";
                this->mMyEnv.abort(this->mMyWorld,-1);
            }
        }
        else
        {
            std::cerr << "\nREMARK: Use --endtime=<real value>.\n\n";
        }

        /*--------------------------------------------------------------------*/
        if (vm.count("tau"))
        {
            this->mTau = vm["tau"].as<double>();

            if (0.0 >= this->tau())
            {
                std::cerr << "\nERROR: Time step width tau <= 0.0.   "
                             " Exit...\n\n";
                this->mMyEnv.abort(this->mMyWorld,-1);
            }
        }
        else if (vm.count("endtime") && vm.count("nTimesteps"))
        {
            if (0 == this->nTimesteps())
            {
                this->mTau = 0.0;
            }
            else
            {
                this->mTau =
                    (this->timeIntervalEnd() - this->timeIntervalStart()) /
                    this->nTimesteps();
            }
        }
        else
        {
            this->mTau = 1.0;
            std::cout << "REMARK: No time step width tau was specified.\n"
                         " Set tau=1.0 and proceeded. \n\n";
            std::cerr << "\nREMARK: Use --tau=<real value>.\n\n";
        }

        if (vm.count("tau") && vm.count("endTime") && vm.count("nTimesteps"))
        {
            if (this->mTau !=
                (this->timeIntervalEnd() - this->timeIntervalStart()) /
                    this->nTimesteps())
            {
                std::cerr
                    << "\nERROR: Inconsistency between time step width, "
                       "number of time steps and time interval. Exit...\n\n";
                exit(-1);
            }
        }

        /*--------------------------------------------------------------------*/
        if (vm.count("nSnapshots"))
        {
            this->mNsnapshots = vm["nSnapshots"].as<int>();

            if (0 >= nSnapshots())
            {
                std::cerr << "\nERROR: Number of snap shots is less or equal"
                             " to zero. It has to be greater or equal to one."
                             " Exit...\n\n";
                this->mMyEnv.abort(this->mMyWorld,-1);
            }
        }
        else
        {
            std::cerr << "\nREMARK: Use --nSnapshots=<integral value>.\n\n";
        }
        /*--------------------------------------------------------------------*/
        if (vm.count("nLayersAtBorder"))
        {
            this->mNlayersAtBorder = vm["nLayersAtBorder"].as<int>();

            for (std::size_t ii = 0; ii < this->nNodes().size(); ++ii)
            {
                if (this->nLayersAtBorder() >= (this->nNodes()[ii] / 2))
                {
                    std::cerr << "\nERROR: Number of layers at global border"
                                 " is too large (nLayers >= nNodes(ii)/2)."
                                 " Exit...\n\n";
                    this->mMyEnv.abort(this->mMyWorld,-1);
                }
            }
            if (0 >= this->nLayersAtBorder())
            {
                std::cerr << "\nERROR: Number of layers is less or equal"
                             " to zero. It has to be greater or equal to one."
                             " Exit...\n\n";
                this->mMyEnv.abort(this->mMyWorld,-1);
            }
        }
        else
        {
            std::cerr
                << "\nREMARK: Use --nLayersAtBorder=<integral value>.\n\n";
        }

        /*--------------------------------------------------------------------*/
        if (vm.count("outputfile"))
        {
            this->mNameOfDataFile = vm["outputfile"].as<std::string>();
            this->mWriteDataFile = true;
        }
        else
        {
            std::cerr << "\nREMARK: Use --outputfile=<name>.\n\n";
        }

        /*--------------------------------------------------------------------*/
        if (vm.count("configfile"))
        {
            this->mNameOfConfigFile = vm["configfile"].as<std::string>();
        }
        else
        {
            std::cerr << "\nREMARK: Use --configfile=<name>.\n\n";
        }
        /*--------------------------------------------------------------------*/
        if (vm.count("enabledAdolc"))
        {
            int tmpEnabledAdolc = vm["enabledAdolc"].as<int>();
            if (0 == tmpEnabledAdolc)
            {
                this->mEnabledAdolc = false;
            }
            else
            {
                this->mEnabledAdolc = true;
#ifndef SCAFES_HAVE_ADOLC
                // Enforce disabling of ADOL-C mode if software was not
                // built with ADOL-C.
                std::cerr << "\nWARNING: Framework was NOT built with ADOL-C "
                             "support!\n";
                std::cerr << "           ===> Enforced enabledAdolc=0!\n\n";

                this->mEnabledAdolc = false;
#endif
            }
        }
        else
        {
            std::cerr << "\nREMARK: Use --enabledAdolc=<0/1>.\n\n";
        }
        /*--------------------------------------------------------------------*/
        if (vm.count("asynchronMode"))
        {
            int tmpAsynchronMode = vm["asynchronMode"].as<int>();
            if (0 == tmpAsynchronMode)
            {
                this->mAsynchronMode = false;
            }
            else
            {
                this->mAsynchronMode = true;
            }
        }
        else
        {
            std::cerr << "\nREMARK: Use --asynchronMode=<0/1>.\n\n";
        }
        /*--------------------------------------------------------------------*/
        if (vm.count("computeGradients"))
        {
            int tmpComputeGradients = vm["computeGradients"].as<int>();
            if (0 == tmpComputeGradients)
            {
                this->mComputeGradients = false;
            }
            else
            {
                this->mComputeGradients = true;
            }
        }
        else
        {
            std::cerr << "\nREMARK: Use --computeGradients=<0/1>.\n\n";
        }
        /*--------------------------------------------------------------------*/
        if (vm.count("checkGradients"))
        {
            int tmpCheckGradients = vm["checkGradients"].as<int>();
            if (0 == tmpCheckGradients)
            {
                this->mCheckGradients = false;
            }
            else
            {
                this->mCheckGradients = true;
            }
        }
        else
        {
            std::cerr << "\nREMARK: Use --checkGradients=<0/1>.\n\n";
        }
        /*--------------------------------------------------------------------*/
        if (vm.count("rankOutput"))
        {
            this->mRankOutput = vm["rankOutput"].as<int>();
            if (static_cast<int>(this->rankOutput()) >
                this->myWorld().size() - 1)
            {
                throw std::runtime_error("rankOutput > rankMax.");
            }
        }
        else
        {
            std::cerr << "\nREMARK: Use --rankOutput=<integral value>.\n\n";
        }

        /*--------------------------------------------------------------------*/
        if (vm.count("outputlevel"))
        {
            gOutputLevelMax = vm["outputlevel"].as<int>();
            gOutputLevelCurr = gOutputLevelMax;
        }
        else
        {
            std::cerr << "\nREMARK: Use --outputlevel=<integral value>.\n\n";
        }
        } // End of if block related to "help".
#else
        for (int ii = 0; ii < this->dim(); ++ii)
        {
            this->mNnodes.at(ii) = 32;
            this->mCoordNodeFirst.at(ii) = 0.0;
            this->mCoordNodeLast.at(ii) = 1.0;
            this->mDivideGrid.at(ii) = false;
            this->mNpartitions.at(ii) = 1;
        }
        this->mDivideGrid.at(0) = true;
        this->mNpartitions.at(0) = this->mMyWorld.size();
#endif
//     }
//
//     // Communicate result.
//     // REMARK: Broadcasting all members using one call to broadcast() does
//     // not work at Juqueen!
//     // [KF, 04/22/2015]
//     //ScaFES::broadcast(this->mMyWorld, *this, 0);
//
//     ScaFES::broadcast(this->mMyWorld, this->mDim, 0);
//     ScaFES::broadcast(this->mMyWorld, this->mNnodes, 0);
//     ScaFES::broadcast(this->mMyWorld, this->mCoordNodeFirst, 0);
//     ScaFES::broadcast(this->mMyWorld, this->mCoordNodeLast, 0);
//     ScaFES::broadcast(this->mMyWorld, this->mNpartitions, 0);
//     ScaFES::broadcast(this->mMyWorld, this->mDivideGrid, 0);
//     ScaFES::broadcast(this->mMyWorld, this->mReadKindFile, 0);
//     ScaFES::broadcast(this->mMyWorld, this->mWriteKindFile, 0);
//     ScaFES::broadcast(this->mMyWorld, this->mTypeDomainDecomp, 0);
//     ScaFES::broadcast(this->mMyWorld, this->mWritePartitionFile, 0);
//     ScaFES::broadcast(this->mMyWorld, this->mNameWritePartitionFile, 0);
//     ScaFES::broadcast(this->mMyWorld, this->mReadPartitionFile, 0);
//     ScaFES::broadcast(this->mMyWorld, this->mNameReadPartitionFile, 0);
//     ScaFES::broadcast(this->mMyWorld, this->mTimeIntervalEnd, 0);
//     ScaFES::broadcast(this->mMyWorld, this->mNtimesteps, 0);
//     ScaFES::broadcast(this->mMyWorld, this->mTau, 0);
//     ScaFES::broadcast(this->mMyWorld, this->mNsnapshots, 0);
//     ScaFES::broadcast(this->mMyWorld, this->mNlayersAtBorder, 0);
//     ScaFES::broadcast(this->mMyWorld, this->mWriteDataFile, 0);
//     ScaFES::broadcast(this->mMyWorld, this->mNameOfDataFile, 0);
//     ScaFES::broadcast(this->mMyWorld, this->mNameOfConfigFile, 0);
//     ScaFES::broadcast(this->mMyWorld, this->mEnabledAdolc, 0);
//     ScaFES::broadcast(this->mMyWorld, this->mAsynchronMode, 0);
//     ScaFES::broadcast(this->mMyWorld, this->mUseBoostMpiSkeletonConcept, 0);
//     ScaFES::broadcast(this->mMyWorld, this->mComputeGradients, 0);
//     ScaFES::broadcast(this->mMyWorld, this->mCheckGradients, 0);
//     ScaFES::broadcast(this->mMyWorld, this->mRankOutput, 0);
//     ScaFES::broadcast(this->mMyWorld, this->mMyWorld, 0);
//     ScaFES::broadcast(this->mMyWorld, this->gOutputLevelMax, 0);
//     ScaFES::broadcast(this->mMyWorld, this->gOutputLevelCurr, 0);
}

/*----------------------------------------------------------------------------*/
#ifdef SCAFES_HAVE_BOOST_MPI
/** REMARK: Adding mMyEnv(rhs.myEnv()) results in the following compile error:
 * /usr/include/boost/mpi/environment.hpp:85:22: note: 'boost::mpi::environment::environment(const boost::mpi::environment&)' is implicitly deleted because the default definition would be ill-formed:
 class BOOST_MPI_DECL environment : noncopyable {
                      ^
/usr/include/boost/mpi/environment.hpp:85:22: error: use of deleted function 'boost::noncopyable_::noncopyable::noncopyable(const boost::noncopyable_::noncopyable&)'
*/
inline Parameters::Parameters(const Parameters& rhs)
: mDim(rhs.dim())
, mNnodes(rhs.nNodes())
, mCoordNodeFirst(rhs.coordNodeFirst())
, mCoordNodeLast(rhs.coordNodeLast())
, mNpartitions(rhs.nPartitions())
, mDivideGrid(rhs.divideGrid())
, mReadKindFile(rhs.readKindFile())
, mWriteKindFile(rhs.writeKindFile())
, mKindFile(rhs.kindFile())
, mTypeDomainDecomp(rhs.typeDomainDecomp())
, mWritePartitionFile(rhs.writePartitionFile())
, mNameWritePartitionFile(rhs.nameWritePartitionFile())
, mReadPartitionFile(rhs.readPartitionFile())
, mNameReadPartitionFile(rhs.nameReadPartitionFile())
, mTimeIntervalStart(rhs.timeIntervalStart())
, mTimeIntervalEnd(rhs.timeIntervalEnd())
, mNtimesteps(rhs.nTimesteps())
, mTau(rhs.tau())
, mNsnapshots(rhs.nSnapshots())
, mNlayersAtBorder(rhs.nLayersAtBorder())
, mWriteDataFile(rhs.writeDataFile())
, mNameOfDataFile(rhs.nameDataFile())
, mNameOfConfigFile(rhs.nameConfigFile())
, mEnabledAdolc(rhs.enabledAdolc())
, mAsynchronMode(rhs.asynchronMode())
, mUseBoostMpiSkeletonConcept(rhs.useBoostMpiSkeletonConcept())
, mComputeGradients(rhs.computeGradients())
, mCheckGradients(rhs.checkGradients())
, mRankOutput(rhs.rankOutput())
, mMyWorld(rhs.myWorld())
, gOutputLevelMax(rhs.outputLevelMax())
, gOutputLevelCurr(rhs.outputLevelCurr())
{
}
#else
inline Parameters::Parameters(const Parameters& rhs)
: mDim(rhs.dim())
, mNnodes(rhs.nNodes())
, mCoordNodeFirst(rhs.coordNodeFirst())
, mCoordNodeLast(rhs.coordNodeLast())
, mNpartitions(rhs.nPartitions())
, mDivideGrid(rhs.divideGrid())
, mReadKindFile(rhs.readKindFile())
, mWriteKindFile(rhs.writeKindFile())
, mKindFile(rhs.kindFile())
, mTypeDomainDecomp(rhs.typeDomainDecomp())
, mWritePartitionFile(rhs.writePartitionFile())
, mNameWritePartitionFile(rhs.nameWritePartitionFile())
, mReadPartitionFile(rhs.readPartitionFile())
, mNameReadPartitionFile(rhs.nameReadPartitionFile())
, mTimeIntervalStart(rhs.timeIntervalStart())
, mTimeIntervalEnd(rhs.timeIntervalEnd())
, mNtimesteps(rhs.nTimesteps())
, mTau(rhs.tau())
, mNsnapshots(rhs.nSnapshots())
, mNlayersAtBorder(rhs.nLayersAtBorder())
, mWriteDataFile(rhs.writeDataFile())
, mNameOfDataFile(rhs.nameDataFile())
, mNameOfConfigFile(rhs.nameConfigFile())
, mEnabledAdolc(rhs.enabledAdolc())
, mAsynchronMode(rhs.asynchronMode())
, mUseBoostMpiSkeletonConcept(rhs.useBoostMpiSkeletonConcept())
, mComputeGradients(rhs.computeGradients())
, mCheckGradients(rhs.checkGradients())
, mRankOutput(rhs.rankOutput())
, mMyEnv(rhs.myEnv())
, mMyWorld(rhs.myWorld())
, gOutputLevelMax(rhs.outputLevelMax())
, gOutputLevelCurr(rhs.outputLevelCurr())
{
}
#endif
/*----------------------------------------------------------------------------*/
inline Parameters& Parameters::operator=(Parameters rhs)
{
    swap(*this, rhs);
    return *this;
}
/*----------------------------------------------------------------------------*/
inline Parameters::~Parameters()
{
}

/*******************************************************************************
 * GETTER METHODS.
 ******************************************************************************/
inline const int& Parameters::dim() const
{
    return this->mDim;
}
/*----------------------------------------------------------------------------*/
inline const std::vector<int>& Parameters::nNodes() const
{
    return this->mNnodes;
}
/*----------------------------------------------------------------------------*/
inline const std::vector<double>& Parameters::coordNodeFirst() const
{
    return this->mCoordNodeFirst;
}
/*----------------------------------------------------------------------------*/
inline const std::vector<double>& Parameters::coordNodeLast() const
{
    return this->mCoordNodeLast;
}
/*----------------------------------------------------------------------------*/
inline const std::vector<int>& Parameters::nPartitions() const
{
    return this->mNpartitions;
}
/*----------------------------------------------------------------------------*/
inline const std::vector<bool>& Parameters::divideGrid() const
{
    return this->mDivideGrid;
}
/*----------------------------------------------------------------------------*/
inline const bool& Parameters::readKindFile() const
{
    return this->mReadKindFile;
}
/*----------------------------------------------------------------------------*/
inline const bool& Parameters::writeKindFile() const
{
    return this->mWriteKindFile;
}
/*----------------------------------------------------------------------------*/
inline const std::string& Parameters::kindFile() const
{
    return this->mKindFile;
}
/*----------------------------------------------------------------------------*/
inline const std::string& Parameters::typeDomainDecomp() const
{
    return this->mTypeDomainDecomp;
}
/*----------------------------------------------------------------------------*/
inline const bool& Parameters::writePartitionFile() const
{
    return this->mWritePartitionFile;
}
/*----------------------------------------------------------------------------*/
inline const std::string& Parameters::nameWritePartitionFile() const
{
    return this->mNameWritePartitionFile;
}
/*----------------------------------------------------------------------------*/
inline const bool& Parameters::readPartitionFile() const
{
    return this->mReadPartitionFile;
}
/*----------------------------------------------------------------------------*/
inline const std::string& Parameters::nameReadPartitionFile() const
{
    return this->mNameReadPartitionFile;
}
/*----------------------------------------------------------------------------*/
inline const double& Parameters::timeIntervalStart() const
{
    return this->mTimeIntervalStart;
}
/*----------------------------------------------------------------------------*/
inline const double& Parameters::timeIntervalEnd() const
{
    return this->mTimeIntervalEnd;
}
/*----------------------------------------------------------------------------*/
inline const int& Parameters::nTimesteps() const
{
    return this->mNtimesteps;
}
/*----------------------------------------------------------------------------*/
inline const double& Parameters::tau() const
{
    return this->mTau;
}
/*----------------------------------------------------------------------------*/
inline const int& Parameters::nSnapshots() const
{
    return this->mNsnapshots;
}
/*----------------------------------------------------------------------------*/
inline const int& Parameters::nLayersAtBorder() const
{
    return this->mNlayersAtBorder;
}
/*----------------------------------------------------------------------------*/
inline const bool& Parameters::writeDataFile() const
{
    return this->mWriteDataFile;
}
/*----------------------------------------------------------------------------*/
inline const std::string& Parameters::nameDataFile() const
{
    return this->mNameOfDataFile;
}
/*----------------------------------------------------------------------------*/
inline const std::string& Parameters::nameConfigFile() const
{
    return this->mNameOfConfigFile;
}
/*----------------------------------------------------------------------------*/
inline const bool& Parameters::enabledAdolc() const
{
    return this->mEnabledAdolc;
}
/*----------------------------------------------------------------------------*/
inline const bool& Parameters::asynchronMode() const
{
    return this->mAsynchronMode;
}
/*----------------------------------------------------------------------------*/
inline const bool& Parameters::useBoostMpiSkeletonConcept() const
{
    return this->mUseBoostMpiSkeletonConcept;
}
/*----------------------------------------------------------------------------*/
inline const bool& Parameters::computeGradients() const
{
    return this->mComputeGradients;
}
/*----------------------------------------------------------------------------*/
inline const bool& Parameters::checkGradients() const
{
    return this->mCheckGradients;
}
/*----------------------------------------------------------------------------*/
inline int Parameters::rankOutput() const
{
    return this->mRankOutput;
}
/*----------------------------------------------------------------------------*/
inline const ScaFES::Communicator& Parameters::myWorld() const
{
    return this->mMyWorld;
}
/*----------------------------------------------------------------------------*/
inline const ScaFES::Environment& Parameters::myEnv() const
{
    return this->mMyEnv;
}
/*----------------------------------------------------------------------------*/
inline int Parameters::rank() const
{
    return this->myWorld().rank();
}
/*----------------------------------------------------------------------------*/
inline int Parameters::outputLevelMax() const
{
    return this->gOutputLevelMax;
}
/*----------------------------------------------------------------------------*/
inline int Parameters::outputLevelCurr() const
{
    return this->gOutputLevelCurr;
}
/*----------------------------------------------------------------------------*/
inline const int& Parameters::indentDepth() const
{
    return this->gOutputLevelCurr;
}

/*******************************************************************************
 * SETTER METHODS.
 ******************************************************************************/
inline void Parameters::increaseLevel()
{
    ++(this->gOutputLevelCurr);
}
/*----------------------------------------------------------------------------*/
inline void Parameters::decreaseLevel()
{
    --(this->gOutputLevelCurr);
}

/*******************************************************************************
 * COMPARISON METHODS.
 ******************************************************************************/
inline bool Parameters::operator==(const Parameters& rhs) const
{
    bool isEqual = true;

    if (this->dim() != rhs.dim())
    {
        isEqual = false;
    }
    if (this->nNodes() != rhs.nNodes())
    {
        isEqual = false;
    }
    if (this->coordNodeFirst() != rhs.coordNodeFirst())
    {
        isEqual = false;
    }
    if (this->coordNodeLast() != rhs.coordNodeLast())
    {
        isEqual = false;
    }
    if (this->nPartitions() != rhs.nPartitions())
    {
        isEqual = false;
    }
    if (this->divideGrid() != rhs.divideGrid())
    {
        isEqual = false;
    }
    if (this->readKindFile() != rhs.readKindFile())
    {
        isEqual = false;
    }
    if (this->kindFile() != rhs.kindFile())
    {
        isEqual = false;
    }
    if (this->writePartitionFile() != rhs.writePartitionFile())
    {
        isEqual = false;
    }
    if (this->nameWritePartitionFile() != rhs.nameWritePartitionFile())
    {
        isEqual = false;
    }
    if (this->nameReadPartitionFile() != rhs.nameReadPartitionFile())
    {
        isEqual = false;
    }
    if (this->timeIntervalStart() != rhs.timeIntervalStart())
    {
        isEqual = false;
    }
    if (this->timeIntervalEnd() != rhs.timeIntervalEnd())
    {
        isEqual = false;
    }
    if (this->nTimesteps() != rhs.nTimesteps())
    {
        isEqual = false;
    }
    if (this->tau() != rhs.tau())
    {
        isEqual = false;
    }
    if (this->nSnapshots() != rhs.nSnapshots())
    {
        isEqual = false;
    }
    if (this->writeDataFile() != rhs.writeDataFile())
    {
        isEqual = false;
    }
    if (this->nameDataFile() != rhs.nameDataFile())
    {
        isEqual = false;
    }
    if (this->nameConfigFile() != rhs.nameConfigFile())
    {
        isEqual = false;
    }
    if (this->enabledAdolc() != rhs.enabledAdolc())
    {
        isEqual = false;
    }
    if (this->asynchronMode() != rhs.asynchronMode())
    {
        isEqual = false;
    }
    if (this->useBoostMpiSkeletonConcept() != rhs.useBoostMpiSkeletonConcept())
    {
        isEqual = false;
    }
    if (this->computeGradients() != rhs.computeGradients())
    {
        isEqual = false;
    }
    if (this->checkGradients() != rhs.checkGradients())
    {
        isEqual = false;
    }
    if (this->rankOutput() != rhs.rankOutput())
    {
        isEqual = false;
    }
//     if (this->myWorld() != rhs.myWorld())
//     {
//         isEqual = false;
//     }

    return isEqual;
}
/*----------------------------------------------------------------------------*/
inline bool Parameters::operator!=(const Parameters& rhs) const
{
    return !(*this == rhs);
}

/*******************************************************************************
 * WORK METHODS.
 ******************************************************************************/
#ifdef SCAFES_HAVE_BOOST
template <typename TT>
inline void Parameters::checkParameter(const po::variables_map& vm,
                                       const std::string& nameParam,
                                       std::vector<TT>& param, bool isNecessary)
{
    if (vm.count(nameParam))
    {
        std::ostringstream ostrRegex;

        for (int ii = 0; ii < dim(); ++ii)
        {
            ostrRegex << "(-?\\d+)";

            if ((dim() - 1) > ii)
            {
                ostrRegex << "x";
            }
        }

        std::string strRegex = ostrRegex.str();
        boost::regex re(strRegex);
        std::string strParam;
        boost::cmatch matches;

        strParam = vm[nameParam].as<std::string>();
        if (boost::regex_match(strParam.c_str(), matches, re))
        {
            for (std::size_t i = 1; i < matches.size(); ++i)
            {
                std::string match(matches[i].first, matches[i].second);
                param.at(i - 1) = atoi(match.c_str());
            }
        }
        else
        {
            std::ostringstream ostrErr;

            for (int ii = 0; ii < dim(); ++ii)
            {
                ostrErr << "n" << ii;

                if ((dim() - 1) > ii)
                {
                    ostrErr << "x";
                }
            }

            std::string strErr = ostrErr.str();
            std::cerr << "\nERROR: The parameter --" << nameParam << "="
                      << strParam << " is not of a format like"
                                     " --" << nameParam << "=" << strErr
                      << ". Exit... \n\n";
            this->mMyEnv.abort(this->mMyWorld,-1);
        }
    }
    else
    {
        std::cerr << "\nREMARK: Use --" << nameParam << "=<a1xa2xa3>\n\n";

        if (isNecessary)
        {
            this->mMyEnv.abort(this->mMyWorld,-1);
        }
    }
}
template <typename TT>
inline void Parameters::checkParameterDouble(const po::variables_map& vm,
                                       const std::string& nameParam,
                                       std::vector<TT>& param, bool isNecessary)
{
    if (vm.count(nameParam))
    {
        std::ostringstream ostrRegex;

        for (int ii = 0; ii < dim(); ++ii)
        {
            ostrRegex << "([-+]?[0-9]*\\.?[0-9]+)";

            if ((dim() - 1) > ii)
            {
                ostrRegex << "x";
            }
        }

        std::string strRegex = ostrRegex.str();
        boost::regex re(strRegex);
        std::string strParam;
        boost::cmatch matches;

        strParam = vm[nameParam].as<std::string>();
        if (boost::regex_match(strParam.c_str(), matches, re))
        {
            for (std::size_t i = 1; i < matches.size(); ++i)
            {
                std::string match(matches[i].first, matches[i].second);
                param.at(i - 1) = atof(match.c_str());
            }
        }
        else
        {
            std::ostringstream ostrErr;

            for (int ii = 0; ii < dim(); ++ii)
            {
                ostrErr << "n" << ii;

                if ((dim() - 1) > ii)
                {
                    ostrErr << "x";
                }
            }

            std::string strErr = ostrErr.str();
            std::cerr << "\nERROR: The parameter --" << nameParam << "="
                      << strParam << " is not of a format like"
                                     " --" << nameParam << "=" << strErr
                      << ". Exit... \n\n";
            this->mMyEnv.abort(this->mMyWorld,-1);
        }
    }
    else
    {
        std::cerr << "\nREMARK: Use --" << nameParam << "=<a1xa2xa3>\n\n";

        if (isNecessary)
        {
            this->mMyEnv.abort(this->mMyWorld,-1);
        }
    }
}
#endif

/*----------------------------------------------------------------------------*/
/** Free function which sets a prefix to the output of process with
 * identifier \c rank. */
inline std::string Parameters::getPrefix(const std::string& nameFramework) const
{
    std::string blanks("     ");
    std::ostringstream blankStream;
    blankStream << nameFramework << ": Rank=" << this->rankOutput() << ": ";
    int nIndents = this->gOutputLevelMax - this->gOutputLevelCurr;
    for (int jj = 0; jj < nIndents; ++jj)
    {
        blankStream << blanks;
    }
    return blankStream.str();
}

/*----------------------------------------------------------------------------*/
#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
template <class Archive>
void Parameters::serialize(Archive& ar, const unsigned int version)
{
    if (1 <= version)
    {
        ar&(this->mDim);
        ar&(this->mNnodes);
        ar&(this->mCoordNodeFirst);
        ar&(this->mCoordNodeLast);
        ar&(this->mNpartitions);
        ar&(this->mDivideGrid);
        ar&(this->mReadKindFile);
        ar&(this->mWriteKindFile);
        ar&(this->mKindFile);
        ar&(this->mTypeDomainDecomp);
        ar&(this->mWritePartitionFile);
        ar&(this->mNameWritePartitionFile);
        ar&(this->mReadPartitionFile);
        ar&(this->mNameReadPartitionFile);
        ar&(this->mTimeIntervalStart);
        ar&(this->mTimeIntervalEnd);
        ar&(this->mNtimesteps);
        ar&(this->mTau);
        ar&(this->mNsnapshots);
        ar&(this->mNlayersAtBorder);
        ar&(this->mWriteDataFile);
        ar&(this->mNameOfDataFile);
        ar&(this->mNameOfConfigFile);
        ar&(this->mEnabledAdolc);
        ar&(this->mAsynchronMode);
        ar&(this->mUseBoostMpiSkeletonConcept);
        ar&(this->mComputeGradients);
        ar&(this->mCheckGradients);
        ar&(this->mRankOutput);
        ar&(this->gOutputLevelMax);
        ar&(this->gOutputLevelCurr);
    }
}
#endif

/*******************************************************************************
 * FREE METHODS WHICH ARE FRIENDS OF THIS CLASS.
 ******************************************************************************/
inline void swap(Parameters& first, Parameters& second)
{
    std::swap(first.mDim, second.mDim);
    std::swap(first.mNnodes, second.mNnodes);
    std::swap(first.mCoordNodeFirst, second.mCoordNodeFirst);
    std::swap(first.mCoordNodeLast, second.mCoordNodeLast);
    std::swap(first.mNpartitions, second.mNpartitions);
    std::swap(first.mDivideGrid, second.mDivideGrid);
    std::swap(first.mReadKindFile, second.mReadKindFile);
    std::swap(first.mWriteKindFile, second.mWriteKindFile);
    std::swap(first.mKindFile, second.mKindFile);
    std::swap(first.mTypeDomainDecomp, second.mTypeDomainDecomp);
    std::swap(first.mWritePartitionFile, second.mWritePartitionFile);
    std::swap(first.mNameWritePartitionFile, second.mNameWritePartitionFile);
    std::swap(first.mReadPartitionFile, second.mReadPartitionFile);
    std::swap(first.mNameReadPartitionFile, second.mNameReadPartitionFile);
    std::swap(first.mTimeIntervalStart, second.mTimeIntervalStart);
    std::swap(first.mTimeIntervalEnd, second.mTimeIntervalEnd);
    std::swap(first.mNtimesteps, second.mNtimesteps);
    std::swap(first.mTau, second.mTau);
    std::swap(first.mNsnapshots, second.mNsnapshots);
    std::swap(first.mNlayersAtBorder, second.mNlayersAtBorder);
    std::swap(first.mWriteDataFile, second.mWriteDataFile);
    std::swap(first.mNameOfDataFile, second.mNameOfDataFile);
    std::swap(first.mNameOfConfigFile, second.mNameOfConfigFile);
    std::swap(first.mEnabledAdolc, second.mEnabledAdolc);
    std::swap(first.mAsynchronMode, second.mAsynchronMode);
    std::swap(first.mUseBoostMpiSkeletonConcept, second.mUseBoostMpiSkeletonConcept);
    std::swap(first.mComputeGradients, second.mComputeGradients);
    std::swap(first.mCheckGradients, second.mCheckGradients);
    std::swap(first.mRankOutput, second.mRankOutput);
    std::swap(first.mMyWorld, second.mMyWorld);
    std::swap(first.gOutputLevelMax, second.gOutputLevelMax);
    std::swap(first.gOutputLevelCurr, second.gOutputLevelCurr);
}

} // End of namespace. //

/*******************************************************************************
 ******************************************************************************/
#ifdef SCAFES_HAVE_BOOST_SERIALIZATION
/** Boost serialization version of class \c Parameter. */
BOOST_CLASS_VERSION(ScaFES::Parameters, 2);
#endif

#endif
