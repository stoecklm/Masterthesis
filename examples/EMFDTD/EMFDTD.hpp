/* ScaFES
 * Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file EMFDTD.hpp
 *
 *  @brief Implementation of EMFDTD.
 */

#ifndef EMFDTD_HPP_
#define EMFDTD_HPP_

#include <fstream>
#include <string>
#include <vector>

#include <boost/property_tree/ini_parser.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/serialization/vector.hpp>

#include <ScaFES.hpp>

#include "physical_constants.hpp"
#include "geometry/IndicatorEllipsoid.hpp"
#include "CPMLParams.hpp"
#include "stencils/CompleteCPMLCurlStencil.hpp"

#include "sources/ModulatedGaussianSource.hpp"
#include "conversion.hpp"

/*******************************************************************************
 ******************************************************************************/
/**
 * \class EMFDTDng
 *  @brief Class for discretized EM problem.
 *
 * \section EMFDTDng3D EMFDTD at ellipsoid.
 *
 * The EM problem is discretized using the FDTD method (Finite Differences
 * at the Time Domain).
 *
 * The workflow for running the simulation is as follows:
 *
 *      - Set up the computational domain within in the source code.
 *      - Recompile source code.
 *      - Set up EM parameters within the ini file.
 *      - Run simulation.
 *
 * Parameters of the PML (perfectly matched layers):
 *
 *      [PML]
        d=10
        m_a=1.
        m=3
        k_max=1.
        s_max=0.01222
        a_max=1.e-3

 * Location and parameters of the source.
 *
        [source]
        ; Location within computational domain
        location=25x25x25
        ; modulation frequency:
        f_mod = 0
        Gauss pulse: time step with maximal stimulation (=1)
        t0 = 40.0
        ; Width
        spread = 8.0

 * Location and target file of the sink.
 *
        [sink]
        ; Location within computational domain
        location=30x30x30
        ; Target file
        data_file = reference.txt

 * Domain initialization
 *
        [initialization]
        ; for ellipsoid: (px,py,pz,rx,ry,rz)
        parameters = 100 100 100 50 50 50
 */

/** data field access by symbolic name. */
enum class DFs: std::vector<std::string>::size_type {
    SIGMA=0,
    Ex, Ey, Ez,
    Hx, Hy, Hz,
    Dx, Dy, Dz,
    Ix, Iy, Iz,
    Psi_ExDy, Psi_ExDz, Psi_EyDx,
    Psi_EyDz, Psi_EzDx, Psi_EzDy,
    Psi_HxDy, Psi_HxDz, Psi_HyDx,
    Psi_HyDz, Psi_HzDx, Psi_HzDy
};

/** when to write which data field. */
typedef ScaFES::WriteHowOften W;
std::vector<W> initWritePolicy()
{
    std::vector<W> v;
    v.push_back(W::AT_START);
    v.push_back(W::LIKE_GIVEN_AT_CL); v.push_back(W::LIKE_GIVEN_AT_CL); v.push_back(W::LIKE_GIVEN_AT_CL);
    v.push_back(W::NEVER); v.push_back(W::NEVER); v.push_back(W::NEVER);
    v.push_back(W::NEVER); v.push_back(W::NEVER); v.push_back(W::NEVER);
    v.push_back(W::NEVER); v.push_back(W::NEVER); v.push_back(W::NEVER);
    v.push_back(W::NEVER); v.push_back(W::NEVER); v.push_back(W::NEVER);
    v.push_back(W::NEVER); v.push_back(W::NEVER); v.push_back(W::NEVER);
    v.push_back(W::NEVER); v.push_back(W::NEVER); v.push_back(W::NEVER);
    v.push_back(W::NEVER); v.push_back(W::NEVER); v.push_back(W::NEVER);
    return v;
}
static const std::vector<W> writePolicy = initWritePolicy();

/** data field names. */
std::vector<std::string> initDFNames()
{
    std::vector<std::string> v;
    v.push_back("SIGMA");
    v.push_back("Ex"); v.push_back("Ey"); v.push_back("Ez");
    v.push_back("Hx"); v.push_back("Hy"); v.push_back("Hz");
    v.push_back("Dx"); v.push_back("Dy"); v.push_back("Dz");
    v.push_back("Ix"); v.push_back("Iy"); v.push_back("Iz");
    v.push_back("Psi_ExDy"); v.push_back("Psi_ExDz"); v.push_back("Psi_EyDx");
    v.push_back("Psi_EYDz"); v.push_back("Psi_EzDx"); v.push_back("Psi_EzDy");
    v.push_back("Psi_HxDy"); v.push_back("Psi_HxDz"); v.push_back("Psi_HyDx");
    v.push_back("Psi_HyDz"); v.push_back("Psi_HzDx"); v.push_back("Psi_HzDy");
    return v;
}
static const std::vector<std::string> DFNames = initDFNames();

/** MPI halos. */
std::vector<int> initHaloWidths()
{
    std::vector<int> v;
    v.push_back(0);
    v.push_back(1); v.push_back(1); v.push_back(1);
    v.push_back(1); v.push_back(1); v.push_back(1);
    v.push_back(0); v.push_back(0); v.push_back(0);
    v.push_back(0); v.push_back(0); v.push_back(0);
    v.push_back(0); v.push_back(0); v.push_back(0);
    v.push_back(0); v.push_back(0); v.push_back(0);
    v.push_back(0); v.push_back(0); v.push_back(0);
    v.push_back(0); v.push_back(0); v.push_back(0);
    return v;
}
static const std::vector<int> HaloWidths = initHaloWidths();

template<typename CT, bool ForwardOnly>
class EMFDTDng: public ScaFES::Problem< EMFDTDng<CT, ForwardOnly>, CT, 3> {
public:
    static const size_t DIM=3;
private:
    typedef ScaFES::Problem< EMFDTDng<CT, ForwardOnly>, CT, 3> Super;
    typedef std::vector<std::string>::size_type size_type;
    typedef ScaFES::Ntuple<int, DIM> Index;
    typedef typename Index::value_type index_type;
    using PTree=boost::property_tree::ptree;

    // Stencils
    typedef CompleteCurlCPMLStencil<CT, Index, 0, 0, 1, 0, 1, 0> HxStencil;
    typedef CompleteCurlCPMLStencil<CT, Index, 1, 0, 0, 0, 0, 1> HyStencil;
    typedef CompleteCurlCPMLStencil<CT, Index, 0, 1, 0, 1, 0, 0> HzStencil;
    typedef CompleteCurlCPMLStencil<CT, Index, 0, -1, 0, 0, 0, -1> DxStencil;
    typedef CompleteCurlCPMLStencil<CT, Index, 0, 0, -1, -1, 0, 0> DyStencil;
    typedef CompleteCurlCPMLStencil<CT, Index, -1, 0, 0, 0, -1, 0> DzStencil;
public:

    EMFDTDng(
        ScaFES::Parameters const& cl,
        ScaFES::GridGlobal<DIM> const& gg,
        const PTree& ptree
        )
    : Super(cl, gg, true
        , DFNames // DF names
        , HaloWidths // MPI halo width
        , std::vector<bool>(DFNames.size(), false) // is known? --> false
        , std::vector<int>(DFNames.size(), 0) // border layers: 0 -> complete volume
        , std::vector<CT>(DFNames.size(), CT(0.)) // initial values for DFs
        , writePolicy
        , std::vector<bool>(DFNames.size(), false) // compute error? --> false
        , parse<CT>(ptree.get<std::string>("initialization.parameters"), " ") // params, vOld in initialize()
        )
    , ptree_(ptree)
    , cfln_(.5)
    , ds_(1.)
    , dt_((cfln_*ds_)/C_0)
    , cpmlp_(dt_,
        ptree_.get<std::size_t>("PML.d"),
        ptree_.get<CT>("PML.m_a"),
        ptree_.get<CT>("PML.m"),
        ptree_.get<CT>("PML.k_max"),
        ptree_.get<CT>("PML.s_max"),
        ptree_.get<CT>("PML.a_max")
        )
    , first_idx_(0, 0, 0)
    , last_idx_(cl.nNodes().at(0)-1, cl.nNodes().at(2)-1, cl.nNodes().at(2)-1)
    , loc_src_(vec2idx(parse<index_type>(ptree_.get<std::string>("source.location"))))
    , loc_sink_(vec2idx(parse<index_type>(ptree_.get<std::string>("sink.location"))))
    , reference_(load_reference(ptree_.get<std::string>("optimization.reference_data_file"), cl.nTimesteps()))
    , src_(dt_,ptree_)
    , hxstencil_(cpmlp_, Index(1, 0, 0), Index(0, -1, -1), last_idx_, cfln_)
    , hystencil_(cpmlp_, Index(0, 1, 0), Index(-1, 0, -1), last_idx_, cfln_)
    , hzstencil_(cpmlp_, Index(0, 0, 1), Index(-1, -1, 0), last_idx_, cfln_)
    , dxstencil_(cpmlp_, Index(0, 0, 1), Index(-1, 0, 0), last_idx_, -cfln_)
    , dystencil_(cpmlp_, Index(1, 0, 1), Index(0, -1, 0), last_idx_, -cfln_)
    , dzstencil_(cpmlp_, Index(1, 1, 0), Index(0, 0, -1), last_idx_, -cfln_)
    {
        boost::mpi::communicator world;

        if (world.rank()==0)
            std::cout<<"Source at "<<loc_src_
            <<"\nSink at "<<loc_sink_
            <<"\nInitialization parameters: "
            <<ptree.get<std::string>("initialization.parameters")
            <<std::endl;

        if (loc_src_==loc_sink_)
            throw std::runtime_error("Source and sink cannot be at the same position.");

        if ((loc_src_<first_idx_)||(loc_src_>last_idx_))
            throw std::runtime_error("Source not within simulation area.");
        if ((loc_sink_<first_idx_)||(loc_sink_>last_idx_))
            throw std::runtime_error("Sink not within simulation area.");

        if (cl.tau()!=1.0)
            throw std::runtime_error("Please choose \"--starttime\" and \"--endtime\" so that tau()=1.0.");
    }

    /**
     * Empty.
     */
    template<typename TT>
    void initInner(std::vector< ScaFES::DataField<TT, DIM> >& /*vNew*/,
        std::vector< ScaFES::DataField<TT, DIM> > const& /*vOld*/,
        Index const& /*idxNode*/,
        int const& /*timestep*/) { }

    /**
     * Empty.
     */
    template<typename TT>
    void initBorder(std::vector< ScaFES::DataField<TT, DIM> >& /*vNew*/,
            std::vector< ScaFES::DataField<TT, DIM> > const& /*vOld*/,
            Index const& /*idxNode*/,
            int const& /*timestep*/) { }

    /**
     * Initialize from p.
     */
    template <typename TT>
    void initInner
    (
     std::vector< ScaFES::DataField<TT, DIM> >& dfs,
     std::vector<TT> const& p,
     Index const& idx,
     int const& /*timestep*/
    )
    {
            typedef ScaFES::Ntuple<TT, DIM> TT3;
            TT3 fidx(TT(idx.elem(0)), TT(idx.elem(1)), TT(idx.elem(2)));

            const TT inside=3.546e7;
            const TT outside=0.;

            IndicatorEllipsoid<TT> e(p.at(0), p.at(1), p.at(2), p.at(3), p.at(4), p.at(5), inside, outside);

            const TT sigma=e(fidx);

            dfs[0](idx)=sigma;
    }

    /**
     * Wire to initInner.
     */
    template<typename TT>
    void initBorder(std::vector< ScaFES::DataField<TT, DIM> >& dfs,
            std::vector<TT> const& p,
            Index const& idx,
            int const& timestep)
    {
        this->initInner(dfs, p, idx, timestep);
    }

    /**
     * Compute H.
     */
    template<typename TT>
    void updateInner(std::vector< ScaFES::DataField<TT, DIM> >& n,
            std::vector< ScaFES::DataField<TT, DIM> > const& o,
            Index const& idx,
            int const& /*timestep*/)
    {
        // start by updating H from E
        n[sc(DFs::Hx)](idx)=hxstencil_(idx, o[sc(DFs::Hx)],
            o[sc(DFs::Ey)], o[sc(DFs::Psi_EyDz)], n[sc(DFs::Psi_EyDz)],
            o[sc(DFs::Ez)], o[sc(DFs::Psi_EzDy)], n[sc(DFs::Psi_EzDy)]);

        n[sc(DFs::Hy)](idx)=hystencil_(idx, o[sc(DFs::Hy)],
            o[sc(DFs::Ez)], o[sc(DFs::Psi_EzDx)], n[sc(DFs::Psi_EzDx)],
            o[sc(DFs::Ex)], o[sc(DFs::Psi_ExDz)], n[sc(DFs::Psi_ExDz)]);

        n[sc(DFs::Hz)](idx)=hzstencil_(idx, o[sc(DFs::Hz)],
            o[sc(DFs::Ex)], n[sc(DFs::Psi_ExDy)], n[sc(DFs::Psi_ExDy)],
            o[sc(DFs::Ey)], n[sc(DFs::Psi_EyDx)], n[sc(DFs::Psi_EyDx)]);
    }

    /**
     * Wire to updateInner.
     */
    template<typename TT>
    void updateBorder(std::vector< ScaFES::DataField<TT, DIM> >& vNew,
            std::vector< ScaFES::DataField<TT, DIM>>const& vOld,
            Index const& idx,
            int const& timestep)
    {
        updateInner(vNew, vOld, idx, timestep);
    }

    /**
     * Compute E.
     */
    template<typename TT>
    void updateInner2(std::vector< ScaFES::DataField<TT, DIM> >& n,
            std::vector< ScaFES::DataField<TT, DIM> > const& o,
            Index const& idx,
            int const& /*timestep*/)
    {
        // update D from H
        // n.b.:
        // 1) although comfortable, DO NOT USE "auto" for type deduction. ADOL-C yiels obscure and hard to debug compilation failures otherwise ("auto"-deduced rhs type is not assignable).
        // 2) although comfortable, NEVER violate the "compute -> use -> persist" order for (AD-)active variables. ADOL-C will simply not work if you persist a datum (i.e., write it to memory) and re-use it from there. The use of temporaries is a workaround.
        TT Dx=dxstencil_(idx, o[sc(DFs::Dx)],
            o[sc(DFs::Hz)], o[sc(DFs::Psi_HzDy)], n[sc(DFs::Psi_HzDy)],
            o[sc(DFs::Hy)], o[sc(DFs::Psi_HyDz)], n[sc(DFs::Psi_HyDz)]);

        TT Dy=dystencil_(idx, o[sc(DFs::Dy)],
            o[sc(DFs::Hx)], o[sc(DFs::Psi_HxDz)], n[sc(DFs::Psi_HxDz)],
            o[sc(DFs::Hz)], o[sc(DFs::Psi_HzDx)], n[sc(DFs::Psi_HzDx)]);

        TT Dz=dzstencil_(idx, o[sc(DFs::Dz)],
            o[sc(DFs::Hy)], o[sc(DFs::Psi_HyDx)], n[sc(DFs::Psi_HyDx)],
            o[sc(DFs::Hx)], o[sc(DFs::Psi_HxDy)], n[sc(DFs::Psi_HxDy)]);

        TT ga=1./(1.+(o[sc(DFs::SIGMA)](idx))*dt_/EPSILON_0);
        TT gb=(o[sc(DFs::SIGMA)](idx))*dt_/EPSILON_0;

        TT Ex=ga*(Dx-o[sc(DFs::Ix)](idx));
        TT Ey=ga*(Dy-o[sc(DFs::Iy)](idx));
        TT Ez=ga*(Dz-o[sc(DFs::Iz)](idx));

        n[sc(DFs::Ix)](idx)=o[sc(DFs::Ix)](idx)+(gb*Ex);
        n[sc(DFs::Iy)](idx)=o[sc(DFs::Iy)](idx)+(gb*Ey);
        n[sc(DFs::Iz)](idx)=o[sc(DFs::Iz)](idx)+(gb*Ez);

        // write temporaries to dependent fields, see above, item 2
        n[sc(DFs::Dx)](idx)=Dx;
        n[sc(DFs::Dy)](idx)=Dy;
        n[sc(DFs::Dz)](idx)=Dz;
        n[sc(DFs::Ex)](idx)=Ex;
        n[sc(DFs::Ey)](idx)=Ey;

        // inject source, use old time
        if (idx==loc_src_)
            n[sc(DFs::Ez)](idx)=src_((o[sc(DFs::Ez)]).time());
        else
            n[sc(DFs::Ez)](idx)=Ez;

        // collect sink
        if (idx==loc_sink_) {
            convert<CT, TT> c;
            CT val=c(n[sc(DFs::Ez)](idx));
            sink_.push_back(val);
        }
    }

    /**
     * Wire to updateInner2.
     */
    template<typename TT>
    void updateBorder2(std::vector< ScaFES::DataField<TT, DIM> >& vNew,
            std::vector< ScaFES::DataField<TT, DIM>>const& vOld,
            Index const& idx,
            int const& timestep)
    {
        updateInner2(vNew, vOld, idx, timestep);
    }

    /**
     * Empty.
     */
    template<typename TT>
    void evalInner(std::vector< ScaFES::DataField<TT, DIM> >& /*vNew*/,
        ScaFES::Ntuple<int, DIM> const& /*idxNode*/,
        int const& /*timestep*/) { }

    /**
     * Empty.
     */
    template<typename TT>
    void evalBorder(std::vector< ScaFES::DataField<TT, DIM> >& /*vNew*/,
            ScaFES::Ntuple<int, DIM> const& /*idxNode*/,
            int const& /*timestep*/) { }

    /**
     * Write sink data to file.
     */
    void writeSinkData()
    {
        boost::mpi::communicator world;

        if (world.size()!=1) {
            if (world.rank()==0) {
                if (sink_.size()==0)
                    world.recv(boost::mpi::any_source, 0, sink_);
            } else {
                if (sink_.size()!=0)
                    world.send(0, 0, sink_);
            }
        }

        if (world.rank()!=0)
            return;

        const std::string fname=ptree_.get<std::string>("sink.data_file");

        if (fname.size()==0)
            throw std::runtime_error("Sink file name empty.");

        std::cout<<"Writing sink data to file \""
            <<fname<<"\"..."
            <<std::flush;

        std::ofstream of(fname.c_str());
        of.setf(of.scientific);
        of.precision(20);

        for (auto it=sink_.begin(); it!=sink_.end(); ++it)
            of<<(*it)<<std::endl;

        std::cout<<"done."<<std::endl;
    }

private:
    /**
     * Layer structure for datafields with holes.
     */
    std::vector<int> layers(const PTree& p)
    {
        std::vector<int> ret(13, 0);
        const int d=p.get<int>("PML.d")+2;

        for (int i=0; i<12; ++i)
            ret.push_back(d);

        return ret;
    };

    /**
     * Parse into vector.
     * @param s
     * @param sep
     * @return
     */
    template <typename T>
    static std::vector<T> parse(const std::string& s, const char* sep="x")
    {
        boost::char_separator<char> csep(sep);
        boost::tokenizer<boost::char_separator<char> > tok(s, csep);

        std::vector<T> ret;

        for (auto ti=tok.begin(); ti!=tok.end(); ++ti) {
            try {
            ret.push_back(boost::lexical_cast<T>(*ti));
            }
            catch (...) {
                std::cout << "While parsing string \"" << s << "\":" << std::endl;
                throw;
            }
        }
        return ret;
    }

    /**
     * Make Index from vector.
     * @param v
     * @return
     */
    static Index vec2idx(const std::vector<typename Index::value_type>& v)
    {
        Index ret;

        for (unsigned long int i=0; i<Index::DIM; ++i)
            ret[i]=v.at(i);

        return ret;
    }

    /**
     * Cast "nice names" into accessible ints
     * @param df
     * @return
     */
    static size_type sc(const DFs & df)
    {
        return static_cast<size_type> (df);
    }

    /**
     * Load reference data from file.
     */
    static std::vector<CT> load_reference(const std::string& fname, const unsigned int& n)
    {
        boost::mpi::communicator world;
        std::vector<CT> ret;

        if (ForwardOnly==true) {
            if (world.rank()==0)
                std::cout<<"NOT reading reference from file. Calling TF evaluation routines WILL FAIL!"<<std::endl;
            return ret;
        }

        /** TODO: This statement is unreachable in case ForwardOnly=true. */
        if (world.rank()==0) {
            std::ifstream ifs(fname.c_str());

            if (!ifs.is_open())
                throw std::runtime_error("Failed opening file \""+fname+"\".");

            while ((ifs.good()) && (ret.size()<n)) {
                std::string line;
                std::getline(ifs, line);
                try {
                    ret.push_back(boost::lexical_cast<CT, std::string>(line));
                } catch (...) {
                    std::cout<<"While parsing reference file \"" << fname << "\": " <<std::endl;
                    throw;
                }
            }

            if (ret.size()<n)
                throw std::runtime_error("Reference file \""+fname+"\" contains too few items.");
        }

        if (world.size()>1)
            boost::mpi::broadcast(world, ret, 0);

        return ret;
    }

    const PTree ptree_;
    const CT cfln_, ds_, dt_;
    const CPMLParams<CT> cpmlp_;
    const Index first_idx_, last_idx_, loc_src_, loc_sink_;
    const std::vector<CT> reference_;
    const ModulatedGaussianSource<CT> src_;
    std::vector<CT> sink_;

    const HxStencil hxstencil_;
    const HyStencil hystencil_;
    const HzStencil hzstencil_;
    const DxStencil dxstencil_;
    const DyStencil dystencil_;
    const DzStencil dzstencil_;
};

#endif // EMFDTD_HPP_

