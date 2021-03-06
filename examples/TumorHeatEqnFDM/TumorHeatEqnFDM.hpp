/* ScaFES
 * Copyright (c) 2017-2018, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file TumorHeatEqnFDM.hpp
 *
 *  @brief Implementation of a n-dimensional Pennes bioheat equation problem
 *         with tumor model inside the grid.
 */

#include <iostream>
#include "ScaFES.hpp"

#include <boost/property_tree/ini_parser.hpp>

/*******************************************************************************
 ******************************************************************************/
/**
 * \class TumorHeatEqnFDM
 *  @brief Class for discretized Pennes bioheat equation problem.
 *
*/
template<typename CT, std::size_t DIM>
class TumorHeatEqnFDM : public ScaFES::Problem<TumorHeatEqnFDM<CT,DIM>, CT, DIM> {
  private:
    /** Parser for ini files. */
    using PTree = boost::property_tree::ptree;
    const PTree ptree;

  public:

    /* Based on Bousselham et al. (2017). */

    /************************************************************************/
    /** constant rho. Density. */
    const CT RHO; /* kg/m^3 */

    /** constant C. Specific heat capacity. */
    const CT C; /* J/(kg K) */

    /** constant K. Thermal conductivity). */
    const CT K; /* W/(m K) */

    /** constant Q_brain. Metabolism heat generation of the brain. */
    const CT Q_BRAIN; /* W/(m^3) */

    /** constant Q_tumor. Metabolism heat generation of the tumor. */
    const CT Q_TUMOR; /* W/(m^3) */

    /** constant rho_b. Density of the blood. */
    const CT RHO_B; /* kg/m^3 */

    /** constant C_Pb. Specific heat of the blood. */
    const CT C_PB; /* J/(kg K) */

    /** constant omega_b. Blood perfusion rate (normal brain tissue). */
    const CT OMEGA_B_BRAIN; /* 1/s */

    /** constant omega_b. Blood perfusion rate (Astrocytoma brain tumor). */
    const CT OMEGA_B_TUMOR; /* 1/s */

    /** constant T_i. Initial condition for T. */
    const CT T_I; /* K */

    /** constant h. Ambient convetion. */
    const CT H; /* W/(m^2 K) */

    /** constant T_inf. Air temperature. */
    const CT T_INF; /* K */

    /** constant q_bc. Heat flux at surface. */
    const CT Q_BC; /* W/(m^2) */

    /** constant diameter. diameter of the tumor. */
    const CT DIAMETER; /* m */

    /** constant depth. depth of the tumor. */
    const CT DEPTH; /* m */

    /* constant T_a. Temperature of artery. */
    const CT T_A; /* K */

    /************************************************************************/
    /* Auxiliary values. */
    /* Radius of tumor. */
    CT RADIUS;

    /* Center of tumor. */
    ScaFES::Ntuple<CT,DIM> tumorCenter;

    /************************************************************************/

    /** All fields which are related to the underlying problem
     * are added in terms of an entry of the parameters of
     * type \c std::vector.
     * @param params Set of ScaFES parameters.
     * @param gg Global grid.
     * @param useLeapfrog Should the leap frog scheme be used?
     * @param nameDatafield Name of the fields.
     * @param stencilWidth Stencil width of the fields.
     * @param isKnownDf Is the data field are known or unknown one?
     * @param ptree_ Config file parser.
     * @param nLayers Number of layers at the global boundary.
     * @param defaultValue Default value of fields.
     * @param writeToFile How often should the data field be written to file.
     * @param computeError Should the Linf error between the numerical
     *                     and exact solution be computed?
     * @param geomparamsInit Initial guess of geometrical parameters.
     * @param checkConvergence Should convergence be checked?
     */
    TumorHeatEqnFDM(ScaFES::Parameters const& params,
                    ScaFES::GridGlobal<DIM> const& gg,
                    bool useLeapfrog,
                    std::vector<std::string> const& nameDatafield,
                    std::vector<int> const& stencilWidth,
                    std::vector<bool> const& isKnownDf,
                    PTree const& ptree_,
                    std::vector<int> const& nLayers = std::vector<int>(),
                    std::vector<CT> const& defaultValue = std::vector<CT>(),
                    std::vector<ScaFES::WriteHowOften> const& writeToFile
                      = std::vector<ScaFES::WriteHowOften>(),
                    std::vector<bool> const& computeError = std::vector<bool>(),
                    std::vector<CT> const& geomparamsInit = std::vector<CT>(),
                    std::vector<bool> const& checkConvergence = std::vector<bool>() )
        : ScaFES::Problem<TumorHeatEqnFDM<CT, DIM>, CT, DIM>(params, gg, useLeapfrog,
                                                             nameDatafield, stencilWidth,
                                                             isKnownDf, nLayers,
                                                             defaultValue, writeToFile,
                                                             computeError, geomparamsInit,
                                                             checkConvergence),
        ptree(ptree_),
        RHO(ptree.get<CT>("Parameters.RHO")),
        C(ptree.get<CT>("Parameters.C")),
        K(ptree.get<CT>("Parameters.K")),
        Q_BRAIN(ptree.get<CT>("Parameters.Q_BRAIN")),
        Q_TUMOR(ptree.get<CT>("Parameters.Q_TUMOR")),
        RHO_B(ptree.get<CT>("Parameters.RHO_B")),
        C_PB(ptree.get<CT>("Parameters.C_PB")),
        OMEGA_B_BRAIN(ptree.get<CT>("Parameters.OMEGA_B_BRAIN")),
        OMEGA_B_TUMOR(ptree.get<CT>("Parameters.OMEGA_B_TUMOR")),
        T_I(ptree.get<CT>("Parameters.T_I")),
        H(ptree.get<CT>("Parameters.H")),
        T_INF(ptree.get<CT>("Parameters.T_INF")),
        Q_BC(ptree.get<CT>("Parameters.Q_BC")),
        DIAMETER(ptree.get<CT>("Parameters.DIAMETER")),
        DEPTH(ptree.get<CT>("Parameters.DEPTH")),
        T_A(ptree.get<CT>("Parameters.T_A"))
        {
            RADIUS = DIAMETER/2.0;

            /* Calculate center of tumor. */
            for (std::size_t pp = 0; pp < DIM; ++pp) {
                if (pp == (DIM-1)) {
                    /* In the highest dimension parameter 'depth' will be used. */
                    if (this->params().coordNodeLast()[pp] <= DIAMETER) {
                        std::cerr << "WARNING: Diameter of tumor is bigger than grid in dimension: "
                                  << pp << "." << std::endl;
                    }
                    if (this->params().coordNodeLast()[pp] < std::fabs(DEPTH) ||
                        DEPTH < 0.0) {
                        std::cerr << "WARNING: Center of tumor is outside of grid in dimension: "
                                  << pp << "." << std::endl;
                    }
                    if ((this->params().coordNodeLast()[pp] + RADIUS) < std::fabs(DEPTH)) {
                        std::cerr << "WARNING: Tumor is completely outside of grid." << std::endl;
                    }
                    if (std::fabs(this->params().coordNodeLast()[pp] - DEPTH) < RADIUS ||
                        std::fabs(DEPTH) < RADIUS) {
                        std::cerr << "WARNING: Part of tumor is outside of grid in dimension: "
                                  << pp << "." << std::endl;
                    }
                    tumorCenter[pp] = this->params().coordNodeLast()[pp] - DEPTH;
                } else {
                    /* In every other dimension the tumor will be located in the center. */
                    if (this->params().coordNodeLast()[pp] <= DIAMETER) {
                        std::cerr << "WARNING: Diameter of tumor is bigger than grid in dimension: "
                                  << pp << "." << std::endl;
                    }
                    tumorCenter[pp] = this->params().coordNodeLast()[pp]/2.0;
                }
            }
        }

    /** Evaluates all fields at one given global inner grid node.
     */
    void evalInner(std::vector< ScaFES::DataField<CT, DIM> >& /*vNew*/,
                   ScaFES::Ntuple<int,DIM> const& /*idxNode*/,
                   int const& /*timestep*/) {
    }

    /** Evaluates all fields at one given global border grid node.
     */
    void evalBorder(std::vector< ScaFES::DataField<CT, DIM> >& /*vNew*/,
                    ScaFES::Ntuple<int,DIM> const& /*idxNode*/,
                    int const& /*timestep*/) {
    }

    /** Initializes all unknown fields at one given global inner grid node.
     *  @param vNew Set of all unknown fields (return value).
     *  @param idxNode Index of given grid node.
     */
    template<typename TT>
    void initInner(std::vector< ScaFES::DataField<TT, DIM> >& vNew,
                   std::vector<TT> const& /*vOld*/,
                   ScaFES::Ntuple<int,DIM> const& idxNode,
                   int const& /*timestep*/) {
        /* Init all nodes with the same temperature. */
        vNew[0](idxNode) = T_I;
    }

    /** Initializes all unknown fields at one given global border grid node.
     *  @param vNew Set of all unknown fields (return value).
     *  @param vOld Set of all unknown fields at old time step.
     *  @param idxNode Index of given grid node.
     *  @param timestep Given time step.
     */
    template<typename TT>
    void initBorder(std::vector< ScaFES::DataField<TT, DIM> >& vNew,
                   std::vector<TT> const& vOld,
                   ScaFES::Ntuple<int,DIM> const& idxNode,
                   int const& timestep) {
        this->template initInner<TT>(vNew, vOld, idxNode, timestep);
    }

    /** Updates all unknown fields at one given global inner grid node.
     *  @param vNew Set of all unknown fields at new time step (return value).
     *  @param vOld Set of all unknown fields at old time step.
     *  @param idxNode Index of given grid node.
     */
    template<typename TT>
    void updateInner(std::vector<ScaFES::DataField<TT,DIM>>& vNew,
                     std::vector<ScaFES::DataField<TT,DIM>> const& vOld,
                     ScaFES::Ntuple<int,DIM> const& idxNode,
                     int const& /*timestep*/) {
        CT omega = 0.0;
        CT Q = 0.0;
        /* Get coordinates of current node. */
        ScaFES::Ntuple<double,DIM> x = this->coordinates(idxNode);
        /* Calculate distance to tumor center. */
        CT distance = 0.0;
        for (std::size_t pp = 0; pp < DIM; ++pp) {
            distance += (x[pp] - tumorCenter[pp]) * (x[pp] - tumorCenter[pp]);
        }
        /* Check if current point is inside tumor. */
        if (distance <= (this->RADIUS*this->RADIUS)) {
            /* Inside tumor. */
            omega = OMEGA_B_TUMOR;
            Q = Q_TUMOR;
        } else {
            /* Outside tumor, i.e. normal healthy brain tissue. */
            omega = OMEGA_B_BRAIN;
            Q = Q_BRAIN;
        }
        /* Discrete Pennes Bioheat Equation for updating inner nodes. */
        vNew[0](idxNode) = vOld[0](idxNode);
        for (std::size_t pp = 0; pp < DIM; ++pp) {
            vNew[0](idxNode) += this->tau() * (K/(RHO*C))
                                * (vOld[0](this->connect(idxNode, 2*pp))
                                   + vOld[0](this->connect(idxNode, 2*pp+1))
                                   - 2.0 * vOld[0](idxNode))
                                / (this->gridsize(pp) * this->gridsize(pp));
        }
        vNew[0](idxNode) += this->tau() * ((RHO_B*C_PB)/(RHO*C)) * omega
                            * (T_A - vOld[0](idxNode));
        vNew[0](idxNode) += this->tau() * (1.0/(RHO*C)) * Q;
    }

    /** Updates all unknown fields at one given global border grid node.
     *  @param vNew Set of all unknown fields at new time step (return value).
     *  @param vOld Set of all unknown fields at old time step.
     *  @param idxNode Index of given grid node.
     */
    template<typename TT>
    void updateBorder(std::vector<ScaFES::DataField<TT,DIM>>& vNew,
                      std::vector<ScaFES::DataField<TT,DIM>>const& vOld,
                      ScaFES::Ntuple<int,DIM> const& idxNode,
                      int const& /*timestep*/) {
        CT omega = 0.0;
        CT Q = 0.0;
        /* Get coordinates of current node. */
        ScaFES::Ntuple<double,DIM> x = this->coordinates(idxNode);
        /* Calculate distance to tumor center. */
        CT distance = 0.0;
        for (std::size_t pp = 0; pp < DIM; ++pp) {
            distance += (x[pp] - tumorCenter[pp]) * (x[pp] - tumorCenter[pp]);
        }
        /* Check if current point is inside tumor. */
        if (distance <= (this->RADIUS*this->RADIUS)) {
            /* Inside tumor. */
            omega = OMEGA_B_TUMOR;
            Q = Q_TUMOR;
        } else {
            /* Outside tumor, i.e. normal healthy brain tissue. */
            omega = OMEGA_B_BRAIN;
            Q = Q_BRAIN;
        }
        /* Discrete Pennes Bioheat Equation modified for updating border nodes. */
        vNew[0](idxNode) = vOld[0](idxNode);
        for (std::size_t pp = 0; pp < DIM; ++pp) {
            /* Last node/edge/surface in highest dimension will be convection
             * boundary condition. */
            if (idxNode.elem(pp) == (this->nNodes(pp)-1)) {
            /* vOld[0](this->connect(idxNode, 2*pp+1) needs to be replaced. */
                if (pp == (DIM-1)) {
                /* Cauchy boundary condition. */
                    vNew[0](idxNode) += this->tau() * (K/(RHO*C))
                                        * (vOld[0](this->connect(idxNode, 2*pp))
                                        /* + vOld[0](this->connect(idxNode, 2*pp+1) */
                                           + vOld[0](this->connect(idxNode, 2*pp))
                                           - ((2.0*this->gridsize(pp)/K)
                                              * H * (vOld[0](idxNode) - T_INF))
                                        /********************************************/
                                           - 2.0 * vOld[0](idxNode))
                                        / (this->gridsize(pp) * this->gridsize(pp));
                } else {
                /* Neumann boundary condition. */
                    vNew[0](idxNode) += this->tau() * (K/(RHO*C))
                                        * (vOld[0](this->connect(idxNode, 2*pp))
                                        /* + vOld[0](this->connect(idxNode, 2*pp+1)) */
                                           + vOld[0](this->connect(idxNode, 2*pp))
                                           - ((2.0*this->gridsize(pp)/K) *
                                              (-1.0 * Q_BC))
                                        /*********************************************/
                                           - 2.0 * vOld[0](idxNode))
                                        / (this->gridsize(pp) * this->gridsize(pp));
                }
            } else if (idxNode.elem(pp) == 0){
            /* vOld[0](this->connect(idxNode, 2*pp) needs to be replaced.
             * Neumann boundary condition. */
                vNew[0](idxNode) += this->tau() * (K/(RHO*C))
                                    /* + vOld[0](this->connect(idxNode, 2*pp)) */
                                    * (vOld[0](this->connect(idxNode, 2*pp+1))
                                       + ((2.0*this->gridsize(pp)/K) * Q_BC)
                                    /*******************************************/
                                       + vOld[0](this->connect(idxNode, 2*pp+1))
                                       - 2.0 * vOld[0](idxNode))
                                    / (this->gridsize(pp) * this->gridsize(pp));
            } else {
            /* No value needs to be replaced.
             * Use central differencing scheme. */
                vNew[0](idxNode) += this->tau() * (K/(RHO*C))
                                    * (vOld[0](this->connect(idxNode, 2*pp))
                                       + vOld[0](this->connect(idxNode, 2*pp+1))
                                       - 2.0 * vOld[0](idxNode))
                                    / (this->gridsize(pp) * this->gridsize(pp));
            }
        }

        /* These terms are independet of the boundary condition. */
        vNew[0](idxNode) += this->tau() * ((RHO_B*C_PB)/(RHO*C)) * omega
                            * (T_A - vOld[0](idxNode));
        vNew[0](idxNode) += this->tau() * (1.0/(RHO*C)) * Q;
    }

    /** Updates (2nd cycle) all unknown fields at one given global inner grid node.
     *  \remarks Only important if leap frog scheme is used.
     */
    template<typename TT>
    void updateInner2(std::vector<ScaFES::DataField<TT,DIM>>&,
                      std::vector<ScaFES::DataField<TT,DIM>> const&,
                      ScaFES::Ntuple<int,DIM> const&,
                      int const&) { }

    /** Updates (2nd cycle) all unknown fields at one given global border
     *  grid node.
     *  \remarks Only important if leap frog scheme is used.
     */
    template<typename TT>
    void updateBorder2(std::vector<ScaFES::DataField<TT,DIM>>&,
                       std::vector<ScaFES::DataField<TT,DIM>>const&,
                       ScaFES::Ntuple<int,DIM> const&,
                       int const&) { }
};
