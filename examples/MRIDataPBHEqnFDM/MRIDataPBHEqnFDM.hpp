/* ScaFES
 * Copyright (c) 2018, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file MRIDataPBHEqnFDM.hpp
 *
 *  @brief Implementation of a n-dimensional Pennes bioheat equation problem
 *         with data from MRI.
 */

#include <iostream>
#include "ScaFES.hpp"

#include <boost/property_tree/ini_parser.hpp>

/******************************************************************************
 *****************************************************************************/
/**
 * \class MRIDataPBHEqnFDM
 *  @brief Class for discretized Pennes bioheat equation problem.
 *
*/
template<typename CT, std::size_t DIM>
class MRIDataPBHEqnFDM : public ScaFES::Problem<MRIDataPBHEqnFDM<CT,DIM>, CT, DIM> {
  private:
    /** Parser for ini files. */
    using PTree = boost::property_tree::ptree;
    const PTree ptree;

  public:
    /** constant h. Ambient convetion. */
    const CT H; /* W/(m^2 K) */

    /** constant T_amb. Ambient temperature. */
    const CT T_AMB; /* K */

    /** constant q_bc. Heat flux at surface inside the brain. */
    const CT Q_BC; /* W/(m^2) */

    /** constant q_skull. Heat flux at upper surface under skull. */
    const CT Q_SKULL; /* W/(m^2) */

    /** constant epsilon. Emissivity. */
    const CT EPSILON; /* - */

    /** constant sigma. Stefan–Boltzmann constant. */
    const CT SIGMA = 5.670373e-8; /* W/(m^2 K^4) */

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
    MRIDataPBHEqnFDM(ScaFES::Parameters const& params,
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
                     std::vector<bool> const& computeError
                         = std::vector<bool>(),
                     std::vector<CT> const& geomparamsInit = std::vector<CT>(),
                     std::vector<bool> const& checkConvergence
                         = std::vector<bool>() )
        : ScaFES::Problem<MRIDataPBHEqnFDM<CT, DIM>, CT, DIM>(params, gg,
                                                              useLeapfrog,
                                                              nameDatafield,
                                                              stencilWidth,
                                                              isKnownDf,
                                                              nLayers,
                                                              defaultValue,
                                                              writeToFile,
                                                              computeError,
                                                              geomparamsInit,
                                                              checkConvergence),
        ptree(ptree_),
        H(ptree.get<CT>("Parameters.H")),
        T_AMB(ptree.get<CT>("Parameters.T_INF")),
        Q_BC(ptree.get<CT>("Parameters.Q_BC")),
        Q_SKULL(ptree.get<CT>("Parameters.Q_SKULL")),
        EPSILON(ptree.get<CT>("Parameters.EPSILON"))
        { }

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
     */
    template<typename TT>
    void initInner(std::vector< ScaFES::DataField<TT, DIM> >& /*vNew*/,
                   std::vector<TT> const& /*vOld*/,
                   ScaFES::Ntuple<int,DIM> const& /*idxNode*/,
                   int const& /*timestep*/) {
    }

    /** Initializes all unknown fields at one given global border grid node.
     */
    template<typename TT>
    void initBorder(std::vector< ScaFES::DataField<TT, DIM> >& /*vNew*/,
                   std::vector<TT> const& /*vOld*/,
                   ScaFES::Ntuple<int,DIM> const& /*idxNode*/,
                   int const& /*timestep*/) {
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
        CT rho = this->knownDf(0, idxNode);
        CT c = this->knownDf(1, idxNode);
        CT k = this->knownDf(2, idxNode);
        CT rho_blood = this->knownDf(3, idxNode);
        CT c_blood = this->knownDf(4, idxNode);
        CT omega = this->knownDf(5, idxNode);
        CT T_blood = this->knownDf(6, idxNode);
        CT q = this->knownDf(7, idxNode);

        /* Discrete Pennes Bioheat Equation for updating inner nodes. */
        vNew[0](idxNode) = vOld[0](idxNode);
        for (std::size_t pp = 0; pp < DIM; ++pp) {
            vNew[0](idxNode) += this->tau() * (k/(rho*c))
                                * (vOld[0](this->connect(idxNode, 2*pp))
                                   + vOld[0](this->connect(idxNode, 2*pp+1))
                                   - 2.0 * vOld[0](idxNode))
                                / (this->gridsize(pp) * this->gridsize(pp));
        }
        vNew[0](idxNode) += this->tau() * ((rho_blood*c_blood)/(rho*c)) * omega
                            * (T_blood - vOld[0](idxNode));
        vNew[0](idxNode) += this->tau() * (1.0/(rho*c)) * q;
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
        CT rho = this->knownDf(0, idxNode);
        CT c = this->knownDf(1, idxNode);
        CT k = this->knownDf(2, idxNode);
        CT rho_blood = this->knownDf(3, idxNode);
        CT c_blood = this->knownDf(4, idxNode);
        CT omega = this->knownDf(5, idxNode);
        CT T_blood = this->knownDf(6, idxNode);
        CT q = this->knownDf(7, idxNode);
        int trepanationArea = this->knownDf(8, idxNode);

        /* Discrete Pennes Bioheat Equation with boundary conditions. */
        vNew[0](idxNode) = vOld[0](idxNode);
        for (std::size_t pp = 0; pp < DIM; ++pp) {
            if (idxNode.elem(pp) == (this->nNodes(pp)-1)) {
            /* vOld[0](this->connect(idxNode, 2*pp+1) needs to be replaced. */
                if (pp == (DIM-1)) {
                /* Last node/edge/surface in highest dimension will be brain surface. */
                    if (trepanationArea == 1) {
                    /* Open skull: Cauchy boundary condition and
                     * thermal radiation boundary condition. */
                        CT tempOld = vNew[0](idxNode) + 273.15;
                        CT tempOldPow4 = tempOld * tempOld * tempOld * tempOld;
                        CT tempAmb = T_AMB + 273.15;
                        CT tempAmbPow4 = tempAmb * tempAmb * tempAmb * tempAmb;

                        vNew[0](idxNode) += this->tau() * (k/(rho*c))
                                            * (vOld[0](this->connect(idxNode, 2*pp))
                                            /* vOld[0](this->connect(idxNode, 2*pp+1)
                                             * is replaced by                         */
                                               + vOld[0](this->connect(idxNode, 2*pp))
                                               - ((2.0*this->gridsize(pp)/k)
                                                  * H * (vOld[0](idxNode) - T_AMB))
                                               - ((2.0*this->gridsize(pp)/k)
                                                  * EPSILON * SIGMA
                                                  * (tempOldPow4 - tempAmbPow4))
                                            /********************************************/
                                               - 2.0 * vOld[0](idxNode))
                                            / (this->gridsize(pp) * this->gridsize(pp));
                    } else {
                    /* Closed skull: Neumann boundary condition. */
                        vNew[0](idxNode) += this->tau() * (k/(rho*c))
                                            * (vOld[0](this->connect(idxNode, 2*pp))
                                            /* vOld[0](this->connect(idxNode, 2*pp+1))
                                               is replaced by                          */
                                               + vOld[0](this->connect(idxNode, 2*pp))
                                               - ((2.0*this->gridsize(pp)/k)
                                                  * (-1.0 * Q_SKULL))
                                            /*********************************************/
                                               - 2.0 * vOld[0](idxNode))
                                            / (this->gridsize(pp) * this->gridsize(pp));
                    }
                } else {
                /* Neumann boundary condition. */
                    vNew[0](idxNode) += this->tau() * (k/(rho*c))
                                        * (vOld[0](this->connect(idxNode, 2*pp))
                                        /* vOld[0](this->connect(idxNode, 2*pp+1))
                                           is replaced by                            */
                                           + vOld[0](this->connect(idxNode, 2*pp))
                                           - ((2.0*this->gridsize(pp)/k)
                                              * (-1.0 * Q_BC))
                                        /*********************************************/
                                           - 2.0 * vOld[0](idxNode))
                                        / (this->gridsize(pp) * this->gridsize(pp));
                }
            } else if (idxNode.elem(pp) == 0){
            /* vOld[0](this->connect(idxNode, 2*pp) needs to be replaced.
             * Neumann boundary condition. */
                vNew[0](idxNode) += this->tau() * (k/(rho*c))
                                    /* vOld[0](this->connect(idxNode, 2*pp))
                                       is replaced by                          */
                                    * (vOld[0](this->connect(idxNode, 2*pp+1))
                                       + ((2.0*this->gridsize(pp)/k) * Q_BC)
                                    /*******************************************/
                                       + vOld[0](this->connect(idxNode, 2*pp+1))
                                       - 2.0 * vOld[0](idxNode))
                                    / (this->gridsize(pp) * this->gridsize(pp));
            } else {
            /* No value needs to be replaced.
             * Use central differencing scheme. */
                vNew[0](idxNode) += this->tau() * (k/(rho*c))
                                    * (vOld[0](this->connect(idxNode, 2*pp))
                                       + vOld[0](this->connect(idxNode, 2*pp+1))
                                       - 2.0 * vOld[0](idxNode))
                                    / (this->gridsize(pp) * this->gridsize(pp));
            }
        }

        /* These terms are independet of the boundary condition. */
        vNew[0](idxNode) += this->tau() * ((rho_blood*c_blood)/(rho*c)) * omega
                            * (T_blood - vOld[0](idxNode));
        vNew[0](idxNode) += this->tau() * (1.0/(rho*c)) * q;
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
