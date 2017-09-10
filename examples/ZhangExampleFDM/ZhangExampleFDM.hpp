/* ScaFES
 * Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file ZhangExampleFDM.hpp
 *
 *  @brief Implementation of n-dimensional Pennes bio heat equation problem on unit cube.
 */

#include <iostream>
#include "ScaFES.hpp"

/*******************************************************************************
 ******************************************************************************/
/**
 * \class ZhangExampleFDM
 *  @brief Class for discretized heat equation problem.
 *
 * \section ZhangExampleFDM 1D Pennes Bioheat Equation Problem based on Zhang (2008)
 *
 *
 * \subsection mathdescr Mathematical Description
 * \subsection mathdiscretization Discretization of the problem
 * \subsubsection discrTD Discretization of the time interval and the domain
 * \subsubsection discrHeatEqn Discretization of Pennes bioheat equation
*/
template<typename CT, std::size_t DIM>
class ZhangExampleFDM : public ScaFES::Problem<ZhangExampleFDM<CT,DIM>, CT, DIM> {
  public:

    /** constant k. Material parameter (thermal conductivity) for tissue. */
    const double K = 0.5; /* W m^-1 K^-1 */

    /** constant eta_b. Material parameter (perfusion) for tissue. */
    const double ETA_B = 10e-4; /* s^-1 */

    /** constant rho_b. Material parameter (density) for blood. */
    const double RHO_B = 1052.0; /* kg m^-3 */

    /** constant c_pb. Material parameter (specific heat capacity) for blood. */
    const double C_PB = 3800.0; /* J kg^-1 K^-1 */

    /** constant T_a. Material parameter (temperature) for blood. */
    const double T_A = 37.0; /* degree Celcius */

    /** constant Q_m. Material parameter (metabolic heat) for tissue. */
    const double Q_M = 400.0; /* W m^-3 */

    /** constant Q_s. Parameter (spatial heating). */
    const double Q_S = 0.0; /* W m^-3 */

    /** constant rho. Material parameter (density) for tissue. */
    const double RHO = 1052.0; /* kg m^-3 */

    /** constant c_p. Material parameter (specific heat capacity) for tissue. */
    const double C_P = 3800.0; /* J kg^-1 K^-1 */

    /** constant T_L. Boundary condition (temperature) for problem. */
    const double T_L = 30.0; /* degree Celcius */

    /** constant L. Geometry parameter (length) for tissue. */
    const double L = 0.04; /* m */

    /** constant T_e. Initial condition (temperature) for problem. */
    const double T_E = T_A + ((Q_M + Q_S)/(ETA_B * RHO_B * C_PB));


    /** All fields which are related to the underlying problem
     * are added in terms of an entry of the parameters of
     * type \c std::vector.
     * @param params Set of ScaFES parameters.
     * @param gg Global grid.
     * @param useLeapfrog Should the leap frog scheme be used?
     * @param nameDatafield Name of the fields.
     * @param stencilWidth Stencil width of the fields.
     * @param isKnownDf Is the data field are known or unknown one?
     * @param nLayers Number of layers at the global boundary.
     * @param defaultValue Default value of fields.
     * @param writeToFile How often should the data field be written to file.
     * @param computeError Should the Linf error between the numerical
     *                     and exact solution be computed?
     * @param geomparamsInit Initial guess of geometrical parameters.
     */
    ZhangExampleFDM(ScaFES::Parameters const& params,
               ScaFES::GridGlobal<DIM> const& gg,
               bool useLeapfrog,
               std::vector<std::string> const& nameDatafield,
               std::vector<int> const& stencilWidth,
               std::vector<bool> const& isKnownDf,
               std::vector<int> const& nLayers = std::vector<int>(),
               std::vector<CT> const& defaultValue = std::vector<CT>(),
               std::vector<ScaFES::WriteHowOften> const& writeToFile
                 = std::vector<ScaFES::WriteHowOften>(),
               std::vector<bool> const& computeError = std::vector<bool>(),
               std::vector<CT> const& geomparamsInit = std::vector<CT>() )
        : ScaFES::Problem<ZhangExampleFDM<CT, DIM>, CT, DIM>(params, gg, useLeapfrog,
                                                        nameDatafield, stencilWidth,
                                                        isKnownDf, nLayers,
                                                        defaultValue, writeToFile,
                                                        computeError, geomparamsInit)
        { }

    /** Evaluates all fields at one given global inner grid node.
     *  @param vNew Set of all fields.
     *  @param idxNode Index of given grid node.
     */
    void evalInner(std::vector< ScaFES::DataField<CT, DIM> >& vNew,
                   ScaFES::Ntuple<int,DIM> const& idxNode,
                   int const& /*timestep*/) {
    }

    /** Evaluates all fields at one given global border grid node.
     *  @param vNew Set of all fields.
     *  @param idxNode Index of given grid node.
     *  @param timestep Given time step.
     */
    void evalBorder(std::vector< ScaFES::DataField<CT, DIM> >& vNew,
                    ScaFES::Ntuple<int,DIM> const& idxNode,
                    int const& timestep) {
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
        /* Uniform distribution as initial condition. */
        vNew[0](idxNode) = T_E; /* knownDf(0, idxNode). T. */
    }

    /** Initializes all unknown fields at one given global border grid node.
     *  @param vNew Set of all unknown fields (return value).
     *  @param vOld Set of all given fields.
     *  @param idxNode Index of given grid node.
     *  @param timestep Given time step.
     */
    template<typename TT>
    void initBorder(std::vector< ScaFES::DataField<TT, DIM> >& vNew,
                   std::vector<TT> const& /*vOld*/,
                   ScaFES::Ntuple<int,DIM> const& idxNode,
                   int const& /*timestep*/) {
        /* Uniform distribution as initial condition. */
        vNew[0](idxNode) = T_E; /* knownDf(0, idxNode). T. */
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
        vNew[0](idxNode) = vOld[0](idxNode);
        for (std::size_t pp = 0; pp < DIM; ++pp) {
            vNew[0](idxNode) += this->tau()
                             * (K/(RHO * C_P))
                             * ((vOld[0](this->connect(idxNode, 2*pp))
                                 + vOld[0](this->connect(idxNode, 2*pp+1))
                                 - 2.0 * vOld[0](idxNode))
                             /(this->gridsize(pp) * this->gridsize(pp)));
        }
        vNew[0](idxNode) += this->tau()
                         * ((ETA_B * RHO_B * C_PB)/(RHO * C_P))
                         * (T_A - vOld[0](idxNode));
        vNew[0](idxNode) += this->tau()
                         * (Q_M/(RHO * C_P));
        vNew[0](idxNode) += this->tau()
                         * (Q_S/(RHO * C_P));
    }

    /** Updates all unknown fields at one given global border grid node.
     *  @param vNew Set of all unknown fields at new time step (return value).
     *  @param idxNode Index of given grid node.
     */
    template<typename TT>
    void updateBorder(std::vector<ScaFES::DataField<TT,DIM>>& vNew,
                      std::vector<ScaFES::DataField<TT,DIM>>const& vOld,
                      ScaFES::Ntuple<int,DIM> const& idxNode,
                      int const& /*timestep*/) {
        /* First node in 1D.
         * Adiabatic boundary conditon. */
        if (idxNode.elem(0) == 0) {
            vNew[0](idxNode) = vOld[0](idxNode);
            for (std::size_t pp = 0; pp < DIM; ++pp) {
                vNew[0](idxNode) += this->tau()
                                 * (K/(RHO * C_P))
                                 * ((2.0 * vOld[0](this->connect(idxNode, 2*pp))
                                     - 2.0 * vOld[0](idxNode))
                                 /(this->gridsize(pp) * this->gridsize(pp)));
            }
            vNew[0](idxNode) += this->tau()
                             * ((ETA_B * RHO_B * C_PB)/(RHO * C_P))
                             * (T_A - vOld[0](idxNode));
            vNew[0](idxNode) += this->tau()
                             * (Q_M/(RHO * C_P));
            vNew[0](idxNode) += this->tau()
                             * (Q_S/(RHO * C_P));
        }
        /* Last node in 1D.
         * Isothermal boundary conditon. */
        if (idxNode.elem(0) == (this->nNodes(0)-1)) {
            vNew[0](idxNode) = T_L;
        }
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
