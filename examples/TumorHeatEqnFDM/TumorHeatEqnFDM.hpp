/* ScaFES
 * Copyright (c) 2017, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file TumorHeatEqnFDM.hpp
 *
 *  @brief Implementation of a 3-dimensional Pennes bioheat equation problem.
 */

#include <iostream>
#include "ScaFES.hpp"

/*******************************************************************************
 ******************************************************************************/
/**
 * \class TumorHeatEqnFDM
 *  @brief Class for discretized Pennes bioheat equation problem.
 *
*/
template<typename CT, std::size_t DIM>
class TumorHeatEqnFDM : public ScaFES::Problem<TumorHeatEqnFDM<CT,DIM>, CT, DIM> {
  public:

    /* Values from Bousselham et al. (2017). */

    /************************************************************************/
    /* Values from Table 1: Thermophysical properties. */
    /** constant rho. Density. */
    const double RHO = 1040.0; /* kg/m^3 */

    /** constant C. Specific heat capacity. */
    const double C = 3650.0; /* J/(kg K) */

    /** constant K. Thermal conductivity). */
    const double K = 0.6; /* W/(m K) */

    /** constant Q. Metabolism heat generation. */
    const double Q = 25000.0; /* W/(m^3) */

    /** constant rho_b. Density of the blood. */
    const double RHO_B = 1052.0; /* kg/m^3 */

    /** constant C_Pb. Specific heat of the blood. */
    const double C_PB = 3800.0; /* J/(kg K) */

    /** constant omega_b. Blood perfusion rate (normal brain tissue). */
    const double OMEGA_B_BRAIN = 0.004; /* 1/s */

    /** constant omega_b. Blood perfusion rate (Astrocytoma brain tumor).
      * Table 1 presents an interval. Exact value is given in text.
      */
    const double OMEGA_B_TUMOR = 0.0007; /* 1/s */

    /************************************************************************/
    /* Values given in text. */
    /** constant T_i. Initial condition for T. */
    const double T_I = 37.0; /* K */

    /** constant h. Ambient convetion. */
    const double H = 10.0; /* W/(m^2 K) */

    /** constant T_inf. Air temperature. */
    const double T_inf = 22.4; /* K */

    /************************************************************************/
    /* Values varied in text. */
    /** constant diameter. diameter of the tumor. */
    const double DIAMETER = 0.02; /* m */

    /** constant depth. depth of the tumor. */
    const double DEPTH = 0.008; /* m */

    /************************************************************************/
    /* Values missing in text. */
    /* constant T_a. Temperature of artery.
     * Value from Das et al.: Numerical estimation... (2013). */
    const double T_a = 37.0; /* K */

    /************************************************************************/
    /* Auxiliary values. */
    /* Center of tumor = (x0, y0, z0). */
    const double x0 = this->params().coordNodeLast()[0]/2.0;
    const double y0 = this->params().coordNodeLast()[1]/2.0;
    const double z0 = this->params().coordNodeLast()[2]-DEPTH;

    /* Radius of tumor. */
    const double RADIUS = DIAMETER/2.0;

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
     * @param nLayers Number of layers at the global boundary.
     * @param defaultValue Default value of fields.
     * @param writeToFile How often should the data field be written to file.
     * @param computeError Should the Linf error between the numerical
     *                     and exact solution be computed?
     * @param geomparamsInit Initial guess of geometrical parameters.
     */
    TumorHeatEqnFDM(ScaFES::Parameters const& params,
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
        : ScaFES::Problem<TumorHeatEqnFDM<CT, DIM>, CT, DIM>(params, gg, useLeapfrog,
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
        /* Get coordinates of current node. */
        ScaFES::Ntuple<double,DIM> x = this->coordinates(idxNode);

        /* Calculate distance to tumor center. */
        double distance = ((x[0] - this->x0) * (x[0] - this->x0))
                          + ((x[1] - this->y0) * (x[1] - this->y0))
                          + ((x[2] - this->z0) * (x[2] - this->z0));

        /* Check if current point is inside tumor. */
        if (distance <= (this->RADIUS*this->RADIUS)) {
            /* Inside tumor. */
            vNew[0](idxNode) = OMEGA_B_TUMOR;
        } else {
            /* Outside tumor, i.e. normal healthy brain tissue. */
            vNew[0](idxNode) = OMEGA_B_BRAIN;
        }
    }

    /** Evaluates all fields at one given global border grid node.
     *  @param vNew Set of all fields.
     *  @param idxNode Index of given grid node.
     *  @param timestep Given time step.
     */
    void evalBorder(std::vector< ScaFES::DataField<CT, DIM> >& vNew,
                    ScaFES::Ntuple<int,DIM> const& idxNode,
                    int const& timestep) {
        /* Assuming border is always healthy brain tissue. */
        vNew[0](idxNode) = OMEGA_B_BRAIN;
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
    }

    /** Initializes all unknown fields at one given global border grid node.
     *  @param vNew Set of all unknown fields (return value).
     *  @param idxNode Index of given grid node.
     */
    template<typename TT>
    void initBorder(std::vector< ScaFES::DataField<TT, DIM> >& vNew,
                   std::vector<TT> const& /*vOld*/,
                   ScaFES::Ntuple<int,DIM> const& idxNode,
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
