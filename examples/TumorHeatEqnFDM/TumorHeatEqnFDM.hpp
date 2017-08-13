/* ScaFES
 * Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file TumorHeatEqnFDM.hpp
 *
 *  @brief Implementation of n-dimensional Pennes bio heat equation problem on unit cube.
 */

#include <iostream>
#include "ScaFES.hpp"

/*******************************************************************************
 ******************************************************************************/
/**
 * \class TumorHeatEqnFDM
 *  @brief Class for discretized heat equation problem.
 *
 * \section tumorHeatEqnFDM 3D Pennes Bioheat Equation Problem on unit cube
 *
 *
 * \subsection mathdescr Mathematical Description
 * \subsection mathdiscretization Discretization of the problem
 * \subsubsection discrTD Discretization of the time interval and the domain
 * \subsubsection discrHeatEqn Discretization of Pennes bioheat equation
*/
template<typename CT, std::size_t DIM>
class TumorHeatEqnFDM : public ScaFES::Problem<TumorHeatEqnFDM<CT,DIM>, CT, DIM> {
  public:

    /** constant rho_brain. Material parameter (density) for the brain. */
    const double RHO_BRAIN = 1040.0; /* kg/m^3 */

    /** constant c_brain. Material parameter (specific heat capacity)
     * for the brain. */
    const double C_BRAIN = 3650.0; /* J/(kg K) */

    /** constant lambda_brain. Material parameter (thermal conductivity)
     * for the brain. */
    const double LAMBDA_BRAIN = 0.6; /* W/(m K) */

    /** constant rho_blood. Material parameter (density) for blood. */
    const double RHO_BLOOD = 1047.0; /* kg/m^3 */

    /** constant c_blood. Material parameter (specific heat capacity)
     * for blood. */
    const double C_BLOOD = 3600.0; /* J/(kg K) */

    /** constant w_brain. Material parameter (perfusion) for the brain. */
    const double W_BRAIN = 0.0399; /* 1/s */

    /** constant T_brain_init. Parameter (temperature) for
     * initalization of brain temperature. */
    const double T_BRAIN_INIT = 306.25; /* K */

    /** constant T_tumor_init. Parameter (temperature) for
     * initalization of tumor temperature.
     * Tumor = Anaplastic strocytoma. */
    const double T_TUMOR_INIT = 312.35; /* K */

    /** constant T_blood. Parameter (temperature) for blood. */
    const double T_BLOOD = 310.15; /* K */

    /** constant T_amb. Parameter (temperature) for ambient. */
    const double T_AMB = 291.15; /* K */

    /** constant alpha_amb. Parameter (heat transfer coefficent)
     for ambient. */
    const double ALPHA_AMB = 10.0; /* W/(m^2 K) */

    /** constant q_M_brain. Material parameter (metabolism) for the brain. */
    const double Q_DOT_M_BRAIN = 49937.0; /* W/(m^3) */

    /** constant sphereRadius. Geometry parameter (radius) for tumor. */
    const double sphereRadius = 0.025;

    /** constant x0. Geometry parameter (x coordiante) for tumor. */
    const double x0 = 0.05;

    /** constant y0. Geometry parameter (y coordiante) for tumor. */
    const double y0 = 0.05;

    /** constant z0. Geometry parameter (z coordiante) for tumor. */
    const double z0 = 0.05;

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
        vNew[0](idxNode) = RHO_BRAIN;     /* knownDf(0, idxNode). rho.     */
        vNew[1](idxNode) = C_BRAIN;       /* knownDf(1, idxNode). c.       */
        vNew[2](idxNode) = LAMBDA_BRAIN;  /* knownDf(2, idxNode). lambda.  */
        vNew[3](idxNode) = W_BRAIN;       /* knownDf(3, idxNode). w.       */
        vNew[4](idxNode) = Q_DOT_M_BRAIN; /* knownDf(4, idxNode). q_dot_m. */
        /* vNew[5](idxNode) = 0.0;           knownDf(5, idxNode). T.       */
    }

    /** Evaluates all fields at one given global border grid node.
     *  @param vNew Set of all fields.
     *  @param idxNode Index of given grid node.
     *  @param timestep Given time step.
     */
    void evalBorder(std::vector< ScaFES::DataField<CT, DIM> >& vNew,
                    ScaFES::Ntuple<int,DIM> const& idxNode,
                    int const& timestep) {
        vNew[0](idxNode) = RHO_BRAIN;     /* knownDf(0, idxNode). rho.     */
        vNew[1](idxNode) = C_BRAIN;       /* knownDf(1, idxNode). c.       */
        vNew[2](idxNode) = LAMBDA_BRAIN;  /* knownDf(2, idxNode). lambda.  */
        vNew[3](idxNode) = W_BRAIN;       /* knownDf(3, idxNode). w.       */
        vNew[4](idxNode) = Q_DOT_M_BRAIN; /* knownDf(4, idxNode). q_dot_m. */
        /* vNew[5](idxNode) = 0.0;           knownDf(5, idxNode). T.       */
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
        double sphereCo =
            (this->coordinates(idxNode).elem(0) - x0)
            * (this->coordinates(idxNode).elem(0) - x0)
            + (this->coordinates(idxNode).elem(1) - y0)
            * (this->coordinates(idxNode).elem(1) - y0)
            + (this->coordinates(idxNode).elem(2) - z0)
            * (this->coordinates(idxNode).elem(2) - z0);
        if (sphereCo <= (sphereRadius*sphereRadius)) {
            vNew[0](idxNode) = T_TUMOR_INIT; /* knownDf(0, idxNode). T. */
        } else {
            vNew[0](idxNode) = T_BRAIN_INIT; /* knownDf(0, idxNode). T. */
        }
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
        vNew[0](idxNode) = T_BRAIN_INIT; /* knownDf(0, idxNode). T. */
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
        /* knownDf(0, idxNode). rho.     */
        /* knownDf(1, idxNode). c.       */
        /* knownDf(2, idxNode). lambda.  */
        /* knownDf(3, idxNode). w.       */
        /* knownDf(4, idxNode). q_dot_m. */
        /* knownDf(5, idxNode). T.       */
        double sphereCo =
            (this->coordinates(idxNode).elem(0) - x0)
            * (this->coordinates(idxNode).elem(0) - x0)
            + (this->coordinates(idxNode).elem(1) - y0)
            * (this->coordinates(idxNode).elem(1) - y0)
            + (this->coordinates(idxNode).elem(2) - z0)
            * (this->coordinates(idxNode).elem(2) - z0);
        if (sphereCo <= (sphereRadius*sphereRadius)) {
            /* Do nothing. */
        } else {
            vNew[0](idxNode) = vOld[0](idxNode);
            for (std::size_t pp = 0; pp < DIM; ++pp) {
                vNew[0](idxNode) += this->tau()
                * (this->knownDf(2, idxNode)
                    / (this->knownDf(0, idxNode) * this->knownDf(1, idxNode)))
                * ((vOld[0](this->connect(idxNode, 2*pp))
                   + vOld[0](this->connect(idxNode, 2*pp+1))
                   - 2.0 * vOld[0](idxNode))
                    / (this->gridsize(pp) * this->gridsize(pp)));
            }
            vNew[0](idxNode) += this->tau() * (
                (RHO_BLOOD * C_BLOOD)
                / (this->knownDf(0, idxNode) * this->knownDf(1, idxNode)))
            * (this->knownDf(3, idxNode)
            * (T_BLOOD - vOld[0](idxNode)));
            vNew[0](idxNode) += this->tau() * (
                (1.0/(this->knownDf(0, idxNode) * this->knownDf(1, idxNode))))
            * this->knownDf(4, idxNode);
        }
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
        if (idxNode.elem(2) == (this->nNodes(2)-1)) {
            /* knownDf(0, idxNode). rho.    */
            /* knownDf(1, idxNode). c.      */
            /* knownDf(2, idxNode). lambda. */
            /* knownDf(3, idxNode). w.      */
            /* knownDf(4, idxNode). q_dot_m.    */
            /* knownDf(5, idxNode). T.      */
            vNew[0](idxNode) = vOld[0](idxNode);
            /* First and second dimension. */
            for (std::size_t pp = 0; pp < 2; ++pp) {
                vNew[0](idxNode) += this->tau()
                * (this->knownDf(2, idxNode)
                    / (this->knownDf(0, idxNode) * this->knownDf(1, idxNode)))
                * ((vOld[0](this->connect(idxNode, 2*pp))
                   + vOld[0](this->connect(idxNode, 2*pp+1))
                   - 2.0 * vOld[0](idxNode))
                    / (this->gridsize(pp) * this->gridsize(pp)));
            }
            /* Third dimension. */
            vNew[0](idxNode) += this->tau()
            * (this->knownDf(2, idxNode)
                / (this->knownDf(0, idxNode) * this->knownDf(1, idxNode)))
            * (((2.0 * this->gridsize(2) * (ALPHA_AMB/this->knownDf(2, idxNode)) * T_AMB)
                + (2.0 * vOld[0](this->connect(idxNode, 2*2)))
                - (2.0 * vOld[0](idxNode) * (this->gridsize(2) * (ALPHA_AMB/this->knownDf(2, idxNode)) + 1.0)))
                / (this->gridsize(2) * this->gridsize(2)));
            vNew[0](idxNode) += this->tau() * (
                (RHO_BLOOD * C_BLOOD)
                / (this->knownDf(0, idxNode) * this->knownDf(1, idxNode)))
            * (this->knownDf(3, idxNode)
            * (T_BLOOD - vOld[0](idxNode)));
            vNew[0](idxNode) += this->tau() * (
                (1.0/(this->knownDf(0, idxNode) * this->knownDf(1, idxNode))))
            * this->knownDf(4, idxNode);
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
