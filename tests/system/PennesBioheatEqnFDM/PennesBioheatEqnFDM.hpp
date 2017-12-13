/* ScaFES
 * Copyright (c) 2017, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file PennesBioheatEqnFDM.hpp
 *
 *  @brief Implementation of n-dimensional heat equation problem on unit hybercube.
 */

#ifndef PENNESBIOHEATEQNFDM_HPP_
#define PENNESBIOHEATEQNFDM_HPP_

#include "ScaFES.hpp"

/*******************************************************************************
 ******************************************************************************/
/**
 * \class PennesBioheatEqnFDM
 *  @brief Class for discretized Pennes bioheat equation problem.
*/
template<typename CT, std::size_t DIM, typename Class>
class PennesBioheatEqnFDM : public ScaFES::Problem<PennesBioheatEqnFDM<CT,DIM, Class>, CT, DIM> {

   public:
    /** Coefficient lambda. */
    const double LAMBDA = 1.0;

    /** Coefficient rho. */
    const double RHO = 1.0;

    /** Coefficient c. */
    const double C = 1.0;

    /** Coefficient rho_blood. */
    const double RHO_BLOOD = 1.0;

    /** Coefficient c_blood. */
    const double C_BLOOD = 1.0;

    /** Coefficient w. */
    const double W = 1.0;

    /** Coefficient a. */
    const double COEFF_A = LAMBDA/(RHO * C);

   public:
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
    PennesBioheatEqnFDM(ScaFES::Parameters const& params,
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
        : ScaFES::Problem<PennesBioheatEqnFDM<CT, DIM, Class>, CT, DIM>(params, gg, useLeapfrog,
                                                                        nameDatafield, stencilWidth,
                                                                        isKnownDf, nLayers,
                                                                        defaultValue, writeToFile,
                                                                        computeError, geomparamsInit)
        { }

    /** Evaluates all fields at one given global inner grid node.
     *  @param vNew Set of all fields.
     *  @param idxNode Index of given grid node.
     *
     */
    void evalInner(std::vector< ScaFES::DataField<CT, DIM> >& vNew,
                   ScaFES::Ntuple<int,DIM> const& idxNode,
                   int const& timestep) {
         static_cast<Class*>(this)->evalInner(vNew, idxNode, timestep);
   }

    /** Evaluates all fields at one given global border grid node.
     *  @param vNew Set of all fields.
     *  @param idxNode Index of given grid node.
     *  @param timestep Given time step.
     */
    void evalBorder(std::vector< ScaFES::DataField<CT, DIM> >& vNew,
                    ScaFES::Ntuple<int,DIM> const& idxNode,
                    int const& timestep) {
         static_cast<Class*>(this)->evalBorder(vNew, idxNode, timestep);
    }

    /** Initializes all unknown fields at one given global inner grid node.
     *  @param vNew Set of all unknown fields (return value).
     *  @param idxNode Index of given grid node.
     */
    template<typename TT>
    void initInner(std::vector< ScaFES::DataField<TT, DIM> >& vNew,
                   std::vector<TT> const& vOld,
                   ScaFES::Ntuple<int,DIM> const& idxNode,
                   int const& timestep) {
        static_cast<Class*>(this)->initInner(vNew, vOld, idxNode, timestep);
    }

    /** Initializes all unknown fields at one given global border grid node.
     *  @param vNew Set of all unknown fields (return value).
     *  @param vOld Set of all given fields.
     *  @param idxNode Index of given grid node.
     *  @param timestep Given time step.
     */
    template<typename TT>
    void initBorder(std::vector< ScaFES::DataField<TT, DIM> >& vNew,
                    std::vector<TT> const& vOld,
                    ScaFES::Ntuple<int,DIM> const& idxNode,
                    int const& timestep) {
        static_cast<Class*>(this)->initBorder(vNew, vOld, idxNode, timestep);
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
            vNew[0](idxNode) += this->tau() * COEFF_A * (
                     vOld[0](this->connect(idxNode, 2*pp))
                     + vOld[0](this->connect(idxNode, 2*pp+1))
                     - 2.0 * vOld[0](idxNode) )
                     / (this->gridsize(pp) * this->gridsize(pp));
        }
        vNew[0](idxNode) -= this->tau() * ((RHO_BLOOD*C_BLOOD)/(RHO*C)) * W
                            * vOld[0](idxNode);
        vNew[0](idxNode) += this->tau() * (1.0/(RHO*C)) * this->knownDf(0, idxNode);
    }

    /** Updates all unknown fields at one given global border grid node.
     *  @param vNew Set of all unknown fields at new time step (return value).
     *  @param idxNode Index of given grid node.
     */
    template<typename TT>
    void updateBorder(std::vector<ScaFES::DataField<TT,DIM>>& vNew,
                      std::vector<ScaFES::DataField<TT,DIM>>const& vOld,
                      ScaFES::Ntuple<int,DIM> const& idxNode,
                      int const& timestep) {
        static_cast<Class*>(this)->updateBorder(vNew, vOld, idxNode, timestep);
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
#endif
