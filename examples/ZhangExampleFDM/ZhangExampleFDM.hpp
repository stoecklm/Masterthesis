/* ScaFES
 * Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file ZhangExampleFDM.hpp
 *
 *  @brief Implementation of a 1-dimensional Pennes bioheat equation problem
 *  based on Zhang (2008): "Lattice Boltzmann method for solving the bioheat equation".
 */

#define _USE_MATH_DEFINES

#include <cmath>
#include "ScaFES.hpp"

/*******************************************************************************
 ******************************************************************************/
/**
 * \class ZhangExampleFDM
 *  @brief Class for discretized example of a 1D Pennes bioheat equation problem.
 *
 * \section zhangExampleFDM 1D Pennes Bioheat Equation Problem based on Zhang (2008)
 *
 * \subsection mathdescr Mathematical Description
 * Given:
 * <ul>
 * <li> time interval \f[[t_S; t_E] \mbox{ with } 0 \le t_S < t_E,\f] </li>
 * <li> domain \f[\Omega := (0,0.04)^1,\f] </li>
 * <li> source metabolism heat \f[Q_m: \bar{\Omega} \times (t_S;t_E] \to R, \quad
        Q_m(x,t) \neq 0,\f] </li>
 * <li> source spatial heating \f[Q_s: \bar{\Omega} \times (t_S;t_E] \to R, \quad
        Q_s(x,t) := 0,\f] </li>
 * <li> boundary condition (left node) \f[g: \partial\Omega \times (t_S;t_E] \to R, \quad
         \mathrm{adiabatic},\f] </li>
 * <li> boundary condition (right node) \f[g: \partial\Omega \times (t_S;t_E] \to R, \quad
         \mathrm{isothermal},\f] </li>
 * <li> initial condition \f[\widetilde{y}: \bar{\Omega} \to R, \quad
        \widetilde{y}(x) := T_a + \frac{Q_m + Q_s}{\eta_b \rho_b c_{pb}}.\f]
 * </ul>
 * Find \f[y: \bar{\Omega} \times [t_S; t_E] \to {R}\f] such that
 * \f{eqnarray*}{
 *  k \nabla^2 y + \eta_b \rho_b c_{pb} (T_a - y) + Q_m + Q_s
 * &=& \rho c_p \frac{\partial y}{\partial t}\quad \mbox{in }
 *          \Omega \times (t_S;t_E], \\
 *                      y & =&  g \quad \quad \quad \mbox{on }
 *          \partial\Omega \times (t_S;t_E], \\
 *             y(\cdot,t_S) & =& \widetilde{y} \quad \quad  \quad \mbox{in } \bar{\Omega}.
 * \f}
 *
 * \subsection mathdiscretization Discretization of the problem
 * \subsubsection discrTD Discretization of the time interval and the domain
 * Let the time interval [t_S;t_E] be uniformly discretised with
 * \f[t_l := \{ t_S + l \cdot \tau\}_l \f]
 * with time step size tau > 0.
 *
 * Let the domain Omega be uniformly discretised with
 * \f[ x_{(i)} :=  \{ (i \cdot h_0) \}_{(i)}
 * \f]
 * with grid sizes h_p>0 for all p.
 * <ul>
 * <li>G_I: Set of all interior grid nodes, </li>
 * <li>G_BL: Left boundary grid node,</li>
 * <li>G_BR: Right boundary grid node,</li>
 * <li>G := G_I with G_BL and G_BR: Set of all grid nodes</li>
 * </ul>
 *
 * \subsubsection discrZhangEx Discretization of Pennes bioheat equation
 * Define the following vectors:
 * \f{eqnarray*}{
 * T^{(l)}_{(i)} &:=& y(x_{(i)},t_l), \\
 * \widetilde{T}_{(i)} &:=& \widetilde{y}(x_{(i)}) \\
 * && \quad \mbox{for all } t_l \in \tau_h, x_{(i)} \in \Omega_h.
 * \f}
 * Discretize derivative in time using explicit Euler scheme and
 * discretize derivative in space using central difference scheme:
 * \f{eqnarray*}{
 * T^{(l+1)}_{i} &=& T_{i}^{l}
 * + \tau \cdot \frac{k}{\rho c_p}
 * \frac{T^{l}_{i+1} + T^{l}_{i-1} - 2 \cdot T^{l}_{i}}{h_0^2} \\
 * && + \tau \cdot \frac{\eta_b \rho_b c_{pb}}{\rho c_p} (T_a - T^l_i) \\
 * && + \tau \cdot \frac{Q_m}{\rho c_p} \\
 * && + \tau \cdot \frac{Q_s}{\rho c_p}
 * \quad \quad \forall l, \, \forall {(i)} \in \cal{G}_I, \\
 * \f}
 * Left boundary (first node):
 * \f{eqnarray*}{
 * T^{(l+1)}_{i} &=& T_{i}^{l}
 * + \tau \cdot \frac{k}{\rho c_p}
 * \frac{2 \cdot T^{l}_{i+1} - 2 \cdot T^{l}_{i}}{h_0^2} \\
 * && + \tau \cdot \frac{\eta_b \rho_b c_{pb}}{\rho c_p} (T_a - T^l_i) \\
 * && + \tau \cdot \frac{Q_m}{\rho c_p} \\
 * && + \tau \cdot \frac{Q_s}{\rho c_p}
 * \quad \quad \forall l, \, \forall {(i)} \in \cal{G}_{BL}, \\
 * \f}
 * Right boundary (last node):
 * \f[  T^{l+1} = T_L
 * \quad \quad \forall l, \, \forall {(i)} \in \cal{G}_{BR},\f]
 * Initial condition:
 * \f[  T^{(0)}_{(i)}  =  \widetilde{T}_{(i)}
 * \quad \quad \forall {(i)} \in \cal{G}.\f]
 *
 * \subsection analyitcalsol Analytical solution of the problem
 * \f[  T(x,t) = T_e + \frac{2 \cdot \alpha}{L} \cdot (T_L - T_e)
 *      \cdot \sum^{\infty}_{m=1} \bigg((-1)^{m-1} \cdot \beta_m
 *      \cdot \cos(\beta_m \cdot x)
 *      \cdot \frac{1 - e^{-(\alpha \cdot \beta^2_m + \eta^*) \cdot t}}
 *      {\alpha \cdot \beta^2_m + \eta^*} \bigg) \f]
 * with
 * \f[  \alpha = \frac{k}{\rho \cdot c_p}, \f]
 * \f[  \beta_m = \frac{(m - 0.5) \cdot \pi}{L} \ \mathrm{and} \f]
 * \f[  \eta^* = \frac{\eta_b \cdot \rho_b \cdot c_{pb}}{\rho \cdot c_p}. \f]
*/
template<typename CT, std::size_t DIM>
class ZhangExampleFDM : public ScaFES::Problem<ZhangExampleFDM<CT,DIM>, CT, DIM> {
  public:

    /** constant k. Material parameter (thermal conductivity) for tissue. */
    const double K = 0.5; /* W m^-1 K^-1 */

    /** constant eta_b. Material parameter (perfusion) for tissue. */
    const double ETA_B = 1.0e-4; /* s^-1 */

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

    /** constant alpha. Material parameter (thermal diffusivity) for problem. */
    const double ALPHA = K/(RHO * C_P);

    /** constant eta_star. Material parameter for problem. */
    const double ETA_STAR = (ETA_B * RHO_B * C_PB)/(RHO * C_P);

    /** m_max. Mathematical parameter. */
    /* Use M_MAX = 1000000 for SCAFESRUN_END_TIME=200. */
    //const int M_MAX = 1000000;
    /* Use M_MAX = 1000000 for SCAFESRUN_END_TIME=1000. */
    const int M_MAX = 1000000;
    /* Try M_MAX = 1000 for SCAFESRUN_END_TIME=5000. */
    //const int M_MAX = 1000;


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
                   int const& timestep) {
        double beta_m = 0.0;
        double sum = 0.0; /* Value of sum from m=1 to m=m_max-1. */
        double valueCurrIter = 0.0; /* Value for current iteration of loop. */
        double time = timestep * this->tau();
        double x = idxNode.elem(0) * this->gridsize(0);
        if (x == L) {
            x = L - 0.00001;
        }
        for (int m = 1; m < M_MAX; ++m) {
            beta_m = ((m - 0.5) * M_PI)/L;
            valueCurrIter = std::pow(-1, (m-1));
            valueCurrIter *= beta_m;
            valueCurrIter *= std::cos(beta_m * x);
            valueCurrIter *= 1.0 - std::exp(-1.0*(ALPHA * beta_m * beta_m + ETA_STAR) * time);
            valueCurrIter /= (ALPHA * beta_m * beta_m) + ETA_STAR;
            sum += valueCurrIter;
        }
        vNew[0](idxNode) = sum;
        vNew[0](idxNode) *= T_L - T_E;
        vNew[0](idxNode) *= (2.0 * ALPHA)/L;
        vNew[0](idxNode) += T_E;
    }

    /** Evaluates all fields at one given global border grid node.
     *  @param vNew Set of all fields.
     *  @param idxNode Index of given grid node.
     *  @param timestep Given time step.
     */
    void evalBorder(std::vector< ScaFES::DataField<CT, DIM> >& vNew,
                    ScaFES::Ntuple<int,DIM> const& idxNode,
                    int const& timestep) {
        this->evalInner(vNew, idxNode, timestep);
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
                   std::vector<TT> const& vOld,
                   ScaFES::Ntuple<int,DIM> const& idxNode,
                   int const& timestep) {
        /* Uniform distribution as initial condition. */
        this->initInner(vNew, vOld, idxNode, timestep);
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
        vNew[0](idxNode) += this->tau()
                         * (K/(RHO * C_P))
                         * ((vOld[0](this->connect(idxNode, 2*0))
                             + vOld[0](this->connect(idxNode, 2*0+1))
                             - 2.0 * vOld[0](idxNode))
                         /(this->gridsize(0) * this->gridsize(0)));
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
         * Adiabatic boundary conditon.
         * 2*pp+1 is current idxNode+1, i.e. right-hand side neighbour. */
        if (idxNode.elem(0) == 0) {
            vNew[0](idxNode) = vOld[0](idxNode);
            vNew[0](idxNode) += this->tau()
                             * (K/(RHO * C_P))
                             * ((2.0 * vOld[0](this->connect(idxNode, 2*0+1))
                                 - 2.0 * vOld[0](idxNode))
                             /(this->gridsize(0) * this->gridsize(0)));
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
