/* ScaFES
 * Copyright (c) 2017, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file HeatEqnFDM.hpp
 *
 *  @brief Implementation of n-dimensional heat equation problem on unit hybercube.
 */

#include "ScaFES.hpp"
#include "analyticalSolutions.hpp"

/*******************************************************************************
 ******************************************************************************/
/**
 * \class HeatEqnFDM
 *  @brief Class for discretized heat equation problem.
 *
 * \section heatEqnFDM n-dimensional Heat Equation Problem on Unit Hypercube
 *
 * \subsection mathdescr Mathematical Description
 * Given:
 * <ul>
 * <li> Time interval \f[[t_S; t_E] \mbox{ with } 0 \le t_S < t_E,\f] </li>
 * <li> domain \f[\Omega := (0,1)^3,\f] </li>
 * <li> source \f[f: \bar{\Omega} \times (t_S;t_E] \to R, \quad
        f(x,t) := 1 + \sum_p=^d x_p,\f] </li>
 * <li> boundary condition \f[g: \partial\Omega \times (t_S;t_E] \to R, \quad
         g(x,t) := t \cdot (1 + \sum_p=^d x_p),\f] </li>
 * <li> initial condition \f[\widetilde{y}: \bar{\Omega} \to R, \quad
        \widetilde{y}(x) := t \cdot (1 + \sum_p=^d x_p),\f] </li>
 * </ul>
 * Find \f[y: \bar{\Omega} \times [t_S; t_E] \to {R}\f] such that
 * \f{eqnarray*}{
 * \partial_t y - \Delta y & =&  f  \quad \mbox{in }
 *          \Omega \times (t_S;t_E], \\
 *                      y & =&  g \quad \mbox{on }
 *          \partial\Omega \times (t_S;t_E], \\
 *             y(\cdot,t_S) & =& \widetilde{y}  \quad \mbox{in } \bar{\Omega}.
 * \f}
 *
 * \subsection mathdiscretization Discretization of the problem
 * \subsubsection discrTD Discretization of the time interval and the domain
 * Let the time interval [t_S;t_E] be uniformly discretised with
 * \f[t_l := \{ t_S + l \cdot \tau\}_l \f]
 * with time step size tau > 0.
 *
 * Let the domain Omega be uniformly discretised with
 * \f[ x_{(i,j,k)} :=  \{ (i \cdot h_0, j \cdot h_1, k \cdot h_2) \}_{(i,j,k)}
 * \f]
 * with grid sizes h_p>0 for all p.
 * <ul>
 * <li>G_I: Set of all interior grid nodes, </li>
 * <li>G_B: Set of all boundary grid nodes,</li>
 * <li>G := G_I with G_B: Set of all grid nodes</li>
 * </ul>
 *
 * \subsubsection discrHeatEqn Discretization of the heat equation
 * Define the following vectors:
 * \f{eqnarray*}{
 * Y^{(l)}_{(i,j,k)} &:=& y(x_{(i,j,k)},t_l), \\
 * F^{(l)}_{(i,j,k)} &:=& f(x_{(i,j,k)},t_l), \\
 * G^{(l)}_{(i,j,k)} &:=& g(x_{(i,j,k)},t_l), \\
 * \widetilde{Y}_{(i,j,k)} &:=& \widetilde{y}(x_{(i,j,k)}) \\
 * && \quad \mbox{for all } t_l \in \tau_h, x_{(i,j,k)} \in \Omega_h.
 * \f}
 * Discretise derivatives in space using the symmetric difference quotients
 * (Finite Difference Method with 7-point stencil) and in time using the
 * forward difference quotient (explicit Euler scheme):
 * Discretise in space (Finite Difference Method with 7-point stencil)
 * and in time (explicit Euler scheme):
 * \f{eqnarray*}{
 * Y^{(l+1)}_{(i,j,k)} &=&  Y^{(l)}_{(i,j,k)}
 *                   + \big( Y^{(l)}_{(i-1,j,k)} - 2 \cdot Y^{(l)}_{(i,j,k)}
 *               + Y^{(l)}_{(i+1,j,k)} \big) \cdot \tau / h_0^2   \\
 * &&
 *               \quad  \quad \quad  \quad  + \big(  Y^{(l)}_{(i,j-1,k)} - 2 \cdot Y^{(l)}_{(i,j,k)}
 *               + Y^{(l)}_{(i,j+1,k)} \big) \cdot \tau / h_1^2   \\
 * &&            \quad \quad  \quad \quad  + \big( Y^{(l)}_{(i,j,k-1)}  - 2 \cdot Y^{(l)}_{(i,j,k)}
 *               + Y^{(l)}_{(i,j,k+1)} \big) \cdot \tau / h_2^2  \\
 * &&  \quad \quad \quad\quad + \tau \cdot F^{(l)}_{(i,j,k)}  \quad\quad
 * \forall l, \, \forall {(i,j,k)} \in \cal{G}_I, \\
 * Y^{(l+1)}_{(i,j,k)} &=&  G^{(l+1)}_{(i,j,k)}
 *               \quad   \quad \forall l, \, \forall {(i,j,k)} \in \cal{G}_B, \\
 *            Y^{(0)}_{(i,j,k)} & =&  \widetilde{Y}_{(i,j,k)}\quad \quad
 * \forall {(i,j,k)} \in \cal{G}.
 * \f}
*/
template<typename CT, std::size_t DIM>
class HeatEqnFDM : public ScaFES::Problem<HeatEqnFDM<CT,DIM>, CT, DIM> {

    /* Defines types of equation which can be used for validation. */
    enum typesOfEqn {constant = 0, linear = 1, quadratic = 2, cubic = 3};

    /* Defines types of boundary conditions. */
    enum typesOfBCs {dirichlet = 1, neumann = 2, cauchy = 3};

    /** constant rho. Material parameter (density) for brain. */
    const double RHO = 1.0; /* kg/m^3 */

    /** constant c. Material parameter (specific heat capacity) for brain. */
    const double C = 1.0; /* J/(kg K) */

    /** constant lambda. Material parameter (thermal conductivity)
     * for brain. */
    const double LAMBDA = 1.0; /* W/(m K) */

    /** constant alpha. Material parameter (heat transfer coefficient)
     * for brain. */
    const double ALPHA = 1.0; /* W/(m^2 K) */

    /** constant T_inf. Parameter (ambient temperature)
     * for brain. */
    const double T_INF = 1.0; /* K */

    /** constant q_dot. Material parameter (heat flux)
     * for brain. */
    const double Q_DOT = 1.0; /* W/(m^2) */

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
     * @param eqnDegree_ Type of degree used for test.
     * @param boundaryCond_ Type of boundary condition used for test.
     */
    HeatEqnFDM(ScaFES::Parameters const& params,
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
               std::vector<CT> const& geomparamsInit = std::vector<CT>(),
               int const& eqnDegree_ = 1,
               int const& boundaryCond_ = 1)
        : ScaFES::Problem<HeatEqnFDM<CT, DIM>, CT, DIM>(params, gg, useLeapfrog,
                                                        nameDatafield, stencilWidth,
                                                        isKnownDf, nLayers,
                                                        defaultValue, writeToFile,
                                                        computeError, geomparamsInit),
          eqnDegree(eqnDegree_), boundaryCond(boundaryCond_)
        { }

    /* Member function for type of analytical solution which will be used. */
    int const eqnDegree;

    /* Member function for type of boundary condition which will be used. */
    int const boundaryCond;

    /** Evaluates all fields at one given global inner grid node.
     *  @param vNew Set of all fields.
     *  @param idxNode Index of given grid node.
     */
    void evalInner(std::vector< ScaFES::DataField<CT, DIM> >& vNew,
                   ScaFES::Ntuple<int,DIM> const& idxNode,
                   int const& timestep) {
        ScaFES::Ntuple<double,DIM> x = this->coordinates(idxNode);
        double t = this->time(timestep);
        int eq = this->eqnDegree;

        /* Vector for F. */
        if (eq == constant) {
            vNew[0](idxNode) = RHO * C * consFuncTimeDerivative<CT,DIM>(x);
            vNew[0](idxNode) -= LAMBDA * consFuncSumOfSpaceDerivatives2ndOrder<CT,DIM>(x, t);
        } else if (eq == linear) {
            vNew[0](idxNode) = RHO * C * linFuncTimeDerivative<CT,DIM>(x);
            vNew[0](idxNode) -= LAMBDA * linFuncSumOfSpaceDerivatives2ndOrder<CT,DIM>(x, t);
        } else if (eq == quadratic) {
            vNew[0](idxNode) = RHO * C * quadFuncTimeDerivative<CT,DIM>(x);
            vNew[0](idxNode) -= LAMBDA * quadFuncSumOfSpaceDerivatives2ndOrder<CT,DIM>(x, t);
        } else if (eq == cubic) {
            vNew[0](idxNode) = RHO * C * cubicFuncTimeDerivative<CT,DIM>(x);
            vNew[0](idxNode) -= LAMBDA * cubicFuncSumOfSpaceDerivatives2ndOrder<CT,DIM>(x, t);
        } else {
            vNew[0](idxNode) = 0.0;
        }

        /* Vector for G. */
        /* Analytical solution for G. */
        vNew[1](idxNode) = 0.0;

        /* Vector for U. */
        /* Analytical solution for U. */
        if (eq == constant) {
            vNew[2](idxNode) = consFunc<CT,DIM>(x, t);
        } else if (eq == linear) {
            vNew[2](idxNode) = linFunc<CT,DIM>(x, t);
        } else if (eq == quadratic) {
            vNew[2](idxNode) = quadFunc<CT,DIM>(x, t);
        } else if (eq == cubic) {
            vNew[2](idxNode) = cubicFunc<CT,DIM>(x, t);
        } else {
            vNew[2](idxNode) = 0.0;
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
        ScaFES::Ntuple<double,DIM> x = this->coordinates(idxNode);
        double t = this->time(timestep);
        int eq = this->eqnDegree;
        int bc = this->boundaryCond;

        /* Vector for F. */
        if (bc == dirichlet) {
            if (eq == constant) {
                vNew[0](idxNode) = RHO * C * consFuncTimeDerivative<CT,DIM>(x);
                vNew[0](idxNode) -= LAMBDA * consFuncSumOfSpaceDerivatives2ndOrder<CT,DIM>(x, t);
            } else if (eq == linear) {
                vNew[0](idxNode) = RHO * C * linFuncTimeDerivative<CT,DIM>(x);
                vNew[0](idxNode) -= LAMBDA * linFuncSumOfSpaceDerivatives2ndOrder<CT,DIM>(x, t);
            } else if (eq == quadratic) {
                vNew[0](idxNode) = RHO * C * quadFuncTimeDerivative<CT,DIM>(x);
                vNew[0](idxNode) -= LAMBDA * quadFuncSumOfSpaceDerivatives2ndOrder<CT,DIM>(x, t);
            } else if (eq == cubic) {
                vNew[0](idxNode) = RHO * C * cubicFuncTimeDerivative<CT,DIM>(x);
                vNew[0](idxNode) -= LAMBDA * cubicFuncSumOfSpaceDerivatives2ndOrder<CT,DIM>(x, t);
            } else {
                vNew[0](idxNode) = 0.0;
            }
        } else if (bc == neumann) {
            if (eq == constant) {
                vNew[0](idxNode) = LAMBDA * consFuncSpaceDerivative1stOrder<CT,DIM>(x, t, 0);
                vNew[0](idxNode) += Q_DOT;
            } else if (eq == linear) {
                vNew[0](idxNode) = LAMBDA * linFuncSpaceDerivative1stOrder<CT,DIM>(x, t, 0);
                vNew[0](idxNode) += Q_DOT;
            } else if (eq == quadratic) {
                /* to be implemented. */
                vNew[0](idxNode) = 0.0;
            } else if (eq == cubic) {
                /* to be implemented. */
                vNew[0](idxNode) = 0.0;
            } else {
                vNew[0](idxNode) = 0.0;
            }
        } else if (bc == cauchy) {
            if (eq == constant) {
                vNew[0](idxNode) = LAMBDA * consFuncSpaceDerivative1stOrder<CT,DIM>(x, t, 0);
                vNew[0](idxNode) += ALPHA * (consFunc<CT,DIM>(x, t) - T_INF);
            } else if (eq == linear) {
                vNew[0](idxNode) = LAMBDA * linFuncSpaceDerivative1stOrder<CT,DIM>(x, t, 0);
                vNew[0](idxNode) += ALPHA * (linFunc<CT,DIM>(x, t) - T_INF);
            } else if (eq == quadratic) {
                /* to be implemented. */
                vNew[0](idxNode) = 0.0;
            } else if (eq == cubic) {
                /* to be implemented. */
                vNew[0](idxNode) = 0.0;
            } else {
                vNew[0](idxNode) = 0.0;
            }
        } else {
            vNew[0](idxNode) = 0.0;
        }

        /* Vector for G. */
        /* Analytical solution for G. */
        if (eq == constant) {
            vNew[1](idxNode) = consFunc<CT,DIM>(x, t);
        } else if (eq == linear) {
            vNew[1](idxNode) = linFunc<CT,DIM>(x, t);
        } else if (eq == quadratic) {
            vNew[1](idxNode) = quadFunc<CT,DIM>(x, t);
        } else if (eq == cubic) {
            vNew[1](idxNode) = cubicFunc<CT,DIM>(x, t);
        } else {
            vNew[1](idxNode) = 0.0;
        }

        /* Vector for U. */
        /* Analytical solution for U. */
        if (eq == constant) {
            vNew[2](idxNode) = consFunc<CT,DIM>(x, t);
        } else if (eq == linear) {
            vNew[2](idxNode) = linFunc<CT,DIM>(x, t);
        } else if (eq == quadratic) {
            vNew[2](idxNode) = quadFunc<CT,DIM>(x, t);
        } else if (eq == cubic) {
            vNew[2](idxNode) = cubicFunc<CT,DIM>(x, t);
        } else {
            vNew[2](idxNode) = 0.0;
        }
    }

    /** Initializes all unknown fields at one given global inner grid node.
     *  @param vNew Set of all unknown fields (return value).
     *  @param idxNode Index of given grid node.
     */
    template<typename TT>
    void initInner(std::vector< ScaFES::DataField<TT, DIM> >& vNew,
                   std::vector<TT> const& /*vOld*/,
                   ScaFES::Ntuple<int,DIM> const& idxNode,
                   int const& timestep) {
        ScaFES::Ntuple<double,DIM> const x = this->coordinates(idxNode);
        double t_s = this->time(timestep);
        int eq = this->eqnDegree;

        /* Vector for U. */
        /* Initial condition for U. */
        if (eq == constant) {
            vNew[0](idxNode) = consFunc<CT,DIM>(x, t_s);
        } else if (eq == linear) {
            vNew[0](idxNode) = linFunc<CT,DIM>(x, t_s);
        } else if (eq == quadratic) {
            vNew[0](idxNode) = quadFunc<CT,DIM>(x, t_s);
        } else if (eq == cubic) {
            vNew[0](idxNode) = cubicFunc<CT,DIM>(x, t_s);
        } else {
            vNew[0](idxNode) = 0.0;
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
        vNew[0](idxNode) = vOld[0](idxNode)
                            + this->tau() * (1.0/(RHO*C))
                            * this->knownDf(0, idxNode);
        for (std::size_t pp = 0; pp < DIM; ++pp) {
            vNew[0](idxNode) += this->tau() * (LAMBDA/(RHO*C)) * (
                     vOld[0](this->connect(idxNode, 2*pp))
                     + vOld[0](this->connect(idxNode, 2*pp+1))
                     - 2.0 * vOld[0](idxNode) )
                     / (this->gridsize(pp) * this->gridsize(pp));
        }
    }

    /** Updates all unknown fields at one given global border grid node.
     *  @param vNew Set of all unknown fields at new time step (return value).
     *  @param idxNode Index of given grid node.
     */
    template<typename TT>
    void updateBorder(std::vector<ScaFES::DataField<TT,DIM>>& vNew,
                      std::vector<ScaFES::DataField<TT,DIM>>const& /*vOld*/,
                      ScaFES::Ntuple<int,DIM> const& idxNode,
                      int const& /*timestep*/) {
        vNew[0](idxNode) = this->knownDf(1, idxNode);
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
