/* ScaFES
 * Copyright (c) 2017-2018, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file InitFile.hpp
 *
 *  @brief Test implementation of n-dimensional heat equation problem on unit
 *         hybercube in combination with init data from netCDF file.
 */

#include "ScaFES.hpp"

#include <boost/property_tree/ini_parser.hpp>

/*******************************************************************************
 ******************************************************************************/
/**
 * \class InitFile
 *  @brief Class for discretized heat equation problem.
 *
 * \section heatEqnFDM 3D Heat Equation Problem on Unit Cube
 *
 *
 * \subsection mathdescr Mathematical Description
 * Given:
 * <ul>
 * <li> Time interval \f[[t_S; t_E] \mbox{ with } 0 \le t_S < t_E,\f] </li>
 * <li> domain \f[\Omega := (0,1)^3,\f] </li>
 * <li> source \f[f: \bar{\Omega} \times (t_S;t_E] \to R, \quad
        f(x,t) := 0,\f] </li>
 * <li> boundary condition \f[g: \partial\Omega \times (t_S;t_E] \to R, \quad
         g(x,t) := 0,\f] </li>
 * <li> initial condition \f[\widetilde{y}: \bar{\Omega} \to R, \quad
        \widetilde{y}(x) :=  \prod_{i=0}^{d-1} x_i \cdot
                    (x_i - 1/2)^2 \cdot (1 - x_i),\f] </li>
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
class InitFile : public ScaFES::Problem<InitFile<CT,DIM>, CT, DIM> {
   private:
   using PTree = boost::property_tree::ptree;
   const PTree ptree;

   public:
    /** Coefficient a. */
    const double A;

    /** Coefficient alpha. */
    const double ALPHA;

    /** Coefficient c. */
    const double C;

    /** Coefficient lambda. */
    const double LAMBDA;

    /** Coefficient rho. */
    const double RHO;

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
     * @param ptree_ Config file parser.
     * @param nLayers Number of layers at the global boundary.
     * @param defaultValue Default value of fields.
     * @param writeToFile How often should the data field be written to file.
     * @param computeError Should the Linf error between the numerical
     *                     and exact solution be computed?
     * @param geomparamsInit Initial guess of geometrical parameters.
     */
    InitFile(ScaFES::Parameters const& params,
             ScaFES::GridGlobal<DIM> const& gg,
             bool useLeapfrog,
             std::vector<std::string> const& nameDatafield,
             std::vector<int> const& stencilWidth,
             std::vector<bool> const& isKnownDf,
             const PTree& ptree_,
             std::vector<int> const& nLayers = std::vector<int>(),
             std::vector<CT> const& defaultValue = std::vector<CT>(),
             std::vector<ScaFES::WriteHowOften> const& writeToFile
               = std::vector<ScaFES::WriteHowOften>(),
             std::vector<bool> const& computeError = std::vector<bool>(),
             std::vector<CT> const& geomparamsInit = std::vector<CT>() )
        : ScaFES::Problem<InitFile<CT, DIM>, CT, DIM>(params, gg, useLeapfrog,
                                                      nameDatafield, stencilWidth,
                                                      isKnownDf, nLayers,
                                                      defaultValue, writeToFile,
                                                      computeError, geomparamsInit),
        ptree(ptree_),
        A(ptree.get<CT>("Parameters.a")),
        ALPHA(ptree.get<CT>("Parameters.alpha")),
        C(ptree.get<CT>("Parameters.c")),
        LAMBDA(ptree.get<CT>("Parameters.lambda")),
        RHO(ptree.get<CT>("Parameters.rho"))
        { }

    /** Evaluates all fields at one given global inner grid node.
     *  @param vNew Set of all fields.
     *  @param idxNode Index of given grid node.
     */
    void evalInner(std::vector< ScaFES::DataField<CT, DIM> >& vNew,
                   ScaFES::Ntuple<int,DIM> const& idxNode,
                   int const& /*timestep*/) {
        vNew[0](idxNode) = 1.0;
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
        vNew[0](idxNode) = 1.0;
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
        vNew[0](idxNode) = 100.0;
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
            vNew[0](idxNode) += this->tau() * A * (
                     vOld[0](this->connect(idxNode, 2*pp))
                     + vOld[0](this->connect(idxNode, 2*pp+1))
                     - 2.0 * vOld[0](idxNode) )
                     / (this->gridsize(pp) * this->gridsize(pp));
        }
        vNew[0](idxNode) += this->tau() * (1.0/(RHO*C)) * this->knownDf(0, idxNode);
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
        vNew[0](idxNode) = 100.0;
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
