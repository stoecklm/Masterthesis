
/* ScaFES
 * Copyright (c) 2016-2018, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file LameNavierEqnFDMFirstOrder.hpp
 *
 *  @brief Implementation of n-dimensional Lame Navier equations problem on unit hybercube.
 */
#define _USE_MATH_DEFINES
#include <cmath>

#include "ScaFES.hpp"

/*******************************************************************************
 ******************************************************************************/
/**
 * \class LameNavierEqnFDMFirstOrder
 *  @brief Class for discretized Lame Navier problem Represented as a 1st Order System.
 *
 * \section lameNavierEqnFDMFirstOrder 3D Lame Navier Problem Represented as a 1st Order System on Unit Cube
 *
 * \warning Currently, ScaFES is not capable to use other values than the ones of
 * the last time step. However, it is necessary for the underlying FD discretization
 * to access values of the last two time steps. Thus, this example will currently not work. [KF], 2016-08-29.
 *
 * \subsection mathdescr Mathematical Description
 * Given:
 * <ul>
 * <li> Space dimension \f$d \in N,\f$ </li>
 * <li> time interval \f$[t_S; t_E] \mbox{ with } 0 \le t_S < t_E,\f$ </li>
 * <li> domain \f$\Omega := (0,1)^d,\f$ </li>
 * <li> source \f$f: \bar{\Omega} \times (t_S;t_E] \to R^d, \quad
        f(x,t) := 0,\f$ </li>
 * <li> boundary condition \f$g: \partial\Omega \times (t_S;t_E] \to R^d, \quad
         g(x,t) := 0,\f$ </li>
 * <li> initial condition \f$\widetilde{y}: \bar{\Omega} \to R^d, \quad
        \widetilde{y}(x) :=  \mbox{TBD},\f$ </li>
 * <li> initial condition \f$\widehat{y}: \bar{\Omega} \to R, \quad
        \widehat{y}(x) := \mbox{TBD},\f$ </li>
 * <li> constant density \f$\rho \in R_{\neq 0},\f$ </li>
 * <li> constant material parameters \f$\mu, \lambda \in R_{\neq 0}.\f$ </li>
 * </ul>
 * Find \f$y: \bar{\Omega} \times [t_S; t_E] \to R^d\f$ such that
 * \f{eqnarray*}{
 * \rho \cdot \partial^2_{tt} (y) + \mu \cdot \mbox{div} ( \nabla (y) ) + (\mu + \lambda) \cdot (\nabla (\mbox{div} (y) ) )^{\dagger} & =&  f  \quad \mbox{in }
 *          \Omega \times (t_S;t_E], \\
 *                      y & =&  g \quad \mbox{on }
 *          \partial\Omega \times (t_S;t_E], \\
 *             y(\cdot,t_S) & =& \widetilde{y}  \quad \mbox{in } \bar{\Omega}, \\
 *            \partial_t (y)(\cdot,t_S) & =& \widehat{y}  \quad \mbox{in } \bar{\Omega}.
 * \f}
 * Transform into system of first order wrt. time:
 * Introduce new functions: \f$ u:=y, \f$ \f$ v:=\partial_t (y). \f$
 *
 * Find \f$[u;v]: \bar{\Omega} \times [t_S; t_E] \to R^{2d}\f$ such that
 * \f{eqnarray*}{
 *             \partial_t (u) & =& v  \quad \mbox{in }
 *          \bar{\Omega} \times (t_S;t_E], \\
 * \rho \cdot \partial_{t} (v) + \mu \cdot \mbox{div} ( \nabla (u) ) + (\mu + \lambda) \cdot (\nabla (\mbox{div} (u) ) )^{\dagger} & =&  f  \quad \mbox{in }
 *          \Omega \times (t_S;t_E], \\
 *                      u & =&  g \quad \mbox{on }
 *          \partial\Omega \times (t_S;t_E], \\
 *             u(\cdot,t_S) & =& \widetilde{y}  \quad \mbox{in } \bar{\Omega}, \\
 *            \partial_t (u)(\cdot,t_S) & =& \widehat{y}  \quad \mbox{in } \bar{\Omega}.
 * \f}
 * Find \f$[u;v]: \bar{\Omega} \times [t_S; t_E] \to R^{2d}\f$ such that
 * \f{eqnarray*}{
 *             \partial_t (u) & =& v  \quad \mbox{in }
 *          \bar{\Omega} \times (t_S;t_E], \\
 * \partial_{t} (v) + 1/\rho \cdot \bigl( \mu \cdot \mbox{div} ( \nabla ) + (\mu + \lambda) \cdot (\nabla (\mbox{div}  ) )^{\dagger}\bigr) (u) & =&  f / \rho  \quad \mbox{in }
 *          \Omega \times (t_S;t_E], \\
 *                      u & =&  g \quad \mbox{on }
 *          \partial\Omega \times (t_S;t_E], \\
 *             u(\cdot,t_S) & =& \widetilde{y}  \quad \mbox{in } \bar{\Omega}, \\
 *            \partial_t (u)(\cdot,t_S) & =& \widehat{y}  \quad \mbox{in } \bar{\Omega}.
 * \f}
 * Find \f$[u;v]: \bar{\Omega} \times [t_S; t_E] \to R^{2d}\f$ such that
 * \f{eqnarray*}{
  \partial_t
\left[
\begin{array}{c}
  u \\
  v \\
\end{array}
\right]
 +
\left[
\begin{array}{cc}
  0 & -id \\
  1/\rho \cdot \bigl( \mu \cdot \mbox{div} ( \nabla ) + (\mu + \lambda) \cdot (\nabla (\mbox{div}  ) )^{\dagger}\bigr) &  0  \\
\end{array}
\right]
\left[
\begin{array}{c}
  u \\
  v \\
\end{array}
\right]
&=&
\left[
\begin{array}{c}
  0 \\
  f/\rho \\
\end{array}
\right]
\quad \mbox{in }  \Omega \times (t_S;t_E], \\
 *                      \partial_t u & =&  v \quad \mbox{on }
 *          \partial\Omega \times (t_S;t_E), \\
 *                      u & =&  g \quad \mbox{on }
 *          \partial\Omega \times (t_S;t_E), \\
\left[
\begin{array}{c}
  u \\
  v \\
\end{array}
\right](\cdot,t_S)
 & =&
 \left[
\begin{array}{c}
  \widetilde{y} \\
  \widehat{y} \\
\end{array}
\right]
\quad \mbox{in } \bar{\Omega}.
 * \f}
*
*
 * \subsection mathdiscretization Discretization of the problem
 * \subsubsection discrTD Discretization of the time interval and the domain
 * Let the time interval \f$[t_S;t_E]\f$ be uniformly discretized in each direction with
 * \f[t_l := \{ t_S + l \cdot \tau\}_l \f]
 * with time step size \f$\tau > 0\f$.
 *
 * Let the domain \f$\Omega\f$ be uniformly discretized with
 * \f[ x_{(i_0,i_1,\ldots, i_{d-1})} :=  \{ (i_0 \cdot h_0, i_1 \cdot h_1,\ldots, i_{d-1} \cdot h_{d-1}) \}_{(i_0,\ldots,i_{d-1})}
 * \f]
 * with grid sizes \f$h_p>0\f$ for all \f$p\in \{0,1,\ldots,d-1\}\f$.
 * <ul>
 * <li>\f$\mathcal{G}_I\f$: Set of all interior grid nodes, </li>
 * <li>\f$\mathcal{G}_B\f$: Set of all boundary grid nodes,</li>
 * <li>\f$\mathcal{G} := \mathcal{G}_I \cup \mathcal{G}_B\f$: Set of all grid nodes.</li>
 * </ul>
 *
 * \subsubsection discrSpaceLameNavierEqn Semi-discretization wrt. space of the Lame Navier equations
 * Define the following vectors:
 * \f{eqnarray*}{
 * Y^{(p)}_{(i,j,k)}(t) &:=& y_p(x_{(i,j,k)},t), \\
 * F^{(p)}_{(i,j,k)}(t) &:=& f_p(x_{(i,j,k)},t), \\
 * G^{(p)}_{(i,j,k)}(t) &:=& g_p(x_{(i,j,k)},t), \\
 * \widetilde{Y}^{(p)}_{(i,j,k)} &:=& \widetilde{y}_p(x_{(i,j,k)}), \\
 * \widehat{Y}^{(p)}_{(i,j,k)} &:=& \widehat{y}_p(x_{(i,j,k)}) \\
 * && \quad \mbox{for all } p \in \{0,1,\ldots,d-1\}, \mbox{for all } x_{(i,j,k)} \in \Omega_h.
 * \f}
 * Partial derivatives of 2nd order wrt.\ same direction will be discretized using the symmetric difference quotients, resulting in the well-known 7-point stencil in 3-D.
 * Mixed partial derivatives of 2nd order wrt.\ space will be discretized as follows:
 * Firstly, the forward difference quotient will be applied, secondly, the backward difference quotient.
 *
 * Find \f$Y^{(\cdot)}: [t_S;t_E] \to R^{d \cdot N}\f$ such that
 * \f{eqnarray*}{
 * \frac{d}{dt} Y^{(0)}_{(i,j,k)}(t) &=& \big( Y^{(0)}_{(i-1,j,k)}(t) - 2 \cdot Y^{(0)}_{(i,j,k)}(t)
 *               + Y^{(0)}_{(i+1,j,k)}(t) \big) \cdot (2\cdot \mu + \lambda ) /(\rho \cdot  h_0^2)   \\
 * &&
 *               \quad  \quad \quad  \quad  + \big(  Y^{(0)}_{(i,j-1,k)}(t) - 2 \cdot Y^{(0)}_{(i,j,k)}(t)
 *               + Y^{(0)}_{(i,j+1,k)}(t) \big) \cdot \mu  / (\rho \cdot h_1^2)   \\
 * &&            \quad \quad  \quad \quad  + \big(  Y^{(0)}_{(i,j,k-1)}(t) - 2 \cdot Y^{(0)}_{(i,j,k)}(t)
 *               + Y^{(0)}_{(i,j,k+1)}(t) \big)  \cdot \mu  / (\rho \cdot h_2^2)  \\
 * &&            \quad \quad  \quad \quad  + \big(  Y^{(0)}_{(i,j+1,k)}(t) -  Y^{(0)}_{(i-1,j+1,k)}(t)
 *               - Y^{(0)}_{(i,j,k)}(t) +  Y^{(0)}_{(i-1,j,k)}  \big)  \cdot (\mu + \lambda)   / (\rho \cdot h_0 \cdot h_1)  \\
 * &&            \quad \quad  \quad \quad  + \big(  Y^{(0)}_{(i,j,k+1)}(t) -  Y^{(0)}_{(i-1,j,k+1)}(t)
 *               - Y^{(0)}_{(i,j,k)}(t) +  Y^{(0)}_{(i-1,j,k)}(t)  \big)  \cdot (\mu + \lambda)  / (\rho \cdot h_0 \cdot h_2)  \\
 * &&  \quad \quad \quad\quad + 1 / \rho  \cdot F^{(0;l)}_{(i,j,k)}(t)  \quad\quad
 * \forall t, \, \forall {(i,j,k)} \in \mathcal{G}_I, \\
 * \frac{d}{dt} Y^{(1)}_{(i,j,k)}(t) &=&  \big( Y^{(1)}_{(i-1,j,k)}(t) - 2 \cdot Y^{(1)}_{(i,j,k)}(t)
 *               + Y^{(1)}_{(i+1,j,k)}(t) \big) \cdot \mu  /(\rho \cdot  h_0^2)   \\
 * &&
 *               \quad  \quad \quad  \quad  + \big(  Y^{(1)}_{(i,j-1,k)}(t) - 2 \cdot Y^{(1)}_{(i,j,k)}(t)
 *               + Y^{(1)}_{(i,j+1,k)}(t) \big) \cdot  (2\cdot \mu + \lambda )  / (\rho \cdot h_1^2)   \\
 * &&            \quad \quad  \quad \quad  + \big(  Y^{(1)}_{(i,j,k-1)}(t) - 2 \cdot Y^{(1)}_{(i,j,k)}(t)
 *               + Y^{(1)}_{(i,j,k+1)} \big)  \cdot \mu  / (\rho \cdot h_2^2)  \\
 * &&            \quad \quad  \quad \quad  + \big(  Y^{(1)}_{(i+1,j,k)}(t) -  Y^{(1)}_{(i+1,j-1,k)}(t)
 *               - Y^{(1)}_{(i,j,k)}(t) +  Y^{(1)}_{(i,j-1,k)}(t)  \big)  \cdot (\mu + \lambda)   / (\rho \cdot h_0 \cdot h_1)  \\
 * &&            \quad \quad  \quad \quad  + \big(  Y^{(1)}_{(i,j+1,k)}(t) -  Y^{(1)}_{(i,j+1,k-1)}(t)
 *               - Y^{(1)}_{(i,j,k)}(t) +  Y^{(1)}_{(i,j,k-1)}(t)  \big)  \cdot (\mu + \lambda) / (\rho \cdot h_0 \cdot h_2)  \\
 * &&  \quad \quad \quad\quad + 1 / \rho  \cdot F^{(1)}_{(i,j,k)}(t)  \quad\quad
 * \forall t, \, \forall {(i,j,k)} \in \mathcal{G}_I, \\
 * \frac{d}{dt} Y^{(2)}_{(i,j,k)}(t) &=&  \big( Y^{(2)}_{(i-1,j,k)}(t) - 2 \cdot Y^{(2)}_{(i,j,k)}(t)
 *               + Y^{(2)}_{(i+1,j,k)} (t)\big) \cdot \mu  /(\rho \cdot  h_0^2)   \\
 * &&
 *               \quad  \quad \quad  \quad  + \big(  Y^{(2)}_{(i,j-1,k)}(t) - 2 \cdot Y^{(2)}_{(i,j,k)}(t)
 *               + Y^{(2)}_{(i,j+1,k)}(t) \big) \cdot  \mu  / (\rho \cdot h_1^2)   \\
 * &&            \quad \quad  \quad \quad  + \big(  Y^{(2)}_{(i,j,k-1)}(t) - 2 \cdot Y^{(2)}_{(i,j,k)}(t)
 *               + Y^{(2)}_{(i,j,k+1)}(t) \big)  \cdot (2 \cdot \mu + \lambda)   / (\rho \cdot h_2^2)  \\
 * &&            \quad \quad  \quad \quad  + \big(  Y^{(2)}_{(i+1,j,k)}(t) -  Y^{(2)}_{(i+1,j,k-1)}(t)
 *               - Y^{(2)}_{(i,j,k)}(t) +  Y^{(2)}_{(i,j,k-1)}(t)  \big)  \cdot (\mu + \lambda)   / (\rho \cdot h_0 \cdot h_1)  \\
 * &&            \quad \quad  \quad \quad  + \big(  Y^{(2)}_{(i,j+1,k)}(t) -  Y^{(2)}_{(i,j+1,k-1)}
 *               - Y^{(2)}_{(i,j,k)}(t) +  Y^{(2)}_{(i,j,k-1)}(t)  \big)  \cdot (\mu + \lambda)  / (\rho \cdot h_0 \cdot h_2)  \\
 * &&  \quad \quad \quad\quad + 1 / \rho  \cdot F^{(2)}_{(i,j,k)}(t)  \quad\quad
 * \forall t, \, \forall {(i,j,k)} \in \mathcal{G}_I, \\
 * Y^{(p)}_{(i,j,k)}(t) &=&  G^{(p;l)}_{(i,j,k)}
 *               \quad   \quad \mbox{for all } p \in \{0,1,\ldots,d-1\}, \forall l, \, \forall {(i,j,k)} \in \mathcal{G}_B, \\
 *            Y^{(p;0)}_{(i,j,k)}(t_S) & =&  \widetilde{Y}^{(p)}_{(i,j,k)}\quad \quad
 * \mbox{for all } p \in \{0,1,\ldots,d-1\}, \forall {(i,j,k)} \in \mathcal{G}, \\
 *            \frac{d}{dt}  Y^{(p)}_{(i,j,k)}(t)|_{t_S} & =& \widehat{Y}^{(p)}_{(i,j,k)}\quad \quad
 * \mbox{for all } p \in \{0,1,\ldots,d-1\}, \forall {(i,j,k)} \in \mathcal{G}.
 * \f}
 * Rewrite system in more compact form:
 *
 * Find \f$Y^{(\cdot)}: [t_S;t_E] \to R^{d \cdot N}\f$ such that
 * \f{eqnarray*}{
 * \frac{d}{dt} Y(t)|_{(i,j,k)} &=&  [S  Y]|_{(i,j,k)}
 *               \quad   \quad  \forall {(i,j,k)} \in \mathcal{G}_I, \mbox{ in } [t_S;t_E], \\
 *          Y(t)|_{(i,j,k)} &=&  G|_{(i,j,k)}
 *               \quad   \quad  \forall {(i,j,k)} \in \mathcal{G}_B \mbox{ in } [t_S;t_E], \\
 *            Y_{(i,j,k)}(t_S) & =&  \widetilde{Y}_{(i,j,k)}\quad \quad
 * \forall {(i,j,k)} \in \mathcal{G}, \\
 *          \frac{d}{dt}  Y_{(i,j,k)}(t)|_{t_S} & =&  \widehat{Y}_{(i,j,k)} \quad \quad
 * \forall {(i,j,k)} \in \mathcal{G}. \\
 * \f}
 * Transform system of ordinary differential equations of second order wrt. time into a system of ordinary differential equations.
 * Introduce new vector-valued functions for each component and the total vector: \f$ U^{(p)}: [t_S;t_E] \to R^{N}, t \mapsto U^{(p)}_{(i,j,k)}(t) :=Y^{(p)}_{(i,j,k)}(t),\f$
\f$ V^{(p)}: [t_S;t_E] \to R^{N}: t \mapsto V^{(p)}_{(i,j,k)}(t) :=d/dt Y^{(p)}_{(i,j,k)}(t)\f$ for all \f$ p \f$,
\f$ U: [t_S;t_E] \to R^{d\cdot N}, t \mapsto U(t) :=[U^{(0)}(t);\ldots U^{(p-1)}(t)], \f$
\f$ V: [t_S;t_E] \to R^{d\cdot N}, t \mapsto V(t) :=[V^{(0)}(t);\ldots V^{(p-1)}(t)]. \f$
 *
 * Find \f$ [U;V]: [t_S;t_E] \to R^{2 \cdot d \cdot N}\f$ such that
 * \f{eqnarray*}{
 * (\frac{d}{dt} U|)_{(i,j,k)}(t) &=&  V_{(i,j,k)}(t)
 *               \quad   \quad  \forall {(i,j,k)} \in \mathcal{G}  \mbox{ in } [t_S;t_E], \\
 * (\frac{d}{dt} V|)_{(i,j,k)}(t) &=& - [T U]_{(i,j,k)}(t)
 *               \quad   \quad  \forall {(i,j,k)} \in \mathcal{G}_I  \mbox{ in } [t_S;t_E], \\
 * Y_{(i,j,k)} &=&  G_{(i,j,k)}
 *               \quad   \quad  \forall {(i,j,k)} \in \mathcal{G}_B  \mbox{ in } [t_S;t_E], \\
 *            Y_{(i,j,k)}|_{t_S} & =&  \widetilde{Y}_{(i,j,k)}\quad \quad
 * \forall {(i,j,k)} \in \mathcal{G}, \\
 *          \frac{d}{dt}  Y_{(i,j,k)}|_{t_S} & =&  \widehat{Y}_{(i,j,k)} \quad \quad
 * \forall {(i,j,k)} \in \mathcal{G} \\
 * \f}
 * with matrix \f$ T \in R^{d \cdot N, d \cdot N} \f$. Rewrite system:
 * Introduce new vectors: \f$ U^{(p,l)}_{(i,j,k)}:=Y^{(p)}_{(i,j,k)}(t_l) \in R^{d \cdot N} \f$ for each time step \f$ t_l \f$ and for all \f$ p \f$.
 * Apply explicit Euler scheme (forward difference quotient).
 *
 * Find \f$[U^{(\cdot,l+1)}; V^{(\cdot,l+1)}]: R^{2 \cdot d\cdot N} \to R^{2 \cdot d \cdot N}\f$ such that
 * \f{eqnarray*}{
 * U^{(\cdot, l+1)} &=& U^{(\cdot ,l)}  + \tau \cdot V^{(\cdot, l)}
 *               \quad   \quad \forall l\ge1, \,  \forall {(i,j,k)} \in \mathcal{G}, \\
 * V^{(\cdot ,l+1)} &=& V^{(\cdot ,l)}  - \tau \cdot T U^{(\cdot, l)}
 *               \quad   \quad  \forall l\ge1, \, \forall {(i,j,k)} \in \mathcal{G}_I, \\
 * Y^{(\cdot, l+1)}_{(i,j,k)} &=&  G_{(i,j,k)}
 *               \quad   \quad \forall l\ge1, \,  \forall {(i,j,k)} \in \mathcal{G}_B, \\
 *            Y^{(\cdot, 0)}_{(i,j,k)} & =&  \widetilde{Y}_{(i,j,k)}\quad \quad
 * \forall {(i,j,k)} \in \mathcal{G}, \\
 *            Y^{(\cdot, 1)}_{(i,j,k)}& =&  Y^{(\cdot, 0)}_{(i,j,k)} + \tau \cdot \widehat{Y}_{(i,j,k)} \quad \quad
 * \forall {(i,j,k)} \in \mathcal{G}. \\
* \f}.
*/
template<typename CT, std::size_t DIM>
class LameNavierEqnFDMFirstOrder : public ScaFES::Problem<LameNavierEqnFDMFirstOrder<CT,DIM>, CT, DIM> {
  public:

    /** Constant mu. Material parameter. */
    const double ELAST_MODULE = 1.0e6;

    /** Constant NU. Material parameter. */
    const double NU = 0.4999;

    /** Constant lambda. Material parameter. */
    const double LAMBDA = NU * ELAST_MODULE /( (1.0 + NU) * (1.0 - 2.0 * NU) );

    /** Constant mu. Material parameter. */
    const double MU = ELAST_MODULE /( 2.0 * (1.0 + NU) );

    /** Constant rho. Parameter. */
    const double RHO = 1.0;

    /** Constant FREQUENCY. Parameter. */
    const double FREQUENCY = 1.0;

    /** Constant OMEGA. Parameter. */
    const double OMEGA = 2.0 * M_PI * FREQUENCY;

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
    LameNavierEqnFDMFirstOrder(ScaFES::Parameters const& params,
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
        : ScaFES::Problem<LameNavierEqnFDMFirstOrder<CT, DIM>, CT, DIM>(params, gg, useLeapfrog,
                                                        nameDatafield, stencilWidth,
                                                        isKnownDf, nLayers,
                                                        defaultValue, writeToFile,
                                                        computeError, geomparamsInit)
        { }

    /** Evaluates all fields at one given global inner grid node.
     *  @param vNew Set of all fields (return value), \f$ F^{(p;0)}_{(i,j,k)} \f$, \f$ G^{(p;0)}_{(i,j,k)} \f$, \f$ Y^{(p;0)}_{(i,j,k)} \f$.
     *  @param idxNode Index \f$(i,j,k)\f$ of given inner grid node \f$x_{(i,j,k)}\f$.
     *  @param timestep Given first time step.
     */
    void evalInner(std::vector< ScaFES::DataField<CT, DIM> >& vNew,
                   ScaFES::Ntuple<int,DIM> const& idxNode,
                   int const& timestep) {
        double t = this->time(timestep);
        vNew[0](idxNode) = 0.0;
        vNew[1](idxNode) = 0.0;
        vNew[2](idxNode) = 0.0;
        vNew[3](idxNode) = 0.0;
        vNew[4](idxNode) = 0.0;
        vNew[5](idxNode) = 0.0;
        vNew[6](idxNode) = 0.0;
        vNew[7](idxNode) = 0.0;
        vNew[8](idxNode) = sin(OMEGA * t);
        vNew[6](idxNode) = 0.0;
        vNew[7](idxNode) = 0.0;
        vNew[8](idxNode) = OMEGA * cos(OMEGA * t);
    }

    /** Evaluates all fields at one given global border grid node.
     *  @param vNew Set of all fields (return value), \f$ F^{(p;0)}_{(i,j,k)} \f$, \f$ G^{(p;0)}_{(i,j,k)} \f$, \f$ Y^{(p;0)}_{(i,j,k)} \f$.
     *  @param idxNode Index \f$(i,j,k)\f$ of current inner grid node \f$ x_{(i,j,k)}\f$.
     *  @param timestep Given first time step.
     */
    void evalBorder(std::vector< ScaFES::DataField<CT, DIM> >& vNew,
                    ScaFES::Ntuple<int,DIM> const& idxNode,
                    int const& timestep) {
            this->evalInner(vNew, idxNode, timestep);
    }

    /** Initializes all unknown fields at one given global inner grid node.
     *  @param vNew Set of all unknown fields (return value), \f$ Y^{(p;0)}_{(i,j,k)} \f$.
     *  @param idxNode Index \f$(i,j,k)\f$ of current inner grid node \f$ x_{(i,j,k)}\f$.
     *  @param timestep Given time step 0.
     */
    template<typename TT>
    void initInner(std::vector< ScaFES::DataField<TT, DIM> >& vNew,
                   std::vector<TT> const& /*vOld*/,
                   ScaFES::Ntuple<int,DIM> const& idxNode,
                   int const& timestep) {
        double t = this->time(timestep);
        vNew[0](idxNode) = 0.0;
        vNew[1](idxNode) = 0.0;
        vNew[2](idxNode) = sin(OMEGA * t);
        vNew[3](idxNode) = 0.0;
        vNew[4](idxNode) = 0.0;
        vNew[5](idxNode) = OMEGA * cos(OMEGA * t);
    }

    /** Initializes all unknown fields at one given global border grid node.
     *  @param vNew Set of all unknown fields (return value), \f$ Y^{(p;0)}_{(i,j,k)} \f$.
     *  @param vOld Set of all given fields, \f$ Y^{(p;-1}_{(i,j,k)} \f$.
     *  @param idxNode Index \f$(i,j,k)\f$ of current border grid node \f$x_{(i,j,k)}\f$.
     *  @param timestep Given time step 0.
     */
    template<typename TT>
    void initBorder(std::vector< ScaFES::DataField<TT, DIM> >& vNew,
                    std::vector<TT> const& vOld,
                    ScaFES::Ntuple<int,DIM> const& idxNode,
                    int const& timestep) {
        this->template initInner<TT>(vNew, vOld, idxNode, timestep);
    }

    /** Updates all unknown fields at one given global inner grid node.
     *  @param vNew Set of all unknown fields (return value), \f$ [U^{(p;l+1)}_{(i,j,k)}; V^{(p;l+1)}_{(i,j,k)}]\f$.
     *  @param vOld Set of all given fields, \f$ [U^{(p;l}_{(i,j,k)}; V^{(p;l}_{(i,j,k)}] \f$.
     *  @param idxNode Index \f$(i,j,k)\f$ of current inner grid node \f$x_{(i,j,k)}\f$.
     */
    template<typename TT>
    void updateInner(std::vector<ScaFES::DataField<TT,DIM>>& vNew,
                     std::vector<ScaFES::DataField<TT,DIM>> const& vOld,
                     ScaFES::Ntuple<int,DIM> const& idxNode,
                     int const& /*timestep*/) {
        vNew[0](idxNode) = vOld[0](idxNode) + this->tau() * vOld[3](idxNode);
        vNew[1](idxNode) = vOld[1](idxNode) + this->tau() * vOld[4](idxNode);
        vNew[2](idxNode) = vOld[2](idxNode) + this->tau() * vOld[5](idxNode);

        vNew[3](idxNode) += vOld[3](idxNode);
        for (std::size_t pp = 0; pp < DIM; ++pp) {
            vNew[0](idxNode) -= MU * this->tau() *
                  (
                        2.0 * vOld[0](idxNode)
                         - vOld[0](this->connect(idxNode, 2*pp))
                         - vOld[0](this->connect(idxNode, 2*pp+1))
                  )
                  / (this->gridsize(pp) * this->gridsize(pp) * RHO);
        }
        vNew[3](idxNode) -= (MU + LAMBDA) * this->tau() *
                      (
                         vOld[0](this->connect(idxNode, 1))
                       - vOld[0](idxNode)
                       - vOld[0](idxNode)
                       + vOld[0](this->connect(idxNode, 0))
                      )
                      / (this->gridsize(0) * this->gridsize(0));
        vNew[3](idxNode) -= (MU + LAMBDA) * this->tau() *
                      (
                         vOld[0](this->connect(idxNode, 3))
                       - vOld[0](this->connect(this->connect(idxNode, 3), 0))
                       - vOld[0](idxNode)
                       + vOld[0](this->connect(idxNode, 0))
                      )
                      / (this->gridsize(0) * this->gridsize(1));
        vNew[3](idxNode) -= (MU + LAMBDA) * this->tau() *
                      (
                         vOld[0](this->connect(idxNode, 5))
                       - vOld[0](this->connect(this->connect(idxNode, 5), 0))
                       - vOld[0](idxNode)
                       + vOld[0](this->connect(idxNode, 0))
                      )
                      / (this->gridsize(0) * this->gridsize(2));
        vNew[3](idxNode) -= this->tau() * this->knownDf(0, idxNode);

        vNew[4](idxNode) = vOld[4](idxNode);
        for (std::size_t pp = 0; pp < DIM; ++pp) {
            vNew[1](idxNode) -= MU * this->tau() *
                  ( 2.0 * vOld[1](idxNode)
                          - vOld[1](this->connect(idxNode, 2*pp))
                          - vOld[1](this->connect(idxNode, 2*pp+1))
                  )
                  / (this->gridsize(pp) * this->gridsize(pp));
        }
        vNew[4](idxNode) -= (MU + LAMBDA) * this->tau() *
                      (
                         vOld[1](this->connect(idxNode, 3))
                       - vOld[1](idxNode)
                       - vOld[1](idxNode)
                       + vOld[1](this->connect(idxNode, 2))
                      )
                      / (this->gridsize(1) * this->gridsize(1));
        vNew[4](idxNode) -= (MU + LAMBDA) * this->tau() *
                      (
                         vOld[1](this->connect(idxNode, 1))
                       - vOld[1](this->connect(this->connect(idxNode, 1), 2))
                       - vOld[1](idxNode)
                       + vOld[1](this->connect(idxNode, 2))
                      )
                      / (this->gridsize(1) * this->gridsize(0));
        vNew[4](idxNode) -= (MU + LAMBDA) * this->tau() *
                      (
                         vOld[1](this->connect(idxNode, 5))
                       - vOld[1](this->connect(this->connect(idxNode, 5), 2))
                       - vOld[1](idxNode)
                       + vOld[1](this->connect(idxNode, 2))
                      )
                      / (this->gridsize(1) * this->gridsize(2));
        vNew[4](idxNode) -= this->tau() * this->knownDf(1, idxNode);


        vNew[5](idxNode) = vOld[5](idxNode);
        for (std::size_t pp = 0; pp < DIM; ++pp) {
            vNew[0](idxNode) -= MU * this->tau() *
                  ( 2.0 * vOld[2](idxNode)
                          - vOld[2](this->connect(idxNode, 2*pp))
                          - vOld[2](this->connect(idxNode, 2*pp+1))
                  )
                  / (this->gridsize(pp) * this->gridsize(pp));
        }
        vNew[5](idxNode) -= (MU + LAMBDA) * this->tau() *
                      (
                         vOld[2](this->connect(idxNode, 1))
                       - vOld[2](idxNode)
                       - vOld[2](idxNode)
                       + vOld[2](this->connect(idxNode, 4))
                      )
                      / (this->gridsize(2) * this->gridsize(2));
        vNew[5](idxNode) -= (MU + LAMBDA) * this->tau() *
                      (
                         vOld[2](this->connect(idxNode, 1))
                       - vOld[2](this->connect(this->connect(idxNode, 1), 4))
                       - vOld[2](idxNode)
                       + vOld[2](this->connect(idxNode, 4))
                      )
                      / (this->gridsize(2) * this->gridsize(1));
        vNew[5](idxNode) -= (MU + LAMBDA) * this->tau() *
                      (
                         vOld[2](this->connect(idxNode, 3))
                       - vOld[2](this->connect(this->connect(idxNode, 3), 4))
                       - vOld[2](idxNode)
                       + vOld[2](this->connect(idxNode, 4))
                      )
                      / (this->gridsize(2) * this->gridsize(2));
        vNew[5](idxNode) -= this->tau() * this->knownDf(2, idxNode);
    }

    /** Updates all unknown fields at one given global border grid node.
     *  @param vNew Set of all unknown fields (return value), \f$ Y^{(p;l+1)}_{(i,j,k)} \f$.
     *  @param vOld Set of all given fields.
     *  @param idxNode Index \f$(i,j,k)\f$ of current border grid node \f$x_{(i,j,k)}\f$.
     */
    template<typename TT>
    void updateBorder(std::vector<ScaFES::DataField<TT,DIM>>& vNew,
                      std::vector<ScaFES::DataField<TT,DIM>>const& vOld,
                      ScaFES::Ntuple<int,DIM> const& idxNode,
                      int const& /*timestep*/) {
        vNew[0](idxNode) = vOld[0](idxNode) + this->tau() * vOld[3](idxNode);
        vNew[1](idxNode) = vOld[1](idxNode) + this->tau() * vOld[4](idxNode);
        vNew[2](idxNode) = vOld[2](idxNode) + this->tau() * vOld[5](idxNode);
        vNew[3](idxNode) = this->knownDf(3, idxNode);
        vNew[4](idxNode) = this->knownDf(4, idxNode);
        vNew[5](idxNode) = this->knownDf(5, idxNode);
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
