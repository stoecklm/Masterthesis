/* ScaFES
 * Copyright (c) 2017, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file PennesBioheatEqnFDM.hpp
 *
 *  @brief Implementation of n-dimensional Pennes bioheat equation problem on unit hybercube.
 */

#include "ScaFES.hpp"
#include "analyticalSolutions.hpp"
#include "boundaryConditions.hpp"

/*******************************************************************************
 ******************************************************************************/
/**
 * \class PennesBioheatEqnFDM
 *  @brief Class for discretized Pennes bioheat equation problem.
*/
template<typename CT, std::size_t DIM>
class PennesBioheatEqnFDM : public ScaFES::Problem<PennesBioheatEqnFDM<CT,DIM>, CT, DIM> {

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

    /** constant rho_blood. Material parameter (density) for blood. */
    const double RHO_BLOOD = 1.0; /* kg/m^3 */

    /** constant c_blood. Material parameter (specific heat capacity)
     * for blood. */
    const double C_BLOOD = 1.0; /* J/(kg K) */

    /** constant w. Material parameter (perfusion) for the brain. */
    const double W = 1.0; /* */

    /** constant T_blood. Parameter (temperature) for blood. */
    const double T_BLOOD = 1.0; /* K */

    /** constant q_M_dot. Material parameter (metabolism heat rate)
     *  for the brain. */
    const double Q_M_DOT = 1.0; /* W/(m^2) */

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
                        std::vector<CT> const& geomparamsInit = std::vector<CT>(),
                        int const& eqnDegree_ = 1,
                        int const& boundaryCond_ = 1)
        : ScaFES::Problem<PennesBioheatEqnFDM<CT, DIM>, CT, DIM>(params, gg, useLeapfrog,
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
            vNew[0](idxNode) += RHO_BLOOD * C_BLOOD * W * consFunc<CT,DIM>(x, t);
        } else if (eq == linear) {
            vNew[0](idxNode) = RHO * C * linFuncTimeDerivative<CT,DIM>(x);
            vNew[0](idxNode) -= LAMBDA * linFuncSumOfSpaceDerivatives2ndOrder<CT,DIM>(x, t);
            vNew[0](idxNode) += RHO_BLOOD * C_BLOOD * W * linFunc<CT,DIM>(x, t);
        } else if (eq == quadratic) {
            vNew[0](idxNode) = RHO * C * quadFuncTimeDerivative<CT,DIM>(x);
            vNew[0](idxNode) -= LAMBDA * quadFuncSumOfSpaceDerivatives2ndOrder<CT,DIM>(x, t);
            vNew[0](idxNode) += RHO_BLOOD * C_BLOOD * W * quadFunc<CT,DIM>(x, t);
        } else if (eq == cubic) {
            vNew[0](idxNode) = RHO * C * cubicFuncTimeDerivative<CT,DIM>(x);
            vNew[0](idxNode) -= LAMBDA * cubicFuncSumOfSpaceDerivatives2ndOrder<CT,DIM>(x, t);
            vNew[0](idxNode) += RHO_BLOOD * C_BLOOD * W * cubicFunc<CT,DIM>(x, t);
        } else {
            std::cerr << "ERROR in evalInner: Degree of equation does not have a valid value." << std::endl;
            vNew[0](idxNode) = -1.0;
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
            std::cerr << "ERROR in evalInner: Degree of equation does not have a valid value." << std::endl;
            vNew[2](idxNode) = -1.0;
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
            vNew[0](idxNode) = 0.0;
        } else {
            std::cerr << "ERROR in evalBorder: Type of boundary condition does not have a valid value." << std::endl;
            vNew[0](idxNode) = -1.0;
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
            std::cerr << "ERROR in evalBorder: Degree of equation does not have a valid value." << std::endl;
            vNew[1](idxNode) = -1.0;
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
            std::cerr << "ERROR in evalBorder: Degree of equation does not have a valid value." << std::endl;
            vNew[2](idxNode) = -1.0;
        }
    }

    /** Initializes all unknown fields at one given global inner grid node.
     *  @param vNew Set of all unknown fields (return value).
     *  @param idxNode Index of given grid node.
     *  @param timestep Given time step.
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
            std::cerr << "ERROR in initInner: Degree of equation does not have a valid value." << std::endl;
            vNew[0](idxNode) = -1.0;
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
        vNew[0](idxNode) = vOld[0](idxNode);
        for (std::size_t pp = 0; pp < DIM; ++pp) {
            vNew[0](idxNode) += this->tau() * (LAMBDA/(RHO*C)) * (
                                vOld[0](this->connect(idxNode, 2*pp))
                                + vOld[0](this->connect(idxNode, 2*pp+1))
                                - 2.0 * vOld[0](idxNode) )
                                / (this->gridsize(pp) * this->gridsize(pp));
        }
        vNew[0](idxNode) -= this->tau() * ((RHO_BLOOD*C_BLOOD)/(RHO*C))
                                        * W * vOld[0](idxNode);
        vNew[0](idxNode) += this->tau() * (1.0/(RHO*C)) * this->knownDf(0, idxNode);
    }

    /** Updates all unknown fields at one given global border grid node.
     *  @param vNew Set of all unknown fields at new time step (return value).
     *  @param vOld Set of all given fields.
     *  @param idxNode Index of given grid node.
     *  @param timestep Given time step.
     */
    template<typename TT>
    void updateBorder(std::vector<ScaFES::DataField<TT,DIM>>& vNew,
                      std::vector<ScaFES::DataField<TT,DIM>>const& vOld,
                      ScaFES::Ntuple<int,DIM> const& idxNode,
                      int const& /*timestep*/) {
        int bc = this->boundaryCond;

        /* Vector for U. */
        if (bc == dirichlet) {
            vNew[0](idxNode) = this->knownDf(1, idxNode);
        } else {
            std::cerr << "ERROR: Type of boundary condition does not have a valid value." << std::endl;
            vNew[0](idxNode) = -1.0;
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

    /** Updates all unknown fields at one given global border grid node
     *  by using either Neumann or Cauchy boundary condition.
     *  @param lhsNeighbour Value of lhs neighbour of current node in specified direction.
     *  @param rhsNeighbour Value of rhs neighbour of current node in specified direction.
     *  @param vOld Set of all given fields.
     *  @param idxNode Index of given grid node.
     *  @param pp Specified direction (current element of dimension).
     */
    template<typename TT>
    void updateBorderHelper(TT& lhsNeighbour, TT& rhsNeighbour,
                            std::vector<ScaFES::DataField<TT,DIM>>const& vOld,
                            ScaFES::Ntuple<int,DIM> const& idxNode,
                            std::size_t const& pp) {
        ScaFES::Ntuple<int,DIM> const& lhsIdxNode = this->connect(idxNode, 2*pp);
        ScaFES::Ntuple<int,DIM> const& rhsIdxNode = this->connect(idxNode, 2*pp+1);
        int bc = this->boundaryCond;

        if (this->insideGrid(lhsIdxNode) == true && this->insideGrid(rhsIdxNode) == true) {
            lhsNeighbour = vOld[0](lhsIdxNode);
            rhsNeighbour = vOld[0](rhsIdxNode);
        } else if (this->insideGrid(rhsIdxNode) == false) {
            lhsNeighbour = vOld[0](lhsIdxNode);
            if (bc == neumann) {
                rhsNeighbour = boundaryCondition2ndTypeIndexPlusOne(lhsNeighbour,
                                   this->gridsize(pp), LAMBDA, Q_DOT);
            } else { /*(bc == cauchy)*/
                rhsNeighbour = boundaryCondition3rdTypeIndexPlusOne(lhsNeighbour, vOld[0](idxNode), T_INF,
                                   this->gridsize(pp), ALPHA, LAMBDA);
            }
        } else { /*(this->insideGrid(lhsIdxNode) == false)*/
            rhsNeighbour = vOld[0](rhsIdxNode);
            if (bc == neumann) {
                lhsNeighbour = boundaryCondition2ndTypeIndexMinusOne(rhsNeighbour,
                                   this->gridsize(pp), LAMBDA, Q_DOT);
            } else { /*(bc == cauchy)*/
                lhsNeighbour = boundaryCondition3rdTypeIndexMinusOne(rhsNeighbour, vOld[0](idxNode), T_INF,
                                   this->gridsize(pp), ALPHA, LAMBDA);
            }
        }
    }

    /** Checks if given node is inside grid.
     *  @param idxNode Index of given grid node.
     */
    bool insideGrid(ScaFES::Ntuple<int,DIM> const& idxNode) {
        for (std::size_t pp = 0; pp < DIM; pp++) {
            if (idxNode.elem(pp) < 0) {
                return false;
            }
        }
        return true;
    }
};
