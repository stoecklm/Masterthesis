#ifndef PROBLEM_MAXWELLEQNS_FDTD_HPP_
#define PROBLEM_MAXWELLEQNS_FDTD_HPP_

#include "ScaFES.hpp"

/*******************************************************************************
*******************************************************************************/
/**
 * @file Problem_MaxwellEqnsFDTD.hpp
 *
 * @brief Maxwell equations
 *
 * \ingroup
 *
 * @section Maxwell_equations
 *
 * First, consider the Maxwell equations in differential form:
 * -----------------------------------------------------------
 *  Let be
 *    * t_0 \in R_{>0} starting time,
 *    * T \in R_{>0} ending time
 *    * rho: \OmegaBar x [t_0;T] --> \R^3 (Ladungsdichte),
 *    * j_l: \OmegaBar x [t_0;T] --> \R^3 (Leitungsstromdichte),
 *    * D: \OmegaBar x [t_0;T] --> \R^3 (elektrische Flussdichte),
 *    * B: \OmegaBar x [t_0;T] --> \R^3 (magnetische Flussdichte),
 *    * E: \OmegaBar x [t_0;T] --> \R^3 (electric field),
 *    * H: \OmegaBar x [t_0;T] --> \R^3 (magnetic field),
 *  then the Maxwell equations are given by:
 *       div(D)                 = \rho     in \Omega x [t_0;T],
 *       div(B)                 = 0        in \Omega x [t_0;T],
 *       rot(E) + \partial_t(B) = 0        in \Omega x [t_0;T],
 *       rot(H) + \partial_t(D) = j_l      in \Omega x [t_0;T],
 *
 *
 * Next, consider the material equations:
 * --------------------------------------
 *  Let be
 *    * \epsilon_0 \in \R (permittivity of the vacuum),
 *    * \epsilon_r \in \R (relative permittivity),
 *    * \mu_0 \in \R (permeability of the vacuum),
 *    * \mu_r \in \R (relative permeability),
 *  Furthermore, let be the material be linear, homogeneous, isotropic,
 *  then the permittivities and the permeabilities are scalar and constant:
 *    * \epsilon_0  = const. \in \R,
 *    * \epsilon_r  = const. \in \R,
 *    * \mu_0       = const. \in \R,
 *    * \mu_r       = const. \in \R.
 *  Thus, define:
 *    * \epsilon := \epsilon_0 \epsilon_r (permittivity=dielektrische Leitf.),
 *    * \mu      := \mu_0 \mu_r (permeability=magnetische Leitfaehigkeit).
 *  and the material equations are given by:
 *       D = \epsilon E    in \Omega x [t_0;T],
 *       B = \mu H         in \Omega x [t_0;T],
 *
 *
 *
 * =============================================================================
 * Now, consider the following initial value problem:
 * (with pure Dirichlet boundary conditions)
 * -------------------------------------------------
 *  Given:
 *    * t_0 \in R_{>0} starting time,
 *    * T \in R_{>0} ending time,
 *    * domain \Omega \subset \R^d3 with boundary \partial\Omega
 *    * \epsilon \in R,
 *    * \mu \in R,
 *    * \sigma \in R (eletric Leitfaehigkeit),
 *    * rho: \OmegaBar x [t_0;T] --> \R^3 (Ladungsdichte),
 *    * j_l: \OmegaBar x [t_0;T] --> R^3,
 *    * E_0: \OmegaBar --> R^3,
 *    * H_0: \OmegaBar --> R^3,
 *    * E_D: \partial\Omega x (t_0;T) --> R^3,
 *    * H_D: \partial\Omega x (t_0;T) --> R^3,
 *    * g: \partial\Omega x [t_0;T] --> R.
 *  Seek:
 *    * E: \OmegaBar x [t_0;T] --> \R^3 (electric field),
 *    * H: \OmegaBar x [t_0;T] --> \R^3 (magnetic field),
 *    such that
 *       \epsilon \partial_t(E) = rot(H) - \sigma E  in \Omega x [t_0;T],
 *       \mu \partial_t(H)      = -rot(E)            in \Omega x [t_0;T],
 *       \epsilon div(E)        = \rho               in \Omega x [t_0;T],
 *       \mu div(H)             = 0                  in \Omega x [t_0;T],
 *                      E(.,t0) = E_0                in \OmegaBar,
 *                      H(.,t0) = H_0                in \OmegaBar,
 *                      E       = E_D          on \partial\Omega x (t_0;T).
 *                      H       = H_D          on \partial\Omega x (t_0;T).
 *
 * Aim:
 * ----
 * Solve the initial value problem (i.v.p.) using the finite difference
 * method at time domain (FDTD).
 */
template<class CT = double>
class Problem_MaxwellEqnsFDTD
    : public ScaFES::Problem<Problem_MaxwellEqnsFDTD<CT>, CT, 3>

{
    public:
        /*----------------------------------------------------------------------
        | TYPE DEFINITIONS.
        ----------------------------------------------------------------------*/
        typedef ScaFES::Ntuple<int, 3> Index;

        /*----------------------------------------------------------------------
        | ENUMERATIONS.
        ----------------------------------------------------------------------*/
        /** Geometry parameters. */
        enum geomparams { COORDX=0, COORDY=1, COORDZ=2,
                          RADIUSX=3, RADIUSY=4, RADIUSZ=5 };

        /** Data field names. */
        enum dfnames { EPS_R=0, MUE_R=1, SIGMA=2,
                       BX=3, BY=4, BZ=5, DX=6, DY=7, DZ=8,
                       EX=9, EY=10, EZ=11, HX=12, HY=13, HZ=14 };

        /*----------------------------------------------------------------------
        | LIFE CYCLE METHODS.
        ----------------------------------------------------------------------*/
        /** Creates own constructor. */
        Problem_MaxwellEqnsFDTD(ScaFES::Parameters const& cl,
                                ScaFES::GridGlobal<3> const& gg);

        /*----------------------------------------------------------------------
        | WORK METHODS.
        ----------------------------------------------------------------------*/
        /** Initialises the source in -x and y-direction. */
        template<typename TT>
        void funcZero(TT& fx, ScaFES::Ntuple<TT,3> const& x, TT const& t);

        /** Initialises the source in z-direction. */
        template<typename TT>
        void funcSourceZ(TT& fx, ScaFES::Ntuple<TT,3> const& x, TT const& t);

        /** Computes new iteration at interior grid nodes. */
        template<typename TT>
        void updateInner(std::vector< ScaFES::DataField<TT, 3> >& vectNew,
                         std::vector< ScaFES::DataField<TT, 3> > const& vectOld,
                         Index const& idxNode,
                         int const& timestep);

        /** Computes new iteration at boundary grid nodes. */
        template<typename TT>
        void updateBorder(std::vector< ScaFES::DataField<TT, 3> >& vectNew,
                          std::vector< ScaFES::DataField<TT, 3> > const& vectOld,
                          Index const& idxNode,
                          int const& timestep);

    private:
        /*----------------------------------------------------------------------
        | STATIC CONSTANTS.
        ----------------------------------------------------------------------*/
        /** EpsR_0. */
        static const double EPSILON_0;

        /** Constant permeability 'mu_0'. */
        static const double MUE_0;

        /** Light speed = const. */
        static const double SPEED_LIGHT;

        /** Number of perfectly matched layers. */
        static const int N_PML;
}; // End of class. //

/*******************************************************************************
 ******************************************************************************/
template<class CT>
double const Problem_MaxwellEqnsFDTD<CT>::EPSILON_0 = 8.854187817620389850536563031e-12;
/*----------------------------------------------------------------------------*/
template<class CT>
double const Problem_MaxwellEqnsFDTD<CT>::MUE_0 = 1.256637061435917295385057353e-6;
/*----------------------------------------------------------------------------*/
template<class CT>
double const Problem_MaxwellEqnsFDTD<CT>::SPEED_LIGHT = 299792458.0;
/*----------------------------------------------------------------------------*/
template<class CT>
int const Problem_MaxwellEqnsFDTD<CT>::N_PML = 10;

/*******************************************************************************
 ******************************************************************************/
template<class CT>
inline Problem_MaxwellEqnsFDTD<CT>::Problem_MaxwellEqnsFDTD(ScaFES::Parameters const& cl,
                                               ScaFES::GridGlobal<3> const& gg)
    : ScaFES::Problem<Problem_MaxwellEqnsFDTD, double, 3>(cl, gg)
{
    this->addDataField("SourceX", 0, *this, &Problem_MaxwellEqnsFDTD<CT>::funcZero<double>, true);
    this->addDataField("SourceY", 0, *this, &Problem_MaxwellEqnsFDTD<CT>::funcZero<double>, true);
    this->addDataField("SourceZ", 0, *this, &Problem_MaxwellEqnsFDTD<CT>::funcSourceZ<double>, true);
    this->addDataField("EpsR", 1, *this, &Problem_MaxwellEqnsFDTD<CT>::funcZero<double>, false);
    this->addDataField("MueR", 1, *this, &Problem_MaxwellEqnsFDTD<CT>::funcZero<double>, false);
    this->addDataField("Sigma", 1, *this, &Problem_MaxwellEqnsFDTD<CT>::funcZero<double>, false);
    this->addDataField("Bx", 1, *this, &Problem_MaxwellEqnsFDTD<CT>::funcZero<double>, false);
    this->addDataField("By", 1, *this, &Problem_MaxwellEqnsFDTD<CT>::funcZero<double>, false);
    this->addDataField("Bz", 1, *this, &Problem_MaxwellEqnsFDTD<CT>::funcZero<double>, false);
    this->addDataField("Dx", 1, *this, &Problem_MaxwellEqnsFDTD<CT>::funcZero<double>, false);
    this->addDataField("Dy", 1, *this, &Problem_MaxwellEqnsFDTD<CT>::funcZero<double>, false);
    this->addDataField("Dz", 1, *this, &Problem_MaxwellEqnsFDTD<CT>::funcZero<double>, false);
    this->addDataField("Ex", 1, *this, &Problem_MaxwellEqnsFDTD<CT>::funcZero<double>, false);
    this->addDataField("Ey", 1, *this, &Problem_MaxwellEqnsFDTD<CT>::funcZero<double>, false);
    this->addDataField("Ez", 1, *this, &Problem_MaxwellEqnsFDTD<CT>::funcZero<double>, false);
    this->addDataField("Hx", 1, *this, &Problem_MaxwellEqnsFDTD<CT>::funcZero<double>, false);
    this->addDataField("Hy", 1, *this, &Problem_MaxwellEqnsFDTD<CT>::funcZero<double>, false);
    this->addDataField("Hz", 1, *this, &Problem_MaxwellEqnsFDTD<CT>::funcZero<double>, false);
}
/*----------------------------------------------------------------------------*/
template<class CT>
template<typename TT>
inline void Problem_MaxwellEqnsFDTD<CT>::funcZero(TT & fx,
                                   ScaFES::Ntuple<TT,3> const& x, TT const& t)
{
    fx = 0.0;
}
/*----------------------------------------------------------------------------*/
template<class CT>
template<typename TT>
inline void Problem_MaxwellEqnsFDTD<CT>::funcSourceZ(TT& fx,
                                                 ScaFES::Ntuple<TT,3> const& x,
                                                 TT const& t)
{
    TT nDecay = 10.0 * this->tau();
    TT n0     = 3.0 * nDecay;
    if ( (fabs(x[0] - 1.0) < 2.2e-15)
      && (fabs(x[1] - 1.0) < 2.2e-15)
      && (fabs(x[2] - 1.0) < 2.2e-15) ) {
        fx = 20.0 * exp(-1.0 * ((this->tau() * this->tau()
                                 / this->nTimesteps() - n0)
                        * (this->tau() * this->tau() / this->nTimesteps() - n0)
                          / (nDecay * nDecay)));
    } else {
        fx = 0.0;
    }
}
/*----------------------------------------------------------------------------*/
template<class CT>
template<typename TT>
inline void Problem_MaxwellEqnsFDTD<CT>::updateInner(
                      std::vector< ScaFES::DataField<TT, 3> >& vectNew,
                      std::vector< ScaFES::DataField<TT, 3> > const& vectOld,
                      Index const& idxNode,
                      int const& timestep)
{
    int i = idxNode[0];
    int j = idxNode[1];
    int k = idxNode[2];
    TT tau_epsilon = this->tau() / vectOld[EPS_R](idxNode);
    TT tau_epsilon_hx = tau_epsilon / this->gridsize(0);
    TT tau_epsilon_hy = tau_epsilon / this->gridsize(1);
    TT tau_epsilon_hz = tau_epsilon / this->gridsize(2);

    vectNew[EX](idxNode)
            = vectOld[DZ](idxNode)
            + tau_epsilon_hy
            * (vectNew[HZ](idxNode) - vectNew[HZ](i,j+1,k))
            - tau_epsilon_hz
            * (vectNew[HY](idxNode) - vectNew[HY](i,j,k-1));
    vectNew[EY](idxNode)
            = vectOld[EY](idxNode)
            + tau_epsilon_hz
            * (vectNew[HX](idxNode) - vectNew[HX](i,j,k-1))
            - tau_epsilon_hx
            * (vectNew[HZ](idxNode) - vectNew[HZ](i-1,j,k));
    vectNew[EZ](idxNode)
            = vectOld[EY](idxNode)
            + tau_epsilon_hx
            * (vectNew[HY](idxNode) - vectNew[HY](i-1,j,k))
            - tau_epsilon_hy
            * (vectNew[HX](idxNode) - vectNew[HX](i,j+1,k));

    // Set source. //
    vectNew[EX](idxNode) += this->knownDf(0, idxNode);
    vectNew[EY](idxNode) += this->knownDf(1, idxNode);
    vectNew[EZ](idxNode) += this->knownDf(2, idxNode);

    // Update of H
    TT tau_mu = this->tau() / vectOld[MUE_R](idxNode);
    TT tau_mu_hx = tau_mu / this->gridsize(0);
    TT tau_mu_hy = tau_mu / this->gridsize(1);
    TT tau_mu_hz = tau_mu / this->gridsize(2);
    vectNew[HX](idxNode)
        = vectOld[HX](idxNode)
        + tau_mu_hz * (vectOld[EY](i,j,k+1) - vectOld[EY](idxNode))
        - tau_mu_hy * (vectOld[EZ](i,j+1,k) - vectOld[EY](idxNode));
    vectNew[HY](idxNode)
        = vectOld[HY](idxNode)
          + tau_mu_hx * (vectOld[EZ](i+1,j,k) - vectOld[EZ](idxNode))
          - tau_mu_hz * (vectOld[EX](i,j,k+1) - vectOld[EX](idxNode));
    vectNew[HZ](idxNode)
        = vectOld[HZ](idxNode)
          + tau_mu_hy * (vectOld[EX](i,j+1,k) - vectOld[EX](idxNode))
          - tau_mu_hx * (vectOld[EY](i+1,j,k) - vectOld[EY](idxNode));

    double c0 = 2.0 * EPSILON_0;
    double c1 = c0 * this->tau();
    TT koeffPx = c0 + vectOld[SIGMA](idxNode) * this->tau();
    TT koeffPy = c0 + vectOld[SIGMA](idxNode) * this->tau();
    TT koeffPz = c0 + vectOld[SIGMA](idxNode) * this->tau();
    TT koeffMx = c0 - vectOld[SIGMA](idxNode) * this->tau();
    TT koeffMy = c0 - vectOld[SIGMA](idxNode) * this->tau();
    TT koeffMz = c0 - vectOld[SIGMA](idxNode) * this->tau();
    TT divx = koeffMx / koeffPx;
    TT divy = koeffMy / koeffPy;
    TT divz = koeffMz / koeffPz;

    vectNew[BX](idxNode)
            = divy * vectOld[BX](idxNode)
            - c1 / koeffPy * (
              (vectOld[EZ](i,j-1,k) - vectOld[EY](idxNode)) / this->gridsize(1)
            - (vectOld[EY](i,j,k+1) - vectOld[EY](idxNode)) / this->gridsize(2)
            );
    vectNew[HX](idxNode)
            = divz * vectOld[HX](idxNode)
            + 1.0 / koeffPz
            * (koeffPx * vectNew[BX](idxNode) - koeffMx * vectOld[BX](idxNode))
            / vectOld[MUE_R](idxNode);

    vectNew[BY](idxNode)
            = divy * vectOld[BY](idxNode)
            - c1 / koeffPz * (
              (vectOld[EX](i,j-1,k) - vectOld[DZ](idxNode)) / this->gridsize(2)
            - (vectOld[EZ](i,j,k+1) - vectOld[EY](idxNode)) / this->gridsize(0)
            );
    vectNew[HY](idxNode)
            = divx * vectOld[HY](idxNode)
            + 1.0 / koeffPx
            * (koeffPy * vectNew[BY](idxNode) - koeffMy * vectOld[BY](idxNode))
            / vectOld[MUE_R](idxNode);

    vectNew[BZ](idxNode)
            = divx * vectOld[BZ](idxNode)
            - c1 / koeffPx * (
              (vectOld[EY](i+1,j,k) - vectOld[EY](idxNode)) / this->gridsize(0)
            - (vectOld[EX](i,j-1,k) - vectOld[DZ](idxNode)) / this->gridsize(1)
            );
    vectNew[HZ](idxNode)
            = divy * vectOld[HZ](idxNode)
            + 1.0 / koeffPy
            * (koeffPz * vectNew[BZ](idxNode) - koeffMz * vectOld[BZ](idxNode))
            / vectOld[MUE_R](idxNode);
}
/*----------------------------------------------------------------------------*/
template<class CT>
template<typename TT>
inline void Problem_MaxwellEqnsFDTD<CT>::updateBorder(
                      std::vector< ScaFES::DataField<TT, 3> >& vectNew,
                      std::vector< ScaFES::DataField<TT, 3> > const& vectOld,
                      Index const& idxNode,
                      int const& timestep)
{
    int i = idxNode[0];
    int j = idxNode[1];
    int k = idxNode[2];
    double c0 = 2.0 * EPSILON_0;
    double c1 = c0 * this->tau();
    TT koeffPx = c0 + vectOld[SIGMA](idxNode) * this->tau();
    TT koeffPy = c0 + vectOld[SIGMA](idxNode) * this->tau();
    TT koeffPz = c0 + vectOld[SIGMA](idxNode) * this->tau();
    TT koeffMx = c0 - vectOld[SIGMA](idxNode) * this->tau();
    TT koeffMy = c0 - vectOld[SIGMA](idxNode) * this->tau();
    TT koeffMz = c0 - vectOld[SIGMA](idxNode) * this->tau();
    TT divx = koeffMx / koeffPx;
    TT divy = koeffMy / koeffPy;
    TT divz = koeffMz / koeffPz;

    vectNew[DX](idxNode)
            = divy * vectOld[DX](idxNode)
            + c1 / koeffPy * (
            (vectOld[HZ](idxNode) - vectOld[HZ](i,j+1,k)) / this->gridsize(1)
          - (vectOld[HY](idxNode) - vectOld[HY](i,j,k-1)) / this->gridsize(2)
            );
    vectNew[EX](idxNode)
            = divz * vectOld[DZ](idxNode)
            + 1.0 / koeffPz
            * (koeffPx * vectNew[DX](idxNode) - koeffMx * vectOld[DX](idxNode))
            / vectOld[EPS_R](idxNode);

    vectNew[DY](idxNode)
            = divz * vectOld[DY](idxNode)
            + c1 / koeffPz * (
            (vectOld[HX](idxNode) - vectOld[HX](i,j,k-1)) / this->gridsize(2)
          - (vectOld[HZ](idxNode) - vectOld[HZ](i-1,j,k)) / this->gridsize(0)
            );
    vectNew[EY](idxNode)
            = divx * vectOld[EY](idxNode)
            + 1.0 / koeffPx
            * (koeffPy * vectNew[DY](idxNode) - koeffMy * vectOld[DY](idxNode))
            / vectOld[EPS_R](idxNode);

    vectNew[DZ](idxNode)
            = divx * vectOld[DZ](idxNode)
            + c1 / koeffPx * (
            (vectOld[HY](idxNode) - vectOld[HY](i-1,j,k)) / this->gridsize(0)
          - (vectOld[HX](idxNode) - vectOld[HX](i,j+1,k)) / this->gridsize(1)
            );
    vectNew[EZ](idxNode)
            = divy * vectOld[EY](idxNode)
            + 1.0 / koeffPy
            * (koeffPz * vectNew[DZ](idxNode) - koeffMz * vectOld[DZ](idxNode))
            / vectOld[EPS_R](idxNode);


    vectNew[BX](idxNode)
            = divy * vectOld[BX](idxNode)
            - c1 / koeffPy * (
              (vectOld[EZ](i,j-1,k) - vectOld[EY](idxNode)) / this->gridsize(1)
            - (vectOld[EY](i,j,k+1) - vectOld[EY](idxNode)) / this->gridsize(2)
            );
    vectNew[HX](idxNode)
            = divz * vectOld[HX](idxNode)
            + 1.0 / koeffPz
            * (koeffPx * vectNew[BX](idxNode) - koeffMx * vectOld[BX](idxNode))
            / vectOld[MUE_R](idxNode);

    vectNew[BY](idxNode)
            = divy * vectOld[BY](idxNode)
            - c1 / koeffPz * (
              (vectOld[EX](i,j-1,k) - vectOld[DZ](idxNode)) / this->gridsize(2)
            - (vectOld[EZ](i,j,k+1) - vectOld[EY](idxNode)) / this->gridsize(0)
            );
    vectNew[HY](idxNode)
            = divx * vectOld[HY](idxNode)
            + 1.0 / koeffPx
            * (koeffPy * vectNew[BY](idxNode) - koeffMy * vectOld[BY](idxNode))
            / vectOld[MUE_R](idxNode);

    vectNew[BZ](idxNode)
            = divx * vectOld[BZ](idxNode)
            - c1 / koeffPx * (
              (vectOld[EY](i+1,j,k) - vectOld[EY](idxNode)) / this->gridsize(0)
            - (vectOld[EX](i,j-1,k) - vectOld[DZ](idxNode)) / this->gridsize(1)
            );
    vectNew[HZ](idxNode)
            = divy * vectOld[HZ](idxNode)
            + 1.0 / koeffPy
            * (koeffPz * vectNew[BZ](idxNode) - koeffMz * vectOld[BZ](idxNode))
            / vectOld[MUE_R](idxNode);
}

#endif
