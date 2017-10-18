/* ScaFES
 * Copyright (c) 2017, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file boundaryConditions.hpp
 *  @brief Contains functions which can be used as boundary conditions.
 */

#ifndef BOUNDARY_CONDITIONS_HPP_
#define BOUNDARY_CONDITIONS_HPP_

/** Calculates value for node which is outside of grid
 *  by using Neumann boundary condition.
 *  To be used if rhs neighbour (index+1) is outside of grid.
 *  @param valueIndexMinusOne Value of lhs neighbour (index-1).
 *  @param gridsize Distance between two grid points.
 *  @param lambda Material parameter (thermal conductivity).
 *  @param q_dot Boundary condition (heat flux).
 */
template<typename TT>
TT boundaryCondition2ndTypeIndexPlusOne(TT valueIndexMinusOne,
                                        double const& gridsize,
                                        double const& lambda,
                                        double const& q_dot) {
    return (valueIndexMinusOne - ((2.0*gridsize)/lambda)*q_dot);
}

/** Calculates value for node which is outside of grid
 *  by using Neumann boundary condition.
 *  To be used if lhs neighbour (index-1) is outside of grid.
 *  @param valueIndexMinusOne Value of rhs neighbour (index+1).
 *  @param gridsize Distance between two grid points.
 *  @param lambda Material parameter (thermal conductivity).
 *  @param q_dot Boundary condition (heat flux).
 */
template<typename TT>
TT boundaryCondition2ndTypeIndexMinusOne(TT valueIndexplusOne,
                                         double const& gridsize,
                                         double const& lambda,
                                         double const& q_dot) {
    return (valueIndexplusOne + ((2.0*gridsize)/lambda)*q_dot);
}

/** Calculates value for node which is outside of grid
 *  by using Cauchy boundary condition.
 *  To be used if rhs neighbour (index+1) is outside of grid.
 *  @param valueIndexMinusOne Value of lhs neighbour (index-1).
 *  @param valueIndex Value of node on boundary (index).
 *  @param T_inf Ambient temperature.
 *  @param gridsize Distance between two grid points.
 *  @param alpha Material parameter (heat transer coefficient).
 *  @param lambda Material parameter (thermal conductivity).
 */
template<typename TT>
TT boundaryCondition3rdTypeIndexPlusOne(TT valueIndexMinusOne,
                                        TT valueIndex,
                                        double const& T_inf,
                                        double const& gridsize,
                                        double const& alpha,
                                        double const& lambda) {
    return (valueIndexMinusOne
            - (2.0*gridsize*(alpha/lambda)*(valueIndex - T_inf)));
}

/** Calculates value for node which is outside of grid
 *  by using Cauchy boundary condition.
 *  To be used if lhs neighbour (index-1) is outside of grid.
 *  @param valueIndexMinusOne Value of rhs neighbour (index+1).
 *  @param valueIndex Value of node on boundary (index).
 *  @param T_inf Ambient temperature.
 *  @param gridsize Distance between two grid points.
 *  @param alpha Material parameter (heat transer coefficient).
 *  @param lambda Material parameter (thermal conductivity).
 */
template<typename TT>
TT boundaryCondition3rdTypeIndexMinusOne(TT valueIndexPlusOne,
                                         TT valueIndex,
                                         double const& T_inf,
                                         double const& gridsize,
                                         double const& alpha,
                                         double const& lambda) {
    return (valueIndexPlusOne
            - (2.0*gridsize*(alpha/lambda)*(valueIndex - T_inf)));
}

#endif
