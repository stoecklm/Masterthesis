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

template<typename TT>
TT boundaryCondition2ndTypeIndexPlusOne(TT valueIndexMinusOne,
                                        double const& gridsize,
                                        double const& lambda,
                                        double const& q_dot) {
    return (valueIndexMinusOne - ((2.0*gridsize)/lambda)*q_dot);
}

template<typename TT>
TT boundaryCondition2ndTypeIndexMinusOne(TT valueIndexplusOne,
                                         double const& gridsize,
                                         double const& lambda,
                                         double const& q_dot) {
    return (valueIndexplusOne + ((2.0*gridsize)/lambda)*q_dot);
}

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
