/* ScaFES
 * Copyright (c) 2011-2015, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file ScaFES_Buffer.hpp
 *  @brief Wrapper include file: Includes dummy or real version depending
 *         on availability of Boost MPI.
 */

#ifndef SCAFES_BUFFER_HPP_
#define SCAFES_BUFFER_HPP_

#include "ScaFES_Config.hpp"

#ifdef SCAFES_HAVE_BOOST_MPI
#include "ScaFES_BufferApi.hpp"
#else
#include "ScaFES_BufferDummy.hpp"
#endif

#endif

