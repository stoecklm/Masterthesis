/* ScaFES
 * Copyright (c) 2011-2015, 2018, ZIH, TU Dresden, Federal Republic of Germany.
 * For details, see the files COPYING and LICENSE in the base directory
 * of the package.
 */

/**
 *  @file ScaFES_DataFile.hpp
 *  @brief Wrapper include file: Includes dummy or real version depending
 *         on availability of NetCDF.
 */

#ifndef SCAFES_DATAFILE_HPP_
#define SCAFES_DATAFILE_HPP_

#include "ScaFES_Config.hpp"

namespace ScaFES
{
/*------------------------------------------------------------------------------
| ENUMERATIONS.
------------------------------------------------------------------------------*/
/** Possible types of writing data to file. */
enum class WriteHowOften
{
    NEVER = 1,
    AT_START = 2,
    AT_END = 3,
    AT_START_AND_END = 4,
    LIKE_GIVEN_AT_CL = 8,
    ALWAYS = 16,
};
}

#ifdef SCAFES_HAVE_NETCDF
#include "ScaFES_DataFileApi.hpp"
#else
#include "ScaFES_DataFileDummy.hpp"
#endif

#endif
