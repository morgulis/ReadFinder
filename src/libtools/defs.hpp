/*  $Id$
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*
* Authors:  Aleksandr Morgulis
*
*/
//------------------------------------------------------------------------------
/*! \file libtools/defs.hpp
    \brief High level definitions for tools library.
*/

#ifndef LIBTOOLS_DEFS_HPP
#define LIBTOOLS_DEFS_HPP

#include <cassert>
#include <iostream>
#include <mutex>
#include <type_traits>

#include <config.h>

//------------------------------------------------------------------------------
/*! \brief Internal namespace name for libtools entities.

    If macro \c OUTER_NS is defined, then its value becomes the enclosing
    namespace for \c ATB_NS::TOOLS_NS_LCL.
*/
#define TOOLS_NS_LCL tools

//------------------------------------------------------------------------------
/*! \def TOOLS_NS
    \brief Fully qualified namespace name for libtools entities.

    The value should be used when a qualifier for a libtools entity is needed.
*/
/*! \def TOOLS_NS_BEGIN
    \brief Open libtools namespace.
*/
/*! \def TOOLS_NS_END
    \brief Close libtools namespace.
*/
#ifdef OUTER_NS
#   define TOOLS_NS OUTER_NS::ATB_NS::TOOLS_NS_LCL
#   define TOOLS_NS_BEGIN namespace OUTER_NS { \
                          namespace ATB_NS { \
                          namespace TOOLS_NS_LCL {
#   define TOOLS_NS_END }}}
#else
#   define TOOLS_NS ATB_NS::TOOLS_NS_LCL
#   define TOOLS_NS_BEGIN namespace ATB_NS { namespace TOOLS_NS_LCL {
#   define TOOLS_NS_END }}
#endif

TOOLS_NS_BEGIN

/// Number of bits in a byte.
static size_t const BYTE_BITS = (size_t)8;

///\var FLIP_BYTES
///\brief Flag set to indicate that byte order needs to be flipped when
///       converting between words of different sizes.

#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
static bool const FLIP_BYTES = false;
#else
static bool const FLIP_BYTES = true;
#endif

//------------------------------------------------------------------------------
// Commonly used types.
//
typedef std::lock_guard< std::mutex > TLock;    ///< Type for RAII mutex locks.

//------------------------------------------------------------------------------
template< typename T_Int1, typename T_Int2 >
auto Max( T_Int1 x, T_Int2 y )
    -> typename std::common_type< T_Int1, T_Int2 >::type
{
    typedef typename std::common_type< T_Int1, T_Int2 >::type Result;
    return std::max( (Result)x, (Result)y );
}

template< typename T_Int1, typename T_Int2 >
auto Min( T_Int1 x, T_Int2 y )
    -> typename std::common_type< T_Int1, T_Int2 >::type
{
    typedef typename std::common_type< T_Int1, T_Int2 >::type Result;
    return std::min( (Result)x, (Result)y );
}

TOOLS_NS_END

#endif

