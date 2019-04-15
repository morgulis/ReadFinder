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
/*! \file libtools/exception.hpp
    \brief Tools to assist with throwing exceptions.
*/

#ifndef LIBTOOLS_EXCEPTION_HPP
#define LIBTOOLS_EXCEPTION_HPP

#include <cerrno>
#include <cstring>
#include <sstream>

#include <libtools/defs.hpp>

//------------------------------------------------------------------------------
/*! \brief Throw a runtime exception.

    Adds exception origin (source file and line number) to the error message.
    The message itself can be formatted in C++ stream like manner.

    \b Example:
    \code{.cpp}
        M_THROW( "failed to read from file " << file_name );
    \endcode
*/
#define M_THROW(_m) { \
        std::ostringstream _os; \
        _os << _m << " [" << __FILE__ << ':' << __LINE__ << ']'; \
        throw std::runtime_error( _os.str().c_str() ); \
    }

//------------------------------------------------------------------------------
/*! \brief Throw a system exception.
*/
#define M_THROW_SYSTEM(_m,_e) {\
        M_THROW( _m << '(' << _e << " : " << strerror( _e ) << ')' ) \
    }

//------------------------------------------------------------------------------
/*! \brief Throw a system exception based on errno.
*/
#define M_THROW_ERRNO(_m) M_THROW_SYSTEM( _m, errno )

//------------------------------------------------------------------------------
/** \brief Replacement for a body of unimplemented function/method.
*/
#define M_TODO M_THROW( "FUNCTION IS NOT IMPLEMENTED" )

#endif

