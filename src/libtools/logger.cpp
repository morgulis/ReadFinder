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
/*! \file libtools/logger.cpp Implementation of CLogger class. */

#include <mutex>
#include <string>

#include <libtools/logger.hpp>

TOOLS_NS_BEGIN

//------------------------------------------------------------------------------
static std::mutex Log_Mtx;  ///< Global mutex to sync all logging.

//------------------------------------------------------------------------------
void CLogger::Log( CLogHandler::Severity l, std::string const & msg, 
                   char const * file, int lineno, bool force )
{
    if( handler_ ) 
    {
        if( l != CLogHandler::QUIET && (l <= l_ || force) )
        {
            TLock lock( Log_Mtx );
            (*handler_)( l, id_ + ": " + msg + "[" + file + ":" + 
                            std::to_string( lineno ) + "]" );
        }
    }
}

TOOLS_NS_END

