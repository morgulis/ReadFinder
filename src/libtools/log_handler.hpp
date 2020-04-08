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
/*! \file libtools/log_handler.hpp 
    Classes for handling pre-formatted log messages. */

#ifndef LIBTOOLS_LOG_HANDLER_HPP
#define LIBTOOLS_LOG_HANDLER_HPP

#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include <libtools/exception.hpp>
#include <libtools/defs.hpp>

TOOLS_NS_BEGIN

//------------------------------------------------------------------------------
/*! Base class for log message handlers.
*/
class CLogHandler
{
private:

    /** Default text representations of different severity levels. */
    static char const * SeverityStr[];

public:

    /** Severity levels.  */
    enum Severity : int
    {
        START_LEVEL = 0,      ///< Initial level.
        QUIET = START_LEVEL,  /**< No logging.  
                                   This should not appear in actual logs. */
        ERROR,      ///< Error level message.
        WARNING,    ///< Warning level message.
        INFO,       ///< Informational message.
        N_LEVELS    ///< Limit on the number of severity levels.
    };

    /** Convert a severity level to text representation.

        \param [in] s Severity level.

        \return Default textual representation of s.
    */
    static char const * Severity2Str( Severity s );

    /** Virtual destructor. */
    virtual ~CLogHandler() = 0;

    /** Message handling operator.

        \param [in] l   Severity level.
        \param [in] msg Log message.
    */
    virtual void operator()( Severity l, std::string const & msg ) = 0;
};

//------------------------------------------------------------------------------
/** Parse a string representation of severity level into enumerated value. */
inline std::istream & operator>>( 
        std::istream & is, CLogHandler::Severity & v )
{
    std::string vstr;
    is >> vstr;

    if( vstr == "quiet" )
    {
        v = CLogHandler::QUIET;
    }
    else if( vstr == "error" )
    {
        v = CLogHandler::ERROR;
    }
    else if( vstr == "warning" )
    {
        v = CLogHandler::WARNING;
    }
    else if( vstr == "info" )
    {
        v = CLogHandler::INFO;
    }
    else
    {
        M_THROW( "wrong log severity name: " << vstr <<
                 "; allowed values: quiet, error, warning, info" );
    }

    return is;
}

//==============================================================================
// IMPLEMENTATION
//==============================================================================
inline char const * CLogHandler::Severity2Str( Severity s )
{
    return SeverityStr[s];
}

//------------------------------------------------------------------------------
inline CLogHandler::~CLogHandler() {}
//==============================================================================

//------------------------------------------------------------------------------
/** Simple log handler that redirects log messages to a given 
    C++ output stream.
*/
class CStreamLogHandler : public CLogHandler
{
public:

    /** Type for holding output stream pointer.
    */
    typedef std::shared_ptr< std::ostream > TStreamHandle;

    /** Instance constructor.

        If \c \b h is not initialized with a valid pointer to a C++ output
        stream, then no logging is performed.

        Assumes shared ownership of the output stream.

        \param [in] id Logger identifier (arbitrary text).
        \param [in] h  Points to the output stream to use for logging.

        \throws std::runtime_error if output stream is in failed state
                                   after construction.
    */
    CStreamLogHandler( std::string const & id, 
                       TStreamHandle h = TStreamHandle() );

    /** Instance constructor.

        If \c \b osp is a null pointer, then no logging is performed.

        \param [in] id  Logger identifier (arbitrary text).
        \param [in] osp Points to the output stream to use for logging.

        \throws std::runtime_error if output stream is in failed state
                                   after construction.
    */
    CStreamLogHandler( std::string const & id, std::ostream * osp );

    /** Log the message with at given severity level.

        Prepends the textual representation of the severity level to
        the message text.

        \param [in] l   Severity level.
        \param [in] msg Log message text.

        \throws std::runtime_error if output stream is in failed state
                                   before or after logging.
    */
    void operator()( Severity l, std::string const & msg ) override;

private:

    /** Check the output stream for failure.

        \throws std::runtime_error if the output stream is in the failed 
                                   state.
    */
    void CheckFail();

    std::string id_;        ///< Log handler id.
    TStreamHandle os_;      ///< Holds a reference to the output stream.

    std::ostream * osp_ = nullptr; ///< Actual pointer to the output stream.
};

//==============================================================================
// IMPLEMENTATION
//==============================================================================
inline void CStreamLogHandler::CheckFail()
{
    if( osp_ && osp_->fail() ) 
    {
        M_THROW( "stream log handler " << id_ << " failed" );
    }
}

//------------------------------------------------------------------------------
inline CStreamLogHandler::CStreamLogHandler( 
        std::string const & id, TStreamHandle h )
    : id_( id ), os_( h ), osp_( os_.get() )
{
    CheckFail();
}

//------------------------------------------------------------------------------
inline CStreamLogHandler::CStreamLogHandler( 
        std::string const & id, std::ostream * osp )
    : id_( id ), osp_( osp )
{
    CheckFail();
}

//------------------------------------------------------------------------------
inline void CStreamLogHandler::operator()( Severity l, std::string const & msg )
{
    if( osp_ )
    {
        CheckFail();
        *osp_ << Severity2Str( l ) << ": " << msg << std::endl;
        CheckFail();
    }
}
//==============================================================================

//------------------------------------------------------------------------------
/** Log handler for logging messages to a text file. 

    Derives from a stream log handler. Creates an output stream corresponding
    to a given file, and delegates logging to CStreamLogHandler.

    The output file is truncated on open.
*/
class CFileLogHandler : public CStreamLogHandler
{
public:

    /** Instance constructor.

        \param [in] fname Name of the log file. This is used as an identifier
                          for the stream log handler.
    */
    CFileLogHandler( std::string const & fname );
};

//==============================================================================
// IMPLEMENTATION
//==============================================================================
inline CFileLogHandler::CFileLogHandler( std::string const & fname )
    : CStreamLogHandler( fname, 
                         TStreamHandle( new std::ofstream( 
                                 fname.c_str(), std::ios_base::app ) ) )
{}

TOOLS_NS_END

#endif

