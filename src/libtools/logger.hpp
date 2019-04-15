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
/*! \file libtools/logger.hpp Logger class definition. */

#ifndef LIBTOOLS_LOGGER_HPP
#define LIBTOOLS_LOGGER_HPP

#include <memory>
#include <sstream>

#include <libtools/log_handler.hpp>
#include <libtools/defs.hpp>

/*! Log a message via a given logger object, with given severity level.

    The macro is forwards the message to the logger \c Log() method 
    providing source file name and line information from the point
    of macro invocation.

    The message can be formatted in C++ stream like manner.

    \b Example:
    \code{.cpp}
        M_LOG( logger_, CLogHandler::INFO, "number of sequences: " << n_seq );
    \endcode
*/
#define M_LOG(_o,_l,_m) { \
        std::ostringstream _os; _os << _m; \
        (_o).Log( (_l), _os.str(), __FILE__, __LINE__ ); \
    }

#define M_FORCE_LOG(_o,_m) { \
        std::ostringstream _os; _os << _m; \
        (_o).Log( TOOLS_NS::CLogHandler::INFO, _os.str(), \
                  __FILE__, __LINE__, true ); \
    }

/*!
    \def M_INFO(o,m)
    Logs an informational message.
*/

/*!
    \def M_WARN(o,m)
    Logs a warning message.
*/

/*!
    \def M_ERR(o,m)
    Logs an error message.
*/

#define M_INFO(_o,_m) M_LOG( _o, TOOLS_NS::CLogHandler::INFO, _m )
#define M_WARN(_o,_m) M_LOG( _o, TOOLS_NS::CLogHandler::WARNING, _m )
#define M_ERR(_o,_m) M_LOG( _o, TOOLS_NS::CLogHandler::ERROR, _m )

TOOLS_NS_BEGIN

//------------------------------------------------------------------------------
/** No-op logger class.

    This class just discards all log requests. Its intended use is as
    tracer, when the tracing is disabled.

    All methods of the class do nothing.
*/
class CNopLogger
{
public:

    /** Type used to hold a reference to the installed log handler.

        Loggers assume shared ownership of their log handlers.
    */
    typedef std::shared_ptr< CLogHandler > LogHandler;

    CNopLogger( std::string const & id = "", 
                LogHandler handler = LogHandler() );

    void Log( CLogHandler::Severity l, std::string const & msg,
              char const * file, int lineno, bool = false );

    void SetSeverity( CLogHandler::Severity l );
};

//==============================================================================
// IMPLEMENTATION
//==============================================================================
inline CNopLogger::CNopLogger( std::string const & id, LogHandler handler ) {}

//------------------------------------------------------------------------------
inline void CNopLogger::SetSeverity( CLogHandler::Severity l ) {}

//------------------------------------------------------------------------------
inline void CNopLogger::Log( 
        CLogHandler::Severity l, std::string const & msg, 
        char const * file, int lineno, bool )
{}
//==============================================================================

//------------------------------------------------------------------------------
/** Generic logger.

    Logger objects are responsible for filtering log requests based on set
    severity level and pre-formatting log messages with log origin in the
    source code.

    A logger object has a string identifier that is prepended to the log
    message. This can help to identify log origin when, e.g., multiple 
    loggers send messages to the same log file.

    Derives from CNopLogger only to inherit some common type definitions.
*/
class CLogger : public CNopLogger
{
public:

    /** Instance constructor.

        Default severity level threshold is set to CLogHandler::QUIET,
        so no logging will take place.

        \param [in] id Logger identification string.
        \param [in] handler Log handler to use for actual message output.
    */
    CLogger( std::string const & id, LogHandler handler = LogHandler() );

    /** Process a log request.

        If no log handler is installed, or if severity level threshold is
        set to CLogHandler::QUIET then nothing is logged.

        Otherwise, if \c \b b is at most the severity level threshold, then
        the msg is formatted and processed via the installed log handler.

        \param [in] l       Log message severity level.
        \param [in] msg     Log message text.
        \param [in] file    Name of the source file that is the message origin.
        \param [in] lineno  Line number in \c \b file of the message origin.
    */
    void Log( CLogHandler::Severity l, std::string const & msg,
              char const * file, int lineno, bool force = false );

    /** Set the severity level threshold.

        \param [in] l New severity level.
    */
    void SetSeverity( CLogHandler::Severity l );

private:

    CLogHandler::Severity l_;   ///< Severity level threshold.
    std::string const id_;      ///< Logger identification string.
    LogHandler handler_;        ///< Reference to the installed log handler.
};

//==============================================================================
// IMPLEMENTATION
//==============================================================================
inline CLogger::CLogger( std::string const & id, LogHandler handler )
    : l_( CLogHandler::QUIET ), id_( id ), handler_( handler )
{}

//------------------------------------------------------------------------------
inline void CLogger::SetSeverity( CLogHandler::Severity l )
{
    l_ = l;
}

TOOLS_NS_END

#endif

