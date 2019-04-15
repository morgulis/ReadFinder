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

/** \file libtools/stopwatch.hpp
    \brief Simple stop watch class used to measure tasks.
*/

#ifndef LIBTOOLS_STOPWATCH_HPP
#define LIBTOOLS_STOPWATCH_HPP

#include <chrono>

#include <libtools/logger.hpp>
#include <libtools/defs.hpp>

TOOLS_NS_BEGIN

//==============================================================================
class Timer
{
public:

    void Start()
    {
        start_time = std::chrono::system_clock::now();
    }

    void Stop()
    {
        end_time = std::chrono::system_clock::now();
    }

    std::chrono::microseconds GetElapsed() const
    {
        using namespace std::chrono;
        return duration_cast< microseconds >( end_time - start_time );
    }

private:

    std::chrono::system_clock::time_point start_time,
                                          end_time;
};

//==============================================================================
class StopWatch : public Timer
{
public:

    StopWatch( CLogger & logger, std::string const & hdr, bool start = true )
        : logger_( &logger ), hdr_( hdr )
    {
        if( start )
        {
            Start();
        }
    }

    StopWatch( std::ostream & os, std::string const & hdr, std::mutex & mtx,
               bool start = true )
        : os_( &os ), hdr_( hdr ), mtx_( &mtx )
    {
        if( start )
        {
            Start();
        }
    }

    StopWatch( CLogger & logger, bool start = true )
        : StopWatch( logger, "", start )
    {}

    ~StopWatch()
    {
        Stop();
        std::ostringstream oss;
        oss << hdr_ << (hdr_.empty() ? "" : ": ") << "complete in "
            << GetElapsed().count() << " microseconds";

        if( logger_ != nullptr )
        {
            M_INFO( *logger_, oss.str() );
        }
        else if( os_ != nullptr )
        {
            TLock lock( *mtx_ );
            *os_ << oss.str() << std::endl;
        }
    }

private:

    CLogger * logger_ = nullptr;
    std::ostream * os_ = nullptr;
    std::string hdr_;
    std::mutex * mtx_ = nullptr;
};

TOOLS_NS_END

#endif

