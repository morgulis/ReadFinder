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

/** \file libtools/progress.cpp
    \brief Simple progress bar.
*/

#include <chrono>

#include <boost/format.hpp>

#include <libtools/progress.hpp>

TOOLS_NS_BEGIN

//==============================================================================
//------------------------------------------------------------------------------
CProgressBase::CProgressBase(
        std::string const & title, std::string const & units, int flags )
    : title_( title ), units_( units ), top_( 1 ),
      quiet_( (flags&QUIET) != 0 ), done_( false )
{
    if( (flags&START) != 0 )
    {
        Start();
    }
}

//------------------------------------------------------------------------------
auto CProgressBase::GetTop() -> ProgressHandle
{
    return ProgressHandle( this, &this->top_ );
}

//------------------------------------------------------------------------------
void CProgressBase::MonitorProc( CProgressBase & self )
{
    self.Monitor();
}

//------------------------------------------------------------------------------
void CProgressBase::Start()
{
    if( quiet_ || done_.load() || monitor_.joinable() )
    {
        return;
    }

    timer_.Start();
    stopped_ = false;
    monitor_ = std::thread( MonitorProc, std::ref( *this ) );
}

//------------------------------------------------------------------------------
void CProgressBase::Stop()
{
    if( quiet_ || done_.load() || !monitor_.joinable() )
    {
        return;
    }

    done_.store( true );
    monitor_.join();
    timer_.Stop();
    stopped_ = true;
}

//==============================================================================
//------------------------------------------------------------------------------
void CProgress::Monitor()
{
    static std::chrono::milliseconds const POLL_INTERVAL( 500 );
    current_santipct_ = 0;
    bool start( true );
    std::string progress_bar;

    while( !done_.load() )
    {
        int64_t santipct( *top_[0] );

        if( start || santipct != current_santipct_ )
        {
            start = false;
            current_santipct_ = santipct;
            progress_bar = title_ + ": [";
            auto plen( std::min( (size_t)(santipct/200),
                                 (size_t)PROGRESS_BAR_LEN ) );
            progress_bar += std::string( plen, '=' );
            progress_bar += std::string( PROGRESS_BAR_LEN - plen, ' ' );
            progress_bar += "] ";
            progress_bar +=
                (boost::format( "%6.2f%%" )%(santipct/100.0)).str();
            std::cerr << progress_bar << "\r" << std::flush;
        }

        std::this_thread::sleep_for( POLL_INTERVAL );
    }
}

//------------------------------------------------------------------------------
void CProgress::Stop()
{
    if( !stopped_ )
    {
        CProgressBase::Stop();
        int64_t total( top_[0].GetTotal() );
        std::cerr << std::string( title_.size() + PROGRESS_BAR_LEN + 12, ' ' )
                  << '\r' << std::flush;
        auto microseconds( timer_.GetElapsed().count() );
        double rate( microseconds == 0 ? (double)total
                                       : ((double)total)/microseconds );
        std::cerr << title_ + ": done in " << microseconds/1000 << " ms; "
                  << boost::format( "%.2f " )%(rate*1000000.0)
                  << units_ << "/sec" << std::endl;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
void CCounterProgress::Monitor()
{
    static std::chrono::milliseconds const POLL_INTERVAL( 500 );

    while( !done_.load() )
    {
        std::cerr << title_ + ": " << top_[0].current << ' ' << units_
                  << "\r" << std::flush;
        std::this_thread::sleep_for( POLL_INTERVAL );
    }
}

//------------------------------------------------------------------------------
void CCounterProgress::Stop()
{
    if( !stopped_ )
    {
        CProgressBase::Stop();
        auto microseconds( timer_.GetElapsed().count() );
        auto current( GetTop().GetCurrent() );
        double rate( microseconds == 0 ? (double)current
                                       : (double)current/microseconds );
        std::cerr << title_ + ": done in " << microseconds/1000 << " ms; "
                  << boost::format( "%.2f " )%(rate*1000000.0)
                  << units_ << "/sec" << std::endl;
    }
}

TOOLS_NS_END

