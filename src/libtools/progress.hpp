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

/** \file libtools/progress.hpp
    \brief Simple progress bar.
*/

#ifndef LIBTOOLS_PROGRESS_HPP
#define LIBTOOLS_PROGRESS_HPP

#include <atomic>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

#include <libtools/stopwatch.hpp>
#include <libtools/defs.hpp>

TOOLS_NS_BEGIN

//==============================================================================
class CProgressBase
{
private:

    static int64_t const FACTOR = 10000;

    struct SimpleProgress;

    typedef std::vector< SimpleProgress > Parts;
    typedef Parts::iterator PartsIter;

    struct SimpleProgress
    {
        SimpleProgress( int64_t t = 1, int64_t c = 0 )
        {
            if( t <= 0 ) t = 1;
            current.store( c );
            total.store( t );
        }

        SimpleProgress( SimpleProgress const & x )
        {
            current.store( x.current.load() );
            total.store( x.total.load() );
        }

        int64_t operator*() const
        {
            int64_t res( 0 );
            int64_t t( total.load() );

            if( t == 0 )
            {
                if( parts && !parts->empty() )
                {
                    auto denom( parts->size() );

                    for( auto const & part : *parts )
                    {
                        res += (*part)/denom;
                    }
                }
                else
                {
                    return 0;
                }
            }
            else
            {
                int64_t c( current.load() );
                res = c*FACTOR/t;
            }

            return res;
        }

        int64_t GetTotal() const
        {
            int64_t t( total.load() );

            if( t == 0 && parts != nullptr )
            {
                int64_t res( 0 );

                for( auto const & part : *parts )
                {
                    res += part.GetTotal();
                }

                return res;
            }

            return t;
        }

        std::atomic< int64_t > current,
                               total;
        std::shared_ptr< Parts > parts;
    };

public:

    class ProgressHandle;

    enum : int {
        QUIET = 1,
        START = 2,
    };

    CProgressBase( std::string const & title,
               std::string const & units,
               int flags = START );
    virtual ~CProgressBase() {}
    bool IsQuiet() const { return quiet_; }

    void Start();
    virtual void Stop();

    ProgressHandle GetTop();

    std::chrono::microseconds GetElapsed() const
    {
        return timer_.GetElapsed();
    }

protected:

    static size_t const PROGRESS_BAR_LEN = 50;

    static void MonitorProc( CProgressBase & self );

    virtual void Monitor() = 0;

    std::string const title_,
                      units_;
    Parts top_;
    int64_t current_santipct_ = 0;
    std::thread monitor_;
    Timer timer_;
    bool quiet_;
    bool stopped_ = true;
    std::atomic< bool > done_;
    std::mutex mtx_;
};

//==============================================================================
class CProgressBase::ProgressHandle
{
public:

    ProgressHandle() : o_( nullptr ), p_( nullptr ) {}

    ProgressHandle( CProgressBase * o, Parts * p, size_t idx = 0 )
        : o_( o ), p_( p ), c_( p_->begin() + idx )
    {}

    void SetTotal( int64_t new_total )
    {
        c_->total.store( new_total );
    }

    int64_t GetCurrent() const { return c_->current.load(); }

    void SetCurrent( int64_t new_current )
    {
        c_->current.store( new_current );
    }

    void Increment( int64_t inc = 1 )
    {
        c_->current.fetch_add( inc );
    }

    operator bool() const { return c_ != p_->end(); }

    void operator++()
    {
        assert( c_ != p_->end() );
        ++c_;
    }

    ProgressHandle Split( int64_t n_parts ) const
    {
        TLock lock( o_->mtx_ );
        assert( n_parts > 0 );
        assert( c_ != p_->end() );
        c_->parts.reset( new Parts( n_parts ) );
        c_->total.store( 0 );
        return ProgressHandle( o_, c_->parts.get() );
    }

    ProgressHandle operator[]( size_t idx ) const
    {
        return ProgressHandle( o_, p_, idx );
    }

private:

    CProgressBase * o_;
    Parts * p_;
    PartsIter c_;
};

//==============================================================================
class CProgress : public CProgressBase
{
public:

    CProgress( std::string const & title,
               std::string const & units,
               int flags = START )
        : CProgressBase( title, units, flags )
    {}

    virtual ~CProgress() { Stop(); }
    virtual void Stop() override;

protected:

    virtual void Monitor() override;
};

//==============================================================================
class CCounterProgress : public CProgressBase
{
public:

    CCounterProgress( std::string const & title,
                      std::string const & units,
                      int flags = START )
        : CProgressBase( title, units, flags )
    {}

    virtual ~CCounterProgress() { Stop(); }
    virtual void Stop() override;

protected:

    virtual void Monitor() override;
};

TOOLS_NS_END

#endif

