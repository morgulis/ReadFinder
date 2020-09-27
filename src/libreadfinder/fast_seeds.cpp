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

#include <algorithm>
#include <atomic>
#include <chrono>
#include <list>
#include <tuple>

#include <boost/format.hpp>

#include <libreadfinder/fast_seeds.hpp>

#include <libtools/taskarray.hpp>

READFINDER_NS_BEGIN

//==============================================================================
//------------------------------------------------------------------------------
CFastSeeds::CFastSeeds( CBatch & bctx )
    : bctx_( bctx ),
      fsidx_( *bctx.GetSearchCtx().fsidx ),
      anchor_use_map_( ANCHOR_TBL_SIZE, false )
{
    bctx_.GetSearchCtx().refs->LoadAll();
    prescreen_ = bctx_.GetSearchCtx().pre_screen;
}

//==============================================================================
class IndexScanner
{
public:

    IndexScanner( CSearchContext & ctx, CFastSeedsIndex & fsidx );
    void AcquireAnchor( uint32_t anchor, size_t thread_idx );
    void Finalize( size_t thread_idx );

private:

    struct ChunkInfo
    {
        uint32_t start_anchor,
                 end_anchor;
        std::atomic_flag requested_;
        std::atomic< bool > loaded_;
    };

    typedef std::vector< ChunkInfo > ChunkData;

    static void LoadProc( IndexScanner & self, size_t chunk )
    {
        self.Load( chunk );
    }

    void LoadAsync( size_t chunk, size_t thread_idx )
    {
        if( loading_threads_[thread_idx] != nullptr )
        {
            loading_threads_[thread_idx]->join();
            delete loading_threads_[thread_idx];
            loading_threads_[thread_idx] = nullptr;
        }

        loading_threads_[thread_idx] = new std::thread(
                IndexScanner::LoadProc, std::ref( *this ), chunk );
    }

    void Load( size_t chunk );

    CSearchContext & ctx_;
    CFastSeedsIndex & fsidx_;
    ChunkData chunks_;
    std::vector< size_t > used_chunks_;
    std::vector< std::thread * > loading_threads_;
    std::mutex load_mtx_;
};

//------------------------------------------------------------------------------
IndexScanner::IndexScanner( CSearchContext & ctx, CFastSeedsIndex & fsidx )
    : ctx_( ctx ), fsidx_( fsidx ), chunks_( fsidx_.GetNumChunks() + 1 ),
      used_chunks_( ctx_.n_threads, fsidx_.GetNumChunks() ),
      loading_threads_( ctx_.n_threads, nullptr )
{
    for( size_t i( 0 ); i < chunks_.size() - 1; ++i )
    {
        auto & ci( chunks_[i] );
        std::tie( ci.start_anchor, ci.end_anchor ) =
            fsidx_.GetChunkAnchorRange( i );
        ci.requested_.clear();
        ci.loaded_.store( false );
    }

    auto & ci( chunks_.back() );
    ci.start_anchor = ci.end_anchor = 0;
    ci.requested_.clear();
    ci.loaded_.store( false );
}

//------------------------------------------------------------------------------
void IndexScanner::Load( size_t chunk )
{
    TLock lock( load_mtx_ );
    auto & c( chunks_[chunk] );

    if( !c.loaded_.load() )
    {
        try
        {
            fsidx_.Load( chunk );
        }
        catch( std::exception const & e )
        {
            M_ERR( ctx_.logger_,
                   "error loading chunk " << chunk <<
                   ": " << e.what() );
        }
        catch( ... )
        {
            M_ERR( ctx_.logger_, "unknown error loading chunk " << chunk );
        }

        c.loaded_.store( true );
    }
}

//------------------------------------------------------------------------------
inline void IndexScanner::AcquireAnchor( uint32_t anchor, size_t thread_idx )
{
    auto * ci( &chunks_[used_chunks_[thread_idx]] );
    assert( anchor >= ci->start_anchor );

    if( anchor >= ci->end_anchor ) // entered new chunk
    {
        // update new current chunk
        //
        size_t n_chunks( fsidx_.GetNumChunks() );
        auto ccidx( used_chunks_[thread_idx] == n_chunks ?
                        0 : used_chunks_[thread_idx] );
        for( ; ccidx < chunks_.size() - 1
               && anchor >= chunks_[ccidx].end_anchor;
               ++ccidx );
        assert( ccidx < chunks_.size() - 1 );
        ci = &chunks_[ccidx];

        {
            TLock lock( load_mtx_ );
            used_chunks_[thread_idx] = ccidx;
            size_t min_used_chunk( n_chunks );

            for( size_t uc : used_chunks_ )
            {
                if( uc == n_chunks )
                {
                    min_used_chunk = 0;
                    break;
                }
                else if( uc < min_used_chunk )
                {
                    min_used_chunk = uc;
                }
            }

            for( size_t i( 0 ); i < min_used_chunk; ++i )
            {
                fsidx_.Unload( i );
            }
        }

        // request current chunk load if needed
        //
        {
            bool requested( ci->requested_.test_and_set() );

            if( !requested )
            {
                LoadAsync( ccidx, thread_idx );
            }
        }

        // wait until current chunk is loaded
        //
        while( !(ci->loaded_.load()) ); // spin

        if( !fsidx_.ChunkIsLoaded( ccidx ) )
        {
            M_THROW( "failed to load index chunk " << ccidx );
        }

        // request next chunk load if needed
        //
        if( ++ccidx < fsidx_.GetNumChunks() )
        {
            auto * next_chunk( &chunks_[ccidx] );
            bool requested( next_chunk->requested_.test_and_set() );

            if( !requested )
            {
                LoadAsync( ccidx, thread_idx );
            }
        }
    }
}

//------------------------------------------------------------------------------
void IndexScanner::Finalize( size_t thread_idx )
{
    assert( thread_idx < loading_threads_.size() );
    auto lt( loading_threads_[thread_idx] );

    if( lt != nullptr )
    {
        lt->join();
        delete lt;
        loading_threads_[thread_idx] = nullptr;
    }
}

//==============================================================================
using TOOLS_NS::StopWatch;

typedef CFastSeedsIndex::IndexEntry IndexEntry;

//==============================================================================
struct CFastSeeds::WordCountingJobData
{
    WordMap wmap;
    std::vector< uint32_t > tasks;
};

struct CFastSeeds::WordCountingJob
{
    typedef std::atomic< uint32_t > TaskIdx;
    typedef std::vector< WordCountingJobData > WCJData;

    WordCountingJob( CFastSeeds & o, WCJData & wcj_data,
                     TaskIdx & task_idx, size_t n_reads_per_job,
                     size_t & job_idx, CProgress::ProgressHandle & ph );
    void operator()();

    CFastSeeds & o_;
    WordCountingJobData & wcj_data_;
    TaskIdx & task_idx_;
    size_t n_reads_per_job_;
    CProgress::ProgressHandle & ph_;
};

//------------------------------------------------------------------------------
inline CFastSeeds::WordCountingJob::WordCountingJob(
        CFastSeeds & o, WCJData & wcj_data, TaskIdx & task_idx,
        size_t n_reads_per_job, size_t & job_idx,
        CProgress::ProgressHandle & ph )
    : o_( o ), wcj_data_( wcj_data[job_idx++] ), task_idx_( task_idx ),
      n_reads_per_job_( n_reads_per_job ), ph_( ph )
{
    wcj_data_.wmap.resize( WMAP_SIZE, 0 );
}

//------------------------------------------------------------------------------
inline void CFastSeeds::WordCountingJob::operator()()
{
    auto const & reads( o_.bctx_.GetReads() );
    auto n_reads( reads.GetNReads() );

    while( true )
    {
        uint32_t task_idx( task_idx_.fetch_add( 1 ) );
        OrdId start_read( task_idx*n_reads_per_job_ ),
              end_read( start_read + n_reads_per_job_ );

        if( start_read >= n_reads )
        {
            return;
        }

        wcj_data_.tasks.push_back( task_idx );

        if( end_read > n_reads )
        {
            end_read = n_reads;
        }

        for( ; start_read < end_read; ++start_read )
        {
            auto const & read( reads[start_read] );

            for( EStrand s : { eFWD, eREV } )
            {
                auto const & sd( reads.GetSeqData( start_read, eFIRST, s ) ),
                           & md( reads.GetMaskData( start_read, eFIRST, s ) );
                HashWordSource hws( sd.GetBuf(), sd.size(), sd.GetStart() );
                HashMaskSource hms( md.GetBuf(), md.size(), md.GetStart() );

                while( hws )
                {
                    assert( hms );

                    if( GetNSet( hms.GetNMer() ) == 0 )
                    {
                        ++wcj_data_.wmap[hws.GetAnchor()];
                    }

                    ++hws;
                    ++hms;
                }
            }

            if( read.IsPaired() )
            {
                for( EStrand s : { eFWD, eREV } )
                {
                    auto const & sd( reads.GetSeqData( 
                                        start_read, eSECOND, s ) ),
                               & md( reads.GetMaskData( 
                                        start_read, eSECOND, s ) );
                    HashWordSource hws( 
                            sd.GetBuf(), sd.size(), sd.GetStart() );
                    HashMaskSource hms( 
                            md.GetBuf(), md.size(), md.GetStart() );

                    while( hws )
                    {
                        assert( hms );

                        if( GetNSet( hms.GetNMer() ) == 0 )
                        {
                            ++wcj_data_.wmap[hws.GetAnchor()];
                        }

                        ++hws;
                        ++hms;
                    }
                }
            }

            ph_.Increment();
        }
    }
}

//==============================================================================
struct CFastSeeds::WordTableGeneratorJob
{
    typedef std::vector< WordCountingJobData > WCJData;

    WordTableGeneratorJob( 
            CFastSeeds & o, WCJData & wcj_data, size_t n_reads_per_job,
            size_t & job_idx, CProgress::ProgressHandle const & ph );
    void operator()();

    CFastSeeds & o_;
    WordCountingJobData & wcj_data_;
    size_t n_reads_per_job_;
    CProgress::ProgressHandle ph_;
};

//------------------------------------------------------------------------------
inline CFastSeeds::WordTableGeneratorJob::WordTableGeneratorJob(
        CFastSeeds & o, WCJData & wcj_data, size_t n_reads_per_job,
        size_t & job_idx, CProgress::ProgressHandle const & ph )
    : o_( o ), wcj_data_( wcj_data[job_idx] ),
      n_reads_per_job_( n_reads_per_job ), ph_( ph[job_idx++] )
{
}

//------------------------------------------------------------------------------
inline void CFastSeeds::WordTableGeneratorJob::operator()()
{
    auto const & reads( o_.bctx_.GetReads() );
    auto n_reads( reads.GetNReads() );
    auto & wmap( wcj_data_.wmap );
    auto & wt( o_.wt_ );
    ph_.SetTotal( wcj_data_.tasks.size() );

    for( auto task_idx : wcj_data_.tasks )
    {
        OrdId start_read( task_idx*n_reads_per_job_ ),
              end_read( start_read + n_reads_per_job_ );
        assert( start_read < n_reads );

        if( end_read > n_reads )
        {
            end_read = n_reads;
        }

        for( ; start_read < end_read; ++start_read )
        {
            auto const & read( reads[start_read] );

            for( EStrand s : { eFWD, eREV } )
            {
                auto const & sd( reads.GetSeqData( start_read, eFIRST, s ) ),
                           & md( reads.GetMaskData( start_read, eFIRST, s ) );
                HashWordSource hws( sd.GetBuf(), sd.size(), sd.GetStart() );
                HashMaskSource hms( md.GetBuf(), md.size(), md.GetStart() );

                while( hws )
                {
                    assert( hms );

                    if( GetNSet( hms.GetNMer() ) == 0 )
                    {
                        auto & w( wt[wmap[hws.GetAnchor()]++] );
                        w.wd.w.word = hws.GetWord();
                        w.readid = start_read;
                        w.hashoff = hws.GetHashOff();
                        w.strand = s;
                        w.mate = eFIRST;
                    }

                    ++hws;
                    ++hms;
                }
            }

            if( read.IsPaired() )
            {
                for( EStrand s : { eFWD, eREV } )
                {
                    auto const & sd( reads.GetSeqData( 
                                        start_read, eSECOND, s ) ),
                               & md( reads.GetMaskData( 
                                        start_read, eSECOND, s ) );
                    HashWordSource hws(
                            sd.GetBuf(), sd.size(), sd.GetStart() );
                    HashMaskSource hms(
                            md.GetBuf(), md.size(), md.GetStart() );

                    while( hws )
                    {
                        assert( hms );

                        if( GetNSet( hms.GetNMer() ) == 0 )
                        {
                            auto & w( wt[wmap[hws.GetAnchor()]++] );
                            w.wd.w.word = hws.GetWord();
                            w.readid = start_read;
                            w.hashoff = hws.GetHashOff();
                            w.strand = s;
                            w.mate = eSECOND;
                        }

                        ++hws;
                        ++hms;
                    }
                }
            }
        }

        ph_.Increment();
    }
}

//==============================================================================
struct CFastSeeds::WordTableSortJob
{
    typedef std::atomic< uint32_t > TaskIdx;

    WordTableSortJob( CFastSeeds & o, TaskIdx & task_idx,
                      CProgress::ProgressHandle & ph );
    void operator()();

    CFastSeeds & o_;
    TaskIdx & task_idx_;
    CProgress::ProgressHandle & ph_;
};

//------------------------------------------------------------------------------
inline CFastSeeds::WordTableSortJob::WordTableSortJob( 
        CFastSeeds & o, TaskIdx & task_idx, CProgress::ProgressHandle & ph )
    : o_( o ), task_idx_( task_idx ), ph_( ph )
{}

//------------------------------------------------------------------------------
inline void CFastSeeds::WordTableSortJob::operator()()
{
    auto const & wmap( o_.wmap_ );
    auto & wt( o_.wt_ );

    while( true )
    {
        uint32_t wa( task_idx_.fetch_add( 1 ) );

        if( wa >= WMAP_SIZE - 1 )
        {
            break;
        }

        auto b( wt.begin() + wmap[wa] ),
             e( wt.begin() + wmap[wa + 1] );
        std::sort( b, e,
                   []( HashWord const & x, HashWord const & y )
                   { return x.wd.w.word < y.wd.w.word; } );
        ph_.Increment();
    }
}

//==============================================================================
struct CFastSeeds::SeedSearchJobData
{
    size_t hits = 0,
           phits = 0,
           pmisses = 0,
           good_anchors = 0;
    std::vector< HashWord > words;
    std::vector< ExtHashWord > ewords;
    std::vector< IndexEntry > idxwords;
    Hits results;

    SeedSearchJobData & operator+=( SeedSearchJobData const & jd )
    {
        hits += jd.hits;
        phits += jd.phits;
        pmisses += jd.pmisses;
        good_anchors += jd.good_anchors;
        return *this;
    }
};

struct CFastSeeds::SeedSearchJob
{
    static ssize_t const GOOD_ANCHOR_SIZE = 32*1024ULL;

    static size_t const PARRAY_SIZE = 2*1024ULL;
    static size_t const SUBINDEX_SIZE = 2*1024ULL;

    static size_t const PWBITS = 64;
    static uint32_t const PWMASK = 0x1FFFFUL;
    static uint32_t const SIMASK = 0xFFE000UL;
    static size_t const SISHIFT = 13;

    typedef std::atomic< uint32_t > TaskIdx;
    typedef std::vector< TaskEntry > WordList;
    typedef std::vector< uint64_t > PArray;
    typedef std::vector< uint32_t > SubIndex;

    SeedSearchJob( CFastSeeds & o, IndexScanner & scanner,
                   TaskIdx & task_idx,
                   SeedSearchJobData * jd, size_t & job_idx,
                   CProgress::ProgressHandle & ph );
    void operator()();
    void PrepareIndex( uint32_t anchor, bool good = true );
    void MatchSubst( uint32_t anchor );
    void MatchScan( uint32_t anchor );

    void SaveHit( IndexEntry const & iw,
                  HashWord const & hw );

    static uint32_t MkSI( uint32_t w )
    {
        return (w&SIMASK)>>SISHIFT;
    }

    CFastSeeds & o_;
    IndexScanner & scanner_;
    TaskIdx & task_idx_;
    size_t job_idx_;
    SeedSearchJobData & jd_;
    PArray pa_ = PArray( PARRAY_SIZE );
    SubIndex sidx_ = SubIndex( SUBINDEX_SIZE + 1 );
    CProgress::ProgressHandle & ph_;
};

//------------------------------------------------------------------------------
inline CFastSeeds::SeedSearchJob::SeedSearchJob( 
        CFastSeeds & o, IndexScanner & scanner,
        TaskIdx & task_idx, SeedSearchJobData * jd,
        size_t & job_idx, CProgress::ProgressHandle & ph )
    : o_( o ), scanner_( scanner ), task_idx_( task_idx ), job_idx_( job_idx ),
      jd_( jd[job_idx++] ), ph_( ph )
{
}

//------------------------------------------------------------------------------
inline void CFastSeeds::SeedSearchJob::PrepareIndex( 
        uint32_t anchor, bool good )
{
    std::fill( pa_.begin(), pa_.end(), 0 );

    if( good )
    {
        std::fill( sidx_.begin(), sidx_.end(), 0 );
    }

    for( auto ib( o_.fsidx_.cbegin( anchor ) ),
              ie( o_.fsidx_.cend( anchor ) );
         ib != ie; ++ib )
    {
        auto w( (ib->wd.w.word)&PWMASK );
        SetBit( pa_[w/PWBITS], w%PWBITS );
    }

    if( good )
    {
        uint32_t si_off( 0 );
        auto ib( o_.fsidx_.cbegin( anchor ) ),
             ie( o_.fsidx_.cend( anchor ) );

        for( auto ic( ib ); ic != ie; )
        {
            auto w( MkSI( ic->wd.w.word ) );

            while( si_off < w )
            {
                sidx_[si_off++] = ic - ib;
            }

            assert( si_off == w );
            sidx_[si_off++] = ic++ - ib;
            for( ; ic != ie && MkSI( ic->wd.w.word ) == w; ++ic );
        }

        for( ; si_off <= SUBINDEX_SIZE; ++si_off )
        {
            sidx_[si_off] = ie - ib;
        }
    }
}

//------------------------------------------------------------------------------
inline void CFastSeeds::SeedSearchJob::SaveHit( 
        IndexEntry const & iw, HashWord const & hw )
{
    ++jd_.hits;
    TReadOff readoff( hw.hashoff );
    Hit h { iw.pos, hw.readid, readoff,
            (uint8_t)(hw.strand - 1), (uint8_t)(hw.mate - 1) };
    jd_.results.push_back( h );
}

//------------------------------------------------------------------------------
inline void CFastSeeds::SeedSearchJob::MatchSubst( uint32_t anchor )
{
    auto const & wt( o_.wt_ );
    auto const & wmap( o_.wmap_ );
    auto ib( o_.fsidx_.cbegin( anchor ) );

    for( WordTable::const_iterator wtb( wt.begin() + wmap[anchor] ),
                                   wte( wt.begin() + wmap[anchor + 1] );
         wtb != wte; ++wtb )
    {
        auto ww( wtb->wd.w.word );
        auto w( ww&PWMASK );

        if( GetBit( pa_[w/PWBITS], w%PWBITS ) == 0 )
        {
            ++jd_.pmisses;
        }
        else
        {
            ++jd_.phits;
            auto ii( ib + sidx_[MkSI( ww )] ),
                 iie( ib + sidx_[MkSI( ww ) + 1] );
            for( ; ii != iie && ii->wd.w.word < ww; ++ii );

            for( ; ii != iie && ii->wd.w.word == ww; ++ii )
            {
                SaveHit( *ii, *wtb );
            }
        }
    }
}

//------------------------------------------------------------------------------
inline void CFastSeeds::SeedSearchJob::MatchScan( uint32_t anchor )
{
    std::vector< ExtHashWord > & words( jd_.ewords );
    words.clear();
    auto const & wt( o_.wt_ );
    auto const & wmap( o_.wmap_ );

    for( WordTable::const_iterator wtb( wt.begin() + wmap[anchor] ),
                                   wte( wt.begin() + wmap[anchor + 1] );
         wtb != wte; ++wtb )
    {
        ExtHashWord hs { *wtb };
        auto w( hs.hw.wd.w.word&PWMASK );

        if( GetBit( pa_[w/PWBITS], w%PWBITS ) == 0 )
        {
            ++jd_.pmisses;
        }
        else
        {
            ++jd_.phits;
            words.push_back( hs );
        }
    }

    std::sort( words.begin(), words.end(),
               []( ExtHashWord const & x, ExtHashWord const & y )
               { return x.hw.wd.w.word < y.hw.wd.w.word; } );

    auto ib( o_.fsidx_.cbegin( anchor ) ),
         ie( o_.fsidx_.cend( anchor ) );
    std::vector< ExtHashWord >::const_iterator wb( words.begin() ),
                                               we( words.end() );

    while( ib != ie && wb != we )
    {
        for( ; ib != ie && ib->wd.w.word < wb->hw.wd.w.word; ++ib );
        for( ; wb != we && ib->wd.w.word > wb->hw.wd.w.word; ++wb );

        if( ib != ie && wb != we && ib->wd.w.word == wb->hw.wd.w.word )
        {
            auto iie( ib );
            for( ; iie != ie && ib->wd.w.word == iie->wd.w.word; ++iie );
            auto wwe( wb );
            for( ; wwe != we && wb->hw.wd.w.word == wwe->hw.wd.w.word; ++wwe );

            for( ; ib != iie; ++ib )
            {
                for( auto wc( wb ); wc != wwe; ++wc ) SaveHit( *ib, wc->hw );
            }

            wb = wwe;
        }
    }
}

//------------------------------------------------------------------------------
inline void CFastSeeds::SeedSearchJob::operator()()
{
    auto const & anchor_use_map( o_.anchor_use_map_ );

    while( true )
    {
        uint32_t task_idx( task_idx_.fetch_add( 1 ) );

        if( task_idx >= ANCHOR_TBL_SIZE )
        {
            break;
        }

        if( anchor_use_map[task_idx] )
        {
            scanner_.AcquireAnchor( task_idx, job_idx_ );
            auto ib( o_.fsidx_.cbegin( task_idx ) ),
                 ie( o_.fsidx_.cend( task_idx ) );

            if( ie - ib <= GOOD_ANCHOR_SIZE )
            {
                ++jd_.good_anchors;
                PrepareIndex( task_idx );
                MatchSubst( task_idx );
            }
            else
            {
                PrepareIndex( task_idx, false );
                MatchScan( task_idx );
            }
        }

        ph_.Increment();
    }

    scanner_.Finalize( job_idx_ );
}

//==============================================================================
struct CFastSeeds::ReadMarkingJob
{
    typedef std::atomic< uint32_t > JobIdx;
    typedef std::tuple< OrdId, float, size_t, uint8_t > MarkedRead;
    typedef std::vector< MarkedRead > MarkedReads;
    typedef CRefData::TRefOId TRefOId;
    typedef std::pair< TRefOId, Hit > ExtHit;
    typedef std::vector< ExtHit > ExtHits;

    ReadMarkingJob( CFastSeeds & o, size_t & thread_idx, JobIdx & job_idx,
                    std::vector< Hits > & hits,
                    std::vector< MarkedReads > & marked_reads,
                    CProgress::ProgressHandle & ph )
        : o_( o ), job_idx_( job_idx ), hits_( hits ),
          marked_reads_( marked_reads[thread_idx++] ),
          per_mate_( o.bctx_.GetSearchCtx().per_mate_marks ),
          ph_( ph )
    {}

    void operator()();
    void ProcessRead( ExtHits::const_iterator b, ExtHits::const_iterator e );
    float ReadIsCovered( ExtHits::const_iterator b, ExtHits::const_iterator e );

    CFastSeeds & o_;
    JobIdx & job_idx_;
    JobIdx curr_job_idx_;
    std::vector< Hits > & hits_;
    MarkedReads & marked_reads_;
    BitSet cover_;
    bool per_mate_ = false;
    CProgress::ProgressHandle & ph_;
};

//------------------------------------------------------------------------------
inline float CFastSeeds::ReadMarkingJob::ReadIsCovered(
        ExtHits::const_iterator b, ExtHits::const_iterator e )
{
    TSeqOff diag_delta( o_.bctx_.GetSearchCtx().max_diag_delta );
    assert( b != e );
    auto cth( o_.bctx_.GetSearchCtx().coverage_th );
    auto clenth( o_.bctx_.GetSearchCtx().covered_bases );
    TSeqOff d1( -diag_delta - cover_.size() - 1 ), dc( d1 );
    cth = cover_.size()*cth;

    for( auto c( b ); b != e; ++b )
    {
        auto d2( b->second.GetDiag() );

        if( dc == d1 && d2 > d1 )
        {
            c = b;
            dc = d2;
        }

        if( d2 - d1 > diag_delta )
        {
            if( cth <= 1.0f*cover_.count() &&
                (size_t)clenth <= cover_.count() )
            {
                return 1.0f*cover_.count()/cover_.size();
            }

            cover_.reset();
            b = c;
            d1 = dc;
        }

        auto const & h( b->second );

        for( ssize_t pos( h.readpos ), i( 0L );
             i < NMER_BASES && (size_t)pos < cover_.size(); ++pos, ++i )
        {
            cover_.set( pos );
        }
    }

    if( cth <= 1.0f*cover_.count() && (size_t)clenth <= cover_.count() )
    {
        return 1.0f*cover_.count()/cover_.size();
    }

    return -1.0f;
}

//------------------------------------------------------------------------------
inline void CFastSeeds::ReadMarkingJob::ProcessRead(
        ExtHits::const_iterator b, ExtHits::const_iterator e )
{
    assert( b != e );

    auto const & read( o_.GetReads()[b->second.read] );
    auto mate_idx( b->second.mate );
    cover_.resize( read.mates_[mate_idx].len );
    cover_.reset();

    while( b != e )
    {
        auto r( std::equal_range(
                    b, e, *b,
                    []( ExtHit const & x, ExtHit const & y )
                    {
                        return x.second.mate == y.second.mate ?
                               x.second.strand == y.second.strand ?
                               x.first < y.first :
                               x.second.strand < y.second.strand :
                               x.second.mate < y.second.mate;
                    } ) );

        if( mate_idx != b->second.mate )
        {
            mate_idx = b->second.mate;
            cover_.resize( read.mates_[mate_idx].len );
            cover_.reset();
        }

        float c( ReadIsCovered( r.first, r.second ) );

        if( c >= 0.0f )
        {
            marked_reads_.push_back( std::make_tuple(
                        b->second.read, c, cover_.count(), mate_idx ) );
            return;
        }

        b = r.second;
    }
}

//------------------------------------------------------------------------------
inline void CFastSeeds::ReadMarkingJob::operator()()
{
    ExtHits ehits;

    while( true )
    {
        curr_job_idx_ = job_idx_.fetch_add( 1 );

        if( curr_job_idx_ >= hits_.size() )
        {
            break;
        }

        // convert hits to extended hits
        //
        auto & hits( hits_[curr_job_idx_] );
        std::sort( hits.begin(), hits.end(),
                   []( Hit const & x, Hit const & y )
                   { return x.refpos < y.refpos; } );
        auto const & refs( o_.GetRefs() );
        Hits::const_iterator hi( hits.begin() ),
                             hie( hits.end() ),
                             hc( hi );

        for( TRefOId i( 0 ), ie( refs.GetSize() ); i < ie; ++i )
        {
            auto refpos_range( o_.fsidx_.GetAbsPosRange( i ) );
            for( ; hc != hie && hc->refpos < refpos_range.second; ++hc );

            for( ; hi != hc; ++hi )
            {
                ehits.push_back( std::make_pair( i, *hi ) );
                ehits.back().second.refpos -= refpos_range.first;
            }
        }

        // free memory used by original hits
        //
        { Hits t; t.swap( hits ); }

        // sort extended hits by read, ref, mate, strand, diag
        //
        std::sort( ehits.begin(), ehits.end(),
                   []( ExtHit const & x, ExtHit const & y )
                   {
                        return x.second.read == y.second.read ?
                               x.second.mate == y.second.mate ?
                               x.second.strand == y.second.strand ?
                               x.first == y.first ?
                               x.second.GetDiag() < y.second.GetDiag() :
                               x.first < y.first :
                               x.second.strand < y.second.strand :
                               x.second.mate < y.second.mate :
                               x.second.read < y.second.read;
                   } );

        // process extended hits on a per-read basis
        //
        for( auto ehi( ehits.cbegin() ), ehie( ehits.cend() ); ehi != ehie; )
        {
            if( per_mate_ )
            {
                auto r( std::equal_range(
                            ehi, ehie, *ehi,
                            []( ExtHit const & x, ExtHit const & y )
                            {
                                return x.second.read == y.second.read ?
                                       x.second.mate < y.second.mate :
                                       x.second.read < y.second.read;
                            } ) );
                ProcessRead( r.first, r.second );
                ehi = r.second;
            }
            else
            {
                auto r( std::equal_range(
                            ehi, ehie, *ehi,
                            []( ExtHit const & x, ExtHit const & y )
                            { return x.second.read < y.second.read; } ) );
                ProcessRead( r.first, r.second );
                ehi = r.second;
            }
        }

        ehits.clear();
        ph_.Increment();
    }
}

//------------------------------------------------------------------------------
inline void CFastSeeds::ComputeSeeds()
{
    auto & ctx( bctx_.GetSearchCtx() );
    auto const & refs( *ctx.refs );
    std::vector< SeedSearchJobData > jd( ctx.n_threads );

    if( refs.GetSize() > 0 )
    {
        IndexScanner scanner( ctx, fsidx_ );
        StopWatch w( ctx.logger_ );
        std::atomic< uint32_t > task_idx( 0 );
        size_t job_idx( 0 );
        CProgress p( "looking for initial seeds", "anchors",
                     ctx.progress_flags_ );
        auto ph( p.GetTop() );
        ph.SetTotal( ANCHOR_TBL_SIZE );
        CTaskArray< SeedSearchJob > jobs( 
                ctx.n_threads, *this, scanner,
                task_idx, &jd[0], job_idx, ph );
        p.Start();
        jobs.Start();
        p.Stop();
        M_INFO( ctx.logger_, "hits computed" );
    }

    SeedSearchJobData jdsum;

    for( auto const & j : jd )
    {
        jdsum += j;
    }

    M_INFO( ctx.logger_, "hits: " << jdsum.hits <<
                          "; phits: " << jdsum.phits <<
                          "; pmisses: " << jdsum.pmisses <<
                          "; good anchors: " << jdsum.good_anchors );

    // free memory used by reference and read index
    //
    // fsidx_.Unload();
    { WordTable t; wt_.swap( t ); }
    { WordMap t; t.swap( wmap_ ); }

    // hit filtering
    //
    {
        size_t reads_per_job( bctx_.GetJobSize() );
        assert( reads_per_job > 0 );
        size_t n_jobs( bctx_.GetReads().GetNReads()/reads_per_job + 1 );
        std::vector< Hits > filter_job_hits( n_jobs );

        // distribute hits accross jobs
        //
        {
            StopWatch w( ctx.logger_ );

            {
                size_t total_hits( 0 );

                for( auto const & job_data : jd )
                {
                    total_hits += job_data.results.size();
                }

                M_INFO( ctx.logger_, "TOTAL HITS: " << total_hits );
            }

            M_INFO( ctx.logger_, "distributing hits..." );

            for( auto & job_data : jd )
            {
                for( auto const & h : job_data.results )
                {
                    size_t job_num( h.read/reads_per_job );
                    filter_job_hits[job_num].push_back( h );
                }

                // clean up search job results
                //
                {
                    Hits t;
                    t.swap( job_data.results );
                }
            }

            M_INFO( ctx.logger_, "distributed hits for filtering" );
        }

        // hit fitering jobs
        //
        {
            StopWatch w( ctx.logger_ );

            std::atomic< uint32_t > job_idx( 0 );
            CProgress p( "filtering initial seeds", "tasks",
                         ctx.progress_flags_ );
            auto ph( p.GetTop() );
            ph.SetTotal( n_jobs );

            {
                bool per_mate( bctx_.GetSearchCtx().per_mate_marks );
                std::vector< ReadMarkingJob::MarkedReads > marked_reads(
                        ctx.n_threads );
                size_t thread_idx( 0 );
                CTaskArray< ReadMarkingJob > jobs(
                        ctx.n_threads, *this, thread_idx, job_idx,
                        filter_job_hits, marked_reads, ph );
                p.Start();
                jobs.Start();
                ReadMarkingJob::MarkedReads result;

                for( auto & mr : marked_reads )
                {
                    std::copy( mr.begin(), mr.end(),
                               std::back_inserter( result ) );
                    ReadMarkingJob::MarkedReads t;
                    t.swap( mr );
                }

                std::sort( result.begin(), result.end(),
                           []( ReadMarkingJob::MarkedRead const & x,
                               ReadMarkingJob::MarkedRead const & y )
                           {
                                return std::get< 0 >( x )
                                            == std::get< 0 >( y ) ?
                                       std::get< 3 >( x )
                                            < std::get< 3 >( y ) :
                                       std::get< 0 >( x )
                                            < std::get< 0 >( y );
                           } );

                if( per_mate )
                {
                    M_INFO( ctx.logger_,
                            "number of marked mates: " << result.size() );
                }
                else
                {
                    M_INFO( ctx.logger_,
                            "number of marked reads: " << result.size() );
                }

                auto const & reads( bctx_.GetReads() );
                auto & os( bctx_.GetSearchCtx().GetOutStream() );
                std::vector< char > seq;
                typedef SEQ_NS::CRecoder< eIUPACNA, eNCBI2NA > R;
                typedef CReadData::SeqConstView SeqView;
                bctx_.GetSearchCtx().n_mapped_reads += result.size();

                for( auto i : result )
                {
                    auto const & read( reads[std::get< 0 >( i )] );
                    std::string id( reads.GetReadId( std::get< 0 >( i ) ) );

                    if( per_mate )
                    {
                        auto mate_idx( std::get< 3 >( i ) );
                        auto len( read.mates_[mate_idx].len );
                        os << ">" << id;

                        if( read.ReadIsPaired() )
                        {
                            os << '.' << mate_idx + 1;
                        }

                        os << " " << std::get< 1 >( i )
                           << " " << std::get< 2 >( i )
                           << '\n';
                        seq.clear();
                        seq.resize( len + 1, 0 );
                        auto sd( reads.GetSeqData( std::get< 0 >( i ),
                                                   ToMate( mate_idx ),
                                                   eFWD ) ),
                             md( reads.GetMaskData( std::get< 0 >( i ),
                                                    ToMate( mate_idx ),
                                                    eFWD ) );

                        for( TSeqLen j( 0 ); j < len; ++j )
                        {
                            if( md.GetLetter( j ) == SeqView::Code::MASK_BASE )
                            {
                                seq[j] = 'N';
                            }
                            else
                            {
                                seq[j] = R::Recode( sd.GetLetter( j ) );
                            }
                        }

                        os << &seq[0] << '\n';
                    }
                    else
                    {
                        if( read.mates_[0].len > 0 )
                        {
                            auto len( read.mates_[0].len );
                            os << ">" << id;

                            if( read.ReadIsPaired() )
                            {
                                os << ".1";
                            }

                            os << " " << std::get< 1 >( i )
                               << " " << std::get< 2 >( i );
                            os << '\n';
                            seq.clear();
                            seq.resize( len + 1, 0 );
                            auto sd( reads.GetSeqData(
                                        std::get< 0 >( i ), eFIRST, eFWD ) ),
                                 md( reads.GetMaskData(
                                        std::get< 0 >( i ), eFIRST, eFWD ) );

                            for( TSeqLen j( 0 ); j < len; ++j )
                            {
                                if( md.GetLetter( j ) ==
                                        SeqView::Code::MASK_BASE )
                                {
                                    seq[j] = 'N';
                                }
                                else
                                {
                                    seq[j] = R::Recode( sd.GetLetter( j ) );
                                }
                            }

                            os << &seq[0] << '\n';
                        }

                        if( read.mates_[1].len > 0 )
                        {
                            auto len( read.mates_[1].len );
                            os << ">" << id;

                            if( read.ReadIsPaired() )
                            {
                                os << ".2";
                            }

                            os << " " << std::get< 1 >( i )
                               << " " << std::get< 2 >( i );
                            os << '\n';
                            seq.clear();
                            seq.resize( len + 1, 0 );
                            auto sd( reads.GetSeqData(
                                        std::get< 0 >( i ), eSECOND, eFWD ) ),
                                md( reads.GetMaskData(
                                        std::get< 0 >( i ), eSECOND, eFWD ) );

                            for( TSeqLen j( 0 ); j < len; ++j )
                            {
                                if( md.GetLetter( j ) ==
                                        SeqView::Code::MASK_BASE )
                                {
                                    seq[j] = 'N';
                                }
                                else
                                {
                                    seq[j] = R::Recode( sd.GetLetter( j ) );
                                }
                            }

                            os << &seq[0] << '\n';
                        }
                    }
                }
            }

            p.Stop();
        }
    }
}

//------------------------------------------------------------------------------
void CFastSeeds::CreateWordTable()
{
    static size_t const DEFAULT_READS_PER_JOB = 10000ULL;
    auto & ctx( bctx_.GetSearchCtx() );
    auto const & reads( bctx_.GetReads() );
    M_INFO( ctx.logger_, "generating read index" );
    auto n_threads( ctx.n_threads );
    size_t n_reads_per_job( DEFAULT_READS_PER_JOB );

    // compute number of reads per job
    //
    {
        size_t n_reads( reads.GetNReads() );
        n_reads_per_job = std::min( n_reads_per_job, n_reads/n_threads + 1U );
        M_INFO( ctx.logger_,
                "using " << n_reads_per_job << " reads per job" );
    }

    // count words in reads
    //
    WordCountingJob::WCJData wcj_data( n_threads );

    {
        StopWatch w( ctx.logger_ );
        size_t job_idx( 0 );
        std::atomic< uint32_t > task_idx( 0 );
        CProgress p( "counting read words", "tasks", ctx.progress_flags_ );
        auto ph( p.GetTop() );
        ph.SetTotal( reads.GetNReads() );
        CTaskArray< WordCountingJob > jobs(
                n_threads, *this, wcj_data, task_idx, n_reads_per_job, 
                job_idx, ph );
        p.Start();
        jobs.Start();
        p.Stop();
        M_INFO( ctx.logger_, "word counting complete" );
    }

    // count number of words per anchor value in total and per thread
    //
    {
        StopWatch w( ctx.logger_ );
        BitSet wa_map( WMAP_SIZE, false );
        uint32_t wmap_off( 0 );

        for( size_t i( 0 ); i < WMAP_SIZE; ++i )
        {
            wmap_[i] = wmap_off;

            for( size_t j( 0 ); j < wcj_data.size(); ++j )
            {
                auto t( wcj_data[j].wmap[i] );
                wcj_data[j].wmap[i] = wmap_off;
                wmap_off += t;
            }

            if( wmap_off > wmap_[i] )
            {
                wa_map[i] = true;
            }
        }

        M_INFO( ctx.logger_, wmap_off << " words in the word table" );
        wt_.resize( wmap_off );
        M_INFO( ctx.logger_, "word counts aggregation complete" );

        for( BitSet::size_type wa( wa_map.find_first() );
             wa != BitSet::npos; wa = wa_map.find_next( wa ) )
        {
            UpdateAnchorUseMap( (uint32_t)wa );
        }
    }

    // fill in word tables
    //
    {
        StopWatch w( ctx.logger_ );
        size_t job_idx( 0 );
        CProgress p( "generating read word index", "tasks",
                     ctx.progress_flags_ );
        auto ph( p.GetTop().Split( n_threads ) );
        CTaskArray< WordTableGeneratorJob > jobs(
                n_threads, *this, wcj_data, n_reads_per_job, job_idx, ph );
        p.Start();
        jobs.Start();
        p.Stop();
        M_INFO( ctx.logger_, "read word index generation complete" );
    }

    // sort word table
    //
    {
        StopWatch w( ctx.logger_ );
        std::atomic< uint32_t > task_idx( 0 );
        CProgress p( "sorting read word index", "anchors",
                     ctx.progress_flags_ );
        auto ph( p.GetTop() );
        ph.SetTotal( WMAP_SIZE - 1 );
        CTaskArray< WordTableSortJob > jobs( n_threads, *this, task_idx, ph );
        p.Start();
        jobs.Start();
        p.Stop();
        M_INFO( ctx.logger_, "sorting of word table complete" );
    }

#ifdef TRACE_WT
    PrintWordTable();
#endif
}

//------------------------------------------------------------------------------
inline void CFastSeeds::UpdateAnchorUseMap( uint32_t wa )
{
    anchor_use_map_[wa] = true;
}

//------------------------------------------------------------------------------
void CFastSeeds::Run()
{
    CreateWordTable();
    ComputeSeeds();
}

READFINDER_NS_END

