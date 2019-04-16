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
CFastSeeds::CFastSeeds( CBatch & bctx, bool seeder_mode )
    : bctx_( bctx ),
      fsidx_( CommonCtxP( new CCommonContext(
                      static_cast< CCommonContext & >( 
                          bctx.GetSearchCtx() ) ) ),
              *bctx.GetSearchCtx().refs ),
      anchor_use_map_( ANCHOR_TBL_SIZE, false ),
      seeder_mode_( seeder_mode )
{
    bctx_.GetSearchCtx().refs->LoadAll();
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
class CFastSeeds::HashWordSource
{
public:

    HashWordSource( TWord const * seq_data, TSeqLen len, TSeqOff off,
                    TSeqOff strode = WORD_STRIDE );
    void operator++();

    uint32_t GetAnchor() const { return data_.f.anchor; }
    uint32_t GetWord() const { return data_.f.word; }
    uint32_t GetSfx() const { return data_.f.sfx; }
    TSeqOff GetHashOff() const { return off_ + 1 - off_adj_; }
    operator bool() const { return !done_; }

private:

    struct
    {
        union
        {
            struct
            {
                uint64_t word   : WORD_BITS;
                uint64_t anchor : ANCHOR_BITS + LB;
                uint64_t sfx    : SFX_BITS;
            } f;

            TWord w;
        };
    } data_;

    static_assert( sizeof( data_ ) == 8, "" );

    TWord nw_;
    TWord const * seq_data_;
    TSeqLen len_;
    TSeqOff off_,
            off_adj_,
            stride_;
    bool done_ = true;
};

//------------------------------------------------------------------------------
inline CFastSeeds::HashWordSource::HashWordSource(
        TWord const * seq_data, TSeqLen len, TSeqOff off, TSeqOff stride )
    : nw_( 0 ), seq_data_( seq_data ), len_( off + len ),
      off_( off ), off_adj_( off ), stride_( stride )
{
    if( off_ >= len_ )
    {
        return;
    }

    done_ = false;
    data_.w = *seq_data_++;
    nw_ = *seq_data_++;

    auto shift( off*LB );
    TWord mask( (1ULL<<shift) - 1 );
    data_.w >>= shift;
    data_.w += ((nw_&mask)<<(WL*LB - shift));
    nw_ >>= shift;
}

//------------------------------------------------------------------------------
inline void CFastSeeds::HashWordSource::operator++()
{
    TSeqOff woff( off_%WL );
    off_ += stride_;
    
    if( off_ >= len_ )
    {
        done_ = true;
        return;
    }

    woff += stride_;
    data_.w >>= stride_*LB;
    data_.w += (nw_<<((WL - stride_)*LB));

    if( woff >= WL )
    {
        woff -= WL;
        nw_ = *seq_data_++;

        if( woff > 0 )
        {
            data_.w += (nw_<<((WL - woff)*LB));
            nw_ >>= woff*LB;
        }
    }
    else
    {
        nw_ >>= stride_*LB;
    }
}

//==============================================================================
class CFastSeeds::HashMaskSource
{
public:

    HashMaskSource( TWord const * mask_data, TSeqLen len, TSeqOff off,
                    TSeqOff stride = WORD_STRIDE );
    void operator++();

    uint64_t GetNMer() const { return data_.f.nmer; }
    uint64_t GetSfx() const { return data_.f.sfx; }
    operator bool() const { return !done_; }

private:

    struct
    {
        union
        {
            struct
            {
                uint64_t pfx  : LB;
                uint64_t nmer : NMER_BITS;
                uint64_t sfx  : SFX_BITS;
            } f;

            TWord w;
        };
    } data_;

    static_assert( sizeof( data_ ) == 8, "" );

    TWord nw_ = 0;
    TWord const * mask_data_;
    TSeqLen len_;
    TSeqOff off_ = 0,
            stride_;
    bool done_ = true;
};

//------------------------------------------------------------------------------
inline CFastSeeds::HashMaskSource::HashMaskSource(
        TWord const * mask_data, TSeqLen len, TSeqOff off, TSeqOff stride )
    : mask_data_( mask_data ), len_( off + len ), off_( off ),
      stride_( stride )
{
    if( off_ >= len )
    {
        return;
    }

    done_ = false;
    data_.w = *mask_data_++;
    nw_ = *mask_data_++;

    auto shift( off*LB );
    TWord mask( (1ULL<<shift) - 1 );
    data_.w >>= shift;
    data_.w += ((nw_&mask)<<(WL*LB - shift));
    nw_ >>= shift;
}

//------------------------------------------------------------------------------
inline void CFastSeeds::HashMaskSource::operator++()
{
    TSeqOff woff( off_%WL );
    off_ += stride_;
    
    if( off_ >= len_ )
    {
        done_ = true;
        return;
    }

    woff += stride_;
    data_.w >>= stride_*LB;
    data_.w += (nw_<<((WL - stride_)*LB));

    if( woff >= WL )
    {
        woff -= WL;
        nw_ = *mask_data_++;

        if( woff > 0 )
        {
            data_.w += (nw_<<((WL - woff)*LB));
            nw_ >>= woff*LB;
        }
    }
    else
    {
        nw_ >>= stride_*LB;
    }
}

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
    auto stride( o_.seeder_mode_ ? 1 : WORD_STRIDE );

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
                HashWordSource hws(
                        sd.GetBuf(), sd.size(), sd.GetStart(), stride );
                HashMaskSource hms(
                        md.GetBuf(), md.size(), md.GetStart(), stride );

                while( hws )
                {
                    assert( hms );

                    if( GetNSet( hms.GetNMer() ) <= LB )
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
                            sd.GetBuf(), sd.size(), sd.GetStart(), stride );
                    HashMaskSource hms( 
                            md.GetBuf(), md.size(), md.GetStart(), stride );

                    while( hws )
                    {
                        assert( hms );

                        if( GetNSet( hms.GetNMer() ) <= LB )
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
    auto stride( o_.seeder_mode_ ? 1 : WORD_STRIDE );

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
                HashWordSource hws(
                        sd.GetBuf(), sd.size(), sd.GetStart(), stride );
                HashMaskSource hms(
                        md.GetBuf(), md.size(), md.GetStart(), stride );

                while( hws )
                {
                    assert( hms );

                    if( GetNSet( hms.GetNMer() ) <= LB )
                    {
                        auto & w( wt[wmap[hws.GetAnchor()]++] );
                        w.wd.w.word = hws.GetWord();
                        w.wd.w.sfx = hws.GetSfx();
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
                            sd.GetBuf(), sd.size(), sd.GetStart(), stride );
                    HashMaskSource hms(
                            md.GetBuf(), md.size(), md.GetStart(), stride );

                    while( hws )
                    {
                        assert( hms );

                        if( GetNSet( hms.GetNMer() ) <= LB )
                        {
                            auto & w( wt[wmap[hws.GetAnchor()]++] );
                            w.wd.w.word = hws.GetWord();
                            w.wd.w.sfx = hws.GetSfx();
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
    void MatchIns( uint32_t anchor );
    void MatchSubst( uint32_t anchor );
    void MatchDel( uint32_t anchor );
    void MatchScan( uint32_t anchor );
    void MatchPfx_1( uint32_t anchor );
    void MatchPfx_2( uint32_t anchor );
    void MatchPfx_2I( uint32_t anchor );
    void MatchPfx_2E( uint32_t anchor );
    void MatchPfx_2D( uint32_t anchor );

    void SaveHit( IndexEntry const & iw,
                  HashWord const & hw,
                  bool exact, int ins, int del,
                  TReadOff posadj = 0 );

    static uint32_t MkSI( uint32_t w )
    {
        return (w&SIMASK)>>SISHIFT;
    }

    static bool CheckIns( uint32_t iw, uint32_t ww, size_t & n_init_matches );
    static bool CheckDel( uint32_t iw, uint32_t ww, size_t & n_init_matches );
    static bool CheckSubst(
            uint32_t iw, uint32_t ww, size_t & n_init_matches );
    static int Check1Err( uint32_t w1, uint32_t w2 );

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
inline bool CFastSeeds::SeedSearchJob::CheckIns(
        uint32_t iw, uint32_t ww, size_t & n_init_matches )
{
    uint32_t xw( iw^ww );
    n_init_matches = (GetNFirstClear( xw )>>1);
    size_t nc( n_init_matches<<1 );
    iw >>= nc;
    ww >>= (nc + LB);
    xw = (iw^ww);
    return xw == 0;
}

//------------------------------------------------------------------------------
inline bool CFastSeeds::SeedSearchJob::CheckDel(
        uint32_t iw, uint32_t ww, size_t & n_init_matches )
{
    uint32_t xw( iw^ww );
    n_init_matches = (GetNFirstClear( xw )>>1);
    size_t nc( n_init_matches<<1 );
    iw >>= (nc + LB);
    ww >>= nc;
    xw = (iw^ww);
    return xw == 0;
}

//------------------------------------------------------------------------------
inline bool CFastSeeds::SeedSearchJob::CheckSubst(
        uint32_t iw, uint32_t ww, size_t & n_init_matches )
{
    uint32_t xw( iw^ww );
    xw = ((xw|(xw>>1))&0x55555555UL);
    n_init_matches = (GetNFirstClear( xw )>>1);
    return GetNSet( xw ) <= 1;
}

//------------------------------------------------------------------------------
inline int CFastSeeds::SeedSearchJob::Check1Err( uint32_t w1, uint32_t w2 )
{
    static uint32_t const MASK = 0x3FFUL;

    uint32_t xw( w1^w2 );
    xw = ((xw|(xw>>1))&0x55555555UL);

    if( GetNSet( xw ) <= 1 )
    {
        return 0;
    }

    size_t nc( GetNFirstClear( xw ) );
    w1 >>= nc;
    w2 >>= nc;
    xw = (w1>>LB);
    xw ^= (w2&MASK);

    if( xw == 0 )
    {
        return 2;
    }

    xw = (w2>>LB);
    xw ^= (w1&MASK);

    if( xw == 0 )
    {
        return 1;
    }

    return 3;
}

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
        IndexEntry const & iw, HashWord const & hw, 
        bool exact, int ins, int del, TReadOff posadj )
{
    ++jd_.hits;
    TReadOff readoff( hw.hashoff + posadj + del - ins );
    Hit h { iw.pos + posadj, hw.readid, readoff,
            (uint8_t)(hw.strand - 1), (uint8_t)(hw.mate - 1), exact };
    
    if( o_.seeder_mode_ )
    {
        h.mid = (posadj != 0);
    }

    jd_.results.push_back( h );
}

//------------------------------------------------------------------------------
inline void CFastSeeds::SeedSearchJob::MatchIns( uint32_t anchor )
{
    auto & ctx( o_.bctx_.GetSearchCtx() );

    if( ctx.exact_seeds )
    {
        return;
    }

    auto const & wt( o_.wt_ );
    auto const & wmap( o_.wmap_ );
    auto const & al( o_.atbl_[anchor] );
    auto const & d( al.data );

    auto ib( o_.fsidx_.cbegin( anchor ) );

    for( size_t i( 0 ); i < al.n_ins; ++i )
    {
        auto wa( d[i].f.wanchor );

        for( WordTable::const_iterator wtb( wt.begin() + wmap[wa] ),
                                       wte( wt.begin() + wmap[wa + 1] );
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
                    if( !o_.seeder_mode_ && ii->wd.w.repeat )
                    {
                        break;
                    }

                    SaveHit( *ii, *wtb, false, 1, 0 );
                }
            }
        }
    }
}

//------------------------------------------------------------------------------
inline void CFastSeeds::SeedSearchJob::MatchSubst( uint32_t anchor )
{
    auto & ctx( o_.bctx_.GetSearchCtx() );
    auto const & wt( o_.wt_ );
    auto const & wmap( o_.wmap_ );
    auto const & al( o_.atbl_[anchor] );
    auto const & d( al.data );
    bool es( ctx.exact_seeds );

    auto ib( o_.fsidx_.cbegin( anchor ) );

    for( size_t i( al.n_ins ), ie( al.n_ins + al.n_subst ); i < ie; ++i )
    {
        bool exact( d[i].f.exact );

        if( es && !exact )
        {
            continue;
        }

        auto wa( d[i].f.wanchor );

        for( WordTable::const_iterator wtb( wt.begin() + wmap[wa] ),
                                       wte( wt.begin() + wmap[wa + 1] );
             wtb != wte; ++wtb )
        {
            auto ww( wtb->wd.w.word );
            ww >>= LB;
            ww += ((wa&LMASK)<<(WORD_BITS - LB));
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
                    if( !o_.seeder_mode_ )
                    {
                        if( (!exact && ii->wd.w.repeat) || ii->wd.w.erepeat )
                        {
                            break;
                        }
                    }

                    SaveHit( *ii, *wtb, exact, 0, 0 );
                }
            }
        }
    }
}

//------------------------------------------------------------------------------
inline void CFastSeeds::SeedSearchJob::MatchDel( uint32_t anchor )
{
    auto & ctx( o_.bctx_.GetSearchCtx() );

    if( ctx.exact_seeds )
    {
        return;
    }

    static auto const MASK = (LMASK<<LB) + LMASK;

    auto const & wt( o_.wt_ );
    auto const & wmap( o_.wmap_ );
    auto const & al( o_.atbl_[anchor] );
    auto const & d( al.data );

    auto ib( o_.fsidx_.cbegin( anchor ) );

    for( size_t i( al.n_ins + al.n_subst ); i < al.len; ++i )
    {
        auto wa( d[i].f.wanchor );

        for( WordTable::const_iterator wtb( wt.begin() + wmap[wa] ),
                                       wte( wt.begin() + wmap[wa + 1] );
             wtb != wte; ++wtb )
        {
            auto ww( wtb->wd.w.word );
            ww >>= (LB + LB);
            ww += ((wa&MASK)<<(WORD_BITS - LB - LB));
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
                    if( !o_.seeder_mode_ && ii->wd.w.repeat )
                    {
                        break;
                    }

                    SaveHit( *ii, *wtb, false, 0, 1 );
                }
            }
        }
    }
}

//------------------------------------------------------------------------------
inline void CFastSeeds::SeedSearchJob::MatchScan( uint32_t anchor )
{
    static auto const MASK = (LMASK<<LB) + LMASK;

    auto & ctx( o_.bctx_.GetSearchCtx() );
    std::vector< ExtHashWord > & words( jd_.ewords );
    words.clear();
    auto const & wt( o_.wt_ );
    auto const & wmap( o_.wmap_ );
    auto const & al( o_.atbl_[anchor] );
    auto const & d( al.data );
    bool es( ctx.exact_seeds );

    for( size_t i( 0 ); i < al.len; ++i )
    {
        bool exact( d[i].f.exact );

        if( es && !exact )
        {
            continue;
        }

        auto wa( d[i].f.wanchor );

        for( WordTable::const_iterator wtb( wt.begin() + wmap[wa] ),
                                       wte( wt.begin() + wmap[wa + 1] );
             wtb != wte; ++wtb )
        {
            ExtHashWord hs { *wtb, exact };

            if( d[i].f.ins != 0 )
            {
                --hs.hw.hashoff;
            }
            else if( d[i].f.del != 0 )
            {
                ++hs.hw.hashoff;
                auto w( hs.hw.wd.w.word );
                w >>= (LB + LB);
                w += ((wa&MASK)<<(WORD_BITS - LB - LB));
                hs.hw.wd.w.word = w;
            }
            else
            {
                auto w( hs.hw.wd.w.word );
                w >>= LB;
                w += ((wa&LMASK)<<(WORD_BITS - LB));
                hs.hw.wd.w.word = w;
            }

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

            if( !o_.seeder_mode_ && ib->wd.w.erepeat )
            {
                ib = iie;
                wb = wwe;
                continue;
            }

            bool repeat( !o_.seeder_mode_ && ib->wd.w.repeat );

            for( ; ib != iie; ++ib )
            {
                for( auto wc( wb ); wc != wwe; ++wc )
                {
                    if( wc->exact || !repeat )
                    {
                        SaveHit( *ib, wc->hw, wc->exact, 0, 0 );
                    }
                }
            }

            wb = wwe;
        }
    }
}

//------------------------------------------------------------------------------
inline void CFastSeeds::SeedSearchJob::MatchPfx_2I( uint32_t anchor )
{
    static uint32_t const LWMASK = 0xFFFUL;
    static uint32_t const HWMASK = 0xFFF000UL;
    static size_t const WSHIFT = 12;

    auto const & idxwords( jd_.idxwords );
    auto & words( jd_.words );
    auto const & wt( o_.wt_ );
    auto const & wmap( o_.wmap_ );

    for( uint32_t l( 0 ); l < 3; ++l )
    {
        words.clear();
        auto wa( (anchor<<LB) + l );
        WordTable::const_iterator wtb( wt.begin() + wmap[wa] ),
                                  wte( wt.begin() + wmap[wa + 1] );
        std::copy( wtb, wte, std::back_inserter( words ) );

        std::sort( words.begin(), words.end(),
                   []( HashWord const & x, HashWord const & y )
                   { return (x.wd.w.word&LWMASK) < (y.wd.w.word&LWMASK); } );

        std::vector< IndexEntry >::const_iterator ib( idxwords.begin() ),
                                                  ie( idxwords.end() ),
                                                  iie, ic;
        std::vector< HashWord >::const_iterator wb( words.begin() ),
                                                we( words.end() ),
                                                wwe, wc;

        while( ib != ie && wb != we )
        {
            uint32_t ikey( ib->wd.w.word&LWMASK ),
                     wkey( wb->wd.w.word&LWMASK );

            if( ikey < wkey )
            {
                for( ; ib != ie && (ib->wd.w.word&LWMASK) == ikey; ++ib );
                continue;
            }

            if( ikey > wkey )
            {
                for( ; wb != we && (wb->wd.w.word&LWMASK) == wkey; ++wb );
                continue;
            }

            for( iie = ib;
                 iie != ie && (iie->wd.w.word&LWMASK) == ikey;
                 ++iie );
            for( wwe = wb;
                 wwe != we && (wwe->wd.w.word&LWMASK) == wkey;
                 ++wwe );

            for( ic = ib; ic != iie; ++ic )
            {
                uint32_t iword( (ic->wd.w.word&HWMASK)>>WSHIFT );

                for( wc = wb; wc != wwe; ++wc )
                {
                    uint32_t wword( (wc->wd.w.word&HWMASK)>>WSHIFT );
                    wword += (l<<(WORD_BITS - WSHIFT));
                    size_t n_init_matches( 0 );

                    if( CheckIns( iword, wword, n_init_matches ) )
                    {
                        if( n_init_matches >= 4 )
                        {
                            SaveHit( *ic, *wc, false, 1, 0 );
                        }
                        else
                        {
                            SaveHit( *ic, *wc, false, 0, 0, WORD_BASES );
                        }
                    }
                }
            }

            ib = iie;
            wb = wwe;
        }
    }
}

//------------------------------------------------------------------------------
inline void CFastSeeds::SeedSearchJob::MatchPfx_2E( uint32_t anchor )
{
    static uint32_t const LWMASK = 0xFFFUL;
    static uint32_t const HWMASK = 0xFFF000UL;
    static size_t const WSHIFT = 12;

    auto const & idxwords( jd_.idxwords );
    auto & words( jd_.words );
    words.clear();
    auto const & wt( o_.wt_ );
    auto const & wmap( o_.wmap_ );

    for( uint32_t l( 0 ); l < 3; ++l )
    {
        auto wa( (anchor<<LB) + l );
        WordTable::const_iterator wtb( wt.begin() + wmap[wa] ),
                                  wte( wt.begin() + wmap[wa + 1] );

        for( ; wtb != wte; ++wtb )
        {
            words.push_back( *wtb );
            words.back().wd.w.word = 
                (wtb->wd.w.word>>LB) + (l<<(WORD_BITS - LB));
        }
    }

    std::sort( words.begin(), words.end(),
               []( HashWord const & x, HashWord const & y )
               { return (x.wd.w.word&LWMASK) < (y.wd.w.word&LWMASK); } );

    std::vector< IndexEntry >::const_iterator ib( idxwords.begin() ),
                                              ie( idxwords.end() ),
                                              iie, ic;
    std::vector< HashWord >::const_iterator wb( words.begin() ),
                                            we( words.end() ),
                                            wwe, wc;

    while( ib != ie && wb != we )
    {
        uint32_t ikey( ib->wd.w.word&LWMASK ),
                 wkey( wb->wd.w.word&LWMASK );

        if( ikey < wkey )
        {
            for( ; ib != ie && (ib->wd.w.word&LWMASK) == ikey; ++ib );
            continue;
        }

        if( ikey > wkey )
        {
            for( ; wb != we && (wb->wd.w.word&LWMASK) == wkey; ++wb );
            continue;
        }

        for( iie = ib; iie != ie && (iie->wd.w.word&LWMASK) == ikey; ++iie );
        for( wwe = wb; wwe != we && (wwe->wd.w.word&LWMASK) == wkey; ++wwe );

        for( ic = ib; ic != iie; ++ic )
        {
            uint32_t iword( (ic->wd.w.word&HWMASK)>>WSHIFT );

            for( wc = wb; wc != wwe; ++wc )
            {
                uint32_t wword( (wc->wd.w.word&HWMASK)>>WSHIFT );
                size_t n_init_matches( 0 );

                if( CheckSubst( iword, wword, n_init_matches ) )
                {
                    if( n_init_matches >= 4 )
                    {
                        SaveHit( *ic, *wc, false, 0, 0 );
                    }
                    else
                    {
                        SaveHit( *ic, *wc, false, 0, 0, WORD_BASES );
                    }
                }
            }
        }

        ib = iie;
        wb = wwe;
    }
}

//------------------------------------------------------------------------------
inline void CFastSeeds::SeedSearchJob::MatchPfx_2D( uint32_t anchor )
{
    static uint32_t const LWMASK = 0xFFFUL;
    static uint32_t const HWMASK = 0xFFF000UL;
    static size_t const WSHIFT = 12;

    auto const & idxwords( jd_.idxwords );
    auto & words( jd_.words );
    words.clear();
    auto const & wt( o_.wt_ );
    auto const & wmap( o_.wmap_ );

    for( uint32_t l( 0 ); l < 3; ++l )
    {
        auto wa( (anchor<<LB) + l );
        WordTable::const_iterator wtb( wt.begin() + wmap[wa] ),
                                  wte( wt.begin() + wmap[wa + 1] );

        for( ; wtb != wte; ++wtb )
        {
            words.push_back( *wtb );
            words.back().wd.w.word = 
                (wtb->wd.w.word>>(LB + LB)) + (l<<(WORD_BITS - LB - LB));
        }
    }

    std::sort( words.begin(), words.end(),
               []( HashWord const & x, HashWord const & y )
               { return (x.wd.w.word&LWMASK) < (y.wd.w.word&LWMASK); } );

    std::vector< IndexEntry >::const_iterator ib( idxwords.begin() ),
                                              ie( idxwords.end() ),
                                              iie, ic;
    std::vector< HashWord >::const_iterator wb( words.begin() ),
                                            we( words.end() ),
                                            wwe, wc;

    while( ib != ie && wb != we )
    {
        uint32_t ikey( ib->wd.w.word&LWMASK ),
                 wkey( wb->wd.w.word&LWMASK );

        if( ikey < wkey )
        {
            for( ; ib != ie && (ib->wd.w.word&LWMASK) == ikey; ++ib );
            continue;
        }

        if( ikey > wkey )
        {
            for( ; wb != we && (wb->wd.w.word&LWMASK) == wkey; ++wb );
            continue;
        }

        for( iie = ib; iie != ie && (iie->wd.w.word&LWMASK) == ikey; ++iie );
        for( wwe = wb; wwe != we && (wwe->wd.w.word&LWMASK) == wkey; ++wwe );

        for( ic = ib; ic != iie; ++ic )
        {
            uint32_t iword( (ic->wd.w.word&HWMASK)>>WSHIFT );

            for( wc = wb; wc != wwe; ++wc )
            {
                uint32_t wword( (wc->wd.w.word&HWMASK)>>WSHIFT );
                size_t n_init_matches( 0 );

                if( CheckDel( iword, wword, n_init_matches ) )
                {
                    if( n_init_matches >= 4 )
                    {
                        SaveHit( *ic, *wc, false, 0, 1 );
                    }
                    else
                    {
                        SaveHit( *ic, *wc, false, 0, 0, WORD_BASES );
                    }
                }
            }
        }

        ib = iie;
        wb = wwe;
    }
}

//------------------------------------------------------------------------------
inline void CFastSeeds::SeedSearchJob::MatchPfx_2( uint32_t anchor )
{
    static uint32_t const LWMASK = 0xFFFUL;

    auto ib( o_.fsidx_.cbegin( anchor ) ),
         ie( o_.fsidx_.cend( anchor ) );
    auto & idxwords( jd_.idxwords );
    idxwords.clear();

    for( ; ib != ie; ++ib )
    {
        if( o_.seeder_mode_ || !ib->wd.w.repeat )
        {
            idxwords.push_back( *ib );
        }
    }

    std::sort( idxwords.begin(), idxwords.end(),
               []( IndexEntry const & x, IndexEntry const & y )
               { return (x.wd.w.word&LWMASK) < (y.wd.w.word&LWMASK); } );
    MatchPfx_2I( anchor );
    MatchPfx_2E( anchor );
    MatchPfx_2D( anchor );
}

//------------------------------------------------------------------------------
inline void CFastSeeds::SeedSearchJob::MatchPfx_1( uint32_t anchor )
{
    static uint32_t const LWMASK = 0xFFFUL;
    static uint32_t const WMASK = 0xFFF000UL;
    static size_t const WA_SHIFT = WORD_BITS - LB;

    auto & wt( o_.wt_ );
    auto const & wmap( o_.wmap_ );
    auto ib( o_.fsidx_.cbegin( anchor ) ),
         ie( o_.fsidx_.cend( anchor ) ),
         iie( ib );

    for( uint32_t l( 0 ); l < 3; ++l )
    {
        uint32_t ll( l<<WA_SHIFT );
        auto wa( (anchor<<LB) + l );
        WordTable::iterator wtb( wt.begin() + wmap[wa] ),
                            wte( wt.begin() + wmap[wa + 1] );
        for( iie = ib;
             iie != ie && (size_t)(iie->wd.w.word>>WA_SHIFT) < l;
             ++iie );
        ib = iie;
        for( ; iie != ie && (size_t)(iie->wd.w.word>>WA_SHIFT) == l;
               ++iie );

        while( ib != iie && wtb != wte )
        {
            uint32_t iw( ib->wd.w.word&WMASK ),
                     ww( ((wtb->wd.w.word>>LB) + ll)&WMASK );

            if( iw < ww )
            {
                ++ib;
            }
            else if( iw > ww )
            {
                ++wtb;
            }
            else
            {
                auto iiie( ib );
                for( ; iiie != iie && (iiie->wd.w.word&WMASK) == iw; ++iiie );
                auto wtte( wtb );
                for( ; wtte != wte
                       && (((wtte->wd.w.word>>LB) + ll)&WMASK) == ww; ++wtte );

                for( ; ib != iiie; ++ib )
                {
                    if( !o_.seeder_mode_ && ib->wd.w.repeat )
                    {
                        continue;
                    }

                    uint32_t iww( ib->wd.w.word );
                    uint32_t liw( iww&LWMASK );

                    for( auto wtc( wtb ); wtc != wtte; ++wtc )
                    {
                        if( iww != wtc->wd.w.word )
                        {
                            uint32_t lww( (wtc->wd.w.word>>LB)&LMASK );

                            if( Check1Err( liw, lww ) < 3 )
                            {
                                SaveHit( *ib, *wtc, false, 0, 0,
                                         HALF_WORD_BASES );
                            }
                        }
                    }
                }

                wtb = wtte;
            }
        }
    }
}

//------------------------------------------------------------------------------
inline void CFastSeeds::SeedSearchJob::operator()()
{
    auto & ctx( o_.bctx_.GetSearchCtx() );
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
                MatchIns( task_idx );
                MatchSubst( task_idx );
                MatchDel( task_idx );
            }
            else
            {
                PrepareIndex( task_idx, false );
                MatchScan( task_idx );
            }

            if( !ctx.exact_seeds )
            {
                MatchPfx_1( task_idx );
                MatchPfx_2( task_idx );
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
    auto d1( -diag_delta - 1 ), dc( d1 );
    cth = cover_.size()*cth;

    for( auto c( b ); b != e; ++b )
    {
        auto d2( b->second.GetDiag() );

        if( d2 > d1 )
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

        for( auto pos( h.mid ? h.readpos - WORD_BASES : h.readpos ), i( 0L );
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
    fsidx_.Unload();
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
                std::vector< size_t > job_sizes( n_jobs, 0 );

                for( auto const & job_data : jd )
                {
                    for( auto const & h : job_data.results )
                    {
                        size_t job_num( h.read/reads_per_job );
                        ++job_sizes[job_num];
                    }
                }

                for( size_t i( 0 ); i < n_jobs; ++i )
                {
                    filter_job_hits[i].reserve( job_sizes[i] );
                }
            }

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

            // if( seeder_mode_ )
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
inline void CFastSeeds::UpdateAnchorLists( uint32_t wanchor )
{
    static uint32_t const LMASK = (1U<<LB) - 1;

    AnchorListEntry ale;
    ale.f.wanchor = wanchor;

    // exact
    //
    {
        ale.f.ins = ale.f.del = 0;
        ale.f.exact = true;
        uint32_t ww( wanchor>>LB );
        assert( atbl_[ww].len < ANCHOR_LIST_MAX_LEN );
        atbl_[ww].data[atbl_[ww].len++] = ale;
        ++atbl_[ww].n_subst;
    }

    // substitutions
    //
    {
        ale.f.exact = false;
        size_t shift( LB );
        uint32_t mask( ((1<<LB) - 1)<<LB );

        for( size_t i( 0 ); i < ANCHOR_BASES; ++i, shift += LB, mask <<= LB )
        {
            uint32_t l( (wanchor&mask)>>shift ),
                     w( wanchor&~mask );

            for( uint32_t j( 1 ); j < 4; ++j )
            {
                uint32_t ww( (w + (((l + j)%4)<<shift))>>LB );
                assert( ww < ANCHOR_TBL_SIZE );
                assert( atbl_[ww].len < ANCHOR_LIST_MAX_LEN );
                atbl_[ww].data[atbl_[ww].len++] = ale;
                ++atbl_[ww].n_subst;
            }
        }
    }

    // deletions
    //
    {
        ale.f.del = 1;
        size_t shift( LB + LB );
        uint32_t pl( 4 ),
                 mask( (LMASK<<LB) + LMASK );

        for( size_t i( 0 );
             i < ANCHOR_BASES; ++i, shift += LB, mask = (mask<<LB) + LMASK,
                               pl = ((wanchor>>shift)&LMASK) )
        {
            for( uint32_t l( 0 ); l < 4; ++l )
            {
                if( l != pl )
                {
                    uint32_t ww( ((wanchor&mask)>>LB) + (wanchor&~mask) );
                    ww += (l<<(shift - LB));
                    ww >>= LB;
                    assert( ww < ANCHOR_TBL_SIZE );
                    assert( atbl_[ww].len < ANCHOR_LIST_MAX_LEN );
                    atbl_[ww].data[atbl_[ww].len++] = ale;
                    ++atbl_[ww].n_del;
                }
            }
        }
    }

    // insertions
    //
    {
        ale.f.del = 0;
        ale.f.ins = 1;
        size_t shift( LB );
        uint32_t pl( 4 ),
                 mask( LMASK ),
                 l;

        for( size_t i( 0 ); i < ANCHOR_BASES; ++i, pl = l, shift += LB )
        {
            l = ((wanchor>>shift)&LMASK);

            if( l != pl )
            {
                uint32_t ww( (wanchor&mask)<<LB );
                mask = (mask<<LB) + LMASK;
                ww += (wanchor&~mask);
                ww >>= LB;
                assert( ww < ANCHOR_TBL_SIZE );
                assert( atbl_[ww].len < ANCHOR_LIST_MAX_LEN );
                atbl_[ww].data[atbl_[ww].len++] = ale;
                ++atbl_[ww].n_ins;
            }
        }
    }
}

//------------------------------------------------------------------------------
inline void CFastSeeds::UpdateAnchorUseMap( uint32_t wa )
{
    static uint32_t const LMASK = (1U<<LB) - 1;
    auto & ctx( bctx_.GetSearchCtx() );

    // exact
    //
    {
        uint32_t a( wa>>LB );
        anchor_use_map_[a] = true;
    }

    if( !ctx.exact_seeds )
    {
        // substitutions
        //
        {
            size_t shift( LB );
            uint32_t mask( ((1<<LB) - 1)<<LB );

            for( size_t i( 0 );
                 i < ANCHOR_BASES; ++i, shift += LB, mask <<= LB )
            {
                uint32_t l( (wa&mask)>>shift ),
                         w( wa&~mask );

                for( uint32_t j( 1 ); j < 4; ++j )
                {
                    uint32_t a( (w + (((l + j)%4)<<shift))>>LB );
                    anchor_use_map_[a] = true;
                }
            }
        }

        // deletions
        //
        {
            size_t shift( LB + LB );
            uint32_t pl( 4 ),
                     mask( (LMASK<<LB) + LMASK );

            for( size_t i( 0 );
                 i < ANCHOR_BASES; ++i, shift += LB, mask = (mask<<LB) + LMASK,
                                   pl = ((wa>>shift)&LMASK) )
            {
                for( uint32_t l( 0 ); l < 4; ++l )
                {
                    if( l != pl )
                    {
                        uint32_t a( ((wa&mask)>>LB) + (wa&~mask) );
                        a += (l<<(shift - LB));
                        a >>= LB;
                        anchor_use_map_[a] = true;
                    }
                }
            }
        }

        // insertions
        //
        {
            size_t shift( LB );
            uint32_t pl( 4 ),
                     mask( LMASK ),
                     l;

            for( size_t i( 0 ); i < ANCHOR_BASES; ++i, pl = l, shift += LB )
            {
                l = ((wa>>shift)&LMASK);

                if( l != pl )
                {
                    uint32_t a( (wa&mask)<<LB );
                    mask = (mask<<LB) + LMASK;
                    a += (wa&~mask);
                    a >>= LB;
                    anchor_use_map_[a] = true;
                }
            }
        }
    }
}

//------------------------------------------------------------------------------
inline void CFastSeeds::CreateAnchorTable()
{
    auto & ctx( bctx_.GetSearchCtx() );
    StopWatch w( ctx.logger_ );

    for( uint32_t wa( 0 ); wa < WMAP_SIZE - 1; ++wa )
    {
        UpdateAnchorLists( wa );
    }

    for( uint32_t a( 0 ); a < ANCHOR_TBL_SIZE; ++a )
    {
        auto & d( atbl_[a].data );
        size_t first_subst( 0 ),
               first_del( atbl_[a].len ),
               curr( 0 );

        while( curr < first_del )
        {
            auto & e( d[curr] );

            if( e.f.del == 1 )
            {
                std::swap( e, d[--first_del] );
            }
            else if( e.f.ins == 1 )
            {
                std::swap( e, d[first_subst++] );
                ++curr;
            }
            else
            {
                ++curr;
            }
        }
    }

    M_INFO( ctx.logger_, "generated anchor table" );
}

//------------------------------------------------------------------------------
void CFastSeeds::Run()
{
    CreateAnchorTable();
    CreateWordTable();

    auto & ctx( bctx_.GetSearchCtx() );

    if( !ctx.db_name.empty() )
    {
        try { fsidx_.Load( ctx.db_name ); }
        catch( std::exception const & e )
        {
            M_INFO( ctx.logger_, "index load failed: " );
            M_INFO( ctx.logger_, e.what() );
            fsidx_.Create( ctx.n_threads );
        }
    }
    else
    {
        fsidx_.Create( ctx.n_threads );
    }

    ComputeSeeds();
}

READFINDER_NS_END

