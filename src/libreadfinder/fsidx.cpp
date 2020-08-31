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
#include <vector>

#include <libreadfinder/fsidx.hpp>

#include <libtools/progress.hpp>
#include <libtools/stopwatch.hpp>
#include <libtools/taskarray.hpp>

READFINDER_NS_BEGIN

//------------------------------------------------------------------------------
namespace {
    static char const * IDXMAP_EXT = ".idm";
    static char const * IDX_EXT = ".idx";
    static char const * LFM_EXT = ".lfm";
}

//==============================================================================
/// Descriptor for an index population job; each job is handled by one thread.
///
struct CFastSeedsIndex::PopulateIndexJobData
{
    TRefOId start_id,
            end_id;
    TSeqOff start_off,
            end_off;
    WordMap idxmap;

    friend std::ostream & operator<<(
            std::ostream & os, PopulateIndexJobData const & x )
    {
        return os << x.start_id << ':' << x.start_off << "---"
                  << x.end_id << ':' << x.end_off;
    }
};

//==============================================================================
class CFastSeedsIndex::IndexEntrySource
{
public:

    IndexEntrySource() {}
    IndexEntrySource( TWord const * seq_data, TSeqLen len, TSeqOff start,
                      uint32_t init_pos );

    void operator++();

    uint32_t GetAnchor() const { return data_.f.anchor; }
    uint32_t GetWord() const { return data_.f.word; }
    // uint32_t GetSfx() const { return data_.f.sfx; }
    uint32_t GetPos() const { return pos_; }
    operator bool() const { return !done_; }

private:

    struct
    {
        union
        {
            struct
            {
                uint64_t word   : WORD_BITS;
                uint64_t anchor : ANCHOR_BITS;
                // uint64_t sfx    : SFX_BITS;
            } f;

            TWord w;
        };
    } data_;

    TWord nw_;

    static_assert( sizeof( data_ ) == 8, "" );

    TWord const * seq_data_;
    TSeqLen len_;
    TSeqOff off_;
    uint32_t pos_;
    bool done_ = true;
};

//------------------------------------------------------------------------------
inline CFastSeedsIndex::IndexEntrySource::IndexEntrySource(
        TWord const * seq_data, TSeqLen len, TSeqOff start, uint32_t init_pos )
    : nw_( 0 ), seq_data_( seq_data ), len_( start + len ), off_( start ),
      pos_( init_pos + start )
{
    if( off_ >= len_ )
    {
        return;
    }

    done_ = false;
    seq_data_ += off_/WL;
    TSeqOff woff( off_%WL );
    data_.w = *seq_data_++;
    nw_ = *seq_data_++;

    if( woff > 0 )
    {
        auto shift( woff*LB );
        data_.w >>= shift;
        data_.w += (nw_<<(LB*WL - shift));
        nw_ >>= shift;
    }
}

//------------------------------------------------------------------------------
inline void CFastSeedsIndex::IndexEntrySource::operator++()
{
    if( ++off_ >= len_ )
    {
        done_ = true;
        return;
    }

    data_.w >>= LB;
    data_.w += (nw_<<((WL - 1)*LB));
    ++pos_;

    if( off_%WL == 0 )
    {
        nw_ = *seq_data_++;
    }
    else
    {
        nw_ >>= LB;
    }
}

//==============================================================================
class CFastSeedsIndex::IndexMaskSource
{
public:

    IndexMaskSource() {}
    IndexMaskSource( TWord const * mask_data, TSeqLen len, TSeqOff start );

    void operator++();

    uint64_t GetNMer() const { return data_.f.nmer; }
    // uint64_t GetSfx() const { return data_.f.sfx; }

private:

    struct
    {
        union
        {
            struct
            {
                uint64_t nmer : NMER_BITS;
                // uint64_t sfx  : SFX_BITS;
            } f;

            TWord w;
        };
    } data_;

    TWord nw_;

    static_assert( sizeof( data_ ) == 8, "" );

    TWord const * mask_data_;
    TSeqLen len_;
    TSeqOff off_;
    bool done_ = true;
};

//------------------------------------------------------------------------------
inline CFastSeedsIndex::IndexMaskSource::IndexMaskSource(
        TWord const * mask_data, TSeqLen len, TSeqOff start )
    : nw_( 0 ), mask_data_( mask_data ), len_( start + len ), off_( start )
{
    if( off_ >= len_ )
    {
        return;
    }

    done_ = false;
    mask_data_ += off_/WL;
    TSeqOff woff( off_%WL );
    data_.w = *mask_data_++;
    nw_ = *mask_data_++;

    if( woff > 0 )
    {
        auto shift( woff*LB );
        data_.w >>= shift;
        data_.w += (nw_<<(LB*WL - shift));
        nw_ >>= shift;
    }
}

//------------------------------------------------------------------------------
inline void CFastSeedsIndex::IndexMaskSource::operator++()
{
    if( ++off_ >= len_ )
    {
        done_ = true;
        return;
    }

    data_.w >>= LB;
    data_.w += (nw_<<((WL - 1)*LB));

    if( off_%WL == 0 )
    {
        nw_ = *mask_data_++;
    }
    else
    {
        nw_ >>= LB;
    }
}

//==============================================================================
struct CFastSeedsIndex::IndexNMerCountingJob
{
    typedef std::vector< PopulateIndexJobData > PIJData;

    IndexNMerCountingJob(
            CFastSeedsIndex & o, PIJData & pij_data, size_t & job_idx,
            CProgress::ProgressHandle const & ph )
        : o_( o ), pij_data_( pij_data[job_idx] ), ph_( ph[job_idx++] )
    {}

    void operator()();

    CFastSeedsIndex & o_;
    PopulateIndexJobData & pij_data_;
    CProgress::ProgressHandle ph_;
};

//------------------------------------------------------------------------------
inline void CFastSeedsIndex::IndexNMerCountingJob::operator()()
{
    pij_data_.idxmap.resize( IDXMAP_SIZE, 0 );
    typedef std::tuple< TRefOId, TSeqOff, TSeqLen > SubJob;
    auto const & refs( o_.refs_ );
    auto const & roff( o_.roff_ );
    std::vector< SubJob > sub_jobs;

    if( pij_data_.start_id == pij_data_.end_id )
    {
        sub_jobs.push_back(
                SubJob { pij_data_.start_id, pij_data_.start_off,
                         pij_data_.end_off - pij_data_.start_off } );
    }
    else
    {
        sub_jobs.push_back(
                SubJob { pij_data_.start_id, pij_data_.start_off,
                         refs.GetLength( pij_data_.start_id )
                         - pij_data_.start_off } );

        for( auto i( pij_data_.start_id + 1 ); i < pij_data_.end_id; ++i )
        {
            sub_jobs.push_back( SubJob { i, 0, refs.GetLength( i ) } );
        }

        if( pij_data_.end_off > 0 )
        {
            sub_jobs.push_back(
                    SubJob { pij_data_.end_id, 0, pij_data_.end_off } );
        }
    }

    ph_.SetTotal( sub_jobs.size() );

    for( auto const & sj : sub_jobs )
    {
        IndexEntrySource ies( refs.GetSeqData( std::get< 0 >( sj ) ),
                              std::get< 2 >( sj ), std::get< 1 >( sj ),
                              roff[std::get< 0 >( sj )] );
        IndexMaskSource ims( refs.GetMaskData( std::get< 0 >( sj ) ),
                             std::get< 2 >( sj ), std::get< 1 >( sj ) );

        while( ies )
        {
            if( GetNSet( ims.GetNMer() ) <= LB )
            {
                ++pij_data_.idxmap[ies.GetAnchor()];
            }

            ++ies;
            ++ims;
        }

        ph_.Increment();
    }
}

//==============================================================================
struct CFastSeedsIndex::PopulateIndexJob
{
    typedef std::vector< PopulateIndexJobData > PIJData;

    PopulateIndexJob(
            CFastSeedsIndex & o, PIJData & pij_data, size_t & job_idx,
            CProgress::ProgressHandle const & ph )
        : o_( o ), pij_data_( pij_data[job_idx] ), ph_( ph[job_idx++] )
    {}

    void operator()();

    CFastSeedsIndex & o_;
    PopulateIndexJobData & pij_data_;
    CProgress::ProgressHandle ph_;
};

//------------------------------------------------------------------------------
inline void CFastSeedsIndex::PopulateIndexJob::operator()()
{
    typedef std::tuple< TRefOId, TSeqOff, TSeqLen > SubJob;
    auto const & refs( o_.refs_ );
    auto const & roff( o_.roff_ );
    std::vector< SubJob > sub_jobs;

    if( pij_data_.start_id == pij_data_.end_id )
    {
        sub_jobs.push_back(
                SubJob { pij_data_.start_id, pij_data_.start_off,
                         pij_data_.end_off - pij_data_.start_off } );
    }
    else
    {
        sub_jobs.push_back(
                SubJob { pij_data_.start_id, pij_data_.start_off,
                         refs.GetLength( pij_data_.start_id )
                         - pij_data_.start_off } );

        for( auto i( pij_data_.start_id + 1 ); i < pij_data_.end_id; ++i )
        {
            sub_jobs.push_back( SubJob { i, 0, refs.GetLength( i ) } );
        }

        if( pij_data_.end_off > 0 )
        {
            sub_jobs.push_back(
                    SubJob { pij_data_.end_id, 0, pij_data_.end_off } );
        }
    }

    ph_.SetTotal( sub_jobs.size() );
    auto & idx( o_.idx_ );
    auto & chunk_map( o_.chunk_map_ );
    auto idxmap( pij_data_.idxmap );

    for( auto const & sj : sub_jobs )
    {
        IndexEntrySource ies( refs.GetSeqData( std::get< 0 >( sj ) ),
                              std::get< 2 >( sj ), std::get< 1 >( sj ),
                              roff[std::get< 0 >( sj )] );
        IndexMaskSource ims( refs.GetMaskData( std::get< 0 >( sj ) ),
                             std::get< 2 >( sj ), std::get< 1 >( sj ) );

        while( ies )
        {
            if( GetNSet( ims.GetNMer() ) <= LB )
            {
                auto anchor( ies.GetAnchor() );
                IndexEntry ie;
                ie.wd.data = 0;
                ie.wd.w.word = ies.GetWord();
                // ie.wd.w.sfx = ies.GetSfx();
                ie.pos = ies.GetPos();
                idx[chunk_map[anchor]].SetData( idxmap[anchor]++, ie );
            }

            ++ies;
            ++ims;
        }

        ph_.Increment();
    }
}

//==============================================================================
struct CFastSeedsIndex::SortIndexJob
{
    typedef std::atomic< uint32_t > JobIdx;

    SortIndexJob( CFastSeedsIndex & o, JobIdx & job_idx,
                  CProgress::ProgressHandle & ph )
        : o_( o ), job_idx_( job_idx ), ph_( ph )
    {}

    void operator()();

    CFastSeedsIndex & o_;
    JobIdx & job_idx_;
    CProgress::ProgressHandle & ph_;
};

//------------------------------------------------------------------------------
inline void CFastSeedsIndex::SortIndexJob::operator()()
{
    while( true )
    {
        uint32_t anchor( job_idx_.fetch_add( 1 ) );

        if( anchor >= IDXMAP_SIZE - 1 )
        {
            break;
        }

        auto b( o_.begin( anchor ) ),
             e( o_.end( anchor ) );
        std::sort( b, e,
                   []( IndexEntry const & x, IndexEntry const & y )
                   { return x.wd.w.word < y.wd.w.word; } );
        ph_.Increment();
    }
}

//==============================================================================
struct CFastSeedsIndex::RepeatsFinderJob
{
    struct Data
    {
        size_t n_repeats = 0,
               n_rpos = 0,
               n_erepeats = 0,
               n_erpos = 0;
    };

    RepeatsFinderJob(
        CFastSeedsIndex & o,
        std::atomic< uint32_t > & task_idx,
        Data * data, size_t & job_idx,
        CProgress::ProgressHandle & ph )
    : o_( o ), task_idx_( task_idx ), data_( data[job_idx++] ), ph_( ph )
    {}

    void operator()();

    CFastSeedsIndex & o_;
    std::atomic< uint32_t > & task_idx_;
    Data & data_;
    CProgress::ProgressHandle & ph_;
};

//------------------------------------------------------------------------------
inline void CFastSeedsIndex::RepeatsFinderJob::operator()()
{
    static size_t const EXACT_REPEAT_THRESHOLD = 128;
    static size_t const REPEAT_THRESHOLD = 16;

    while( true )
    {
        uint32_t anchor( task_idx_.fetch_add( 1 ) );
        if( anchor >= IDXMAP_SIZE - 1 ) break;
        auto * ib( o_.begin( anchor ) ),
             * ie( o_.end( anchor ) ),
             * iie( ib );

        while( ib != ie )
        {
            uint32_t w( ib->wd.w.word );
            for( iie = ib; iie != ie && iie->wd.w.word == w; ++iie );

            if( (size_t)(iie - ib) > EXACT_REPEAT_THRESHOLD )
            {
                ++data_.n_erepeats;
                ++data_.n_repeats;
                data_.n_erpos += (iie - ib);
                data_.n_rpos += (iie - ib);

                for( ; ib != iie; ++ib )
                {
                    ib->wd.w.erepeat = true;
                    ib->wd.w.repeat = true;
                }
            }
            else if( (size_t)(iie - ib) > REPEAT_THRESHOLD )
            {
                ++data_.n_repeats;
                data_.n_rpos += (iie - ib);

                for( ; ib != iie; ++ib )
                {
                    ib->wd.w.repeat = true;
                    ib->wd.w.erepeat = false;
                }
            }
            else
            {
                for( ; ib != iie; ++ib )
                {
                    ib->wd.w.repeat = false;
                    ib->wd.w.erepeat = false;
                }
            }
        }

        ph_.Increment();
    }
}

//==============================================================================
struct CFastSeedsIndex::CreateFreqHistJob
{
    struct Data
    {
        std::vector< size_t > fhist = std::vector< size_t >( 65, 0 );
    };

    CreateFreqHistJob(
        CFastSeedsIndex & o,
        std::atomic< uint32_t > & task_idx,
        Data * data, size_t & job_idx,
        CProgress::ProgressHandle & ph )
    : o_( o ), task_idx_( task_idx ), data_( data[job_idx++] ), ph_( ph )
    {}

    void operator()();

    CFastSeedsIndex & o_;
    std::atomic< uint32_t > & task_idx_;
    Data & data_;
    CProgress::ProgressHandle & ph_;
};

//------------------------------------------------------------------------------
inline void CFastSeedsIndex::CreateFreqHistJob::operator()()
{
    while( true )
    {
        uint32_t anchor( task_idx_.fetch_add( 1 ) );
        if( anchor >= IDXMAP_SIZE - 1 ) break;
        auto * ib( o_.begin( anchor ) ),
             * ie( o_.end( anchor ) ),
             * iie( ib );

        for( ; ib != ie; ib = iie )
        {
            uint32_t w( ib->wd.w.word );
            for( iie = ib; iie != ie && iie->wd.w.word == w; ++iie );
            uint64_t f( iie - ib );
            assert( f > 0ULL );
            uint64_t i( 0ULL );

            for( uint64_t j( 1ULL ); i < 64ULL; ++i, j <<= 1 )
            {
                if( f <= j )
                {
                    ++data_.fhist[i];
                    break;
                }
            }

            if( i == 64 ) ++data_.fhist[64];
        }

        ph_.Increment();
    }
}

//==============================================================================
struct CFastSeedsIndex::CreateFreqTableJob
{
    struct Data
    {
        std::vector< FreqTableEntry > freq_table_;
    };

    CreateFreqTableJob(
        CFastSeedsIndex & o,
        std::atomic< uint32_t > & task_idx,
        Data * data, size_t & job_idx,
        CProgress::ProgressHandle & ph )
    : o_( o ), task_idx_( task_idx ), data_( data[job_idx++] ), ph_( ph )
    {}

    void operator()();

    CFastSeedsIndex & o_;
    std::atomic< uint32_t > & task_idx_;
    Data & data_;
    CProgress::ProgressHandle & ph_;
};

//------------------------------------------------------------------------------
inline void CFastSeedsIndex::CreateFreqTableJob::operator()()
{
    while( true )
    {
        uint32_t anchor( task_idx_.fetch_add( 1 ) );
        if( anchor >= IDXMAP_SIZE - 1 ) break;
        FreqTableEntry fte;
        fte.data.f.anchor = anchor;
        auto * ib( o_.begin( anchor ) ),
             * ie( o_.end( anchor ) ),
             * iie( ib );

        for( ; ib != ie; ib = iie )
        {
            uint32_t w( ib->wd.w.word );
            for( iie = ib; iie != ie && iie->wd.w.word == w; ++iie );
            uint64_t f( iie - ib );
            assert( f > 0ULL );
            uint64_t i( 0ULL );

            for( uint64_t j( 1ULL ); i < 64ULL; ++i, j <<= 1 )
            {
                if( f <= j ) break;
            }

            if( i > o_.cutoff_idx_ )
            {
                fte.data.f.word = w;
                fte.data.f.freq = i;
                data_.freq_table_.push_back( fte );
            }
        }

        ph_.Increment();
    }
}

//==============================================================================
//------------------------------------------------------------------------------
CFastSeedsIndex::IndexChunk::IndexChunk(
        CCommonContext & ctx, IndexMap const & idxmap,
        size_t start_anchor, size_t end_anchor, size_t sz )
    : start_anchor_( start_anchor ), end_anchor_( end_anchor ),
      idxmap_( &idxmap ),
      data_off_( idxmap[start_anchor] ),
      data_( idxmap[end_anchor] - idxmap[start_anchor], IndexEntry() ),
      loaded_( true )
{
    assert( sz == data_.size()*sizeof( IndexEntry ) );
}

//------------------------------------------------------------------------------
CFastSeedsIndex::IndexChunk::IndexChunk(
        CCommonContext & ctx, IndexMap const & idxmap,
        size_t start_anchor, size_t end_anchor )
    : start_anchor_( start_anchor ), end_anchor_( end_anchor ),
      idxmap_( &idxmap ),
      data_off_( idxmap[start_anchor] )
{
}

//------------------------------------------------------------------------------
void CFastSeedsIndex::IndexChunk::Save( std::ostream & os ) const
{
    assert( loaded_ );
    size_t sz( data_.size()*sizeof( IndexEntry ) );
    os.write( (char const *)&data_[0], sz );

    if( !os )
    {
        M_THROW( "error writing " << sz << " bytes of index data" );
    }
}

//------------------------------------------------------------------------------
void CFastSeedsIndex::IndexChunk::Load( std::ifstream & is, size_t * used_mem )
{
    if( !loaded_ )
    {
        is.seekg( data_off_*sizeof( IndexEntry ), std::ios_base::beg );
        size_t sz( (*idxmap_)[end_anchor_] - data_off_ );
        data_.resize( sz );
        sz *= sizeof( IndexEntry );
        is.read( (char *)&data_[0], sz );

        if( !is )
        {
            M_THROW( "error reading " << sz << " bytes of index data" );
        }

        if( used_mem != nullptr ) *used_mem += sz;
        loaded_ = true;
    }
}

//------------------------------------------------------------------------------
void CFastSeedsIndex::IndexChunk::Unload( CCommonContext & ctx )
{
    Data().swap( data_ );
    loaded_ = false;
}

//==============================================================================
//------------------------------------------------------------------------------
CFastSeedsIndex::CFastSeedsIndex( CCommonContext & ctx, CRefData const & refs )
    : ctx_( ctx ), refs_( refs )
{
}

//------------------------------------------------------------------------------
void CFastSeedsIndex::SetUpChunks( bool allocate )
{
    auto & logger( ctx_.logger_ );

    //
    // chunk data contains (in order) the first anchor, one past
    // the last anchor, and size (in bytes) of the chunk
    //
    typedef std::tuple< size_t, size_t, size_t > ChunkData;
    std::vector< ChunkData > chunk_data;
    chunk_map_.resize( IDXMAP_SIZE, 0 );
    used_mem_ += IDXMAP_SIZE*sizeof( uint32_t );
    size_t j( 0 ),
           s( 0 );

    for( size_t i( 1 ); i < IDXMAP_SIZE; ++i )
    {
        s += (idxmap_[i] - idxmap_[i-1])*sizeof( IndexEntry );

        if( s >= MIN_IDX_CHUNK_SIZE )
        {
            chunk_data.push_back( ChunkData{ j, i, s } );
            s = 0;
            j = i;
        }
    }

    if( s > 0 )
    {
        chunk_data.push_back( ChunkData( j, IDXMAP_SIZE - 1, s ) );
    }

    M_INFO( logger, "number of index chunks: " << chunk_data.size() );

    for( size_t i( 0 ); i < chunk_data.size(); ++i )
    {
        auto const & cd( chunk_data[i] );

        for( auto anchor_start( std::get< 0 >( cd ) ),
                  anchor_end( std::get< 1 >( cd ) );
             anchor_start < anchor_end; ++anchor_start )
        {
            chunk_map_[anchor_start] = i;
        }
    }

    M_INFO( logger, "generated chunk map" );
    idx_.clear();

    for( auto const & cd : chunk_data )
    {
        if( allocate )
        {
            idx_.emplace_back( ctx_, idxmap_,
                               std::get< 0 >( cd ),
                               std::get< 1 >( cd ),
                               std::get< 2 >( cd ) );
            used_mem_ += std::get< 2 >( cd );
        }
        else
        {
            idx_.emplace_back( ctx_, idxmap_,
                               std::get< 0 >( cd ),
                               std::get< 1 >( cd ) );
        }

        used_mem_ += sizeof( IndexChunk );
    }

    M_INFO( logger, "initialized index chunk data" );
}

//------------------------------------------------------------------------------
CFastSeedsIndex & CFastSeedsIndex::Create( size_t n_threads )
{
    auto & logger( ctx_.logger_ );

    if( keep_loaded_ )
    {
        M_INFO( logger, "reference index is already created" )
        return *this;
    }

    M_INFO( logger, "generating reference index" );
    ssize_t n_job_pos( 0 );
    std::vector< PopulateIndexJobData > pij_data;
    pij_data.reserve( n_threads );

    // compute number of reference positions per nmer counting job and
    // generate job descriptors
    //
    {
        StopWatch w( logger );
        roff_.resize( refs_.GetSize() + 1 );
        used_mem_ += roff_.size()*sizeof( uint32_t );

        for( size_t i( 0 ), ie( refs_.GetSize() ); i < ie; ++i )
        {
            roff_[i] = n_job_pos;
            n_job_pos += refs_.GetLength( i );
        }

        roff_[refs_.GetSize()] = n_job_pos;
        M_INFO( logger, "genome length is " << n_job_pos << " bases" );
        n_job_pos = 1 + n_job_pos/n_threads;
        M_INFO( logger, n_job_pos << " bases per job" );
        ssize_t p( 0 );

        for( TRefOId i( 0 ), ie( refs_.GetSize() ); i < ie; )
        {
            TSeqLen l( refs_.GetLength( i ) );

            if( i > 0 || p > 0 )
            {
                pij_data.back().end_id = i;
                pij_data.back().end_off = p;
            }

            while( p + n_job_pos < l )
            {
                pij_data.push_back( PopulateIndexJobData {
                        i, i, (TSeqOff)p, (TSeqOff)(p + n_job_pos),
                        WordMap() } );
                p += n_job_pos;
            }

            pij_data.push_back( PopulateIndexJobData {
                    i, i, (TSeqOff)p, (TSeqOff)p, WordMap() } );
            p += n_job_pos;

            while( p >= l && ++i < ie )
            {
                p -= l;
                l = refs_.GetLength( i );
            }
        }

        pij_data.back().end_id = refs_.GetSize();
        pij_data.back().end_off = 0;
        assert( pij_data.size() <= n_threads );

        if( pij_data.size() < n_threads )
        {
            n_threads = pij_data.size();
        }

        if( n_threads ==0 )
        {
            n_threads = 1;
        }
        M_INFO( logger,
                "created " << pij_data.size() << " position counting jobs" );
    }

    // count nmers for each index creating job
    //
    if( pij_data.size() > 0 )
    {
        StopWatch w( logger );
        size_t job_idx( 0 );
        CProgress p( "counting reference words", "tasks",
                     ctx_.progress_flags_ );
        auto ph( p.GetTop().Split( n_threads ) );
        CTaskArray< IndexNMerCountingJob > jobs(
                n_threads, *this, pij_data, job_idx, ph );
        p.Start();
        jobs.Start();
        p.Stop();
        M_INFO( logger, "nmer counting complete" );
    }

    // count number of positions per anchor value in total and per thread
    //
    {
        StopWatch w( logger );
        idxmap_.resize( IDXMAP_SIZE, 0 );
        used_mem_ += idxmap_.size()*sizeof( uint32_t );
        uint32_t idxmap_off( 0 );

        for( size_t i( 0 ); i < IDXMAP_SIZE; ++i )
        {
            idxmap_[i] = idxmap_off;

            for( size_t j( 0 ); j < pij_data.size(); ++j )
            {
                auto t( pij_data[j].idxmap[i] );
                pij_data[j].idxmap[i] = idxmap_off;
                idxmap_off += t;
            }
        }

        M_INFO( logger, "nmer counts aggregation complete" );
        M_INFO( logger, idxmap_off << " entries in the index" );
    }

    SetUpChunks();

    // fill in index data
    //
    if( pij_data.size() > 0 )
    {
        StopWatch w( logger );
        size_t job_idx( 0 );
        CProgress p( "populating reference index", "tasks",
                     ctx_.progress_flags_ );
        auto ph( p.GetTop().Split( n_threads ) );
        CTaskArray< PopulateIndexJob > jobs(
                n_threads, *this, pij_data, job_idx, ph );
        p.Start();
        jobs.Start();
        p.Stop();
        M_INFO( logger, "reference index data generation complete" );
    }

    // sort index data
    //
    if( pij_data.size() > 0 )
    {
        StopWatch w( logger );
        std::atomic< uint32_t > job_idx( 0 );
        CProgress p( "sorting reference index", "anchors",
                     ctx_.progress_flags_ );
        auto ph( p.GetTop() );
        ph.SetTotal( IDXMAP_SIZE - 1 );
        CTaskArray< SortIndexJob > jobs( n_threads, *this, job_idx, ph );
        p.Start();
        jobs.Start();
        p.Stop();
        M_INFO( logger, "sorting of index data complete" );
    }

    // identify repeats
    //
    if( pij_data.size() > 0 )
    {
        StopWatch w( logger );
        size_t n_repeats( 0 ),
               n_rpos( 0 ),
               n_erepeats( 0 ),
               n_erpos( 0 );
        CProgress p( "identifying repeats", "anchors",
                     ctx_.progress_flags_ );
        auto ph( p.GetTop() );
        ph.SetTotal( IDXMAP_SIZE - 1 );
        std::vector< RepeatsFinderJob::Data > jd( n_threads );
        std::atomic< uint32_t > task_idx( 0 );
        size_t job_idx( 0 );
        CTaskArray< RepeatsFinderJob > jobs(
            n_threads, *this, task_idx, &jd[0], job_idx, ph );
        p.Start();
        jobs.Start();
        p.Stop();
        
        for( auto const & data : jd )
        {
            n_erepeats += data.n_erepeats;
            n_erpos += data.n_erpos;
            n_repeats += data.n_repeats;
            n_rpos += data.n_rpos;
        }

        M_INFO( logger, n_erepeats << " exact repeats marked" );
        M_INFO( logger, n_erpos << " exact repeat positions" );
        M_INFO( logger, n_repeats << " repeats marked" );
        M_INFO( logger, n_rpos << " repeat positions" );
    }

    // create a histogram of logs of nmer frequencies
    //
    size_t fhist[65];
    std::fill( fhist, fhist + 65, 0 );
    cutoff_idx_ = 65;

    if( pij_data.size() > 0 )
    {
        StopWatch w( logger, "building nmer frequency histogram" );
        CProgress p( "building frequency histogram", "anchors",
                     ctx_.progress_flags_ );
        auto ph( p.GetTop() );
        ph.SetTotal( IDXMAP_SIZE - 1 );
        std::vector< CreateFreqHistJob::Data > jd( n_threads );
        std::atomic< uint32_t > task_idx( 0 );
        size_t job_idx( 0 );
        CTaskArray< CreateFreqHistJob > jobs(
            n_threads, *this, task_idx, &jd[0], job_idx, ph );
        p.Start();
        jobs.Start();
        p.Stop();

        for( auto const & data : jd )
        {
            for( size_t i( 0 ); i < 65; ++i )
            {
                fhist[i] += data.fhist[i];
            }
        }

        for( size_t i( 0 ); i < 65; ++i )
        {
            M_INFO( logger, "freq: " << i << "; count: " << fhist[i] );
        }

        constexpr size_t const NUM_WORDS_CUTOFF = 128*1024*1024ULL;
        // constexpr size_t const NUM_WORDS_CUTOFF = 1024*1024ULL;

        for( size_t num_words( 0 ); cutoff_idx_ > 0; )
        {
            if( num_words + fhist[--cutoff_idx_] > NUM_WORDS_CUTOFF ) break;
            num_words += fhist[cutoff_idx_];
        }

        M_INFO( logger, "cutoff frequency is " << cutoff_idx_ );
    }

    // create the log frequency table for repeated words
    //
    if( pij_data.size() > 0 )
    {
        StopWatch w( logger, "building nmer log frequency table" );
        CProgress p( "building log frequency table", "anchors",
                     ctx_.progress_flags_ );
        auto ph( p.GetTop() );
        ph.SetTotal( IDXMAP_SIZE - 1 );
        std::vector< CreateFreqTableJob::Data > jd( n_threads );
        std::atomic< uint32_t > task_idx( 0 );
        size_t job_idx( 0 );
        CTaskArray< CreateFreqTableJob > jobs(
            n_threads, *this, task_idx, &jd[0], job_idx, ph );
        p.Start();
        jobs.Start();

        for( auto & data : jd )
        {
            std::copy(
                data.freq_table_.begin(), data.freq_table_.end(),
                std::back_inserter( freq_table_ ) );
            std::vector< FreqTableEntry >().swap( data.freq_table_ );
        }

        std::sort( freq_table_.begin(), freq_table_.end() );
        p.Stop();
        used_mem_ += freq_table_.size()*sizeof( FreqTableEntry );
        M_INFO( logger, freq_table_.size() << " entries in frequency table" );
    }

    keep_loaded_ = true;
    return *this;
}

//------------------------------------------------------------------------------
CFastSeedsIndex & CFastSeedsIndex::Save( std::string const & basename )
{
    auto & logger( ctx_.logger_ );

    {
        auto fname( basename + LFM_EXT );
        std::ofstream ofs( fname, std::ios_base::binary );

        if( !ofs )
        {
            M_THROW( "error opening frequency map file " <<
                     fname << " for writing" );
        }

        ofs.exceptions( std::ios::badbit );
        uint64_t sz( freq_table_.size() );
        ofs.write( reinterpret_cast< char const * >( &sz ), sizeof( sz ) );
        ofs.write(
            reinterpret_cast< char const * >( freq_table_.data() ),
            sz*sizeof( FreqTableEntry ) );
        ofs.write(
            reinterpret_cast< char const * >( &cutoff_idx_ ),
            sizeof( uint64_t ) );
        M_INFO( logger, "saved log frequency data to " << fname );
    }

    {
        auto fname( basename + IDXMAP_EXT );
        std::ofstream ofs( fname.c_str(), std::ios_base::binary );

        if( !ofs )
        {
            M_THROW( "error opening index map file " <<
                     fname << " for writing" );
        }

        ofs.exceptions( std::ios::badbit );
        ofs.write( (char const *)&idxmap_[0],
                   IDXMAP_SIZE*sizeof( uint32_t ) );

        if( !ofs )
        {
            M_THROW( "error writing to index map file "<< fname );
        }

        M_INFO( logger, "saved index map to " << fname );
    }

    {
        auto fname( basename + IDX_EXT );
        std::ofstream ofs( fname.c_str(), std::ios_base::binary );

        if( !ofs )
        {
            M_THROW( "error opening index file " << fname << " for writing" );
        }

        ofs.exceptions( std::ios::badbit );
        CProgress p( "saving reference index", "chunks",
                     ctx_.progress_flags_ );
        auto ph( p.GetTop() );
        ph.SetTotal( idx_.size() );
        p.Start();

        for( auto const & chunk : idx_ )
        {
            chunk.Save( ofs );
            M_INFO( logger, "saved index chunk of size " <<
                            chunk.GetSize()*sizeof( IndexEntry ) <<
                            " bytes" );
            ph.Increment();
        }

        p.Stop();
    }

    return *this;
}

//------------------------------------------------------------------------------
CFastSeedsIndex & CFastSeedsIndex::Load(
        std::string const & basename, bool load_chunks )
{
    auto & logger( ctx_.logger_ );

    if( keep_loaded_ )
    {
        M_INFO( logger, "reference index is already loaded" )
        return *this;
    }

    // create reference offset data
    {
        ssize_t n_job_pos( 0 );
        roff_.resize( refs_.GetSize() + 1 );
        used_mem_ += roff_.size()*sizeof( uint32_t );

        for( size_t i( 0 ), ie( refs_.GetSize() ); i < ie; ++i )
        {
            roff_[i] = n_job_pos;
            n_job_pos += refs_.GetLength( i );
        }

        roff_[refs_.GetSize()] = n_job_pos;
    }

    // load index map data
    {
        idxmap_.clear();
        idxmap_.resize( IDXMAP_SIZE );
        auto fname( basename + IDXMAP_EXT );
        std::ifstream ifs( fname.c_str(), std::ios_base::binary );

        if( !ifs )
        {
            M_THROW( "error opening index map file " << fname <<
                     " for reading" );
        }

        ifs.read( (char *)&idxmap_[0], IDXMAP_SIZE*sizeof( uint32_t ) );
        used_mem_ += IDXMAP_SIZE*sizeof( uint32_t );

        if( !ifs )
        {
            M_THROW( "error reading from index map file " << fname );
        }

        M_INFO( ctx_.logger_, "loaded index map from " << fname );
    }

    SetUpChunks( false );

    // load index data
    {
        auto fname( basename + IDX_EXT );
        index_stream_.open( fname.c_str(), std::ios_base::binary );

        if( !index_stream_ )
        {
            M_THROW( "error opening index data file " << fname <<
                     " for reading" );
        }

        if( load_chunks )
        {
            for( auto & chunk : idx_ )
            {
                chunk.Load( index_stream_, &used_mem_ );
            }

            M_INFO( ctx_.logger_, "loaded index data from " << fname );
        }
    }

    // load frequency table
    {
        auto fname( basename + LFM_EXT );
        std::ifstream ifs( fname, std::ios_base::binary );

        if( !ifs )
        {
            M_INFO(
                ctx_.logger_,
                "error opening log frequency table file " << fname
                << " for reading; dynamic batch sizing will be disabled" );
        }
        else
        {
            ifs.exceptions( std::ios_base::badbit );
            uint64_t sz;
            ifs.read( reinterpret_cast< char * >( &sz ), sizeof( sz ) );
            freq_table_.clear();
            freq_table_.resize( sz );
            ifs.read(
                reinterpret_cast< char * >( freq_table_.data() ),
                sizeof( FreqTableEntry )*sz );
            used_mem_ += sizeof( FreqTableEntry )*sz;
            ifs.read(
                reinterpret_cast< char * >( &cutoff_idx_ ),
                sizeof( uint64_t ) );
            assert( cutoff_idx_ <= 65 );
            M_INFO( ctx_.logger_, "frequency cutoff idx: " << cutoff_idx_ );
            ctx_.dynamic_batches_ = true;
            M_INFO(
                ctx_.logger_, "loaded word frequency table from " << fname );
        }
    }

    if( load_chunks ) keep_loaded_ = true;
    return *this;
}

READFINDER_NS_END

