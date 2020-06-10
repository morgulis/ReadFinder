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

#ifndef LIBREADFINDER_FSIDX_HPP
#define LIBREADFINDER_FSIDX_HPP

#include <libreadfinder/common_ctx.hpp>
#include <libreadfinder/refdata.hpp>
#include <libreadfinder/fsdefs.hpp>
#include <libreadfinder/defs.hpp>

READFINDER_NS_BEGIN

//==============================================================================
class CFastSeedsIndex : public CFastSeedsDefs
{
public:

    struct IndexEntry
    {
        WordData wd;
        uint32_t pos;

        friend std::ostream & operator<<(
                std::ostream & os, IndexEntry const & ie )
        {
            return os << "{ wd: " << ie.wd << "; pos: " << ie.pos << " }";
        }
    };

    static_assert( sizeof( CFastSeedsIndex::IndexEntry ) == 8, "" );

private:

    static size_t const IDXMAP_SIZE = 1ULL + (1ULL<<ANCHOR_BITS);
    static size_t const MIN_IDX_CHUNK_SIZE = 128*1024*1024ULL;

    typedef std::vector< uint32_t > RefOffsets;
    typedef std::vector< uint32_t > IndexMap;

    class IndexChunk
    {
    private:

        void Check( uint32_t anchor ) const
        {
            assert( loaded_ );
            assert( anchor >= start_anchor_ && anchor < end_anchor_ );
        }

    public:

        IndexChunk() {}

        IndexChunk( CCommonContext & ctx, IndexMap const & idxmap,
                    size_t start_anchor, size_t end_anchor, size_t sz );

        IndexChunk( CCommonContext & ctx, IndexMap const & idxmap,
                    size_t start_anchor, size_t end_anchor );

        void Load( std::ifstream & is );
        void Unload( CCommonContext & ctx );
        void Save( std::ostream & os ) const;

        IndexEntry * begin( uint32_t anchor )
        {
            Check( anchor );
            return &data_[0] + ((*idxmap_)[anchor] - data_off_);
        }

        IndexEntry const * cbegin( uint32_t anchor ) const
        {
            Check( anchor );
            return &data_[0] + ((*idxmap_)[anchor] - data_off_);
        }

        IndexEntry * end( uint32_t anchor )
        {
            Check( anchor );
            return &data_[0] + ((*idxmap_)[anchor + 1] - data_off_);
        }

        IndexEntry const * cend( uint32_t anchor ) const
        {
            Check( anchor );
            return &data_[0] + ((*idxmap_)[anchor + 1] - data_off_);
        }

        void SetData( size_t off, IndexEntry const & ie )
        {
            assert( loaded_ );
            assert( off >= data_off_ );
            assert( off - data_off_ < data_.size() );
            data_[off - data_off_] = ie;
        }

        size_t GetSize() const
        {
            assert( loaded_ );
            return data_.size();
        }

        uint32_t GetStartAnchor() const { return start_anchor_; }
        uint32_t GetEndAnchor() const { return end_anchor_; }
        bool IsLoaded() const { return loaded_; }

    private:

        typedef std::vector< IndexEntry > Data;

        uint32_t start_anchor_,
                 end_anchor_;
        IndexMap const * idxmap_ = nullptr;
        size_t data_off_;
        Data data_;
        bool loaded_ = false;
    };

    typedef std::vector< IndexChunk > Index;
    typedef std::vector< uint32_t > ChunkMap;

    struct PopulateIndexJobData;
    struct IndexNMerCountingJob;
    struct PopulateIndexJob;
    struct SortIndexJob;

    class IndexEntrySource;
    class IndexMaskSource;
    class WordSource;

public:

    CFastSeedsIndex( CCommonContext & ctx, CRefData const & refs );

    CFastSeedsIndex & Create( size_t n_threads );
    CFastSeedsIndex & Save( std::string const & basename );
    CFastSeedsIndex & Load( std::string const & basename,
                            bool load_chunks = false );

    void Unload()
    {
        if( !keep_loaded_ )
        {
            Index().swap( idx_ );
            IndexMap().swap( idxmap_ );
            index_stream_.close();
            M_INFO( ctx_.logger_, "Index data unloaded" );
        }
    }

    IndexEntry * begin( uint32_t anchor )
    {
        assert( anchor < IDXMAP_SIZE - 1 );
        return idx_[chunk_map_[anchor]].begin( anchor );
    }

    IndexEntry const * cbegin( uint32_t anchor ) const
    {
        assert( anchor < IDXMAP_SIZE - 1 );
        return idx_[chunk_map_[anchor]].cbegin( anchor );
    }

    IndexEntry * end( uint32_t anchor )
    {
        assert( anchor < IDXMAP_SIZE - 1 );
        return idx_[chunk_map_[anchor]].end( anchor );
    }

    IndexEntry const * cend( uint32_t anchor ) const
    {
        assert( anchor < IDXMAP_SIZE - 1 );
        return idx_[chunk_map_[anchor]].cend( anchor );
    }

    std::pair< uint32_t, uint32_t > GetAbsPosRange( TRefOId refid ) const
    {
        assert( refid < roff_.size() - 1 );
        return std::make_pair( roff_[refid], roff_[refid + 1] );
    }

    size_t GetNumChunks() const { return idx_.size(); }

    bool ChunkIsLoaded( size_t chunk ) const
    {
        assert( chunk < idx_.size() );
        return idx_[chunk].IsLoaded();
    }

    std::pair< uint32_t, uint32_t > GetChunkAnchorRange( size_t chunk ) const
    {
        assert( chunk < idx_.size() );
        auto & cd( idx_[chunk] );
        return std::make_pair( cd.GetStartAnchor(), cd.GetEndAnchor() );
    }

    void Load( size_t chunk )
    {
        assert( chunk < idx_.size() );
        idx_[chunk].Load( index_stream_ );
    }

    void Unload( size_t chunk )
    {
        assert( chunk < idx_.size() );

        if( !keep_loaded_ )
        {
            idx_[chunk].Unload( ctx_ );
            M_INFO( ctx_.logger_, "chunk " << chunk << " is unloaded" );
        }
    }

private:

    void SetUpChunks( bool allocate = true );

    CCommonContext & ctx_;
    CRefData const & refs_;
    RefOffsets roff_;
    IndexMap idxmap_;
    ChunkMap chunk_map_;
    Index idx_;
    std::ifstream index_stream_;
    std::vector< FreqTableEntry > freq_table_;
    bool keep_loaded_ = false;
};

READFINDER_NS_END

#endif

