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

#ifndef LIBREADFINDER_READDATA_HPP
#define LIBREADFINDER_READDATA_HPP

#include <cstring>
#include <deque>
#include <memory>
#include <set>
#include <vector>

#include <libreadfinder/defs.hpp>
#include <libreadfinder/fsidx.hpp>

#include <libseq/seqbuf.hpp>
#include <libseq/seqinput.hpp>
#include <libseq/sequtil.hpp>
#include <libseq/seqview.hpp>

#include <libtools/logger.hpp>

#include <boost/dynamic_bitset.hpp>

READFINDER_NS_BEGIN

using namespace SEQ_NS;
using namespace TOOLS_NS;

//------------------------------------------------------------------------------
class CReadData
{
private:

    class Reads;

    struct WordData
    {
        static constexpr size_t const ANCHOR_BITS = 14;
        static constexpr size_t const NMER_BITS =
            LB*CFastSeedsDefs::NMER_BASES;
        static constexpr size_t const WORD_BITS = NMER_BITS - ANCHOR_BITS;
        static constexpr TWord const WORD_MASK = (1ULL<<WORD_BITS) - 1;

        static_assert( WORD_BITS <= BYTE_BITS*sizeof( uint32_t ), "" );

        uint32_t word,
                 oid;
    };

    typedef std::deque< WordData > WordSet;
    typedef std::vector< WordSet > Words;

    static constexpr size_t const N_WORD_SETS = (1ULL<<WordData::ANCHOR_BITS);

    static_assert( sizeof( WordData ) == 8, "" );

public:

    class HashWordSource;
    class HashMaskSource;

    struct Read;

    CReadData( CLogger & logger,
               CSeqInput & seqs,
               boost::dynamic_bitset< TWord > const * ws,
               CFastSeedsIndex const & fsidx,
               size_t batch_size,
               size_t mem_limit,
               int progress_flags = 0
             );

    OrdId AddSeqData( std::string const & id,
                      std::string const & iupac,
                      EStrand strand,
                      EMate mate, bool read_is_paired,
                      size_t & bases_read, size_t & mem_used,
                      Words & words,
                      CFastSeedsIndex const & fsidx
                    );

    Read const & operator[]( OrdId oid ) const;

    void Freeze();

    size_t GetNReads() const;
    size_t GetNSeq() const;
    size_t GetNPaired() const;
    size_t GetNSingle() const;
    size_t GetNAmbigSeq() const;

    OrdId GetStartOId() const { return start_oid_; }
    OrdId GetEndOId() const { return end_oid_; }

private:

    typedef std::vector< char > Ids;
    typedef std::vector< Read > ReadData;

    struct IdHandle
    {
        Ids::size_type id;
        OrdId ordid;
    };

    struct IdHandleCompare
    {
        IdHandleCompare( Ids const * ids );
        bool operator()( IdHandle const & l, IdHandle const & r ) const;

        Ids const * ids;
    };

    typedef std::set< IdHandle, IdHandleCompare > IdSet;

    static ECoding const CODING = eNCBI2NA;

    typedef CCode< CODING > Code;
    typedef CWord< TWord > Word;
    typedef CUnit< Code::UNIT > Unit;

    static size_t const WUNITS = Unit::WordUnits< Word >::N_UNITS;
    static size_t const WBITS = WUNITS*Unit::N_BITS;
    static size_t const WBYTES = WBITS/BYTE_BITS;

    typedef CFixedSeqBuf< CODING, TWord > SeqBuf;

public:

    typedef CFixedSeqView< CODING, TWord > SeqView;
    typedef CFixedConstSeqView< CODING, TWord > SeqConstView;

    SeqConstView GetSeqData( OrdId oid, EMate mate, EStrand strand ) const;
    SeqConstView GetMaskData( OrdId oid, EMate mate, EStrand strand ) const;

    TReadLen GetLen( OrdId oid, EMate mate ) const;

    void Update()
    {
        start_oid_ = end_oid_;
    }

    char const * GetReadId( OrdId oid ) const;

    /// \name Interface for 'range-based for' iteration.
    /// @{
    Reads begin();
    Reads const & end() const;
    /// @}

private:

    typedef std::vector< ssize_t > MaskCache;

    size_t AppendData(
        OrdId oid, std::string const & iupac,
        EStrand strand, uint8_t mate_idx,
        CFastSeedsIndex const & fsidx, Words & words );
    void CollectWords( OrdId oid, uint8_t mate_idx, Words & words );
    void ExtendBuffers( size_t len );
    void ExtendBuffer( size_t len );

    bool PreScreen(
        boost::dynamic_bitset< TWord > const & ws, std::string const & iupac );

    size_t EstimateMemory( CFastSeedsIndex const & fsidx, Words & words );

    std::unique_ptr< IdSet > idset_;

    TOOLS_NS::CLogger & logger_;
    Ids ids_;
    ReadData reads_;
    SeqBuf seq_data_,
           mask_data_;
    MaskCache mask_cache_;
    size_t n_seq_ = 0,
           n_ambig_seq_ = 0;
    OrdId start_oid_ = 0,
          end_oid_ = 0;
    std::deque< uint32_t > word_freq_;
};

//==============================================================================
class CReadData::HashWordSource : public CFastSeedsDefs
{
public:

    HashWordSource( TWord const * seq_data, TSeqLen len, TSeqOff off );
    void operator++();

    uint32_t GetAnchor() const { return data_.f.anchor; }
    uint32_t GetWord() const { return data_.f.word; }
    TSeqOff GetHashOff() const { return off_ - off_adj_; }
    operator bool() const { return !done_; }

private:

    struct _Data
    {
        uint64_t GetWord() const { return f.word; }
        uint64_t GetAnchor() const { return f.anchor; }
        uint64_t GetNMer() const { return (GetAnchor()<<WORD_BITS) + GetWord(); }
        TWord GetData() const { return w; }

        union
        {
            struct
            {
                uint64_t word   : WORD_BITS;
                uint64_t anchor : ANCHOR_BITS;
            } f;

            TWord w;
        };

        _Data() {}
        _Data( TWord ww ) { w = ww; }
    } data_;

    static_assert( sizeof( data_ ) == 8, "" );

public:

    typedef _Data Data;

    Data & GetData() { return data_; }
    Data const & GetData() const { return data_; }

private:

    TWord nw_;
    TWord const * seq_data_;
    TSeqLen len_;
    TSeqOff off_,
            off_adj_;
    bool done_ = true;
};

//------------------------------------------------------------------------------
inline CReadData::HashWordSource::HashWordSource(
        TWord const * seq_data, TSeqLen len, TSeqOff off )
    : nw_( 0 ), seq_data_( seq_data ), len_( off + len ),
      off_( off ), off_adj_( off )
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
inline void CReadData::HashWordSource::operator++()
{
    TSeqOff woff( off_++%WL );
    
    if( off_ >= len_ )
    {
        done_ = true;
        return;
    }

    ++woff;
    data_.w >>= LB;
    data_.w += (nw_<<((WL - 1)*LB));

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
        nw_ >>= LB;
    }
}

//==============================================================================
class CReadData::HashMaskSource : public CFastSeedsDefs
{
public:

    HashMaskSource( TWord const * mask_data, TSeqLen len, TSeqOff off );
    void operator++();

    uint64_t GetNMer() const { return data_.f.nmer; }
    operator bool() const { return !done_; }

private:

    struct
    {
        union
        {
            struct
            {
                uint64_t nmer : NMER_BITS;
            } f;

            TWord w;
        };
    } data_;

    static_assert( sizeof( data_ ) == 8, "" );

    TWord nw_ = 0;
    TWord const * mask_data_;
    TSeqLen len_;
    TSeqOff off_ = 0;
    bool done_ = true;
};

//------------------------------------------------------------------------------
inline CReadData::HashMaskSource::HashMaskSource(
        TWord const * mask_data, TSeqLen len, TSeqOff off )
    : mask_data_( mask_data ), len_( off + len ), off_( off )
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
inline void CReadData::HashMaskSource::operator++()
{
    TSeqOff woff( off_++%WL );
    
    if( off_ >= len_ )
    {
        done_ = true;
        return;
    }

    ++woff;
    data_.w >>= LB;
    data_.w += (nw_<<((WL - 1)*LB));

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
        nw_ >>= LB;
    }
}
//------------------------------------------------------------------------------
struct CReadData::Read
{
    Read( OrdId oid = 0 );
    OrdId GetOrdId() const;
    bool IsPaired() const;
    bool ReadIsPaired() const { return read_is_paired_; }

    struct MateInfo
    {
        TReadLen len;
        ssize_t data[MAX_STRANDS - 1];
        ssize_t mask[MAX_STRANDS - 1];
        TReadLen A_tail[MAX_STRANDS - 1];
    };

    OrdId oid_;
    MateInfo mates_[FromMate( MAX_MATES )];
    size_t id_off_;
    bool read_is_paired_;
};

//==============================================================================
class CReadData::Reads
{
public:

    static Reads const End;

    Reads() {}
    Reads( CReadData const & o );

    /// \name Interface for 'range-based for' iteration.
    /// @{
    Read const & operator*() const;
    Reads & operator++();

    friend bool operator!=( Reads const & x, Reads const & y )
    {
        return x.o_ != y.o_;
    }
    /// @}

private:

    CReadData const * o_ = nullptr;
    ReadData::const_iterator curr_;
};

//------------------------------------------------------------------------------
inline CReadData::Reads::Reads( CReadData const & o )
    : o_( &o ), curr_( o_->reads_.begin() )
{
}

//------------------------------------------------------------------------------
inline auto CReadData::Reads::operator*() const -> Read const &
{
    return *curr_;
}

//------------------------------------------------------------------------------
inline auto CReadData::Reads::operator++() -> Reads &
{
    assert( o_ != nullptr );
    
    if( ++curr_ == o_->reads_.end() )
    {
        o_ = nullptr;
    }

    return *this;
}

//==============================================================================
// IMPLEMENTATION
//==============================================================================
//------------------------------------------------------------------------------
inline auto CReadData::begin() -> Reads
{
    return Reads( *this );
}

//------------------------------------------------------------------------------
inline auto CReadData::end() const -> Reads const &
{
    return Reads::End;
}

//------------------------------------------------------------------------------
inline auto CReadData::operator[]( OrdId oid ) const -> Read const &
{
    return reads_[oid];
}

//------------------------------------------------------------------------------
inline void CReadData::Freeze()
{
    idset_.reset();
}

//------------------------------------------------------------------------------
inline auto CReadData::GetSeqData( 
        OrdId oid, EMate mate, EStrand strand ) const -> SeqConstView
{
    auto const & m( reads_[oid].mates_[FromMate( mate )] );
    size_t off( m.data[strand - 1]/SeqConstView::WUNITS );
    return SeqView( seq_data_.GetBuf() + off, 
                    (m.data[strand - 1])%SeqConstView::WUNITS, m.len );
}

//------------------------------------------------------------------------------
inline auto CReadData::GetMaskData( 
        OrdId oid, EMate mate, EStrand strand ) const -> SeqConstView
{
    auto const & m( reads_[oid].mates_[FromMate( mate )] );
    size_t off( m.mask[strand - 1]/SeqConstView::WUNITS );
    return SeqView( mask_data_.GetBuf() + off, 
                    m.mask[strand - 1]%SeqConstView::WUNITS, m.len );
}

//------------------------------------------------------------------------------
inline TReadLen CReadData::GetLen( OrdId oid, EMate mate ) const
{
    return reads_[oid].mates_[FromMate( mate )].len;
}

//------------------------------------------------------------------------------
inline char const * CReadData::GetReadId( OrdId oid ) const
{
    return &ids_[0] + reads_[oid].id_off_;
}

//------------------------------------------------------------------------------
inline size_t CReadData::GetNReads() const 
{
    return reads_.size();
}

//------------------------------------------------------------------------------
inline size_t CReadData::GetNSeq() const 
{
    return n_seq_;
}

//------------------------------------------------------------------------------
inline size_t CReadData::GetNPaired() const 
{
    return GetNSeq() - GetNReads();
}

//------------------------------------------------------------------------------
inline size_t CReadData::GetNSingle() const 
{
    return GetNReads() - GetNPaired();
}

//------------------------------------------------------------------------------
inline size_t CReadData::GetNAmbigSeq() const
{
    return n_ambig_seq_;
}

//------------------------------------------------------------------------------
inline CReadData::IdHandleCompare::IdHandleCompare( Ids const * ids )
    : ids( ids )
{}

//------------------------------------------------------------------------------
inline bool CReadData::IdHandleCompare::operator()( 
        CReadData::IdHandle const & l, CReadData::IdHandle const & r ) const
{
    char const * data( &(*ids)[0] );
    return strcmp( data + l.id, data + r.id ) < 0;
}

//------------------------------------------------------------------------------
inline CReadData::Read::Read( OrdId oid )
    : oid_( oid )
{
    mates_[0] = mates_[1] = MateInfo{ 0, { -1, -1 }, { -1, -1 } };
}

//------------------------------------------------------------------------------
inline auto CReadData::Read::GetOrdId() const -> OrdId
{
    return oid_;
};

//------------------------------------------------------------------------------
inline bool CReadData::Read::IsPaired() const
{
    return mates_[0].data[eFWD - 1] >= 0 && mates_[1].data[eFWD - 1] >= 0;
}

READFINDER_NS_END

#endif

