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
#include <memory>
#include <set>
#include <vector>

#include <libreadfinder/defs.hpp>

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

public:

    struct Read;

    CReadData();
    CReadData( CLogger & logger,
               CSeqInput & seqs,
               boost::dynamic_bitset< TWord > const * ws,
               size_t batch_size,
               int progress_flags = 0
             );

    OrdId AddSeqData( std::string const & id,
                      std::string const & iupac,
                      EStrand strand,
                      EMate mate, bool read_is_paired
                    );

    Read const & operator[]( OrdId oid ) const;

    void Freeze();

    size_t GetNReads() const;
    size_t GetNSeq() const;
    size_t GetNPaired() const;
    size_t GetNSingle() const;
    size_t GetNAmbigSeq() const;

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

    typedef CFixedSeqBuf< CODING, TWord > SeqBuf;

public:

    typedef CFixedSeqView< CODING, TWord > SeqView;
    typedef CFixedConstSeqView< CODING, TWord > SeqConstView;

    SeqConstView GetSeqData( OrdId oid, EMate mate, EStrand strand ) const;
    SeqConstView GetMaskData( OrdId oid, EMate mate, EStrand strand ) const;

    TReadLen GetLen( OrdId oid, EMate mate ) const;
    TReadLen GetATail( OrdId oid, EMate mate, EStrand strand ) const;
    TReadLen GetTHead( OrdId oid, EMate mate, EStrand strand ) const;

    char const * GetReadId( OrdId oid ) const;

    /// \name Interface for 'range-based for' iteration.
    /// @{
    Reads begin();
    Reads const & end() const;
    /// @}

private:

    typedef std::vector< ssize_t > MaskCache;

    void AppendData( OrdId oid, std::string const & iupac, 
                     EStrand strand, uint8_t mate_idx );
    void ExtendBuffers( size_t len );
    void ExtendBuffer( size_t len );

    bool PreScreen(
        boost::dynamic_bitset< TWord > const & ws, std::string const & iupac );

    std::unique_ptr< IdSet > idset_;

    Ids ids_;
    ReadData reads_;
    SeqBuf seq_data_,
           mask_data_;
    MaskCache mask_cache_;
    size_t n_seq_ = 0,
           n_ambig_seq_ = 0;
};

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
inline CReadData::CReadData()
    : mask_cache_( 1 + std::numeric_limits< TReadLen >::max(), -1 )
{
}

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
inline TReadLen CReadData::GetATail( 
        OrdId oid, EMate mate, EStrand strand ) const
{
    return reads_[oid].mates_[FromMate( mate )].A_tail[strand - 1];
}

//------------------------------------------------------------------------------
inline TReadLen CReadData::GetTHead( 
        OrdId oid, EMate mate, EStrand strand ) const
{
    return reads_[oid].mates_[FromMate( mate )].A_tail[2 - strand];
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

