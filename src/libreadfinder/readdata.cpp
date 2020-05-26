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

#include <libreadfinder/readdata.hpp>

#include <libtools/progress.hpp>

READFINDER_NS_BEGIN

//------------------------------------------------------------------------------
CReadData::Reads const CReadData::Reads::End;

//------------------------------------------------------------------------------
inline void CReadData::ExtendBuffer( size_t len )
{
    static TWord const ZERO = 0;
    seq_data_.Append( &ZERO, WUNITS );
    len = WUNITS*( (len + WUNITS - 1)/WUNITS );
    seq_data_.resize( seq_data_.size() + len, ZERO );
    seq_data_.Append( &ZERO, WUNITS );
}

//------------------------------------------------------------------------------
inline void CReadData::ExtendBuffers( size_t len )
{
    static TWord const ZERO = 0;
    static TWord const ONE = CFullMask< TWord >::V;

    // one word padding in front of the sequence
    seq_data_.Append( &ZERO, WUNITS );
    mask_data_.Append( &ONE, WUNITS );

    // reserve enough full words for len bases.
    len = WUNITS*( (len + WUNITS - 1)/WUNITS );
    seq_data_.resize( seq_data_.size() + len, ZERO );
    mask_data_.resize( mask_data_.size() + len, ONE );

    // one word padding after the sequence
    seq_data_.Append( &ZERO, WUNITS );
    mask_data_.Append( &ONE, WUNITS );
}

//------------------------------------------------------------------------------
inline void CReadData::AppendData( 
        OrdId oid, std::string const & iupac, 
        EStrand strand, uint8_t mate_idx )
{
    typedef CFixedConstSeqView< eIUPACNA, uint8_t > SrcView;
    assert( oid < reads_.size() );
    auto & read( reads_[oid] );
    
    if( read.mates_[mate_idx].data[eFWD - 1] < 0 )
    {
        auto & mate( read.mates_[mate_idx] );
        mate.len = iupac.size();
        size_t orig_data_sz( seq_data_.size() ),
               orig_mask_sz( mask_data_.size() );
        assert( orig_data_sz%WUNITS == 0 );
        assert( orig_mask_sz%WUNITS == 0 );
        size_t data_off( orig_data_sz/WUNITS ),
               mask_off( orig_mask_sz/WUNITS );
        ExtendBuffers( mate.len );
        SeqView d1( seq_data_.GetBuf() + data_off, WUNITS, mate.len ),
                m1( mask_data_.GetBuf() + mask_off, WUNITS, mate.len );
        bool ambig( 
            Recode( 
                d1, m1, 
                SrcView( (uint8_t const *)iupac.data(), 0, mate.len ) ) );
        mate.data[strand - 1] = orig_data_sz + WUNITS;

        if( ambig )
        {
            ExtendBuffers( mate.len );
            SeqView d1( seq_data_.GetBuf() + data_off, WUNITS, mate.len ),
                    d2( seq_data_.GetBuf() + data_off, 
                        3*WUNITS + mate.len, mate.len ),
                    m1( mask_data_.GetBuf() + mask_off, WUNITS, mate.len ),
                    m2( mask_data_.GetBuf() + mask_off, 
                        3*WUNITS + mate.len, mate.len );
            Copy( d2, d1 );
            Copy( m2, m1 );
            ReverseComplement( d2 );
            Reverse( m2 );
            mate.mask[strand - 1] = orig_mask_sz + WUNITS;
            mate.mask[2 - strand] = orig_mask_sz + 3*WUNITS + mate.len;
            ++n_ambig_seq_;
        }
        else
        {
            if( mask_cache_[mate.len] < 0 )
            {
                mask_cache_[mate.len] = mate.mask[strand - 1] 
                                      = mate.mask[2 - strand] 
                                      = orig_mask_sz + WUNITS;
            }
            else
            {
                mask_data_.resize( orig_mask_sz );
                mate.mask[strand - 1] = mate.mask[2 - strand] 
                                      = mask_cache_[mate.len];
            }

            ExtendBuffer( mate.len );
            SeqView d1( seq_data_.GetBuf() + data_off, WUNITS, mate.len ),
                    d2( seq_data_.GetBuf() + data_off, 
                        3*WUNITS + mate.len, mate.len );
            Copy( d2, d1 );
            ReverseComplement( d2 );
        }

        mate.data[2 - strand] = orig_data_sz + 3*WUNITS + mate.len;
        ++n_seq_;
    }
    else
    {
        M_THROW( "duplicate mate for the same read " << oid );
    }
}

//------------------------------------------------------------------------------
auto CReadData::AddSeqData( 
        std::string const & id, std::string const & iupac, 
        EStrand strand, EMate mate, bool read_is_paired 
    ) -> OrdId
{
    if( !idset_ )
    {
        idset_.reset( new IdSet( IdSet::key_compare( &ids_ ) ) );
    }

    size_t sz( ids_.size() );
    ids_.resize( sz + id.size() + 1 );
    strcpy( &ids_[0] + sz, id.c_str() );
    IdHandle key{ sz, 0 };
    auto i( idset_->find( key ) );
    OrdId oid( idset_->size() );

    if( i == idset_->end() )
    {
        key.ordid = oid;
        i = idset_->insert( key ).first;
        reads_.push_back( Read( key.ordid ) );
        reads_.back().id_off_ = sz;
        reads_.back().read_is_paired_ = read_is_paired;
    }
    else if( !read_is_paired )
    {
        // M_WARN( logger_, "duplicate id: " << id );
        reads_.push_back( Read( oid ) );
        reads_.back().id_off_ = sz;
        reads_.back().read_is_paired_ = read_is_paired;
    }
    else
    {
        oid = i->ordid;
        ids_.resize( sz );
    }

    AppendData( oid, iupac, strand, FromMate( mate ) );
    return i->ordid;
}

//------------------------------------------------------------------------------
namespace
{
    std::pair< uint32_t, uint32_t > IUPACNA2NCBI2NA( char iupacna_base )
    {
        uint32_t base,
                 ambig = 0;

        switch( iupacna_base )
        {
            case 'A': case 'a': base = 0; break;
            case 'C': case 'c': base = 0x40000000UL; break;
            case 'G': case 'g': base = 0x80000000UL; break;
            case 'T': case 't': base = 0xC0000000UL; break;
            default: base = 0; ambig = 1; break;
        }

        return std::make_pair( base, ambig );
    }
}

bool CReadData::PreScreen(
    boost::dynamic_bitset< TWord > const & ws, std::string const & iupac )
{
    // a window is composed of 6 16-mers within a 21-mer.
    //
    static constexpr size_t const WLEN = 6;
    static constexpr size_t const NMER = 21;
    static constexpr uint8_t WINDOW_MASK = 0x3FU;
    static constexpr uint32_t NMER_MASK = 0x1FFFFFUL;

    if( iupac.size() < NMER ) return false;

    uint8_t window = 0;
    uint32_t ambigs = 0;
    size_t n_present = 0,
           n_ambigs = 0;
    uint32_t word( 0 );

    uint32_t base,
             ambig;
    
    size_t i( 0 );

    for( ; i < NMER - WLEN; ++i )
    {
        std::tie( base, ambig ) = IUPACNA2NCBI2NA( iupac[i] );
        ambigs <<= 1;
        ambigs += ambig;
        n_ambigs += ambig;
        word >>= 2;
        word += base;
    }

    for( ; i < NMER - 1; ++i )
    {
        std::tie( base, ambig ) = IUPACNA2NCBI2NA( iupac[i] );
        ambigs <<= 1;
        ambigs += ambig;
        n_ambigs += ambig;
        word >>= 2;
        word += base;
        unsigned int present( ws[word] );
        n_present += present;
        window <<= 1;
        window += present;
    }

    for( ; i < iupac.size(); ++i )
    {
        std::tie( base, ambig ) = IUPACNA2NCBI2NA( iupac[i] );
        ambigs <<= 1;
        ambigs &= NMER_MASK;
        ambigs += ambig;
        n_ambigs += ambig;
        word >>= 2;
        word += base;
        unsigned int present( ws[word] );
        n_present += present;
        window <<= 1;
        window &= WINDOW_MASK;
        window += present;
        if( n_ambigs == 0 && n_present == WLEN ) return true;
        n_ambigs -= (ambigs>>(NMER - 1));
        n_present -= (window>>(WLEN - 1));
    }

    return false;
}

//------------------------------------------------------------------------------
CReadData::CReadData( 
        CLogger & logger, CSeqInput & seqs,
        boost::dynamic_bitset< TWord > const * ws,
        size_t batch_size, int progress_flags
    )
    : logger_( logger ),
      mask_cache_( 1 + std::numeric_limits< TReadLen >::max(), -1 )
{
    static size_t const MAX_READ_LEN = 32*1024;

    size_t i( 0 ),
           n_skipped( 0 ),
           n_screened( 0 );
    CCounterProgress p( "reading input data", "reads", progress_flags );
    p.Start();

    for( auto && sd : seqs.Iterate() )
    {
        bool long_read( false );

        for( size_t mi( 0 ), mie( std::min( (size_t)2, sd.GetNCols() ) );
             mi < mie; ++mi )
        {
            if( sd.GetData( mi ).size() >= MAX_READ_LEN )
            {
                M_WARN( logger,
                        "mate " << mi << " of read " << sd.GetId() <<
                        " is too long (" << sd.GetData( mi ).size() <<
                        " bp); read will be skipped" );
                long_read = true;
                break;
            }
        }

        if( long_read )
        {
            ++n_skipped;
        }
        else
        {
            for( size_t mi( 0 ), mie( std::min( (size_t)2, sd.GetNCols() ) );
                 mi < mie; ++mi )
            {
                if( ws == nullptr || PreScreen( *ws, sd.GetData( mi ) ) )
                {
                    AddSeqData( sd.GetId(), sd.GetData( mi ),
                                eFWD, ToMate( mi ), (mie > 1) );
                }
                else
                {
                    ++n_screened;
                    AddSeqData( sd.GetId(), "",
                                eFWD, ToMate( mi ), (mie > 1) );
                }
            }
        }

        if( ++i >= batch_size )
        {
            break;
        }

        p.GetTop().Increment();
    }

    M_INFO( logger, n_skipped <<
                    " reads were skipped due to length restriction" );
    M_INFO( logger, n_screened <<
                    " mates were skipped due to pre-screening" );
    p.Stop();
    Freeze();
}

READFINDER_NS_END

