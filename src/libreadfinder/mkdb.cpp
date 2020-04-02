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

#include <map>

#include <boost/dynamic_bitset.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <zlib.h>

#include <libreadfinder/fsidx.hpp>
#include <libreadfinder/rf_ctx.hpp>
#include <libreadfinder/mkdb.hpp>
#include <libreadfinder/refdata.hpp>

#include <libtools/progress.hpp>

READFINDER_NS_BEGIN

//==============================================================================
//------------------------------------------------------------------------------
static char const * IDS_EXT = ".ids";
static char const * OFF_EXT = ".off";
static char const * DATA_EXT = ".dat";

//------------------------------------------------------------------------------
template< typename T_IO > 
static void CheckStream( T_IO & s, std::string const & name ) 
{
    if( !s.good() ) 
    {
        M_THROW( "failed stream: " << name );
    }
}

//==============================================================================
class CReadFinderHitsDBFactory
{

public:

    CReadFinderHitsDBFactory( CMkDBOptions const & opts,
                              CommonCtxP ctx );
    void Run();

private:

    void ProcessData();

    CommonCtxP ctx_;
    std::ifstream gzins_;
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf_;
    std::unique_ptr< std::istream > ref_is_;
    std::ofstream ids_os_,
                  offsets_os_;
    gzFile data_os_ = nullptr;
    std::string fasta_id_;
    std::vector< char > fasta_data_;
    std::vector< uint64_t> data_,
                           mask_;
    uint64_t offset_ = 0;
};

//------------------------------------------------------------------------------
CReadFinderHitsDBFactory::CReadFinderHitsDBFactory(
        CMkDBOptions const & opts, CommonCtxP ctx )
    : ctx_( ctx ),
      ids_os_( (opts.output + IDS_EXT).c_str() ),
      offsets_os_( (opts.output + OFF_EXT).c_str() ),
      data_os_( gzopen( (opts.output + DATA_EXT).c_str(), "wb9" ) )
{
    if( data_os_ == nullptr )
    {
        M_THROW( "can not open database compressed data stream " <<
                 opts.output << DATA_EXT << " for writing" );
    }

    if( opts.input.size() > 3 &&
        opts.input.substr( opts.input.size() - 3 ) == ".gz" )
    {
        gzins_.open( opts.input.c_str(),
                     std::ios_base::in | std::ios_base::binary );
        inbuf_.push( boost::iostreams::gzip_decompressor() );
        inbuf_.push( gzins_ );
        ref_is_.reset( new std::istream( &inbuf_ ) );
    }
    else
    {
        ref_is_.reset( new std::ifstream( opts.input.c_str() ) );
    }

    auto & ref_is( *ref_is_ );
    CheckStream( ref_is, opts.input );
    CheckStream( ids_os_, opts.output + IDS_EXT );
    CheckStream( offsets_os_, opts.output + OFF_EXT );
}

//------------------------------------------------------------------------------
void CReadFinderHitsDBFactory::ProcessData()
{
    data_.clear();
    mask_.clear();
    uint64_t len( 0 );
    data_.push_back( 0ULL ); data_.push_back( 0ULL );
    mask_.push_back( ~0ULL ); mask_.push_back( ~0ULL );

    for( std::string::size_type i( 0 ); i < fasta_data_.size(); )
    {
        uint64_t word( 0 ), mword( ~0ULL ), letter, mletter = 3;
        int j( 0 );

        for( ; j < 32 && i < fasta_data_.size(); ++j, ++i )
        {
            ++len;
            bool valid( true );
            letter = 0; mletter = 3;

            switch( fasta_data_[i] )
            {
                case 'A': case 'a':

                    letter = 0;
                    break;

                case 'C': case 'c':

                    letter = 1;
                    break;

                case 'G': case 'g':

                    letter = 2;
                    break;

                case 'T': case 't': case 'U': case 'u':

                    letter = 3;
                    break;

                case 'R': case 'r': case 'Y': case 'y': case 'S': case 's':
                case 'W': case 'w': case 'K': case 'k': case 'M': case 'm':
                case 'B': case 'b': case 'D': case 'd': case 'H': case 'h':
                case 'V': case 'v': case 'N': case 'n':

                    mletter = 0;
                    break;

                default:
                    
                    --j;
                    --len;
                    valid = false;
                    break;
            }

            if( valid )
            {
                letter <<= 2*j;
                mletter <<= 2*j;
                word |= letter;
                mword &= ~mletter;
            }
        }

        data_.push_back( word ); mask_.push_back( mword );
    }

    data_.push_back( 0ULL ); data_.push_back( 0ULL );
    mask_.push_back( ~0ULL ); mask_.push_back( ~0ULL );
    ids_os_ << fasta_id_ << std::endl;
    offsets_os_.write( (char const *)&offset_, sizeof( offset_ ) );
    offsets_os_.write( (char const *)&len, sizeof( len ) );

    if( data_.size() > 0 &&
        gzwrite( data_os_,
                 (char const *)&data_[0],
                 data_.size()*sizeof( uint64_t ) ) == 0 )
    {
        M_THROW( "failed write of " << data_.size() <<
                 " words of sequence data " );
    }

    if( mask_.size() > 0 &&
        gzwrite( data_os_,
                 (char const *)&mask_[0],
                 mask_.size()*sizeof( uint64_t ) ) == 0 )
    {
        M_THROW( "failed write of " << data_.size() <<
                 " words of mask data " );
    }

    offset_ += data_.size();
    offset_ += data_.size();
    fasta_data_.clear();
}

//------------------------------------------------------------------------------
void CReadFinderHitsDBFactory::Run()
{
    std::string line;
    int lineno( 0 );
    offset_ = 0;

    CCounterProgress p( "storing reference data", "sequences",
                        ctx_->progress_flags_ );
    auto ph( p.GetTop() );
    p.Start();
    auto & ref_is( *ref_is_ );

    while( ref_is )
    {
        std::getline( ref_is, line );
        ++lineno;

        if( line.empty() || line[0] == '#' )
        {
            continue;
        }

        if( line[0] == '>' )
        {
            if( !fasta_data_.empty() )
            {
                ProcessData();
                ph.Increment();
            }

            fasta_id_ = line.substr( 1, line.find_first_of( " \t" ) - 1 );

            if( fasta_id_.empty() )
            {
                M_THROW( "empty id at line " << lineno );
            }
        }
        else
        {
            std::copy( line.begin(), line.end(), 
                       std::back_inserter( fasta_data_ ) );
        }
    }

    if( !fasta_data_.empty() )
    {
        ProcessData();
        ph.Increment();
    }

    p.Stop();

    ids_os_.close();
    offsets_os_.close();
    gzclose( data_os_ );
}

//==============================================================================
namespace
{

//==============================================================================
class WordSource
{
public:

    typedef uint32_t Word;

    WordSource( TWord const * data, TSeqLen len );
    void operator++();
    Word operator*() const { return (Word)(w_&0xFFFFFFFFULL); }
    operator bool() const { return !done_; }

private:

    TWord const * data_;
    TWord w_ = 0,
          nw_ = 0;
    TSeqLen len_;
    TSeqOff off_ = 0;
    bool done_ = true;
};

//------------------------------------------------------------------------------
inline WordSource::WordSource( TWord const * data, TSeqLen len )
    :   data_( data ), len_( len )
{
    if( off_ >= len_ ) return;
    done_ = false;
    w_ = *data_++;
    nw_ = *data_++;
}

//------------------------------------------------------------------------------
inline void WordSource::operator++()
{
    if( ++off_ >= len_ )
    {
        done_ = true;
        return;
    }

    w_ >>= LB;
    w_ += (nw_<<((WL - 1)*LB));

    if( off_%WL == 0 ) nw_ = *data_++;
    else nw_ >>= LB;
}

//==============================================================================
//------------------------------------------------------------------------------
void CreateWorkSet( CMkDBOptions const & opts, CRefData const & refs )
{
    static char const * WS_EXT = ".ws";
    static size_t const WORD_SET_SIZE = 1ULL<<32;
    boost::dynamic_bitset< TWord > ws( WORD_SET_SIZE );

    for( size_t i( 0 ), ie( refs.GetSize() ); i < ie; ++i )
    {
        for( WordSource seq( refs.GetSeqData( i ), refs.GetLength( i ) ),
                        mask( refs.GetMaskData( i ), refs.GetLength( i ) );
                seq; ++seq, ++mask )
        {
            WordSource::Word mw( *mask );

            if( mw == 0 )
            {
                WordSource::Word sw( *seq ),
                                 rsw( sw );
                SEQ_NS::ReverseComplement< eNCBI2NA >( rsw );
                ws.set( sw );
                ws.set( rsw );
            }
        }
    }

    std::vector< TWord > blocks( ws.num_blocks(), 0 );
    boost::to_block_range( ws, blocks.begin() );
    auto fname( opts.output + WS_EXT );
    std::ofstream ofs( fname, std::ios_base::binary );

    if( !ofs )
    {
        M_THROW( "error opening index map file "<< fname << " for writing" );
    }

    ofs.exceptions( std::ios_base::badbit );
    ofs.write( (char const *)blocks.data(), sizeof( TWord )*blocks.size() );
}

}

//==============================================================================
//------------------------------------------------------------------------------
void MakeDB( CMkDBOptions const & opts )
{
    CommonCtxP ctx( new CCommonContext( opts ) );
    CReadFinderHitsDBFactory( opts, ctx ).Run();

    if( opts.mkidx || opts.mkws )
    {
        CRefData refs( opts.output );
        refs.LoadAll();
        if( opts.mkws ) CreateWorkSet( opts, refs );

        if( opts.mkidx )
        {
            CFastSeedsIndex( ctx, refs )
                .Create( opts.n_threads )
                .Save( opts.output );
        }
    }
}

READFINDER_NS_END

