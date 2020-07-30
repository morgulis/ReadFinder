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

#include <cassert>
#include <string>

#include <libreadfinder/refdata.hpp>

READFINDER_NS_BEGIN

static char const * IDMAP_EXT = ".ids";
static char const * OFFSETS_EXT = ".off";
static char const * DATA_EXT = ".dat";

//------------------------------------------------------------------------------
CRefData::CRefData( std::string const & db_name ) 
    : dbname_( db_name )
{
    // load reference sequence id map
    {
        std::string fname( dbname_ + IDMAP_EXT );
        std::ifstream is( fname.c_str() );
        if( !is.good() ) M_THROW( "could not open file " << fname );
        std::string line;

        while( is ) {
            getline( is, line );

            if( !line.empty() ) {
                id_map_.push_back( line );
            }
        }
    }

    // load offsets
    {
        std::string fname( dbname_ + OFFSETS_EXT );
        std::ifstream is( fname.c_str() );
        if( !is.good() ) M_THROW( "could not open file " << fname );
        uint64_t off, len;

        for( size_t i( 0 ), e( id_map_.size() ); i < e; ++i ) {
            is.read( (char *)&off, sizeof( uint64_t ) );
            is.read( (char *)&len, sizeof( uint64_t ) );

            if( !is.good() ) {
                M_THROW( "failed to read offsets at index " << i );
            }

            offsets_.push_back( off );
            lengths_.push_back( len );
        }
    }

    // open sequence data file
    std::string fname( dbname_ + DATA_EXT );

    if( (data_is_ = gzopen( fname.c_str(), "rb" )) == nullptr )
    {
        M_THROW( "can not open database compressed data stream " << fname <<
                 " for reading" );
    }

    data_.resize( GetSize() );
}

//------------------------------------------------------------------------------
void CRefData::Load( TRefOId refid )
{
    assert( refid < GetSize() );
    assert( data_[refid].get() == nullptr );
    ssize_t byteoff( offsets_[refid]*sizeof( TWord ) );

    if( gztell( data_is_ ) != byteoff )
    {
        auto coff( gzseek( data_is_, byteoff, SEEK_SET ) );

        if( coff != byteoff || coff < 0 )
        {
            M_THROW( "could not set position in compressed database for "
                     "sequence " << refid << " at offset " << byteoff );
        }
    }

    size_t n_words( 
            4 + (lengths_[refid] + WL - 1)/WL );
    data_[refid].reset( new SeqInfo );
    TData & seqdata( data_[refid]->seqdata ),
          & maskdata( data_[refid]->maskdata );
    seqdata.resize( n_words );
    maskdata.resize( n_words );

    auto bytes_read( gzread( data_is_,
                             (char *)&seqdata[0],
                             n_words*sizeof( TWord ) ) );

    if( (size_t)bytes_read != n_words*sizeof( TWord ) || bytes_read < 0 )
    {
        M_THROW( "read failed for sequence " << refid << " data" );
    }

    bytes_read = gzread( data_is_,
                         (char *)&maskdata[0],
                         n_words*sizeof( TWord ) );

    if( (size_t)bytes_read != n_words*sizeof( TWord ) || bytes_read < 0 )
    {
        M_THROW( "read failed for sequence " << refid << " mask data" );
    }
}

//------------------------------------------------------------------------------
size_t CRefData::GetMemUsage() const
{
    size_t res( 0 );
    res += id_map_.size()*sizeof( std::string );
    for( auto const & id : id_map_ ) res += id.size();
    res += offsets_.size()*sizeof( uint64_t );
    res += lengths_.size()*sizeof( TSeqLen );
    res += data_.size()*sizeof( SeqInfoHolder );

    for( auto const & d : data_ )
    {
        if( d != nullptr )
        {
            res += (d->seqdata.size() + d->maskdata.size());
        }
    }

    return res;
}

//------------------------------------------------------------------------------
#ifdef TEST_ENABLED

bool CRefData::Check( int test_id, void * data )
{
    switch( test_id )
    {
        case 0:
        {
            if( data_[0].get() != nullptr ) return false;
            Hold( 0 );
            if( data_[0].get() == nullptr ) return false;
            auto * ptr( data_[0].get() );
            if( data_[0]->count != 1 ) return false;
            auto const * sdata( &(data_[0]->seqdata)[0] ),
                       * mdata( &(data_[0]->maskdata)[0] );
            Hold( 0 );
            if( data_[0].get() == nullptr ) return false;
            if( data_[0].get() != ptr ) return false;
            if( data_[0]->count != 2 ) return false;
            if( &(data_[0]->seqdata)[0] != sdata ) return false;
            if( &(data_[0]->maskdata)[0] != mdata ) return false;
            Release( 0 );
            if( data_[0].get() == nullptr ) return false;
            if( data_[0].get() != ptr ) return false;
            if( data_[0]->count != 1 ) return false;
            if( &(data_[0]->seqdata)[0] != sdata ) return false;
            if( &(data_[0]->maskdata)[0] != mdata ) return false;
            Release( 0 );
            if( data_[0].get() != nullptr ) return false;
            return true;
        }

        default: M_THROW( "" );
    }

    return false;
}

#endif

READFINDER_NS_END

