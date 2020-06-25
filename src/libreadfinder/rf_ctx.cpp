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

#include <libreadfinder/rf_ctx.hpp>

#include <unistd.h>

READFINDER_NS_BEGIN

//==============================================================================
//------------------------------------------------------------------------------
CCommonContext::CCommonContext( CommonOptions const & opts )
    : logger_( "readfinder", CLogger::LogHandler( 
                    opts.log_fname.empty() ? 
                        new CStreamLogHandler( "<stderr>", &std::cerr ) :
                        new CFileLogHandler( opts.log_fname ) ) )
{
    logger_.SetSeverity( opts.trace_level );

    if( opts.quiet )
    {
        progress_flags_ |= CProgress::QUIET;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
CSearchContext::CSearchContext( CSearchOptions const & opts )
    : CCommonContext( opts ), CSearchOptions( opts )
{
    Check();
    ResetOutStream();

    if( !output.empty() )
    {
        ResetOutStream( output );
    }

    std::vector< std::string > fnames;

    if( !opts.input_1.empty() )
    {
        fnames.push_back( opts.input_1 );
    }

    if( !opts.input_2.empty() )
    {
        fnames.push_back( opts.input_2 );
    }

    seqs.reset( MkSeqInput( fnames,
                            opts.input_format,
                            opts.start_read,
                            opts.end_read - opts.start_read,
                            opts.force_paired ) );

    if( max_mem == 0 )
    { // compute memory limit based on physical memory size
        long pages( sysconf( _SC_PHYS_PAGES ) ),
             page_size( sysconf( _SC_PAGE_SIZE ) );
        max_mem_bytes = pages*page_size;
    }
    else max_mem_bytes = max_mem*1024ULL*1024ULL;

    refs.reset( new CRefData( opts.db_name ) );
    refs->LoadAll();

    if( pre_screen )
    {
        static char const * WS_EXT = ".ws";
        static size_t const WORD_SET_SIZE = 1ULL<<32;
        ws.resize( WORD_SET_SIZE );
        std::vector< TWord > blocks( ws.num_blocks(), 0 );
        auto fname( db_name + WS_EXT );
        std::ifstream ifs( fname, std::ios_base::binary );

        if( !ifs )
        {
            M_THROW( "error opening word set file " << fname <<
                     " for reading" );
        }

        ifs.exceptions( std::ios_base::badbit );
        ifs.read( (char *)blocks.data(), sizeof( TWord )*blocks.size() );
        boost::from_block_range( blocks.begin(), blocks.end(), ws );
    }

    fsidx.reset( new CFastSeedsIndex( *this, *refs ) );
}

//------------------------------------------------------------------------------
void CSearchContext::Check()
{
}

READFINDER_NS_END

