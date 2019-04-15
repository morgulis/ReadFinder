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

READFINDER_NS_BEGIN

//==============================================================================
//------------------------------------------------------------------------------
CCommonContext::CCommonContext( CommonOptions const & opts )
    : logger_( "readfinder", CLogger::LogHandler( 
                    opts.log_fname.empty() ? 
                        new CStreamLogHandler( "<stderr>", &std::cerr ) :
                        new CFileLogHandler( opts.log_fname ) ) )// ,
      // memmgr_( new CMemMonitor( opts.max_mem ) )
{
    logger_.SetSeverity( opts.trace_level );

    if( opts.quiet )
    {
        progress_flags_ |= CProgress::QUIET;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
CSearchContext::CSearchContext( CSearchOptions const & opts/* , bool init_output  */)
    : CCommonContext( opts ), CSearchOptions( opts )// ,
      // run_second_pass( !no_splice && !one_pass )
{
    /*
    penalties.Finalize();

    if( penalties.intron_len_penalty_exp_scale <= 0.0 )
    {
        penalties.intron_len_penalty_exp_scale =
            log( 2 )/max_internal_region;
        assert( penalties.intron_len_penalty_exp_scale > 0.0 );
    }
    */

    Check();
    ResetOutStream();

    if( !output.empty() )
    {
        ResetOutStream( output );
    }

    /*
    auto & os( *osp );

    if( init_output )
    {
        os << "@HD\tVN:1.0\tGO:query\n";
        os << "@PG\tID:readfinder\tPN:readfinder\tVN:" << opts.version_string
           << "\tCL:" << opts.cmd_line << '\n';
    }
    */

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
                            opts.end_read - opts.start_read ) );
    // stat.reset( new SearchStat );
    // refs.reset( new CRefData( opts.db_name ) );
    // (*stat)[MHStat::N_REFS] = refs->GetSize();

    /*
    if( init_output )
    {
        for( size_t i( 0 ); i < refs->GetSize(); ++i )
        {
            os << "@SQ\tSN:" << refs->GetRefId( i )
               << "\tLN:" << refs->GetLength( i ) << '\n';
        }
    }
    */

    refs->LoadAll();
}

//------------------------------------------------------------------------------
void CSearchContext::Check()
{
}

READFINDER_NS_END

