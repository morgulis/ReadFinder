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
#include <iostream>
#include <memory>
#include <thread>

#include <boost/program_options.hpp>

#include <libreadfinder/fast_seeds.hpp>
#include <libreadfinder/mkdb.hpp>
#include <libreadfinder/rf_options.hpp>
#include <libreadfinder/search.hpp>

#include <libseq/seqinput.hpp>

//==============================================================================

using namespace TOOLS_NS;
using namespace SEQ_NS;
using namespace READFINDER_NS;

//==============================================================================
namespace
{

static std::string const VERSION_STRING = "0.9.3";

namespace po = boost::program_options;

//------------------------------------------------------------------------------
std::map< std::string, int > INFMT_MAP = {
    { "fasta", CSeqInput::FASTA },
    { "zfasta", CSeqInput::ZFASTA },
    { "fastq", CSeqInput::FASTQ },
    { "zfastq", CSeqInput::ZFASTQ },
#ifdef USE_NGS
    { "sra", CSeqInput::SRA },
#endif
};

//==============================================================================
class Action
{
public:

    enum
    {
        MKDB = 0,
        SEARCH,

        N_ACTIONS
    };

    static std::string const NAMES[];

    virtual int Run() { return 0; }
};

//------------------------------------------------------------------------------
std::string const Action::NAMES[] = { "mkdb", "search" };

//==============================================================================
class MkDBAction : public Action
{
public:

    MkDBAction( CMkDBOptions const & opts ) : opts_( opts ) {}

    int Run() override;

private:

    CMkDBOptions opts_;
};

//------------------------------------------------------------------------------
int MkDBAction::Run()
{
    MakeDB( opts_ );
    return 0;
}

//==============================================================================
class SearchAction : public Action
{
public:

    SearchAction( CSearchOptions const & opts ) : opts_( opts ) {}

    int Run() override;

    CSearchOptions opts_;
};

//------------------------------------------------------------------------------
int SearchAction::Run()
{
    SearchSeeds( opts_ );
    return 0;
}

//==============================================================================
//------------------------------------------------------------------------------
std::shared_ptr< Action > ParseOptions( int argc, char ** argv )
{
    po::options_description global_options( "Global options" ),
                            common_options( "Options common to all actions" ),
                            mkdb_options( "Database creation options" ),
                            search_options( "Search options" );

    //==========================================================================
    // global options description
    //
    std::string requested_action;

    global_options.add_options()
        ( "help,h", "print usage information" )
        ( "version,v", "print application version and exit" )
        ( "action,A",
          po::value< std::string >( &requested_action )
          ->default_value( Action::NAMES[Action::SEARCH] ),
          "action to perform (possible values: mkdb, search) "
          "[default: search]" );

    //==========================================================================
    // common options description
    //
    CommonOptions copts;

    common_options.add_options()
        ( "log-file",
          po::value< std::string >( &copts.log_fname ),
          "log file name [default: <stdout>]" )
        ( "threads,t",
          po::value< size_t >( &copts.n_threads )->default_value(
              std::thread::hardware_concurrency() ),
          "number of threads to use [default: max available by the OS]" )
        ( "trace-level",
          po::value< TOOLS_NS::CLogHandler::Severity >( &copts.trace_level ),
          "minimum log severity level;\n"
          "one of { quiet, info, warning, error } "
          "[default: warning]")
        ( "quiet",
          po::bool_switch( &copts.quiet ),
          "do not report progress" );

    //==========================================================================
    // mkdb options description
    //
    CMkDBOptions mkdb_opts;

    mkdb_options.add_options()
        ( "input,i",
          po::value< std::string >( &mkdb_opts.input )->required(),
          "name of FASTA file containing reference sequence data; "
          "if the name ends in .gz, the file is assumed gzip-compressed" )
        ( "output,o",
          po::value< std::string >( &mkdb_opts.output )->required(),
          "base name of output database files" )
        ( "mkidx",
          po::bool_switch( &mkdb_opts.mkidx ),
          "enable creation of the database index" )
        ( "mkws",
          po::bool_switch( &mkdb_opts.mkws ),
          "enable creation of word set for pre-screening" );

    //==========================================================================
    // search options description
    //
    CSearchOptions search_opts;

    search_options.add_options()
        ( "db,d",
          po::value< std::string >( &search_opts.db_name )->required(),
          "base name of the reference database" )
        ( "keep-db-loaded",
          po::bool_switch( &search_opts.keep_loaded ),
          "keep the whole database index loaded in RAM" )
        ( "input,i",
          po::value< std::string >( &search_opts.input_1 ),
#ifdef USE_NGS
          "FASTA/FASTQ (possibly gzip-compressed) file name or SRA "
          "accession of the input reads source "
          "(or of the first mate in the case of paired search with "
#else
          "FASTA/FASTQ (possibly gzip-compressed) file name "
          "of the input reads source "
          "(or of the first mate in the case of paired search with "
#endif
          "FASTA/FASTQ input)" )
        ( "first-mate,1",
          po::value< std::string >( &search_opts.input_1 ),
          "alias for --input, -i" )
        ( "second-mate,2",
          po::value< std::string >( &search_opts.input_2 ),
          "Name of the FASTA/FASTQ file containing the second mates of the "
          "reads" )
        ( "in-fmt,F",
          po::value< std::string >()->notifier(
              [&search_opts]( std::string const & v )
              {
                    if( INFMT_MAP.count( v ) > 0 )
                    {
                        search_opts.input_format = INFMT_MAP[v];
                    }
                    else
                    {
                        M_THROW( "unknown input format: " << v );
                    }
              } ),
#ifdef USE_NGS
          "input file format {fasta,zfasta,fastq,zfastq,sra}; "
#else
          "input file format {fasta,zfasta,fastq,zfastq}; "
#endif
          "[default: fasta or zfasta if name ends in .gz]" )
        ( "memory,M",
          po::value< size_t >( &search_opts.max_mem ),
          "memory limit in megabytes "
          "[default: 90% of available machine memory]" )
        ( "batch-size,b",
          po::value< size_t >( &search_opts.batch ),
          "process reads in batches of this size "
          "(to prevent memory overflow) [default: 5000000]" )
        ( "first-read",
          po::value< size_t >( &search_opts.start_read ),
          "first read to process (0-based) [default: 0]" )
        ( "last-read",
          po::value< size_t >( &search_opts.end_read ),
          "last read to process (0-based) "
          "[default: number of reads - 1]" )
        ( "output,o",
          po::value< std::string >( &search_opts.output ),
          "output file name [default: <stdout>]" )
        ( "continuous-bases,C",
          po::value< TReadLen >( &search_opts.continuous_bases ),
          "length of min exact match needed to mark read as matched "
          "(must be >= 21) [default: 21]" )
        ( "covered-bases,B",
          po::value< TReadLen >( &search_opts.covered_bases ),
          "min number of matched bases needed to mark read as matched "
          "[default: 0]" )
        ( "max-diag-delta,D",
          po::value< TReadLen >( &search_opts.max_diag_delta ),
          "max difference in diagonals between hits within a chain"
          "[default: 4]" )
        ( "coverage-ratio,c",
          po::value< float >( &search_opts.coverage_th ),
          "coverage ratio threshold to mark read as matched "
          "[default: 0.5]" )
        ( "per-mate-marks,m",
          po::bool_switch( &search_opts.per_mate_marks ),
          "reporting on a per-mate basis" )
        ( "paired",
          po::bool_switch( &search_opts.force_paired ),
          "interpret single fasta ot fastq input as interleaved paired reads" )
        ( "prescreen",
          po::bool_switch( &search_opts.pre_screen ),
          "pre-screen reads for presense of words from reference" );

    //==========================================================================
    // parse and process global options
    //

    po::variables_map cmd_line_args;
    po::store( po::command_line_parser( argc, argv ).
                    options( global_options ).
                    style( po::command_line_style::default_style &
                           ~po::command_line_style::allow_guessing ).
                    allow_unregistered().
                    run(),
               cmd_line_args );
    po::notify( cmd_line_args );

    if( cmd_line_args.count( "version" ) > 0 )
    {
        std::cout << "readfinder version "
                  << VERSION_STRING << std::endl;
        return std::make_shared< Action >();
    }

    if( cmd_line_args.count( "help" ) > 0 )
    {
        std::cout << "readfinder [options] \n" << std::endl;
        po::options_description all_opts( "Available Program Options are" );
        all_opts.add( global_options )
                .add( common_options )
                .add( mkdb_options )
                .add( search_options );
        std::cout << all_opts << std::endl;
        return std::make_shared< Action >();
    }

    //==========================================================================
    common_options.add_options()
        ( "action,A", po::value< std::string >() );

    if( requested_action == Action::NAMES[Action::MKDB] )
    {
        common_options.add( mkdb_options );
        po::variables_map mkdb_args;
        po::store( po::command_line_parser( argc, argv ).
                        options( common_options ).
                        style( po::command_line_style::default_style &
                               ~po::command_line_style::allow_guessing ).
                        run(),
                   mkdb_args );
        po::notify( mkdb_args );
        (CommonOptions &)mkdb_opts = copts;
        return std::make_shared< MkDBAction >( mkdb_opts );
    }
    else if( requested_action == Action::NAMES[Action::SEARCH] )
    {
        common_options.add( search_options );
        po::variables_map search_args;
        po::store( po::command_line_parser( argc, argv ).
                        options( common_options ).
                        style( po::command_line_style::default_style &
                               ~po::command_line_style::allow_guessing ).
                        run(),
                   search_args );
        po::notify( search_args );
        search_opts.version_string = VERSION_STRING;
        (CommonOptions &)search_opts = copts;
        return std::make_shared< SearchAction >( search_opts );
    }

    M_THROW( "unknown action: " << requested_action );
}

}

//==============================================================================
//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    try
    {
        return ParseOptions( argc, argv )->Run();
    }
    catch( std::exception const & e )
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        std::cerr << "Run \"readfinder -h\" for usage information."
                  << std::endl;
        return EXIT_FAILURE;
    }

}

