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

#ifndef LIBREADFINDER_RF_OPTIONS_HPP
#define LIBREADFINDER_RF_OPTIONS_HPP

#include <vector>

#include <libreadfinder/defs.hpp>

#include <libseq/seqinput.hpp>

#include <libtools/log_handler.hpp>

READFINDER_NS_BEGIN

using namespace TOOLS_NS;
using namespace SEQ_NS;

//==============================================================================
/// Parameters common to all readfinder actions
///
struct CommonOptions
{
    std::string log_fname;  ///< Log file name.

    /// Severity threshold for logging.
    ///
    CLogHandler::Severity trace_level = CLogHandler::WARNING;

    size_t n_threads = 1;       ///< Number of threads to use.
    bool quiet = false;         ///< Disable progress reporting.

    friend std::ostream & operator<<(
            std::ostream & os, CommonOptions const & x )
    {
        return os << '\n'
                  <<    "COMMON OPTIONS: \n\n"
                  <<    "   log file name: " << x.log_fname << '\n'
                  <<    "   trace_level: " << x.trace_level << '\n'
                  <<    "   number of threads: "
                        << (x.n_threads > 0 ? std::to_string( x.n_threads )
                                            : "default")
                        << '\n'
                  <<    "   suppress progress reporting: " << x.quiet << '\n'
                  << std::endl;
    }
};

//==============================================================================
/// Parameters to use for readfinder database generation
///
struct CMkDBOptions : public CommonOptions
{
    std::string input;      ///< File name of the input FASTA file.
    std::string output;     ///< Base name of the database files.
    bool mkidx = false;     ///< Generate seeder index.
    bool mkws = false;      ///< Generate word set for pre-screeening.

    friend std::ostream & operator<<(
            std::ostream & os, CMkDBOptions const & x )
    {
        return os << (CommonOptions const &)x
                  <<    "DATABASE OPTIONS: \n\n"
                  <<    "   input: " << x.input << '\n'
                  <<    "   output: " << x.output << '\n'
                  <<    "   generate seeder index: " << x.mkidx << '\n'
                  <<    "   generate word set: " << x.mkws << '\n'
                  << std::endl;
    }
};

//==============================================================================
/// Parameters to use for readfinder search object construction.
///
struct CSearchOptions : public CommonOptions
{
    /// Application command line to use in SAM header.
    std::string cmd_line;

    /// Application version information to use in SAM header.
    std::string version_string;

    std::string db_name;        ///< Sequence database base name.

    /** Input data.

        Depending on the input format, this can be either FASTA file,
        or SRA accession.

        In the case of FASTA file, if both input_1 and input_2 are specified,
        then paired search is performed with input_1 and input_2 taken as the
        sources of the first and the second mates of the input reads.

        If only input_1 is specified and paired mode is not forced, then single
        search is performed.

        Otherwise, input_1 is the source of both mates. In this case the
        sequences in the file should alternate exactly between the first and the
        second mate (starting with the first mate). Mates of the same read must
        have either identical ids, or the should have the common prefix <id> and
        ids <id>[_./][12]. In the resulting SAM file, the reads are identified
        by the common prefix <id>.
    */
    std::string input_1;

    /// Mate input data (if applicable, see input_1).
    std::string input_2;

    std::string output;     ///< Name of the output (stdout, if empty).

    int input_format = CSeqInput::FASTA;    ///< Input format.

    size_t batch = 5000000, ///< Input batch size.
           start_read = 0,                                   ///< First read.
           end_read = std::numeric_limits< ssize_t >::max(); ///< Last read.

    /** Seeds that extend to  within this many bases from the start or end
        of the mate (excluding A-tails and T-heads) are considered for
        generating flanking regions for short exon search.
    */

    ///< Seeder min matched bases threshold
    TReadLen covered_bases = 0;

    ///< Max difference in diagonals between hits within a chain.
    TReadLen max_diag_delta = 4;

    ///< Seeder read coverage threshold.
    float coverage_th = 0.5f;

    bool per_mate_marks = false; ///< Report individual mates in seeder.
    bool pre_screen = false;     ///< Pre-screen reads.
    bool force_paired = false;   ///< Force paired for single fasta/fastq input.

    friend std::ostream & operator<<(
            std::ostream & os, CSearchOptions const & x )
    {
        return os << (CommonOptions const &)x
                  <<    "SEARCH OPTIONS: \n\n"
                  <<    "   application command line: " << x.cmd_line << '\n'
                  <<    "   database name: " << x.db_name << '\n'
                  <<    "   input: " << x.input_1 << '\n'
                  <<    "   mate input: " << x.input_2 << '\n'
                  <<    "   output: " << x.output << '\n'
                  <<    "   input format: " << x.input_format << std::endl;
    }
};

READFINDER_NS_END

#endif

