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

/** \file libseq/seqinput.cpp
    \brief Base class and factory for reading sequence data represented in
           different formats.
*/

#include <libseq/fasta_input.hpp>
#include <libseq/fastq_input.hpp>
#include <libseq/paired_input.hpp>
#include <libseq/seqinput.hpp>

#ifdef USE_NGS
#   include <libseq/sra_input.hpp>
#endif

#include <libtools/exception.hpp>

SEQ_NS_BEGIN

//------------------------------------------------------------------------------
CSeqInput * MkSeqInput( std::vector< std::string > const & fnames, int fmt,
                        size_t start, ssize_t n_seq, bool paired )
{
    bool compressed( fmt == CSeqInput::ZFASTA || fmt == CSeqInput::ZFASTQ );

    if( fnames.size() < 2 )
    {
        if( fmt == CSeqInput::FASTA || fmt == CSeqInput::ZFASTA )
        {
            return new CFastaInput(
                fnames[0], compressed, start, n_seq, paired );
        }

        if( fmt == CSeqInput::FASTQ || fmt == CSeqInput::ZFASTQ )
        {
            return new CFastqInput(
                fnames[0], compressed, start, n_seq, paired );
        }

#ifdef USE_NGS
        if( fmt == CSeqInput::SRA )
        {
            return new CSRAInput( fnames[0], start, n_seq );
        }
#endif
    }
    else
    {
        if( fmt == CSeqInput::FASTA || fmt == CSeqInput::ZFASTA )
        {
            return CreatePairedInput( 
                    new CFastaInput( fnames[0], compressed, start, n_seq ),
                    new CFastaInput( fnames[1], compressed, start, n_seq ) );
        }

        if( fmt == CSeqInput::FASTQ || fmt == CSeqInput::ZFASTQ )
        {
            return CreatePairedInput(
                    new CFastqInput( fnames[0], compressed, start, n_seq ),
                    new CFastqInput( fnames[1], compressed, start, n_seq ) );
        }
    }

    M_THROW( "failed to create input" );
}

SEQ_NS_END

