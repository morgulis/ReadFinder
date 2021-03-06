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

/** \file libseq/fastq_input.hpp
    \brief Single stream fastq input.
*/

#ifndef LIBSEQ_FASTQ_INPUT_HPP
#define LIBSEQ_FASTQ_INPUT_HPP

#include <fstream>

#include <libseq/seqinput.hpp>
#include <libseq/defs.hpp>

#include <libtools/exception.hpp>
#include <libtools/textfile_input.hpp>

SEQ_NS_BEGIN

using namespace TOOLS_NS;

//==============================================================================
class CFastqInput : public CSeqInput
{
public:

    CFastqInput( std::string const & fname, bool compressed,
                 size_t start, ssize_t n_seq, bool paired = false );
    virtual ~CFastqInput() override {}

    virtual bool Next() override;
    virtual SeqData const & GetSeqData() const override;

private:

    bool NextPriv( size_t col );

    std::shared_ptr< CTextInput > in_;
    SeqData sd_;
    std::string line_;
    ssize_t n_left_;
};

//------------------------------------------------------------------------------
inline auto CFastqInput::GetSeqData() const -> SeqData const &
{
    return sd_;
}

SEQ_NS_END

#endif

