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

/** \file libseq/fasta_input.hpp
    \brief Single stream fasta input.
*/

#ifndef LIBSEQ_FASTA_INPUT_HPP
#define LIBSEQ_FASTA_INPUT_HPP

#include <fstream>

#include <libseq/seqinput.hpp>
#include <libseq/defs.hpp>

#include <libtools/exception.hpp>
#include <libtools/textfile_input.hpp>
#include <libtools/textinput.hpp>

SEQ_NS_BEGIN

using namespace TOOLS_NS;

//==============================================================================
class CFastaInput : public CSeqInput
{
public:

    CFastaInput( std::string const & fname, bool compressed,
                 size_t start = 0, ssize_t n_seq = -1 );
    virtual ~CFastaInput() override {}

    virtual bool Next() override;
    virtual SeqData const & GetSeqData() const override;

private:

    std::shared_ptr< CTextInput > in_;
    SeqData sd_;
    std::string line_;
    ssize_t n_left_;
};

//------------------------------------------------------------------------------
inline auto CFastaInput::GetSeqData() const -> SeqData const &
{
    return sd_;
}

SEQ_NS_END

#endif

