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

/** \file libseq/fastq_input.cpp
    \brief Single stream fastq input.
*/

#include <libseq/fastq_input.hpp>

#include <libtools/gzinfile.hpp>

SEQ_NS_BEGIN

//------------------------------------------------------------------------------
CFastqInput::CFastqInput(
        std::string const & fname, bool compressed,
        size_t start, ssize_t n_seq )
    : in_( (compressed || CheckCompressed( fname )) ?
           (CTextInput *)(new CGZipInput( fname )) :
           (CTextInput *)(new CTextFileInput( fname )) ),
      sd_( 1 ), n_left_( n_seq )
{
    for( ; start > 0 && Next(); --start );
}

//------------------------------------------------------------------------------
bool CFastqInput::Next()
{
    CTextInput & in( *in_ );

    if( !in || n_left_ == 0 )
    {
        return false;
    }

    line_ = *in;
    assert( !line_.empty() );

    if( line_[0] != '@' )
    {
        M_THROW( in.GetId() <<
                 ": expected @defline at line " << in.GetLineNo() );
    }

    {
        auto pos( line_.find_first_of( " \t" ) );

        if( pos == std::string::npos )
        {
            sd_.SetId( line_.substr( 1 ) );
        }
        else
        {
            sd_.SetId( line_.substr( 1, pos - 1 ) );
        }
    }

    std::vector< char > data;

    while( in )
    {
        ++in;

        if( in )
        {
            line_ = *in;
            assert( !line_.empty() );

            if( line_[0] == '+' )
            {
                break;
            }

            data.insert( data.end(), line_.begin(), line_.end() );
        }
    }

    if( !in )
    {
        M_THROW( in.GetId() <<
                 ": expected +defline at line " << in.GetLineNo() );
    }

    size_t len( 0 );

    while( in )
    {
        ++in;

        if( in )
        {
            line_ = *in;
            assert( !line_.empty() );

            if( line_[0] == '@' && len == data.size() )
            {
                break;
            }
            else
            {
                len += line_.size();
            }
        }
    }

    sd_.SetData( 0, std::string( data.begin(), data.end() ) );
    --n_left_;
    return true;
}

SEQ_NS_END

