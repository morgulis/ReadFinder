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

/** \file libtools/gzinfile.hpp
    \brief GZip encrypted input stream based on boost::iostreams.
*/

#include <fstream>
#include <iostream>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <libtools/textinput.hpp>
#include <libtools/defs.hpp>

TOOLS_NS_BEGIN

//==============================================================================
class CGZipInput : public CTextInput
{
public:

    CGZipInput( std::string const & fname );

private:

    std::ifstream gzins_;
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf_;
    std::shared_ptr< std::istream > ins_;
};

//------------------------------------------------------------------------------
inline CGZipInput::CGZipInput( std::string const & fname )
    : CTextInput( fname ),
      gzins_( fname.c_str(), std::ios_base::in | std::ios_base::binary )
{
    inbuf_.push( boost::iostreams::gzip_decompressor() );
    inbuf_.push( gzins_ );
    ins_.reset( new std::istream( &inbuf_ ) );
    Reset( ins_ );
}

TOOLS_NS_END
