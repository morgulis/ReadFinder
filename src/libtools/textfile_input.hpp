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
/*! \file libtools/textfile_input.hpp CTextFileInput class definition. */

#ifndef LIBTOOLS_TEXTFILE_INPUT_HPP
#define LIBTOOLS_TEXTFILE_INPUT_HPP

#include <fstream>
#include <string>

#include <libtools/textinput.hpp>
#include <libtools/defs.hpp>

TOOLS_NS_BEGIN

//------------------------------------------------------------------------------
/*! Line based text input handler for files. */
class CTextFileInput : public CTextInput
{
public:

    /*! Instance constructor.
        
        Creates a C++ input file stream from the given name and constructs
        the base CTextInput object from that stream. \c \b fname is also
        used as an identifier for CTextInput base.

        \param [in] fname Input file name.
    */
    CTextFileInput( std::string const & fname );
};

//==============================================================================
// IMPLEMENTATION
//==============================================================================
inline CTextFileInput::CTextFileInput( std::string const & fname )
    : CTextInput( new std::ifstream( fname.c_str() ), fname )
{}

TOOLS_NS_END

#endif

