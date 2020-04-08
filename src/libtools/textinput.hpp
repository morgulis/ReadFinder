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
/*! \file libtools/textinput.hpp CTextInput class definition. */

#ifndef LIBTOOLS_TEXTINPUT_HPP
#define LIBTOOLS_TEXTINPUT_HPP

#include <iostream>
#include <iterator>
#include <memory>
#include <string>

#include <libtools/exception.hpp>
#include <libtools/defs.hpp>

TOOLS_NS_BEGIN

//------------------------------------------------------------------------------
/*! Simple text input handler.

    Provides checked line based input from a C++ input stream. Provides
    access to current line number. Conforms to input iterator concept
    and provides interface to support 'ranged for' construct.

    A text id can be defined for the stream to simplify referencing it
    in human readable messages.
*/
class CTextInput 
    : public std::iterator< std::input_iterator_tag, std::string, void >
{
private:

    /*! Temporary structure used to hold value for postfix operator++. */
    struct ValueHolder
    {
        /*! Instance constructor.

            \param [in] v Value to hold.
        */
        ValueHolder( value_type const & v );

        /*! Get the held value.

            \return The value being held.
        */
        value_type operator*() const;

        value_type const v_;    ///< The value being held.
    };

    /*! Check the object for consistency.

        \throws std::runtime_error if no input stream is set up or if
                                   the input stream is in failed state,
                                   but not at the end of file.
    */
    void Check();

    /*! Read next line from the input stream.

        If the input stream is not set up or is at the end of file, then
        makes sets boolean state of the object to \c \b false.
        Otherwise, if the input stream is not in the failed state, then 
        reads the next line from it. The state of the input stream is
        checked before and after reading, and if it fails, the boolean
        state of the object is set to \c \b false.
    */
    void Next();

public:

    /*! Type used to hold shared ownership of underlying input stream.
    */
    typedef std::shared_ptr< std::istream > Stream;

    /*! CTextInput object that is always in the boolean state of
        \c \b false. Used as a 'one past the end' iterator value.
    */
    static CTextInput const END;

    /*! Instance constructor.

        Shares the input stream ownership with \c \b s. Reads the first
        line of the stream, if possible.

        \param [in] s Holds the input stream that should be used by this instance.
        \param [in] id Text identifier.
    */
    CTextInput( Stream s, std::string const & id = "" );

    /*! Instance constructor.

        Assumes ownership of its argument. Reads the first
        line of the stream, if possible.

        \param [in] s Pointer to an input stream.
        \param [in] id Text identifier.
    */
    CTextInput( std::istream * s = nullptr, std::string const & id = "" );

    /*! Instance constructor.

        Creates an empty object (no actual stream is set up) with a given
        identifier.

        \param [in] id Text identifier.
    */
    CTextInput( std::string const & id );

    /*! Replace the input stream.

        Use the new input stream. Shared the new input stream ownership 
        with \c \b s.

        \param [in] s New input stream.

        \return The old input stream.
    */
    Stream Reset( Stream s );

    /*! Replace the input stream.

        Assumes ownership of \c \b s.

        \param [in] s New input stream.

        \return The old input stream.
    */
    Stream Reset( std::istream * s = nullptr );

    /*! Get the input stream that is currently installed.

        \return Currently installed input stream.
    */
    Stream GetStream() const;

    /*! Get the current line number in the input stream.

        The line number correponds to the number of lines read by
        this instance. It may not correspond to the actual position
        in the input stream if it was not at the beginning, when
        it was installed.

        \return The current line number.
    */
    size_t GetLineNo() const;

    /*! Skip all blank lines and lines starting with a given character.

        \param [in] pfx Character starting a comment line.
    */
    void SkipBlanksAndComments( char pfx );

    /*! Get the stream identifier.

        \return The stream identifier.
    */
    std::string const & GetId() const;

    /*! Set the stream identifier.

        \param [in] id The new identifier.
    */
    void SetId( std::string const & id );

    /*! Get the boolean status of the object.

        \return \c \b true if the stream is installed, good, and not
                at the end of file; \c \b false otherwise.
    */
    operator bool() const;

    /*! Get the line at the current position. 

        \return Current input line.
    */
    reference operator*();

    /*! Get the pointer to the line at the current position. 

        \return Pointer to the current input line.
    */
    pointer operator->();

    /*! Prefix increment. Reads the next line of input.

        \return This object.
    */
    CTextInput & operator++();

    /*! Postfix increment. Reads the next line of input.

        \return An object holding the previous line of input.
    */
    ValueHolder operator++( int );

    /*! No operation.

        This method does not rewind to the beginning of the stream. It
        is definied solely to support 'range-based for' C++ construct.

        \return This object.
    */
    CTextInput & begin();

    /*! Get a value that compares equal to the input stream objects at
        \c \b false boolean state.

        This method is defined to support 'range-based for' C++ construct.

        \return END
    */
    CTextInput const & end();

    /*! Comparison operator for CTextInput objects.

        Two CTextInput instances compare equal only if both are in the
        \c \b false boolean state.

        \param [in] a Left object.
        \param [in] b Right object.

        \return \c \b true, if the two object compare equal, \c \b false
                otherwise.
    */
    friend bool operator==( CTextInput const & a, CTextInput const & b );

protected:

    Stream s_;          ///< Holds the pointer to the installed input stream.
    std::string id_;    ///< Text identifier.
    value_type v_;      ///< Current input line.
    size_t lineno_ = 0; ///< Current line number.
    bool good_ = true;  ///< Current boolean state.
};

//------------------------------------------------------------------------------
/*! Compare CTextInput objects for inequality.

    \param [in] a Left object.
    \param [in] b Right object.

    \return \c \b false if the two objects are equal, \c \b true otherwise.
*/
inline bool operator!=( CTextInput const & a, CTextInput const & b )
{
    return !(a == b);
}

//==============================================================================
// IMPLEMENTATION
//==============================================================================
inline CTextInput::ValueHolder::ValueHolder( value_type const & v )
    : v_( v )
{}

//------------------------------------------------------------------------------
inline CTextInput::value_type CTextInput::ValueHolder::operator*() const
{
    return v_;
}

//------------------------------------------------------------------------------
inline void CTextInput::Check()
{
    if( s_ && s_->fail() && !s_->eof() )
    {
        good_ = false;
        M_THROW( "text input stream " << (id_.empty() ? "" : id_ + " ") << 
                 "failure at line " << lineno_ );
    }
}

//------------------------------------------------------------------------------
inline void CTextInput::Next()
{
    if( s_ && s_->eof() )
    {
        good_ = false;
    }
    else if( s_ && good_ )
    {
        Check();
        std::getline( *s_, v_ );
        ++lineno_;
        Check();
    }
}

//------------------------------------------------------------------------------
inline CTextInput::CTextInput( Stream s, std::string const & id ) 
    : s_( s ), id_( id )
{
    Next();
}

//------------------------------------------------------------------------------
inline CTextInput::CTextInput( std::istream * s, std::string const & id ) 
    : CTextInput( Stream( s ), id )
{}

//------------------------------------------------------------------------------
inline CTextInput::CTextInput( std::string const & id )
    : CTextInput( nullptr, id )
{}

//------------------------------------------------------------------------------
inline CTextInput::Stream CTextInput::Reset( Stream s )
{
    good_ = true;
    auto res( s_ );
    s_ = s;
    lineno_ = 0;
    Next();
    return res;
}

//------------------------------------------------------------------------------
inline CTextInput::Stream CTextInput::Reset( std::istream * s )
{
    good_ = true;
    auto res( s_ );
    s_.reset( s );
    lineno_ = 0;
    Next();
    return res;
}

//------------------------------------------------------------------------------
inline CTextInput::Stream CTextInput::GetStream() const 
{ 
    return s_; 
}

//------------------------------------------------------------------------------
inline size_t CTextInput::GetLineNo() const
{
    return lineno_;
}

//------------------------------------------------------------------------------
inline std::string const & CTextInput::GetId() const
{
    return id_;
}

//------------------------------------------------------------------------------
inline void CTextInput::SetId( std::string const & id )
{
    id_ = id;
}

//------------------------------------------------------------------------------
inline CTextInput::operator bool() const
{
    return s_ && good_ && !s_->fail();
}

//------------------------------------------------------------------------------
inline CTextInput::reference CTextInput::operator*()
{
    return v_;
}

//------------------------------------------------------------------------------
inline CTextInput::pointer CTextInput::operator->()
{
    return &v_;
}

//------------------------------------------------------------------------------
inline CTextInput & CTextInput::operator++()
{
    Next();
    return *this;
}

//------------------------------------------------------------------------------
inline CTextInput::ValueHolder CTextInput::operator++( int )
{
    ValueHolder res( v_ );
    Next();
    return res;
}

//------------------------------------------------------------------------------
inline CTextInput & CTextInput::begin()
{
    return *this;
}

//------------------------------------------------------------------------------
inline CTextInput const & CTextInput::end()
{
    return END;
}

//------------------------------------------------------------------------------
inline void CTextInput::SkipBlanksAndComments( char pfx )
{
    for( ; *this && (v_.empty() || v_[0] == pfx); Next() );
}

//------------------------------------------------------------------------------
inline bool operator==( CTextInput const & a, CTextInput const & b )
{
    if( !a && !b ) 
    {
        return true;
    }

    return false;
}

TOOLS_NS_END

#endif

