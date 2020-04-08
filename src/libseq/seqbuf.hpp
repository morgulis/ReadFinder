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
/** \file libseq/seqbuf.hpp
    \brief Sequence buffers and operations on them.
*/

#ifndef LIBSEQ_SEQBUF_HPP
#define LIBSEQ_SEQBUF_HPP

#include <iterator>

#include <libseq/defs.hpp>

#include <libtools/databuf.hpp>

SEQ_NS_BEGIN

using namespace TOOLS_NS;

//==============================================================================
template< ECoding T_CODE, typename T_Word, typename T_Alloc >
class CFixedSeqBuf;

//==============================================================================
template< typename T_Alloc = std::allocator< char > >
class CSeqBuf
{
public:

    typedef T_Alloc Alloc;

    template< ECoding T_CODE, typename T_Word >
    CSeqBuf( CFixedSeqBuf< T_CODE, T_Word, T_Alloc > const & other );

    template< ECoding T_CODE, typename T_Word, typename T_OtherAlloc >
    CSeqBuf( CFixedSeqBuf< T_CODE, T_Word, T_OtherAlloc > const & other,
             T_OtherAlloc const & other_alloc = T_OtherAlloc() );

private:

    typedef CDataBuf< Alloc > DataBuf;

    DataBuf databuf_;
    ECoding coding_;
};

//==============================================================================
/** Sequence buffers with encoding and data word fixed at compile time.

    \tparam T_CODE  Sequence encoding.
    \tparam T_Word  Underlying word type (unsigned integral).
    \tparam T_Alloc Sequence data allocator.
*/
template< ECoding T_CODE,
          typename T_Word, 
          typename T_Alloc = std::allocator< T_Word > >
class CFixedSeqBuf
{
public:

    typedef T_Alloc Alloc;

    /// Code traits;
    typedef CCode< T_CODE > Code;

    static ECoding const CODING = T_CODE;

private:

    /// Data unit conrresponding to T_CODE.
    typedef CUnit< Code::LBITS > Unit;  

    /// Underlying data buffer type.
    typedef CFixedDataBuf< CWord< T_Word >, Unit::UNIT, Alloc > DataBuf;

public:

    /// Word traits.
    typedef typename DataBuf::Word Word;

    static EDataWord const WORD = Word::WORD;

    /// Instance constructor. Creates an empty buffer.
    CFixedSeqBuf( Alloc const & alloc = Alloc() );

    /// Clear and deallocate.
    void Reset();

    /// Clear without deallocation.
    void clear();

    /** Append data to the buffer.

        \param [in] words   Start of the source data.
        \param [in] n_bases Number of bases to append.
    */
    void Append( T_Word const * words, size_t n_bases );

    /** Resize the buffer, new bases undefined.

        \param [in] n_bases New size of the buffer.
    */
    void resize( size_t n_bases );

    /** Resize and set new bases to the given value.

        \param [in] n_bases New size of the buffer.
        \param [in] v       Value for new bases.
    */
    void resize( size_t n_bases, T_Word v );

    /// Get the size of the buffer in bases.
    size_t size() const;

    /// Get the minimum number of letters needed to bring data to a whole
    /// number of words.
    size_t GetPadding() const;

    /// Get the underlying data buffer.
    DataBuf const & GetDataBuf() const;

    T_Word * GetBuf() const;  ///< Get the underlying raw data.

private:

    DataBuf databuf_;   ///< Underlying data buffer.
};

//==============================================================================
//  IMPLEMENTATION
//==============================================================================
template< typename T_Alloc >
template< ECoding T_CODE, typename T_Word >
CSeqBuf< T_Alloc >::CSeqBuf( 
        CFixedSeqBuf< T_CODE, T_Word, T_Alloc > const & other )
    : databuf_( other.GetDataBuf() ),
      coding_( other.CODING )
{
}

//------------------------------------------------------------------------------
template< typename T_Alloc >
template< ECoding T_CODE, typename T_Word, typename T_OtherAlloc >
CSeqBuf< T_Alloc >::CSeqBuf( 
        CFixedSeqBuf< T_CODE, T_Word, T_OtherAlloc > const & other,
        T_OtherAlloc const & other_alloc )
    : databuf_( other.GetDataBuf(), other_alloc ),
      coding_( other.CODING )
{
}

//------------------------------------------------------------------------------
template< ECoding T_CODE, typename T_Word, typename T_Alloc >
inline CFixedSeqBuf< T_CODE, T_Word, T_Alloc >::CFixedSeqBuf(
        Alloc const & alloc )
    : databuf_( alloc )
{
}

//------------------------------------------------------------------------------
template< ECoding T_CODE, typename T_Word, typename T_Alloc >
inline void CFixedSeqBuf< T_CODE, T_Word, T_Alloc >::clear()
{
    databuf_.clear();
}

//------------------------------------------------------------------------------
template< ECoding T_CODE, typename T_Word, typename T_Alloc >
inline void CFixedSeqBuf< T_CODE, T_Word, T_Alloc >::Append( 
        T_Word const * words, size_t n_bases )
{
    databuf_.Append( words, n_bases );
}

//------------------------------------------------------------------------------
template< ECoding T_CODE, typename T_Word, typename T_Alloc >
void CFixedSeqBuf< T_CODE, T_Word, T_Alloc >::Reset()
{
    databuf_.Reset();
}

//------------------------------------------------------------------------------
template< ECoding T_CODE, typename T_Word, typename T_Alloc >
inline auto CFixedSeqBuf< T_CODE, T_Word, T_Alloc >::GetBuf() const -> T_Word *
{
    return databuf_.GetBuf();
}

//------------------------------------------------------------------------------
template< ECoding T_CODE, typename T_Word, typename T_Alloc >
inline auto CFixedSeqBuf< T_CODE, T_Word, T_Alloc >::GetDataBuf() const 
    -> DataBuf const &
{
    return databuf_;
}

//------------------------------------------------------------------------------
template< ECoding T_CODE, typename T_Word, typename T_Alloc >
inline void CFixedSeqBuf< T_CODE, T_Word, T_Alloc >::resize( size_t n_bases )
{
    databuf_.resize( n_bases );
}

//------------------------------------------------------------------------------
template< ECoding T_CODE, typename T_Word, typename T_Alloc >
inline void CFixedSeqBuf< T_CODE, T_Word, T_Alloc >::resize( 
        size_t n_bases, T_Word v )
{
    databuf_.resize( n_bases, v );
}

//------------------------------------------------------------------------------
template< ECoding T_CODE, typename T_Word, typename T_Alloc >
inline size_t CFixedSeqBuf< T_CODE, T_Word, T_Alloc >::size() const
{
    return databuf_.size();
}

//------------------------------------------------------------------------------
template< ECoding T_CODE, typename T_Word, typename T_Alloc >
inline size_t CFixedSeqBuf< T_CODE, T_Word, T_Alloc >::GetPadding() const
{
    return databuf_.GetPadding();
}

SEQ_NS_END

#endif

