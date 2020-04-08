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
/** \file libtools/databuf.hpp
    \brief Data buffers and operations on them.
*/

#ifndef LIBTOOLS_DATABUF_HPP
#define LIBTOOLS_DATABUF_HPP

#include <memory>
#include <vector>

#include <libtools/bitutils.hpp>
#include <libtools/defs.hpp>

TOOLS_NS_BEGIN

//------------------------------------------------------------------------------
/** Raw data holder.

    Holds a vector of 64-bit words allocated via a specified allocator.

    \tparam T_Alloc Allocator type satisfying allocator concept requirements of 
                    C++-11.
*/
template< typename T_Alloc = std::allocator< char > >
class CRawDataBuf
{
public:

    typedef uint64_t Word;  ///< Buffer word type.

private:

    /// Allocator type adjusted to allocate data in words.
    typedef typename T_Alloc::template rebind< Word >::other Alloc;

    typedef std::vector< Word > Data;   ///< Container type to hold the data.

public:

    /** Instance constructor.

        Allocates a new buffer of a given size. Contents of the buffer
        is undefined.
        
        \param [in] sz  Initial size of the buffer in words.
        \param [in] a   Allocator instance.
    */
    CRawDataBuf( size_t sz = 0, T_Alloc const & a = T_Alloc() );

    /** Instance constructor.

        Copy constructor. The new instance shares data with \c other.

        \param [in] other   Source buffer.
    */
    CRawDataBuf( CRawDataBuf const & other );

    /** Instance constructor.

        Copy constructor for buffers using different allocators. The data
        is copied over, since sharing across different allocators is not
        feasible.

        \tparam T_OtherAlloc    The allocator type used by the source buffer.

        \param [in] other   The source buffer.
        \param [in] a       Allocator instance.
    */
    template< typename T_OtherAlloc >
    CRawDataBuf( CRawDataBuf< T_OtherAlloc > const & other,
                 T_Alloc const & a = T_Alloc() );

    void resize( size_t new_sz );   ///< Resize to the new number of words.

    /// Resize with value of extra words derived from the given value.
    template< typename T_Word > void resize( size_t new_sz, T_Word v );

    size_t size() const;    ///< Get the number of words in the buffer.
    void clear();   ///< Clear the buffer.

    Word * GetBuf() const;  ///< Get the actual data buffer.
    Data const & GetData() const;   ///< Get the actual data as a vector.

    void Reset();   ///< Clear the buffer and deallocate.

private:

    std::shared_ptr< Data > data_;  ///< Actual data holder.
};

//------------------------------------------------------------------------------
template< typename T_Word, EDataUnit T_UNIT, typename T_Alloc >
class CFixedDataBuf;

template< typename T_Alloc = std::allocator< char > >
class CDataBuf
{
public:

    template< typename T_Word, EDataUnit T_UNIT >
    CDataBuf( CFixedDataBuf< T_Word, T_UNIT, T_Alloc > const & other );

    template< typename T_Word, EDataUnit T_UNIT, typename T_OtherAlloc >
    CDataBuf( CFixedDataBuf< T_Word, T_UNIT, T_OtherAlloc > const & other,
              T_OtherAlloc const & other_alloc = T_OtherAlloc() );

private:

    typedef CRawDataBuf< T_Alloc > RawBuf;

    EDataWord word_;
    EDataUnit unit_;
    RawBuf rawbuf_;
    size_t sz_ = 0;
};

//------------------------------------------------------------------------------
/** Data buffer for compile-time known word and unit sizes.
  
    \tparam T_Word  Word traits type.
    \tparam T_UNIT  Data unit type name.
    \tparam T_Alloc Allocator type satisfying allocator concept requirements of 
                    C++-11.
*/
template< typename T_Word, 
          EDataUnit T_UNIT, 
          typename T_Alloc = std::allocator< typename T_Word::T > >
class CFixedDataBuf
{
private:

    typedef CRawDataBuf< T_Alloc > RawBuf;  ///< Underlying raw buffer type.

public:

    typedef T_Alloc Alloc;
    typedef T_Word Word;
    typedef CUnit< T_UNIT > Unit;

    /** Instance constructor.

        Creates an empty buffer.

        \param [in] alloc   Allocator instance.
    */
    CFixedDataBuf( Alloc const & alloc = Alloc() );

    void Reset();   ///< Clear buffer and deallocate.

    void clear();   ///< Clear buffer without deallocation.

    /// Resize the buffer to new number of units.
    void resize( size_t n_units );

    /// Resize the buffer to new number of units and initialize
    /// extra words with a given value.
    void resize( size_t n_units, typename Word::T v );

    size_t size() const;    ///< Get the size of the buffer in units.

    /// Get the minimum number of units needed to bring data to a whole
    /// number of words.
    size_t GetPadding() const;

    /** Append content from raw storage.

        \param [in] words   Raw data.
        \param [in] n_units Number of units to append.
    */
    void Append( typename Word::T const * words, size_t n_units );

    typename Word::T * GetBuf() const;  ///< Get the underlying raw data.

    RawBuf GetRawBuf() const;   ///< Get the underlying raw buffer.

private:

    static size_t const WORD_RATIO = sizeof( typename RawBuf::Word )/Word::WORD;
    static size_t const UNIT_RATIO = Word::N_BITS/Unit::N_BITS;
    static size_t const RAW_UNIT_RATIO = WORD_RATIO*UNIT_RATIO;

    RawBuf rawbuf_; ///< Underlying raw data buffer.
    size_t sz_ = 0; ///< Size of this buffer in units.
};

//==============================================================================
//  IMPLEMENTATION
//==============================================================================
template< typename T_Alloc >
CRawDataBuf< T_Alloc >::CRawDataBuf( size_t sz, T_Alloc const & a )
    : data_( new Data( sz, a ) )
{
}

//------------------------------------------------------------------------------
template< typename T_Alloc >
CRawDataBuf< T_Alloc >::CRawDataBuf( CRawDataBuf const & other )
    : data_( other.data_ )
{
}

//------------------------------------------------------------------------------
template< typename T_Alloc >
template< typename T_OtherAlloc >
CRawDataBuf< T_Alloc >::CRawDataBuf( 
        CRawDataBuf< T_OtherAlloc > const & other, T_Alloc const & a )
    : data_( new Data( other.size(), a ) )
{
    memcpy( GetBuf(), other.GetBuf(), size()*sizeof( Word ) );
}

//------------------------------------------------------------------------------
template< typename T_Alloc >
inline void CRawDataBuf< T_Alloc >::resize( size_t new_sz )
{
    data_->resize( new_sz );
}

//------------------------------------------------------------------------------
template< typename T_Alloc >
template< typename T_Word > 
inline void CRawDataBuf< T_Alloc >::resize( size_t new_sz, T_Word v )
{
    static size_t const N_WORDS = sizeof( Word )/sizeof( T_Word );

    if( new_sz > data_->size() )
    {
        Word w;
        T_Word * ww( reinterpret_cast< T_Word * >( &w ) );

        for( size_t i( 0 ); i < N_WORDS; ++i )
        {
            ww[i] = v;
        }

        data_->resize( new_sz, w );
    }
    else
    {
        data_->resize( new_sz );
    }
}

//------------------------------------------------------------------------------
template< typename T_Alloc >
inline size_t CRawDataBuf< T_Alloc >::size() const
{
    return data_->size();
}

//------------------------------------------------------------------------------
template< typename T_Alloc >
inline auto CRawDataBuf< T_Alloc >::GetBuf() const -> Word *
{
    return &((*data_)[0]);
}

//------------------------------------------------------------------------------
template< typename T_Alloc >
inline auto CRawDataBuf< T_Alloc >::GetData() const -> Data const &
{
    return *data_.get();
}

//------------------------------------------------------------------------------
template< typename T_Alloc >
inline void CRawDataBuf< T_Alloc >::Reset()
{
    data_.reset( new Data );
}

//------------------------------------------------------------------------------
template< typename T_Alloc >
inline void CRawDataBuf< T_Alloc >::clear()
{
    data_->clear();
}

//------------------------------------------------------------------------------
template< typename T_Alloc >
template< typename T_Word, EDataUnit T_UNIT >
inline CDataBuf< T_Alloc >::CDataBuf( 
        CFixedDataBuf< T_Word, T_UNIT, T_Alloc > const & other )
    : word_( T_Word::WORD ), 
      unit_( T_UNIT ),
      rawbuf_( other.GetRawBuf() ),
      sz_( other.size() )
{
}

//------------------------------------------------------------------------------
template< typename T_Alloc >
template< typename T_Word, EDataUnit T_UNIT, typename T_OtherAlloc >
inline CDataBuf< T_Alloc >::CDataBuf( 
        CFixedDataBuf< T_Word, T_UNIT, T_OtherAlloc > const & other,
        T_OtherAlloc const & other_alloc )
    : word_( T_Word::WORD ),
      unit_( T_UNIT ),
      rawbuf_( other.GetRawBuf(), other_alloc ),
      sz_( other.size() )
{
}

//------------------------------------------------------------------------------
template< typename T_Word, EDataUnit T_UNIT, typename T_Alloc >
inline CFixedDataBuf< T_Word, T_UNIT, T_Alloc >::CFixedDataBuf( 
        Alloc const & alloc )
    : rawbuf_( 0, alloc )
{
}

//------------------------------------------------------------------------------
template< typename T_Word, EDataUnit T_UNIT, typename T_Alloc >
void CFixedDataBuf< T_Word, T_UNIT, T_Alloc >::Reset()
{
    rawbuf_.Reset();
    sz_ = 0;
}

//------------------------------------------------------------------------------
template< typename T_Word, EDataUnit T_UNIT, typename T_Alloc >
void CFixedDataBuf< T_Word, T_UNIT, T_Alloc >::clear()
{
    rawbuf_.clear();
    sz_ = 0;
}

//------------------------------------------------------------------------------
template< typename T_Word, EDataUnit T_UNIT, typename T_Alloc >
inline void CFixedDataBuf< T_Word, T_UNIT, T_Alloc >::resize( size_t n_units )
{
    rawbuf_.resize( (n_units + RAW_UNIT_RATIO - 1)/RAW_UNIT_RATIO );
    sz_ = n_units;
}

//------------------------------------------------------------------------------
template< typename T_Word, EDataUnit T_UNIT, typename T_Alloc >
inline void CFixedDataBuf< T_Word, T_UNIT, T_Alloc >::resize( 
        size_t n_units, typename Word::T v )
{
    rawbuf_.resize( (n_units + RAW_UNIT_RATIO - 1)/RAW_UNIT_RATIO, v );
    sz_ = n_units;
}

//------------------------------------------------------------------------------
template< typename T_Word, EDataUnit T_UNIT, typename T_Alloc >
inline auto CFixedDataBuf< T_Word, T_UNIT, T_Alloc >::GetBuf() const 
    -> typename Word::T *
{
    return (typename Word::T *)rawbuf_.GetBuf();
}

//------------------------------------------------------------------------------
template< typename T_Word, EDataUnit T_UNIT, typename T_Alloc >
inline auto CFixedDataBuf< T_Word, T_UNIT, T_Alloc >::GetRawBuf() const 
    -> RawBuf
{
    return rawbuf_;
}

//------------------------------------------------------------------------------
template< typename T_Word, EDataUnit T_UNIT, typename T_Alloc >
inline size_t CFixedDataBuf< T_Word, T_UNIT, T_Alloc >::size() const
{
    return sz_;
}

//------------------------------------------------------------------------------
template< typename T_Word, EDataUnit T_UNIT, typename T_Alloc >
inline size_t CFixedDataBuf< T_Word, T_UNIT, T_Alloc >::GetPadding() const
{
    size_t rem( sz_%UNIT_RATIO );
    return rem == 0 ? 0 : UNIT_RATIO - rem;
}

//------------------------------------------------------------------------------
template< typename T_Word, EDataUnit T_UNIT, typename T_Alloc >
inline void CFixedDataBuf< T_Word, T_UNIT, T_Alloc >::Append( 
        typename Word::T const * words, size_t n_units )
{
    size_t old_sz( sz_ );
    resize( sz_ + n_units );

    if( old_sz%UNIT_RATIO == 0 )
    {
        CopyData< T_UNIT >( GetBuf() + old_sz/UNIT_RATIO, words, n_units );
    }
    else
    {
        CopyData< T_UNIT >( GetBuf(), words, old_sz, n_units );
    }
}

TOOLS_NS_END

#endif

