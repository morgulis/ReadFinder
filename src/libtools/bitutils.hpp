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
/** \file libtools/bitutils.hpp
    \brief Bit manipulation utilities.
*/

#ifndef LIBTOOLS_BITUTILS_HPP
#define LIBTOOLS_BITUTILS_HPP

#include <cstring>

#include <libtools/exception.hpp>
#include <libtools/defs.hpp>

TOOLS_NS_BEGIN

//------------------------------------------------------------------------------
/** Accepted word type names. */
typedef size_t EDataWord;

static EDataWord const eW8  = 1;   ///< One byte.
static EDataWord const eW16 = 2;   ///< Two bytes.
static EDataWord const eW32 = 4;   ///< Four bytes.
static EDataWord const eW64 = 8;   ///< Eight bytes.

//------------------------------------------------------------------------------
/** Accepted unit (sub-word) sizes. */
typedef size_t EDataUnit;

static EDataUnit const eU1 = 1;    ///< 1 bit units.
static EDataUnit const eU2 = 2;    ///< 2 bit units.
static EDataUnit const eU4 = 4;    ///< 4 bit units.
static EDataUnit const eU8 = 8;    ///< 8 bit units.

//------------------------------------------------------------------------------
template< EDataWord T_SZ > struct TypeBySize {};

template<> struct TypeBySize< eW32 >
{
    typedef uint32_t UNSIGNED_INT_TYPE;
};

//------------------------------------------------------------------------------
/** Traits of each word type. */
template< typename T_Word > struct CWord
{
    /// The corresponding unsigned integer type.
    typedef T_Word T;

    /// Word type name and the number of bytes.
    static EDataWord const WORD = sizeof( T );

    /// Number of bits in a word.
    static size_t const N_BITS = WORD*BYTE_BITS;
};

//------------------------------------------------------------------------------
/** Word value with all bits set to 1. */
template< typename T_Word > struct CFullMask
{
    static T_Word const V = ~((T_Word)0ULL);
};

//------------------------------------------------------------------------------
/** Sub-word unit traits. */
template< EDataUnit T_UNIT > struct CUnit
{
    static_assert( T_UNIT == eU1 ||
                   T_UNIT == eU2 ||
                   T_UNIT == eU4 ||
                   T_UNIT == eU8,
                   "" );

    static EDataUnit const UNIT = T_UNIT;   ///< Unit name.

    /// Number of bits in a unit.
    static size_t const N_BITS = UNIT;

    /// Number of units in a word.
    template< typename T_Word > struct WordUnits
    {
        static size_t const N_UNITS = T_Word::N_BITS/N_BITS;
    };
};

//------------------------------------------------------------------------------
/** Get a word value with specified bits set to 1.

    \tparam T_Word  Unsigned integer word type.

    \param [in] begin   The first set bit.
    \param [in] end     One past the last set bit.

    \returns The corresponding mask value of type T_Word.
*/
template< typename T_Word >
inline T_Word MaskBits( size_t begin, size_t end )
{
    typedef CWord< T_Word > Word;

    assert( begin < end );
    assert( end <= Word::N_BITS );

    size_t l( end - begin );
    T_Word res( CFullMask< T_Word >::V );
    
    if( l < Word::N_BITS )
    {
        res >>= (Word::N_BITS - l);
        res <<= begin;
    }

    return res;
}

//------------------------------------------------------------------------------
/** Select specified bits from a word.

    All other bits are set to 0.

    \tparam T_Word  Unsigned integer word type.

    \param [in] src     Source value.
    \param [in] begin   The first bit.
    \param [in] end     One past the last bit.

    \returns Value of type T_Word with bits <tt>[begin, end)</tt> equal
             to the corresponding bits in \c src; all other bits set to 0.
*/
template< typename T_Word >
inline T_Word SelectBits( T_Word src, size_t begin, size_t end )
{
    return src&MaskBits< T_Word >( begin, end );
}

//------------------------------------------------------------------------------
/** Set the specified bits in a destination word from the source word.

    The rest of the bits in the destination are left unchanged.
    Source must have all bits outside of <tt>[begin,end)</tt> cleared.

    \tparam T_Word  Unsigned integer word type.

    \param [in,out] dst     Destination value.
    \param [in]     src     Source value.
    \param [in]     begin   The first bit.
    \param [in]     end     One past the last bit.
*/
template< typename T_Word >
inline void SetBits( T_Word & dst, T_Word src, size_t begin, size_t end )
{
    dst &= ~MaskBits< T_Word >( begin, end );
    dst |= src;
}

//------------------------------------------------------------------------------
/** Set the specified bits in a destination word from the corresponding bits
    in the source word.

    The rest of the bits in the destination are left unchanged.
    Source bits outside of <tt>[begin, end)</tt> are ignored and can have any
    value.

    \tparam T_Word  Unsigned integer word type.

    \param [in,out] dst     Destination value.
    \param [in]     src     Source value.
    \param [in]     begin   The first bit.
    \param [in]     end     One past the last bit.
*/
template< typename T_Word >
inline void CopyBits( T_Word & dst, T_Word src, size_t begin, size_t end )
{
    SetBits( dst, SelectBits( src, begin, end ), begin, end );
}

//------------------------------------------------------------------------------
/** Set a particular bit in a word.

    \tparam T_Word  Unsigned integer word type.

    \param [in,out] dst     Destination value.
    \param [in]     bit     Bit to set.
*/
template< typename T_Word > inline void SetBit( T_Word & dst, size_t bit )
{
    assert( bit < CWord< T_Word >::N_BITS );
    dst |= (((T_Word)1)<<bit);
}

//------------------------------------------------------------------------------
/** Clear a particular bit in a word.

    \tparam T_Word  Unsigned integer word type.

    \param [in,out] dst     Destination value.
    \param [in]     bit     Bit to clear.
*/
template< typename T_Word > inline void ClearBit( T_Word & dst, size_t bit )
{
    assert( bit < CWord< T_Word >::N_BITS );
    dst &= (~(((T_Word)1)<<bit));
}

//------------------------------------------------------------------------------
/** Get the value of a particular bit in a word.

    \tparam T_Word  Unsigned integer word type.

    \param [in] dst Source value.
    \param [in] bit Bit to read.
*/
template< typename T_Word > inline T_Word GetBit( T_Word src, size_t bit )
{
    assert( bit < CWord< T_Word >::N_BITS );
    return ((src>>bit)&1);
}

//------------------------------------------------------------------------------
/** Get a bit field and return its value as an unsigned integer.

    The value of the bitfield is returned in the lower significant
    <tt>end - begin</tt> bits of the result, other bits being cleared.

    \tparam T_Word  Unsigned integer word type.

    \param [in]     src     Source value.
    \param [in]     begin   The first bit.
    \param [in]     end     One past the last bit.

    \return The value of the specified bit field.
*/
template< typename T_Word >
inline T_Word GetField( T_Word src, size_t begin, size_t end )
{
    return SelectBits( src, begin, end )>>begin;
}

//------------------------------------------------------------------------------
/** Set a bit field.

    The input is the value of the bit field that is set into the 
    specified bits of the destination word.

    \tparam T_Word  Unsigned integer word type.

    \param [in,out] dst     Destination value.
    \param [in]     begin   The first bit.
    \param [in]     end     One past the last bit.
*/
template< typename T_Word >
inline void SetField( T_Word & dst, T_Word src, size_t begin, size_t end )
{
    SetBits( dst, (T_Word)(src<<begin), begin, end );
}

//------------------------------------------------------------------------------
/** Reverse 1-bit, 2-bit, or 4-bit units within a byte value.
   
    \tparam T_UNIT  Type of unit (number of bits).

    \param [in,out] byte    Byte in which units are to be reversed.
*/
template< EDataUnit T_UNIT > void ReverseByteUnits( uint8_t & byte );

extern uint8_t const REVERSE_BYTE_UNITS_TABLE_1[];
extern uint8_t const REVERSE_BYTE_UNITS_TABLE_2[];

template<> 
inline void ReverseByteUnits< 1 >( uint8_t & byte )
{
    byte = REVERSE_BYTE_UNITS_TABLE_1[byte];
}

template<> 
inline void ReverseByteUnits< 2 >( uint8_t & byte )
{
    byte = REVERSE_BYTE_UNITS_TABLE_2[byte];
}

template<> 
inline void ReverseByteUnits< 4 >( uint8_t & byte )
{
    byte = (byte>>4) + ((byte&0xF)<<4);
}

//------------------------------------------------------------------------------
/** Reverse bytes in a word.

    \tparam T_Word  Unsigned integer word type.

    \param [in.out] word    The value in which bytes should be reversed.
*/
template< typename T_Word >
void ReverseBytes( T_Word & word )
{
    uint8_t * bytes( (uint8_t *)&word );

    for( size_t i( 0 ); i < sizeof( T_Word )/2; ++i )
    {
        std::swap( bytes[i], bytes[sizeof( T_Word ) - i - 1] );
    }
}

//------------------------------------------------------------------------------
template< EDataUnit T_UNIT, typename T_Word > struct CWordReverser
{
    static_assert( T_UNIT == eU1 ||
                   T_UNIT == eU2 ||
                   T_UNIT == eU4,
                   "" );

    static void Do( T_Word & word )
    {
        ReverseBytes( word );
        uint8_t * bytes( (uint8_t *)&word );

        for( size_t i( 0 ); i < sizeof( T_Word ); ++i )
        {
            ReverseByteUnits< T_UNIT >( bytes[i] );
        }
    }
};

template< typename T_Word > struct CWordReverser< 8, T_Word >
{
    static void Do( T_Word & word )
    {
        ReverseBytes( word );
    }
};

/** Reverse units in a word.

    \tparam T_UNIT  Unit type (number of bits).
    \tparam T_Word  Word type (unsigned integer).

    \param [in,out] word    Value to reverse.
*/
template< EDataUnit T_UNIT, typename T_Word > void ReverseWord( T_Word & word )
{
    CWordReverser< T_UNIT, T_Word >::Do( word );
}

//------------------------------------------------------------------------------
/** Reverse several units within a word.

    \tparam T_UNIT  Unit type (number of bits).
    \tparam T_Word  Word type (unsigned integer).

    \param [in,out] word    Value to reverse.
    \param [in]     begin   Index of the first unit.
    \param [out]    end     Index of the one past the last unit.
*/
template< EDataUnit T_UNIT, typename T_Word >
void ReverseField( T_Word & w, size_t begin, size_t end )
{
    static_assert( T_UNIT == eU1 ||
                   T_UNIT == eU2 ||
                   T_UNIT == eU4 ||
                   T_UNIT == eU8,
                   "" );

    typedef CWord< T_Word > Word;

    size_t bb( begin*T_UNIT ),
           ee( end*T_UNIT );
    T_Word ww( w );
    ww >>= bb;
    ReverseWord< T_UNIT >( ww );
    ww >>= (Word::N_BITS - ee);
    T_Word mask( MaskBits< T_Word >( bb, ee ) );
    ww &= mask;
    w &= ~mask;
    w |= ww;
}

//------------------------------------------------------------------------------
/** Copy specified number of units of data from a source buffer to a 
    destination buffer.

    Data is copied from source offset 0 to destination offset 0. 
    All units in destination beyond \c n_units are left unchanged.

    \tparam T_UNIT  Data unit name.
    \tparam T_Word  Unsigned integer word type.

    \param [in,out] dst     Destination buffer.
    \param [in]     src     Source buffer.
    \param [in]     n_units Number of units to copy.
*/
template< EDataUnit T_UNIT, typename T_Word >
inline void CopyData( T_Word * dst, T_Word const * src, size_t n_units )
{
    static_assert( T_UNIT == eU1 ||
                   T_UNIT == eU2 ||
                   T_UNIT == eU4 ||
                   T_UNIT == eU8,
                   "" );

    typedef CWord< T_Word > Word;
    static size_t const UNIT_RATIO = Word::N_BITS/T_UNIT;

    size_t n_words( n_units/UNIT_RATIO );
    memmove( dst, src, n_words*sizeof( T_Word ) );
    size_t tail_units( n_units%UNIT_RATIO );

    if( tail_units > 0 )
    {
        CopyBits( dst[n_words], src[n_words], 0, tail_units*T_UNIT );
    }
}

//------------------------------------------------------------------------------
/** Get a word from a word stream at the specified bit position.

    \tparam T_Word  Unsigned integer word type.

    \param [in] src Source buffer.
    \param [in] off Bit offset from the start of the buffer.

    \return Word that starts at the specified bit offset.
*/
template< typename T_Word >
inline T_Word GetWord( T_Word const * src, ssize_t off )
{
    static size_t const WBITS = CWord< T_Word >::N_BITS;

    if( off%WBITS == 0 ) 
    {
        return src[off/(ssize_t)WBITS];
    }

    src += off/(ssize_t)WBITS;

    if( off < 0 )
    {
        --src;
    }

    off %= WBITS;
    return GetField( *src, off, WBITS ) + 
           (SelectBits( *(src + 1), 0, off )<<(WBITS - off));
}

//------------------------------------------------------------------------------
/** Set a word in a word stream at specified bit position.

    \tparam T_Word  Unsigned integer word type.

    \param [in,out] dst Destination buffer.
    \param [in]     src Source word.
    \param [in]     off Bit offset from the start of the destination buffer.
*/
template< typename T_Word >
inline void SetWord( T_Word * dst, T_Word src, ssize_t off )
{
    static size_t const WBITS = CWord< T_Word >::N_BITS;

    dst += off/(ssize_t)WBITS;

    if( off%WBITS == 0 )
    {
        *dst = src;
        return;
    }

    if( off < 0 )
    {
        --dst;
    }

    off %= WBITS;
    SetField( *dst++, SelectBits( src, 0, WBITS - off ), off, WBITS );
    SetBits( *dst, GetField( src, WBITS - off, WBITS ), 0, off );
}

//------------------------------------------------------------------------------
/** Copy a specified number of data units from source buffer into the
    destination buffer.

    Buffers must not overlap. The bits ouside the of the copied-to area
    of the destination are unaffected.

    \tparam T_UNIT  Unit size in bits.
    \tparam T_Word  Unsigned integer word type.

    \param [in,out] dst     Destination buffer.
    \param [in]     src     Source buffer.
    \param [in]     dst_off Offset (in units) into destination.
    \param [in]     src_off Offset (in units) into source.
    \param [in]     n_units Number of units to copy.
*/
template< EDataUnit T_UNIT, typename T_Word >
inline void CopyData(
        T_Word * dst, T_Word const * src, 
        ssize_t dst_off, ssize_t src_off, size_t n_units )
{
    typedef CWord< T_Word > Word;
    typedef CUnit< T_UNIT > Unit;
    static size_t const WUNITS = Unit::template WordUnits< Word >::N_UNITS;
    static size_t const WBITS = Word::N_BITS;

    if( dst_off < 0 )
    {
        dst += (dst_off - (ssize_t)WUNITS + 1)/(ssize_t)WUNITS;
    }
    else
    {
        dst += dst_off/(ssize_t)WUNITS;
    }

    if( src_off < 0 )
    {
        src += (src_off - (ssize_t)WUNITS + 1)/(ssize_t)WUNITS;
    }
    else
    {
        src += src_off/(ssize_t)WUNITS;
    }

    dst_off %= WUNITS;
    src_off %= WUNITS;

    if( src_off == 0 && dst_off == 0 )
    {
        CopyData< T_UNIT >( dst, src, n_units );
        return;
    }

    size_t n_words( n_units/WUNITS );
    ssize_t si( src_off*T_UNIT ),
            di( dst_off*T_UNIT );

    for( size_t i( 0 ); i < n_words; ++i, ++dst, ++src )
    {
        SetWord( dst, GetWord( src, si ), di );
    }

    size_t n_bits( (n_units%WUNITS)*T_UNIT ),
           low_bits( std::min( n_bits, WBITS - di ) );
    T_Word w( GetWord( src, si ) );

    if( n_bits > 0 )
    {
        SetField( *dst, SelectBits( w, 0, low_bits ), di, di + low_bits );

        if( low_bits < n_bits )
        {
            ++dst;
            SetBits( *dst, GetField( w, low_bits, n_bits ), 
                    0, n_bits - low_bits );
        }
    }
}

//------------------------------------------------------------------------------
/** Clear a range of bits in a buffer.

    \tparam T_Word  Unsigned integer word type.

    \param [in,out] buf     Buffer to work on.
    \param [in]     off     Offset of the first bit to clear.
    \param [in]     n_bits  Number of bits to clear.
*/
template< typename T_Word >
void ClearBits( T_Word * buf, ssize_t off, size_t n_bits )
{
    static size_t const WBITS = CWord< T_Word >::N_BITS;

    if( off < 0 )
    {
        buf += (off - (ssize_t)WBITS + 1)/(ssize_t)WBITS;
    }
    else
    {
        buf += off/WBITS;
    }

    off %= WBITS;

    if( off + n_bits <= WBITS )
    {
        *buf &= ~MaskBits< T_Word >( off, off + n_bits );
        return;
    }

    if( off > 0 )
    {
        *buf++ &= ~MaskBits< T_Word >( off, WBITS );
        n_bits -= (WBITS - off);
    }

    assert( n_bits > 0 );

    for( ; n_bits >= WBITS; n_bits -= WBITS )
    {
        *buf++ = 0;
    }

    if( n_bits > 0 )
    {
        *buf &= ~MaskBits< T_Word >( 0, n_bits );
    }
}

//------------------------------------------------------------------------------
/** Shift bits in a buffer to the left (towards higher significance).

    \tparam T_Word  Unsigned integer word type.

    \param [in,out] buf     Buffer to work on.
    \param [in]     off     Offset of the first bit to shift.
    \param [in]     n_bits  Number of bits to shift.
    \param [in]     shift   Shift value.
*/
template< typename T_Word >
void ShiftBitsLeft( T_Word * buf, ssize_t off, size_t n_bits, size_t shift )
{
    static ssize_t const WBITS = (ssize_t)CWord< T_Word >::N_BITS;

    if( shift == 0 )
    {
        return;
    }

    if( shift >= n_bits )
    {
        ClearBits( buf, off, n_bits );
        return;
    }

    ssize_t off_w( off + n_bits - shift - WBITS );

    for( ; off_w >= off; off_w -= WBITS )
    {
        SetWord( buf, GetWord( buf, off_w ), off_w + shift );
    }

    if( off_w + WBITS > off )
    {
        ssize_t d( off_w + WBITS - off );
        assert( d > 0 && d < WBITS );
        T_Word * b( buf + (off < 0 ? (off - WBITS + 1)/WBITS : off/WBITS ) );
        ssize_t o( off%(size_t)WBITS );
        T_Word w;

        if( o + d < WBITS )
        {
            w = GetField( *b, o, o + d );
        }
        else
        {
            w = (GetWord( b, o )&MaskBits< T_Word >( 0, d ));
        }

        CopyData< 1 >( buf, &w, off + shift, 0, d );
    }

    ClearBits( buf, off, shift );
}

//------------------------------------------------------------------------------
/** Shift bits in a buffer to the right (towards lower significance).

    \tparam T_Word  Unsigned integer word type.

    \param [in,out] buf     Buffer to work on.
    \param [in]     off     Offset of the first bit to shift.
    \param [in]     n_bits  Number of bits to shift.
    \param [in]     shift   Shift value.
*/
template< typename T_Word >
void ShiftBitsRight( T_Word * buf, ssize_t off, size_t n_bits, size_t shift )
{
    static ssize_t const WBITS = (ssize_t)CWord< T_Word >::N_BITS;

    if( shift == 0 ) return;

    if( shift >= n_bits )
    {
        ClearBits( buf, off, n_bits );
        return;
    }

    ssize_t off_e( off + n_bits ),
            off_s( off + shift );

    for( ; off_s + WBITS <= off_e; off_s += WBITS )
    {
        SetWord( buf, GetWord( buf, off_s ), off_s - shift );
    }

    if( off_s < off_e )
    {
        ssize_t d( off_e - off_s );
        assert( d > 0 && d < WBITS );
        T_Word * b( 
                buf + (off_s < 0 ? (off_s - WBITS + 1)/WBITS : off_s/WBITS ) );
        ssize_t o( off_s%(size_t)WBITS );
        T_Word w;

        if( o + d < WBITS )
        {
            w = GetField( *b, o, o + d );
        }
        else
        {
            w = (GetWord( b, o )&MaskBits< T_Word >( 0, d ));
        }

        CopyData< 1 >( buf, &w, off_s - shift, 0, d );
    }

    ClearBits( buf, off_e - shift, shift );
}

//------------------------------------------------------------------------------
/** Copy a specified number of data units from source buffer into the
    destination buffer.

    Buffers must not overlap. The bits ouside the of the copied-to area
    of the destination are unaffected. The data is copied from the start
    of the source buffer.

    \tparam T_UNIT  Unit size in bits.
    \tparam T_Word  Unsigned integer word type.

    \param [in,out] dst     Destination buffer.
    \param [in]     src     Source buffer.
    \param [in]     dst_off Offset (in units) into destination.
    \param [in]     n_units Number of units to copy.
*/
template< EDataUnit T_UNIT, typename T_Word >
inline void CopyData( 
        T_Word * dst, T_Word const * src, ssize_t dst_off, size_t n_units )
{
    CopyData< T_UNIT >( dst, src, dst_off, 0, n_units );
}

//------------------------------------------------------------------------------
/** Reverse units in a buffer for a whole number of words.

    \tparam T_UNIT  Unit size in bits.
    \tparam T_Word  Unsigned integer word type.

    \param [in,out] buf     Buffer to reverse.
    \param [in]     n_words Number of words.
*/
template< EDataUnit T_UNIT, typename T_Word >
inline void ReverseUnits( T_Word * buf, size_t n_words )
{
    for( size_t i( 0 ), ie( n_words/2 ); i < ie; ++i )
    {
        ReverseWord< T_UNIT >( buf[i] );
        ReverseWord< T_UNIT >( buf[n_words - i - 1] );
        std::swap( buf[i], buf[n_words - i - 1] );
    }

    if( n_words%2 != 0 )
    {
        ReverseWord< T_UNIT >( buf[n_words/2] );
    }
}

//------------------------------------------------------------------------------
/** Reverse units in a buffer for any number of units.

    \tparam T_UNIT  Unit size in bits.
    \tparam T_Word  Unsigned integer word type.

    \param [in,out] buf     Start of the buffer.
    \param [in]     off     Offset (in units) of the start of the region being
                            reversed
    \param [in]     n_units Number of units to reverse.
*/
template< EDataUnit T_UNIT, typename T_Word >
inline void ReverseUnits( T_Word * buf, ssize_t off, size_t n_units )
{
    typedef CWord< T_Word > Word;
    typedef CUnit< T_UNIT > Unit;
    static size_t const WUNITS = Unit::template WordUnits< Word >::N_UNITS;

    if( off < 0 )
    {
        buf += (off - (ssize_t)WUNITS + 1)/(ssize_t)WUNITS;
    }
    else
    {
        buf += off/WUNITS;
    }

    off %= WUNITS;

    if( off == 0 && n_units%WUNITS == 0 )
    {
        ReverseUnits< T_UNIT >( buf, n_units/WUNITS );
        return;
    }

    size_t n_bits( n_units*T_UNIT );

    if( off == 0 )
    {
        size_t n_words( n_units/WUNITS ),
               shift_units( n_units%WUNITS );
        size_t shift_bits( T_UNIT*shift_units );
        T_Word last( buf[n_words] );
        ReverseUnits< T_UNIT >( buf, n_words );
        ShiftBitsLeft( buf, 0, n_bits, shift_bits );
        ReverseField< T_UNIT >( last, 0, shift_units );
        SetBits( *buf, SelectBits( last, 0, shift_bits ), 0, shift_bits );
        return;
    }

    if( off + n_units < WUNITS )
    {
        ReverseField< T_UNIT >( *buf, off, off + n_units );
        return;
    }

    size_t n_units_adj( n_units - (WUNITS - off) );
    size_t n_words( n_units_adj/WUNITS );
    ReverseUnits< T_UNIT >( buf + 1, n_words );
    size_t right_flank( WUNITS - off ),
           left_flank( n_units_adj%WUNITS );
    T_Word first( *buf ),
           last( left_flank > 0 ? *(buf + n_words + 1) : 0 );
    ReverseField< T_UNIT >( first, off, WUNITS );

    if( left_flank > 0 )
    {
        ReverseField< T_UNIT >( last, 0, left_flank );
    }

    if( left_flank > right_flank )
    {
        ShiftBitsLeft( buf, off*T_UNIT, n_bits, 
                       (left_flank - right_flank)*T_UNIT );
    }
    else if( left_flank < right_flank )
    {
        ShiftBitsRight( buf, off*T_UNIT, n_bits,
                        (right_flank - left_flank)*T_UNIT );
    }

    if( left_flank > 0 )
    {
        CopyData< T_UNIT >( buf, &last, off, 0, left_flank );
    }

    CopyData< T_UNIT >( 
            buf, &first, off + left_flank + n_words*WUNITS, off, right_flank );
}

//------------------------------------------------------------------------------
/** Get the number of first clear bits in a word. */
template< typename T_Word > inline size_t GetNFirstClear( T_Word w )
{
    size_t res( __builtin_ffsll( w ) );
    return res == 0 ? CWord< T_Word >::N_BITS : res - 1;
}

//------------------------------------------------------------------------------
/** Get the number of first set bits in a word. */
template< typename T_Word > inline size_t GetNFirstSet( T_Word w )
{
    return GetNFirstClear( ~w );
}

//------------------------------------------------------------------------------
/** Get the number of set bits in a word. */
template< typename T_Word > inline size_t GetNSet( T_Word w )
{
    return __builtin_popcountll( (uint64_t )w );
}

//------------------------------------------------------------------------------
/** Get the number of clear bits in a word. */
template< typename T_Word > inline size_t GetNClear( T_Word w )
{
    return CWord< T_Word >::N_BITS - GetNSet( w );
}

//------------------------------------------------------------------------------
template< EDataUnit T_UNIT, typename T_Word > struct CMatchMask;

template< typename T_Word > struct CMatchMask< 1, T_Word >
{
    static T_Word const V = CFullMask< T_Word >::V;
};

template<> struct CMatchMask< 2, uint8_t >  { static uint8_t  const V = 0x55U; };
template<> struct CMatchMask< 2, uint16_t > { static uint16_t const V = 0x5555U; };
template<> struct CMatchMask< 2, uint32_t > { static uint32_t const V = 0x55555555UL; };
template<> struct CMatchMask< 2, uint64_t > { static uint64_t const V = 0x5555555555555555ULL; };

template<> struct CMatchMask< 4, uint8_t >  { static uint8_t  const V = 0x11U; };
template<> struct CMatchMask< 4, uint16_t > { static uint16_t const V = 0x1111U; };
template<> struct CMatchMask< 4, uint32_t > { static uint32_t const V = 0x11111111UL; };
template<> struct CMatchMask< 4, uint64_t > { static uint64_t const V = 0x1111111111111111ULL; };

template<> struct CMatchMask< 8, uint8_t >  { static uint8_t  const V = 0x1U; };
template<> struct CMatchMask< 8, uint16_t > { static uint16_t const V = 0x101U; };
template<> struct CMatchMask< 8, uint32_t > { static uint32_t const V = 0x1010101UL; };
template<> struct CMatchMask< 8, uint64_t > { static uint64_t const V = 0x101010101010101ULL; };

//------------------------------------------------------------------------------
/** Get the length of the initial exact match (in units) for two words. 

    Positions masked by \c mask are counted as mismatches.
*/
template< EDataUnit T_UNIT, typename T_Word >
inline size_t GetExactMatch( T_Word w1, T_Word w2, T_Word mask = 0 )
{
    T_Word w( (w1^w2)|mask ), w_shift( w );

    for( size_t i( 1 ); i < T_UNIT; ++i )
    {
        w_shift >>= 1;
        w |= w_shift;
    }

    w &= CMatchMask< T_UNIT, T_Word >::V;
    return GetNFirstClear( w )/T_UNIT;
}

TOOLS_NS_END

#endif

