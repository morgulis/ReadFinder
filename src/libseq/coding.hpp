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
/** \file libseq/coding.hpp
    \brief Definitions related to sequence coding.
*/

#ifndef LIBSEQ_CODING_HPP
#define LIBSEQ_CODING_HPP

#include <libseq/defs.hpp>

#include <libtools/bitutils.hpp>

SEQ_NS_BEGIN

using namespace TOOLS_NS;

//------------------------------------------------------------------------------
/** Supported nucleotide sequence encodings. */
enum ECoding : unsigned int
{
    eNCBI2NA,   ///< 2-bits per base encoding (no ambiguities).
    eIUPACNA,   ///< IUPAC nucleotide encoding (via ASCII text).
};

//------------------------------------------------------------------------------
/** Static computation of bases per byte

    \tparam T_LBITS Bits per base.
*/
template< size_t T_LBITS > struct CBL
{
    /// Bases per byte
    static TSeqLen const V = (TSeqLen)BYTE_BITS/T_LBITS;
};

//------------------------------------------------------------------------------
/** Common (non-specialized) definitions for base coding functions. 

    \tparam T_CODE  Sequence encoding.
*/
template< ECoding T_CODE > struct CCodeBaseCommon
{
    static ECoding const CODING = T_CODE;   ///< Encoding identifier.
};

//------------------------------------------------------------------------------
/** Base (specialized) encoding functions. 

    \tparam T_CODE  Sequence encoding.
*/
template< ECoding T_CODE > struct CCodeBase;

/** Specialization for NCBI2NA encoding. */
template<> struct CCodeBase< eNCBI2NA > : public CCodeBaseCommon< eNCBI2NA >
{
    static size_t const LBITS = (size_t)2;  ///< Number of bits per base.
};

/** Specialization for IUPACNA encoding. */
template<> struct CCodeBase< eIUPACNA > : public CCodeBaseCommon< eIUPACNA >
{
    static size_t const LBITS = (size_t)8;  ///< Number of bits per base.
};

//------------------------------------------------------------------------------
/** Common (non-specialized) encoding related definitions. 

    \tparam T_CODE  Sequence encoding.
*/
template< ECoding T_CODE > struct CCodeCommon : public CCodeBase< T_CODE >
{
private:

    typedef CCodeBase< T_CODE > TBase;

public:

    /// Data unit used to hold a base in this encoding.
    static EDataUnit const UNIT = TBase::LBITS;

    /// Value representing a non-masked base in mask data.
    static uint64_t const CLEAR_BASE = 0;

    /// Value representing a masked base in mask data.
    static uint64_t const MASK_BASE = (1<<TBase::LBITS) - 1;

    /// Unit trait types corresponding to UNIT.
    typedef CUnit< UNIT > Unit;

    static bool const AMBIG[];      ///< Ambiguity table.
    static uint8_t const COMPL[];   ///< Complement base table.

    template< typename T_Word > struct WordBases
    {
        static TSeqLen const V = Unit::template WordUnits< T_Word >::N_UNITS;
    };
};

//------------------------------------------------------------------------------
/** Specialized encoding related definitions.

    \tparam T_CODE  Sequence encoding.
*/
template< ECoding T_CODE > struct CCode;

//------------------------------------------------------------------------------
/// Specializations.
/// @{
template<> struct CCode< eNCBI2NA > : public CCodeCommon< eNCBI2NA > {};
template<> struct CCode< eIUPACNA > : public CCodeCommon< eIUPACNA > {};
/// @}

//------------------------------------------------------------------------------
/** Get a single letter from coded word.

    \tparam T_CODE  Encoding used.
    \tparam T_Word  Word type.

    \param [in] src Source word.
    \param [in] off Offset into the sequence coded into the word.

    \return The selected base.
*/
template< ECoding T_CODE, typename T_Word >
inline T_Word GetLetter( T_Word src, TSeqOff off )
{
    static size_t const LB = CCode< T_CODE >::LBITS;
    off *= LB;
    return GetField( src, off, off + LB );
}

//------------------------------------------------------------------------------
/** Get a single letter from a sequence buffer.

    \tparam T_CODE  Encoding used.
    \tparam T_Word  Word type.

    \param [in] src Source sequence data.
    \param [in] off Offset into the sequence.

    \return The selected base.
*/
template< ECoding T_CODE, typename T_Word >
inline T_Word GetLetter( T_Word const * src, TSeqOff off )
{
    static TSeqLen const WL = 
        CCode< T_CODE >::template WordBases< CWord< T_Word > >::V;
    return GetLetter< T_CODE >( 
            src[off < 0 ? (off - WL + 1)/WL : off/WL], 
            off%(size_t)WL );
}

template< ECoding T_CODE, typename T_Word >
inline T_Word GetLetter( T_Word * src, TSeqOff off )
{
    return GetLetter< T_CODE >( (T_Word const *)src, off );
}

//------------------------------------------------------------------------------
/** Set a single letter of a coded word.

    \tparam T_CODE  Encoding used.
    \tparam T_Word  Word type.

    \param [in,out] dst     Destination word.
    \param [in]     off     Offset into the sequence coded into the word.
    \param [in]     letter  The new value of the letter.
*/
template< ECoding T_CODE, typename T_Word >
inline void SetLetter( T_Word & dst, TSeqLen off, T_Word letter )
{
    static size_t const LB = CCode< T_CODE >::LBITS;
    off *= LB;
    SetField( dst, letter, off, off + LB );
}

//------------------------------------------------------------------------------
/** Set a single letter in sequence data.

    \tparam T_CODE  Encoding used.
    \tparam T_Word  Word type.

    \param [in,out] dst     Destination sequence.
    \param [in]     off     Offset into the sequence coded into the word.
    \param [in]     letter  The new value of the letter.
*/
template< ECoding T_CODE, typename T_Word >
inline void SetLetter( T_Word * dst, TSeqLen off, T_Word letter )
{
    static TSeqLen const WL = 
        CCode< T_CODE >::template WordBases< CWord< T_Word > >::V;
    SetLetter< T_CODE >(
            dst[off < 0 ? (off - WL + 1)/WL : off/WL],
            off%(size_t)WL, letter );
}

//------------------------------------------------------------------------------
/** Converting individual letters between encodings.

    \tparam T_DST_CODE  Destination encoding.
    \tparam T_SRC_CODE  Source encoding.
*/
template< ECoding T_DST_CODE, ECoding T_SRC_CODE >
struct CRecoder
{
    static_assert( T_DST_CODE != T_SRC_CODE, "" );
    static uint8_t const TBL[]; ///< Conversion table.

    /** Convert a single letter from source to destination encoding. */
    static uint8_t Recode( uint8_t l ) { return TBL[l]; }
};

/** Specialization for the case when source and destination encodings are
    the same.

    In this case the conversion function is identity.
*/
template< ECoding T_CODE >
struct CRecoder< T_CODE, T_CODE >
{
    static uint8_t Recode( uint8_t l ) { return l; }
};

//------------------------------------------------------------------------------
/** Set a single letter in a destination word (and mask in the mask word)
    from a letter in a different encoding.

    The function sets/clears the mask if the source letter is/is not
    ambiguous in the source encoding.

    \tparam T_DST_CODE  Destination encoding.
    \tparam T_SRC_CODE  Source encoding.
    \tparam T_DstWord   Destination word type.
    \tparam T_SrcWord   Source word type.

    \param [in,out] dst         Destination word.
    \param [in]     doff        Destination offset (in bases).
    \param [in,out] dmask       Destination mask word.
    \param [in]     dmask_off   Destination mask offset (in bases).
    \param [in]     letter      Source letter in the source encoding.

    \return \c true if \c letter is ambiguous in the source encoding;
            \c false otherwise.
*/
template< ECoding T_DST_CODE, ECoding T_SRC_CODE, 
          typename T_DstWord, typename T_SrcWord >
bool RecodeLetter( T_DstWord & dst, TSeqLen doff,
                   T_DstWord & dmask, TSeqLen dmask_off,
                   T_SrcWord letter )
{
    SetLetter< T_DST_CODE >( 
            dst, doff, 
            (T_DstWord)CRecoder< T_DST_CODE, T_SRC_CODE >::Recode( letter ) );
    typedef CCode< T_SRC_CODE > SCode;
    typedef CCode< T_DST_CODE > DCode;
    bool ambig( SCode::AMBIG[(int)letter] );
    SetLetter< T_DST_CODE >( 
            dmask, dmask_off,
            (T_DstWord)(ambig ? DCode::MASK_BASE : DCode::CLEAR_BASE) );
    return ambig;
}

//------------------------------------------------------------------------------
/** Convert sequence data from one encoding to another.

    \tparam T_DST_CODE  Destination encoding.
    \tparam T_SRC_CODE  Source encoding.
    \tparam T_DstWord   Destination word type.
    \tparam T_SrcWord   Source word type.

    \param [in,out] dst         Desination sequence buffer.
    \param [in,out] dmask       Destination mask buffer.
    \param [in]     src         Source buffer.
    \param [in]     dst_off     Destination sequence offset.
    \param [in]     dmask_off   Destination mask offset.
    \param [in]     src_off     Source offset.
    \param [in]     len         Source sequence length.

    \return \c true if the source sequence contained an ambiguous letter;
            \c false otherwise.
*/
template< ECoding T_DST_CODE, ECoding T_SRC_CODE,
          typename T_DstWord, typename T_SrcWord >
bool Recode( T_DstWord * dst, T_DstWord * dmask, T_SrcWord const * src,
             TSeqOff dst_off, TSeqOff dmask_off, TSeqOff src_off, TSeqLen len )
{
    typedef CCode< T_DST_CODE > DCode;
    typedef CCode< T_SRC_CODE > SCode;
    typedef CWord< T_DstWord > DWord;
    typedef CWord< T_SrcWord > SWord;
    static TSeqLen const DWL = DCode::template WordBases< DWord >::V;
    static TSeqLen const SWL = SCode::template WordBases< SWord >::V;

    bool res( false );

    for( TSeqLen i( 0 ); i < len; ++i, ++dst_off, ++dmask_off, ++src_off )
    {
        dst += (dst_off < 0 ? (dst_off - DWL + 1)/DWL : dst_off/DWL);
        dmask += (dmask_off < 0 ? (dmask_off - DWL + 1)/DWL : dmask_off/DWL);
        src += (src_off < 0 ? (src_off - SWL + 1)/SWL : src_off/SWL);
        dst_off %= (size_t)DWL;
        dmask_off %= (size_t)DWL;
        src_off %= (size_t)SWL;
        bool r( RecodeLetter< T_DST_CODE, T_SRC_CODE >(
                    *dst, dst_off, *dmask, dmask_off,
                    GetLetter< T_SRC_CODE >( *src, src_off ) ) );
        res = (res || r);
    }

    return res;
}

//------------------------------------------------------------------------------
/** Reverse bases in a word.

    \tparam T_CODE  Encoding.
    \tparam T_Word  Word type.

    \param [in,out] word    Value to reverse.
*/
template< ECoding T_CODE, typename T_Word > inline void Reverse( T_Word & w )
{
    ReverseWord< CCode< T_CODE >::UNIT >( w );
}

//------------------------------------------------------------------------------
/** Reverse bases in a sequence.

    \tparam T_CODE  Encoding.
    \tparam T_Word  Word type.

    \param [in,out] buf Buffer containing the sequence data.
    \param [in]     off Offset into the buffer.
    \param [in]     len Length of the sequence to reverse.
*/
template< ECoding T_CODE, typename T_Word >
inline void Reverse( T_Word * buf, TSeqOff off, TSeqLen len )
{
    ReverseUnits< CCode< T_CODE >::UNIT >( buf, off, len );
}

//------------------------------------------------------------------------------
/** Class for computing complemented (non-reversed) sequences in a single
    data word.

    \tparam T_CODE  Sequence encoding.
    \tparam T_Word  Word type.
*/
template< ECoding T_CODE, typename T_Word > struct CComplementerBase
{
    typedef CCode< T_CODE > Code;
    typedef typename Code::Unit Unit;

    static size_t const WL = 
        Unit::template WordUnits< CWord< T_Word > >::N_UNITS;

    static void ComplementWord( T_Word & w )
    {
        T_Word s( w );

        for( size_t off( 0 ); off < WL; ++off )
        {
            SetLetter< T_CODE >( 
                    w, off, 
                    (T_Word)Code::COMPL[GetLetter< T_CODE >( s, off )] );
        }
    }
};

/** Specialization for NCBI2NA encoding, where complementing is trivial.
*/
template< typename T_Word > struct CComplementerBase< eNCBI2NA, T_Word >
{
    static void ComplementWord( T_Word & w )
    {
        w = ~w;
    }
};

/** Computing complemented sequences for multiple data words.

    \tparam T_CODE  Sequence encoding.
    \tparam T_Word  Word type.
*/
template< ECoding T_CODE, typename T_Word > struct CComplementer 
    : public CComplementerBase< T_CODE, T_Word >
{
    typedef CComplementerBase< T_CODE, T_Word > Base;

    static void ComplementWords( T_Word * buf, size_t n_words )
    {
        for( size_t i( 0 ); i < n_words; ++i, ++buf )
        {
            Base::ComplementWord( *buf );
        }
    }
};

//------------------------------------------------------------------------------
/** Complement a sequence encoded in a single word.

    \tparam T_CODE  Sequence encoding.
    \tparam T_Word  Word type.

    \param [in,out] w   Word to complement.
*/
template< ECoding T_CODE, typename T_Word >
inline void Complement( T_Word & w )
{
    CComplementer< T_CODE, T_Word >::ComplementWord( w );
}

//------------------------------------------------------------------------------
/** Complement a sequence encoded in multiple words.

    \tparam T_CODE  Sequence encoding.
    \tparam T_Word  Word type.

    \param [in,out] buf     Buffer containing the sequence buffer.
    \param [in]     n_words Number of words to complement.
*/
template< ECoding T_CODE, typename T_Word >
inline void ComplementWords( T_Word * buf, size_t n_words )
{
    CComplementer< T_CODE, T_Word >::ComplementWords( buf, n_words );
}

//------------------------------------------------------------------------------
/** Complement a sequence.

    \tparam T_CODE  Sequence encoding.
    \tparam T_Word  Word type.

    \param [in,out] buf Buffer containing the sequence buffer.
    \param [in]     off Start offset.
    \param [in]     len Length (in bases).
*/
template< ECoding T_CODE, typename T_Word >
void Complement( T_Word * buf, TSeqOff off, TSeqLen len )
{
    typedef CCode< T_CODE > Code;
    typedef CWord< T_Word > Word;
    static TSeqLen const WL = Code::template WordBases< Word >::V;

    buf += (off < 0 ? (off - WL + 1)/WL : off/WL);
    off %= (size_t)WL;

    for( ; len > 0 && off < WL; ++off, --len )
    {
        SetLetter< T_CODE >( 
                *buf, off, 
                (T_Word)Code::COMPL[GetLetter< T_CODE >( *buf, off )] );
    }

    ++buf;

    if( len >= WL )
    {
        size_t n_words( len/WL );
        ComplementWords< T_CODE >( buf, n_words );
        buf += n_words;
        len -= n_words*WL;
    }

    for( off = 0; len > 0; --len, ++off )
    {
        SetLetter< T_CODE >(
                *buf, off,
                (T_Word)Code::COMPL[GetLetter< T_CODE >( *buf, off )] );
    }
}

//------------------------------------------------------------------------------
/** Reverse and complement a single word.

    \tparam T_CODE  Sequence encoding.
    \tparam T_Word  Word type.

    \param [in,out] w   The word to reverse-complement.
*/
template< ECoding T_CODE, typename T_Word >
void ReverseComplement( T_Word & w )
{
    Reverse< T_CODE >( w );
    Complement< T_CODE >( w );
}

//------------------------------------------------------------------------------
template< ECoding T_CODE, typename T_Word > struct CReverseComplementer
{
    static void ReverseComplement( T_Word * buf, TSeqOff off, TSeqLen len )
    {
        Reverse< T_CODE >( buf, off, len );
        Complement< T_CODE >( buf, off, len );
    }
};

/** Reverse and complement a nucleotide sequence.

    \tparam T_CODE  Sequence encoding.
    \tparam T_Word  Word type.

    \param [in,out] buf Buffer containing the sequence.
    \param [in]     off Offset of the start of the sequence.
    \param [in]     len Length of the sequence (in bases).
*/
template< ECoding T_CODE, typename T_Word >
inline void ReverseComplement( T_Word * buf, TSeqOff off, TSeqLen len )
{
    CReverseComplementer< T_CODE, T_Word >::ReverseComplement(
            buf, off, len );
}

SEQ_NS_END

#endif

