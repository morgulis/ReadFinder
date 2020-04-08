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
/** \file libseq/sequtil.hpp
    \brief Various nucleotide sequence utilities.
*/

#ifndef LIBSEQ_SEQUTIL_HPP
#define LIBSEQ_SEQUTIL_HPP

// #include <libseq/align.hpp>
#include <libseq/coding.hpp>
#include <libseq/defs.hpp>

#include <libtools/exception.hpp>

SEQ_NS_BEGIN

using namespace TOOLS_NS;

//------------------------------------------------------------------------------
/** Parse a string into a strand value.

    \param [in]  s       Input string.
    \param [out] val     Parsed strand value.
    \param [in]  fwdchar Character parsed as forward strand.
    \param [in]  revchar Character parsed as reverse strand.

    \throws std::runtime_error if the input string is different from 
    \c \b fwdchar and \c \b revchar.
*/
inline void FromString( 
        std::string const & s, EStrand & val, char fwdchar, char revchar );

//------------------------------------------------------------------------------
/** Convert mate index (0 or 1) to the corresponding EMate value. */
inline EMate ToMate( uint8_t mate_idx );

/** Convert EMate value to a mate index (0 or 1). */
inline constexpr uint8_t FromMate( EMate mate );

//==============================================================================
// IMPLEMENTATION
//==============================================================================
inline void FromString( 
        std::string const & s, EStrand & val, char fwdchar, char revchar )
{
    if( s.size() == 1 && s[0] == fwdchar )
    {
        val = eFWD;
    }
    else if( s.size() == 1 && s[0] == revchar )
    {
        val = eREV;
    }
    else
    {
        M_THROW( "conversion from " << s << " to EStrand failed; " <<
                 "allowed values: " << fwdchar << ',' << revchar );
    }
}

//------------------------------------------------------------------------------
inline EMate ToMate( uint8_t mate_idx )
{
    return mate_idx + 1;
}

//------------------------------------------------------------------------------
inline EStrand ToStrand( uint8_t strand_idx )
{
    return strand_idx + 1;
}

//------------------------------------------------------------------------------
inline constexpr uint8_t FromMate( EMate mate )
{
    return mate - 1;
}

//------------------------------------------------------------------------------
inline constexpr uint8_t StrandIdx( EStrand s ) { return s - 1; }

//==============================================================================
/** Compute a poly A tail of a sequence.

    \tparam T_CODE  Sequence encoding.
    \tparam T_Word  Word type of sequence storage.
*/
template< ECoding T_CODE, typename T_Word > struct CPolyATailComputer;

/** Specialization for NCBI2NA encoding. */
template< typename T_Word > struct CPolyATailComputer< eNCBI2NA, T_Word >
{
    /** Compute poly A tail of a sequence.

        \param [in] seq     Pointer to the sequence data.
        \param [in] mask    Pointer to the mask data.
        \param [in] start   Offset of the start of the sequence (in bases).
        \param [in] len     Length of the sequence (in bases).
        \param [in] w       Minimum tail of consider (in bases).
        \param [in] r       Threshold fraction of A bases in the tail.

        \return Length of the poly A tail of the sequence.
    */
    static TSeqLen GetPolyATail( 
            T_Word const * seq, T_Word const * mask,
            TSeqOff start, TSeqLen len, TSeqLen w, float r );
};

template< typename T_Word >
TSeqLen CPolyATailComputer< eNCBI2NA, T_Word >::GetPolyATail(
        T_Word const * seq, T_Word const * mask,
        TSeqOff start, TSeqLen len, TSeqLen w, float r )
{
    typedef CWord< T_Word > Word;
    typedef CCode< eNCBI2NA > Code;
    static size_t WBITS = Word::N_BITS;
    static TSeqLen const WUNITS = Code::template WordBases< Word >::V;

    assert( w <= len );

    TSeqLen res( w );
    TSeqOff off( start + len - w );
    T_Word wseq, wmask;
    size_t n_A( 0 );

    {
        TSeqOff off_1( off );

        for( ; off_1 <= start + len - WUNITS; off_1 += WUNITS )
        {
            wseq = GetWord( seq, (off_1<<1) ),
            wmask = GetWord( mask, (off_1<<1) );
            wseq |= wmask;
            wseq |= (wseq>>1);
            wseq &= CMatchMask< Code::UNIT, T_Word >::V;
            n_A += WUNITS - GetNSet( wseq );
        }

        if( off_1 < start + len )
        {
            wseq = GetWord( seq, (off_1<<1) ),
            wmask = GetWord( mask, (off_1<<1) );
            off_1 %= WUNITS;
            off_1 <<= 1;
            wseq = GetField( wseq, off_1, WBITS );
            wmask = GetField( wmask, off_1, WBITS );
            wmask |= MaskBits< T_Word >( 0, off_1 );
            wseq |= wmask;
            wseq |= (wseq>>1);
            wseq &= CMatchMask< Code::UNIT, T_Word >::V;
            n_A += WUNITS - GetNSet( wseq );
        }
    }

    if( n_A < res*r )
    {
        TSeqOff e( start + len );

        for( ; off < e && n_A < res*r; ++off, --res )
        {
            if( GetLetter< eNCBI2NA >( seq, off ) == 0 &&
                GetLetter< eNCBI2NA >( mask, off ) == Code::CLEAR_BASE )
            {
                --n_A;
            }
        }

        return std::min( res, len );
    }

    while( off > start )
    {
        wseq = GetWord( seq, ((off - WUNITS)*2) );
        wmask = GetWord( mask, ((off - WUNITS)*2) );
        wseq |= wmask;
        wseq |= (wseq>>1);
        wseq &= CMatchMask< Code::UNIT, T_Word >::V;
        ReverseWord< 1 >( wseq );
        size_t s( GetNFirstClear( wseq )>>1 );
        n_A += s;
        res += s;
        off -= s;
        wseq = GetWord( seq, ((off - WUNITS)*2) );
        wmask = GetWord( mask, ((off - WUNITS)*2) );
        wseq |= wmask;
        wseq |= (wseq>>1);
        wseq &= CMatchMask< Code::UNIT, T_Word >::V;
        wseq |= (~CMatchMask< Code::UNIT, T_Word >::V);
        ReverseWord< 1 >( wseq );
        s = (GetNFirstSet( wseq )>>1);
        res += s;

        if( n_A < res*r )
        {
            for( ; n_A < res*r; --res );
            return std::min( res, len );
        }

        off -= s;
    }

    return len;
}

//------------------------------------------------------------------------------
/** Compute poly A tail of a sequence.

    \tparam T_CODE  Sequence encoding.
    \tparam T_Word  Word type of sequence storage.

    \param [in] seq     Pointer to the sequence data.
    \param [in] mask    Pointer to the mask data.
    \param [in] start   Offset of the start of the sequence (in bases).
    \param [in] len     Length of the sequence (in bases).
    \param [in] w       Minimum tail of consider (in bases).
    \param [in] r       Threshold fraction of A bases in the tail.

    In the first \c w bases find the longest tail of length \c len
    such that number \c n_A of A bases satisfies <tt>n_A >= len*r</tt>
    If \c len is less than \c w then return len. Otherwise look at the
    tails of length \c len \c w bases and higher until the relation
    <tt>n_A >= len*r</tt> holds and return the length of maximum such 
    tail.

    \return Length of the poly A tail of the sequence.
*/
template< ECoding T_CODE, typename T_Word >
inline TSeqLen GetPolyATail( T_Word const * seq, T_Word const * mask, 
                      TSeqOff start, TSeqLen len, TSeqLen w, float r )
{
    return CPolyATailComputer< T_CODE, T_Word >::GetPolyATail( 
            seq, mask, start, len, w, r );
}

SEQ_NS_END

#endif

