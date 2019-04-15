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

#if 0
//------------------------------------------------------------------------------
/** Iterate over alignment events encoded as CIGAR string. */
class CCigarIterator
{
private:

    static CCigarIterator const End;

public:

    /// Default instance constructor.
    CCigarIterator() {}

    /** Instance constructor.

        Copies its argument into internal string.

        \param [in] cigar   CIGAR string to parse.
    */
    CCigarIterator( std::string const & cigar );

    /** Instance constructor.

        Moves its argument into internal string.

        \param [in] cigar   CIGAR string to parse.
    */
    CCigarIterator( std::string && cigar );

    /** Iterate to the next event.

        \return false if at the end of iteration; 
                true otherwise.
    */
    bool Next();

    /// \name Interface for 'range based for' type iteration.
    /// @{
    CCigarIterator & begin();
    CCigarIterator const & end();
    CCigarIterator & operator++();
    CAlignBase::Event & operator*();

    friend bool operator!=( 
            CCigarIterator const & x, CCigarIterator const & y )
    {
        if( x.e_.type == CAlignBase::Event::eNone ||
            y.e_.type == CAlignBase::Event::eNone )
        {
            return x.e_.type != y.e_.type;
        }
        else
        {
            return &x != &y;
        }
    }

    friend bool operator==( 
            CCigarIterator const & x, CCigarIterator const & y )
    {
        return !(x != y);
    }
    /// @}

private:

    void Init();    ///< Common initialization for both constructors.

    std::string cigar_;                 ///< CIGAR string to parse.
    std::string::size_type pos_ = 0;    ///< Current position in cigar_.

    /// Current event.
    CAlignBase::Event e_ = CAlignBase::Event{ 0, CAlignBase::Event::eNone },
                      pe_;  ///< Previous event.
};
#endif

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

/*
//------------------------------------------------------------------------------
inline void CCigarIterator::Init()
{
    Next();
}

//------------------------------------------------------------------------------
inline CCigarIterator::CCigarIterator( std::string const & cigar )
    : cigar_( cigar )
{
    Init();
}

//------------------------------------------------------------------------------
inline CCigarIterator::CCigarIterator( std::string && cigar )
    : cigar_( std::move( cigar ) )
{
    Init();
}

//------------------------------------------------------------------------------
inline CCigarIterator & CCigarIterator::begin()
{
    return *this;
}

//------------------------------------------------------------------------------
inline CCigarIterator const & CCigarIterator::end()
{
    return End;
}

//------------------------------------------------------------------------------
inline CCigarIterator & CCigarIterator::operator++()
{
    Next();
    return *this;
}

//------------------------------------------------------------------------------
inline CAlignBase::Event & CCigarIterator::operator*()
{
    return e_;
}

//------------------------------------------------------------------------------
inline bool CCigarIterator::Next()
{
    if( pos_ >= cigar_.size() )
    {
        e_.type = CAlignBase::Event::eNone;
        return false;
    }

    pe_ = e_;

    if( !isdigit( cigar_[pos_] ) )
    {
        M_THROW( "error parsing CIGAR string: " << cigar_ <<
                 ": digit expected at position " << pos_ );
    }

    e_.len = 0;

    for( ; pos_ < cigar_.size() && isdigit( cigar_[pos_] ); ++pos_ )
    {
        e_.len = e_.len*10 + (cigar_[pos_] - '0');
    }

    if( e_.len == 0 )
    {
        M_THROW( "error parsing CIGAR string: " << cigar_ <<
                 ": event length must be positive (at position " << pos_ <<
                 ")" );
    }

    if( pos_ >= cigar_.size() )
    {
        M_THROW( "error parsing CIGAR string: " << cigar_ <<
                 ": unexpected end of string" );
    }

    switch( cigar_[pos_++] )
    {
        case '=': e_.type = CAlignBase::Event::eMatch; break;
        case 'X':
        case 'R': e_.type = CAlignBase::Event::eMismatch; break;
        case 'I': e_.type = CAlignBase::Event::eInsert; break;
        case 'D': e_.type = CAlignBase::Event::eDelete; break;
        case 'S': e_.type = CAlignBase::Event::eSoftClip; break;
        case 'H': e_.type = CAlignBase::Event::eHardClip; break;
        case 'N': e_.type = CAlignBase::Event::eSkip; break;
        default:
            M_THROW( "error parsing CIGAR string: " << cigar_ <<
                     ": unsupported event type tag " << cigar_[pos_-1] <<
                     " at position " << pos_ - 1 );
    }

    if( pe_.type == e_.type )
    {
        M_THROW( "error parsing CIGAR string: "<< cigar_ << 
                 ": neighbor events of the same type: " <<
                 "at position " << pos_ );
    }

    return true;
}

//==============================================================================
class CCigarGenerator
{
public:

    CCigarGenerator & operator<<( CAlignBase::Event e );
    std::string ToString() const;

    template< typename T_Iter >
    static std::string MakeCigar( T_Iter start, T_Iter end );

private:

    std::string str_;
    CAlignBase::Event pe_;
    bool start_ = true;
};

#define M_EVENT_LETTERS "=XIDSHN-";

//------------------------------------------------------------------------------
inline CCigarGenerator & CCigarGenerator::operator<<( CAlignBase::Event e )
{
    static char const EVENT_LETTERS[] = M_EVENT_LETTERS;

    if( start_ )
    {
        pe_ = e;
        start_ = false;
    }
    else if( pe_.type != e.type )
    {
        str_ += std::to_string( pe_.len );
        str_.append( 1, EVENT_LETTERS[pe_.type] );
        pe_ = e;
    }
    else
    {
        pe_.len += e.len;
    }

    return *this;
}

//------------------------------------------------------------------------------
inline std::string CCigarGenerator::ToString() const
{
    static char const EVENT_LETTERS[] = M_EVENT_LETTERS;
    std::string res( str_ );

    if( !start_ )
    {
        res += std::to_string( pe_.len );
        res.append( 1, EVENT_LETTERS[pe_.type] );
    }

    return res;
}

#undef M_EVENT_LETTERS

//------------------------------------------------------------------------------
template< typename T_Iter >
inline std::string CCigarGenerator::MakeCigar( T_Iter start, T_Iter end )
{
    CCigarGenerator gen;

    for( ; start != end; ++start )
    {
        gen << *start;
    }

    return gen.ToString();
}

//------------------------------------------------------------------------------
template< typename T_Iter >
inline std::string MakeBTOP( 
        T_Iter start, T_Iter end, unsigned char const * letters )
{
    static const char EVENT_LETTERS[] = "=XIDSHN";
    static const bool EVENT_ADD_LETTERS[] = { 
        false, true, false, true, false, false, true };

    std::string res;
    uint8_t pe( CAlignBase::Event::eNone );
    TSeqLen len( 0 );

    for( ; start != end; ++start )
    {
        if( start->type != pe )
        {
            if( len > 0 )
            {
                res += std::to_string( len ) + EVENT_LETTERS[pe];

                if( EVENT_ADD_LETTERS[pe] )
                {
                    TSeqLen l( pe == CAlignBase::Event::eSkip ? 4 : len );

                    for( TSeqLen ll( 0 ); ll < l; ++ll )
                    {
                        res.append( 1, *letters++ );
                    }
                }

                len = 0;
            }

            pe = start->type;
        }

        len += start->len;
    }

    if( len > 0 )
    {
        res += std::to_string( len ) + EVENT_LETTERS[pe];

        if( EVENT_ADD_LETTERS[pe] )
        {
            TSeqLen l( pe == CAlignBase::Event::eSkip ? 4 : len );

            for( TSeqLen ll( 0 ); ll < l; ++ll )
            {
                res.append( 1, *letters++ );
            }
        }
    }

    return res;
}

//------------------------------------------------------------------------------
template< typename T_Iter >
inline std::string MakeMD( 
        T_Iter start, T_Iter end, unsigned char const * letters )
{
    size_t len( 0 );
    std::string res;

    for( ; start != end; ++start )
    {
        switch( start->type )
        {
            case CAlignBase::Event::eMatch:

                len += start->len;
                break;

            case CAlignBase::Event::eMismatch:

                if( len > 0 )
                {
                    res += std::to_string( len );
                    len = 0;
                }

                for( size_t i( 0 ); i < start->len; ++i )
                {
                    res.append( 1, *letters++ );
                }

                break;

            case CAlignBase::Event::eDelete:

                if( len > 0 )
                {
                    res += std::to_string( len );
                    len = 0;
                }

                res.append( 1, '^' );

                for( size_t i( 0 ); i < start->len; ++i )
                {
                    res.append( 1, *letters++ );
                }

                break;

            case CAlignBase::Event::eSkip:

                if( len > 0 )
                {
                    res += std::to_string( len );
                    len = 0;
                }

                res.append( 1, '|' );
                letters += 4;
                break;

            default: break;
        }
    }

    if( len > 0 )
    {
        res += std::to_string( len );
        len = 0;
    }

    return res;
}

//------------------------------------------------------------------------------
template< typename T_Iter >
inline size_t GetEditDistance( T_Iter start, T_Iter end )
{
    size_t res( 0 );

    for( ; start != end; ++start )
    {
        switch( start->type )
        {
            case CAlignBase::Event::eMismatch:
            case CAlignBase::Event::eInsert:
            case CAlignBase::Event::eDelete:

                res += start->len;
                break;

            default: break;
        }
    }

    return res;
}
*/

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

