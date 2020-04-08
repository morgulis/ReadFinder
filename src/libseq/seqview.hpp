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

#ifndef LIBSEQ_SEQVIEW_HPP
#define LIBSEQ_SEQVIEW_HPP

#include <libseq/coding.hpp>
#include <libseq/sequtil.hpp>
#include <libseq/defs.hpp>

SEQ_NS_BEGIN

//------------------------------------------------------------------------------
/** Base class for sequence view classes with compile time defined coding
    and word types.

    \tparam T_CODE  Sequence encoding.
    \tparam T_Word  Data word type.
*/
template< ECoding T_CODE, typename T_Word >
class CFixedSeqViewBase
{
public:

    typedef CWord< T_Word > Word;   ///< Word traits.
    typedef CCode< T_CODE > Code;   ///< Encoding traits.

    static EDataUnit const UNIT = Code::UNIT;   ///< Data unit size.

    typedef CUnit< UNIT > Unit; ///< Data unit traits.

    /// Type connecting unit and word traits.
    typedef typename Unit::template WordUnits< Word > WordUnits;

    /// Number of units (bases) per word.
    static size_t const WUNITS = WordUnits::N_UNITS;

    /// Alias for sequence encoding.
    static ECoding const CODING = Code::CODING;

    /// Default instance constructor.
    CFixedSeqViewBase(){}

    /** Instance constructor.

        \param [in] sz  Sequence length in bases.
        \param [in] off Offset in bases from the start of the buffer.
    */
    CFixedSeqViewBase( TSeqLen sz, TSeqOff off );

    TSeqLen size() const;       ///< Get the sequence view length in bases.
    TSeqOff GetStart() const;   ///< Get the starting offset in bases.

protected:

    TSeqLen sz_ = 0;    ///< Sequence view length.
    TSeqOff off_ = 0;   ///< Starting offset.
};

//------------------------------------------------------------------------------
template< ECoding T_CODE, typename T_Word > class CFixedConstSeqView;

/** Mutable sequence view with compile time coding and word type. 

    \tparam T_CODE  Sequence encoding.
    \tparam T_Word  Data word type.
*/
template< ECoding T_CODE, typename T_Word > 
class CFixedSeqView : public CFixedSeqViewBase< T_CODE, T_Word >
{
    friend class CFixedConstSeqView< T_CODE, T_Word >;

private:

    typedef CFixedSeqViewBase< T_CODE, T_Word > TBase;

public:

    /// Default instance constructor.
    CFixedSeqView(){}

    /** Instance constructor.

        \param [in] buf Underlying raw sequence data.
        \param [in] off Offset in bases from the start of the buffer.
        \param [in] len Sequence length in bases.
    */
    CFixedSeqView( T_Word * buf, ssize_t off, TSeqLen len );

    T_Word * GetBuf() const;    ///< Get the underlying buffer.

    /** Advance the start of the view by the given number of bases.

        \param [in] off Increment value.
    */
    CFixedSeqView & operator+=( TSeqOff off )
    {
        this->off_ += off;
        this->sz_ -= off;
        return *this;
    }

    /** Return a new view obtained by adding a given increment (in bases)
        to the existing view.

        \param [in] v   Source view.
        \param [in] off Increment value.

        \return Resulting view.
    */
    friend CFixedSeqView operator+( CFixedSeqView const & v, TSeqOff off )
    {
        CFixedSeqView res( v );
        res += off;
        return res;
    }

private:

    T_Word * buf_ = nullptr;    ///< Raw sequence data.
};

//------------------------------------------------------------------------------
/** Immutable sequence view with compile time coding and word type. 

    \tparam T_CODE  Sequence encoding.
    \tparam T_Word  Data word type.
*/
template< ECoding T_CODE, typename T_Word > 
class CFixedConstSeqView : public CFixedSeqViewBase< T_CODE, T_Word >
{
private:

    typedef CFixedSeqViewBase< T_CODE, T_Word > TBase;

public:

    /// Default instance constructor.
    CFixedConstSeqView(){}

    /** Instance constructor.

        \param [in] buf Underlying raw sequence data.
        \param [in] off Offset in bases from the start of the buffer.
        \param [in] len Sequence length in bases.
    */
    CFixedConstSeqView( T_Word const * buf, ssize_t off, TSeqLen len );

    /** Conversion from the corresponding mutable view. */
    CFixedConstSeqView( CFixedSeqView< T_CODE, T_Word > const & other );

    T_Word const * GetBuf() const;  ///< Get the underlying sequence data.

    /** Advance the start of the view by the given number of bases.

        \param [in] off Increment value.
    */
    CFixedConstSeqView & operator+=( TSeqOff off )
    {
        this->off_ += off;
        this->sz_ -= off;
        return *this;
    }

    /** Return the base value at a given offset in the view.

        \param [in] off Offset relative to the start of the view.

        \return The base value at the requested offset.
    */
    T_Word GetLetter( TSeqOff off ) const;

    /** Return the base value at a given offset in the view in the
        specified encoding.

        \tparam T_OTHER_CODE    Target encoding.

        \param [in] off Offset relative to the start of the view.

        \return The base value at the requested offset.
    */
    template< ECoding T_OTHER_CODE > T_Word GetLetter( TSeqOff off ) const;

    /** Return a new view obtained by adding a given increment (in bases)
        to the existing view.

        \param [in] v   Source view.
        \param [in] off Increment value.

        \return Resulting view.
    */
    friend CFixedConstSeqView operator+( 
            CFixedConstSeqView const & v, TSeqOff off )
    {
        CFixedConstSeqView res( v );
        res += off;
        return res;
    }

private:

    T_Word const * buf_ = nullptr;  ///< Raw sequence data.
};

//------------------------------------------------------------------------------
/** Copy data from source view to destination view with possible 
    conversion between encodings.

    For bases not representable in destination encoding the corresponding
    positions are masked in mask view.
*/
template< typename T_DstView, typename T_SrcView >
bool Recode( 
        T_DstView const & dview, T_DstView const & dmaskview, 
        T_SrcView const & sview, TSeqLen len = -1 );

//------------------------------------------------------------------------------
/** Copy data from the source view to the destination view. */
template< typename T_DstView, typename T_SrcView >
void Copy( T_DstView && dview, T_SrcView const & sview );

//------------------------------------------------------------------------------
/** Reverse bases in the view. */
template< typename T_View > void Reverse( T_View const & view );

/** Reverse and complement bases in the view. */
//------------------------------------------------------------------------------
template< typename T_View > void ReverseComplement( T_View const & view );

//------------------------------------------------------------------------------
/** Compute the length polyA tail in the sequence data.

    \tparam T_View  Sequence view type.

    \param [in] seq     Sequence data.
    \param [in] mask    Sequence mask.
    \param [in] w       Minimum tail to consider in bases.
    \param [in] r       Threshold fraction of A bases in the tail.

    \return The length of the polyA tail.
*/
template< typename T_View > TSeqLen GetPolyATail( 
        T_View const & seq, T_View const & mask, TSeqLen w, float r );

//------------------------------------------------------------------------------
/** Get the word of sequence data at the given offset (in bases) from the
    start of the view.
*/
template< typename T_View > 
typename T_View::Word::T GetWord( T_View const & view, TSeqOff off );

//------------------------------------------------------------------------------
/** Find the length of the initial matching segment of two sequences */
template< typename T_View > TSeqLen AlignExact(
        T_View const & seq1, T_View const & mask1,
        T_View const & seq2, T_View const & mask2 );

//------------------------------------------------------------------------------
/** Find the length of the initial matching segment of the first sequence and
    the reverse complement of the second sequence. 
*/
template< typename T_View > TSeqLen AlignExactReverse(
        T_View const & seq1, T_View const & mask1,
        T_View const & seq2, T_View const & mask2 );

//==============================================================================
// IMPLEMENTATION
//==============================================================================
template< ECoding T_CODE, typename T_Word >
inline CFixedSeqViewBase< T_CODE, T_Word >::CFixedSeqViewBase( 
        TSeqLen sz, TSeqOff off )
    : sz_( sz ), off_( off )
{}

//------------------------------------------------------------------------------
template< ECoding T_CODE, typename T_Word >
inline TSeqLen CFixedSeqViewBase< T_CODE, T_Word >::size() const
{
    return sz_;
}

//------------------------------------------------------------------------------
template< ECoding T_CODE, typename T_Word >
inline TSeqLen CFixedSeqViewBase< T_CODE, T_Word >::GetStart() const
{
    return off_;
}

//------------------------------------------------------------------------------
template< ECoding T_CODE, typename T_Word >
inline CFixedSeqView< T_CODE, T_Word >::CFixedSeqView( 
        T_Word * buf, ssize_t off, TSeqLen len )
    : TBase( len, off ),
      buf_( buf )
{
}

//------------------------------------------------------------------------------
template< ECoding T_CODE, typename T_Word >
inline auto CFixedSeqView< T_CODE, T_Word >::GetBuf() const -> T_Word *
{
    return buf_;
}

//------------------------------------------------------------------------------
template< ECoding T_CODE, typename T_Word >
inline CFixedConstSeqView< T_CODE, T_Word >::CFixedConstSeqView( 
        T_Word const * buf, ssize_t off, TSeqLen len )
    : TBase( len, off ),
      buf_( buf )
{
}
 
//------------------------------------------------------------------------------
template< ECoding T_CODE, typename T_Word >
inline CFixedConstSeqView< T_CODE, T_Word >::CFixedConstSeqView( 
        CFixedSeqView< T_CODE, T_Word > const & other )
    : TBase( other ),
      buf_( other.buf_ )
{
}

//------------------------------------------------------------------------------
template< ECoding T_CODE, typename T_Word >
inline auto CFixedConstSeqView< T_CODE, T_Word >::GetBuf() const 
    -> T_Word const *
{
    return buf_;
}

//------------------------------------------------------------------------------
template< ECoding T_CODE, typename T_Word >
inline auto CFixedConstSeqView< T_CODE, T_Word >::GetLetter( 
        TSeqOff off ) const -> T_Word
{
    return SEQ_NS::GetLetter< T_CODE >( buf_, this->off_ + off );
}

//------------------------------------------------------------------------------
template< ECoding T_CODE, typename T_Word >
template< ECoding T_OTHER_CODE >
inline auto CFixedConstSeqView< T_CODE, T_Word >::GetLetter( 
        TSeqOff off ) const -> T_Word
{
    return CRecoder< T_OTHER_CODE, T_CODE >::Recode( GetLetter( off ) );
}

//------------------------------------------------------------------------------
template< typename T_DstView, typename T_SrcView >
inline bool Recode( 
        T_DstView const & dview, T_DstView const & dmaskview, 
        T_SrcView const & sview, TSeqLen len )
{
    if( len < 0 ) 
    {
        len = sview.size();
    }

    return Recode< T_DstView::CODING, T_SrcView::CODING >(
            dview.GetBuf(), dmaskview.GetBuf(), sview.GetBuf(),
            dview.GetStart(), dmaskview.GetStart(), sview.GetStart(), len );
}

//------------------------------------------------------------------------------
template< typename T_DstView, typename T_SrcView >
inline void Copy( T_DstView && dview, T_SrcView const & sview )
{
    if( dview.GetStart() == 0 && sview.GetStart() == 0 )
    {
        CopyData< T_SrcView::Code::UNIT >( 
                dview.GetBuf(), sview.GetBuf(), sview.size() );
    }
    else
    {
        CopyData< T_SrcView::Code::UNIT >(
                dview.GetBuf(), sview.GetBuf(), 
                dview.GetStart(), sview.GetStart(), sview.size() );
    }
}

//------------------------------------------------------------------------------
template< typename T_View > inline void Reverse( T_View const & view )
{
    Reverse< T_View::CODING >( view.GetBuf(), view.GetStart(), view.size() );
}

//------------------------------------------------------------------------------
template< typename T_View > 
inline void ReverseComplement( T_View const & view )
{
    ReverseComplement< T_View::CODING >( 
            view.GetBuf(), view.GetStart(), view.size() );
}

//------------------------------------------------------------------------------
template< typename T_View > 
inline TSeqLen GetPolyATail( 
        T_View const & seq, T_View const & mask, TSeqLen w, float r )
{
    typedef typename T_View::Word::T T_Word;
    return GetPolyATail< T_View::CODING, T_Word >( 
            seq.GetBuf(), mask.GetBuf(), seq.GetStart(), seq.size(), w, r );
}

//------------------------------------------------------------------------------
template< typename T_View > 
inline typename T_View::Word::T GetWord( T_View const & view, TSeqLen off )
{
    typedef typename T_View::Code Code;
    return TOOLS_NS::GetWord( 
            view.GetBuf(), (view.GetStart() + off)*Code::LBITS );
}

//------------------------------------------------------------------------------
template< typename T_View > TSeqLen AlignExact(
        T_View const & seq1, T_View const & mask1,
        T_View const & seq2, T_View const & mask2 )
{
    typedef typename T_View::Word::T T_Word;
    static TSeqLen const WL = T_View::WUNITS;

    TSeqLen n_matches( 0 ),
            off( 0 ),
            n_m;

    do 
    {
        T_Word w1( GetWord( seq1, off ) ),
               m1( GetWord( mask1, off ) ),
               w2( GetWord( seq2, off ) ),
               m2( GetWord( mask2, off ) );
        n_m = GetExactMatch< T_View::UNIT >( w1, w2, (T_Word)(m1|m2) );
        n_matches += n_m;
        off += WL;
    }
    while( n_m == WL );

    return n_matches;
}

//------------------------------------------------------------------------------
template< typename T_View > TSeqLen AlignExactReverse(
        T_View const & seq1, T_View const & mask1,
        T_View const & seq2, T_View const & mask2 )
{
    typedef typename T_View::Word::T T_Word;
    typedef typename T_View::Code Code;
    static TSeqLen const WL = T_View::WUNITS;
    static ECoding const CODE = Code::CODING;

    TSeqLen n_matches( 0 ),
            off( 0 ),
            roff( -WL ),
            n_m;

    do
    {
        T_Word w1( GetWord( seq1, off ) ),
               m1( GetWord( mask1, off ) ),
               w2( GetWord( seq2, roff ) ),
               m2( GetWord( mask2, roff ) );
        ReverseComplement< CODE >( w2 );
        Reverse< CODE >( m2 );
        n_m = GetExactMatch< T_View::UNIT >( w1, w2, (T_Word)(m1|m2) );
        n_matches += n_m;
        off += WL;
        roff -= WL;
    }
    while( n_m == WL );

    return n_matches;
}

SEQ_NS_END

#endif

