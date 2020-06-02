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

/** \file libseq/seqinput.hpp
    \brief Base class and factory for reading sequence data represented in
           different formats.
*/

#ifndef LIBSEQ_SEQINPUT_HPP
#define LIBSEQ_SEQINPUT_HPP

#include <string>
#include <vector>

#include <libseq/defs.hpp>

#include <libtools/exception.hpp>

SEQ_NS_BEGIN

//==============================================================================
class CSeqInput
{
public:

    enum : int
    {
        FASTA = 0,
        ZFASTA,
        FASTQ,
        ZFASTQ,
        SRA
    };

protected:

    class SeqData
    {
    private:

        static size_t const MAX_COLS = 2;

    public:

        SeqData( size_t n_cols = 1 );

        size_t GetNCols() const;
        std::string const & GetId() const;
        std::string const & GetData( size_t col ) const;

        void SetId( std::string const & id );
        void SetData( size_t col, std::string const & data );

        void SetNCols( size_t n_cols )
        {
            assert( n_cols <= MAX_COLS );
            n_cols_ = n_cols;
        }

    private:

        size_t n_cols_;
        std::string id_;
        std::string data_[MAX_COLS];
    };

    static bool CheckCompressed( std::string const & fname )
    {
        return fname.size() > 3 && fname.substr( fname.size() - 3 ) == ".gz";
    }

    size_t n_failed_reads_ = 0;

private:

    class Iterator
    {
    public:

        Iterator() {}
        Iterator( CSeqInput * o );
        Iterator begin();
        Iterator end() const;
        Iterator & operator++();

        SeqData const & operator*() const;

        friend bool operator!=( Iterator const & x, Iterator const & y )
        {
            return x.o_ != y.o_;
        }

    private:

        CSeqInput * o_ = nullptr;
    };

public:

    virtual ~CSeqInput() {}

    virtual bool Next() = 0;
    virtual SeqData const & GetSeqData() const = 0;
    size_t GetNumFailedReads() const { return n_failed_reads_; }

    Iterator Iterate();
};

//------------------------------------------------------------------------------
inline CSeqInput::SeqData::SeqData( size_t n_cols )
    : n_cols_( n_cols )
{
    assert( n_cols_ > 0 );
    assert( n_cols_ <= MAX_COLS );
}

//------------------------------------------------------------------------------
inline size_t CSeqInput::SeqData::GetNCols() const
{
    return n_cols_;
}

//------------------------------------------------------------------------------
inline std::string const & CSeqInput::SeqData::GetId() const
{
    return id_;
}

//------------------------------------------------------------------------------
inline void CSeqInput::SeqData::SetId( std::string const & id )
{
    id_ = id;
}

//------------------------------------------------------------------------------
inline std::string const & CSeqInput::SeqData::GetData( size_t col ) const
{
    assert( col < n_cols_ );
    return data_[col];
}

//------------------------------------------------------------------------------
inline void CSeqInput::SeqData::SetData( size_t col, std::string const & data )
{
    assert( col < n_cols_ );
    data_[col] = data;
}

//------------------------------------------------------------------------------
inline CSeqInput::Iterator::Iterator( CSeqInput * o ) : o_( o )
{
    assert( o_ != nullptr );

    if( !(o_->Next()) )
    {
        o_ = nullptr;
    }
}

//------------------------------------------------------------------------------
inline auto CSeqInput::Iterator::begin() -> Iterator
{
    return *this;
}

//------------------------------------------------------------------------------
inline auto CSeqInput::Iterator::end() const -> Iterator
{
    return Iterator();
}

//------------------------------------------------------------------------------
inline auto CSeqInput::Iterator::operator++() -> Iterator &
{
    if( *this != end() )
    {
        if( !(o_->Next()) )
        {
            o_ = nullptr;
        }
    }

    return *this;
}

//------------------------------------------------------------------------------
inline auto CSeqInput::Iterator::operator*() const -> SeqData const &
{
    assert( *this != end() );
    return o_->GetSeqData();
}

//------------------------------------------------------------------------------
inline auto CSeqInput::Iterate() -> Iterator
{
    return Iterator( this );
}

//==============================================================================
//------------------------------------------------------------------------------
CSeqInput * MkSeqInput( 
        std::vector< std::string > const & fnames, int fmt,
        size_t start, ssize_t n_seq, bool paired = false );

//------------------------------------------------------------------------------
inline std::string CombineIds( std::string const & id1, std::string const id2 )
{
    if( id1.size() == id2.size() )
    {
        size_t i( 0 );
        for( ; i < id1.size() && id1[i] == id2[i]; ++i );

        if( i == id1.size() )
        {
            return id1;
        }

        if( i > 0 &&
            i == id1.size() - 1 &&
            id1[i] == '1' &&
            id2[i] == '2' &&
            (id1[i-1] == '.' || id1[i-1] == '_' || id1[i-1] == '/') )
        {
            return id1.substr( 0, i - 1 );
        }
    }

    M_THROW( "id mismatch: " << id1 << " and " << id2 );
}


SEQ_NS_END

#endif

