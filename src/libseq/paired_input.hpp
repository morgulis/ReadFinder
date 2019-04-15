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

/** \file libseq/paired_input.hpp
    \brief Combines two single mate input objects into a paired input.
*/

#ifndef LIBSEQ_PAIRED_INPUT_HPP
#define LIBSEQ_PAIRED_INPUT_HPP

#include <libseq/seqinput.hpp>
#include <libseq/defs.hpp>

SEQ_NS_BEGIN

//==============================================================================
template< typename T_SeqInput_1, typename T_SeqInput_2 >
class CPairedInput : public CSeqInput
{
public:

    typedef T_SeqInput_1 SeqInput_1;
    typedef T_SeqInput_2 SeqInput_2;

    CPairedInput( SeqInput_1 * seq_1, SeqInput_2 * seq_2 );
    virtual ~CPairedInput() {}

    virtual bool Next() override;
    virtual SeqData const & GetSeqData() const override;

private:

    std::string CombinedIds( 
            std::string const & id1, std::string const id2 ) const;

    std::unique_ptr< SeqInput_1 > seq_1_;
    std::unique_ptr< SeqInput_2 > seq_2_;
    SeqData sd_;
};

//------------------------------------------------------------------------------
namespace
{

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

}

//------------------------------------------------------------------------------
template< typename T_SeqInput_1, typename T_SeqInput_2 >
inline CPairedInput< T_SeqInput_1, T_SeqInput_2 > * CreatePairedInput(
        T_SeqInput_1 * seq_1, T_SeqInput_2 * seq_2 )
{
    return new CPairedInput< T_SeqInput_1, T_SeqInput_2 >( seq_1, seq_2 );
}

//------------------------------------------------------------------------------
template< typename T_SeqInput_1, typename T_SeqInput_2 >
inline CPairedInput< T_SeqInput_1, T_SeqInput_2 >::CPairedInput(
        SeqInput_1 * seq_1, SeqInput_2 * seq_2 )
    : seq_1_( seq_1 ), seq_2_( seq_2 ), sd_( 2 )
{
}

//------------------------------------------------------------------------------
template< typename T_SeqInput_1, typename T_SeqInput_2 >
inline auto CPairedInput< T_SeqInput_1, T_SeqInput_2 >::GetSeqData() const 
    -> SeqData const &
{
    return sd_;
}

//------------------------------------------------------------------------------
template< typename T_SeqInput_1, typename T_SeqInput_2 >
inline bool CPairedInput< T_SeqInput_1, T_SeqInput_2 >::Next()
{
    bool r1( seq_1_->Next() ),
         r2( seq_2_->Next() );

    if( r1 != r2 )
    {
        M_THROW( "mate streams have different number of sequences" );
    }

    if( !r1 )
    {
        return false;
    }

    sd_.SetId( CombineIds( seq_1_->GetSeqData().GetId(),
                           seq_2_->GetSeqData().GetId() ) );
    sd_.SetData( 0, seq_1_->GetSeqData().GetData( 0 ) );
    sd_.SetData( 1, seq_2_->GetSeqData().GetData( 0 ) );
    return true;
}

SEQ_NS_END

#endif

