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

#ifndef LIBREADFINDER_REFDATA_HPP
#define LIBREADFINDER_REFDATA_HPP

#include <vector>
#include <map>
#include <memory>
#include <string>
#include <fstream>

#include <zlib.h>

#include <libreadfinder/defs.hpp>

#include <libtools/exception.hpp>

READFINDER_NS_BEGIN

using namespace SEQ_NS;

//==============================================================================
/** Simple reference database for Readfinder.
*/
class CRefData 
{
public:

    typedef uint32_t TRefOId;   ///< Type to use for ordinal sequence ids.

private:

    /// List of string ids of reference sequences in the input order.
    typedef std::vector< std::string > TIdMap;

    /// Type of mapping from string ids to ordinal ids.
    typedef std::map< std::string, TRefOId > TOIdMap;

    /// Offsets of the start of each sequence data relative to the start
    /// of the database (indexed by ordinal id).
    typedef std::vector< uint64_t > TOffsets;

    /// Lengths of reference sequences (indexed by ordinal id).
    typedef std::vector< TSeqLen > TLengths;

    /// Type to store sequence and ambiguity data.
    typedef std::vector< TWord > TData;

    /// Descriptor of a reference sequence.
    struct SeqInfo
    {
        TData seqdata,      ///< Storage for sequence data.
              maskdata;     ///< Storage for ambiguity data.
        size_t count = 0;   ///< Reference count.
    };

    /** Type used to hold sequence data holders.

        This is a vector of pointers to sequence descriptors (SeqInfo) indexed
        by ordinal id. If the data is not referenced, the corresponding
        pointer is NULL. Otherwise it points to the valid sequence and
        ambuguity data loaded from the reference database.
    */
    typedef std::unique_ptr< SeqInfo > SeqInfoHolder;

public:

    /** Instance constructor.

        Opens the database and loads meta data: sequence ids and offsets.
        The real data need to be loaded by clients via Hold()/Release()
        reference counting methods.

        \param [in] db_name Base name of the database.
    */
    CRefData( std::string const & db_name );

    /// Instance destructor.
    ~CRefData();

    /** Get the length in bases of the reference sequence by its ordinal id. */
    TSeqLen GetLength( TRefOId refid ) const;

    /** Get the string id of the reference sequence by its ordinal id. */
    std::string const & GetRefId( TRefOId oid ) const;

    /** Get the number of reference sequences in the database. */
    TIdMap::size_type GetSize() const;

    /** Acquire access to the sequence data.

        If the requested sequence is not loaded, loads it.
        Increments the reference counter.

        \param [in] refid   Ordinal id of the requested sequence.
    */
    void Hold( TRefOId refid );

    /** Release the sequence data.

        Decrement the reference counter for the requested sequence.
        If the counter reaches 0, unload the sequence data.

        \param [in] refid   Ordinal id of the requested sequence.
    */
    void Release( TRefOId refid );

    /** Load and hold all sequences.
    */
    void LoadAll();

    /** Release all sequences.
    */
    void ReleaseAll();

    /** Get the sequence data in NCBI2NA encoding.

        \param [in] refid   Ordinal id of the requested sequence.

        \return Pointer to the start of the sequence data.
    */
    TWord const * GetSeqData( TRefOId refid ) const;

    /** Get the ambiguity data.

        \param [in] refid   Ordinal id of the requested sequence.

        \return Pointer to the start of the ambiguity data. The data is
                encoded with 2-bits per base, with non-ambiguous bases
                set to 0 and ambiguous bases set to 3.
    */
    TWord const * GetMaskData( TRefOId refid ) const;

    /** Return an estimate of used memory. */
    size_t GetMemUsage() const;

#ifdef TEST_ENABLED
    bool Check( int test_id, void * data = nullptr );
#endif

private:

    /** Load the sequence data from the reference database.

        \param [in] refid   Ordinal id of the requested sequence.
    */
    void Load( TRefOId refid );

    std::string dbname_;    ///< Name (prefix) of the reference database.
    TIdMap id_map_;         ///< String ids indexed by ordinal id.
    TOffsets offsets_;      /**< Sequence offsets within the database 
                                    indexed by ordinal id. */
    TLengths lengths_;      ///< Sequence lengths indexed by ordinal id.
    gzFile data_is_ = nullptr; ///< Sequence data stream.

    /// Sequence data for sequences with positive reference counts.
    std::vector< SeqInfoHolder > data_;
};

//------------------------------------------------------------------------------
inline CRefData::~CRefData() 
{ 
    gzclose( data_is_ );
}

//------------------------------------------------------------------------------
inline TSeqLen CRefData::GetLength( TRefOId refid ) const 
{ 
    return lengths_[refid]; 
}

//------------------------------------------------------------------------------
inline std::string const & CRefData::GetRefId( TRefOId oid ) const 
{ 
    return id_map_[oid]; 
}

//------------------------------------------------------------------------------
inline auto CRefData::GetSize() const -> TIdMap::size_type
{
    return id_map_.size();
}

//------------------------------------------------------------------------------
inline void CRefData::Hold( TRefOId refid )
{
    assert( refid < GetSize() );

    if( data_[refid].get() == nullptr )
    {
        Load( refid );
    }

    ++data_[refid]->count;
}

//------------------------------------------------------------------------------
inline void CRefData::LoadAll()
{
    for( size_t i( 0 ); i < GetSize(); ++i )
    {
        Hold( i );
    }
}

//------------------------------------------------------------------------------
inline void CRefData::Release( TRefOId refid )
{
    assert( refid < GetSize() );
    assert( data_[refid].get() != nullptr );
    assert( data_[refid]->count > 0 );

    if( --data_[refid]->count == 0 )
    {
        data_[refid].reset();
    }
}

//------------------------------------------------------------------------------
inline void CRefData::ReleaseAll()
{
    for( size_t i( 0 ); i < GetSize(); ++i )
    {
        Release( i );
    }
}

//------------------------------------------------------------------------------
inline TWord const * CRefData::GetSeqData( TRefOId refid ) const
{
    assert( refid < GetSize() );
    assert( data_[refid].get() != nullptr );
    return &(data_[refid]->seqdata[0]) + 2;
}

//------------------------------------------------------------------------------
inline TWord const * CRefData::GetMaskData( TRefOId refid ) const
{
    assert( refid < GetSize() );
    assert( data_[refid].get() != nullptr );
    return &(data_[refid]->maskdata[0]) + 2;
}

READFINDER_NS_END

#endif

