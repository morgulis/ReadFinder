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

#ifndef LIBREADFINDER_FAST_SEEDS_HPP
#define LIBREADFINDER_FAST_SEEDS_HPP

#include <deque>

#include <boost/dynamic_bitset.hpp>

#include <config.h>

#include <libreadfinder/batch.hpp>
#include <libreadfinder/fsdefs.hpp>
#include <libreadfinder/fsidx.hpp>
#include <libreadfinder/defs.hpp>

READFINDER_NS_BEGIN

//==============================================================================
class CFastSeeds : public CFastSeedsDefs
{
private:

    static size_t const ANCHOR_TBL_SIZE = (1ULL<<ANCHOR_BITS);
    static size_t const WMAP_SIZE = 1ULL + (1ULL<<ANCHOR_BITS);

    //--------------------------------------------------------------------------
    struct WordCountingJobData;
    struct WordCountingJob;

    //--------------------------------------------------------------------------
    typedef boost::dynamic_bitset< uint64_t > BitSet;

    //--------------------------------------------------------------------------
    struct WordTableGeneratorJob;
    struct WordTableSortJob;

    static_assert( sizeof( HashWord ) == 12, "" );

    struct ExtHashWord
    {
        HashWord hw;
    };

    typedef std::vector< HashWord > WordTable;

    typedef CReadData::HashWordSource HashWordSource;
    typedef CReadData::HashMaskSource HashMaskSource;

    struct SeedSearchJobData;
    struct SeedSearchJob;

    //--------------------------------------------------------------------------
    struct Hit
    {
        uint32_t refpos;

        union {
            OrdId read;
            TRefOId refseq;
            uint32_t len;
        };

        TReadOff readpos;
        uint8_t strand : 1;
        uint8_t mate : 1;

        TSeqOff GetDiag() const { return refpos - readpos; }

        friend std::ostream & operator<<( std::ostream & os, Hit const & h )
        {
            os << h.refpos << ':' << h.read << ':' << h.readpos << ':'
               << (int)h.strand << ':' << (int)h.mate;
            return os;
        }
    };

    static_assert( sizeof( Hit ) == 12, "" );

    struct Hits
    {
        static constexpr size_t const MAX_BUCKETS = 1024ULL;
        typedef std::deque< Hit > Bucket;
        typedef std::vector< Bucket > Data;

        Data data = Data( MAX_BUCKETS );
    };

    typedef HashWord TaskEntry;

    struct HitFilteringJob;
    struct ReadMarkingJob;

public:

    CFastSeeds( CBatch & bctx );
    void Run();

private:

    void UpdateAnchorUseMap( uint32_t wa );
    void CreateWordTable();
    void ComputeSeeds();

    CRefData const & GetRefs() const { return *bctx_.GetSearchCtx().refs; }
    CReadData const & GetReads() const { return bctx_.GetReads(); }

    CBatch & bctx_;
    CFastSeedsIndex & fsidx_;
    WordTable wt_;
    WordMap wmap_ = WordMap( WMAP_SIZE, 0 );
    BitSet anchor_use_map_;
    size_t reads_per_bucket_ = 0;
    OrdId start_oid_ = 0;
    bool prescreen_ = false;
};

READFINDER_NS_END

#endif

