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

/** \file libreadfinder/batch.hpp
    \brief Single batch of reads.
*/

#ifndef LIBREADFINDER_BATCH_HPP
#define LIBREADFINDER_BATCH_HPP

#include <condition_variable>
#include <memory>
#include <mutex>

#include <libreadfinder/rf_ctx.hpp>
#include <libreadfinder/readdata.hpp>
// #include <libreadfinder/seeddata.hpp>
#include <libreadfinder/defs.hpp>

READFINDER_NS_BEGIN

//==============================================================================
class CBatch//  : private BatchDefs
{
public:

    struct StatParams
    {
        enum : size_t
        {
            N_READS = 0,
            N_PAIRED_READS,
            N_JOBS,
            MAX_JOB_READS,

            N_SEEDER_HITS,
            N_EXTENDED_SEEDS,
            N_SEEDER_EEXONS,
            N_SEEDER_EXONS,
            N_PRIMARY_SEEDS,
            N_MAPPED_READS,

            N_PARAMS
        };
    };

    // typedef std::vector< CSeedData > SeedData;

    struct SearchThread;    ///< Executes jobs sequentially in a thread.
    struct SearchJob;       ///< Executes one job (one CSeedData instance).
    struct TaskGroup;       ///< Process seeds for one read.
    struct TaskDescriptor;  ///< Base of Task.
    struct Task;            /**< Process seeds for a combination of
                                 (readid, refid, fwd mate). */

public:

    // typedef Stat< StatParams::N_PARAMS > BatchStat;

    CBatch( CSearchContext & ctx, size_t batch_num = 0 );
    bool Run();
    bool RunSeeder();

    CSearchContext const & GetSearchCtx() const { return ctx_; }
    CSearchContext & GetSearchCtx() { return ctx_; }
    CReadData const & GetReads() const { return *reads_; }
    // SeedData & GetSeeds() { return seeds_; }
    // SeedData const & GetSeeds() const { return seeds_; }
    size_t GetJobSize() const { return reads_per_job_; }

    // BatchStat & GetStat() { return stat_; }

private:

    typedef std::vector< std::string * > OutStr;

    // static BatchStat::UpdateMap const BatchStatMap;
    static std::vector< std::string > const BatchStatDescriptions;

    void OutputThread();

    CSearchContext & ctx_;
    std::unique_ptr< CReadData > reads_;
    // SeedData seeds_;
    OutStr out_str_;
    std::mutex out_mtx_;
    std::condition_variable out_cvar_;
    // BatchStat stat_;
    size_t reads_per_job_;
    size_t batch_num_ = 0;
};

//------------------------------------------------------------------------------
// inline std::unique_ptr< CBatch > MakeBatch( CSearchContext & ctx )
inline std::unique_ptr< CBatch > MakeBatch(
        CSearchContext & ctx, size_t batch_num )
{
    /*
    return std::unique_ptr< CBatch >(
            new CBatch( ctx, (*ctx.stat )[MHStat::N_BATCHES] ) );
    */
    return std::unique_ptr< CBatch >( new CBatch( ctx, batch_num ) );
}

READFINDER_NS_END

#endif

