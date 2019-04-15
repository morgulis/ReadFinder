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

#include <libreadfinder/batch.hpp>
#include <libreadfinder/fast_seeds.hpp>
// #include <libmatchhits/mh_stat.hpp>
// #include <libmatchhits/search_thread.hpp>

#include <libtools/taskarray.hpp>

READFINDER_NS_BEGIN

/*
//==============================================================================
//------------------------------------------------------------------------------
CBatch::BatchStat::UpdateMap const CBatch::BatchStatMap = {
    MHStat::N_READS,
    MHStat::N_PAIRED_READS,
    MHStat::N_JOBS,
    MHStat::N_PARAMS,
    MHStat::N_SEEDER_HITS,
    MHStat::N_EXTENDED_SEEDS,
    MHStat::N_SEEDER_EEXONS,
    MHStat::N_SEEDER_EXONS,
    MHStat::N_PRIMARY_SEEDS,
    MHStat::N_MAPPED_READS,
};

//------------------------------------------------------------------------------
std::vector< std::string > const CBatch::BatchStatDescriptions = {
    MHSTAT_N_READS_DESCR,
    MHSTAT_N_PAIRED_READS_DESCR,
    MHSTAT_N_JOBS_DESCR,
    MHSTAT_MAX_JOB_READS_DESCR,
    MHSTAT_N_SEEDER_HITS_DESCR,
    MHSTAT_N_EXTENDED_SEEDS_DESCR,
    MHSTAT_N_SEEDER_EEXONS_DESCR,
    MHSTAT_N_SEEDER_EXONS_DESCR,
    MHSTAT_N_PRIMARY_SEEDS_DESCR,
    MHSTAT_N_MAPPED_READS_DESCR,
};
*/

//------------------------------------------------------------------------------
CBatch::CBatch( CSearchContext & ctx, size_t batch_num )
    : ctx_( ctx ),
      reads_( new CReadData( /*ctx_.memmgr_, */ctx_.logger_, *ctx_.seqs,
                             ctx_.batch,
                             ctx_.progress_flags_ ) ),
      batch_num_( batch_num )
{
    /*
    stat_[StatParams::N_READS] = reads_->GetNReads();
    stat_[StatParams::N_PAIRED_READS] = reads_->GetNPaired();
    */

    ctx_.n_reads += reads_->GetNReads();

    // estimate the number of sub-batches
    //
    {
        static size_t const CFACTOR = 100;

        size_t n_jobs( 0 );
        auto n_threads( ctx_.n_threads );
        auto n_reads( reads_->GetNReads() );

        if( n_reads > n_threads*CFACTOR )
        {
            reads_per_job_ = n_reads/(n_threads*CFACTOR) + 1;
            n_jobs = n_threads*CFACTOR;
        }
        else
        {
            reads_per_job_ = n_reads/n_threads + 1;
            n_jobs = n_threads;
        }

        out_str_.resize( n_jobs, nullptr );
        /*
        stat_[StatParams::N_JOBS] = n_jobs;
        stat_[StatParams::MAX_JOB_READS] = reads_per_job_;
        seeds_.resize( n_jobs );
        */
    }
}

//------------------------------------------------------------------------------
void CBatch::OutputThread()
{
    auto & os( ctx_.GetOutStream() );

    for( size_t start_job( 0 ), end_job( 0 ); end_job < out_str_.size(); )
    {
        {
            std::unique_lock< std::mutex > lock( out_mtx_ );
            for( ; end_job < out_str_.size() && out_str_[end_job] != nullptr;
                   ++end_job );

            if( end_job == start_job )
            {
                out_cvar_.wait_for( lock, std::chrono::seconds( 1 ) );
            }
        }

        for( ; start_job < end_job; ++start_job )
        {
            os << *out_str_[start_job];
            delete out_str_[start_job];
            out_str_[start_job] = nullptr;
        }
    }

    os << std::flush;
}

/*
//------------------------------------------------------------------------------
static std::string GetSeedsFileNamePrefix()
{ return std::string( "./seeds." ); }

//------------------------------------------------------------------------------
bool CBatch::Run()
{
    if( reads_->GetNReads() != 0 )
    {
        M_TRACE( "batch " << batch_num_ << " started" );
        auto batch_acc( stat_.Accumulate( *ctx_.stat, BatchStatMap ) );

        if( ctx_.load_seeds )
        {
            std::string fname( GetSeedsFileNamePrefix() +
                               std::to_string( batch_num_ ) );
            std::ifstream is( fname.c_str(), std::ios::binary );

            for( auto & seeds : seeds_ )
            {
                seeds.Load( is );
            }
        }
        else
        {
            CFastSeeds( *this ).Run();

            if( ctx_.save_seeds )
            {
                std::string fname(
                        GetSeedsFileNamePrefix() +
                        std::to_string( batch_num_ ) );
                std::ofstream os( fname.c_str(), std::ios::binary );

                for( auto const & seeds : seeds_ )
                {
                    seeds.Save( os );
                }

                return true;
            }
        }

        std::thread out_thread( [this](){ OutputThread(); } );

        if( !seeds_.empty() )
        {
            std::atomic< size_t > job_idx( 0 );
            CProgress p( "aligning", "tasks", ctx_.progress_flags_ );
            auto ph( p.GetTop().Split( seeds_.size() ) );
            CTaskArray< SearchThread > jobs(
                    ctx_.n_threads, *this, job_idx, ph );
            p.Start();
            jobs.Start();
            p.Stop();
        }

        out_thread.join();
        M_INFO( ctx_.logger_,
                "batch statistics:\n" <<
                stat_.Format( BatchStatDescriptions ) );
        M_TRACE( "batch " << batch_num_ << " finished" );
        return true;
    }

    return false;
}
*/

//------------------------------------------------------------------------------
bool CBatch::RunSeeder()
{
    if( reads_->GetNReads() != 0 )
    {
        // auto batch_acc( stat_.Accumulate( *ctx_.stat, BatchStatMap ) );
        CFastSeeds( *this, true ).Run();
        return true;
    }

    return false;
}

READFINDER_NS_END

