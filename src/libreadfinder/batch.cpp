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

#include <libtools/taskarray.hpp>

READFINDER_NS_BEGIN

//------------------------------------------------------------------------------
CBatch::CBatch( CSearchContext & ctx, size_t batch_num )
    : ctx_( ctx ),
      reads_( new CReadData( /*ctx_.memmgr_, */ctx_.logger_, *ctx_.seqs,
                             ctx_.batch,
                             ctx_.progress_flags_ ) ),
      batch_num_( batch_num )
{
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

//------------------------------------------------------------------------------
bool CBatch::RunSeeder()
{
    if( reads_->GetNReads() != 0 )
    {
        CFastSeeds( *this, true ).Run();
        return true;
    }

    return false;
}

READFINDER_NS_END

