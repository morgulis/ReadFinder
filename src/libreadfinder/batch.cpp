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
      reads_( new CReadData( ctx_.logger_, *ctx_.seqs,
                             (ctx_.pre_screen ? &ctx.ws : nullptr),
                             *ctx_.fsidx,
                             ctx_.batch,
                             (ctx_.used_mem_bytes < ctx_.max_mem_bytes ?
                                ctx_.max_mem_bytes - ctx_.used_mem_bytes : 0),
                             ctx_.n_threads, ctx_.progress_flags_ ) ),
      batch_num_( batch_num )
{
    ctx_.n_reads += reads_->GetNReads();

    // estimate the number of sub-batches
    //
    {
        static size_t const CFACTOR = 100;

        auto n_threads( ctx_.n_threads );
        auto n_reads( reads_->GetNReads() );

        if( n_reads > n_threads*CFACTOR )
        {
            reads_per_job_ = n_reads/(n_threads*CFACTOR) + 1;
        }
        else
        {
            reads_per_job_ = n_reads/n_threads + 1;
        }
    }
}

//------------------------------------------------------------------------------
bool CBatch::RunSeeder()
{
    if( reads_->GetNReads() == 0 ) return false;

    for( reads_->Update();
            reads_->GetStartOId() < reads_->GetNReads(); reads_->Update() )
    {
        CFastSeeds( *this ).Run();
    }

    return true;
}

READFINDER_NS_END

