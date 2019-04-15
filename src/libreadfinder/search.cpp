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

#include <config.h>

#include <libreadfinder/batch.hpp>
#include <libreadfinder/search.hpp>
#include <libreadfinder/defs.hpp>

READFINDER_NS_BEGIN

//==============================================================================
/*
//------------------------------------------------------------------------------
void Search( CSearchOptions const & opts )
{
    CSearchContext ctx( opts );

    while( MakeBatch( ctx )->Run() )
    {
        ++(*ctx.stat)[MHStat::N_BATCHES];
    }

    M_INFO( ctx.logger_,
            "search statistics:\n" <<
            ctx.stat->Format( MHStat::Descriptions ) );
    M_FORCE_LOG( ctx.logger_, "matchhits search: complete" );
}
*/

//------------------------------------------------------------------------------
void SearchSeeds( CSearchOptions const & opts )
{
    CSearchContext ctx( opts/* , false */ );
    size_t batch_num( 0 );
    while( MakeBatch( ctx, batch_num++ )->RunSeeder() );
    // while( MakeBatch( ctx )->RunSeeder() );
    M_FORCE_LOG( ctx.logger_, "total reads: " << ctx.n_reads );
                 // "total reads: " << ctx.stat->at( MHStat::N_READS ) );
    M_FORCE_LOG( ctx.logger_,
                 "failed reads: " << ctx.seqs->GetNumFailedReads() );
    M_FORCE_LOG( ctx.logger_, "mapped reads: " << ctx.n_mapped_reads );
                 // "mapped reads: " << ctx.stat->at( MHStat::N_MAPPED_READS ) );
    M_FORCE_LOG( ctx.logger_, "matchhits seeder: complete" );
}

READFINDER_NS_END

