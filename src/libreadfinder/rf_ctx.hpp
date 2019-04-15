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

#ifndef LIBREADFINDER_MH_CTX_HPP
#define LIBREADFINDER_MH_CTX_HPP

#include <functional>
#include <memory>

#include <libreadfinder/rf_options.hpp>
// #include <libreadfinder/mh_stat.hpp>
#include <libreadfinder/refdata.hpp>
#include <libreadfinder/defs.hpp>

#include <libtools/logger.hpp>
#include <libtools/progress.hpp>
// #include <libtools/stat.hpp>

READFINDER_NS_BEGIN

//==============================================================================
struct CCommonContext
{
    CCommonContext( CommonOptions const & opts );
    CLogger logger_;
    // std::shared_ptr< CMemMonitor > memmgr_;
    int progress_flags_ = 0;
};

typedef std::shared_ptr< CCommonContext > CommonCtxP;

//==============================================================================
struct CSearchContext : public CCommonContext, public CSearchOptions
{
private:

    typedef std::unique_ptr<
        std::ostream,
        std::function< void( std::ostream * ) > > OutStream;

    // typedef TOOLS_NS::Stat< MHStat::N_PARAMS > SearchStat;

public:

    CSearchContext( CSearchOptions const & opts/* , bool init_output = true */ );

    void Check();
    std::ostream & GetOutStream() { return *osp; }

    void ResetOutStream( std::string const & name )
    {
        OutStream( new std::ofstream( name.c_str() ),
                   []( std::ostream * os ){ delete os; } ).swap( osp );
    }

    void ResetOutStream( std::ostream & os = std::cout )
    {
        OutStream( &os, []( std::ostream * osp ){} ).swap( osp );
    }

    /*
    TScore GetScoreThreshold( TScore s ) const
    {
        return std::min( penalties.top_cutoff*s, 0.0f );
    }
    */

    std::unique_ptr< CSeqInput > seqs;
    std::unique_ptr< CRefData > refs;
    // std::shared_ptr< SearchStat > stat;
    // bool run_second_pass = false;
    size_t n_reads = 0,
           n_mapped_reads = 0;

private:

    OutStream osp;
};

READFINDER_NS_END

#endif

