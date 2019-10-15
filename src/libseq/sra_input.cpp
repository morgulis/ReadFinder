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

/** \file libseq/sra_input.hpp
    \brief Getting reads from SRA.
*/

#include <libseq/sra_input.hpp>

#ifdef USE_NGS

#include <ncbi-vdb/NGS.hpp>

SEQ_NS_BEGIN

//------------------------------------------------------------------------------
CSRAInput::CSRAInput( std::string const & acc, size_t start, ssize_t n_seq )
    : acc_( acc ), start_( ++start ),
      run_( ncbi::NGS::openReadCollection( acc ) ),
      it_( run_.getReadRange( start_, n_seq, Read::all ) )
{
}

//------------------------------------------------------------------------------
bool CSRAInput::Next()
{
    while( true )
    {
        if( !it_.nextRead() )
        {
            return false;
        }

        sd_.SetId( acc_ + "." + std::to_string( start_++ ) );

        try
        {
            if( it_.getNumFragments() == 2 )
            {
                sd_.SetNCols( 2 );

                {
                    it_.nextFragment();
                    StringRef sr( it_.getFragmentBases() );
                    sd_.SetData( 0, std::string( sr.data(), sr.size() ) );
                }

                {
                    it_.nextFragment();
                    StringRef sr( it_.getFragmentBases() );
                    sd_.SetData( 1, std::string( sr.data(), sr.size() ) );
                }
            }
            else
            {
                sd_.SetNCols( 1 );
                it_.nextFragment();
                StringRef sr( it_.getFragmentBases() );
                sd_.SetData( 0, std::string( sr.data(), sr.size() ) );
            }
        }
        catch( std::exception const & e )
        {
            ++n_failed_reads_;
            continue;
        }
        catch( ... )
        {
            throw;
        }

        break;
    }

    return true;
}

SEQ_NS_END

#endif

