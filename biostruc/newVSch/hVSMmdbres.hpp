#ifndef _VSMmdbRES_HPP
#define _VSMmdbRES_HPP

/*
 * $Id: hVSMmdbres.hpp,v 1.1 2005/07/26 17:11:46 chenj Exp $
 *
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *            National Center for Biotechnology Information (NCBI)
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government do not place any restriction on its use or reproduction.
 *  We would, however, appreciate having the NCBI and the author cited in
 *  any work or product based on this material
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 * ===========================================================================
 *
 * Author: Jie Chen
 * Source file for VSMmdb.cgi
 *
 * $Log: hVSMmdbres.hpp,v $
 * Revision 1.1  2005/07/26 17:11:46  chenj
 * Making linux VSMmdb.cgi
 *
 *
 *
 */



#include <cgi/ncbires.hpp>

using namespace ncbi;

BEGIN_NCBI_SCOPE

// class CVSMmdbResource

class CVSMmdbResource : public CNcbiResource
{
 public:

	CVSMmdbResource (CNcbiRegistry& config);
	virtual ~CVSMmdbResource() {};

	// define the command to be executed when no other command matches

	virtual CNcbiCommand* GetDefaultCommand(void) const;

// 	virtual void HandleRequest(CCgiContext& ctx);

};

END_NCBI_SCOPE

#endif
