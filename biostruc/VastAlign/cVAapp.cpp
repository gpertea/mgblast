/*
 * $Id: cVAapp.cpp,v 1.1 2005/07/26 16:25:22 chenj Exp $
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
 *
 * $Log: cVAapp.cpp,v $
 * Revision 1.1  2005/07/26 16:25:22  chenj
 * source for vastalign.cgi
 *
 *
 */


#include <corelib/ncbiexec.hpp>
#include <util/util_exception.hpp>
#include <serial/exception.hpp>

#include "hVAapp.hpp"
#include "hVAcmd.hpp"
#include "hUtilib.hpp"
#include "PubVastApi.hpp"
#include "SHGlobal.hpp"

#undef _DEBUG
#define _DEBUG

USING_NCBI_SCOPE;

static CVAApp theVAApp;

int main(int argc, const char* argv[])
{

  CCgiResponse	res;

  try {

    // Execute main application function
    	return theVAApp.AppMain(argc, argv, 0, eDS_Default, 0);
  }
  catch (CUtilException&) {
	cout << "VAApp.AppMain failed: CUtilException\n";
  }
  catch (CSerialException&) {
	cout << "VAApp.AppMain failed CSerialException\n";
  }

}


BEGIN_NCBI_SCOPE

void CVAApp::Init()
{

   if ( !SHDBApiInitialize(true,"vastalign") ) {

	PrtMes::PrintErrorHeader((CCgiResponse*)0, "VastAlign", 1);
        PrtMes::PrintMsgWithTail((CCgiResponse*)0,
                        "Can't Initialize SHDBApi database");
	return;
   }

   if (!VastSrvInitialize()) {

	PrtMes::PrintErrorHeader((CCgiResponse*)0, "VastAlign", 1);
	PrtMes::PrintMsgWithTail((CCgiResponse*)0, 
			"Can't Initialize VastSrv database");
	exit(1);
   }

   CCgiApplication :: Init();

}



CNcbiResource* CVAApp::LoadResource(void)
{
    auto_ptr<CVAResource> resource (new CVAResource (GetConfig()));

    // add commands to the resource class

    resource->AddCommand (new CVACommand(*resource));

    return resource.release();

}


int CVAApp::ProcessRequest(CCgiContext& ctx)
{

    // execute request
    ctx.GetResource().HandleRequest(ctx);
    return(0);
}


void CVAApp::Exit()
{
   VastSrvFinish();

   CCgiApplication :: Exit();
}


END_NCBI_SCOPE
