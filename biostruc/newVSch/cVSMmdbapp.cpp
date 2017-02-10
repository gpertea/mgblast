/*
 * $Id: cVSMmdbapp.cpp,v 1.1 2005/07/26 17:11:46 chenj Exp $
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
 * $Log: cVSMmdbapp.cpp,v $
 * Revision 1.1  2005/07/26 17:11:46  chenj
 * Making linux VSMmdb.cgi
 *
 *
 *
 */

#include <corelib/ncbiexec.hpp>
#include <util/util_exception.hpp>
#include <serial/exception.hpp>

#include "hVSMmdbapp.hpp"
#include "hVSMmdbcmd.hpp"
#include "VastSrchUti.hpp"
#include "mmdbuti.hpp"

#undef _DEBUG
#define _DEBUG

USING_NCBI_SCOPE;

static CVSMmdbApp theVSMmdbApp;
ofstream 		LogStr;	


// Initialization

string VSMmdbConf::BaseUrl;
string VSMmdbConf::CgiName;
string VSMmdbConf::HeadFName;
string VSMmdbConf::HtmlDir;
string VSMmdbConf::VSNbrUrl;
string VSMmdbConf::VSNbrCgi;
string VSMmdbConf::HomePage;
string VSMmdbConf::QMqueStrName;
string VSMmdbConf::QMqueNbrName;
string VSMmdbConf::MonitUrl;
string VSMmdbConf::MonitCgi;
string VSMmdbConf::IbutHelpF;


int main(int argc, const char* argv[])
{


// string curtime = CTime(CTime::eCurrent).AsString();
// fprintf(stderr, "Main begins at curtime %s\n", curtime.c_str());

  CCgiResponse	res;

  try {

    // Execute main application function
    	return theVSMmdbApp.AppMain(argc, argv);
  }
  catch (CUtilException&) {
	cout << "VSMmdbApp.AppMain failed: CUtilException\n";
  }
  catch (CSerialException&) {
	cout << "VSMmdbApp.AppMain failed CSerialException\n";
  }

}


BEGIN_NCBI_SCOPE

void CVSMmdbApp::Init()
{

   const CNcbiRegistry& reg = GetConfig();

   if (!VastSrchInitialize()) {

	PrtMes::PrintErrorHeader((CCgiResponse*)0, "VastSrch", 1);
	PrtMes::PrintMsgWithTail((CCgiResponse*)0, 
			"Can't Initialize VastSrch database");
	exit(1);
   }

   if (!OpenMMDBAPI(0, NULL)) {
         PrtMes::PrintErrorHeader((CCgiResponse*)0, "VastSrch", 1);
         PrtMes::PrintMsgWithTail((CCgiResponse*)0,
		"Unable to Open MMDB Api on Server.<p>\nPlease inform info@ncbi.nlm.nih.gov\n");

          exit(1);
   }

   // Read Conf. file

   static VSMmdbConf thisConf;
   static ImgConf    thisImgConf;

   thisConf.BaseUrl =      reg.Get("VASTSEARCH", "BaseUrl");
   thisConf.CgiName =      reg.Get("VASTSEARCH", "CgiName");
   thisConf.HeadFName =    reg.Get("VASTSEARCH", "HeadFName");
   thisConf.HtmlDir =      reg.Get("VASTSEARCH", "HtmlDir");
   thisConf.VSNbrUrl =     reg.Get("VASTSEARCH", "VSNbrUrl");
   thisConf.VSNbrCgi =     reg.Get("VASTSEARCH", "VSNbrCgi");
   thisConf.HomePage =     reg.Get("VASTSEARCH", "HomePage");
   thisConf.QMqueStrName = reg.Get("VASTSEARCH", "QMqueStrName");
   thisConf.QMqueNbrName = reg.Get("VASTSEARCH", "QMqueNbrName");
   thisConf.MonitUrl = 	   reg.Get("VASTSEARCH", "MonitUrl");
   thisConf.MonitCgi = 	   reg.Get("VASTSEARCH", "MonitCgi");
   thisConf.IbutHelpF =    reg.Get("VASTSEARCH", "IbutHelpF");

   thisImgConf.VastUrl =   thisConf.VSNbrUrl;
   thisImgConf.VastCgi =   thisConf.VSNbrCgi;
   thisImgConf.CddUrl =    reg.Get("VASTSEARCH", "CddUrl");
   thisImgConf.CddCgi =    reg.Get("VASTSEARCH", "CddCgi");
   thisImgConf.EntrezUrl = reg.Get("VASTSEARCH", "EntrezUrl");
   thisImgConf.EntrezCgi = reg.Get("VASTSEARCH", "EntrezCgi");

   if ( thisConf.BaseUrl.empty() || thisConf.CgiName.empty()
	|| thisConf.HeadFName.empty() || thisConf.HtmlDir.empty()
	|| thisConf.VSNbrUrl.empty() || thisConf.VSNbrCgi.empty()
	|| thisConf.HomePage.empty() || thisConf.QMqueStrName.empty()
	|| thisConf.QMqueNbrName.empty() || thisConf.MonitUrl.empty()
	|| thisConf.MonitCgi.empty() || thisConf.IbutHelpF.empty() )
   {
	 PrtMes::PrintErrorHeader((CCgiResponse*)0, "VastSrch", 1);
         PrtMes::PrintMsgWithTail((CCgiResponse*)0,
		"Missing Configuratio file or the Configuration file is imcomplete.\n");

   }



   CCgiApplication :: Init();

}



CNcbiResource* CVSMmdbApp::LoadResource(void)
{
    auto_ptr<CVSMmdbResource> resource (new CVSMmdbResource (GetConfig()));

    // add commands to the resource class

    resource->AddCommand (new CVSMmdbSubmitCommand(*resource));

    resource->AddCommand (new CVSMmdbStrTextCommand(*resource));
    resource->AddCommand (new CVSMmdbStrImgCommand(*resource));

    resource->AddCommand (new CVSMmdbViewCommand(*resource));

    resource->AddCommand (new CVSMmdbSearchCommand(*resource));

    resource->AddCommand (new CVSMmdbChkErrCommand(*resource));

    return resource.release();

}


int CVSMmdbApp::ProcessRequest(CCgiContext& ctx)
{

    // execute request
    ctx.GetResource().HandleRequest(ctx);
    return(0);
}


void CVSMmdbApp::Exit()
{
   VastSrchFinish();

   CCgiApplication :: Exit();
}


END_NCBI_SCOPE
