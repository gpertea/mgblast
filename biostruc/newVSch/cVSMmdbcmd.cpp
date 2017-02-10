/*
 * $Id: cVSMmdbcmd.cpp,v 1.1 2005/07/26 17:11:46 chenj Exp $
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
 * $Log: cVSMmdbcmd.cpp,v $
 * Revision 1.1  2005/07/26 17:11:46  chenj
 * Making linux VSMmdb.cgi
 *
 *
 *
 */


#include <corelib/ncbiexec.hpp>
#include <corelib/ncbifile.hpp>
#include <serial/objostr.hpp>

#include "hVSMmdbapp.hpp"
#include "hVSMmdbcmd.hpp"
#include "VastSrchUti.hpp"
#include "SHGlobal.hpp"
#include "mmdbuti.hpp"
#include "rdpdb.hpp"
#include "qmuti.hpp"
#include "SendSummary.hpp"

#include <sstream>

#include <objmime.h>

#include <string.h>
#include <sys/time.h>
#include <math.h>


// #include "../../newVSch/mmdbdep.h"
// #include "../../newVSch/mmdbapi.h"


// Initialization

string 		DataInfo::JobID;
string 		DataInfo::JobType;
string 		DataInfo::subset;
stringstream 	DataInfo::PDB_ios;
stringstream 	DataInfo::PDB_chk;
QmReqID 	DataInfo::ReqID = 0;
QmGrpID		DataInfo::GrpID = 0;


using namespace SHProjNS;
using namespace ncbi;

static ostringstream output;

CVSMmdbCommand::CVSMmdbCommand(CNcbiResource& resource) : CNcbiCommand(resource)
{};

string CVSMmdbCommand::GetEntry() const
{
	return string("cmdVSMmdb"); 
}



// CVSMmdbSubmitCommand

CVSMmdbSubmitCommand :: CVSMmdbSubmitCommand(CNcbiResource& resource) : CVSMmdbCommand(resource)
{};

CNcbiCommand* CVSMmdbSubmitCommand :: Clone(void) const
{
	return new CVSMmdbSubmitCommand(GetVSMmdbResource());
}

string CVSMmdbSubmitCommand :: GetName() const
{
	return string("Submit");
}


string CVSMmdbSubmitCommand :: ReadReqAndPDBSubmission(CCgiContext& ctx)
{

  string HeadValue = "HEADER";
  string Remark = "REMARK  1";
  string molname = "QUERY PROTEIN";

  CCgiRequest& req = ctx.GetRequest();
  CCgiResponse& resp = ctx.GetResponse();

  static DataInfo Dinfo;
  static VSMmdbConf thisConf;

  string date = GetTodaysDate();

  /* Write the contents of the pdb text area to buffer */

  string pdbdata = req.GetEntry("pdbfile").GetValue();
  if (pdbdata.empty()) {
      PrtMes::PrintErrorHeader(&resp, "VastSrch", 1);
      PrtMes::PrintMsgWithoutTail(&resp,
            "<h3>Missing Required PDB File:</h3><p>\nPlease submit a file.\n");

      
	resp.out() << "<HR align=left width=800><br>\n"
		<< "<a href=\"" << thisConf.HomePage
                << "\">Back</a> to the HomePage\n" 
		<< "</body>\n</html>";
	exit(1);
  }

  if (pdbdata.size() > PDBMAXLEN) {

      PrtMes::PrintErrorHeader(&resp, "VastSrch", 1);
      PrtMes::PrintMsgWithTail(&resp,
           "<h2>You cannot submit a file greater than 5MB!</h2><p>\nIf this is a problem, please notify info@ncbi.nlm.nih.gov\n");

             THROWS("input file too bigger");
  }

  unsigned pos = pdbdata.find(HeadValue);
  if (pos == string::npos) {

                Dinfo.PDB_ios.setf(ios::left);
                Dinfo.PDB_ios << setw(10) << HeadValue << setw(40) << molname;
                Dinfo.PDB_ios.unsetf(ios::left);
                Dinfo.PDB_ios << date;
                Dinfo.PDB_ios << setw(7) << Dinfo.JobID.substr(0,4)<< "\n";
                Dinfo.PDB_ios.setf(ios::left);
                Dinfo.PDB_ios << pdbdata;
  }
  else {
		if (pos) pdbdata = pdbdata.substr(pos);
		char data_line[100];
		istringstream  strm_tmp(pdbdata);

		strm_tmp.getline(data_line, 100);
		string strtmp = data_line;
		unsigned siz = strtmp.size();
		if (siz < 59)
		{
		     for (unsigned i= siz; i < 59; i++) strtmp += ' ';
		     strtmp.replace(50, 9, date);
		     Dinfo.PDB_ios << strtmp << "\n"
			    << pdbdata.substr(siz+1);
		}
		else {
	
                	pdbdata.replace(50, 9, date);
                	Dinfo.PDB_ios << pdbdata; 
		}
  }

  Dinfo.subset = req.GetEntry("dataset").GetValue();
  if (Dinfo.subset.empty()) {

        PrtMes::PrintErrorHeader(&resp, "VastSrch", 1);
        PrtMes::PrintMsgWithoutTail(&resp,
                "<h3>No search set provided. ");
        PrtMes::PrintMsgWithoutTail(&resp, Dinfo.subset);
        PrtMes::PrintMsgWithTail(&resp,
                "Please select 'All' or 'Nonredundant' search set\n");
        exit(1);
     }

  return date;

} // end of ReadReqAndPDBSubmission


   
void CVSMmdbSubmitCommand :: Execute(CCgiContext& ctx)
{

   	CCgiResponse& resp = ctx.GetResponse();
   	static DataInfo Dinfo;
   
   	int iJobId = GetJobID_DB();
   	Dinfo.JobID = "VS" + ToString(iJobId);
   
   	string date = ReadReqAndPDBSubmission(ctx);

	string headerStr = Dinfo.PDB_ios.str().substr(10, 39);
	if (headerStr.find_first_not_of(" ") == string::npos) 
		headerStr = "NULL";
        UploadReqInfo(Dinfo.subset, iJobId, date, headerStr);


        //if (!UploadAllPdbOnly(Dinfo.JobID, Dinfo.PDB_ios.str(), Dinfo.PDB_chk.str())) {
        if (!UploadOriPdbOnly(Dinfo.JobID, Dinfo.PDB_ios.str())) {

             THROWS("UploadAllPdb failed");
        }
// return;

	static VSMmdbConf thisConf;

        if (!sqmInit())  {

	     resp.WriteHeader();
	     resp.out() << "<html>\n<body>\n"
                 << "<h1>VastSch Error: </h1>\n<br><br>\n"
                 << "Can't Initialize QM\n"
                 << "</body>\n</html>\n\n";

            ERR_POST("Can't Initialize QM.");
            throw exception();
        }

   	Dinfo.ReqID = sqmSubmit(thisConf.QMqueStrName.c_str(), 0, 0);
	if (!Dinfo.ReqID) {


	      resp.WriteHeader();

              resp.out() << "<html>\n<body>\n"
                 << "<h1>VSchSrv Erro:</h1>\n<br><br>\n"
                 << "ReqId1 = 0\n"
                 << "</body>\n</html>\n\n";

	      ERR_POST("Got 0 as reqId from QM");
	      THROWS("Got 0 as reqId from QM");
//	      throw exception();

        }

// 	sqmSetAction(Dinfo.ReqID, qmReqActionRun);

/*
	char *IP = getenv("REMOTE_ADDR");
        if (!IP) {

                resp.WriteHeader();

                resp.out() << "<html>\n<body>\n"
                     << "<h1>VachSrv Error: </h1>\n<br><br>\n"
                     << "No IP information.\n"
                     << "</body>\n</html>\n";

                ERR_POST("No IP information.\n");
                throw exception();
        }

	sqmSetIP(Dinfo.ReqID, IP);
*/

	Dinfo.GrpID = sqmSetGroup(Dinfo.ReqID, 0);
	ChkGrpId(ctx, 0, Dinfo.GrpID);

	UploadQMIdsAndSetPdbLoaded(Dinfo.ReqID, Dinfo.GrpID, iJobId);

	if (!sqmFree()) {

        	resp.WriteHeader();
        	resp.out() << "<html>\n<body>\n"
                        << "<h1>VastSch Error: </h1>\n<br><br>\n"
                        << "Can't close QM\n"
                        << "</body>\n</html>\n\n";

                ERR_POST("Can't close QM.");
                throw exception();

    	};

	string UrlStr = thisConf.MonitUrl + thisConf.MonitCgi
			+ "?cmdSch=Status&ReqId=" + ReqId2Str(Dinfo.ReqID)
			+ "&GrpId=" + ReqId2Str(Dinfo.GrpID)
			+ "&QMQue=" + thisConf.QMqueStrName;

	resp.WriteHeader();
    	resp.out() << "<html>\n"
                << "<meta http-equiv=refresh content=0;url=\""
                << UrlStr
                << "\">\n</html>\n";


} 	// CVSMmdbSubmitCommand :: Execute()




// CVSMmdbStrTextCommand

CVSMmdbStrTextCommand :: CVSMmdbStrTextCommand(CNcbiResource& resource)
 : CVSMmdbCommand(resource)
{};

string CVSMmdbStrTextCommand :: GetName() const
{

	return string("StrText"); 	// cmdVSMmdb=StrText,
};


CNcbiCommand* CVSMmdbStrTextCommand :: Clone(void)  const
{
     return new CVSMmdbStrTextCommand(GetVSMmdbResource());
};


void CVSMmdbStrTextCommand :: Execute(CCgiContext& ctx)
{

     CCgiResponse& resp = ctx.GetResponse();
     CCgiRequest&  req = ctx.GetRequest();
     static DataInfo Dinfo;
     static VSMmdbConf thisConf;
     static ImgInfo Iinfo;
     string strtmp;

     string value = req.GetEntry("grpid").GetValue();
     Dinfo.GrpID = atoll(value.c_str());
     if (!Dinfo.GrpID) {
	
	output.str("");
	output << "<h3>Missing request_id the request_id provided is 0.</h3><p>\n"
		<< "Please input a request_id.\n";
	PrtMes::PrintErrorHeader(&resp, "VastSrch", 1);
        PrtMes::PrintMsgWithoutTail(&resp, output.str());

        resp.out() << "<HR align=left width=800><br>\n"
                << "<a href=\"" << thisConf.HomePage
                << "\">Back</a> to the HomePage\n"
                << "</body>\n</html>";
        exit(1);

     }
     Iinfo.id = GrpId2JobId(Dinfo.GrpID);
     if (!Iinfo.id) {
	
	output.str("");
	output << "<h3>Invalid request_id.</h3><p>\n"
	   << "You provided either an invalid request_id, "
	   << "or your request_id has expired since the search results are kept only for 10 days.<p>\n"
  	   << "Please check your request_id carefully, or you may resubmit your file to start a new search.\n";

	PrtMes::PrintErrorHeader(&resp, "VastSrch", 1);
	PrtMes::PrintMsgWithoutTail(&resp, output.str());

	resp.out() << "<HR align=left width=800><br>\n"
                << "<a href=\"" << thisConf.HomePage
                << "\">Back</a> to the HomePage\n"
                << "</body>\n</html>";
        exit(1);


     }
     Dinfo.JobID = "VS" + ToString(Iinfo.id);
       
     Dinfo.ReqID = JobId2ReqId(Iinfo.id);

     if (!VSMmdbJobDone(Iinfo.id)) {

	string errmsg, wmsg, wmsgH;
	DownloadErrWarnMsgs( (long long)Dinfo.GrpID, errmsg, wmsg, wmsgH);
	if ( errmsg.empty() && wmsg.empty() && wmsgH.empty() ) {
	
	// redirect  

		string UrlStr = thisConf.MonitUrl + thisConf.MonitCgi
                        + "?cmdSch=Status&ReqId=" + ReqId2Str(Dinfo.ReqID)
			+ "&GrpId=" + ReqId2Str(Dinfo.GrpID)
			+ "&QMQue=" + thisConf.QMqueStrName;

        	resp.WriteHeader();
        	resp.out() << "<html>\n"
                	<< "<meta http-equiv=refresh content=0;url=\""
                	<< UrlStr
                	<< "\">\n</html>\n";
	}
	else {
	
	     resp.WriteHeader();
     	     resp.out() << "<html>\n<head>\n"
                << "<title>VAST Search Error Report</title>\n";
     	     PrintFileData(resp, thisConf.HtmlDir + thisConf.HeadFName);

     	     resp.out() << "<pre><HR align=left width=800><br>\n";
     	     if  ( !errmsg.empty() ) {

        	resp.out()
                  << "<font class=H2>Errors in the submitted file:</font><br>\n"
                  << errmsg << "<br><br>\n\n";
     	     }

     	     if ( !wmsg.empty() ) {

        	resp.out() 
		  << "\n<font class=H2>Warning: the following chains may have no structure neighbors due to missing Ca atoms.</font><br>\n"
                  << wmsg << "<br><br>\n\n";
     	     }

     	     if ( !wmsgH.empty() ) {

        	resp.out() << "\n<font class=H2>Warning: the following HETATM records interpreted as modified residues, changed to ATOM records in the submitted file.</font><br>\n"
                	<< wmsgH << "<br><br>\n\n";
     	     }

	     resp.out() << "<a href=\"" << thisConf.HomePage
                << "\">Back</a> to the HomePage\n"; 

	}

        return;
     }


     Dinfo.subset = GetSubsetName4Job(Iinfo.id);
     if ( Dinfo.subset == "NULL" )  {
		ERR_POST("GetSubsetName4Job() failed");
		throw exception();
     }
     
     resp.WriteHeader();
     resp.out() << "<html>\n<head>\n"
                << "<title>VAST Search Structure Summary:&nbsp;"
                << GetTodaysDate()
                << "; ReqId: "
                << Dinfo.GrpID
                << "</title>\n";
  
     PrintFileData(resp, thisConf.HtmlDir + thisConf.HeadFName);

string curtime;
     value= req.GetEntry("ViewNbr").GetValue();
     if (value.empty()) {

	if (!VSNbrStarted(Iinfo.id)) Dinfo.JobType = "Biostr";
	else Dinfo.JobType = "Nbr";
     }
     else if (value == "No") Dinfo.JobType = "Biostr";
     else if (value=="Yes" || value=="Show") Dinfo.JobType = "Nbr";

curtime = CTime(CTime::eCurrent).AsString();
fprintf(stderr, "SendSummaryPageText curtime %s\n", curtime.c_str());

     if (Dinfo.JobType == "Biostr") {

      	SendSummaryPageText(ctx, 0, 0);

curtime = CTime(CTime::eCurrent).AsString();
fprintf(stderr, "after SendSummaryPageText curtime %s\n", curtime.c_str());

     }
     else if (Dinfo.JobType == "Nbr") {
	
	if (VSJobDone(Iinfo.id)) {
	    SendSummaryPageText(ctx, 1);
	}
	else {
	 
		SendSummaryPageText(ctx, 1, 0);
	  	exit(0);
	}
     }
     else {
	 THROWS("ViewNbr input incorrect");
     }

	
     ctx.GetResponse().out() << "</body>\n</html>\n\n";

     return;

};  // end of CVSMmdbStrCommand :: Execute()





// CVSMmdbStrImgCommand
CVSMmdbStrImgCommand :: CVSMmdbStrImgCommand(CNcbiResource& resource)
 : CVSMmdbCommand(resource)
{};


string CVSMmdbStrImgCommand :: GetName() const
{
        return string("StrImg"); 
};


CNcbiCommand* CVSMmdbStrImgCommand :: Clone(void)  const
{
     return new CVSMmdbStrImgCommand(GetVSMmdbResource());
};



void CVSMmdbStrImgCommand :: Execute(CCgiContext& ctx)
{
    CCgiRequest& req = ctx.GetRequest();
    static DataInfo Dinfo;
    static ImgInfo Iinfo;

    string value = req.GetEntry("grpid").GetValue();
    if (value.empty()) {

	ERR_POST("No grpid");
	throw exception();

    }
    Dinfo.GrpID = atoll(value.c_str());
    if (!Dinfo.GrpID) {
	ERR_POST("Grp ID is 0, from CVSMmdbStrImgCommand :: Execute");
	throw exception();
    }

    Iinfo.id = GrpId2JobId(Dinfo.GrpID);
    Dinfo.JobID = "VS" + ToString(Iinfo.id);

    value = req.GetEntry("chbeg").GetValue();
    if (value.empty()) THROWS("NO chbeg");
    Iinfo.chbeg = atoi(value.c_str());

    value = req.GetEntry("chend").GetValue();
    if (value.empty()) THROWS("NO chend");
    Iinfo.chend = atoi(value.c_str());

    value = req.GetEntry("imgsize").GetValue();
    if (value.empty()) THROWS("NO chend");
    Iinfo.imgsize = atoi(value.c_str());

/*
    value = req.GetEntry("ViewNbr").GetValue();  
    if (value == "Yes" || value=="Show") 
		DrawImg(ctx, 1);  // 1 or 0 VastLink
    else DrawImg(ctx, 0);
*/

    DrawImg(ctx);

} // CVSMmdbSubmitStrImgCommand :: Execute



// CVSMmdbViewCommand
	
CVSMmdbViewCommand :: CVSMmdbViewCommand(CNcbiResource& resource) : CVSMmdbCommand(resource)
{};

CNcbiCommand* CVSMmdbViewCommand :: Clone(void) const
{
     return new CVSMmdbViewCommand(GetVSMmdbResource());
}

string CVSMmdbViewCommand :: GetName() const
{
     return string("View3D");
}



void CVSMmdbViewCommand :: SendStructureMIME(char Filetype, int Mime, 
			int Complexity, int Models, CCgiContext& ctx)
{

  static DataInfo Dinfo;
  CCgiResponse& resp = ctx.GetResponse();

  int iJobId = atoi(Dinfo.JobID.substr(2).c_str());

  BiostrucPtr bsp;
  SeqEntryPtr sep;
  NcbiMimeAsn1Ptr mime;
  BiostrucSeqPtr bssp;

  /* Save or View  = add appropriate MIME header */

// bsp = FetchBS((char*)bFile.c_str(), 1, Complexity, Models, POWER_VIEW);

   bsp = VSOpenBSP(iJobId, Complexity, Models);
   if (bsp == NULL) {

   	string strtmp = "Missing Biostruc data for reqid = "
			+ ReqId2Str(Dinfo.GrpID) 
			+ ". Please inform chenj@ncbi.nlm.nih.gov for this problem. Thanks!";
        PrtMes::PrintErrorHeader(&resp, "VastSrch", 1);
        PrtMes::PrintMsgWithTail(&resp,  strtmp);

        exit(1);
   }

   bssp = BiostrucSeqNew();
   bssp->structure = bsp;

   sep = GetSeqEntryForJobId(iJobId, 0, NULL);

   if (sep == NULL) {

	string strtmp = "Missing SeqEntry data for reqid = "
                        + ReqId2Str(Dinfo.GrpID)
                        + ". Please inform chenj@ncbi.nlm.nih.gov for this problem. Thanks!";

        PrtMes::PrintErrorHeader(&resp, "VastSrch", 1);
	PrtMes::PrintMsgWithTail(&resp,  strtmp);

        exit(1);
   }

   if (sep->choice == 2) {  // BioseqSet, multiple chains

       BioseqSetPtr seps = (BioseqSetPtr) sep->data.ptrvalue;
       SeqEntryPtr subsep = (SeqEntryPtr) seps->seq_set;
       SeqIdPtr sid;
       int ori_gi=0;

       while (subsep != NULL) {

          sid = (SeqIdPtr) ((BioseqPtr)subsep->data.ptrvalue)->id;
          while (sid != NULL) {

                if (sid->choice == SEQID_GI) {

                    sid->data.intvalue = ori_gi++;
                    break;
                }
                sid=sid->next;
          }

          subsep = subsep->next;
       }

   }

   ValNodeLink(&(bssp->sequences), sep);

   mime = (NcbiMimeAsn1Ptr) ValNodeNew(NULL);    /* yanli */
   mime->choice = NcbiMimeAsn1_strucseq;
   mime->data.ptrvalue = bssp;

   /* the following headers are format-independent */

   ESerialDataFormat dataformat;
  // if (Filetype == 'i') {   BiostrucAsnWrite()
   if (Filetype == 'j') {
        /* Cn3D asn.1 format */
        if (Mime == LAUNCH_VIEWER || Mime == SAVE_FILE)
                         dataformat = eSerial_AsnBinary;
        else dataformat = eSerial_AsnText;

        auto_ptr <CObjectOStream> oos
                (CObjectOStream::Open(dataformat, resp.out()));
        CObjectOStream::AsnIo aip(*oos, "Ncbi-mime-asn1");

        switch (Mime) {
           case LAUNCH_VIEWER:
                resp.SetContentType("chemical/ncbi-asn1-binary");
                resp.WriteHeader();
                NcbiMimeAsn1AsnWrite(mime, aip, NULL);
                break;

           case SAVE_FILE:
                resp.SetContentType("application/octet-stream");
                resp.WriteHeader();
                NcbiMimeAsn1AsnWrite(mime, aip, NULL);
                break;

           case SEE_FILE:
                resp.SetContentType("text/html");
                resp.WriteHeader();
                resp.out() << "<HTML><PRE>\r\n";
                NcbiMimeAsn1AsnWrite(mime, aip, NULL);
                break;
            }

        aip.End();
   }
   if (Mime == SEE_FILE)
                resp.out() << "</PRE></HTML>\r\n";

//   CExec::System((string("rm -r ") + dir).c_str());

} // end of SendStructureMIME()



void CVSMmdbViewCommand :: Execute(CCgiContext& ctx)
{

fprintf(stderr, "CVSMmdbViewCommand :: Execute()\n");


    CCgiRequest& req = ctx.GetRequest();
    CCgiResponse& resp = ctx.GetResponse();
    static DataInfo Dinfo;

    string value = req.GetEntry("grpid").GetValue();
    if (value.empty()) {
	THROWS("No GrpId");
    }
    Dinfo.GrpID = atoll(value.c_str());

    Dinfo.JobID = "VS" + ToString(GrpId2JobId(Dinfo.GrpID));
    if (Dinfo.JobID.empty()) {
		ERR_POST("No JobID");
		throw exception();
    }

    string View = req.GetEntry("dopt").GetValue();   
    if (View.empty()) {

      	PrtMes::PrintErrorHeader(&resp, "VastSrch", 1);
      	PrtMes::PrintMsgWithoutTail(&resp,
			"		Wrong Display Option Type: dopt= ");
      	PrtMes::PrintMsgWithoutTail(&resp, View);
      	PrtMes::PrintMsgWithTail(&resp, ".<p>\n<pre>dopt supports: 'i' or 'a' ASN.1</pre><p>\n");

      	ERR_POST("No dopt provided");
	throw exception();
      
    }
   
    char DispOpt = '\0'; 
    int  Save = 0;
    switch (View[0])
    {
      case 'i':
      	DispOpt = 'i';
      	break;
      case 'a':
      	DispOpt = 'a';
      	Save = 1;
      	break;
      case 'j':
	DispOpt = 'j';
	break;
    }

    string com_value = req.GetEntry("Complexity").GetValue();
    if (com_value.empty()) {
	PrtMes::PrintErrorHeader(&resp, "VastSrch", 1);
	PrtMes::PrintMsgWithTail(&resp, "Can't find value of Complexity.\n");
	ERR_POST("No complexity provided");
	throw exception();
    }

    int iCount = -1;
    vector<string> ComplexityDescriptions;
    ComplexityDescriptions.push_back("Virtual Bond Model");
    ComplexityDescriptions.push_back("All Atom Model");
    ComplexityDescriptions.push_back("Up to 5 Models");
    ComplexityDescriptions.push_back("Up to 10 Models");
    ComplexityDescriptions.push_back("All Models");
    ComplexityDescriptions.push_back("Cn3D Subset");

    for (int i=0; i< (int)ComplexityDescriptions.size(); i++) 
	if (com_value.find(ComplexityDescriptions[i]) !=string::npos) {
		iCount = i;
		break;
        }

    if (iCount < 0) {

      string strtmp = "Wrong Complexity Type\nPlease select one from: "
			+ ComplexityDescriptions[5] + ", "
			+ ComplexityDescriptions[0] + ", and "
			+ ComplexityDescriptions[4] + ".";
      PrtMes::PrintErrorHeader(&resp, "VastSrch", 1);
      PrtMes::PrintMsgWithTail(&resp, strtmp);

      THROWS("Wrong Complexity type");
   }

    int MaxModels = 1;
    int Complex = ALLMDL;
    switch (iCount)
    {
      case 0 :
        Complex = ONECOORDRES;
        break;
      case 1:
        Complex = ONECOORDATOM;
        break;
      case 2:
        MaxModels = 5;
        break;
      case 3 :
        MaxModels = 10;
        break;
      case 4 :
        MaxModels = INT2_MAX;
        break;
      case 5:
        Complex = ONECOORDATOM;
        break;
      }

    string SaveChoice = req.GetEntry("save").GetValue();

    if (SaveChoice.find("See") != string::npos)
      	Save = 1;  /* sends MIMED-type files to view in ascii*/
    else if (SaveChoice.find("Save") != string::npos)
      	Save = 2;  /* sends MIMED-type files in raw save form */

    SendStructureMIME(DispOpt, Save, Complex, MaxModels, ctx);


} // end of CVSMmdbViewCommand :: Execute()





// CVSMmdbSearchCommand

CVSMmdbSearchCommand::CVSMmdbSearchCommand(CNcbiResource& resource):CVSMmdbCommand(resource)
{};


CNcbiCommand* CVSMmdbSearchCommand :: Clone(void) const
{
     return new CVSMmdbSearchCommand(GetVSMmdbResource());
}

string CVSMmdbSearchCommand :: GetName() const
{
     return string("Search");
}

void CVSMmdbSearchCommand :: Execute(CCgiContext& ctx)
{
    static DataInfo Dinfo;
    static VSMmdbConf thisConf;
    CCgiRequest& req = ctx.GetRequest();
    CCgiResponse& resp= ctx.GetResponse();


    string value = req.GetEntry("grpid").GetValue();
    if (value.empty())
		ERR_POST("No reqid in CVSMmdbSearchCommand :: Execute");
    Dinfo.GrpID = atoll(value.c_str());

    if ( SetupDoitInDB(Dinfo.GrpID) ) {
	
    	if (!sqmInit())  {

		resp.WriteHeader();

  		resp.out()  << "<html>\n<body>\n"
             		<< "<h1>VastSch Error: </h1>\n<br><br>\n"
	     		<< "Can't Initialize QM\n"
             		<< "</body>\n</html>\n\n";

            	ERR_POST("Can't Initialize QM.");
            	throw exception();
    	}
    
    	Dinfo.ReqID = sqmSubmit(thisConf.QMqueNbrName.c_str(), 0, 0);
    	if (!Dinfo.ReqID) {

	 	resp.WriteHeader();

              	resp.out() << "<html>\n<body>\n"
                 	<< "<h1>VSchSrv Erro:</h1>\n<br><br>\n"
                 	<< "ReqId2 = 0\n"
                 	<< "</body>\n</html>\n\n";

              	ERR_POST("Got 0 as reqId from QM");
              	throw exception();

	}

//    sqmSetAction(Dinfo.ReqID, qmReqActionRun);

    	QmReqID grpid = sqmSetGroup(Dinfo.ReqID, Dinfo.GrpID);
    	ChkGrpId(ctx, Dinfo.GrpID, grpid);

    	UpdateReqIdOfGrpId(Dinfo.GrpID, Dinfo.ReqID);

/*
    	char *IP = getenv("REMOTE_ADDR");
    	if (!IP) {

		resp.WriteHeader();

                resp.out() << "<html>\n<body>\n"
                     << "<h1>VachSrv Error: </h1>\n<br><br>\n"
                     << "No IP information.\n"
                     << "</body>\n</html>\n";

                ERR_POST("No IP information.\n");
                throw exception();
    	}

    	sqmSetIP(Dinfo.ReqID, IP);
*/


    	if (!sqmFree()) {

		resp.WriteHeader();
		resp.out() << "<html>\n<body>\n"
                        << "<h1>VastSch Error: </h1>\n<br><br>\n"
                        << "Can't close QM\n"
                        << "</body>\n</html>\n\n";

                ERR_POST("Can't close QM.");
                throw exception();

    	};

    }

/*

    output << thisConf.MonitUrl << thisConf.MonitCgi
		<< "?cmdSch=Status&ReqId=" << Dinfo.ReqID
		<< "&GrpId=" << Dinfo.GrpID
 		<< "&QMQue=" << thisConf.QMqueNbrName;
*/

    output.str("");
    output << thisConf.BaseUrl << thisConf.CgiName
		<< "?cmdVSMmdb=StrText&grpid="  << Dinfo.GrpID
		<< "&ViewNbr=Yes";

    resp.WriteHeader();
    resp.out() << "<html>\n"
                << "<meta http-equiv=refresh content=0;url=\""
                << output.str()
                << "\">\n</html>\n";



} // CVSMmdbSearchCommand :: Execute()



// CVSMmdbChkErrCommand

CVSMmdbChkErrCommand :: CVSMmdbChkErrCommand(CNcbiResource& resource):CVSMmdbCommand(resource)
{};

CNcbiCommand* CVSMmdbChkErrCommand :: Clone(void) const
{
     return new CVSMmdbChkErrCommand(GetVSMmdbResource());
}

string CVSMmdbChkErrCommand :: GetName() const
{
     return string("ChkErr");
}

void CVSMmdbChkErrCommand :: Execute(CCgiContext& ctx)
{

     CCgiResponse& resp = ctx.GetResponse();
     CCgiRequest&  req = ctx.GetRequest();
     static DataInfo Dinfo;
     static VSMmdbConf thisConf;

     string errmsg, wmsg, wmsgH;
     string value = req.GetEntry("Jobid").GetValue();
     if ( !value.empty() )
	DownloadErrWarnMsgs((unsigned)atoi(value.c_str()), errmsg, wmsg, wmsgH);
     else {
	
	value = req.GetEntry("grpid").GetValue();
	if ( !value.empty() )
	    DownloadErrWarnMsgs((long long)atoll(value.c_str()), errmsg, 
								wmsg, wmsgH);
     	else throw exception();
     }

     resp.WriteHeader();
     resp.out() << "<html>\n<head>\n"
		<< "<title>VAST Search Error Report</title>\n";
     PrintFileData(resp, thisConf.HtmlDir + thisConf.HeadFName);

     resp.out() << "<pre><HR align=left width=800><br>\n";
     if  ( !errmsg.empty() ) {

	resp.out()
		<< "<font class=H2>Errors in the submitted file:</font><br>\n"
		<< errmsg << "<br><br>\n\n";
     }

     if ( !wmsg.empty() ) {
	
	resp.out() << "\n<font class=H2>Warning: the following chains may have no structure neighbors due to missing Ca atoms.</font><br>\n"
		<< wmsg << "<br><br>\n\n";
     }

     if ( !wmsgH.empty() ) {

        resp.out() << "\n<font class=H2>Warning: the following HETATM records interpreted as modified residues, changed to ATOM records in the submitted file.</font><br>\n"
                << wmsgH << "<br><br>\n\n";
     }
     
     resp.out() << "<HR align=left width=800><br>\n"
	 	<< "<a href=# onclick=\"self.close();return false;\">"
		<< "Close Window</a>"
		<< "<br>\n</pre>\n</body>\n</html>\n\n";

/*
		<< "<a href=\"" << thisConf.HomePage
		<< "\">Back</a> to the HomePage\n" 
*/

}

