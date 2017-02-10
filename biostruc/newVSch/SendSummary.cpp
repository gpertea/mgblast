/*
 * $Id: SendSummary.cpp,v 1.1 2005/07/26 17:11:46 chenj Exp $
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
 * $Log: SendSummary.cpp,v $
 * Revision 1.1  2005/07/26 17:11:46  chenj
 * Making linux VSMmdb.cgi
 *
 *
 *
 */

#include "SendSummary.hpp"
#include "SHGlobal.hpp"

static string tmp_dir = "/tmp/";

using namespace ncbi;
using namespace SHProjNS;


static void CreateIbutn(string& urlstr, string& dispOp, 
		unsigned refT, CCgiContext& ctx)
{
     static VSMmdbConf thisConf;

     ctx.GetResponse().out()
	<< "<a href=# onclick=\"PopWin(" << urlstr
       	<< ", " << dispOp << ", " << refT*1000 << ");\">\n"
       	<<  "<img src=\"" << thisConf.HtmlDir << "infosmall.gif\" border=0>"
        << "</a>\n";
	
}

static void PrintActionButtons(CCgiContext& ctx)
{

  CCgiResponse& resp = ctx.GetResponse();
  static DataInfo 	Dinfo;
  static VSMmdbConf 	thisConf;
  static ImgInfo	Iinfo;
  string urlstr, dispOp;

  string altstr= "Click to view 3D structure";
  resp.out() <<"<Form method=get action="
	<< thisConf.BaseUrl << thisConf.CgiName << "/VS.c3d>\n"
	<< "<input name=grpid type=hidden value=" << Dinfo.GrpID << ">\n"
	<< "  <tr>    <td><br></td>  </tr>\n  <tr>\n    <td>\n"
	<< "      <table>\n        <tr>\n          <td>\n"
	<< "<input type=hidden name=cmdVSMmdb value=View3D>\n"
	<<"<input type=image name=cmdVSMmdb value=View3D alt=\"" << altstr
	<< "\" title=\"" << altstr << "\" src=" 
	<< thisConf.HtmlDir << "v3ds.gif>\n"
	<< "          </td>\n\n          <td class=H4 valign=TOP>\n"
	<< "&nbsp;&nbsp;of&nbsp;\n"
	<< "<select name=Complexity>\n"
	<< "<option value=\"Cn3D Subset\">All Atom Model\n"
	<<"<option value=\"Virtual Bond Model\">Backbone only Model"
//	<< "\n<option value=\"All Models\">All Models\n"
	<< "</select>\n"
	<< "&nbsp;with \n"
	<< "<a href=/Structure/CN3D/cn3d.shtml>Cn3D</a>"
	<< "<img src=\"/Structure/new.gif\" alt=New title=New"
	<< " valign=top>\n"
	<< "<input type=hidden name=dopt value=j>\n"
	<< "      <SELECT NAME=save>\n"
	<< "    	<OPTION VALUE=Viewer>Display\n"
	<< "    	<OPTION VALUE=See>See File\n"
	<< "    	<OPTION VALUE=Save>Save File\n"
	<< "      </SELECT>\n";

  urlstr = "'" + thisConf.IbutHelpF + "_Cn3D.html'"; 
  dispOp = "'resizable=yes, scrollbars=yes, width=420, height=200'";
  CreateIbutn(urlstr, dispOp, 15, ctx);

  resp.out()	<< "          </td>\n        </tr>\n      </table>\n"
	<< "    </td>\n  </tr>\n\n"
	<< "</Form>\n\n";


  if ( Dinfo.JobType == "Biostr" ) {

     resp.out() << "  <tr><td width=800><br><HR></td></tr>\n"
		<< "  <tr><td><br>\n";


     if ( !Iinfo.prot_ch_cnt) return;
     string errmsg, wmsg, wmsgH;
     DownloadErrWarnMsgs((long long)Dinfo.GrpID, errmsg, wmsg, wmsgH);
     if ( !wmsg.empty()  || !wmsgH.empty() ) {
	resp.out() << "<a href=# onclick=\"PopWin('" 
		<< thisConf.BaseUrl  << thisConf.CgiName
		<< "?cmdVSMmdb=ChkErr&grpid=" << Dinfo.GrpID << "', "
	        << "'resizable=yes, scrollbars=yes, width=800, height=500', "
		<< 300*1000 << ");\">\n"
		<< "<font class=TEXT color=red>Warning message from data parsing</font></a>";
     }

     resp.out() << "    </td>\n  </tr>\n\n";

     altstr = "Click to launch VAST neighbor searching";
     resp.out()<< "  <tr>  <td><br></td>  </tr>\n\n"
	<< "  <tr>\n    <td class=SMALL1>"
	<<"Click the \"Start\" button to launch the VAST search.</td></tr>\n\n";

     resp.out() << "<Form method=get action=\"" 
	<< thisConf.BaseUrl << thisConf.CgiName << "\">\n"
	<< "<input type=hidden name=grpid value=" << Dinfo.GrpID <<">\n"
        << "  <tr>\n    <td>\n      <table>\n        <tr>\n          <td>\n"
	<< "<input type=hidden name=cmdVSMmdb value=Search>\n"
        << "<input type=Image name=cmdVSMmdb value=Search alt=\"" << altstr
	<< "\" title=\"" << altstr << "\" src=" 
	<< thisConf.HtmlDir << "start2.gif>\n"
	<< "          </td>\n\n          <td class=H4 valign=TOP>\n&nbsp;&nbsp;"
        << "the VAST Calculation</font>\n"
	<< "          </td>\n        </tr>\n      </table>\n"
        << "    </td>\n  </tr>\n"
        << "</Form>\n";

  }

/*
  resp.out() << "  <tr>\n    <td>&nbsp;</td>  </tr>\n"
	     << "  <tr>\n    <td class=SMALL1>\n<a href=\"" << thisConf.HomePage
	     << "\">Back</a> to the HomePage\n    </td>\n  </tr>\n";

*/

} // PrintActionButtons



static void
PrintChainsMap (BiostrucPtr bsp, PDNMS ModelStruc, unsigned short chain1, 
		unsigned short cnt, CCgiContext& ctx, bool VastLink=0)
{

        BiostrucFeaturePtr 	domain_bfp = 0;
        bool         		hasDomain = FALSE;
        unsigned            	*DomIdx;
        unsigned short          imgsize=0,  i, chain2;
        PMSD            	pmsd;
        PDNMM           	pdnmm, pdnmm1;
        IntervalHead    	**DomHead;

	CCgiResponse& resp = ctx.GetResponse();
	static DataInfo Dinfo;
	static VSMmdbConf thisConf;

        pmsd = (PMSD) ModelStruc->data.ptrvalue;
        pdnmm = pmsd->pdnmmHead;
        if (pdnmm == 0) return;

        DomHead = NewDataType<IntervalHead *>(cnt+1);
        for (i=0; i< cnt+1; i++) DomHead[i] = NULL;

        CalDomIdx(pdnmm, &DomIdx);
        GetDomFeaPtr(bsp, &hasDomain, &domain_bfp);
        CheckDomColIdx(domain_bfp, pdnmm, DomHead, DomIdx, hasDomain, FALSE);

//        GroupingChains(pmsd);

        pdnmm = pmsd->pdnmmHead;
        pdnmm1 = PdnmmforChainX(pdnmm, chain1);
        pdnmm = pdnmm1;

        chain2 = MIN(cnt, chain1+ChainPerImg-1);

        resp.out() << "<map name=chain_map>\n";
       
     	static ImgInfo Iinfo; 
	Iinfo.ismap = true;
        Iinfo.chbeg = chain1;
        Iinfo.chend = chain2;
	Iinfo.subset = Dinfo.subset;
        unsigned short y =ModelMapOrImg(0,pmsd, DomHead, DomIdx, ctx, VastLink);

//        imgsize = y + FontBH + 80;
	imgsize = y + 50;

        resp.out() << "</map>\n"
		<< "<img src=\"" << thisConf.BaseUrl << thisConf.CgiName
		<< "?cmdVSMmdb=StrImg&grpid=" << Dinfo.GrpID
		<< "&chbeg=" << chain1
		<< "&chend=" << chain2
		<< "&imgsize=" << imgsize << "\""
		<< " usemap=#chain_map border=0 ismap>\n";

        delete [] DomIdx;

}       /* end PrintChainsMap */



static void SumTab(CCgiContext& ctx, vector< vector<string> > vvstr, bool isStr)
{

	CCgiResponse& resp = ctx.GetResponse();
        static DataInfo         Dinfo;
        static ImgConf          thisConf;
        static ImgInfo          Iinfo;

/*
	vector < vector<string> > vvstr;
        GetChnInfo(Iinfo.id, vvstr);
	if ( !vvstr.size() ) {
		Iinfo.prot_ch_cnt = 0;
		return;
	}
*/
        DomIntvlData dip[100];

	resp.out() << "      <table border=0>\n       <tr bgcolor=silver>\n"
                << "          <td align=middle width=%%25>"
		<< "<b>Chain ID</b></td>\n"
                << "          <td align=middle width=%%25>"
		<< "<b>Domains</b></td>\n"
                << "          <td align=middle width=%%25>"
		<< "<b>Residue Range</b></td>\n";

	if (isStr) {
	}
	else 
	    resp.out() << "          <td align=middle>"
			<< "<b>No. of Neighbors</b></td>\n";

  	resp.out() << "       </tr>\n\n";

	unsigned sdi, line_cnt=1, total_nbr=0;

	for (unsigned i=0; i< vvstr.size(); i++)
        {
             unsigned chnNo = ToUnsigned(vvstr[i][0]);
             unsigned domintvl_no =
                        GetVSDomIntvlData(Dinfo.JobID, chnNo, dip, 100);

	     for (unsigned j=0; j< domintvl_no; j++)
             {

//		 if (!j && isStr && domintvl_no > 1) continue;

                 if (dip[j].domNo >= 10000) continue;
                 sdi = DomNo2Sdi(Iinfo.id, chnNo, dip[j].domNo);
		 total_nbr = GetVSNumOfAllNbrs(sdi);

                 resp.out() << "        <tr";

                 if ( 0 == line_cnt%2 )
                        resp.out() << " bgcolor=#e5e5e5";

                 resp.out() << ">\n          <td align=middle>["
                    		<< vvstr[i][1] << "]</td>\n";
		 if (!isStr && total_nbr) {
		    resp.out() << "          <td align=middle>\n"
                    	<< "<a href=" <<thisConf.VastUrl << thisConf.VastCgi
                    	<< "?reqid=" << Dinfo.GrpID
                    	<< "&subsetstr=" << Dinfo.subset.substr(0,3)
                    	<< "&sdid=" << sdi <<">";
		 }
		 else resp.out() << "          <td align=middle>\n";
                 if ( !dip[j].domNo )		//  dip[0].domNo = j = 0;
                    resp.out() << "entire chain</a></td>\n";
                 else
                    resp.out() << "domain " <<dip[j].domNo << "</a></td>\n";

                 resp.out()
                        << "          <td align=middle>" << dip[j].frm
                        << " - " << dip[j].tu;

                 if (j) {
                    for (unsigned k=j+1; k< domintvl_no; k++) {

                        if (dip[j].domNo == dip[k].domNo) {

                                resp.out() << ", " << dip[k].frm << " - " << dip
[k].tu;
                                dip[k].domNo += 10000;
                        }
                    }
		}

                resp.out() << "</td>\n";

		if (!isStr) {

		     resp.out()
                        << "          <td align=middle>"<< total_nbr
                        << "</td>\n"
                        << "        </tr>\n\n";
		}
                line_cnt ++;
             }
        }

        resp.out() << "      </table>\n\n";


}	 // SumTab()



void SendSummaryPageText(CCgiContext& ctx, bool VastLink, bool JobDone)
{

     CCgiRequest&	req = ctx.GetRequest();
     CCgiResponse& 	response = ctx.GetResponse();
     static DataInfo 	Dinfo;
     static VSMmdbConf 	thisConf;
     static ImgInfo	Iinfo;

     BiostrucPtr bsp;
     PDNMS  ModelStruc;
     PMSD pmsdThis = 0;
     BiostrucSourcePtr pbssThis = 0;
     ValNodePtr pvnThis = 0;
     BiostrucHistoryPtr pbshThis = 0;
//     char* pcAuthors = 0;

 // bsp = FetchBS((char*)bFile.c_str(), 1, 2, 1000, POWER_VIEW);

     bsp = VSOpenBSP(Iinfo.id, 2, 1000);
       
     if (bsp == 0) {

	string strtmp = "Missing Biostruc data for reqid = "
                        + ReqId2Str(Dinfo.ReqID)
                        + ". Please inform chenj@ncbi.nlm.nih.gov for this problem. Thanks!";

     	PrtMes::PrintErrorHeader(&response, "Vast Search", 1);
	PrtMes::PrintMsgWithTail(&response, strtmp);

    	exit(1);
     }

     ModelStruc = MakeAModelstruc(bsp);

     if (ModelStruc == 0) {
	PrtMes::PrintErrorHeader(&response, "Vast Search", 1);
	PrtMes::PrintMsgWithTail(&response, "ModelStruc==NULL.\nPlease inform chenj@ncbi.nlm.nih.gov.\n");
	return;
     }
	
     pmsdThis = (PMSD) ModelStruc->data.ptrvalue;
     pvnThis
	=ValNodeFindNext(pmsdThis->pbsBS->descr,NULL,BiostrucDescr_history);
     if (pvnThis) {
    	pbshThis = (BiostrucHistoryPtr) pvnThis->data.ptrvalue;
    	pbssThis = pbshThis->data_source;

// how to delete pbshThis and pbssThis 

     }


  string urlstr = "" ;
  string dispOp = "";
  unsigned chaincnt = ChainCount(pmsdThis);
  unsigned chbeg = 1;  // (CurPage-1)*ChainPerImg +1;
  response.out() 
      << "<TABLE width=800 BORDER=0 CELLPADDING=3 CELLSPACING=0>\n\n"
	<< "  <tr>\n    <td VALIGN=TOP NOWRAP>&nbsp;</td>\n  </tr>\n\n";


  vector< vector<string> > vvstr;
  GetChnInfo(Iinfo.id, vvstr);
  Iinfo.prot_ch_cnt = 1;
  if ( !vvstr.size() ) Iinfo.prot_ch_cnt = 0;

  if (!VastLink) {
      response.out()
	<< "  <tr>\n    <td valign=TOP class=H1>"
	<< "Data Parsing Done</td>\n"
	<< "  </tr>\n\n"
	<< "  <tr> <td>&nbsp;</td></tr>\n\n"
	<< "  <tr>\n    <td valign=TOP class=SMALL1>\n"
	<< "Your structure data has been uploaded. Request ID: "
	<< "<b>" << Dinfo.GrpID << "</b>\n" 
	<< "    </td>\n  </tr>\n\n";

      if (Iinfo.prot_ch_cnt) 
	  response.out()<< "  <tr>\n    <td class=SMALL1>\n"
      	  	<< "<b>Bookmark</b> this page to return later.\n";
      response.out() << "  <tr><td width=800><br><HR></td></tr>\n\n";


      if (Iinfo.prot_ch_cnt) {
      	  response.out() << "  <tr>\n    <td class=H4>\n"
		<< "Query Chains/Domains Summary<br><br></td></tr>\n\n"
		<< "  <tr>\n    <td>\n";
      	  SumTab(ctx, vvstr, 1);
      	  response.out() << "    </td>\n  </tr>\n\n";
      }
      else {
	  response.out() << "  <tr>\n    <td class=H4>\n"
		<< "VAST Chains Summary: "
		<< "No protein chains submitted for VAST Search.\n";

	  urlstr= "'" + thisConf.IbutHelpF + "_NulCh.html'";
	  dispOp = "'resizable=yes, scrollbars=yes, width=400, height=200'";
 	  CreateIbutn(urlstr, dispOp, 15, ctx);

	  response.out() << "    </td>\n  </tr>\n\n";
      }
	
      response.out() << "  <tr>\n    <td>\n";
      PrintChainsMap(bsp, ModelStruc, chbeg, chaincnt, ctx, VastLink);
      response.out() << "    </td>\n  </tr>\n\n";

      response.out() << "  <tr>\n    <td class=SMALL1>The \n";
      if (Iinfo.prot_ch_cnt)
      	  response.out() << "table and ";

      response.out() << "graphics above indicate the individual chains and 3D domains identified in your structure. You may view them in "
                << "<a href=/Structure/CN3D/cn3d.shtml>Cn3D</a>\n"
                <<" to verify whether the data parsing is correct by clicking \"View 3D Structure\" below.\n";

        urlstr = "'" + thisConf.IbutHelpF + "_Graph.html'";
        dispOp = "'resizable=yes, scrollbars=yes, width=420, height=200'";
        CreateIbutn(urlstr, dispOp, 15, ctx);
        response.out() << "    </td>\n  </tr>\n";

	PrintActionButtons(ctx);

	response.out() << "  <tr>\n    <td>&nbsp;</td>  </tr>\n"
             << "  <tr>\n    <td class=SMALL1>\n<a href=\"" << thisConf.HomePage
             << "\">Back</a> to the HomePage\n    </td>\n  </tr>\n";


  	ctx.GetResponse().out() << "</TABLE>\n" << "<br>\n";

  }
  else {

    string errmsg;

    if (HaveErrmsgFromVSCal(Dinfo.GrpID, errmsg)) {

	response.out() << "  <tr>\n"
	    << "    <td valign=TOP class=H1>VAST Search Error Report</td>\n"
	    << "  </tr>\n\n"
            << "  <tr> <td>&nbsp;</td>\n  </tr>\n\n"
	    << "  <tr>\n    <td valign=TOP class=SMALL1>\n"
	    << "Your VAST search has stopped due to the following error.\n"
	    << "    </td></tr>\n"
	    << "  <tr><td class=TEXT><b>"
	    << errmsg << "</b>"
	    << "    </td></tr>\n"
	    << "  <tr>\n    <td>&nbsp;</td>  </tr>\n"
	    << " <tr>\n    <td width=800><HR><br></td>  </tr>\n"
             << "  <tr>\n    <td class=SMALL1>\n<a href=\"" << thisConf.HomePage
             << "\">Back</a> to the HomePage\n    </td>\n  </tr>\n";



    }
    else if (JobDone) {
     
      response.out()
	<< "  <tr>\n"
	<< "    <td valign=TOP class=H1>VAST Search Done</td>\n  </tr>\n\n"
	<< "  <tr> <td>&nbsp;</td>\n  </tr>\n\n"
	<< "  <tr>\n    <td valign=TOP class=SMALL1>\n"
	<< "Your VAST search has finished successfully. Follow links on the Query Structure Summary to see results.\n"
	<< "    </td></tr>\n";

/*
      urlstr = "'" + thisConf.IbutHelpF + "_VSlink.html'";
      dispOp = "'resizable=yes, scrollbars=yes, width=420, height=200'";
      CreateIbutn(urlstr, dispOp, 15, ctx);
*/

      response.out() << "    <tr><td width=800><br><HR></td></tr>\n\n"
		<< "  <tr>\n    <td class=H4>\n"
                << "Query Structure Summary<br><br>"
                << "    </td>\n  </tr>\n\n  <tr>\n    <td>\n";

      SumTab(ctx, vvstr, 0);
      response.out() << "    </td>\n  </tr>\n\n";
/*
                << "    <tr><td width=800><br><HR></td></tr>\n\n"
		<< "    <tr><td><br></td></tr>\n\n";
*/

      string showImg = req.GetEntry("showImg_fin").GetValue();
      if ("Yes" == showImg) {
	
/*
	  response.out() << "  <tr><td class=H4>"
			<< "Query Structure Summary<br></td>\n"
*/
	  response.out() << "    <tr><td><br></td></tr>\n\n"
			<< "  </tr>\n\n  <tr>\n    <td>\n";
      	  PrintChainsMap(bsp, ModelStruc, chbeg, chaincnt, ctx, VastLink);
      	  response.out() << "    </td>\n  </tr>\n\n"
                << "  <tr>\n    <td class=SMALL1>\n"
                << "The graphics above indicate the individual chains and 3D domains identified in your structure. You may view them in "
                << "<a href=/Structure/CN3D/cn3d.shtml>Cn3D</a>\n"
                <<" by clicking \"View 3D Structure\" below. "
		<< "You may also click each of them to look at their VAST neighbors.\n";

        urlstr = "'" + thisConf.IbutHelpF + "_Graph.html'";
        dispOp = "'resizable=yes, scrollbars=yes, width=420, height=200'";
        CreateIbutn(urlstr, dispOp, 15, ctx);

/*
        response.out() << "    </td>\n  </tr>\n"
                << "    <tr><td width=800><HR></td></tr>\n\n";
*/

	VastLink=false;
        PrintActionButtons(ctx);

	response.out() << "    <tr><td width=800><HR></td></tr>\n\n"
		<< "    <tr><td>&nbsp;</td></tr>\n\n"
		<< "    <tr><td class=SMALL1>\n"
		<< "<a href=\"" << thisConf.BaseUrl << thisConf.CgiName << "?"
                << "cmdVSMmdb=StrText&grpid=" << Dinfo.GrpID
		<< "\">Close</a> the graphical Summary.\n    </td></tr>\n"
		<< "  <tr>\n    <td>&nbsp;</td>  </tr>\n"
                << "  <tr>\n    <td class=SMALL1>\n"
                << "<a href=\"" << thisConf.HomePage
                << "\">Back</a> to the HomePage\n    </td>\n  </tr>\n"
		<< "</TABLE>\n\n";
 
	
      }
      else {
	
	    response.out() 
                << "    <tr><td width=800><br><HR></td></tr>\n\n"
                << "    <tr><td><br></td></tr>\n\n"
		<< "  <tr>\n    <td class=SMALL1>\n"
		<< "<a href=\"" << thisConf.BaseUrl << thisConf.CgiName << "?"
		<< "cmdVSMmdb=StrText&grpid=" << Dinfo.GrpID
		<< "&showImg_fin=Yes\">View</a>"
		<< " the Summary (graphical).\n    </td>\n  </tr>\n"
		<< "  <tr>\n    <td>&nbsp;</td>  </tr>\n"
             	<< "  <tr>\n    <td class=SMALL1>\n"
		<< "<a href=\"" << thisConf.HomePage
             	<< "\">Back</a> to the HomePage\n    </td>\n  </tr>\n"
		<< "</TABLE>\n\n";
      }
	
    }
    else { 

       unsigned DomCnt = TotalDom(Dinfo.JobID);
    
       sqmInit();
       string subTime = sqmGetReqTime(GrpId2ReqId(Dinfo.GrpID));
       sqmFree();

       string curTime  = CTime(CTime::eCurrent).AsString();

       string subsetstr = GetSubsetName4Grpid(Dinfo.GrpID);
       unsigned estT = 3 + (int)((DomCnt-1)*2.5);
       if (subsetstr == "All") estT *= 10;

       unsigned subH = atoi(subTime.substr(11,2).c_str());
       unsigned subM = atoi(subTime.substr(14,2).c_str());
       unsigned curH = atoi(curTime.substr(11,2).c_str());
       unsigned curM = atoi(curTime.substr(14,2).c_str());
       if (curH < subH) curH += 24;

       unsigned waiT = (curH-subH)*60 + (curM-subM); 

       unsigned refT;
       ostringstream output;
       string   refMsg;
       if ( estT <= waiT ) {
	  if ( waiT < 60 ) {

		refT = 5;
		output << "Your job needs more time to finish, its status will be checked automatically after "
			<< "<b>" << refT << "</b> minutes.";
	  }
	  else if (waiT < 120) {

		refT = 10;
                output << "Your job hasn't finised after running for <b>1</b> hour, its status will be checked automatically after "
                        << "<b>" << refT << "</b> minutes.";
	  }
	  else {

		refT = 0;
		output << "Your job has been running longer than <b>2</b> hours, it may need more time to finish, please "
			<< "<a href=" << thisConf.BaseUrl << thisConf.CgiName
			<< "?cmdVSMmdb=StrText&ViewNbr=Yes&grpid=" <<Dinfo.GrpID
			<< ">check</a> your job status later manually.";

	  }
       }
       else {

		refT = estT - waiT;
		output << "\n      <table width=800 border=0 cellpadding=0>\n"
			<< "      <tr><td class=SMALL1 width=30%>"
		 	<< "Estimated completion Time: </td>\n"
			<< "          <td class=SMALL1>"
			<< "<b>" << refT << "</b> minutes.</td>\n      </tr>\n"
			<< "      <tr><td>&nbsp;</td>\n"
			<< "          <td class=SMALL1>"
			<< "(actual time depends on the number of chains, "
			<< "domains, and the size of the query protein; "
			<< "highly populated folds such as a TIM barrel "
			<< "will take a longer time).</td>\n      </tr>\n\n"
			<< "      <tr><td class=SMALL1 colspan=2>"
			<< "The job status will be checked automatically after "
			<< "<b>" <<refT << "</b> minutes.</td>\n      </tr>\n\n" 
			<< "      </table>\n\n";
       }

       if (refT) {

/*
       	   response.out() << "<script>\n"
		<< "  <!--\n"
		<< "       setTimeout(\'document.location.replace(\""
        	<< thisConf.MonitUrl << thisConf.MonitCgi
		<< "?cmdSch=Status&ReqId=" << Dinfo.ReqID
        	<< "&GrpId=" << Dinfo.GrpID
		<< "&QMQue=" << thisConf.QMqueNbrName
      		<< "\");\', "
                << refT*60000 << ")\n"
                << "  //-->\n</script>\n";
*/

	   response.out() << "<script>\n"
                << "  <!--\n"
                << "       setTimeout(\'document.location.replace(\""
		<< thisConf.BaseUrl << thisConf.CgiName
		<< "?cmdVSMmdb=StrText&ViewNbr=Yes&grpid=" << Dinfo.GrpID
                << "\");\', "
                << refT*60000 << ")\n"
                << "  //-->\n</script>\n";
	}


       response.out() 
		<< "  <tr>\n    <td valign=TOP class=H1>"
		<< "VAST Search in Progress<br>\n"
		<< "    </td>\n  </tr>\n\n"
		<< "  <tr><td>&nbsp;</td></tr>\n\n"
		<< "  <tr><td valign=TOP class=SMALL1>\n"
		<< "Your VAST Search job was submitted at "
		<< subTime
		<< "(EDT). Request ID: <b>" <<Dinfo.GrpID <<"</b></td></tr>\n\n"
		<< "  <tr><td class=SMALL1 colspan=2>" << output.str() 
		<< "</td></tr>\n\n"
		<< "  <tr><td>&nbsp;</td></tr>\n\n"
		<< "  <tr>\n    <td class=SMALL1>\n"
		<< "<b>Bookmark</b> this page to return and manually \n"
		<< "<a href=" << thisConf.BaseUrl << thisConf.CgiName
		<< "?cmdVSMmdb=StrText&ViewNbr=Yes&grpid=" << Dinfo.GrpID
		<< ">check</a> later.\n    </td>\n  </tr>\n\n";
/*
		<< "<a href=" << thisConf.MonitUrl << thisConf.MonitCgi
		<< "?cmdSch=Status&ReqId=" << Dinfo.ReqID
		<< "&GrpId=" << Dinfo.GrpID<< "&QMQue=" << thisConf.QMqueNbrName
		<< ">check</a> later.\n Request ID: "
		<< Dinfo.GrpID 
		<< ".\n    </td>\n  </tr>\n\n"
		<< "  <tr>\n    <td valign=TOP class=SMALL1>\nOr click "
		<< "<a href=\"" << thisConf.HomePage << "#ChkRes"
		<< "\">here</a> to check previous jobs.\n    </td>\n  </tr>\n\n"
		<< "  <tr>\n    <td>&nbsp;</td>  </tr>\n";
*/

       string showImg = req.GetEntry("showImg_mid").GetValue();
       if ( "Yes" == showImg ) {
	
	   response.out() << "  <tr><td><HR></td></tr>\n\n"
		<< "  <tr><td class=H4>\n"
		<< "Query Structure Summary<br></td>"
		<< "  </tr>\n\n"
		<< "  <tr>\n    <td>\n";
	   PrintChainsMap(bsp, ModelStruc, chbeg, chaincnt, ctx, VastLink);
           response.out() << "    </td>\n  </tr>\n\n"
                << "  <tr>\n    <td class=SMALL1>\n"
                << "The graphics above indicate the individual chains and 3D domains identified in your structure. You may view them in "
                << "<a href=/Structure/CN3D/cn3d.shtml>Cn3D</a>\n"
                <<" by clicking \"View 3D Structure\" below. ";

           urlstr = "'" + thisConf.IbutHelpF + "_Graph.html'";
           dispOp = "'resizable=yes, scrollbars=yes, width=420, height=200'";
           CreateIbutn(urlstr, dispOp, 15, ctx);
	 
	   VastLink=false;
           PrintActionButtons(ctx);

	   response.out() << "    <tr><td width=800><HR></td></tr>\n\n"
                << "    <tr><td>&nbsp;</td></tr>\n\n"
                << "    <tr><td class=SMALL1>\n"
                << "<a href=\"" << thisConf.BaseUrl << thisConf.CgiName << "?"
                << "cmdVSMmdb=StrText&grpid=" << Dinfo.GrpID
                << "\">Close</a> the graphical Summary .\n    </td></tr>\n"
                << "  <tr>\n    <td>&nbsp;</td>  </tr>\n";
	}
	else {
	   
	   response.out() << "  <tr><td width=800><HR></td></tr>\n\n"
                << "  <tr><td>&nbsp;</td></tr>\n\n"
	   	<< "  <tr>\n    <td class=SMALL1>\n"
                << "<a href=\"" << thisConf.BaseUrl << thisConf.CgiName << "?"
                << "cmdVSMmdb=StrText&grpid=" << Dinfo.GrpID
                << "&showImg_mid=Yes\">View</a>"
                << " the Query Structure Summary.\n    </td>\n  </tr>\n"
		<< "  <tr><td>&nbsp;</td></tr>\n\n";

	}


	response.out() << "  <tr>\n    <td class=SMALL1>\n<a href=\"" 
		<< thisConf.HomePage
                << "\">Back</a> to the HomePage\n    </td>\n  </tr>\n"
		<< "</TABLE>\n\n";

    }

  }   // ViewNbr


  return;

} // SendSummaryPage 


