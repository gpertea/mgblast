/*
* $Id: PrintBanner_comb.cpp,v 1.1 2005/07/26 17:12:50 chenj Exp $
*
*
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
*
* Author:  Jie Chen
*
*
* $Log: PrintBanner_comb.cpp,v $
* Revision 1.1  2005/07/26 17:12:50  chenj
* Making linux VSNbr.cgi
*
* Revision 1.2  2003/01/22 17:12:50  chenj
* Get Cn3D 4.1
*
* Revision 1.1.1.1  2002/12/06 20:17:21  chenj
* Imported Scouces
*
*
*
*
* This file is used to print the html head.
*
* ==========================================================================
*/


#include "hUtilib.hpp"
#include "vastuti.hpp"
#include <sstream>

#define DISP_IN_GRAPH 	'n'
#define DISP_IN_TABLE	'y'

extern	char 		URLcgi[PATH_MAX];
extern	char		CGIname[PATH_MAX];
extern 	char		MMDBCGIname[PATH_MAX], DATApath[PATH_MAX];
extern 	Uint1		numSubsets;
extern 	bool		SubsetInfoLoaded;
extern	unsigned	aSdi, aMmdbId;
extern	string		JobID, Passwd;
extern	long long	ReqId;
	
void PrintHitsSortBanner(FILE *table, SortBy sortby, SubSetID subsetnum,
	unsigned pagenum, unsigned numpages, char cTable)
{
	unsigned i;
        static char subsetname[6][100]={"Low redundancy",
                                        "Medium redundancy",
                                        "High redundancy",
                                        "Non-identical seq.",
                                        "All sequences",
                                        "Selected neighbor(s)"};

        WWWInfoPtr    www_info;
 	if (WWWGetArgs(&www_info) != WWWErrOk);
	int idx = -1;
	char * subsetstr=NULL;
	if ((idx = WWWFindName(www_info, (char *)"subsetstr")) >=0) {
	    subsetstr = WWWGetValueByIndex(www_info, idx);
	}
        //else / error;

	/* subset 1 should be "All of PDB"; put it second last */
        if (! SubsetInfoLoaded) loadSubsetInfo();

	fprintf(table, "<TR>\n<TD>\n");

	vector<string> subsetdisp;
	if (subsetstr) {
	    
             fprintf(table, "<input type=hidden name=subsetstr value=\"%s\">\n", subsetstr);
	     if ((string)subsetstr == "Non") {
	
	     	subsetdisp.reserve(3);
	     	subsetdisp.push_back(GetSubsetName_DB(2));
	     	subsetdisp.push_back(GetSubsetName_DB(3));
	     	subsetdisp.push_back(GetSubsetName_DB(6));
	     }
	     else {
		     subsetdisp.reserve(6);
		     for (i=2; i< numSubsets; i++) 
	                subsetdisp.push_back( GetSubsetName_DB(i));
		     subsetdisp.push_back(GetSubsetName_DB(1));
		     subsetdisp.push_back(GetSubsetName_DB(6)); 
	     }
	}
					

/* Use the Form created in PrintAlignViewBanner */
	fprintf(table, "<input type=hidden name=presubset value=%d>\n",
			subsetnum);

   	string altstr = "Click to display and sort the neighbors";
    	fprintf(table, "  <table cellpadding=0>   ");
	fprintf(table, "<!-- Subtable for display/sort hits -->\n");

/* Insert a blank row */
        fprintf(table, "    <tr><td></td></tr>\n\n");

    	fprintf(table, "    <tr>\n      <td>\n          ");
	fprintf(table, "<input type=image name=dispsub value=List ");
	fprintf(table, "src=\"%slist.gif\" ", DATApath);
	fprintf(table, "alt=\"%s\" title=\"%s\">\n", altstr.c_str(), altstr.c_str());
	fprintf(table, "      </td>\n\n");

	fprintf(table, "      <td class=SMALL1 valign=middle>\n");
	fprintf(table, "          &nbsp;&nbsp;\n"); 

	fprintf(table, "          <select name=subset>\n");
	
/* subset 1 should be "All of PDB"; put it second last */
	for (i=0; i< 2; i++) {
	     fprintf(table, "            <option value=\"%s\"", 
							subsetdisp[i].c_str());
	     if (i+2 == subsetnum) fprintf(table, " selected");
	     fprintf(table, ">%s\n", subsetname[i]); 
	}

	if (subsetdisp.size() - i == 1) {
	
	    fprintf(table, "            <option value=\"%s\"", 
						subsetdisp[2].c_str());
	    if (6 == subsetnum) fprintf(table, " selected");
	    fprintf(table, ">%s\n", subsetname[5]);
	}
	else {
	    for (i=2; i< subsetdisp.size()-2; i++) {

		fprintf(table,"            <option value=\"%s\"",
							subsetdisp[i].c_str());
                if (i+2 == subsetnum) fprintf(table, " selected");
             	fprintf(table, ">%s\n", subsetname[i]);
	    }

	    fprintf(table,"            <option value=\"%s\"", 
							subsetdisp[4].c_str());
	    if ( subsetnum == 1 ) fprintf(table, " selected");
	    fprintf(table, ">%s\n", subsetname[4]); 

	    fprintf(table, "            <option value=\"%s\"", 
							subsetdisp[5].c_str());
	    if (subsetnum ==6) fprintf(table, " selected");
	    fprintf(table, ">%s\n", subsetname[5]); 
	}

	fprintf(table, "          </select>\n\n");

	fprintf(table, "          <b>&nbsp;subset, sorted by&nbsp;\n\n");

	fprintf(table, "          <select name=sort>\n");
	fprintf(table, "            <option value=%d", NRes);
	if (!sortby || sortby== NRes) fprintf(table, " selected");
	fprintf(table, ">Aligned Length\n");
	fprintf(table, "            <option value=%d", VScore);
	if (sortby == VScore) fprintf(table, " selected");
	fprintf(table, ">Vast Score\n");
	fprintf(table, "            <option value=%d", PValue);
	if (sortby== PValue) fprintf(table, " selected");
	fprintf(table, ">Vast P_value\n");
	fprintf(table, "            <option value=%d", Rmsd);
	if (sortby == Rmsd)
		fprintf(table, " selected");
	fprintf(table, ">Rmsd\n");
	fprintf(table, "            <option value=%d", PcntId);
	if (sortby == PcntId) fprintf(table, " selected");
	fprintf(table, ">Seq. Identity\n");
	fprintf(table, "          </select>\n\n");

	if (pagenum && numpages > 1) {
	    fprintf(table, "          <b>&nbsp;page&nbsp;\n\n");
	    fprintf(table, "          <select name=doclistpage>\n");
	    for (i =1; i<= numpages; i++) { 
	      if (i == pagenum)
	       fprintf(table, "            <option value=%d selected>%d\n",i,i);
	      else fprintf(table, "            <option value=%d>%d\n", i,i);
	      }
	    fprintf(table, "          </select>\n\n");
	}

	fprintf(table, "          &nbsp;in&nbsp;\n");
	fprintf(table, "          <select name=table>\n");
	fprintf(table, "            <option value=n");
	if (cTable == DISP_IN_GRAPH) 
		fprintf(table, " selected");
	fprintf(table, ">Graphics\n");
	fprintf(table, "            <option value=y");
	if (cTable == DISP_IN_TABLE)
		fprintf(table, " selected");
	fprintf(table, ">Table\n");
	fprintf(table, "            <option value=s>ASN1\n"); 
	fprintf(table, "          </select>\n\n");

    	string urlstr = (string)"'" + DATApath + "help_List.html'";
    	string dispOpstr= "'resizable=yes, scrollbars=yes, width=420, height=250'";
    	fprintf(table, "         <a href=\"#\" onclick=\"PopWin(%s, %s);\">\n", 
                                        urlstr.c_str(), dispOpstr.c_str());
    	fprintf(table, "         <img src=\"%sinfosmall.gif\" border=0></a>\n",
                                                                DATApath);

	fprintf(table, "      </td>\n    </tr>\n  </table>\n</TD>\n</TR>\n\n");

} /* end of PrintHitsSortBanner */



void 
PrintAlignViewBanner(FILE *table, unsigned iFSID, VastPageDataPtr vpp, 
	unsigned numhitsdisplayed)
{
    Int4  i;
    Char  str[50];
    string AllId;

    fprintf(table, "<TR>\n<TD>\n"); 

/* Create one Form with 3 submits */

    fprintf(table, "<Form method=get action=\"%s%s\">\n", URLcgi, CGIname);
    fprintf(table, "<input type=hidden name=sdid value=%d>\n", aSdi);

    if (JobID != "") {
        if (ReqId) 
	    fprintf(table, "<input type=hidden name=reqid value=%llu>\n",ReqId);
 	else {
	    fprintf(table, "<input type=hidden name=vsid value=\"%s\">\n", 
								JobID.c_str());
	    fprintf(table, "<input type=hidden name=pass value=\"%s\">\n", 
								Passwd.c_str());
	}
    }

    ostringstream output;
    for (i=0; i< numhitsdisplayed; i++) {

	output.str("");
	output << getfid(vpp[i].bSdi, vpp[i].bAlignId) << ",";
	AllId += output.str();
    }
    fprintf(table, "<input type=hidden name=allbfid value=\"%s\">\n", 
							AllId.c_str());

    fprintf(table, "<input type=hidden name=defhit value=%d>\n",
				getfid(vpp[0].bSdi, vpp[0].bAlignId));
     //        fprintf(table, "<input type=hidden name=subsetstr value=\"Non\">\n");

    fprintf(table, "  <table cellpadding=0>");
    fprintf(table, "   <!-- Subtable for View/Save Structure-->\n");

    string altstr = "Click to view 3D superpositions";
    fprintf(table, "    <tr>\n      <td>\n         ");
    fprintf(table, "<input type=image name=viewstr src=\"%sv3d.gif\"",DATApath);
    fprintf(table, " alt=\"%s\" title=\"%s\">\n",altstr.c_str(),altstr.c_str());
    fprintf(table, "      </td>\n\n");

    fprintf(table, "      <td class=SMALL1 valign=middle>\n");
    fprintf(table, "         <b>&nbsp;&nbsp;of&nbsp;\n\n");

    fprintf(table, "         <select name=atm_complexity>\n");
    fprintf(table, "            <option value=0>All Atoms\n");
    fprintf(table, "            <option value=1>Backbone\n");
    fprintf(table, "         </select>\n\n");

    fprintf(table, "         &nbsp;with&nbsp;\n\n");

    fprintf(table, "         <select name=calltype>\n");
    fprintf(table, "            <option value=c>Cn3D\n");
    fprintf(table, "            <option value=a>Cn3D/Cache\n");
    fprintf(table, "         </select>\n\n");

    fprintf(table, "         <select name=action>\n");
    fprintf(table, "            <option value=0>Display\n");
    fprintf(table, "            <option value=1>See File\n");
    fprintf(table, "            <option value=2>Save File\n");
    fprintf(table, "         </select>\n\n");

    string urlstr = (string)"'" + DATApath + "help_V3D.html'";
    string dispOpstr= "'resizable=yes, scrollbars=yes, width=420, height=250'";
    fprintf(table, "	     <a href=\"#\" onclick=\"PopWin(%s, %s);\">\n", 
					urlstr.c_str(), dispOpstr.c_str());
    fprintf(table, "	     <img src=\"%sinfosmall.gif\" border=0></a>\n",
								DATApath);
    fprintf(table, "         <img src=\"/Structure/new.gif\" alt=\"New\" title=\"New\">\n");
    fprintf(table, "         <a href=\"/Structure/CN3D/cn3d.shtml\"><I>Get Cn3D 4.1!</I></B></a>\n");

    fprintf(table, "      </td>\n    </tr>\n  </table>\n");
    fprintf(table, "</TD>\n</TR>\n\n");



    fprintf(table, "<TR>\n<TD>\n");
    fprintf(table, "  <table cellpadding=0>  "); 
    fprintf(table, "<!-- Subtable for View Alignments -->\n");

/* Insert a blank row */
    fprintf(table, "    <tr><td></td></tr>\n\n");


    altstr = "Click to view sequence alignment";
    fprintf(table, "    <tr>\n      <td>\n         ");
    fprintf(table,"<INPUT TYPE=image NAME=viewali src=\"%svseq.gif\"",DATApath);
    fprintf(table, " alt=\"%s\" title=\"%s\">\n",altstr.c_str(),altstr.c_str());
    fprintf(table, "      </td>\n\n");

    fprintf(table, "      <td class=SMALL1 valign=middle>\n");
    fprintf(table, "         <b>&nbsp;&nbsp;using&nbsp;\n\n");
    
    fprintf(table, "         <select name=alitype>\n");
    fprintf(table, "            <option value=h>Hypertext\n");
    fprintf(table, "            <option value=t>Plain text\n");
    fprintf(table, "            <option value=f>mFASTA\n");
    fprintf(table, "         </select>\n\n");
 
    fprintf(table, "         &nbsp;for&nbsp;\n\n");

    fprintf(table, "         <select name=nbr_complexity>\n");
    fprintf(table, "            <option value=1>Selected\n");
    fprintf(table, "            <option value=0>All on page\n");
    fprintf(table, "         </select>\n\n");
    
    fprintf(table, "         &nbsp;VAST neighbors</b>\n");

    urlstr = (string)"'" + DATApath + "help_VSEQ.html'";
    fprintf(table, "         <a href=\"#\" onclick=\"PopWin(%s, %s);\">\n", 
                                        urlstr.c_str(), dispOpstr.c_str());
    fprintf(table, "         <img src=\"%sinfosmall.gif\" border=0></a>\n",
                                                                DATApath);
    fprintf(table, "      </td>\n    </tr>\n  </table>\n");
    fprintf(table, "</TD>\n</TR>\n\n");

} /* end PrintAlignViewBanner */



void
PrintQueryInfo(FILE *table, char *pcPDB, char cChain,int iDomain)
{

	  Char 	title[MAX_TBUFF];
	
	  fprintf(table, "<TR>\n<TD>\n");
	  fprintf(table, "  <table cellpadding=0>  <!-- Subtable for description-->\n");

/* insert two blank rows for spacing */
	  fprintf(table, "    <tr>\n");
	  fprintf(table, "      <td>&nbsp; </td>\n");
	  fprintf(table, "      <td>&nbsp; </td>\n");
	  fprintf(table, "    </tr>\n\n");

	  fprintf(table, "    <tr>\n");
	  fprintf(table, "      <td width=15%% VALIGN=TOP NOWRAP ALIGN=RIGHT class=H4>\n");
	  if (!ReqId) 
			fprintf(table, "Query: </TD>\n\n<TD class=H4>\n");
	  else 
		fprintf(table, "             Request Id: &nbsp;%llu", ReqId);

          if (cChain != ' ' || iDomain >0) fprintf(table, ",&nbsp;\n");
	  else fprintf(table, "\n");
	  fprintf(table, "      </td>\n");
	  if (ReqId) fprintf(table, "      <td class=H4 width=85%%>");
	  else fprintf(table, "      <td class=TEXT>");
	  if (aMmdbId && JobID == "") {
	     fprintf(table, "MMDB&nbsp;");
	     fprintf(table, "<a href=\"%s%s?uid=%s",URLcgi, MMDBCGIname, pcPDB);
	     fprintf(table, "&Dopt=s\">%d</a>", aMmdbId);
	     fprintf(table, ",&nbsp;");
	  }

	  if (!ReqId) 
		if (pcPDB[0] != NULLB) fprintf(table, "%s", pcPDB); 
  
	  if (cChain != ' ')
		fprintf(table, "&nbsp;chain&nbsp;%c", cChain);

	  if (iDomain > 0)
		fprintf(table, "&nbsp;domain&nbsp;%d", (int) iDomain);
	
	  fprintf(table, "      </td>\n    </tr>\n\n");
	  if (JobID != "" && constructPdbDescr(aMmdbId, title, MAX_TBUFF)) {
	      fprintf(table, "    <tr>\n");
              fprintf(table, "      <td VALIGN=TOP NOWRAP ALIGN=RIGHT class=H4>\n");
              fprintf(table, "<b>Description:</b></TD>\n");
              fprintf(table, "      <td class=TEXT>%s", title);
              fprintf(table, "      </td>\n    </tr>\n\n");
	  }


	  fprintf(table, "  </table>\n");
	  fprintf(table, "</TD>\n</TR>\n\n");

} /* end PrintQueryInfo */



void PrintSearchNbr(FILE *table)
{

	fprintf(table, "<TR>\n<TD>\n");
	fprintf(table, "  <table cellpadding=0>       ");
	fprintf(table, "<!-- Subtable for serach of str. neighbor -->\n");

/* Insert a blank row */
        fprintf(table, "    <tr><td></td></tr>\n\n");

/* use the form created in PrintAlignViewBanner */

    	string altstr = "Click to look for neighbors";
	fprintf(table, "    <tr>\n      <td>\n           ");
	fprintf(table, "<input type=image name=schsub value=Find ");
	fprintf(table, "src=\"%sFind.gif\" ", DATApath);
	fprintf(table, "alt=\"%s\" title=\"%s\">\n", altstr.c_str(), altstr.c_str());
	fprintf(table, "      </td>\n\n");

	fprintf(table, "      <td class=SMALL1 valign=middle>\n");
	fprintf(table, "         <b>&nbsp;MMDB or PDB ids:&nbsp;\n");
	fprintf(table, "         <textarea name=selnbr cols=10 rows=1></textarea>\n");
	fprintf(table, "         &nbsp;&nbsp;or 3D-Domain ids:&nbsp;\n");
	fprintf(table, "         <textarea name=selsdid cols=10 rows=1></textarea>");

    	string urlstr = (string)"'" + DATApath + "help_Find.html'";
    	string dispOpstr= "'resizable=yes, scrollbars=yes, width=420, height=150'";
    	fprintf(table, "         <a href=\"#\" onclick=\"PopWin(%s, %s);\">\n", 
                                        urlstr.c_str(), dispOpstr.c_str());
    	fprintf(table, "         <img src=\"%sinfosmall.gif\" border=0></a>\n",
                                                                DATApath);


	fprintf(table, "      </td>\n    </tr>\n\n");
	fprintf(table, "  </table>\n");
	fprintf(table, "</TD>\n</TR>\n\n");


}	/* PrintSearchNbr */

