/*
 * $Id: vastsrv_comb.cpp,v 1.1 2005/07/26 17:12:50 chenj Exp $
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
 * Author: Christopher Hogue, Tom Madej, Siqian He, Jie Chen
 *
 *
 * $Log: vastsrv_comb.cpp,v $
 * Revision 1.1  2005/07/26 17:12:50  chenj
 * Making linux VSNbr.cgi
 *
 * Revision 1.8  2003/01/31 16:11:44  chenj
 * set ODBCINI correctly
 *
 * Revision 1.7  2003/01/15 17:15:39  chenj
 * delete all SelfNbr stuff
 *
 * Revision 1.6  2003/01/15 16:24:52  chenj
 * Change some posting messages
 *
 * Revision 1.5  2003/01/14 21:15:19  chenj
 * change VastToCn3d() to VastToCn3DAndAli()
 *
 *
 *
 * This file with others together produce graphical display of structure 
 * alignment.
 *
 *
 * ===========================================================================
 *
 */

/* Main program for VAST structure neighbor server. */

/* Use vastlocl.c instead of vastdtbs.c because it's faster to read in the
 * existing BiostrucAnnotSet file than to load it from database.
 */

#include <corelib/ncbiexec.hpp>
#include <corelib/ncbifile.hpp>
#include <corelib/ncbitime.hpp>
#include "SHGlobal.hpp"
#include "hUtilib.hpp"
#include "VastSrchUti.hpp"
#include "vastuti.hpp"
#include <cddapi.h>
#include "dartutil.h"

#define CPUTIME_MAX		120
#define DEFAULT_SUBSET_NUM    	NRBlast10e_40 	/* the NR set BLAST 10e-40 */
#define DEFAULT_SORT_BY		NRes 	/* Number of Aligned Residues */
#define NUM_HITS_PER_PAGE	60
#define DEFAULT_PAGE            1
#define DEF_ALL_NBR		1
#define SELECTED_NBR		0
#define VIEW_STR_SUBMIT		0
#define VIEW_ALI_SUBMIT		1
#define DEF_LIST_SUBMIT		2
#define FIND_SUBMIT		3
#define MAX_KEPT_NBR		30
#define DOMID_SIZE              7       /* used to be 6, as length of domId */
#define MaxVPP			150

#undef HEADFILE
#define HEADFILE   "sshead.txt"

#undef TAILFILE
#define TAILFILE   "sstail.txt"

char 		VSPATH[PATH_MAX];
static char 	URLBase[PATH_MAX];
char 		URLcgi[PATH_MAX];
char 		DATApath[PATH_MAX];
static char 	VASTpath[PATH_MAX];
char 		ENTREZurl[PATH_MAX];
static char 	DOCSUMurl[PATH_MAX];
string 		MAILto;
char 		CGIname[PATH_MAX];
char 		MMDBCGIname[PATH_MAX];
static char 	HELPname[PATH_MAX];
char		Database[PATH_MAX];

static char 	DART[PATH_MAX];
static char 	LOGIN[PATH_MAX];
static char 	PASSWD[PATH_MAX];
static char 	INITpath[PATH_MAX];
static char 	LIBpath[PATH_MAX];
char		CDDurl[PATH_MAX];
unsigned	aSdi=0, aMmdbId=0, aChnNo=0;
string		JobID, Passwd="";
long long 	ReqId=0;

Uint1 			numSubsets = 0;
SubsetNameData 	subsetNames[6];
bool 	SubsetInfoLoaded = false;
extern BiostrucAnnotSetPtr LocalGetBiostrucAnnotSet(unsigned mmdbid);
Boolean		psok;

using namespace ncbi;
using namespace SHProjNS;

static Uint1 getSubsetNum(char *subname)
{
  Uint1 i = 0;
  if(! SubsetInfoLoaded) loadSubsetInfo();
  for(i = 0; i < numSubsets; i++) 
    if(strcmp(subsetNames[i].subName,subname)==0)
      return subsetNames[i].subId;

  PrtMesC(MAILto, "VASTSRV", "GetSubsetNumber_DB() failed.", "", false);
  exit(1);
  

} /* getSubsetNbr() */



static void 
VastTableBegin (FILE *table, unsigned iFSID, SortBy sortby, SubSetID subsetnum, 
	unsigned iKept,  char * NonNbr, VastPageDataPtr vpp, 
	unsigned numhitsdisplayed, unsigned numhits, unsigned numpages,
	unsigned pagenum, char cTable, unsigned NbrFlag)
{
	unsigned domNo, i, numNbrs;
	char pdbId[5], chainLett;
	string ThisSubset;
	static char subsetname[6][100]={"All of MMDB",
					"Low redundancy",
                                        "Medium redundancy",
                                        "High redundancy",
                                        "Non-identical seq.",
                                        "Selected neighbor(s)"};

	StrCut(pdbId, vpp[0].aDomName, 1, 4);
	chainLett = vpp[0].aDomName[4];
	domNo = (unsigned) atol(vpp[0].aDomName+5);

        fprintf(table, "Content-type: text/html\n\n");
        WWWPrintFileData((char *)HEADFILE,  stdout);

        fprintf(table, "<TABLE width=800 BORDER=0 CELLPADDING=3 CELLSPACING=0 bgcolor=#FFFFCC>\n\n");
        PrintQueryInfo(table, pdbId, chainLett, domNo);

	fprintf(table, "<TR>\n<TD><br></TD>\n</TR>\n\n");
        PrintAlignViewBanner(table, iFSID, vpp, numhitsdisplayed);

        PrintHitsSortBanner(table, sortby, subsetnum, pagenum, numpages,cTable);

	PrintSearchNbr(table);

        fprintf(table, "</TABLE>\n\n");
        fflush(table);

    fprintf(table, "<br>\n"); 
    fprintf(table, "<font class=SMALL1><strong>");
    if (NbrFlag == DEF_ALL_NBR) {

	for (i=0; i< numSubsets; i++) {
	    if (i == subsetnum-1) {
		ThisSubset =  subsetname[i];
		break;
	    }
	}

	if (i) {
	    if (JobID == "") numNbrs = getNumOfLiveVastNeighbors(aSdi);
	    else numNbrs = GetVSNumOfAllNbrs(aSdi); // New VS
    	    fprintf(table, "%d neighbors found.\n&nbsp;", numNbrs);
	}

	if (numhitsdisplayed > iKept) {
            int numtmp;

            numtmp  = numhits - numhitsdisplayed + iKept;
            if (numtmp)
                fprintf(table, "%d out of %d ", numhitsdisplayed-iKept,numhits);
            else fprintf(table, "%d ", numhits);
 	    if (i) {
	    	if (numhitsdisplayed > 1) fprintf(table, "representatives");
            	else fprintf(table, "representative");  
	    }
	    else {
		if (numhitsdisplayed > 1) fprintf(table, "neighbors");
                else fprintf(table, "neighbor");
            }
	    if (i && i < numSubsets-1) {
		fprintf(table, " from the <a href=\"%s/vasthelp.html#VASTNR\">", VASTpath);
		fprintf(table, "%s</a>&nbsp;subset displayed", ThisSubset.c_str());
	    }
	    if (!i) fprintf(table, " displayed");
            if (iKept >1)
                fprintf(table, " with %d selected neighbors",iKept);
	    else if (iKept == 1)
                fprintf(table, " with 1 selected neighbor");

        }
        else if (iKept > 1)
                fprintf(table,"%d selected neighbors displayed",iKept);
        else fprintf(table, "1 selected neighbor displayed");
    }
    else {
	if (numhitsdisplayed > iKept)  {
            int numtmp;

            numtmp = numhitsdisplayed - iKept;

            if (numtmp>1)
                fprintf(table, "%d neighbors found", numtmp);
            else fprintf(table, "1 neighbor found");
            if (iKept>1)
                fprintf(table, " and displayed with %d selected ones",iKept);
            else if (iKept ==1)
                fprintf(table, " and displayed with 1 selected neighbor");

	}
        else if (iKept >1)
                fprintf(table, "%d selected neighbors displayed", iKept);
        else fprintf(table, "1 selected neighbor displayed");

        if (NonNbr && NonNbr[0] != NULLB) {
	   char str[MAX_TBUFF];
	
	   if (strchr(NonNbr, ','))
		sprintf(str, " are not structure neighbors");
	   else sprintf(str, " is not a structure neighbor");
	   if (StrLen(NonNbr) <= 41) 
	      fprintf(table,", but <font color=#CC6600>%s</font>%s",NonNbr,str);
	   else fprintf(table, ", but others are not structure neighbors");
 	}
    }
    fprintf(table, ".</strong></font><br>\n");

} /* end of VastTableBegin */



static Int2 CalCellHForCDs(Int4 mmdbid, Int4 chnno)
{
  Int4          i, j, gi, seqlen, numseg, from, to, y, maxy=0;
  Int2          CdNum, *iClus;
  Int4Ptr       starts, lens;
  OverLoc       *head, *end;
  DenseSegPtr   *dsp;
  SeqAnnotPtr   sap = NULL;
  SeqAlignPtr   salp= NULL;
  SeqEntryPtr   sep=NULL;
  BioseqPtr     bseqp=NULL;
  unsigned      *PssmId;

  gi = constructGi((unsigned)mmdbid, (unsigned)chnno);
  seqlen = constructChainLength((unsigned)mmdbid, (unsigned)chnno);

  sep = (SeqEntryPtr) constructSeqEntryForGi(gi, true, Database);
  if (!sep) return (90);
  bseqp = (BioseqPtr)sep->data.ptrvalue;

  if (!bseqp) {
        PrtMesC(MAILto, "VastSrv", "SeqEntryPtr not NULL, but BioseqPtr is NULL and gi=",  (char *)ToString(gi).c_str(), false);
        return 0;
  }

  sap = CddSynchronousQuery(bseqp, 0.01, true, true, FALSE, (char *)"", true);
  if (sap) {

    end=head = NewOverLoc((unsigned)seqlen);
    head->y = 90+FontBH;
    head->next = NULL;

    for (salp = (SeqAlignPtr)sap->data, CdNum=0; salp!=NULL;
                salp = salp->next, CdNum++);

    dsp = (DenseSegPtr *) MemNew (CdNum * sizeof(DenseSegPtr));
    PssmId = (unsigned *) MemNew (CdNum * sizeof(unsigned));
    iClus = (Int2Ptr) MemNew (CdNum * sizeof(Int2));

    for (salp= (SeqAlignPtr)sap->data, i=0; salp!=NULL; salp=salp->next,i++) {
      dsp[i] = (DenseSegPtr)salp->segs;
      PssmId[i] = GetPSSMID(dsp[i]);
      iClus[i] = -1;

      if (!PssmId[i])
        PrtMesC(MAILto, "VastSrv", "PssmId=0", "", false);
    }


    for (i=0; i< CdNum; i++) {
        if (iClus[i] >= 0) continue;

        iClus[i] = i;
        for (j=i+1; j< CdNum; j++)
             if (PssmId[i] == PssmId[j]) iClus[j] = i;
    }

    for (i=0; i< CdNum; i++) {

      numseg = dsp[i]->numseg;
      starts = dsp[i]->starts;
      lens = dsp[i]->lens;
      from = starts[0]+1;
      to = starts[(numseg-1)*2] + lens[numseg-1];
      y = GetY_nr_cddsrv(head, &end, from , to, seqlen, 5);

      y += FontBH+2;
      maxy = MAX(maxy, y);
    }

    FreeOverLoc(head);
    MemFree(dsp);
    MemFree(PssmId);
    MemFree(iClus);

  }

  if (!maxy) maxy = 90;
  return(maxy);

}       /* end of CalCellHForCDs */



static void 
VastInfoRows(FILE *table, unsigned iFSID, VastPageDataPtr vpp, 
	unsigned numhitsdisplayed, unsigned iKept,
	char * selnbrstring, char * selsdidstring, SortBy sortby, 
	SubSetID subsetnum, unsigned pagenum)
{
  unsigned 	i, cellh;

   fprintf(table, "<TABLE border=0 CELLPADDING=0 CELLSPACING=0 width=800>\n");
   fprintf(table, "                     <!-- table for checkbox -->\n");

cellh=90;

/*
   if (JobID != "") cellh = 90;
   else cellh = CalCellHForCDs(aMmdbId, aChnNo); 
*/
   fprintf(table, "<tr>\n");
   fprintf(table, "<td height=%d valign=middle><br></td>\n", cellh);
   fprintf(table, "<td width=770 rowspan=%d>\n", numhitsdisplayed+1);
   fprintf(table, "<img src=");
   fprintf(table, ParURL, URLcgi, CGIname, aSdi, sortby); 
	
   if (selnbrstring == NULL && selsdidstring == NULL) {
	if (pagenum) 
	    fprintf(table, PageSubsetURL, pagenum, subsetnum, subsetnum);
	else fprintf(table, "schsub=Find");
   }
   else {
     fprintf(table, "schsub=Find");
     if (selnbrstring) fprintf(table, "&selnbr=%s", selnbrstring);
     if (selsdidstring) fprintf(table, "&selsdid=%s", selsdidstring);
   }

   if (iKept) 
	for (i=0; i< iKept; i++) 
	    fprintf(table, "&hit=%d", getfid(vpp[i].bSdi, vpp[i].bAlignId));

   if (JobID != "") {
	if (ReqId) fprintf(table, "&reqid=%llu", ReqId);
	else fprintf(table, "&vsid=%s&pass=%s", JobID.c_str(), Passwd.c_str());
   }

   fprintf(table,"&cmd=graph&imgsize=%d\" usemap=#chain_map border=0 ismap>",
		(numhitsdisplayed*30 + cellh));
   fprintf(table, "\n</td>\n</tr>\n\n");

   fprintf(table, "<map name=chain_map>\n");

   ImgMapOrDraw(true, iFSID, vpp, numhitsdisplayed, iKept, selnbrstring, 
	selsdidstring,	sortby, subsetnum, pagenum, cellh, table);
   fprintf(table, "</map>\n");

  for (i = 0; i < numhitsdisplayed; i++) {

    fprintf(table, "<tr>\n");
    fprintf(table, "<td VALIGN=TOP height=30 width=30 class=SMALL1>\n");
    fprintf(table, "<INPUT TYPE=checkbox NAME=hit VALUE=%d",
					 getfid(vpp[i].bSdi, vpp[i].bAlignId));
    if (i < iKept) fprintf(table, " checked");
    fprintf(table, "></td>\n");
    fprintf(table, "</tr>\n\n");

    fflush(table);

  } /* end of for */

  fprintf(table, "</TABLE>\n\n");
  fprintf(table, "</FORM>\n"); 	      /* close Form opened in PrintAlignView */

} /* end of VastInfoRows */



#define SetTd	"<td VALIGN=top ALIGN=center>"
#define SetFont "<font color=#CC6600>"

static void 
VastTableRows(FILE *table, VastPageDataPtr vpp, 
	unsigned numhitsdisplayed, unsigned iKept, SortBy sortby)
{
  unsigned 	bMMDBid;
  Uint1 domNo;
  char 	pcSlaveName[2 * DBStrSize + 4];
  float f;
  int i;

   fprintf(table, "<br>\n");
   fprintf(table, "<TABLE cellspacing=3 cellpadding=2 width=800 border=1>\n");
   fprintf(table,"<tr valign=middle>\n");
   fprintf(table,"<th>&nbsp;</th>\n");
   fprintf(table,"<th align=left><pre> <a href=\"%s/vasthelp.html#VASTTable\">PDB</a>", VASTpath);
   fprintf(table," <a href=\"%s/vasthelp.html#VASTTable\">C</a>", VASTpath);
   fprintf(table," <a href=\"%s/vasthelp.html#VASTTable\">D</a></pre></th>\n",VASTpath);

   fprintf(table,"<th><pre><a href=\"%s/vasthelp.html#VASTTable\">", 
		VASTpath);
   if (sortby == NRes) 
	fprintf(table, "%sAli. Len.</font></a></pre></th>\n", SetFont);
   else fprintf(table, "Ali. Len.</a></pre></th>\n");

   fprintf(table,"<th><pre><a href=\"%s/vasthelp.html#VASTTable\">", 
		VASTpath);
   if (sortby == VScore) 
	fprintf(table, "%sSCORE</font></a></pre></th>\n", SetFont);
   else fprintf(table, "SCORE</a></pre></th>\n");

   fprintf(table,"<th><pre><a href=\"%s/vasthelp.html#VASTTable\">", 
		VASTpath);
   if (sortby == PValue)
	fprintf(table, "%sP-VAL</font></a></pre></th>\n", SetFont);
   else fprintf(table, "P-VAL</a></pre></th>\n");

   fprintf(table,"<th><pre><a href=\"%s/vasthelp.html#VASTTable\">", 
		VASTpath);
   if (sortby == Rmsd) 
	fprintf(table, "%sRMSD</font></a></pre></th>\n", SetFont);
   else fprintf(table, "RMSD</a></pre></th>\n");

   fprintf(table,"<th><pre><a href=\"%s/vasthelp.html#VASTTable\">", 
		VASTpath);
   if (sortby == PcntId)
	fprintf(table, "%s%%Id</font></a></pre></th>\n", SetFont);
   else fprintf(table, "%%Id</a></pre></th>\n");

   fprintf(table,"<th><pre><a href=\"%s/vasthelp.html#VASTTable\">Description</pre></th>\n",VASTpath);
   fprintf(table,"</tr><br>\n");
   fflush(table);


  for (i = 0; i < numhitsdisplayed; i++) {
      
    char	pdbId[5];

    bMMDBid = SdiToMmdbId(vpp[i].bSdi);
    domNo = (unsigned) atoi(vpp[i].bDomName+5);

    pcSlaveName[0] = NULLB;
    if(vpp[i].bDescr1[0] != NULLB || vpp[i].bDescr2[0] != NULLB) {
      if(vpp[i].bDescr1[0] != NULLB) strcpy(pcSlaveName,vpp[i].bDescr1);
      if(vpp[i].bDescr2[0] != NULLB) strcat(pcSlaveName,vpp[i].bDescr2);
    }
    else if (JobID != "") 
	constructPdbDescr(bMMDBid, pcSlaveName, 2 * DBStrSize + 4);

    fprintf(table, "<tr>\n");
    fprintf(table, "<td VALIGN=TOP>"); 

    fprintf(table, "<INPUT TYPE=checkbox NAME=hit VALUE=%d", 
				getfid(vpp[i].bSdi, vpp[i].bAlignId));
    if (i < iKept) fprintf(table, " checked");
    fprintf(table, "></td>\n");
    fprintf(table, "<td VALIGN=TOP><pre>");
    StrCut(pdbId, vpp[i].bDomName, 1, 4);
    fprintf(table, "<a href=\"%s%s?uid=%ld&Dopt=s\">%s</a>", 
		URLcgi, MMDBCGIname, (long) bMMDBid, pdbId);

    if (vpp[i].bDomName[4] == ' ') fprintf(table,"&nbsp;");
    else
      fprintf(table," <a href=\"%s%s?uid=%ld&Dopt=s\">%c</a>", 
		URLcgi, MMDBCGIname, (long) bMMDBid, vpp[i].bDomName[4]);

    if (domNo > 0)
      fprintf(table," <a href=\"%s%s?uid=%ld&Dopt=s\">%d</a></pre></td>\n", 
		URLcgi, MMDBCGIname, (long) bMMDBid, domNo);
    else
      fprintf(table," </pre></td>\n");

    if (vpp[i].nres > 0) {
      	fprintf(table,"%s", SetTd);
	if (sortby == NRes)
	    fprintf(table, "%s%d</font></td>\n", SetFont, vpp[i].nres);
	else fprintf(table, "%d</td>\n", vpp[i].nres);
    }
    else fprintf(table,"<td> </td>\n");

    if (vpp[i].vScore > 0) {
        f = (FloatLo) (vpp[i].vScore);
        f = f/(FloatLo) ASP_SCALE_FACTOR;
        fprintf(table, "%s", SetTd);
	if (sortby == VScore)
	    fprintf(table, "%s%.1f</font></td>\n",SetFont, f);
	else fprintf(table, "%0.1f</td>\n", f);
    }
    else fprintf(table, "<td> </td>\n");

    if (vpp[i].mlogp > 0) {
        f = (float) (vpp[i].mlogp);
        f = f/(float) ASP_SCALE_FACTOR;

        /* adjust for database size */
        f -= LOG10_500;

        if (f <= 4.0) {
          f = (float) exp(-LOG_10*f);
          fprintf(table, "%s", SetTd);
	  if (sortby == PValue)
	      fprintf(table, "%s%.4f</font></td>\n", SetFont, f);
	  else fprintf(table, "%.4f</td>\n", f);
        }
        else {
          fprintf(table, "%s", SetTd);
	  if (sortby == PValue) 
	      fprintf(table, "%s10e-%.1f</font></td>\n", SetFont, f);
	  else fprintf(table, "10e-%.1f</td>\n", f);
	}
    }
    else fprintf(table,"<td> </td>\n");

    if (vpp[i].rmsd >= 0) {
      f = (FloatLo) (vpp[i].rmsd);
      f = f/(FloatLo) ASP_SCALE_FACTOR;
      fprintf(table, "%s", SetTd);
      if (sortby == Rmsd) 
	fprintf(table, "%s%.1f</font></td>\n", SetFont, f);
      else fprintf(table,"%.1f</td>\n", f);
    }
    else fprintf(table,"<td> </td>\n");  

    if (vpp[i].pcntId > 0) {
      f = (FloatLo) (vpp[i].pcntId);
      f = f/(FloatLo) ASP_SCALE_FACTOR;
      fprintf(table, "%s", SetTd);
      if (sortby == PcntId) 
      	  fprintf(table, "%s%.1f</font></td>\n", SetFont, f * 100.0); 
      else fprintf(table, "%.1f</td>\n", f*100.0);
    }
    else
      fprintf(table,"<td VALIGN=top ALIGN=center>0.0</td>\n");  

    fprintf(table,"<td VALIGN=top>\n");
    fprintf(table,"%s" , pcSlaveName);
    fprintf(table,"<BR>\n");
    fprintf(table,"</td>\n</tr>\n");
    fflush(table);
  } /* end of for */


  fprintf(table, "</TABLE>\n\n");
  fprintf(table, "</FORM>\n"); 		/* close Form opened in PrintAlignView */
} /* end of VastTableRows */



static void 
SaveBioAnnotSetInFile(unsigned Fsid, SubSetID subsetnum, SortBy sortby)
{
        unsigned                  	numhitsdisplayed, alloc_size;
        BiostrucAnnotSetPtr     pbsa = NULL;
        AsnIoPtr                aip;
	VastPageDataPtr 	vpp=NULL;

/*
  	if (JobID != "") {
	    OpenMMDBAPI((POWER_VIEW), NULL);
	    pbsa = LocalGetFeatureSet(Fsid/10000, Fsid);
	    CloseMMDBAPI();
	}
	else
*/
	{
        	vpp = constructVastPagesForAllNbrs(aSdi, subsetnum, sortby, 
		      (unsigned *)&alloc_size, (unsigned *)&numhitsdisplayed);
        	pbsa = constructBASPFromVastPagePtr(vpp, numhitsdisplayed, 1);
	}

        printf("Content-type: application/octet-stream\r\n\r\n"); 
/*
        printf("Content-type: text/html\r\n\r\n");
        printf("<HTML><PRE>\r\n");
*/
        aip = AsnIoNew(ASNIO_TEXT_OUT, stdout, NULL, NULL, NULL);
        BiostrucAnnotSetAsnWrite(pbsa, aip, NULL);
        AsnIoFlush(aip);
        AsnIoClose(aip);
	free(vpp);
}



static void
MakeVastTableGraph(VastPageDataPtr vpp, unsigned numhitsdisplayed,unsigned FSID,
	unsigned iKept, char * selnbrstring, 
	char * selsdidstring, char * NonNbr, SortBy sortby,SubSetID subsetnum,
        unsigned pagenum, unsigned numhits, char cTable, unsigned NbrFlag,
        short ImgSize, short GraphFlag)
{
  unsigned numpages;

  numpages = numhits / NUM_HITS_PER_PAGE;
  numpages = numhits % NUM_HITS_PER_PAGE == 0 ? numpages : (numpages+1);

  if (NbrFlag == DEF_ALL_NBR) {
     if (numhits < NUM_HITS_PER_PAGE) pagenum = DEFAULT_PAGE;
  }
  else pagenum = 0;


  if (GraphFlag) 
      	DrawStrucNbr(FSID, vpp, numhitsdisplayed, iKept, selnbrstring, 
	     selsdidstring, sortby, subsetnum, pagenum, ImgSize);
  else {
  	VastTableBegin(stdout, FSID, sortby, subsetnum, iKept,NonNbr,
		vpp, numhitsdisplayed, numhits, numpages, 
		pagenum, cTable, NbrFlag);
	switch (cTable) {
	    case 'n': 

		     fprintf(stdout, "<b><font class=SMALL1>In the graphics below the ");
		     fprintf(stdout, "<font color=red>red</font> regions are aligned segments, please click them to get a correspounding sequence alignment.<b></font>\n");
		     fprintf(stdout, "<br><br>\n");	
		     VastInfoRows(stdout, FSID, vpp, numhitsdisplayed, iKept, 
			selnbrstring, 
			selsdidstring, sortby, subsetnum, pagenum);
		     break;
	    case 'y':
		     VastTableRows(stdout, vpp, numhitsdisplayed, iKept, 
				sortby);
		     break;
	    default:
		     PrtMesC("", "VASTSERV", 
                        "Please select \"Graph\" or \"Tabler\" for the display on screen or \"File\" for a local copy.",
			"", false);
	}

	fprintf(stdout,"\n<BR><align=left HR width=800 SIZE=5 NOSHADE><BR>\n");
        WWWPrintFileData((char *)TAILFILE,  stdout);

  }

  /* RemoveTempFiles();   */

  return;       

} /* MakeVastTableGraph() */




/* Extract vastsrv parameters from the config file. */
static bool GetVastParams()
{
  URLBase[0] = URLcgi[0] = ENTREZurl[0] = DOCSUMurl[0] = HELPname[0]
  = CGIname[0] = MMDBCGIname[0] = DART[0] = LOGIN[0] = PASSWD[0] 
  = INITpath[0] = LIBpath[0] = CDDurl[0] = Database[0] = NULLB;

  GetAppParam("vsch", "VASTSRV", "URLBase", "", URLBase, PATH_MAX);
  if (URLBase[0] == NULLB) {
    ErrPostEx(SEV_FATAL,0,0,
		"VAST config file\nVASTSRV section has no URLBase...\n");
    return FALSE;
  }

  GetAppParam("vsch", "VASTSRV", "URLcgi", "", URLcgi, PATH_MAX);
  if (URLcgi[0] == NULLB) {
    ErrPostEx(SEV_FATAL,0,0,
			"VAST config file\nVASTSRV section has no URLcgi...\n");
    return FALSE;
  }

  GetAppParam("vsch", "VASTSRV", "ENTREZurl", "", ENTREZurl, PATH_MAX);
  if (ENTREZurl[0] == NULLB) {
    ErrPostEx(SEV_FATAL,0,0,
		"VAST config file\nVASTSRV section has no ENTREZurl...\n");
    return FALSE;
  }

  GetAppParam("vsch", "VASTSRV", "DOCSUMurl", "", DOCSUMurl, PATH_MAX);
  if (DOCSUMurl[0] == NULLB) {
    ErrPostEx(SEV_FATAL,0,0,
		"VAST config file\nVASTSRV section has no DOCSUMurl...\n");
    return FALSE;
  }

  GetAppParam("vsch", "VASTSRV", "CGIname", "", CGIname, PATH_MAX);
  if (CGIname[0] == NULLB) {
    ErrPostEx(SEV_FATAL,0,0,
		"VAST config file\nVASTSRV section has no CGIname...\n");
    return FALSE;
  }

  GetAppParam("vsch", "VASTSRV", "MMDBCGIname", "", MMDBCGIname, PATH_MAX);
  if (MMDBCGIname[0] == NULLB) {
    ErrPostEx(SEV_FATAL,0,0,
		"VAST config file\nVASTSRV section has no MMDBCGIname...\n");
    return FALSE;
  }

  GetAppParam("vsch", "VASTSRV", "HELPname", "", HELPname, PATH_MAX);
  if (HELPname[0] == NULLB) {
    ErrPostEx(SEV_FATAL,0,0,
		"VAST config file\nVASTSRV section has no HELPname...\n");
    return FALSE;
  }

char mailto[PATH_MAX];
  GetAppParam("vsch", "VASTSRV", "MAILto", "", mailto, PATH_MAX);
  MAILto = string(mailto);
  if (MAILto == "") {
    ErrPostEx(SEV_FATAL,0,0,
		"VAST config file\nVASTSRV section has no MAILto...\n");
    return FALSE;
  }

  GetAppParam("vsch", "VASTSRV", "VSPATH", "", VSPATH, PATH_MAX);
  if (VSPATH[0] == NULLB) {
    ErrPostEx(SEV_FATAL,0,0,
	"VAST config file\nVASTSRV section has no VSPATH...\n");
    return FALSE;
  }

  GetAppParam("vsch", "VASTSRV", "DATApath", "", DATApath, PATH_MAX);
  if (DATApath[0] == NULLB) {
    ErrPostEx(SEV_FATAL, 0, 0, 
		"VAST config file\nVASTSRV section has no VAST DATApath...\n");
    return FALSE;
  }

  GetAppParam("vsch", "VASTSRV", "VASTpath", "", VASTpath, PATH_MAX);
  if (VASTpath[0] == NULLB) {
    ErrPostEx(SEV_FATAL, 0, 0, 
		"VAST config file\nVASTSRV section has no VASTpath...\n");
    return FALSE;
  }

  GetAppParam("mmdb", "MMDB", "Database", "", Database, PATH_MAX);
  if (Database[0] == NULLB) {
      ErrPostEx(SEV_FATAL,0,0, "MMDB config file \nMMDB section no Database\n");      return FALSE;
  }

  GetAppParam("mmdb", "CD", "DARTdb", "", DART, PATH_MAX);
  if (DART[0] == NULLB) {
      ErrPostEx(SEV_FATAL,0,0, "MMDB config file \nCD section no DARTdb\n");
      return FALSE;
  }

  GetAppParam("mmdb", "CD", "LOGINname", "", LOGIN, PATH_MAX);
  if (LOGIN[0] == NULLB) {
      ErrPostEx(SEV_FATAL,0,0, "MMDB config file \nCD section no LOGINname\n");
      return FALSE;
  }

  GetAppParam("mmdb", "CD", "PASSwd", "", PASSWD, PATH_MAX);
  if (PASSWD[0] == NULLB) {
      ErrPostEx(SEV_FATAL,0,0, "MMDB config file \nCD section no PASSwd\n");
      return FALSE;
  }

  GetAppParam("mmdb", "CD", "INITpath", "", INITpath, PATH_MAX);
  if (INITpath[0] == NULLB) {
      ErrPostEx(SEV_FATAL,0,0, "MMDB config file \nCD section no INITpath\n");
      return FALSE;
  }

  GetAppParam("mmdb", "CD", "LIBpath", "", LIBpath, PATH_MAX);
  if (LIBpath[0] == NULLB) {
      ErrPostEx(SEV_FATAL,0,0, "MMDB config file \nCD section no LIBpath\n");
      return FALSE;
  }

  GetAppParam("mmdb", "CD", "CDDurl", "", CDDurl, PATH_MAX);
  if (CDDurl[0] == NULLB) {
      ErrPostEx(SEV_FATAL,0,0, "MMDB config file \nCD section no CDDurl\n");
      return FALSE;
  }

  return true;

} /* GetVastParams() */




/* The following ValidateMMDBID() was changed by J. Chen, May, 2002) */

static bool ValidateMMDBID(char * pdbcode, unsigned iMMDBid)
{
  DocUid 	uid;
  Uint1		i;

  for (i=0; i< StrLen(pdbcode); i++) pdbcode[i] = toupper(pdbcode[i]); 
  uid = constructLiveOrDeadMmdbIdForPdbId(pdbcode, (bool *)&psok, Database);
  
  if ((unsigned)uid == iMMDBid) return true;
  else return false;

}



static BiostrucFeaturePtr FilterHitsByDomainSubset(BiostrucFeaturePtr pbsf, 
	unsigned subsetnum, unsigned *KepId, unsigned iKept)
{
  BiostrucFeaturePtr current, pbsfHead = NULL, pbsfTail;
  unsigned i, n, dombsfid;
  unsigned gn, gr, hcnt, *min_ranks, *group_num, *group_rank;
  char domid[DOMID_SIZE + 1], pdbcode[4+1];
 
/* The next bit of code is used for filtering the hit lists.  When we go through a hit
 * list we skip domains that do not belong to the subset of interest, or which belong to
 * the subset but for which a group representative has already been encountered.
*/
  for (i = 0; i <= DOMID_SIZE; i++)
    domid[i] = NULLB;

// block to pass the C++ comiler on linux. 3/2/05
//   n = GetNumberOfDomains();
  min_ranks = new unsigned [n];
  group_num = new unsigned [n];
  group_rank = new unsigned [n];

  /* use a first pass to "flag" the selected hits in the specified subset 
     If group_num[i] = 0 or group_rank[i] is larger than min_ranks[i], 
     then the ith neighbor should not be included in summary page */

  if (group_rank != NULL) {	      /* check mem alloc */

    for (i = 0; i < n; i++) {	      /* initializations */
         group_num[i] = 0;
         min_ranks[i] = n + 1;
         group_rank[i] = n + 1;
    }

    for (current = pbsf, hcnt = 0; current != NULL; current = current->next, hcnt++) {
       /* copy domain identifier into domid[] */
      if (current->name[6] == ' '){
	 StrCut(domid, current->name, 8, 13);
	 if (current->name[13] != ' ') domid[6] = current->name[13];
      }
      else {
	StrCut(domid, current->name, 9, 14);
	if (current->name[14] != ' ') domid[6] = current->name[14];
      }
      
      if (domid[5] == '0') domid[5] = ' ';

       /* if not in subset then skip over this domain */
// block to pass the C++ compiler on linux. 3/2/05
//      if (BelongsToSubset(domid, subsetnum, &gn, &gr) <= 0) {
//	continue;
 //     }

      /* otherwise record group data for this hit */
      group_num[hcnt] = gn;
      group_rank[hcnt] = gr;

      /* and reset minimum rank for this group */
      if (gr < min_ranks[gn - 1])
	min_ranks[gn - 1] = gr;
    }
  }

  /* Now use the values just set in group_num, group_rank, and min_ranks
     to decide whether or not the current neighbor should be linked
     into the new feature list */

  current = pbsf;
  hcnt = 0;
  while (current) {
    dombsfid = current->id;

    if (group_rank != NULL)		      /* check mem alloc */
    {
      if (group_num[hcnt] == 0)	     /* group_num not set so NOT in subset */
      {					      /* incr hcnt and do NOT link */
        hcnt++;
	current = current->next;
        continue;
      }
   
      gn = group_num[hcnt];

      if (group_rank[hcnt] != min_ranks[gn - 1])  
			/* neighbor is in subset but of lower rank */
      {
        hcnt++;			         /* incr hcnt and do NOT link */
	current = current->next;
        continue;
      }
    }
    
    /* With the new rcsb depositions we now need to validate mmdbids
       in the .bas files. Extract pdbcode and use it with mmdbID to validate */

    StrCut(pdbcode, current->name, 8, 11);
    pdbcode[4]= NULLB;
       
    if (ValidateMMDBID(pdbcode, SdiToMmdbId((current->id)/100)) == false) {
      hcnt++;
      current = current->next;	         /* incr hcnt and do NOT link */
      continue;
    }

/* remove duplicated lists, J. Chen */
    i=0; 
    while (i< iKept && dombsfid != KepId[i]) i++; 
    if (i < iKept) {
	hcnt++;
	current = current->next;
	continue;
    }

    if (pbsfHead == NULL)    /* neighbor has passed all tests - include it!*/
    {
      pbsfHead = current;
      pbsfTail = pbsfHead;
      hcnt++;
      current = current->next;
      pbsfTail->next = NULL;
    }
    else
    {
      pbsfTail->next = current;
      pbsfTail = current;
      hcnt++;
      current = current->next;
      pbsfTail->next = NULL;
    }
  }
  
  delete [] min_ranks;
  delete [] group_num;
  delete [] group_rank;
  return pbsfHead;

}	/* FilterHitsByDomainSubset */

 

static void FilterHitsByPageNew(BiostrucFeatureSetPtr pbsfs, unsigned PageNum, 
	unsigned HitsPerPage, unsigned *numhits, unsigned *numpages)
{
  
  BiostrucFeaturePtr pbsf, pbfHead=NULL, pbfTail=NULL;
  unsigned index, FidCount, RemainFids, CompleteFidSet; 
  unsigned UpperLimit, LowerLimit;
 
  pbsf = pbsfs->features;
  if (pbsf == NULL) {
    *numpages = 1;
    *numhits = 0;
  }
  
  FidCount = 0;
  
  while (pbsf)
  {
    FidCount++;
    pbsf = pbsf->next;
  }
  *numhits = FidCount;

  if (FidCount <= HitsPerPage) *numpages = 1;
  else
  { 
    RemainFids = FidCount % HitsPerPage;
  
    if (RemainFids)
    {
     CompleteFidSet = FidCount - RemainFids;
     *numpages = (CompleteFidSet/HitsPerPage) + 1;
    }
    else *numpages = FidCount/HitsPerPage;
  }
  
  UpperLimit = HitsPerPage * PageNum;
  LowerLimit = UpperLimit - HitsPerPage + 1;
    
  if ((FidCount < HitsPerPage ) || (LowerLimit > FidCount)) {
    UpperLimit = FidCount;
    LowerLimit = 1;
  }
  
  if (UpperLimit > FidCount) UpperLimit = FidCount;
  
  for (index = 1, pbsf = pbsfs->features; index <= FidCount; 
	index++, pbsf=pbsf->next)
  {
    if ((index >= LowerLimit) && (index <= UpperLimit)) {
	if (!pbfHead) pbfHead = pbfTail = pbsf;
	else {
		pbfTail->next = pbsf;
		pbfTail = pbsf;
	}
    }
  }
  
  if (pbfTail) pbfTail->next = NULL;
  pbsfs->features = pbfHead;

}	/* FilterHitsByPageNew */


static 
BiostrucAnnotSetPtr GetValidBiostrucAnnotSet(BiostrucAnnotSetPtr basp,unsigned Fsid)
{
    BiostrucAnnotSetPtr basp2 = NULL;
    BiostrucFeatureSetPtr pbsfs = NULL, pbsfsLast= NULL;
    BiostrucIdPtr               pbsidThis = NULL;

    if (basp == NULL) return NULL;

    pbsfs = basp->features;
    while (pbsfs) {
       if (pbsfs->id == Fsid)
        {
          basp2 = BiostrucAnnotSetNew();
          basp2->id = basp->id;
          basp->id = NULL; /* unlink the id valnode from basp object */
          basp2->descr = basp->descr;
          basp->descr = NULL;  /* unlink the descr from basp object */
          basp2->features = pbsfs;
          if (pbsfsLast) /* relink next to prev */
            pbsfsLast->next = pbsfs->next;
          else
            basp->features = pbsfs->next;
          basp2->features->next = NULL;
	  break;
	}
       pbsfsLast = pbsfs;
       pbsfs = pbsfs->next;
    }
    BiostrucAnnotSetFree(basp);
    
    if (basp2 == NULL) {
	char    str[15];

        sprintf(str, "%d", Fsid);
	PrtMesC("","VASTSRV","Incorrect sdid or chaindom: chaindom = ",str,false);
    }

    pbsidThis = ValNodeFindNext(basp2->id, NULL, BiostrucId_mmdb_id);
    if (!pbsidThis) {
	PrtMesC(MAILto, "VASTSRV", "No MMDB-ID in Data on Server.", "", true);
	BiostrucAnnotSetFree(basp2);
	return NULL;
    }

    return(basp2);

} /* GetValidBiostrucAnnotSet */


static unsigned 
MakeVppByVS(unsigned Fsid, VastPageDataPtr vpp, SubSetID subsetnum,
 	SortBy sortby, unsigned pagenum, unsigned *Keptfid, unsigned iKept)
{

    BiostrucAnnotSetPtr 	basp = NULL, basp2 = NULL;
    BiostrucFeatureSetPtr 	pbsfs = NULL, pbsfs2 = NULL;
    BiostrucIdPtr 	  	pbsidThis = NULL;
    unsigned			numhits, numpages, numhitsdisplayed;
    char                        str[15];

    basp = LocalGetBiostrucAnnotSet(Fsid/10000);
    if (!basp) {
	sprintf(str, "%d", Fsid);
        PrtMesC("","VASTSRV","Incorrect sdid or chaindom: chaindom = ",str,false);
    }

    basp2 = GetValidBiostrucAnnotSet(basp, Fsid);
if (!basp2) fprintf(stderr, "basp2==NULL\n");
    if (!basp2) return 0;

    pbsfs = basp2->features;
    while (pbsfs)
    {
     if (pbsfs->id == Fsid)
      {

	pbsidThis = ValNodeFindNext(pbsfs->descr, NULL,
					BiostrucFeatureSetDescr_name);
        if (pbsfs->features != NULL) {
          pbsfs2 = BiostrucFeatureSetNew();
          pbsfs2->id = pbsfs->id;
          pbsfs->id = NULL;
          pbsfs2->descr = pbsfs->descr;
          pbsfs->descr = NULL;

          pbsfs2->features =
		FilterHitsByDomainSubset(pbsfs->features, subsetnum, 
			Keptfid, iKept);
          pbsfs->features = NULL;
          if (pbsfs2->features != NULL) VastTableSort(pbsfs2, sortby);
	  else 
	    PrtMesC("chenj@ncbi.nlm.nih.gov", "VASTSrv",
		"Error: pbsfs2->features == NULL\n", "", false);

          FilterHitsByPageNew(pbsfs2, pagenum, NUM_HITS_PER_PAGE, 
					&numhits, &numpages);
          if (numhits < NUM_HITS_PER_PAGE) pagenum = DEFAULT_PAGE;
	  
	  numhitsdisplayed = BFSP2VastPageDataPtr(pbsfs2, vpp, numhits);
        }
        else numhitsdisplayed = 0;

	BiostrucAnnotSetFree(basp2);
	return (numhitsdisplayed);
      }
     pbsfs = pbsfs->next;
   }

  BiostrucAnnotSetFree(basp2);
  sprintf(str, "%d", Fsid);
  PrtMesC("", "VASTSRV", "Incorrect sdid or chaindom: chaindom = ", str,false);

  return (-1);

}	/* MakeVppByVS */




static BiostrucFeaturePtr
FilterHitsByNbrs(BiostrucFeaturePtr bfp, char * *SelNames, unsigned *SelId, 
	unsigned issdid, unsigned iSel, unsigned *KepId, unsigned iKept)
{
   BiostrucFeaturePtr pbfHead = NULL, pbfTail=NULL;
   char		*str, sName[10];
   unsigned 	i, j, id, sdi;

   if (bfp == NULL) return NULL;
   
   while (bfp) {
     if (SelNames != NULL) {
	str = bfp->name;
     	if (str[6] == ' ') {
	    StrCut(sName, str, 8, 13); 
	    if (StrLen(str) >= 14 && str[13] != ' ') StringCat(sName, str+13);
        }
     	else {
	    StrCut(sName, str, 9, 14);
	    if (StrLen(str) >=15 && str[14] != ' ') StringCat(sName, str+14);
     	}
     }
     id = bfp->id;
     sdi = DomNameToSdi(sName);

     for (i=0; i< iSel; i++) {
	if (SelNames != NULL) {
     	    if (StringNCmp(SelNames[i], sName, StrLen(SelNames[i]))) continue;
	}
	else {
	     if (issdid) {
		if (sdi != SelId[i]) continue;
	     }
	     else if (id != SelId[i]) continue;
  	}

	j=0; 
	while (j< iKept && id != KepId[j]) j++;
        if (j < iKept) continue;
        if (pbfHead == NULL)  {
		pbfHead = pbfTail = bfp;
	}
        else {
		pbfTail->next = bfp;
		pbfTail = bfp;
	}
	break; 
    }
    bfp = bfp->next;

  }
		
  if (pbfTail) pbfTail->next = NULL;

  return pbfHead;
} 	/* end of FilterHitsByNbrs */
					



static unsigned MakeVppByVSNbr(unsigned Fsid, VastPageDataPtr vpp, 
  	char * *SelNames, unsigned *SelId, unsigned issdid, unsigned iSel, unsigned *KepId, 		unsigned iKept, SortBy sortby)
{
    BiostrucAnnotSetPtr 	basp = NULL, basp2=NULL;
    BiostrucFeatureSetPtr 	pbsfs = NULL, pbsfs2 = NULL;
    unsigned 			numhitsdisplayed;
    char                        str[15];

    basp = LocalGetBiostrucAnnotSet(Fsid/10000);
    if (!basp) {
        sprintf(str, "%d", Fsid);
        PrtMesC("","VASTSRV","Incorrect sdid or chaindom. chaindom  = ",str,false);
    }

    basp2 = GetValidBiostrucAnnotSet(basp, Fsid);
    pbsfs = basp2->features;
    while (pbsfs) {
	if (pbsfs->id == Fsid) {
	    if (pbsfs->features != NULL) {
		pbsfs2 = BiostrucFeatureSetNew();
		pbsfs2->id = pbsfs->id;
		pbsfs->id = NULL;
		pbsfs2->descr = pbsfs->descr;
		pbsfs->descr = NULL;

		pbsfs2->features = 
		   FilterHitsByNbrs(pbsfs->features, SelNames, SelId, issdid, 
				iSel, KepId, iKept);
		pbsfs->features = NULL;
		if (pbsfs2->features) {
		   if (sortby) VastTableSort(pbsfs2, sortby);
		
		   numhitsdisplayed= 
			BFSP2VastPageDataPtr(pbsfs2, vpp, MaxVPP);
	       } else numhitsdisplayed = 0;
	    }
	    else numhitsdisplayed = 0;

	    BiostrucAnnotSetFree(basp2);
	    return (numhitsdisplayed);
	}
	pbsfs = pbsfs->next;
    }

    BiostrucAnnotSetFree(basp2);
    sprintf(str, "%d", Fsid);
    PrtMesC("", "VASTSRV", "Incorrect sdid or chaindom: chaindom = ", str, false);

    return (-1);

} /* MakeVppByVSNbrs */




static void CheckNbrs(char * arg, char * *SelNames, unsigned iSel, VastPageDataPtr vpp, unsigned vppnum)
{
   char		savedchar, name[10];
   char 	*str, *ptr0=NULL, *ptr=NULL, *ptr2=NULL;
   unsigned	i, j, count=0;
  
   arg[0] = NULLB;
   if (!vppnum) {
	for (i=0; i< iSel; i++) {
	   char tmp[10];	

           sprintf(tmp, "%s, ", SelNames[i]);
	   StringCat(arg, tmp);
 	}
	return;
   }

   str  = new char [vppnum*8+1]; 
   ptr0 = new char [vppnum*8+1];

   str[0] = NULLB;
   for (i=0; i< vppnum; i++) {
	StringCat(str, vpp[i].bDomName);
	StringCat(str, ",");
   }

   for (i=0; i< iSel; i++) {
     sprintf(ptr0, str);
     ptr = ptr0;
     while (*ptr) {
	sprintf(name, SelNames[i]);
	if (IsMmdbId(name)) { 
	    j = (unsigned) atol (name);
	    constructPdbIdForMmdbId(j, name);
	    name[4] = NULLB;
	} 
	ptr = SkipToSet(ptr, name);
	ptr2 = SkipToSet(ptr, (char *)",");
	savedchar = *ptr2;
	*ptr2 = NULLB;
 	if (!strncmp(name, ptr, StrLen(name))) break;
	*ptr2 = savedchar;
	ptr = SkipSet(ptr2, (char *)",");
     }
     if (ptr[0] == NULLB) {
	char	tmp[10];
	
	sprintf(tmp, "%s, ", SelNames[i]);
	StringCat(arg, tmp);
	count++;
     }
  }

  delete [] str; 
  delete [] ptr0;
  return;
}




void CheckNbrsbysdid(char * NonNbr, unsigned *SelSds, unsigned iSds,
        VastPageDataPtr vpp, unsigned numhits)
{
   unsigned  i, j;
   char  str[10];

   NonNbr[0] = NULLB;

   for (i=0; i< iSds; i++) {
        j=0;
        while (j < numhits && SelSds[i] != vpp[j].bSdi) j++;
        if (j == numhits) {
            sprintf(str, "%d, ", SelSds[i]);
            StringCat(NonNbr, str);
        }
   }
   

}       /* CheckNbrsbysdid */



static void SaveBioAnnotSetInFile_DBVS(SubSetID subsetnum, SortBy sortby)
{
        unsigned 	numhitsdisplayed, alloc_size;
        BiostrucAnnotSetPtr     pbsa = NULL;
        AsnIoPtr                aip;
        VastPageDataPtr         vpp=NULL;

        if (JobID != "") 
            vpp = GetVSVppForAllNbrs(aSdi, subsetnum, sortby,&numhitsdisplayed);
        else vpp = constructVastPagesForAllNbrs(aSdi, subsetnum, sortby,
                        (unsigned *)&alloc_size, (unsigned *)&numhitsdisplayed);

	pbsa = constructBASPFromVastPagePtr(vpp, numhitsdisplayed, 1);
        printf("Content-type: application/octet-stream\r\n\r\n");
/*
        printf("Content-type: text/html\r\n\r\n");
        printf("<HTML><PRE>\r\n");
*/
        aip = AsnIoNew(ASNIO_TEXT_OUT, stdout, NULL, NULL, NULL);
        BiostrucAnnotSetAsnWrite(pbsa, aip, NULL);
        AsnIoFlush(aip);
        AsnIoClose(aip);
        free(vpp);

} /* SaveBioAnnotSetInFile_newVS */




short Main()
{
  char 		cTable='n', ODBCInitStr[100], LIBPathStr[100], str[15];
  char          NonNbr[MAX_TBUFF];
  char  	*www_arg, *Cmd=NULL, **SelNames=NULL;
  char 		*selnbrstring = NULL, *selsdidstring=NULL;
  short		ret, ImgSize;
  int		indx;
  unsigned 	aDomId=0, aDomNo=0, pagenum;
  unsigned	numhitsdisplayed, totalnumhits;
  unsigned	maxAsdi = 0, iSel=0, iKept=0, i, *Keptfid;
  unsigned	NbrFlag, iSubBut=DEF_LIST_SUBMIT, *SelSds, iSds;
  struct rlimit rl;
  SubSetID 	pre_subsetnum=DEFAULT_SUBSET_NUM; 
  SubSetID	subsetnum=DEFAULT_SUBSET_NUM;
  SortBy 	sortby = DEFAULT_SORT_BY;
  WWWInfoPtr 	www_info;
  VastPageData 	vpp[MaxVPP];


  /* this sets up the unix time limit */
  getrlimit(RLIMIT_CPU, &rl);
  rl.rlim_max = rl.rlim_cur = CPUTIME_MAX;
  setrlimit(RLIMIT_CPU, &rl);

/* Begin processing www information block */
  if (WWWGetArgs(&www_info) != WWWErrOk) 
	PrtMesC("", "VASTSRV", "Failed to process posting - check your get/post syntax.", "", false);

  if (WWWGetNumEntries(www_info) == 0) 
	PrtMesC("", "VASTSERV", "No input - nothing to report.", "", false);

  if (GetVastParams() == false)
        PrtMesC(MAILto, "VASTSRV","Couldn't read the config file \".vschrc\"...", "", false);

  sprintf(ODBCInitStr, "ODBCINI=%s", INITpath);
  putenv(ODBCInitStr);
  sprintf(LIBPathStr, "LD_LIBRARY_PATH=%s", LIBpath);
  putenv(LIBPathStr);


/*
  if (!Dart_InitReal())
        PrtMesC(MAILto, "VASTSRV", "Unable to Initialize dart db. Please use \"refresh\" button or try again later.", "", false);

*/

  if (!SHDBApiInitialize(TRUE, (char *)"mmdbvast")) 
	PrtMesC(MAILto, "VASTSRV", "SHDBInitialize() failed\n", "", false);
  if (!VastSrvInitialize())
        PrtMesC("", "VASTSRV", "Unable to open VastSrv", "", false);
  if (!MmdbSrvInitialize())
	PrtMesC(MAILto,"VASTSRV","Unable to Initialize MmdbSrv","",false);
  if (!VastSrchInitialize()) 
	PrtMesC(MAILto, "VASTSRV", "Unable to Initialize VastSrch", "", false);

  /* load in the chaindom into memory */

  /* an old VAST Search job */

  if ((indx = WWWFindName(www_info, (char *)"vsid")) >= 0) {
      www_arg = WWWGetValueByIndex(www_info, indx);
      JobID = StringSave(www_arg);

      if ((indx = WWWFindName(www_info, (char *)"pass")) < 0)
        PrtMesC("", "VASTSRV", "Password required -- input password of your Vast Search job", "", false);

      www_arg = WWWGetValueByIndex(www_info, indx);
      Passwd = StringSave(www_arg);
      ret = Check_VastSearch_Password();
      if (ret != 1) {
          if (ret == 2) exit(0);
          PrtMesC("", "VASTSRV", "Incorrect password of a Vast Search result. It could be a wrong Vast Search identifier or a wrong password.", "", false);
      }

  }
  else { 	// new Vast Search job

     if ((indx = WWWFindName(www_info, (char *)"reqid")) >=0) {

	www_arg = WWWGetValueByIndex(www_info, indx);
        ReqId = atoll(www_arg);

 	int iJobID = GrpId2JobId(ReqId);
        if (iJobID) JobID = "VS" + ToString(iJobID);
	else PrtMesC("", "VASTSRV", "Invalid ReqId=", www_arg, false);

     }

  }

  if ((indx = WWWFindName(www_info, (char *)"chaindom")) >= 0) { 
						/* only in URL from VS */
  	www_arg = WWWGetValueByIndex(www_info, indx);
  	if (isInt(www_arg)) aDomId = (long) atol(www_arg);
	else 
	   PrtMesC(MAILto, "VASTSRV", "Invalid Biostruc-feature-set-id: chaindom = ", www_arg, false);
        aMmdbId = aDomId / 10000;
        aChnNo = (aDomId % 10000)/100;

/* non-VS
        if (aMmdbId != 99999) 
		aSdi= DomainInfo2Sdi(aMmdbId, aChnNo, aDomId % 100);
*/

  }
  else if ((indx = WWWFindName(www_info, (char *)"sdid")) >=0) {
	int maxsdi;

/* Could be a sdid for new Vast Search Job currently */

	www_arg = WWWGetValueByIndex(www_info, indx);
	if (isInt(www_arg)) aSdi = (long) atol(www_arg);
	else if (www_arg[0] == '-') aSdi = atol(www_arg);
	else 
	   PrtMesC(MAILto, "VASTSRV", "Non-numeric Structure Domain Identifier (sdid) = ", www_arg, false);

	if (JobID == "") {
	    maxsdi = constructMaxASdi();
	    if (!aSdi || aSdi > maxsdi) 
               PrtMesC(NULL, "VASTSRV", "Vast neighbor data for this domain are not yet available. Please try later again.", "", false);

	    if (aSdi < 0)
	       PrtMesC("", "VASTSRV", "This is an obsolete domain, no Vast information.", "", false);
            aMmdbId = SdiToMmdbId(aSdi);
            aChnNo = SdiToChainNo(aSdi);
	}
        else {
	    aMmdbId=0;
 	    aChnNo = VSSdi2ChnNo(aSdi);
        }
  }
  else PrtMesC("", "VASTSRV", "No chaindom or sdid -- Please input chaindom or Structure Domain Identifier.", "", false);

  if ( (indx= WWWFindName(www_info, (char *)"dispsub.x")) >=0
  		|| (indx= WWWFindName(www_info, (char *)"dispsub")) >=0 )
                iSubBut = DEF_LIST_SUBMIT;
  else if ( (indx=WWWFindName(www_info, (char *)"viewstr.x")) >=0
  		|| (indx=WWWFindName(www_info, (char *)"viewstr")) >=0 )
                iSubBut=VIEW_STR_SUBMIT;
  else if ( (indx=WWWFindName(www_info, (char *)"viewali.x")) >=0
  		|| (indx=WWWFindName(www_info, (char *)"viewali")) >=0 )
                iSubBut=VIEW_ALI_SUBMIT;
  else if ( (indx= WWWFindName(www_info, (char *)"schsub.x")) >=0 
	  	|| (indx= WWWFindName(www_info, (char *)"schsub")) >=0  )
                iSubBut=FIND_SUBMIT;

  if (iSubBut == VIEW_STR_SUBMIT) {
                                /* check whether or not to launch a viewer */
    if ((indx = WWWFindName(www_info, (char *)"calltype")) >= 0) {
      www_arg = WWWGetValueByIndex(www_info, indx);

      switch (www_arg[0]) {
        case 'a':                       /* Cn3d 4.0 */
        case 'c':                       /* Cn3d 3.0 */
                  VastToCn3DAndAli(www_info);
                  exit(0);
        default:
           PrtMesC("", "VASTSRV", "Bad calltype -- check you viewer selection.",
"", false);
      }
    }
    else PrtMesC("", "VASTSrv","No calltype -- please select a viewer.","",false)
;
  }
  else if (iSubBut == VIEW_ALI_SUBMIT) {
    if ((indx = WWWFindName(www_info, (char *)"alitype")) >= 0) {
      www_arg = WWWGetValueByIndex(www_info, indx);

      switch (www_arg[0]) {

        case 'h':                       /* HTML */
        case 'f':                       /* FASTA */
        case 't':                       /* TEXT */
                  VastToCn3DAndAli(www_info);
                  exit(0);
        default:
           PrtMesC("", "VASTSRV", "Bad alitype -- please select a correct alignment display format", "", false);


     }
   }
   else PrtMesC("", "VASTSRV", "No alitype -- please select a alignment display format.", "", false);
 }



  /* subset filtering; identify which subset we're working with */
  if ((indx = WWWFindName(www_info, (char *)"subset")) >= 0) {
    www_arg = WWWGetValueByIndex(www_info, indx);
    if (isInt(www_arg)) subsetnum = (SubSetID)atoi(www_arg);
    else subsetnum = (SubSetID)getSubsetNum(www_arg); 
					/* www_arg is the subsetname now */
  }

  if (iSubBut == FIND_SUBMIT) subsetnum = (SubSetID)1;

  if ((indx = WWWFindName(www_info, (char *)"presubset")) >=0) {
	www_arg = WWWGetValueByIndex(www_info, indx);
	if (isInt(www_arg)) pre_subsetnum = (SubSetID)atoi(www_arg);
	else pre_subsetnum = (SubSetID)getSubsetNum(www_arg);
  }


  if (pre_subsetnum != subsetnum) pagenum = DEFAULT_PAGE;
  else {
     if ((indx = WWWFindName(www_info, (char *)"doclistpage")) < 0)
     	pagenum = DEFAULT_PAGE;
     else {
    	www_arg = WWWGetValueByIndex(www_info, indx);
    	if (isInt(www_arg)) pagenum = (unsigned) atoi(www_arg);
    	else pagenum = DEFAULT_PAGE;
	if (!pagenum) pagenum = DEFAULT_PAGE; 
     }
  }

  if ((indx = WWWFindName(www_info, (char *)"sort")) < 0)
    sortby = DEFAULT_SORT_BY;
  else {
    www_arg = WWWGetValueByIndex(www_info, indx);
    if (isInt(www_arg)) sortby = (SortBy) atoi(www_arg);
    else sortby = DEFAULT_SORT_BY;
  }


  if ((indx = WWWFindName(www_info, (char *)"table")) >= 0) {
	www_arg = WWWGetValueByIndex(www_info, indx);
	cTable = www_arg[0];

        if (cTable == 's' && (indx=WWWFindName(www_info, (char *)"hit")) <0 ) {
	    SaveBioAnnotSetInFile_DBVS(subsetnum, sortby);
 	    return 1;
        }
  }

  iKept = 0;
  Keptfid = NULL;
  if ((iSubBut==DEF_LIST_SUBMIT || iSubBut==FIND_SUBMIT) 
		&& (indx = WWWFindName(www_info, (char *)"hit")) >=0) {
    
    unsigned 	NumLabels;
    char * 	Name; 

    NumLabels = WWWGetNumEntries(www_info);

    for (indx=0; indx < NumLabels; indx++) {
	Name = WWWGetNameByIndex(www_info, indx);
	if (!StrICmp(Name, "hit")) iKept++;
    }

    if (iKept > MAX_KEPT_NBR) {
	sprintf(str, "%d", MAX_KEPT_NBR);
	PrtMesC("", "VASTSRV", "Too many selected neighbors. The maximum allowed is ", str, false);
    }

    Keptfid = new unsigned [iKept];

    iKept = 0;
    for (indx = 0; indx< NumLabels; indx++) {
	Name = WWWGetNameByIndex(www_info, indx);
	if (!StrICmp(Name, "hit")) {

	   char * ptr= NULL;

	   www_arg = WWWGetValueByIndex(www_info, indx);
	   if (isInt(www_arg)) 
		Keptfid[iKept++] = (unsigned)atoi(www_arg);
	   else {
               PrtMesC("", "VASTSRV", "The hit format is not correct.", "", true);
               return 0;
           }
	}
    }

    if (JobID != "")
	numhitsdisplayed = 
	    GetVSVppByNbrs(aSdi, vpp, Keptfid, iKept, sortby);
    else numhitsdisplayed = 
	    constructVastPagesByNbrs(vpp, Keptfid, iKept, aSdi, sortby);

  }

  NbrFlag = DEF_ALL_NBR;
  if (iSubBut==FIND_SUBMIT) {

    unsigned	numhitsbyuid = 0, numhitsbysdid = 0;
    char  	*NonNbrbyuid=NULL, *NonNbrbysdid=NULL;

    if ((indx=WWWFindName(www_info, (char *)"selnbr")) >= 0) {

      www_arg = WWWGetValueByIndex(www_info, indx);


    for (i=0; i< strlen(www_arg); i++)
	if (!isdigit(www_arg[i]) && (!isalpha(www_arg[i])) && (www_arg[i]!=','))
           {
		if (i && www_arg[i-1] == ',') www_arg[i] = ' ';
		else www_arg[i] = ',';
	   }


      if (www_arg[0] != NULLB) { 
	MakeNbrList(www_arg, &SelNames, NULL, &iSel, 1);

        if (iSel) {
	   if (JobID != "") 
		numhitsbyuid =
		   GetVSVppByNbrsWithLikeUnlikeFid(aSdi, vpp+iKept, MaxVPP, 
			SelNames, NULL, 0, iSel, Keptfid, iKept, sortby);
  	   else 
		numhitsbyuid =
	    	   constructVastPagesByNbrsWithLikeUnlikeBsfId(vpp+iKept,
			MaxVPP, aSdi, SelNames, iSel, Keptfid, iKept, sortby);

cerr << "numhits = " << numhitsbyuid << endl;
cerr << "vpp[0].bSdi = " << vpp[0].bSdi << " vpp[i].bAlignId = " << vpp[i].bAlignId << endl;
           numhitsdisplayed += numhitsbyuid;
	   NonNbrbyuid = new char [iSel * 10];
	   CheckNbrs(NonNbrbyuid,SelNames, iSel, vpp, numhitsdisplayed);
	   sprintf(NonNbr, NonNbrbyuid);

	   selnbrstring = new char [strlen(www_arg) +3*iSel];
	   selnbrstring[0] = '\0';
           for (i=0; i< iSel; i++) {
		if (SelNames[i][4] == ' ') SelNames[i][4] = '.';
		StringCat(selnbrstring, SelNames[i]);
		StringCat(selnbrstring, "%2C");
           }
           
	}	/* if (iSel) */
      }		/* if (www_arg[0]) */
    }		/* if (selnbr) */

    if ((indx=WWWFindName(www_info, (char *)"selsdid")) >=0) {

        www_arg = WWWGetValueByIndex(www_info, indx);

/* change return carriage. */

	for (i=0; i< strlen(www_arg); i++)
           if (!isdigit(www_arg[i]) && !isalpha(www_arg[i]) && www_arg[i]!=',')
           {
                if (i && www_arg[i-1] == ',') www_arg[i] = ' ';
                else www_arg[i] = ',';
           }

        if (www_arg[0] != NULLB) {
            MakeNbrList(www_arg, NULL, &SelSds, &iSds, 0);
        }

        if (iSds) {
	    sprintf(www_arg, "%d", SelSds[0]);
	    for (i=1; i< iSds-1; i++) {

		sprintf(str, ",%d", SelSds[i]);
		StringCat(www_arg, str);
	    }
	    if (iSds >1) {
		sprintf(str, ",%d", SelSds[iSds-1]);
	    	StringCat(www_arg, str);
	    }

	    if (JobID != "") 
		numhitsdisplayed=
			GetVSVppByNbrsWithLikeUnlikeFid(aSdi, vpp, MaxVPP, NULL,
				www_arg, 1, iSel, Keptfid, iKept, sortby);

            else 
                numhitsbysdid = 
                   constructVPPBySdi(vpp+numhitsdisplayed, MaxVPP, aSdi,
				www_arg, Keptfid,iKept,sortby);

            numhitsdisplayed += numhitsbysdid;
            NonNbrbysdid = new char [iSds *10];
            CheckNbrsbysdid(NonNbrbysdid, SelSds, iSds, vpp, numhitsdisplayed);

            selsdidstring = new char [StrLen(www_arg)+3*iSds];
            for (i=0; i< iSds; i++) {
                sprintf(str, "%d%%2C", SelSds[i]);
                StringCat(selsdidstring, str);
            }
        }

    }

    if (NonNbrbyuid && NonNbrbyuid[0] != NULLB) {
        sprintf(NonNbr, NonNbrbyuid);
        delete [] NonNbrbyuid;

        if (NonNbrbysdid && NonNbrbysdid[0] != NULLB) {
                StringCat(NonNbr, NonNbrbysdid);
                MemFree(NonNbrbysdid);
        }
    }
    else if (NonNbrbysdid && NonNbrbysdid[0] != NULLB) {
        sprintf(NonNbr, NonNbrbysdid);
        delete [] NonNbrbysdid;
    }

    if (NonNbr[0] != NULLB) NonNbr[StrLen(NonNbr)-2] = NULLB;

    if (!numhitsdisplayed && NonNbr[0] == NULLB)
        PrtMesC("", "VASTSRV", "Please select a neighbor by using checkbox or input a name or a number (sdi) in the test box following \"Find\", then click \"Find\" again.", "", false);

    NbrFlag = SELECTED_NBR;
    vpp[0].numHits = numhitsdisplayed;

  }     /* if (Find) */


  if (NbrFlag == DEF_ALL_NBR) {
        if (subsetnum <6) {
	    if (JobID != "")  {
	 	numhitsdisplayed=
		   GetVSVpp(aSdi, vpp+iKept, NUM_HITS_PER_PAGE, Keptfid, iKept,
			subsetnum, sortby, NUM_HITS_PER_PAGE, pagenum);
	    }
	    else 
    	    	numhitsdisplayed=
		    constructVastPages(vpp+iKept,NUM_HITS_PER_PAGE,
			aSdi, (unsigned *)Keptfid, iKept, subsetnum, sortby,
			NUM_HITS_PER_PAGE, pagenum);

	    vpp[0].numHits = vpp[iKept].numHits;
    	    numhitsdisplayed += iKept;

	}
	else vpp[0].numHits = iKept;
  }

  if (numhitsdisplayed <= 0) {
      if (subsetnum == 6) 
	PrtMesC("", "VASTSRV", "Please use checkbox(es) to select neighbor(s).", "", false); 
      if (NbrFlag == DEF_ALL_NBR) {
     	printf("Content-type: text/html\n\n");
     	printf("<body bgcolor = \"#f0f0f0\">\n");
     	printf("<br>\n");
     	printf("<h2><a href=\"%s%s#VASTNonbr\">", URLBase, HELPname);
	printf("VAST did not find any structure neighbors.</a></h2>");
      }
      else {
	printf("Content-type: text/html\n\n");
        printf("<body bgcolor = \"#f0f0f0\">\n");
	printf("<br><h2>VASTSRV:<p>\n");
	if (NonNbr && NonNbr[0] != NULLB) {
	   if (strchr(NonNbr, ',')) 
	      sprintf(str," are not structure neighbors. ");
	   else sprintf(str," is not a structure neighbor. ");
 	   printf("<font color=#CC6600>%s</font>%s\n", NonNbr, str);
	}
      }
 
      exit(0);
  }

  if (iSubBut == DEF_LIST_SUBMIT || iSubBut == FIND_SUBMIT) { 

    totalnumhits = vpp[0].numHits;
    if (NbrFlag == DEF_ALL_NBR && totalnumhits <=0) 
	PrtMesC(MAILto, "VASTSRV", "TotalNumberOfHits == 0, yet TotalNumberOfHitsDisplayed > 0.", "", false);

    if (JobID=="" || (JobID!="" && strlen(JobID.c_str()) <6)) {  
							/* new Vast Search */
    	aSdi = vpp[0].aSdi;
    	if (aSdi <= 0 ) 
	    PrtMesC(MAILto, "VASTSRV", "TotalNumberOfHits > 0, yet aSdi == 0.", "", false);
    }

    if ((indx = WWWFindName(www_info, (char *)"cmd")) >=0) {
	Cmd = WWWGetValueByIndex(www_info, indx);
	if (!StrICmp(Cmd, "graph")) {
	    if ((indx = WWWFindName(www_info, (char *)"imgsize")) <0) 
		PrtMesC("", "VASTSRV", "Missing imgsize -- please use \"imgsize = number\".", "", false);
	    else {
		www_arg = WWWGetValueByIndex(www_info, indx);
		if (isInt(www_arg)) ImgSize = atoi(www_arg);
		else {
		    PrtMesC(MAILto, "VASTSRV", "Missing imgsize: imgsize = ", www_arg, false);
		}  /* sth is wrong */
		MakeVastTableGraph(vpp, numhitsdisplayed, aDomId, 
			iKept, selnbrstring, selsdidstring, NonNbr, 
			sortby, subsetnum, pagenum, totalnumhits, cTable, 
			NbrFlag, ImgSize,1);
	   }
	}
    }
    else if (cTable == 's') {
	BiostrucAnnotSetPtr     pbsa = NULL;
	AsnIoPtr                aip;
	
	if (JobID == "")
		pbsa = constructBASPFromVastPagePtr(vpp, numhitsdisplayed, 0);
	else pbsa = constructBASPFromVastPagePtr(vpp, numhitsdisplayed, 1);
	printf("Content-type: application/octet-stream\r\n\r\n");
 	        aip = AsnIoNew(ASNIO_TEXT_OUT, stdout, NULL, NULL, NULL);
        BiostrucAnnotSetAsnWrite(pbsa, aip, NULL);
        AsnIoFlush(aip);
        AsnIoClose(aip);
    }
    else   
    	MakeVastTableGraph(vpp,numhitsdisplayed,aDomId, iKept, 
		selnbrstring, selsdidstring, NonNbr, sortby, subsetnum, pagenum,
		totalnumhits, cTable, NbrFlag, 0, 0);
  } 


  MMDBFini();
  VastSrchFinish();
 
  if (selnbrstring != NULL) {
 	delete [] selnbrstring;
  	for (i=0; i< iSel; i++) delete [] SelNames[i];
  	delete [] SelNames;
  }

  if (selsdidstring != NULL) {
       delete [] selsdidstring;
       delete [] SelSds;
  }

  // Dart_FiniReal(); 
  
  return 0;

} /* end Main */

