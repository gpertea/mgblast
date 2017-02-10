/* 
 * $Id: vastgraph_comb.cpp,v 1.1 2005/07/26 17:12:50 chenj Exp $
 *
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
 *
 * Author: Jie Chen
 *
 * $Log: vastgraph_comb.cpp,v $
 * Revision 1.1  2005/07/26 17:12:50  chenj
 * Making linux VSNbr.cgi
 *
 * Revision 1.3  2003/01/14 20:49:20  chenj
 * Imported sources
 *
 *
 * This is to show the "red-cloud" graph. 
 * The footprint alignment is broken at the position of insert.
 *
 * ===========================================================================
*/

#include "hUtilib.hpp"
#include "VastSrchUti.hpp"
#include "vastuti.hpp"
#include "SHGlobal.hpp"

#include <sstream>

#include <time.h>
#include <cddapi.h>
#include "dartutil.h"


#define iNcolors 13

static double pix_per_res;
static char *	QueryBit;
static gdImagePtr im;
static short    white, black, red, blue, gray;
static unsigned iDartCol[13];

extern 	Char 	URLcgi[MAX_TBUFF],MMDBCGIname[MAX_TBUFF], CGIname[MAX_TBUFF];
extern 	string 		MAILto, JobID, Passwd;
extern 	char 		ENTREZurl[PATH_MAX], Database[PATH_MAX];
extern 	Char 		VSPATH[PATH_MAX], CDDurl[PATH_MAX];
extern 	unsigned	aSdi, aMmdbId, aChnNo;
extern	long long	ReqId;


using namespace ncbi;
using namespace SHProjNS;

static void
GetLeftAnglePoints(gdPoint *points, unsigned x1, unsigned x2, unsigned dx, unsigned y, unsigned dy)
{
        unsigned    midp;

        midp = (x1+x2)/2;

        points[0].x = MIN(midp, x1+dx); points[0].y = y;
        points[1].x = x1;               points[1].y = y+0.5*dy;
        points[2].x = MIN(midp, x1+dx); points[2].y = y+dy;

}       /* GetLeftAnglePoints */




static void
GetRightAnglePoints(gdPoint *points, unsigned x1, unsigned x2, unsigned dx, unsigned y, unsigned dy)
{
        unsigned midp;

        midp = (x2+x1)/2;

        points[0].x = MAX(midp, x2-dx); points[0].y = y;
        points[1].x = x2;               points[1].y = y+0.5*dy;
        points[2].x = MAX(midp, x2-dx); points[2].y = y+dy;

}       /* GetRightAnglePoints */




static void GetPoints(gdPoint *points, short *n, char * flag, 
			unsigned x1, unsigned x2, unsigned y, unsigned dy)
{
	unsigned midp, dx;

	midp = (x1+x2)/2;

	switch (flag[0]) {
	case 'L':  	/* Left zigzag */
                GetRightAnglePoints(points, x1, x2, 7, y, dy);
                points[3].x = x1;               points[3].y = y+dy;
                points[4].x = MIN(x1+3, midp);  points[4].y = y+dy*0.75;
                points[5].x = x1;               points[5].y = y+dy*0.5;
                points[6].x = MIN(x1+3, midp);  points[6].y = y+dy*0.25;
                points[7].x = x1;               points[7].y = y;
                *n = 8;
		break;
	case 'R':	/* right zigzag */
                GetLeftAnglePoints(points, x1, x2, 7, y, dy);
                points[3].x = x2;               points[3].y = y+dy;
                points[4].x = MAX(x2-3, midp);  points[4].y = y+dy*0.75;
                points[5].x = x2;               points[5].y = y+dy*0.5;
                points[6].x = MAX(x2-3, midp);  points[6].y = y+dy*0.25;
                points[7].x = x2;               points[7].y = y;
                *n = 8;
		break;
	case 'W':	/* two zigzag lines */
                points[0].x = x1;               points[0].y = y;
                points[1].x = x2;               points[1].y = y;
                points[2].x = MAX(x2-3, midp);  points[2].y = y+dy*0.25;
                points[3].x = x2;               points[3].y = y+dy*0.5;
                points[4].x = MAX(x2-3, midp);  points[4].y = y+dy*0.75;
                points[5].x = x2;               points[5].y = y+dy;
                points[6].x = x1;               points[6].y = y+dy;
                points[7].x = MIN(x1+3, midp);  points[7].y = y+dy*0.75;
                points[8].x = x1;               points[8].y = y+dy*0.5;
                points[9].x = MIN(x1+3, midp);  points[9].y = y+dy*0.25;
                *n = 10;
                break;
        case 'A':	/* Aligned, so two angles */
		if (x2-x1 <=8) dx=3;
		else if (x2-x1 <= 18) dx = 4;
		else dx = 7;

                points[0].x = MIN(midp, x1+dx);  points[0].y = y;
                if (x2-x1 < 8) {
			points[1].x = x2;  	/* points[1].y = y; */
			points[3].x = x2;  	/* points[3].y = y + dy; */
		}
		else {
			points[1].x = MAX(midp, x2-dx); 
			points[3].x = MAX(midp, x2-dx);  
		}
                points[2].x = x2;                points[2].y = y+0.5*dy;
                points[4].x = MIN(midp, x1+dx);  points[4].y = y+dy;
                points[5].x = x1;                points[5].y = y+0.5*dy;

		points[1].y = y; 		 points[3].y = y+dy;

                *n = 6;
                break;
        case 'U':	/* Un...,  special for query */
                points[0].x = x1;       points[0].y = y;
                points[1].x = x2;       points[1].y = y;
                points[2].x = x2;       points[2].y = y+dy;
                points[3].x = x1;       points[3].y = y+dy;
                *n = 4;
                break;
	default:
		PrtMesC(MAILto, "VASTSRV (VastGraph)", 
		    "Error in GetPoints -- didn't choose points", "", false);
	}
		
}	/* GetPoints */



#define HTMLurl "href=\"%s%s?sdid=%d&viewali=View&action=0&alitype=h"
#define VSHTMLurl "href=\"%s%s?chaindom=%d&viewali=View&action=0&alitype=h"

static unsigned 
BlockMapOrImg(bool ismap, unsigned x, unsigned y, unsigned x0, unsigned from, 
	unsigned to, unsigned tick_from, unsigned tick_to, short color, 
	unsigned dy, unsigned pre_t, unsigned next_f, bool isaligned, 
	unsigned Fsid, VastPageDataPtr vpp, 
	int iSlv, unsigned numhitsdisplayed, FILE *File)
{
  char		altstr[200], str[MAX_TBUFF], str2[20], strnum[MAX_TBUFF];
  short		i;
  unsigned 	x1=0, x2=0;
  gdPoint  	points[10];
  
  CalCoor(&x1, &x2, x, from, to, pix_per_res, (unsigned) MaxSeqImgSize);
  if (iSlv < 0) if (x1-x0 <2) x1 = x0+2;
  strnum[0] = '\0';

  if (ismap == TRUE) {

    if (iSlv < 0) {
    	if (dy==8) sprintf(altstr, "Unaligned residues ");
    	else sprintf(altstr, "Maximal aligned region: from residues ");
	sprintf(str, "%d to %d, ", tick_from, tick_to);
      	StringCat(altstr, str);
	StringCat(altstr,"click for multiple alignment with neighbors.");
	for (i=0; i< numhitsdisplayed; i++) {
            sprintf(str2, "%d%%2C", getfid(vpp[i].bSdi, vpp[i].bAlignId));
            StringCat(strnum, str2);
	}	
    }
    else {
	sprintf(altstr, 
		"Aligned residues %d to %d, click for sequence alignment.",
		tick_from, tick_to);
    }

    fprintf(File, "<area shape=rect coords=%d,%d,%d,%d ", 
		x1, y, x2, y+dy);

    fprintf(File, HTMLurl, URLcgi, CGIname, aSdi);
    if (iSlv < 0)  {
	fprintf(File, "&nbr_complexity=0&allbfid=%s", strnum);
    }
    else fprintf(File, "&nbr_complexity=1&hit=%d", 
			getfid(vpp[iSlv].bSdi, vpp[iSlv].bAlignId));
    if (JobID != "") {
	if (ReqId) fprintf(File, "&reqid=%llu", ReqId);
	else 
	      fprintf(File, "&vsid=%s&pass=%s", JobID.c_str(), Passwd.c_str());
    }

    fprintf(File, "\"");
    fprintf(File, " alt=\"%s\" title=\"%s\">\n", altstr, altstr);
  }
  else {
     if (isaligned == TRUE) {
  	if (!pre_t)  {
	    short   n;

	    if (!next_f) {
		GetLeftAnglePoints(points, x1, x2, 7, y, dy);
		points[3].x = x2; 	points[3].y = y+dy;
		points[4].x = x2; 	points[4].y = y;
		n = 5;
	    }
	    else GetPoints(points, &n, (char *)"A", x1, x2, y, dy);

	    gdImageFilledPolygon(im, points, n, color);
	
	}
	else if (next_f) {
		GetRightAnglePoints(points, x1, x2, 7, y, dy);
		points[3].x = x1; 	points[3].y = y+dy;
		points[4].x = x1; 	points[4].y = y;
		gdImageFilledPolygon(im, points, 5, color);
	}
	else gdImageFilledRectangle(im, x1, y, x2, y+dy, color);
     }
     else gdImageFilledRectangle(im, x1, y, x2, y+dy, color);
  }
		
  return(x2);

}  /* end of BlockImg */


#define EntrezLink "href=\"%s?cmd=Retrieve&db=%s&list_uids=%d&dopt=%s\""

static void NameMapOrImg(bool ismap, unsigned *x, unsigned y, char * pdbname, 
		Char chainname, unsigned iDomain, unsigned iGi, unsigned seqlen,
		char * pcSlaveName, FILE *File, bool hasLine)
{  

  Char 	cTmp[10], str[10], altstr[MAX_TBUFF];
  unsigned	x0;

  sprintf(cTmp, "%s %c", pdbname, chainname);
  if (iDomain) {
	sprintf(str, " %d", iDomain);
	StringCat(cTmp, str);
  }
  if (ismap == TRUE) {
    sprintf(altstr, "%s", pdbname);
    if (chainname != ' ') {
	sprintf(str, "_%c", chainname);
	StringCat(altstr, str);
    }
    if (iDomain) {
		StringCat(altstr, "_");
		sprintf(str, "%d", iDomain);
		StringCat(altstr, str);
    }
    fprintf(File, "<area shape=rect coords=%d,%d,%d,%d ", 
	    *x, y, *x+StrLen(cTmp)*FontBW, y+FontBH);
    if (pcSlaveName != NULL) {
	char *chp = NULL, strtmp[MAX_TBUFF];

	sprintf(str, ": %d", seqlen);
    	StringCat(altstr, str),
    	StringCat(altstr, " residues");
        if (!iDomain) StringCat(altstr, ", ");
	else {
	   if (chainname != ' ') sprintf(str, " on chain %c, ", chainname);
	   else sprintf(str, " on chain, ");
	   StringCat(altstr, str);
	}
        while (chp = strchr(pcSlaveName, '\''))
	{
	    StrCut(strtmp, pcSlaveName, 1, chp-pcSlaveName);
	    StringCat(altstr, strtmp);
	    StringCat(altstr, "\\'");
	    pcSlaveName = chp+1;
	}
	StringCat(altstr, pcSlaveName);
        fprintf(File, "href=\"%s%s?uid=%s", URLcgi, MMDBCGIname, pdbname);
        fprintf(File, "&Dopt=s\" ");
    }
    else if (iGi) {
	StringCat(altstr, ", click for Entrez sequence summary");
	fprintf(File, EntrezLink, ENTREZurl, "protein", iGi, "GenPept"); 
    }

    fprintf(File, "alt=\"%s\" title=\"%s\">\n", altstr, altstr);
  }
  else {
    short color;
	
    gdImageString(im, gdFont7X13b, *x, y, cTmp, blue); 
    y += FontBH;
    x0 = *x + StrLen(pdbname)*FontBW;
    if (hasLine) gdImageLine(im, *x, y, x0, y, blue);
    x0 += FontBW;
    if (chainname != ' ' && hasLine) gdImageLine(im, x0, y, x0+FontBW, y, blue);
    x0 += 2*FontBW;
    if (iDomain && hasLine) gdImageLine(im, x0, y, x0+FontBW, y, blue);
  }
 
  *x += 11*FontBW;

}  /* end NameMapOrImg */


static void PrintAliMap(unsigned x1, unsigned y1, unsigned x2, unsigned y2, unsigned tick_f, 
	unsigned tick_t, unsigned Fsid, VastPageDataPtr vpp, unsigned iSlv, FILE *File)
{
   char cTmp[MAX_TBUFF];

   sprintf(cTmp, "Aligned residues %d to %d, click for sequence alignment", 
		tick_f, tick_t);
   fprintf(File, "<area shape=rect coords=%d,%d,%d,%d ", 
		x1, y1, x2, y2);
   fprintf(File, HTMLurl, URLcgi, CGIname, aSdi);
   fprintf(File, "&nbr_complexity=1&hit=%d", 
				getfid(vpp[iSlv].bSdi, vpp[iSlv].bAlignId));
   if (JobID != "") {
       	if (ReqId) fprintf(File, "&reqid=%llu", ReqId);
     	else 	 
	     fprintf(File, "&vsid=%s&pass=%s", JobID.c_str(), Passwd.c_str());
   }

   fprintf(File, "\"");
   fprintf(File, " alt=\"%s\" title=\"%s\">\n", cTmp, cTmp);
}



static unsigned
InsertOnSlaveMapOrImg(bool ismap, unsigned x, unsigned y, int x0, unsigned from, 
	unsigned to, unsigned tick_from, unsigned tick_to, unsigned pre_to, unsigned next_from, 
	IntvlPairDataPtr ipdp, unsigned *segNo,unsigned segNum, unsigned Fsid, 
	VastPageDataPtr vpp, unsigned iSlv, FILE *File)
{
  short n=0; 
  unsigned x1, x2, dy,i, idx, next, f_res, t_res;
  gdPoint points[10];

  dy = 10;
  idx = segNo[0];
  t_res = ipdp[idx-1].mstto;
  CalCoor(&x1, &x2, x, from, t_res, pix_per_res, MaxSeqImgSize);
  if (x1-x0 < 2) x1 = x0+2;
  x2 -= 1.0 + MIN(0, (unsigned)(pix_per_res/2.0));
  if (ismap == TRUE) 
	PrintAliMap(x1, y, x2, y+dy, tick_from, ipdp[idx-1].slvto, 
		Fsid, vpp, iSlv, File);
  else {
     if (!pre_to) {

	GetLeftAnglePoints(points, x1, x2, 7, y, dy);
	points[3].x = x2;       points[3].y = y+dy;
        points[4].x = x2;       points[4].y = y;
        n = 5;

	gdImageFilledPolygon(im, points, n, red);
     }
     else 
     	gdImageFilledRectangle(im, x1, y, x2, y+dy, red);
  }
	
  for (i=0; i< segNum; i++) {
    idx = segNo[i];
    f_res = ipdp[idx].mstfrom;
    CalCoor(&x1, NULL, x, f_res, 0, pix_per_res, MaxSeqImgSize);
    x1 += 1.0 + MIN(0, (unsigned)(pix_per_res/2.0));
    if (i< segNum-1) {
      next = segNo[i+1];
      t_res = ipdp[next].mstfrom-1;
      CalCoor(NULL, &x2, x, 0, t_res, pix_per_res, MaxSeqImgSize);
      x2 -= 1.0 + MIN(0, (unsigned)(pix_per_res/2.0));
      if (ismap==TRUE) 
	PrintAliMap(x1, y, x2, y+dy, ipdp[idx].slvfrom, ipdp[next].slvfrom-1, 
		Fsid, vpp, iSlv, File);
      else gdImageFilledRectangle(im, x1, y, x2, y+dy, red);
    }
  }

  CalCoor(NULL, &x2, x, 0, to, pix_per_res, MaxSeqImgSize);
  if (ismap == TRUE) { 
       PrintAliMap(x1, y, x2, y+dy, ipdp[idx].slvfrom, tick_to, Fsid,
		vpp, iSlv, File);
  }
  else {
  	if (next_from) {
  	     GetRightAnglePoints(points, x1, x2, 7, y, dy);
  	     points[3].x = points[4].x = x1;
  	     points[3].y = y+dy; points[4].y = y;
	     n = 5;
	     gdImageFilledPolygon(im, points, 5, red);

        }
  	else gdImageFilledRectangle(im, x1, y, x2, y+dy, red);
  }

  return(x2);

}	/* end InsertOnSlaveMap... */




static void
AlignmentMapOrImg(bool ismap, unsigned x, unsigned y, unsigned Fsid, unsigned iSlv, 
	unsigned seqlen_s, VastPageDataPtr vpp, char * stats, short sortnamelen, 
	FILE *File)
{
  char		pdbId_s[6], chnLett_s;
  char *	pcSlaveName = NULL;
  short		i, j;
  int 		x0 = -1, ali_to_s=-1;
  unsigned	dy, from_m,from_s, to_m, to_s, ali_from_m, ali_to_m=0;
  unsigned	ali_from_s=0, pre_to_s=0, segNo[20], segNum = 0, domNo_s;
  IntvlPairData	*ipdp;

  dy = 10;

  if (ismap == TRUE) {
      pcSlaveName = new char [2*DBStrSize+4];
      pcSlaveName[0] = '\0';
      if(vpp[iSlv].bDescr1[0] != '\0' || vpp[iSlv].bDescr2[0] != '\0') {
          if(vpp[iSlv].bDescr1[0] != '\0')
                strcpy(pcSlaveName,vpp[iSlv].bDescr1);
          if(vpp[iSlv].bDescr2[0] != '\0')
                strcat(pcSlaveName,vpp[iSlv].bDescr2);
      }
  }

  StrCut(pdbId_s, vpp[iSlv].bDomName, 1, 4);
  chnLett_s = (vpp[iSlv].bDomName)[4];
  domNo_s = (unsigned)atoi(vpp[iSlv].bDomName+5);
  ipdp = vpp[iSlv].Ipdp;

  NameMapOrImg(ismap, &x, y, pdbId_s, chnLett_s, domNo_s, 0, seqlen_s, 
		pcSlaveName, File, 1);

  if (true == ismap) delete [] pcSlaveName;

  for (i = 0; i< vpp[iSlv].IpdpLen; i++) {
    from_m = ipdp[i].mstfrom;
    to_m = ipdp[i].mstto;
    from_s = ipdp[i].slvfrom;
    to_s = ipdp[i].slvto;

/*  red piece for aligned region */		   
    if (!i) {
	ali_from_m = from_m;
	ali_from_s = from_s;
    }

    
    if (i && from_m == ali_to_m+1) {
		if (from_s > ali_to_s+1) segNo[segNum++] = i;
    }
    else {
  	if (i) {
           if (!segNum) 
	 	x0 = BlockMapOrImg(ismap, x, y, x0, ali_from_m, ali_to_m, 
			ali_from_s, ali_to_s, red, dy, pre_to_s, 0, TRUE, Fsid,
			vpp, iSlv, 0, File); 
	   else
		x0 = InsertOnSlaveMapOrImg(ismap, x, y, x0, ali_from_m, 
			ali_to_m, ali_from_s, ali_to_s, pre_to_s, 0, ipdp, 
			segNo, segNum, Fsid, vpp, iSlv, File);

	   for (j = ali_from_m-1; j< ali_to_m; j++) QueryBit[j] = '1';

      	   if (!pre_to_s) pre_to_s = ali_to_s;
      	}
	

      ali_from_m = from_m;
      ali_from_s = from_s;
      segNum = 0;
    }

    ali_to_m = to_m;
    ali_to_s = to_s;

  }	/* for (i = 0 ... */	


/* rectangle with red color from from_m to to_m ----  aligned block */

  if (!segNum)  
  	x0 =BlockMapOrImg(ismap, x, y, x0, ali_from_m, ali_to_m, ali_from_s,
		ali_to_s, red, dy, pre_to_s, 1, TRUE, Fsid, vpp, iSlv, 0, File);
  else 
	x0 = InsertOnSlaveMapOrImg(ismap, x, y, x0, ali_from_m, ali_to_m, 
		ali_from_s, ali_to_s, pre_to_s, 1, ipdp, segNo, segNum, 
		Fsid, vpp, iSlv, File);

  for (j = ali_from_m-1; j< ali_to_m; j++) QueryBit[j] = '1';

/* print stats */

  x =GraphWidth - (MAX(sortnamelen, StrLen(stats)) + StrLen(stats))*FontBW/2.0;
  if (ismap==FALSE) {
	gdImageString(im, gdFont7X13b, x, y, stats, blue);
  }

} /* end of AlignImg() */




static void 
GetLeft4Points(unsigned x1, unsigned x2, unsigned y1, unsigned y2, gdPointPtr points)
{
	unsigned i=0;

	points[i].x = x2; 	points[i++].y = y1;
	points[i].x = x1; 	points[i++].y = y1+2;
	points[i].x = x1; 	points[i++].y = y2-2;
	points[i].x = x2; 	points[i++].y = y2;

}




static void 
GetRight4Points(unsigned x1, unsigned x2, unsigned y1, unsigned y2, gdPointPtr points) 
{
	unsigned i=0;

	points[i].x = x1; 	points[i++].y = y2;
	points[i].x = x2; 	points[i++].y = y2-2;
	points[i].x = x2; 	points[i++].y = y1+2;
	points[i].x = x1; 	points[i++].y = y1;
}


#define CDDsrv "href=\"%s?ascbin=2&maxaln=10&seltype=3&uid=%d&aln=%s&querygi=%d&querynm=%s&version=%s\""
#define CDDOpt "db=cdd&term="
#define DomOpt "db=structure&cmd=Display&dopt=structure_domains&from_uid="


static void 
QueryMapOrImg(bool ismap, unsigned x, unsigned y, char * pcPDB, 
	unsigned Fsid,Char cChain, unsigned domNo, unsigned gi, unsigned seqlen,
	DomIntvlDataPtr dip, unsigned domintvl_no, VastPageDataPtr vpp,
	unsigned numhitsdisplayed, FILE *File)
{
  string	cTmp, shortname, defline, def2;
  string	cdalign, querynm, definition;
  char 		**CddName;
  DenseSegPtr	*dsp;
  int 		x0=-1;
  short		i, color, char_num;
  unsigned 	alinumseg,j, y0, x1, x2, y1, y2, len, dont, x_ori, numseg;
  int 		*starts, *lens;
  unsigned	index, ntick, tick, right, dy=0, from, to, len1, len2, minali; 
  OverLoc	*head, *end;
  SeqAnnotPtr	sap = NULL;
  SeqAlignPtr	salp= NULL;
  static unsigned 	ticksteps[9] = {5, 10, 20, 25, 50, 100, 200, 150, 500};

  x_ori = x;

  NameMapOrImg(ismap, &x, y, pcPDB, cChain, 0, gi, seqlen, NULL, File, 0);

  y1 = y+3;

/* sequence ruler */

  j = minali = 0;
  while (j < seqlen) {
     char  	bitj;
     unsigned 	k, l, pre_t, next_f;
     bool	aliflag;

     bitj = QueryBit[j];
     switch (bitj) {
	case '0': color=gray; aliflag = FALSE; dy = 8;
		  pre_t = next_f = 1; y0 = y1;
		  break;
	case '1': color = red; aliflag = TRUE; dy = 14;
		  if (!minali) {
			minali = 1; pre_t = 0;
		  }
		  else pre_t = 1;
		  y0 = y;
		  break;
     }

     k = j+1;
     while (k< seqlen && bitj == QueryBit[k]) k++;

     if (bitj == '1') {
	l = k+2;
     	while (l< seqlen && QueryBit[l] != bitj) l++;
     	if (l >= seqlen) next_f = 1;
 	else next_f = 0;
     }	
     
     x0 = BlockMapOrImg(ismap, x, y0, x0, ++j, k, j, k, color, dy, pre_t, 
		next_f, aliflag, Fsid, vpp, -1, numhitsdisplayed, File); 
     j = k;
  }

  if (ismap == FALSE) {

    /* Add ticks */

    CalCoor(NULL, &right, x, 0, seqlen, pix_per_res, MaxSeqImgSize);

    y2 = y-5;
    if (QueryBit[0] == '1') gdImageLine(im, x, y, x, y2, red);
    else  gdImageLine(im, x, y1, x, y2, red);
    if (QueryBit[seqlen-1] == '1') gdImageLine(im, right, y, right, y2, red);
    else gdImageLine(im, right, y1, right, y2, red);

    y2 -= FontMH;

    gdImageString(im, gdFont6X12, x, y2, (char *)"1", red);
    cTmp = ToString(seqlen);
    gdImageString(im, gdFont6X12, right-(cTmp.size()-1)*FontMW, y2, (char *)cTmp.c_str(), red);

    y2 += FontMH;

    for (i=0; i< 9; i++) {
      ntick = seqlen/ticksteps[i];
      if (ntick < 10) {
	tick = ticksteps[i];
	break;
      }
    }

    x0 = x + 2*FontMW;
    dont = right - cTmp.size()*FontMW;
    for (i=1; i<=ntick; i++){
      x1 = (int)(x+(float)(i*tick-1) * pix_per_res);
      if (x1 > dont) continue;

      if (QueryBit[i*tick-1] == '1') gdImageLine(im, x1, y, x1, y2, red);
      else gdImageLine(im, x1, y1, x1, y2, red);

      cTmp = ToString(i*tick);
      len = cTmp.size()*FontMW/2;
      x2 = x1 + len;
      x1 -= len;
      if (x2<dont && (x2-x0 > len)) { 
	gdImageString(im, gdFont6X12, x1, y2-FontMH, (char *)cTmp.c_str(), red);
	x0 = x2 + FontMW;
      }
    } 
  }

/* "3d Domains" & "CDs" */

  if (domintvl_no >= 1) {
    x = x_ori;  y += FontBH+5;

    if (ismap == FALSE) 
	NameMapOrImg(ismap, &x, y, (char *)"3d Dom.", ' ', 0, 0, 0,NULL,NULL,0);
    else {
	fprintf(File, "<area shape=rect coords=%d,%d,%d,%d ",
		x, y, x+StrLen("3d Dom.")*FontBW, y+FontBH);
      	if (JobID == "") {
	    cTmp = "Click for Entrez 3d Domains";
	    fprintf(File, "href=\"%s?%s%d\" ", ENTREZurl, DomOpt, aMmdbId);
        }
	else cTmp = "3d Dom";
	fprintf(File, "alt=\"%s\" title=\"%s\">\n", cTmp.c_str(), cTmp.c_str());
	
	x += 11*FontBW;
    }
	
    if (!domNo) domintvl_no = 1;
    for (i = 0; i< domintvl_no; i++) {

      if (domNo && !dip[i].domNo) continue; 

      y2 = y+FontBH+2;
      if (dip[i].domNo != domNo) {
	if (ismap == FALSE) color = gray;
	y1 = y+1;
	y2 -= 2;
      }
      else {
	if (ismap == FALSE) color = red;
	y1 = y;
      }

      CalCoor(&x1, &x2, x, dip[i].frm, dip[i].tu, pix_per_res, MaxSeqImgSize);
      if (ismap == TRUE) {
	
	  cTmp="Residue " + ToString(dip[i].frm) + " to " + ToString(dip[i].tu);
	  cTmp += ", click for structure neighbors";

	  fprintf(File, "<area shape=rect coords=%d,%d,%d,%d ",
		  x1, y1, x2, y2);
		
	  int iJobID = atoi(JobID.substr(2).c_str());
	  fprintf(File, "href=\"%s%s?sdid=%d", URLcgi, CGIname, 
					Dom2Sdi(iJobID, cChain, dip[i].domNo));
	  if (JobID != "") fprintf(File,"&reqid=%llu", ReqId);
	  fprintf(File, "&Dopt=s\" alt=\"%s\" title=\"%s\">\n", 
						cTmp.c_str(), cTmp.c_str());
      }
      else {
	  gdPoint points[10];
	  short	  n;	
	  
	  if (dip[i].domNo != domNo) {
	     if (dip[i].frm == 1) {
		points[0].x = x1;   	points[0].y = y1;
		points[1].x = x1;	points[1].y = y2;
		n = 2;
	     }
	     else {
		GetLeft4Points(x1, MIN((x1+x2)/2, x1+2), y1, y2,points);
		n = 4;
	     }

	     if (dip[i].tu == seqlen) {
		points[n].x = x2;	points[n++].y = y2;
		points[n].x = x2;	points[n++].y = y1;
	     }
	     else {
		GetRight4Points(MAX((x1+x2)/2, x2-2), x2, y1, y2, points+n);
		n += 4;
	     }
	     gdImageFilledPolygon(im, points, n, color);
	  }
	  else 
	  	gdImageFilledRectangle(im, x1, y1, x2, y2-2, color); 

	  if (dip[i].domNo) cTmp = ToString(dip[i].domNo);
	  else cTmp = (string)"Chain " + cChain;
	  if (x2-x1> cTmp.size()*FontBW) {
	    x0 = (x1+x2-cTmp.size()*FontBW)/2.0; 
	    gdImageString(im, gdFont7X13b, x0, y+1, (char *)cTmp.c_str(),white);
	  }
      }

    } /* for (i=0; i< domintvl_no  */
  }

  if (!gi) return;

  SeqEntryPtr sep = (SeqEntryPtr) constructSeqEntryForGi(gi, true, Database);
  if (!sep) return;
  BioseqPtr bseqp = (BioseqPtr)sep->data.ptrvalue;

  if (!bseqp) {
        PrtMesC(MAILto, "MMDBSRV", "SeqEntryPtr not NULL, but BioseqPtr is NULL and gi=",  (char *)ToString(gi).c_str(), false);
  }

  sap = CddSynchronousQuery(bseqp, 0.01, true, true, false, (char *)"", true);

  if (sap) {

    unsigned    	*PssmId;
    unsigned short	*iColor, CdNum, thisColor=0;
    short 		 *iClus;
    right = (Int4)(x + (float)(seqlen-1)*pix_per_res); // why give a new value?

    end = head = NewOverLoc(seqlen);
    head->y = y0 = y+7+FontBH;
    head->next = NULL;

    for (salp = (SeqAlignPtr)sap->data, CdNum=0; salp!=NULL; 
		salp=salp->next, CdNum++); 
    
    iColor = new unsigned short [CdNum];
    PssmId = new unsigned [CdNum];
    dsp = new DenseSegPtr [CdNum];
    CddName = new char* [CdNum];
    for (j=0; j< CdNum; j++) {
	CddName[j] = new char [30];
	CddName[j][0] = '\0';
    }
    iClus = new short [CdNum];
    
    string strtmp;
    for (salp = (SeqAlignPtr)sap->data, i=0; salp!=NULL; salp=salp->next, i++) {
      dsp[i] = (DenseSegPtr)salp->segs;
      PssmId[i] = GetPSSMID(dsp[i]); 
      iClus[i] = -1;

      if (PssmId[i]) {
      	Dart_CDGi2Acc_external(PssmId[i], strtmp);
	sprintf(CddName[i], strtmp.c_str());
      }
      else PrtMesC(MAILto, "VastSrv", "PssmId=0", "", false);
    }

    for (i=0; i< CdNum; i++) {

	if (iClus[i] >= 0) continue;

	if (ismap == false) iColor[i] = (thisColor++) % iNcolors;
	iClus[i] = i;

	for (j=i+1; j< CdNum; j++) {
	      if (PssmId[i] == PssmId[j]) {
		    if (ismap == false) iColor[j] = iColor[i];
		    iClus[j] = i;
	      }
	}
    }


    querynm = (string)pcPDB + cChain;
    for (i=0; i< querynm.size(); i++)
             if (isalpha(querynm[i])) querynm[i] = tolower(querynm[i]);
    if (querynm[4] == ' ') querynm[4] = '+';
    querynm + "((query)";

    for (i=0; i< CdNum; i++) {

      numseg = (dsp[i])->numseg;
      starts = (dsp[i])->starts;
      lens = (dsp[i])->lens;
      from  = (unsigned)(starts[0] + 1);
      to = (unsigned)(starts[(numseg-1)*2] + lens[numseg-1]);

      CalCoor(&x1, &x2, x, from, to, pix_per_res, MaxSeqImgSize);
      x1 = MAX(x, x1);
      x2 = MIN(x2, right);

      y = GetY_nr_cddsrv(head, &end, from, to, seqlen, 5);
      y2 = y+FontBH+2;

      Dart_Acc2Info_external(CddName[i], shortname, defline, definition);

      if (ismap == true) {

	if (defline.size() == 254) 
		defline.replace(defline.size()-4, 3, "...");
	if (defline[defline.size()-1] != '.') defline += ".";

	unsigned pos = defline.find(';');
        if (pos != string::npos) def2 = defline.substr(0, pos) + ".";	
	else {
	     
	     pos = defline.find('.');
	     if (pos != string::npos) def2 = defline.substr(0, pos+1);
	     else def2 = defline + ".";
	}	     
	def2 = (string)CddName[i] + ":" + def2 + " Click for the CD alignment.";

	alinumseg = 0;
        cdalign = "";
	ostringstream output;
	for (j=0; j< numseg; j++) {
             index = 2*j;
             if (starts[index]!= -1 && starts[index+1] != -1) {

		 output.str("");
                 output << "," << starts[index+1] << "," << starts[index]
                                << "," << lens[j];
                 cdalign += output.str();

                 alinumseg++;
              }
        }

	output.str("");
        output << alinumseg << cdalign;
        strtmp = output.str();

	fprintf(File, "<area shape=rect coords=%d,%d,%d,%d ",
		x1, y, x2, y2);

	output.str("");
        output << "href=\"" << CDDurl
               << "?ascbin=2&maxaln=10&seltype=3&uid=" << PssmId[i]
               << "&aln=" << strtmp
               << "&querygi=" << gi
               << "&querynm=" << querynm;

	fprintf(File, (char *)output.str().c_str());
	fprintf(File, "\" alt=\"%s\"", def2.c_str());
	fprintf(File, "title=\"%s\">\n", def2.c_str());

	output.str("");
        output << "+OR+" << PssmId[i] << "[uid]";
        cTmp += output.str();

      }
      else {  	/* ismap == FALSE */

	gdImageRoundRectangle(im, x1, y, x2, y2, 8, 5, 
		iDartCol[iColor[i]], 1);
    
	len1 = shortname.size()*FontBW;
	len2 = shortname.size()*FontMW;
	color = white;
	if (iColor[i] == 2 || iColor[i] == 4 || iColor[i] == 0) color=black;
	if ((x2-x1)> len1) 
		gdImageString(im, gdFont7X13b, (x1+x2-len1)/2, y+1, 
				(char *)shortname.c_str(), color);
	else if ((x2-x1) > len2) 
		gdImageString(im, gdFont6X12, (x1+x2-len2)/2, y+1, 
			      (char *)shortname.c_str(), color);
	else {
		char_num = (x2-x1)/FontBW;
	         if (char_num >= 3) 
			shortname.replace(shortname.size() -3, 3, "...");
         	else switch (char_num) {
                	case 1: shortname =  "."; break;
                	case 2: shortname = "..";
              	}
	  	gdImageString(im, gdFont7X13b, x1+2, y+1, 
					(char *)shortname.c_str(), color);
	}
      }

      if (!dy) {
	dy = y0-y;
	y0 = y;
      }
    }

    if (!dy) dy=4;
    else dy = 10;

    x = x_ori;
    if (ismap == false)
     	   NameMapOrImg(ismap,&x,head->y+dy,(char*)"CDs",' ',0,0,0,NULL,NULL,0);
    else {
	   strtmp = "Click for Entrez Conserved Domains.";
	   fprintf(File, "<area shape=rect coords=%d,%d,%d,%d ",
		x, head->y+dy, x+StrLen("CDs")*FontBW, head->y+dy+FontBH);
	   fprintf(File, "href=\"%s?%s%s\" alt=\"%s\" title=\"%s\">\n", 
	     ENTREZurl, CDDOpt, cTmp.c_str()+4, strtmp.c_str(), strtmp.c_str());
    }

    FreeOverLoc(head);
    delete [] iColor;
    delete [] PssmId;
    delete [] dsp;
    for (i=0; i< CdNum; i++) delete [] CddName[i];
    delete [] CddName;
    delete [] iClus;
  }

} /* end of QueryMapOrImg */


#define MaxDipLen 120

void ImgMapOrDraw(bool ismap, unsigned Fsid, VastPageDataPtr vpp, 
	unsigned numhitsdisplayed, unsigned iKept, char * selnbrstring, 
	char * selsdidstring, SortBy sortby,
	SubSetID subsetnum, unsigned pagenum, unsigned StartLoc, FILE *File)
{
  char          stats[100], str[100], a_pdbId[6];
  char		a_chainLett;
  DomIntvlData	dip[MaxDipLen];
  float         f;
  unsigned          i, a_domNo, a_gi=0, a_seqlen, b_seqlen;
  unsigned          x, y, a_domintvl_no=0;
  static char SortBy_name[6][10] ={"Ali_Len", "Score", "P_value",
                                "Rmsd", "Ali_Res.", "%Id."};
  static char altstr[6][100]={"Number of Aligned Residues", 
				"Vast Score",
				"Vast P_value",
				"Room-mean square deviation",
                                "Number of Aligned Residues",
				"Identity"};

 // StrCut(a_pdbId, vpp[0].aDomName, 1, 4);
  sprintf(a_pdbId, "Chain");
  a_chainLett = vpp[0].aDomName[4];
  a_domNo = (unsigned) atol(vpp[0].aDomName+5);

  if (JobID == "") {
  	a_domintvl_no = 
	    constructDomIntvlData(aMmdbId, aChnNo, dip, MaxDipLen);

  	a_seqlen = dip[0].tu;		
  	a_gi = constructGi(aMmdbId, aChnNo);
  }
  else  {		// new VS
	a_domintvl_no = GetVSDomIntvlData(JobID, aChnNo, dip, MaxDipLen);
	a_seqlen = dip[0].tu;
  }

/*
  else {
	AsnIoPtr		aipr=NULL;
	BiostrucPtr		a_bstp = NULL;
        BioseqPtr		bsqp = NULL; 
	BioseqSetPtr		bssp = NULL;
	BiostrucResidueGraphSetPtr    stdDictionary;
        char 			AsnPath[MAX_TBUFF];
	ObjectIdPtr		objidp = NULL;
	SeqEntryPtr		sep = NULL;
 	SeqIdPtr		sip = NULL;
   	BiostrucFeatureSetPtr 	bfsp = NULL;
        ValNodePtr      	vnp = NULL, vnp1, vnp2;
	BiostrucFeaturePtr	bfp = NULL;	
	ResidueIntervalPntrPtr  ripp;
	unsigned 			chnNo, domid;

	sprintf(AsnPath, VSPATH);
	StringCat(AsnPath, JobID.c_str());
	StringCat(AsnPath, "/b");
	StringCat(AsnPath, a_pdbId);

	a_bstp = FetchBS(AsnPath, 0, ONECOORDRES, 1, POWER_VIEW);

	aipr = AsnIoOpen((char *)"bstdt", (char *)"rb");
	stdDictionary = BiostrucResidueGraphSetAsnRead(aipr, NULL);
	AsnIoFlush(aipr);
    	aipr = AsnIoClose(aipr);

	sep = (SeqEntryPtr) MakeBioseqs(a_bstp, stdDictionary);

	if (!sep) 
	     PrtMesC(MAILto, "VASTSRV (VastGraph, VS1)", 
			"Unable to get SeqEntry", "", false);

	if (sep->choice == 2){
	   i=1;

	   bssp = (BioseqSetPtr)sep->data.ptrvalue;
	   sep  = (SeqEntryPtr) bssp->seq_set;
	   while (i< aChnNo) {
		sep = sep->next;
		i++;
	   }
	   if (!sep)  
		PrtMesC(MAILto, "VASTSRV (VastGraph, VS2)",
                        "Unable to get SeqEntry", "", false);

	}
 	bsqp = (BioseqPtr)sep->data.ptrvalue;	
	sip = bsqp->id;
	while (sip) {
	     if (sip->choice == SEQID_LOCAL) {
	     	sprintf(str, "%s %c %d", a_pdbId, a_chainLett, a_domNo);
	     	objidp = (ObjectIdPtr)sip->data.ptrvalue;
	        if (!StringNCmp(objidp->str, str, 6)) {
		    a_seqlen = bsqp->length; 
		    break;
	      	}
	     }
   	     sip = sip->next;
	}

	dip[0].domNo = 0;
	dip[0].frm = 1;
	dip[0].tu = a_seqlen;
	a_domintvl_no = 1;
        domid = 1;
 	for (bfsp = a_bstp->features; bfsp != NULL; bfsp = bfsp->next)  {
	    for (vnp=bfsp->descr; vnp != NULL; vnp=vnp->next)
	        if (vnp->choice == BiostrucFeatureSetDescr_name &&
                        !StringCmp((char *)vnp->data.ptrvalue,(char *)"NCBI Domains")) {
		    for (bfp = bfsp->features; bfp!=NULL; bfp=bfp->next) {
			vnp=bfp->Location_location;
			vnp1 = (ValNodePtr) vnp->data.ptrvalue;
			vnp2 = (ValNodePtr) vnp1->data.ptrvalue;

			ripp = (ResidueIntervalPntrPtr) vnp2->data.ptrvalue;
                	chnNo = ripp->molecule_id;
			if (chnNo != aChnNo) continue;
           
			for (; ripp!= NULL; ripp = ripp->next) {

			      dip[a_domintvl_no].domNo = domid;
			      dip[a_domintvl_no].frm = ripp->from;
			      dip[a_domintvl_no++].tu = ripp->to;
		        }
			domid++;
		    }
	        }	
	}
	
	BiostrucFree(a_bstp);
  }
*/

  pix_per_res= (float)MaxSeqImgSize/(float)a_seqlen;

  QueryBit = new char [a_seqlen+1];
  for (i=0; i< a_seqlen; i++) QueryBit[i] = '0';

  y = StartLoc+3;  /* the start position of slave */

  x = GraphWidth - StrLen(SortBy_name[sortby])*FontBW;
  if (ismap == TRUE) {

    sprintf(str, "%s. Click for a printable table.", altstr[sortby]);
    
    fprintf(File, "<area shape=rect coords=%d,%d,%d,%d ",
	     x, 25, GraphWidth, 25+FontBH);   /* the position of status label */
    fprintf(File, "href=");
    fprintf(File, ParURL, URLcgi, CGIname, aSdi, sortby);
    if (selnbrstring == NULL && selsdidstring == NULL) {
	if (pagenum)
            fprintf(File, PageSubsetURL, pagenum, subsetnum, subsetnum);
	else fprintf(File, "schsub=Find");
    }
    else {
	fprintf(File, "schsub=Find");
	if (selnbrstring) fprintf(File, "&selnbr=%s", selnbrstring);
	if (selsdidstring) fprintf(File, "&selsdid=%s", selsdidstring);
    }

    if (iKept) 
      	for (i=0; i< iKept; i++)
       		fprintf(File, "&hit=%d", getfid(vpp[i].bSdi, vpp[i].bAlignId));
    if (JobID != "") {
	if (ReqId) fprintf(File, "&reqid=%llu", ReqId);
	else fprintf(File, "&vsid=%s&pass=%s", JobID.c_str(), Passwd.c_str());
    }

    fprintf(File, "&table=y\" alt=\"%s\" title=\"%s\">\n", str, str);
  }
  else 
	NameMapOrImg(ismap, &x, 25, SortBy_name[sortby], ' ', 0, 0, 0, 
			NULL, NULL, 1);

  x = 10;
  for (i=0; i< numhitsdisplayed; i++) { 

      b_seqlen = constructChainLength(SdiToMmdbId(vpp[i].bSdi), 
						SdiToChainNo(vpp[i].bSdi));

      if (!sortby || sortby==4) sprintf(stats, "%d", vpp[i].nres);
      else {
	if (sortby == 1) f = (float) vpp[i].vScore;
	else if (sortby == 2)f = (float) vpp[i].mlogp;
	else if (sortby == 3) f = (float) vpp[i].rmsd;
	else f = (float) vpp[i].pcntId;
	f /= (float) ASP_SCALE_FACTOR;
	if (sortby == 1 || sortby == 3) sprintf(stats, "%.1f", f);
	else if (sortby == 5)sprintf(stats, "%.1f", f*100.0);
	else {

          /* adjust for database size */
	  f -= LOG10_500;
	  if (f <= 4.0) {
	    f = (float) exp(-LOG_10*f);
	    sprintf(stats, "%.4f", f);
	  }
	  else sprintf(stats, "10e-%.1f", f);
	}
      } 

    AlignmentMapOrImg(ismap, x, y, Fsid, i, b_seqlen, vpp, stats, 
		StrLen(SortBy_name[sortby]), File);
    y += 30;

  }  

  /* draw query domain as a ruller */
  
  y = 30;     /* the start position of master seq. ruler */
  QueryMapOrImg(ismap, x, y, a_pdbId, Fsid, a_chainLett, a_domNo, 
		a_gi,a_seqlen, dip, a_domintvl_no, vpp, numhitsdisplayed, File);

}  /* ImgMapOrDraw */



void DrawStrucNbr(unsigned Fsid, VastPageDataPtr vpp, unsigned numhitsdisplayed,
		unsigned iKept, char * selnbrstring, char * selsdidstring,
		SortBy sortby, 
		SubSetID subsetnum, unsigned pagenum, short ImageSize)
{

	unsigned StartLoc;

	printf("Content-type: image/gif\n\n");
	im = gdImageCreate(GraphWidth, ImageSize);
	white   = gdImageColorAllocate(im, 255, 255, 255);
	black   = gdImageColorAllocate(im,   0,   0,   0);
        blue    = gdImageColorAllocate(im, 0, 0, 255);
	red     = gdImageColorAllocate(im, 255, 0, 0);
        gray 	= gdImageColorAllocate(im, 102, 102, 102);

/* Dart Color scheme */
	iDartCol[5]  = red;
	iDartCol[1]  = blue;
	iDartCol[2]  = gdImageColorAllocate(im, 0, 255, 0);   /* green */
	iDartCol[3]  = gdImageColorAllocate(im, 204,102,  0); /* orange */
	iDartCol[4]  = gdImageColorAllocate(im, 204,204,  0); /* yeller */
	iDartCol[0]  = gdImageColorAllocate(im, 153,204,255); /* sky b */
	iDartCol[6]  = gdImageColorAllocate(im, 102,153,  0); /* spring */
	iDartCol[7]  = gdImageColorAllocate(im, 204,102,153); /* lavender */
	iDartCol[8]  = gdImageColorAllocate(im,   0,204,204); /* cyan */
	iDartCol[9]  = gdImageColorAllocate(im, 153,153,  0); /* brown */
	iDartCol[10] = gdImageColorAllocate(im, 153, 51,255); /* violet */
	iDartCol[11] = gdImageColorAllocate(im,   0,153,153); /* blue-green */
	iDartCol[12] = gdImageColorAllocate(im,   0,204,102); /* teal */

	StartLoc = ImageSize - numhitsdisplayed*30;
	ImgMapOrDraw(FALSE, Fsid, vpp, numhitsdisplayed, iKept, selnbrstring, 
		selsdidstring, sortby, subsetnum, pagenum, StartLoc, NULL);

	gdImageGif(im, stdout);
        gdImageDestroy(im);

} /* DrawStrucNeig */

