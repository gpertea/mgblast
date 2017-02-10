/*
* $Id: mmdbuti.cpp,v 1.1 2005/07/26 17:11:46 chenj Exp $
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
* Source file for VSMmdb.cgi
*
*
* $Log: mmdbuti.cpp,v $
* Revision 1.1  2005/07/26 17:11:46  chenj
* Making linux VSMmdb.cgi
*
* Revision 1.3  2003/01/15 16:13:42  chenj
* add uid as an argument to ChainScaleMapOrImg()
*
* Revision 1.2  2002/12/26 21:21:20  chenj
* use Dart_CdNum() to get all Cds in database
*
* Revision 1.1.1.1  2002/12/04 21:12:08  chenj
* Imported sources
*
*
*
*=============================================================================
*/



#include <corelib/ncbiexec.hpp>
#include <corelib/ncbifile.hpp>
#include "VastSrchUti.hpp"
#include "SHGlobal.hpp"
#include "mmdbuti.hpp"
#include "qman.hpp"
#include "dart_external.hpp"

#include <cddapi.h>
#include "dartutil.h"



int	*ColorIdx;

static string tmp_dir = "/tmp/";


// Initialization
unsigned 	ImgInfo :: id = 0;
unsigned 	ImgInfo :: chbeg = 0;
unsigned 	ImgInfo :: chend = 0;
unsigned 	ImgInfo :: imgsize = 0;
string 	 	ImgInfo :: subset;
bool	 	ImgInfo :: ismap = true;
int	 	ImgInfo :: white = 0;
int	 	ImgInfo :: black = 0;
int	 	ImgInfo :: blue = 0;
int	 	ImgInfo :: red = 0;
int	 	ImgInfo :: gray = 0;
int	 	ImgInfo :: maxseqlen = 0;
double	 	ImgInfo :: pix_per_res = 0;
vector<int>	ImgInfo :: DartColr;
vector<int>	ImgInfo :: Cn3DColr;
unsigned	ImgInfo :: prot_ch_cnt = 0;
unsigned	ImgInfo :: nucl_ch_cnt = 0;
	


string ImgConf :: VastUrl;
string ImgConf :: VastCgi;
string ImgConf :: CddUrl;
string ImgConf :: CddCgi;
string ImgConf :: EntrezUrl;
string ImgConf :: EntrezCgi;


using namespace SHProjNS;
using namespace ncbi;

void 
ChainScaleMapOrImg(gdImagePtr im, unsigned short x, unsigned short y, 
	PMMD pmmd, unsigned sdid, unsigned  color, unsigned labelcol, 
	CCgiContext& ctx, bool VastLink)
{
    Char    *chain;
    string  cTmp;
    unsigned short 	x0, x1, y1, x2, y2, right, seqlen;
    unsigned    i, len, tick, ntick, dont; 
    static  unsigned   ticksteps[14] ={5, 10, 20, 25, 50, 100, 200, 150, 500,
				  1000, 2000, 2500, 5000, 10000};
    CCgiResponse& resp = ctx.GetResponse();
    static ImgConf thisConf;
    static ImgInfo Iinfo;

    seqlen = pmmd->iResCount;
    chain = pmmd->pcMolName;
    right = x2 = (unsigned)(x+ (seqlen-1)*Iinfo.pix_per_res);
    y2 = y + FontBH + 2;

    if (Iinfo.ismap)  {

 	if (!VastLink) cTmp = ToString(seqlen) + " residues.";
        else 
	   cTmp = ToString(seqlen) + 
			" residues, click to see its structure neighbors.";

        resp.out() << "<area shape=rect coords=" << x
		<< "," << y << "," << x2 << "," << y2;
        if (VastLink) {

	    resp.out()<< " href=\"" << thisConf.VastUrl << thisConf.VastCgi
	    		<< "?reqid=" << JobId2GrpId(Iinfo.id)
			<< "&subsetstr=" << Iinfo.subset.substr(0,3)
	    		<< "&sdid=" << sdid << "\"";
       	}
        resp.out() << " alt=\"" << cTmp << "\" title=\"" << cTmp << "\">\n";

    }
    else {
        gdImageFilledRectangle(im, x, y, x2, y2, color);
        cTmp = (string)"Chain " + chain;

        if ((right - x) >= (int)cTmp.size()*FontBW) {
	    x1 = (x + right - cTmp.size()*FontBW)/2;
	    gdImageString(im, gdFont7X13b,x1,y+1,(char*)cTmp.c_str(), labelcol);
        }
        else if ((right-x) >= (int)cTmp.size()*FontMW) {
	    x1 = (x + right - cTmp.size()*FontMW)/2;
	    gdImageString(im, gdFont6X12,x1,y+1,(char*)cTmp.c_str(), labelcol);
        }
        else if ((right-x) >= (int)strlen(chain)*FontBW) {
	    x1 = (x + right - strlen(chain)*FontBW)/2;
	    gdImageString(im, gdFont7X13b, x1, y+1, chain, labelcol);
        }
        else if ((right-x) >= (int)strlen(chain)*FontMW) {
	    x1 = (x + right - strlen(chain)*FontMW)/2;
	    gdImageString(im, gdFont6X12, x1, y+1, chain, labelcol);
        }

	/* add ticks */
       	gdImageLine(im, x, y, x, y-5, Iinfo.red);    /* for "1"  */
        gdImageLine(im, x2, y, x2, y-5, Iinfo.red);    /* for "len" */

        y1= y- (5 + FontMH);
        gdImageString(im, gdFont6X12, x-2, y1, (char *)"1", Iinfo.red);   // "1"
        x0 = x+FontMW;
        cTmp[0] = '\0';
        cTmp = ToString(seqlen);
        x2 -= (cTmp.size()-1)*FontMW;
        if (x2 > x)
           gdImageString(im, gdFont6X12, x2, y1, (char*)cTmp.c_str(),Iinfo.red);
								   /* "len" */
        if ((float)(seqlen-1)/(float)(Iinfo.maxseqlen-1)< 0.05) return;

        dont = x2 - FontMW;

        for (i=0; i< 14; i++) {
	    ntick = seqlen/ticksteps[i];
	    if (ntick < 10)  {
		tick = ticksteps[i];
		break;
	    }
        }

        for (i=1; i<= ntick; i++) {
             x1 = (unsigned)(x + ((float)(i*tick)-1.0)*Iinfo.pix_per_res);
             if (x1 > right-FontMW)  continue;
             y1 = y; y2 = y-5;
             gdImageLine(im, x1, y1, x1, y2, Iinfo.red);
             cTmp = ToString(i*tick);
             len = cTmp.size() * FontMW/2;
             x2  = x1 + len;
             x1 -= len;
             y1 -=  (5+ FontMH);
             if (x2 < dont && ((x2-x0) > (int)len)) {
                gdImageString(im, gdFont6X12, x1, y1, (char*)cTmp.c_str(), 
								Iinfo.red);
                x0 = x2+FontMW;
            }
        }
     }

} /* end of ChainScaleMapOrImg */



unsigned short ChainDomMapOrImg(gdImagePtr im, unsigned short x, 
		unsigned short y, PMMD pmmd, unsigned * DomIdx, 
		IntervalHead **DomHead, CCgiContext& ctx, bool VastLink)
{
    string                  cTom;
    unsigned                    x1, x2, y1, y2, ytmp = 0, idx, chnNo;
    unsigned                    colidx, from, to, seqlen, color=0; 
    unsigned		    right, labelcol=0;
    ResidueIntervalPntrPtr  ripp;
    IntervalHead	    *tmp_head;
    static ImgInfo Iinfo;
    static ImgConf thisConf;
    CCgiResponse& resp = ctx.GetResponse();
	
    chnNo = (unsigned)pmmd->iChainId;
    idx = DomIdx[chnNo]; 
    seqlen = (unsigned)pmmd->iResCount;
    right = (unsigned)(x + (float)(seqlen-1)*Iinfo.pix_per_res);

    if (DomHead[idx] == NULL) { 	/* no 3d domain */
	if (!Iinfo.ismap) {
	    colidx = ColorIdx[idx];
	    color = Iinfo.Cn3DColr[colidx];

	    labelcol = Iinfo.white;
	    if (colidx==3 || colidx==4 || colidx==5 || colidx==7
		|| colidx==8 || colidx==9)
		labelcol = Iinfo.black;
	}

	ChainScaleMapOrImg(im, x, y, pmmd, Dom2Sdi(Iinfo.id, chnNo, 0), 
					color, labelcol, ctx, VastLink);

	ytmp = y+FontBH+4;
    }
    else {		/* having 3d domain */
	unsigned	i;

	x1 = 10;
	x2 = x1 + strlen("3d Domains")*FontBW;
	y1 = y + FontBH + 6;
	y2 = y1+ FontBH;
	if (!Iinfo.ismap) {
	    gdImageLine(im, x1, y2, x2, y2, Iinfo.blue);
	    gdImageString(im,gdFont7X13b, x1, y1, (char *)"3d Domains", 
								Iinfo.blue);
	}

	for (tmp_head = DomHead[idx], i=1; tmp_head != NULL;
	     tmp_head = tmp_head->next, i++) {
                colidx = tmp_head->colidx;
                ripp = tmp_head->ripp;

                if (i == 1) {
		    unsigned   lablecol;

		    if (Iinfo.ismap) color = lablecol = 0;
		    else {
			color = Iinfo.gray;
			lablecol = Iinfo.white;
		    }
		    ChainScaleMapOrImg(im, x, y, pmmd, 
			Dom2Sdi(Iinfo.id, chnNo, 0), color, lablecol, 
			ctx, VastLink);
                }

                y1 = y + FontBH + 6;
                y2 = y1 + FontBH + 2;
                ytmp = y2;

                for (; ripp != NULL; ripp = ripp->next) {
		    from = ripp->from;   /* domain range */
		    to = ripp->to;
		    CalCoor(&x1, &x2,x,from,to,Iinfo.pix_per_res,MaxSeqImgSize);
		    x1 = MAX(x1, x);
		    x2 = MIN(x2, right);
		    if (Iinfo.ismap) {
			
			cTom="Residues " + ToString(from) +" to " +ToString(to);
  			if (VastLink) 
				cTom += ", click for structure neighbors";

                        resp.out()<<"<area shape=rect coords="
				      <<x1 <<"," <<y1 <<"," <<x2 <<"," <<y2;
			if (VastLink) {
                            resp.out() <<" href=\"" << thisConf.VastUrl
				<< thisConf.VastCgi
				<< "?reqid=" << JobId2GrpId(Iinfo.id)
				<< "&subsetstr=" 
				<< Iinfo.subset.substr(0, 3)
			        << "&sdid="
			        <<Dom2Sdi(Iinfo.id, chnNo, tmp_head->domCumid);
			};
			resp.out() <<"\" alt=\"" << cTom <<"\" title=\""
			              <<cTom <<"\">\n";
		    }
		    else {
                        gdImageFilledRectangle(im, x1, y1, x2, y2,
				Iinfo.Cn3DColr[colidx]);

                        cTom = ToString(tmp_head->thisdomain);

                        if ((x2 - x1 - 4) > cTom.size()*FontBW) {
			    int  x3;

			    x3=(x1+x2 - cTom.size()*FontBW)/2;
			    if (colidx==3 || colidx==4 || colidx==5
				|| colidx==7 || colidx==8 || colidx==9)
				color = Iinfo.black;
			    else color = Iinfo.white;
			    gdImageString(im, gdFont7X13b, x3, y1+1,
				(char*)cTom.c_str(), color);
                        }
		    }
                }  /* ; ripp != NULL */

           }    /* tmp_head */
    }

    return(ytmp);
}


short MaxSeqLenProOrDRna(PMSD pmsd, unsigned short protein)
{
        unsigned short    tmp_maxseqlen = 0;
        PMMD    pmmd;
        PDNMM   pdnmm;

        tmp_maxseqlen = 0;
        for (pdnmm = pmsd->pdnmmHead; pdnmm != NULL; pdnmm = pdnmm->next)
        {
                pmmd = (PMMD) pdnmm->data.ptrvalue;
                if ((protein && (pmmd->bWhat & AM_PROT)) 
			|| (!protein && (pmmd->bWhat & (AM_DNA | AM_RNA)))) {
                        if (tmp_maxseqlen < pmmd->iResCount)
                                tmp_maxseqlen = pmmd->iResCount;
                }
        }

        return(tmp_maxseqlen);

}       /* end MaxSeqLen */




PDNMM PdnmmforChainX(PDNMM pdnmm, unsigned short chainx)
{
        unsigned short    i = 1;
        PDNMM   pdnmmtmp;
        PMMD    pmmd;

        pdnmmtmp = pdnmm;

        while (i <= chainx)
        {
                pmmd = (PMMD)(pdnmmtmp->data.ptrvalue);
                if ((pmmd->bWhat & (AM_PROT |AM_DNA | AM_RNA)) && pdnmm)
                        i++;
                if (i <= chainx) pdnmmtmp = pdnmmtmp->next;
        }

        return(pdnmmtmp);

}       /* end PdnmmforChainX */




void
GetDomFeaPtr(BiostrucPtr bsp, bool *hasDomain, BiostrucFeaturePtr *domain_bfp)
{
        BiostrucFeatureSetPtr bfsp = NULL;
        ValNodePtr      vnp = NULL;

        *domain_bfp = NULL;
        *hasDomain = FALSE;

        for (bfsp = bsp->features; bfsp != NULL; bfsp = bfsp->next)
        {
             for (vnp = bfsp->descr; vnp != NULL; vnp = vnp->next)
                  if (vnp->choice == BiostrucFeatureSetDescr_name &&
                        !StringCmp((char *)vnp->data.ptrvalue, (char *)"NCBI Domains"))
                  {
                        *hasDomain = TRUE;
                        *domain_bfp = bfsp->features;
                        break;
                  }
        }


}	/* end of GetDomFeaPtr */




unsigned short ChainCount(PMSD pmsd)
{
        unsigned short            cnt;
        PMMD            pmmd;
        PDNMM           pdnmm;

        for (pdnmm = pmsd->pdnmmHead, cnt=0; pdnmm != NULL; pdnmm =pdnmm->next)        {
                pmmd = (PMMD) pdnmm->data.ptrvalue;
                if ((pmmd->bWhat) & (AM_PROT |AM_DNA | AM_RNA))
                        cnt++;
        }

        return(cnt);

}       /* end ChainCount */



void GroupingChains(PMSD pmsd)
{

        unsigned short    i, cnt;
        PDNMM   prohead, proend, nuhead, nuend, pdnmm, pdnmmnext;
        PDNMM   elsehead, elseend;
        PMMD    pmmd;


        prohead = proend = nuhead = nuend = elsehead = elseend = NULL;
        pdnmm = pmsd->pdnmmHead;

        cnt = ChainCount(pmsd);

        for (i=0, pdnmm = pmsd->pdnmmHead; i< cnt; i++)
        {
                pmmd = (PMMD) pdnmm->data.ptrvalue;
                pdnmmnext = pdnmm->next;

                if (pmmd->bWhat & AM_PROT)
                {
                    if (!prohead) prohead = proend = pdnmm; 
		    else {
                        proend->next = pdnmm;
                        proend = pdnmm;
                    }
                }
                else if (pmmd->bWhat & (AM_DNA | AM_RNA)) {
                        if (!nuhead) nuhead = nuend = pdnmm;
                        else {
                                nuend->next = pdnmm;
                                nuend = pdnmm;
                        }
                }
                else i--;

                pdnmm->next = NULL;
                pdnmm = pdnmmnext;
        }


        if (!pdnmm)
                if (!elsehead) elsehead = pdnmm;
                else elseend->next = pdnmm;
        elsehead = elseend = NULL;  /* remove non_PROT,non_DNA and non_RNA */
	/* ??? */

        if (!prohead) {
                if (!nuhead) pmsd->pdnmmHead = elsehead;
                else {
                        pmsd->pdnmmHead = nuhead;
                        nuend->next = elsehead;
                }
        }
        else {
                pmsd->pdnmmHead = prohead;
                if (!nuhead) proend->next = elsehead;
                else {
                        proend->next = nuhead;
                        nuend->next = elsehead;
                }
        }

        ((PDNMM)pmsd->pdnmmHead)->last = NULL;

}       /* end GroupingChains */



void CalDomIdx(PDNMM pdnmm, unsigned **DomIdx)
{
	unsigned	cnt, chnNo, maxchnNo=0;
	PMMD	pmmd;
	PDNMM	tmp_pdnmm;

        for (tmp_pdnmm = pdnmm; tmp_pdnmm != NULL; tmp_pdnmm = tmp_pdnmm->next)
        {
              pmmd = (PMMD) tmp_pdnmm->data.ptrvalue;
              if (pmmd->bWhat & (AM_PROT | AM_DNA |AM_RNA)) {
                    chnNo = pmmd->iChainId;
                    if (maxchnNo < chnNo) maxchnNo = chnNo;
              }
        }

        maxchnNo ++;
        *DomIdx = NewDataType <unsigned>(maxchnNo);

        for (tmp_pdnmm = pdnmm, cnt=0; tmp_pdnmm != NULL;
                                        tmp_pdnmm = tmp_pdnmm->next, cnt++)
        {
             pmmd = (PMMD) tmp_pdnmm->data.ptrvalue;
             if (pmmd->bWhat & (AM_PROT | AM_DNA |AM_RNA))
                 (*DomIdx)[pmmd->iChainId] = cnt;
        }

}	/* CalDomIdx */



void
CheckDomColIdx(BiostrucFeaturePtr bfp, PDNMM pdnmm, IntervalHead **DomHead, 
unsigned *DomIdx, bool hasDomain, bool forImg)
{
	BiostrucFeaturePtr	tmp_bfp;
	unsigned short			colidx = 0, thisdomain;
	unsigned 			chnNo;
	IntervalHead		*tmp_head;
	PDNMM			tmp_pdnmm;
	PMMD			pmmd;
	ResidueIntervalPntrPtr	ripp;
	ValNodePtr		vnp, vnp1, vnp2;
	
	if (hasDomain) {
	   for (tmp_bfp = bfp; tmp_bfp != NULL; tmp_bfp = tmp_bfp->next)
	   {
		unsigned  idx;

                vnp = tmp_bfp->Location_location;
                vnp1 = (ValNodePtr) vnp->data.ptrvalue;
                vnp2 = (ValNodePtr) vnp1->data.ptrvalue;

                ripp = (ResidueIntervalPntrPtr) vnp2->data.ptrvalue;
		chnNo = ripp->molecule_id;

		colidx %= ncycle;
		idx = DomIdx[chnNo];
		 if (DomHead[idx] == NULL)  { 
		     thisdomain = 1;
		     if (forImg == TRUE) ColorIdx[idx] = -11; /*just a label */
		}
		else thisdomain++;

		tmp_head = NewDataType <IntervalHead> (1);
		tmp_head->colidx = colidx;
		tmp_head->ripp = ripp;
		tmp_head->next = DomHead[idx];
		tmp_head->thisdomain = thisdomain;
		tmp_head->domCumid = tmp_bfp->id;
		DomHead[idx] = tmp_head;

		colidx++;
	   }	/* tmp_bfp */
	}	/* if (hasDomain) */

	for (tmp_pdnmm = pdnmm; tmp_pdnmm != NULL; tmp_pdnmm = tmp_pdnmm->next)
	{
	     unsigned  idx;

	     pmmd = (PMMD) tmp_pdnmm->data.ptrvalue;
	     if (pmmd->bWhat & (AM_PROT | AM_DNA |AM_RNA)) {

		chnNo = pmmd->iChainId; 
		idx = DomIdx[chnNo];
		if (DomHead[idx] == NULL) { 
			colidx %= ncycle; 
			if (forImg==TRUE) ColorIdx[idx] = colidx;
			colidx++;
		}
	     }

	}	/* for (tmp_pdnmm) */


}	/* CheckDomColIdx */



unsigned short ChainNameMapOrImg(gdImagePtr im, unsigned short x, 
				unsigned short y, unsigned short protein)
{
    unsigned short        x1, x2, y2;
    string	tmpstr;
    static 	ImgInfo Iinfo;

    if (protein) tmpstr = "Protein";
    else tmpstr = "Nucleotide";

    x1 = x+10;
    x2 = x1 + tmpstr.size()*FontBW;
    y2 = y+FontBH;

    if (!Iinfo.ismap) {
        gdImageString(im,gdFont7X13b, x1, y, (char*)tmpstr.c_str(),Iinfo.blue);
//        gdImageLine(im, x1, y2, x2, y2, Iinfo.blue);
    }

    return(x2);
}



static unsigned short CDDMapOrImg(gdImagePtr im, unsigned short x, 
		unsigned short y, PMSD pmsd, PMMD pmmd, CCgiContext& ctx)
{
    char                **CddName;
    string 		cTmp, str;
    OverLoc             *end, *head;
    unsigned            right, i, j, dy=0, y0, x1, x2, y1, y2, gi, seqlen;
    unsigned            numseg, from, to, len1, len2, alinumseg, index;
    int           	*starts, *lens;
    SeqAnnotPtr         sap = NULL;
    SeqAlignPtr         salp = NULL;
    DenseSegPtr         *dsp;
    SeqEntryPtr         sep=NULL;
    BioseqPtr           bseqp=NULL;
    static ImgInfo	Iinfo;
    CCgiResponse& 	resp = ctx.GetResponse();
    static ImgConf	thisConf;

//    unsigned chnNo = pmmd->iChainId;
    gi = GetVSGi(Iinfo.id, pmmd->iChainId);
    
    if (!gi) return (y);

    seqlen = pmmd->iResCount;
    right = (unsigned)(x + (float)(seqlen-1)*Iinfo.pix_per_res);

    sep = (SeqEntryPtr) GetVSSeqEntryForGi(gi, false, NULL);

    if (!sep) return (y);
    bseqp = (BioseqPtr)sep->data.ptrvalue;

    if (!bseqp) {

	str = "<h3>SeqEntryPtr not NULL, but BioseqPtr is NULL and gi=" 
		+ ToString(gi) + ".</h3>\n";
	PrtMes::PrintErrorHeader(&resp, "VastSrch", 1);
      	PrtMes::PrintMsgWithTail(&resp, str);

      	ERR_POST("bseqp = NULL");
      	throw exception();
    }

//curtime = CTime(CTime::eCurrent).AsString();
//cerr << "before CddSynchronousQuery " << curtime << endl;
//    sap = CddSynchronousQuery(bseqp, 0.01, true, true, false, "", false);

sap = NULL;
/*
    if (sap) {

        unsigned        *PssmId;
        unsigned short  *iColor, CdNum, thisColor=0; 
	short		*iClus;
        string          querynm, shortname, defline, definition;
        string          def2, cdalign;

        end = head = NewOverLoc(seqlen);
        head->y = y0 = y+5;

        for (salp = (SeqAlignPtr)sap->data, CdNum=0; salp!=NULL;
             salp=salp->next, CdNum++);

        iColor = NewDataType <unsigned short> (CdNum);
        PssmId = NewDataType <unsigned> (CdNum);
	dsp = NewDataType <DenseSegPtr> (CdNum);
        CddName = NewDataType <CharPtr> (CdNum);
        for (j=0; j< CdNum; j++) {
                CddName[j] = NewDataType <char>(30);
                CddName[j][0] = '\0';
        }
        iClus = NewDataType <short> (CdNum);

	string strtmp;
        for (salp = (SeqAlignPtr)sap->data, i=0; salp!=NULL;
                                                salp=salp->next, i++) {
            dsp[i] = (DenseSegPtr)salp->segs;
            PssmId[i] = GetPSSMID(dsp[i]);
            iClus[i] = -1;

            if (PssmId[i]) {
                Dart_CDGi2Acc_external(PssmId[i], strtmp);
		sprintf(CddName[i], strtmp.c_str());
	    }
            else {
		THROWS("PssmId = 0");
	    }
        }

        for (i=0; i< CdNum; i++) {
             if (iClus[i] >= 0) continue;

             iClus[i] = i;
             if (!Iinfo.ismap) 
			iColor[i] = (thisColor++) % Iinfo.DartColr.size();
             for (j = i+1; j< CdNum; j++) {
                  if (PssmId[i] == PssmId[j]) {
                      if (!Iinfo.ismap) iColor[j] = iColor[i];
                      iClus[j] = i;
                  }
             }
        }

	strtmp = pmmd->pcMolName;
	querynm = "query";
 	if (strtmp != "") querynm += "_" + strtmp;

        for (i=0; i< CdNum; i++) {

            numseg = (dsp[i])->numseg;
            starts = (dsp[i])->starts;
            lens = (dsp[i])->lens;
            from  = starts[0] + 1;
            to = starts[(numseg-1)*2] + lens[numseg-1];

            CalCoor(&x1, &x2, x, from, to, Iinfo.pix_per_res, MaxSeqImgSize);
            x1 = MAX(x, x1);
            x2 = MIN(x2, right);

            y1= GetY_nr_cddsrv(head, &end, from, to, seqlen, 7);
            y2 = y1+FontBH+4;

            Dart_Acc2Info_external(CddName[i], shortname, defline, definition);

            if (Iinfo.ismap) {

                if (defline.size() == 254)
			defline.replace(defline.size()-4, 3, "...");
		if (defline[defline.size()-1] != '.') defline += ".";

        	unsigned pos = defline.find(';');
        	if (pos != string::npos) def2 = defline.substr(0, pos-1);
        	else {

             		pos = defline.find('.');
             		if (pos != string::npos) def2 = defline.substr(0, pos);
             		else def2 = defline + ".";
        	}
        	def2 = (string)CddName[i] + ":" + def2 
						+ "Click for the CD alignment.";

                alinumseg = 0;
                cdalign = "";
                for (j=0; j< numseg; j++) {
                    index = 2*j;
                    if (starts[index]!= -1 && starts[index+1] != -1) {

			cdalign += "," + ToString(starts[index+1]) + ","
                                + ToString(starts[index]) + ","
                                + ToString(lens[j]);
                        alinumseg++;
                    }
                }
		strtmp = ToString(alinumseg) + cdalign;

//
                resp.out() << "<area shape=rect coords="
			<< x1 << ","
			<< y1 << ","
			<< x2 << ","
			<< y2 << " "
			<< "href=\"" << thisConf.CddUrl << thisConf.CddCgi
			<< "?ascbin=2&maxaln=10&seltype=3&uid=" << PssmId[i]
			<< "&aln=" << strtmp 
			<< "&querygi=" 	
			<< "&querynm=" << querynm << "\" "
			<< "alt=\"" << CddName[i] << ": " << defline
			<< " Click to see CD alignment.\" title=\""
			<< CddName[i] << ": " << defline 
			<< " Click to see CD alignment.\">\n";
//


		cTmp += ("+OR+" + ToString(PssmId[i]) + "[uid]");

            }
            else {    // ismap == FALSE 

                gdImageRoundRectangle(im, x1, y1, x2, y2, 8, 5,
						Iinfo.DartColr[iColor[i]], 1);

                len1 = shortname.size()*FontBW;
                len2 = shortname.size()*FontMW;

                int color = Iinfo.white;
                if (iColor[i] ==2 || iColor[i] ==4 || iColor[i]==0) 
							color=Iinfo.black;
                if ((x2-x1-4)> len1)
                    gdImageString(im, gdFont7X13b, (x1+x2-len1)/2, y1+2,
                                  (char *)shortname.c_str(), color);
                else if ((x2-x1-4) > len2)
                    gdImageString(im, gdFont6X12, (x1+x2-len2)/2, y1+2,
                                  (char *)shortname.c_str(), color);
                else {
                    Int4  char_num;

                    char_num = (x2-x1-4)/FontBW;
                    if (char_num >= 3) 
			shortname.replace(shortname.size() -3, 3, "...");
                    else switch (char_num) {
                        case 1: shortname = "."; break;
                        case 2: shortname = "..";
                    }
                    gdImageString(im,gdFont7X13b,x1+4, y1+2, 
					(char*)shortname.c_str(), color);
                }
            }

            if (!dy) {
                dy = y0-y1;
                y0 = y1;
            }
        }       // for 

        y0 = y2;

        if (!dy) dy=6;
        else dy = 15;

        x1 = 10;
        x2 = x1+ StrLen("CDs") * FontBW;
        y1 = y+dy;
        y2 = y1+FontBH;
        if (!Iinfo.ismap) {
            gdImageLine(im, x1, y2, x2, y2, Iinfo.blue);
            gdImageString(im, gdFont7X13b, x1, y1, (char *)"CDs", Iinfo.blue);
        }
        else {
            strtmp = "Click to see Entrez Conserved Domains.";
	    resp.out() << "<area shape=rect coords=" << ToString(x1) << ","
			<< ToString(y1) << "," << ToString(x2) << ","
			<< ToString(y2) << " href=\"" 
			<< thisConf.EntrezUrl << thisConf.EntrezCgi
			<< "?db=cdd&term=" << cTmp.substr(3)
			<< "\" alt=\"" << strtmp 
			<< "\" title=\"" << strtmp << "\">\n";
        }
        y2 = y0;

        for (end = head; end != NULL; end = end->next) delete end;
        delete [] iColor;
        delete [] PssmId;
        delete [] dsp;
        for (i=0; i< CdNum; i++) delete [] CddName[i];
        delete [] CddName;
        delete [] iClus;
    }
    else y2 = y;
*/
	y2=y;  // because the above block.

    return(y2);

}  // end of CDDMapOrImg




unsigned short ModelMapOrImg(gdImagePtr im, PMSD pmsd, IntervalHead **DomHead,
		unsigned * DomIdx, CCgiContext& ctx, bool VastLink)
{
    PDNMM		pdnmm, pdnmm1;
    PMMD 		pmmd;
    unsigned 		i,colidx; 
    unsigned short 	x, y;
    static ImgInfo 	Iinfo;
    
    pdnmm = (PDNMM)pmsd->pdnmmHead;
    pdnmm1 = PdnmmforChainX(pdnmm, Iinfo.chbeg);
    pdnmm = pdnmm1;

// Protein Chains 

    Iinfo.maxseqlen = MaxSeqLenProOrDRna(pmsd, 1);
    Iinfo.pix_per_res = (float)MaxSeqImgSize/(float)(Iinfo.maxseqlen-1);
/*

    if (!Dart_InitReal()) {

        PrtMes::PrintErrorHeader(&(ctx.GetResponse()), "VastSrch", 1);
        PrtMes::PrintMsgWithTail(&(ctx.GetResponse()),"Dart_InitReal() failed");

        ERR_POST("Dart_InitReal() failed.");
        throw exception();
    }
*/

    if (Iinfo.maxseqlen) {
	for (i = Iinfo.chbeg; i <= Iinfo.chend; i++, pdnmm = pdnmm->next) {
   
	    pmmd = (PMMD)pdnmm->data.ptrvalue;
	    
	    if ((pmmd->bWhat) & (AM_DNA || AM_RNA)) break;
	    if ((pmmd->bWhat) & AM_PROT) {
		
		Iinfo.prot_ch_cnt++;

        // Name
		if (i == 1 || i == Iinfo.chbeg) y = 50;
		else y += 60;
		x = ChainNameMapOrImg(im, 0, y, 1);

	// Chain & Domains

	       x += 30; 
	       y=ChainDomMapOrImg(im,x,y, pmmd, DomIdx, DomHead, ctx, VastLink);

	// Cdd Image
		
		y = CDDMapOrImg(im, x, y, pmsd, pmmd, ctx);

	    }  // AM_PROT
	    
	}    // end for (i=chain1

    }  // end Protein Chains


// Nucleotide Chains
    Iinfo.maxseqlen = MaxSeqLenProOrDRna(pmsd, 0);
    Iinfo.pix_per_res = (float)MaxSeqImgSize/(float)(Iinfo.maxseqlen-1);

    if (Iinfo.maxseqlen) {
        for (i = Iinfo.chbeg, pdnmm = pdnmm1; i <= Iinfo.chend; 
						i++, pdnmm =pdnmm->next){

            pmmd = (PMMD)pdnmm->data.ptrvalue;

            if ((pmmd->bWhat) & (AM_DNA | AM_RNA)) {

            	unsigned  labelcol;
	
		Iinfo.nucl_ch_cnt++;

        // Name
                if (i==1 || i == Iinfo.chbeg) y = 50;
                else y += 80;
                x = ChainNameMapOrImg(im, 0, y, 0);

        // Chain Image

                if (!Iinfo.ismap) {
                    x += 9;
                    colidx = ColorIdx[DomIdx[pmmd->iChainId]];
                    labelcol = Iinfo.white;
                    if (colidx==3 || colidx==4 || colidx==5 || colidx==7
                        || colidx==8 || colidx==9)
                        labelcol = Iinfo.black;
                    ChainScaleMapOrImg(im, x, y, pmmd, 0, 
				Iinfo.Cn3DColr[colidx], labelcol, ctx, VastLink);
                }

            } 	// AM_DNA | AM_RNA 

        } 	// end nucleotide Chains

    }  	// end maxseqlen

//    Dart_FiniReal();

    return(y);

}	// end of ModelMapOrImg



void DrawImg(CCgiContext& ctx)
{
	bool		hasDomain = FALSE;
	BiostrucFeaturePtr	domain_bfp;
	BiostrucPtr	bsp=NULL;
	gdImagePtr      im;
	unsigned		i, cnt, *DomIdx;
	IntervalHead	*tmp_head;
	PDNMS		ModelStruc=NULL;
	PDNMM		pdnmm;
	PMSD		pmsd;
	IntervalHead	**DomHead;
	static 		ImgInfo Iinfo;

	
       /* 
	string dir = tmp_dir + "VS" + ToString(Iinfo.id) + "/";
	string bFile = dir + "biostr.txt";

        if (!CDir(dir).Exists()) {
            CExec::System((string("mkdir ") + dir).c_str());
            CExec::System((string("chmod 777 ") + dir).c_str());  // temp.!!!!
            string strtmp = "VS" + ToString(Iinfo.id);
	    DownloadBiostrFromDB(strtmp, bFile);
        }

	bsp = FetchBS((char*)bFile.c_str(), 1, 2, 1000, POWER_VIEW);
*/

	bsp = VSOpenBSP(Iinfo.id, 2, 10000);
 	if (!bsp) return;
 
	ModelStruc = MakeAModelstruc(bsp);
	pmsd = (PMSD) ModelStruc->data.ptrvalue;
        pdnmm = pmsd->pdnmmHead;
	cnt = ChainCount(pmsd);
	if (!cnt) return;

	ColorIdx = NewDataType<int>(cnt);
	DomHead = NewDataType <IntervalHead *> (cnt+1);
	for (i=0; i< cnt; i++) {
		ColorIdx[i] = -1;
		DomHead[i] = NULL;
	}
	CalDomIdx(pdnmm, &DomIdx);

	GetDomFeaPtr(bsp, &hasDomain, &domain_bfp);	
	CheckDomColIdx(domain_bfp, pdnmm, DomHead, DomIdx, hasDomain, TRUE);

        GroupingChains(pmsd);

        im = gdImageCreate(GraphWidth, Iinfo.imgsize);

 	Iinfo.white = gdImageColorAllocate(im, 255, 255, 255);
        Iinfo.black = gdImageColorAllocate(im,   0,   0,   0);
	Iinfo.blue =  gdImageColorAllocate(im,   0,   0, 225);
	Iinfo.red  =  gdImageColorAllocate(im, 255,   0,   0);
	Iinfo.gray = gdImageColorAllocate(im, 102, 102, 102);


	// Dart Color scheme
	Iinfo.DartColr.reserve( 13 );
	Iinfo.DartColr.push_back(gdImageColorAllocate(im,153,204,255)); // sky b
        Iinfo.DartColr.push_back(Iinfo.blue);
	Iinfo.DartColr.push_back(gdImageColorAllocate(im,  0,255,0));  // green 
        Iinfo.DartColr.push_back(gdImageColorAllocate(im,204,102,0));  // orange
        Iinfo.DartColr.push_back(gdImageColorAllocate(im,204,204,0));  // yellow
	Iinfo.DartColr.push_back(Iinfo.red);
        Iinfo.DartColr.push_back(gdImageColorAllocate(im,102,153,  0)); //spring
        Iinfo.DartColr.push_back(gdImageColorAllocate(im,204,102,153)); // lavender 
        Iinfo.DartColr.push_back(gdImageColorAllocate(im,  0,204,204)); // cyan
        Iinfo.DartColr.push_back(gdImageColorAllocate(im,153,153,  0)); // brown
        Iinfo.DartColr.push_back(gdImageColorAllocate(im,153, 51,255));// violet
        Iinfo.DartColr.push_back(gdImageColorAllocate(im,  0,153,153)); // blue-green
        Iinfo.DartColr.push_back(gdImageColorAllocate(im,  0,204,102)); // teal


	// Cn3D Color scheme
	Iinfo.Cn3DColr.reserve( 10 );
	Iinfo.Cn3DColr.push_back(gdImageColorAllocate(im,255, 0, 255));//magenta
	Iinfo.Cn3DColr.push_back(gdImageColorAllocate(im,  0, 0, 255)); // blue
	Iinfo.Cn3DColr.push_back(gdImageColorAllocate(im,139, 87, 66)); // brown
	Iinfo.Cn3DColr.push_back(gdImageColorAllocate(im,0,255,127)); //l. green
	Iinfo.Cn3DColr.push_back(gdImageColorAllocate(im,179,179,179));//l. gray
	Iinfo.Cn3DColr.push_back(gdImageColorAllocate(im,255,165,  0)); // gold
	Iinfo.Cn3DColr.push_back(gdImageColorAllocate(im,255,114, 86)); // pink
	Iinfo.Cn3DColr.push_back(gdImageColorAllocate(im,  0,255,  0)); // green
	Iinfo.Cn3DColr.push_back(gdImageColorAllocate(im,  0,255,255)); // cyan
	Iinfo.Cn3DColr.push_back(gdImageColorAllocate(im,255,236,139));//yel.tint

	
	printf("Content-type: image/gif\n\n");
	Iinfo.ismap = false;
	ModelMapOrImg(im, pmsd, DomHead, DomIdx, ctx, false);

//        int fd = ctx.GetResponse().GetOutputFD();
//        gdImageGif(im, fd);
	
        gdImageGif(im, stdout);
        gdImageDestroy(im);

	delete [] ColorIdx;
	for (i = 0; i< cnt; i++)
	  if (DomHead[i] != NULL)  {
		for (tmp_head = DomHead[i]; tmp_head != NULL; 
					tmp_head=tmp_head->next)
			delete [] tmp_head;
	  }
	delete [] DomHead;
	delete [] DomIdx;
	
//        CExec::System((string("rm -r ") + dir).c_str());

	return;

}	/* end DrawImg */
