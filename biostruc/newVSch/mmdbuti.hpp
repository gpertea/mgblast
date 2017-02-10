/*
* $Id: mmdbuti.hpp,v 1.1 2005/07/26 17:11:46 chenj Exp $
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
* $Log: mmdbuti.hpp,v $
* Revision 1.1  2005/07/26 17:11:46  chenj
* Making linux VSMmdb.cgi
*
* Revision 1.3  2003/01/15 16:15:09  chenj
* change the definition of ChainScaleMapOrImg()
*
* Revision 1.2  2003/01/14 19:57:17  chenj
* Minor changes
*
* Revision 1.1.1.1  2002/12/04 21:12:08  chenj
* Imported sources
*
*
*
*===========================================================================
*/

#ifndef _MMDBUTI_H
#define _MMDBUTI_H

#include <cgi/cgictx.hpp>
#include "hUtilib.hpp"

#include <ncbi.h>
#include <mmdbdata.h>
#include <gifgen.h>
#include <mmdbapi.h>

#define MAX_TBUFF    8192
#define ChainPerImg  10

#define MaxSeqImgSize   700
#define GraphWidth      800

#define ncycle 	10


using namespace ncbi;
using namespace std;

class ImgInfo
{
  public:
		ImgInfo() {};
		~ImgInfo() {};
		
		static 	unsigned 	id;
		static 	unsigned 	chbeg;
		static 	unsigned 	chend;
		static 	unsigned 	imgsize;
		static 	string 		subset;
		static 	bool 		ismap;
		static 	int		white;
		static 	int		black;
		static 	int		blue;
		static 	int		red;
		static 	int		gray;
		static 	int		maxseqlen;
		static 	double		pix_per_res;
		static 	vector<int> 	DartColr;
		static 	vector<int> 	Cn3DColr;
		static	unsigned 	prot_ch_cnt;
		static	unsigned	nucl_ch_cnt;

};


class ImgConf
{
  public:
		ImgConf() {};
		~ImgConf() {};

		static 	string 		VastUrl;
		static 	string 		VastCgi;
		static 	string 		CddUrl;
        	static 	string 		CddCgi;
        	static 	string 		EntrezUrl;
        	static 	string 		EntrezCgi;

};



typedef struct  domainhead {

        unsigned short    colidx, thisdomain;
	unsigned	domCumid;
        ResidueIntervalPntrPtr  ripp;
        struct  domainhead * next;

} IntervalHead;

unsigned short ModelMapOrImg(gdImagePtr im, PMSD pmsd, IntervalHead **DomHead, unsigned *DomIdx, CCgiContext& ctx, bool VastLink);

unsigned short ChainNameMapOrImg(gdImagePtr im, unsigned short x,
			unsigned short y, unsigned short protein);

void DrawImg(CCgiContext& ctx);

void ChainScaleMapOrImg(gdImagePtr im, unsigned short x, unsigned short y, 
		PMMD pmmd, unsigned sdid, unsigned color, unsigned labelcol,
		CCgiContext& ctx, bool VastLink);


unsigned short ChainDomMapOrImg(gdImagePtr im, unsigned short  x, 
		unsigned short  y, PMMD pmmd, unsigned * DomIdx, 
		IntervalHead **DomHead, CCgiContext& ctx, bool VastLink);

short MaxSeqLenProOrDRna(PMSD pmsd, unsigned short protein);

PDNMM PdnmmforChainX(PDNMM pdnmm, unsigned short chainx);

void GetDomFeaPtr(BiostrucPtr bsp, bool *hasDomain, BiostrucFeaturePtr *domain_bfp);

unsigned short ChainCount(PMSD pmsd);

void GroupingChains(PMSD pmsd);

void CalDomIdx(PDNMM pdnmm, unsigned **DomIdx);

void CheckDomColIdx(BiostrucFeaturePtr bfp, PDNMM pdnmm, IntervalHead **DomHead,
unsigned* DomIdx, bool hasDomain, bool forImg);


#endif
