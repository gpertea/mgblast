#ifndef _VASTUTI_
#define _VASTUTI_


/*
 * $Id: vastuti.hpp,v 1.1 2005/07/26 17:12:50 chenj Exp $
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
 * $Log: vastuti.hpp,v $
 * Revision 1.1  2005/07/26 17:12:50  chenj
 * Making linux VSNbr.cgi
 *
 *
 *
 */


#include "vastgenDB.hpp"

#include <string>

#include <ncbi.h>
#include <accentr.h>
#include <netentr.h>
#include <ncbiwww.h>
#include <sys/resource.h>
#include <asn.h>
#include <accutils.h>
#include <mmdbapi.h>
#include "vastlocl.h"
#include "mmdblocl.h"
#include "mmdbdata.h"
#include <objcdd.h>
#include <www.h>

void PrtMesC(string MAILto, CharPtr str1, CharPtr str2, CharPtr str3, bool ret);

void loadSubsetInfo(void);

SeqEntryPtr MakeBioseqs(BiostrucPtr bsp, BiostrucResidueGraphSetPtr stdDictionary);

void PrintAlignViewBanner(FILE *table, unsigned iFSID, VastPageDataPtr vpp,
        unsigned numhitsdisplayed);

void PrintHitsSortBanner(FILE *table, SortBy sortby, SubSetID subsetnum,
       unsigned pagenum, unsigned numpages, char cTable);

void PrintSearchNbr(FILE *table);

void PrintQueryInfo(FILE *table, char * pcPDB, char cChain, int iDomain);

void ImgMapOrDraw(bool ismap, unsigned Fsid, VastPageDataPtr vpp,
        unsigned numhitsdisplayed, unsigned iKept, char * selnbrstring,
        char * selsdidstring, SortBy sortby,
        SubSetID subsetnum, unsigned pagenum, unsigned StartLoc, FILE *File);


void DrawStrucNbr(unsigned Fsid, VastPageDataPtr vpp, unsigned numhitsdisplayed,
                unsigned iKept, char * selnbrstring, char * selsdidstring,
                SortBy sortby,
                SubSetID subsetnum, unsigned pagenum, short ImageSize);

unsigned GetNumberOfDomains(void);

unsigned BelongsToSubset(char *domid, unsigned sub, unsigned *grpnum, unsigned *grank);

void VastToCn3DAndAli(WWWInfoPtr www_info);

void loadSubsetInfo();

string GetSubsetName_DB(unsigned subsetId);

void WWWPrintFileData(char * FName,  FILE *pFile);

BiostrucAnnotSetPtr LocalGetBiostrucAnnotSet(unsigned mmdbid);

BiostrucAnnotSetPtr LocalGetFeatureSet(unsigned mmdbid,unsigned feature_set_id);

short Check_VastSearch_Password();

BiostrucAnnotSetPtr PruneBiostrucAnnotHits(BiostrucAnnotSetPtr basp, unsigned FSID, ValNodePtr pvnFids);

void MakeNbrList(char * nbrString,  char **(*SelNames), unsigned **SelSds,
	unsigned *iCount, unsigned ischar);

void 
OrderCopyVpp(VastPageDataPtr vpptmp, VastPageDataPtr vpp, unsigned iKept, unsigned * KepBsfId);


#endif
