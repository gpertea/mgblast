/* 
 * $Id: vast2cn3dDB_comb.cpp,v 1.1 2005/07/26 17:12:50 chenj Exp $
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
 * Author: Jie Chen, Lewis Geer, Chris Hogue, Siqian He
 * 
 *
 * $Log: vast2cn3dDB_comb.cpp,v $
 * Revision 1.1  2005/07/26 17:12:50  chenj
 * Making linux VSNbr.cgi
 *
 * Revision 1.3  2003/01/15 16:21:47  chenj
 * add construcMaxSdi()
 *
 * Revision 1.2  2003/01/14 20:42:32  chenj
 * change VastToCn3D() to VastToCn3DAndAli()
 *
 * Revision 1.1.1.1  2002/12/06 20:17:21  chenj
 * Imported Scouces
 *
 *
 *
 * This file fetches data from pubvast database, but for Vast Search job, it 
 * still reads data from files.
 *
 * =========================================================================
 *
 */


#define NLM_GENERATED_CODE_PROTO
#undef  DOWNLOAD_DEBUG_MODE

#include <corelib/ncbiexec.hpp>
#include <corelib/ncbifile.hpp>
#include <ctools/asn_converter.hpp>
#include "hUtilib.hpp"
#include "VastSrchUti.hpp"
#include "vastuti.hpp"
#include <accutils.h>
#include <mmdbapi1.h>
#include <objmime.h>
#include "cddalignview.h"
#include <objalign.h>
#include "vast2mage.h"
#include "mkbioseq.h"

#define CPUTIME_MAX 120

#define VIEW_STR_SUBBUT		0
#define VIEW_ALI_SUBBUT		1
#define VIEW_IN_Cn3DCACHE	'a'
#define VIEW_IN_Cn3D_GEN	'c'
#define VIEW_IN_HTML  		'h'
#define VIEW_IN_FASTA 		'f'
#define VIEW_IN_TEXT 		't'
#define VIEW_IN_HTML_PAGE	'H'
#define VIEW_IN_FASTA_PAGE	'F'
#define VIEW_IN_TEXT_PAGE	'T'


#define ALL_NBR 	0
#define CHECKED_NBR  	1

extern  void    WWWPrintFileData(CharPtr FName,  FILE *pFile);
extern 	string	MAILto, JobID, Passwd;
extern 	Char 	Database[PATH_MAX], VSPATH[PATH_MAX];
extern 	unsigned	aSdi, aMmdbId, aChnNo;
extern 	long long 	ReqId;

using namespace ncbi;

static void AddMMDBIdToBioseq(BioseqPtr bsp, unsigned mmdbid)
{
    DbtagPtr    dbtp;
    SeqIdPtr    sip_new = NULL;
    SeqAnnotPtr sap_new = NULL;

    if (bsp == NULL) return;

    sap_new= SeqAnnotNew();
    sap_new->type = 4;
    sip_new = (SeqIdPtr)ValNodeNew(NULL);
    sip_new->choice = SEQID_GENERAL;
    dbtp = (DbtagPtr)DbtagNew();
    dbtp->db = new char [10];
    sprintf(dbtp->db, "mmdb");
    dbtp->tag = (ObjectIdPtr) ObjectIdNew();
    dbtp->tag->id = mmdbid;
    sip_new->data.ptrvalue = dbtp;
    sap_new->data = sip_new;
    if (bsp->annot == NULL) bsp->annot = sap_new;
    else {
            sap_new->next = bsp->annot;
            bsp->annot = sap_new;
    }

    return;

}	/* AddMMDBIdToBioseq */




static void RemoveGi(BioseqPtr bsp) 
{
    SeqIdPtr preSid=NULL, thisSid = (SeqIdPtr)bsp->id;

    while (thisSid != NULL) {
	
 	if (thisSid->choice == SEQID_GI) {
            if (preSid) preSid->next = thisSid->next;
	    else bsp->id = thisSid->next;
	}	

	preSid = thisSid;
        thisSid = thisSid->next;
    }
}	/* RemoveGi */




/* Display a structural alignment in Cn3D and in HTML*/

void LIBCALL VastToCn3DAndAli(WWWInfoPtr www_info)
{
  char 		szName[5], AsnPath[PATH_MAX], AsnName[10]; 
  char		cViewType, pdbname_m[5], chain_m=' ';
  char 		*Name, *www_arg, *szTemp;
  short		iPDB = 0, complexity, indx;
  short		nbr_complexity;
  unsigned 	iFidCount=0, domNo_m, i, hits_num=0, iSubBut, *gi;
  unsigned 	aGi=0, Fid, Fsid, NumLabels, *BsfId;
  /* Fsid = aMmdbId * 10000 + aChnNo * 100 + aDomCumulNo , only used in VS job*/
  /* Fid = bMmdbId * 100000 + bChnNo * 1000 + bDomCumulNo * 10 + bAlignId */
  AsnIoPtr 		paiFile, aipr;
  BioseqPtr		biosp;
  BiostrucAnnotSetPtr 	pbsa = NULL, pbsaShort = NULL;
  BiostrucAlignSeqPtr	basp = NULL;
  BiostrucFeatureSetPtr pbsfs = NULL;
  BiostrucFeaturePtr 	pbsf = NULL;
  BiostrucResidueGraphSetPtr 	stdDictionary;
  BiostrucPtr 		pbsMaster=NULL, pbsSlave=NULL, pbsSlaveHead = NULL;
  BiostrucPtr		pbsSlaveTail;
  BiostrucAlignPtr 	pbsaStruct;
  BundleSeqsAlignsPtr 	bsap = NULL;
  NcbiMimeAsn1Ptr 	pvnNcbi;
  ObjectIdPtr		objid;
  SeqAnnotPtr 		psaAlignHead = NULL, psaAlignTail;
  SeqAlignPtr 		salpHead, salpTail;
  SeqEntryPtr 		sep;
  SeqIdPtr		sid;
  ValNodePtr 		pvnFids = NULL;
  VastPageDataPtr	vpp;

  /* action == 0 indicates MIME; action == 1 is text; action == 2 is save */
  if ((indx = WWWFindName(www_info, (char *)"action")) < 0) iPDB = 0;
  else {
    www_arg = WWWGetValueByIndex(www_info, indx);

    if (isInt(www_arg)) iPDB = (Int4) atoi(www_arg);
    else iPDB = 0;
  }


  if ( (indx=WWWFindName(www_info, (char *)"viewstr.x")) >=0
  		|| (indx=WWWFindName(www_info, (char *)"viewstr")) >=0 ) 
		iSubBut=VIEW_STR_SUBBUT;  
  else if ((indx=WWWFindName(www_info, (char *)"viewali.x")) >=0 
  		|| (indx=WWWFindName(www_info, (char *)"viewali")) >=0 ) 
		iSubBut=VIEW_ALI_SUBBUT;
  else {
	PrtMesC("", "VASTSRV (Vast2Cn3D)", "Please click either \"View 3D Alignment\" or \"View Sequence Alignment\".", "", false); 
  }


  if (iSubBut == VIEW_STR_SUBBUT) {
  	if ((indx = WWWFindName(www_info, (char *)"calltype")) >=0) {
		www_arg = WWWGetValueByIndex(www_info, indx);
		cViewType = www_arg[0];
  	}
	else 
	   PrtMesC(MAILto, "VASTSRV (Vast2Cn3D)", "No calltype code -- Please choose either Cn3D or Cn3D/Cache.\n", "", false); 

  }
  else if ((indx = WWWFindName(www_info, (char *)"alitype")) >=0) {
	CharPtr www_arg2;

	www_arg = WWWGetValueByIndex(www_info, indx);
        
  	if ((indx = WWWFindName(www_info, (char *)"nbr_complexity")) >=0) {
	    www_arg2 = WWWGetValueByIndex(www_info, indx);
	    if (isInt(www_arg2)) nbr_complexity = (Int2)atoi(www_arg2);
            else nbr_complexity = ALL_NBR;
	}
	else nbr_complexity = ALL_NBR;

   	cViewType = www_arg[0];
	if (nbr_complexity == ALL_NBR) {
	    if (cViewType == VIEW_IN_HTML) cViewType = VIEW_IN_HTML_PAGE;
	    else if (cViewType == VIEW_IN_FASTA) 
			cViewType = VIEW_IN_FASTA_PAGE;
	    else if (cViewType == VIEW_IN_TEXT) 
			cViewType = VIEW_IN_TEXT_PAGE;
	}
  }
  else 
    PrtMesC(MAILto, "VASTSRV (Vast2Cn3D)", "No alitype code -- Please choose Hypertext, Plain text or mFASTA.", "", true); 


  if ((indx = WWWFindName(www_info, (char *)"atm_complexity")) < 0)
    	complexity = ONECOORDRES; /* select alpha Carbons only by default */
  else {
    www_arg = WWWGetValueByIndex(www_info, indx);
    if (isInt(www_arg)) complexity = (Int2) atoi(www_arg);
    else complexity = ONECOORDRES;
  }
  if ((complexity != ONECOORDRES) && (complexity != ONECOORDATOM))
    complexity = ONECOORDRES; 
		/* bizarre value, but default to alpha-Carbons only */

  if (cViewType == 'H' || cViewType == 'F' || cViewType == 'T') {

	nbr_complexity = ALL_NBR;
        if ((indx = WWWFindName(www_info, (char *)"allbfid")) < 0)  {
              PrtMesC(MAILto, "VASTSRV (Vast2Cn3D)", "No \"allbfid\" was submitted!", "", true);
        }

        www_arg = WWWGetValueByIndex(www_info, indx);
        MakeNbrList(www_arg, NULL, &BsfId,  &hits_num, 0);
        gi = new unsigned [hits_num];
        for (i=0; i< hits_num; i++) {
	    Fid = BsfId[i];
	    if (JobID != "") ValNodeAddInt(&pvnFids, 0, Fid);
	}
	iFidCount = i-1;
  }
  else {
     nbr_complexity = CHECKED_NBR;

     if ((indx = WWWFindName(www_info, (char *)"hit")) <0) {
	if ((indx = WWWFindName(www_info, (char *)"defhit")) <0) {
	    PrtMesC(MAILto, "VASTSRV (Vast2Cn3D)", "No defhit was submitted!", "", false);
	}
	else {
	   www_arg = WWWGetValueByIndex(www_info, indx);
		
	   hits_num = 1;
           gi = new unsigned;
	   BsfId = new unsigned; 
	
	   if (isInt(www_arg))
		BsfId[0] = Fid = atol(www_arg);
           else {
                   PrtMesC(MAILto, (char *)"VASTSRV (Vast2Cn3D)", 
			(char *)"Non-numeric defhit: defhit = ", www_arg, TRUE);
           }

           if (cViewType == VIEW_IN_Cn3D_GEN || JobID != "") {
                 	ValNodeAddInt(&pvnFids, 0, Fid);
			iFidCount=1;
	   }
	}
     }
     else {	 /* loop over all the "hit" values in the list */

	 NumLabels = WWWGetNumEntries(www_info);

	 for (indx = 0; indx < NumLabels; indx++) {
		Name = WWWGetNameByIndex(www_info, indx);
		if (!StrICmp(Name, "hit")) hits_num++;
	 }

	 gi = new unsigned [hits_num];
	 BsfId = new unsigned [hits_num];
      
	 for (indx = 0; indx < NumLabels; indx++) {

	     Name = WWWGetNameByIndex(www_info, indx);
	     
	     if (!StrICmp(Name, "hit")) {
		 www_arg = WWWGetValueByIndex(www_info, indx);

		 if (isInt(www_arg)) Fid =atoi(www_arg);
		 else {
		   PrtMesC(MAILto,  (char *)"VASTSRV (Vast2Cn3D)",
			(char *)"Non-numeric hit: hit = .", www_arg, TRUE);
		 }
	
		 if (++iFidCount > 10 && cViewType == VIEW_IN_Cn3D_GEN) 
		 {
		     PrtMesC("", "VASTSrv (Vast2Cn3D)", "With general Cn3D, at most only 10 alignments can be displayed simultaneously. You may reselect neighbors or use Cn3D/Cache.", "", false);
			/* up to 10 boxes, ignore the rest. */
		 }
		 ValNodeAddInt(&pvnFids, 0, Fid);
		 BsfId[iFidCount-1] = Fid;
	     }
	 }
     }
  }


  /* Get BiostrucAnnotSet and SeqAnnot */

/*
  if (JobID != "" && !ReqId) {  // old VastSearch 
    OpenMMDBAPI((POWER_VIEW  // ^ FETCH_ENTREZ ), NULL);
    pbsa = LocalGetFeatureSet(aMmdbId, Fsid);
  }
  else
*/

    VastPageDataPtr vpptmp;

    if (hits_num > 1)
	vpptmp = new VastPageData [hits_num];
    vpp = new VastPageData [hits_num];
    if (JobID == "") {
    	if (hits_num ==1) 
    	    i=constructVastPagesByNbrs(vpp,BsfId,hits_num,aSdi, (SortBy)(-1));
    	else 
    	    i=constructVastPagesByNbrs(vpptmp,BsfId,hits_num,aSdi,(SortBy)(-1));
    }
    else {
	if (hits_num ==1) 
	    i=GetVSVppByNbrs(aSdi, vpp, BsfId, hits_num, (SortBy)(-1));
	else i=GetVSVppByNbrs(aSdi, vpptmp, BsfId, hits_num, (SortBy)(-1));

    }
    if (i != hits_num) 
	    PrtMesC(MAILto,"VASTSRV (Vast2Cn3D)","Error in getting VAST data.",
								 "", false); 
    if (hits_num > 1) {
      OrderCopyVpp(vpptmp, vpp, hits_num, BsfId);
      delete [] vpptmp;
    }
 
    if (JobID == "") pbsa = constructBASPFromVastPagePtr(vpp, hits_num, 0);
    else pbsa = constructBASPFromVastPagePtr(vpp, hits_num, 1);

    delete [] vpp;


  if (pbsa == NULL) {
    printf("Content-type: text/html\n\n");
    printf("<br><center><h2>VASTSrv Error (VastToCn3D):</h2><p>\n");
    printf("No alignment record exists for master mmdb_id = %d.</h3>\n",
                        aMmdbId);
    printf("Please alert \"%s\" of this problem.</center>\n",MAILto.c_str());
  }

   pbsaShort = pbsa;

{AsnIoPtr aipr;
aipr=AsnIoOpen("pbsa2.out", "w");
BiostrucAnnotSetAsnWrite(pbsaShort, aipr, NULL);
AsnIoClose(aipr);
}


  if (pbsaShort == NULL) 
     PrtMesC("","VASTSRV (Vast2Cn3D)","Can't find alignment record.","",false);
  pbsfs = pbsaShort->features;
  if (pbsfs) {
      pbsf = pbsfs->features;
      szTemp = pbsf->name;
      StrCut(szName, pbsf->name, 1, 4);
      sprintf(pdbname_m, szName);
      chain_m = szTemp[4];
      domNo_m = (Int4) atol (szTemp+5);
  }

  if (cViewType == VIEW_IN_Cn3DCACHE) bsap = BundleSeqsAlignsNew();
  else if (cViewType != VIEW_IN_Cn3D_GEN) basp = BiostrucAlignSeqNew();

  for (i=0; i< hits_num; i++) {
      unsigned bsdi;
      bsdi = BsfId[i]/100;
      gi[i] = constructGi(SdiToMmdbId(bsdi), SdiToChainNo(bsdi));
  }

  if (JobID == "") aGi = constructGi(aMmdbId, aChnNo);
  else aGi = GetVSGi(atoi(JobID.c_str()+2), aChnNo); // new VS

  /* Get Biostruc for Cn3d, and SeqEntry (Bioseq)  */

  if (cViewType==VIEW_IN_Cn3D_GEN || JobID!="") {

    pbsaStruct = BiostrucAlignNew();

    if (JobID == "") {     /* not a VS job */
	bool psok;
	pbsMaster = OpenBSP(aMmdbId, complexity, 1, TRUE, FALSE, FALSE, psok, Database);
    }
    else {
	string dir = "/tmp/" + string(JobID) + "_w/";
	string bFile = dir + "biostr.txt";
	if (!CDir(dir).Exists()) {
	    CExec::System((string("mkdir ") + dir).c_str());
      	    CExec::System((string("chmod 777 ") + dir).c_str());  // temp.!!!!
  	}
 	
	DownloadBiostrFromDB(JobID, bFile);
	pbsMaster = FetchBS((char*)bFile.c_str(),1, complexity, 1,POWER_VIEW);
  	CExec::System((string("rm -r ") + dir).c_str());
    }

    if (pbsMaster == NULL) 
      PrtMesC(NULL, "VASTSrv (Vast2Cn3D)", "Unable to load master structure from a Vast Search result.", "", false); 

    if (JobID == "") {
	sep = (SeqEntryPtr) constructSeqEntryForGi(aGi, TRUE, Database);
    }
    else if (ReqId) { 	/* new VS job */

	sep = (SeqEntryPtr) GetVSSeqEntryForGi(aGi, false, NULL);
	if (sep->choice == 1) { /* Bioseq */
	    RemoveGi((BioseqPtr)sep->data.ptrvalue);
	}
	else {
	
	     fprintf(stderr, "sep->choice != 1 sth wrong\n");
	     exit(1);
	}
    }

    if (sep == NULL)  
      PrtMesC(MAILto, "VASTSERV (Vast2Cn3D)", "Unable to get SeqEntry.", "", false);

    if (cViewType == VIEW_IN_Cn3D_GEN) {
	ValNodeLink(&(pbsaStruct->sequences), sep);
	if (JobID != "") {
		biosp = (BioseqPtr)sep->data.ptrvalue;
		sid = biosp->id;
		objid = (ObjectIdPtr) sid->data.ptrvalue;
 	}	
	
    }
    else if (cViewType != VIEW_IN_Cn3DCACHE) 
    		ValNodeLink(&(basp->sequences), sep);
    else 
	PrtMesC("", "VASTSRV (Vast2Cn3D)", "Can't work in Cn3D/Cache mode -- please try Cn3D mode.", "", false); 

	
  /* Make a linked list of Biostrucs of the slave structures */

    for (i=0; i< hits_num; i++) {
      Int4 bMmdbId;

      bMmdbId = SdiToMmdbId(BsfId[i]/100);        
      if (!pbsSlaveHead) {
	bool psok; 

	pbsSlaveHead = OpenBSP(bMmdbId, complexity, 1, TRUE, FALSE, FALSE, psok, Database);
	if (!pbsSlaveHead)  
	  PrtMesC("", "VASTSRV (Vast2Cn3D)", "Unable to load slave structure.", "", false);

      /* Make Bioseq for Slaves */
	sep = (SeqEntryPtr) constructSeqEntryForGi(gi[i], TRUE, Database);
	if (sep == NULL) 
	  PrtMesC(MAILto, "VASTSRV (Vast2Cn3D)","Unable to get SeqEntry.","",false);

        if (cViewType == VIEW_IN_Cn3D_GEN) 
                ValNodeLink(&(pbsaStruct->sequences), sep);
    	else if (cViewType != VIEW_IN_Cn3DCACHE) 
                ValNodeLink(&(basp->sequences), sep);
    	else 
           PrtMesC("", "VASTSRV (Vast2Cn3D)", "Can't work in Cn3D/Cache mode -- please try Cn3D mode.", "", false);

	pbsSlaveTail = pbsSlaveHead;
      }
      else { 
	bool psok;
	
	pbsSlave =OpenBSP(bMmdbId, complexity, 1, TRUE, FALSE, FALSE, psok, Database);
	if (!pbsSlave) 
		PrtMesC("", "VASTSERV (Vast2Cn3D)", "Unable to load slave structure.", "", false);

	sep = (SeqEntryPtr) constructSeqEntryForGi(gi[i], TRUE, Database);
	if (sep == NULL) 
		PrtMesC(MAILto, "VASTSRV (Vast2Cn3D)", "Unable to get SeqEntry.", "", false);

        if (cViewType == VIEW_IN_Cn3D_GEN) 
                ValNodeLink(&(pbsaStruct->sequences), sep);
        else if (cViewType != VIEW_IN_Cn3DCACHE) 
                ValNodeLink(&(basp->sequences), sep);
        else 
           PrtMesC("", "VASTSRV (Vast2Cn3D)", "Can't work in Cn3d/Cache mode -- please try Cn3D mode.", "", false);

	pbsSlaveTail->next = pbsSlave;
	pbsSlaveTail = pbsSlaveTail->next;
	pbsSlaveTail->next = NULL;
      }
    }
  }
  else {     /* cViewType != Cn3D_IN_GEN */

    BioseqPtr 		bsp;

/* Get SeqEntryPtr */

    if (JobID == "")  	/* not a VS job */
    	sep = constructSeqEntryForGi(aGi, TRUE, Database);
    else sep = GetVSSeqEntryForGi(aGi, false, NULL);
    if (sep == NULL) 
       PrtMesC(MAILto,"VASTSRV (Vast2Cn3D)","Unable to get Bioseq.", "", false);

    bsp = (BioseqPtr)sep->data.ptrvalue;
    if (JobID == "") AddMMDBIdToBioseq(bsp, aMmdbId);

    if (cViewType == VIEW_IN_Cn3DCACHE) ValNodeLink(&(bsap->sequences), sep);
    else ValNodeLink(&(basp->sequences), sep);

    for (i=0; i< hits_num; i++) {

      sep = constructSeqEntryForGi(gi[i], TRUE, Database);

      if (sep == NULL) 
         PrtMesC(MAILto,"VASTSRV (Vast2Cn3D)","Unable to get Bioseq.","",false);

      bsp = (BioseqPtr)sep->data.ptrvalue;
      AddMMDBIdToBioseq(bsp, SdiToMmdbId(BsfId[i]/100));      

      if (cViewType == VIEW_IN_Cn3DCACHE) ValNodeLink(&(bsap->sequences), sep);
      else ValNodeLink(&(basp->sequences), sep);
    }
  }

  delete [] gi;
  delete [] BsfId;

  /* Make a linked list of sequence alignments of master and slaves */
  pbsf=pbsfs->features;

  while (pbsf) {
    if (!psaAlignHead) {
      psaAlignHead = fnPBSFtoPSA (pbsf); /* get the sequence alignments */

      if (psaAlignHead == NULL || psaAlignHead->data == NULL) 
	PrtMesC(MAILto, "VASTSERV (Vast2Cn3D)", "Unable to create SeqAnnot.", "", false);

      salpHead = (SeqAlignPtr)(psaAlignHead->data);
      salpTail = salpHead;

      if (JobID!="" && cViewType == VIEW_IN_Cn3D_GEN) {
       	DenseDiagPtr dendiag;

      	for (dendiag=(DenseDiagPtr)salpHead->segs;dendiag;dendiag=dendiag->next)
 	{
      		sid = dendiag->id;
      		sid->choice = 1; 
      		sid->data.ptrvalue = (ObjectIdPtr)objid;
        }
      }
    }
    else {
      psaAlignTail = fnPBSFtoPSA (pbsf);

      salpTail->next = (SeqAlignPtr)(psaAlignTail->data);
      if (psaAlignTail == NULL || psaAlignTail->data == NULL) 
	PrtMesC(NULL, "VASTSRV (Vast2Cn3D)", "Unable to create SeqAnnot.","",false);
      
      salpTail = salpTail->next;
      if (JobID != "" && cViewType == VIEW_IN_Cn3D_GEN) {
	DenseDiagPtr dendiag;

	for (dendiag=(DenseDiagPtr)salpTail->segs;dendiag;dendiag=dendiag->next)
	{
        	sid = dendiag->id;
        	sid->choice = 1;
        	sid->data.ptrvalue = (ObjectIdPtr)objid;
	}
     }
      salpTail->next = NULL;
    }
    pbsf = pbsf->next;
  }

  /* assemble pvnNcbi */
  
  pvnNcbi = ValNodeNew(NULL);
  if (cViewType == VIEW_IN_Cn3D_GEN) {

    pbsaStruct->master = pbsMaster;
    pbsaStruct->slaves = pbsSlaveHead;
    pbsaStruct->alignments = pbsaShort;
    pbsaStruct->seqalign = psaAlignHead;

    pvnNcbi->choice =  NcbiMimeAsn1_alignstruc;
    pvnNcbi->data.ptrvalue = pbsaStruct;
/*    pvnNcbi = (NcbiMimeAsn1Ptr) CheckId(pvnNcbi, JobID); */
  	/* to check identity, yanli  */
  }
  else if (cViewType == VIEW_IN_HTML || cViewType == VIEW_IN_FASTA
	|| cViewType == VIEW_IN_HTML_PAGE || cViewType == VIEW_IN_FASTA_PAGE
        || cViewType == VIEW_IN_TEXT || cViewType == VIEW_IN_TEXT_PAGE) {




    if (JobID != "") { 		// VS job

// first modification, could be done when uploading SeqEntry.
	sep=NULL;
	BioseqPtr   bsp=NULL;
	SeqIdPtr    sip=NULL, sip_next=NULL;
	ObjectIdPtr oip=NULL;
	PDBSeqIdPtr psip=NULL;

        sep=basp->sequences;
	bsp = (BioseqPtr)sep->data.ptrvalue;
	sip = bsp->id;
        if (sip->choice == SEQID_LOCAL) {
	    sip_next = sip->next;
	    oip = (ObjectIdPtr)sip->data.ptrvalue;
	    psip = PDBSeqIdNew();
	    psip->mol= (char *) MemNew (6);
	    StrCut(psip->mol, oip->str, 1, 4);
//strcpy(psip->mol, "Query");
	    psip->chain = (Uint1)(oip->str)[5];
            oip = ObjectIdFree(oip);
	    sip->choice = SEQID_PDB;
	    sip->next = sip_next;
	    sip->data.ptrvalue = psip;
	}
	else 
	    PrtMesC(MAILto, "VastSrv(Vast2Cn3D)", 
		"Can't find local id in JobId = ",(char *)JobID.c_str(), false);
	

/*
// second modification: could be done by changing vpp.aDomName
 	salpHead = (SeqAlignPtr) psaAlignHead->data.ptrvalue;		
	DenseDiagPtr  dendiag;
	for (dendiag = (DenseDiagPtr)salpHead->segs;dendiag;dendiag=dendiag->next)
	{ 
	   sid = dendiag->id;
	   psip = (PDBSeqIdPtr) sip->data.ptrvalue;
	   delete [] psip->mol;
	   psip->mol = new char [6];
	   strcpy(psip->mol, "Query");
	}

*/
    }

    basp->seqalign = psaAlignHead;

    pvnNcbi->choice = NcbiMimeAsn1_alignseq;
    pvnNcbi->data.ptrvalue = basp; 

  }
  else { 	/*if (cViewType == VIEW_IN_Cn3dCache) */

    BiostrucSeqsAlignsCddPtr bsacp;
    SeqAlignData_seq_align_dataPtr sadp;

    bsap->seqaligns = psaAlignHead; 
    bsap->strucaligns = pbsaShort;

    sadp = ValNodeNew(NULL);
    sadp->data.ptrvalue = bsap;
    sadp->choice = SeqAlignData_seq_align_data_bundle;
    bsacp = BiostrucSeqsAlignsCddNew();
    bsacp->SeqAlignData_seq_align_data = sadp;
    if (complexity == ONECOORDRES)     
		bsacp->structure_type = 2;	/* ncbi_backbone */
    else bsacp->structure_type = 3;		/* ncbi_all_atoms */

    pvnNcbi->choice = NcbiMimeAsn1_general;
    pvnNcbi->data.ptrvalue = bsacp;

  } 

{AsnIoPtr aipr;
aipr=AsnIoOpen("pvnNcbi.out", "w");
NcbiMimeAsn1AsnWrite(pvnNcbi, aipr, NULL);
AsnIoClose(aipr);
}

  if (cViewType == VIEW_IN_Cn3D_GEN || cViewType == VIEW_IN_Cn3DCACHE) {
    if (iPDB == 0)    /* cn3d MIME */
      printf ("Content-type: chemical/ncbi-asn1-binary\n\n");

    else if (iPDB == 1) {    /* "See File" */
      printf ("Content-type: text/html\n\n");
      printf ("<HTML><body><pre>\n");
    }
    else        /* "Save File" */
      printf ("Content-type: application/octet-stream\n\n");

    if (iPDB != 1) paiFile = AsnIoNew(ASNIO_BIN_OUT, stdout, NULL, NULL, NULL); 
    else paiFile = AsnIoNew(ASNIO_TEXT_OUT, stdout, NULL, NULL, NULL);

    NcbiMimeAsn1AsnWrite(pvnNcbi, paiFile, NULL); 
    AsnIoFlush(paiFile); 
    AsnIoClose(paiFile);
  }
  else {

    CNcbi_mime_asn1 *ptr = new CNcbi_mime_asn1;
    CAsnConverter < CNcbi_mime_asn1, NcbiMimeAsn1>
                CtoCpp((AsnWriteFunc)NcbiMimeAsn1AsnWrite,
                                        (AsnReadFunc)NcbiMimeAsn1AsnRead);

    if ( CtoCpp.FromC(pvnNcbi, ptr) ) {

    	printf("Content-type: text/html\n\n");
    	WWWPrintFileData((char *)"sshead.txt", stdout);
    	printf("<TABLE width=800 BORDER=0 CELLPADDING=3 CELLSPACING=0 bgcolor=#FFFFCC>\n\n");
    	PrintQueryInfo(stdout, pdbname_m, chain_m, domNo_m);
    	printf("</TABLE>\n");
    	printf("<br>\n");
    	if (cViewType == VIEW_IN_HTML || cViewType == VIEW_IN_HTML_PAGE)
      	    CAV_DisplayMultiple(*ptr, CAV_HTML|CAV_SHOW_IDENTITY, 60, 2.0, 
						NULL, 0, NULL, NULL, NULL); 
        else if (cViewType == VIEW_IN_FASTA || cViewType == VIEW_IN_FASTA_PAGE) 
	{
      	    printf("<pre>\n");
            CAV_DisplayMultiple(*ptr, 
		CAV_FASTA|CAV_LEFTTAILS|CAV_RIGHTTAILS|CAV_FASTA_LOWERCASE, 
		60, 2.0, NULL, 0, NULL, NULL, NULL);
      	    printf("</pre>\n");
     	}
    	else if (cViewType == VIEW_IN_TEXT || cViewType == VIEW_IN_TEXT_PAGE) {
		printf("<pre>\n");
		CAV_DisplayMultiple(*ptr, CAV_TEXT, 60, 2.0, NULL, 0, NULL,NULL, NULL);
		printf("</pre>\n");
    	}
    }
    else {
	printf("Content-type: text/html\n\n");
        WWWPrintFileData((char *)"sshead.txt", stdout);
        printf("CtoCpp.FromC failed\n");

    }

    printf("<HR SIZE=5 NOSHADE>\n");
    printf("</body></html>\n");

    delete ptr;		// forgot 6/23/65
  }

  if (cViewType==VIEW_IN_Cn3D_GEN) {
  	CloseMMDBAPI();
  	MMDBFini();
  }

  VastSrchFinish();

  NcbiMimeAsn1Free(pvnNcbi); 

} 	/* end of VastToCn3DAndAli */
