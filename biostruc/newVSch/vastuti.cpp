/*
* $Id: vastuti.cpp,v 1.1 2005/07/26 17:12:50 chenj Exp $
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
* $Log: vastuti.cpp,v $
* Revision 1.1  2005/07/26 17:12:50  chenj
* Making linux VSNbr.cgi
*
* Revision 1.1.1.1  2002/12/06 20:17:21  chenj
* Imported Scouces
*
*
*
*
* This is the utility program file. 
*
* ==========================================================================
*/




#include "hUtilib.hpp"
#include "vastuti.hpp"


extern Uint1 	numSubsets;
extern SubsetNameData  subsetNames[5];
extern string	MAILto;
extern char	DATApath[MAX_TBUFF], VSPATH[MAX_TBUFF];
extern char	URLcgi[MAX_TBUFF], CGIname[MAX_TBUFF];
extern bool	SubsetInfoLoaded;
extern string  	JobID, Passwd;


void PrtMesC(string MAILto, CharPtr str1, CharPtr str2, CharPtr str3, bool ret)
{
        printf("Content-type: text/html\n\n");
        printf("<body bgcolor=#f0f0f0>\n");
        printf("<br><H2>%s</H2><p>\n", str1);
        if (str2) {
                printf("<H2>%s", str2);
                if (str3) printf("\"%s\". </H2><p>\n", str3);
                else printf("</H2><p>\n");
        }
        if (MAILto != "")
                printf("<H2>Please alert \"%s\" of this problem. </H2>\n",
                                                                        MAILto);
        if (!ret) exit(1);
}


void loadSubsetInfo(void)
{

  numSubsets = 0;
  numSubsets = constructSubsetNames(subsetNames, 10);

  if(numSubsets <= 0)
        PrtMesC(MAILto, "VASTSRV", "loadSubsetInfo() failed.", "", false);

  SubsetInfoLoaded = true;

} /* loadSubsetInfo() */



string GetSubsetName_DB(unsigned subsetId)
{
  Uint1 i = 0;
  if(SubsetInfoLoaded == false) loadSubsetInfo();
  for(i = 0; i < numSubsets; i++)
    if(subsetNames[i].subId == (Uint1) subsetId)
      return (string)subsetNames[i].subName;

  PrtMesC(MAILto, "VASTSRV", "GetSubsetName_DB() failed.", "", false);
  return ("");

} /* GetSubsetName_DB() */





/************************** DZ: extracted from mmdbsrv.c **********************
 * WWWPrintFileData looks in the current CGI-BIN directory, or the "data"
 * subdirectory for the data file and prints it out to pFile.
 */
 
void WWWPrintFileData(char *FName,  FILE *pFile)
{
   FILE *f = NULL;
   char fullpath [PATH_MAX], pcBuf[1024];
   char *ptr;  
   
   fullpath[0] = '\0';
   StringCpy(fullpath,  DATApath); /* look in DATApath */
   StringCat(fullpath,  FName);
   f = FileOpen (fullpath, "r");
   if (f == NULL) {
       f = FileOpen (FName, "r");  /* look in curent */
       if (f == NULL)  {  /* look in ./data/ */
         ProgramPath (fullpath, sizeof (fullpath) - 1);
         ptr = StringRChr (fullpath, DIRDELIMCHR);
         if (ptr != NULL) {
	  *ptr = '\0';
         }
         FileBuildPath (fullpath, (char *)"data", FName);  
         f = FileOpen (fullpath, (char *)"r");
         if (f == NULL)  {
           return;  
         } 
       }
   }
   do {
     pcBuf[0] = '\0';
     ptr = fgets(pcBuf, (size_t)1024, f);
     if (ptr) fprintf(pFile, ptr);
   } while (ptr);
   FileClose(f);
} /* WWWPrintFileData() */



/* Note: A lot of this common html text could be put into a header file,
 * read in, and then spewed out.
 * Cf. the way this is done in mmdbsrv with WWWPrintFile.
 * Sometime I'll take the time to write a more general version of the latter
 * function that can be used for this purpose.  T.M.
 */

/* DZ: most of the HTML file URLs are inside the new VAST subdirectory, and
 * references to these are interwoven with many query-specific references.
 * Becasue of this interleaving, using WWWPrintFileData would require several
 * invocations and corresponding text files. Instead, WWWPrintFileData is only
 * used to define the new path to CN3D, and all vast html references are now
 * prefixed with the VASTpath string. This variable has been added to .vastrc
 * and to the list of parameters to be extracted using GetVastParams.
 */




BiostrucAnnotSetPtr LIBCALL
LocalGetBiostrucAnnotSet(unsigned mmdbid)
{
  AsnIoPtr aip = NULL;
  BiostrucAnnotSetPtr pbsa = NULL;
  char path[PATH_MAX];
  char pcId[20];
  unsigned iFileExists = 0;

  sprintf(pcId, "/%ld", (long) mmdbid);
  StringCpy(path, VSPATH);
  StringCat(path, JobID.c_str());
  StringCat(path, pcId);
  StringCat(path, ".bas");
  iFileExists = FileLength(path);
  if (iFileExists == 0)
    return NULL;

  aip = AsnIoOpen(path, (char *)"r");
  pbsa = BiostrucAnnotSetAsnRead(aip, NULL);
  AsnIoClose (aip);
  if (!pbsa) return NULL;
  return pbsa;

} /* LocalGetBiostrucAnnotSet() */




/*
BiostrucAnnotSetPtr LIBCALL
LocalGetFeatureSet(unsigned mmdbid, unsigned feature_set_id)
{
    BiostrucAnnotSetPtr basp2 = NULL;
    BiostrucFeatureSetPtr pbsfs = NULL;
    BiostrucAnnotSetPtr basp = NULL;
    BiostrucFeatureSetPtr pbsfsLast = NULL;

    if (IsVASTData(mmdbid))
       basp = VASTBsAnnotSetGet(mmdbid);
    else if (IsVASTData(feature_set_id)) {
       basp = VASTBsAnnotSetGet(feature_set_id);
       if (basp != NULL) return basp;
    }
    else
    if (JobID != "")
       basp = LocalGetBiostrucAnnotSet(mmdbid);

    if (basp == NULL)
        return NULL;

    pbsfs = basp->features;
    pbsfsLast  = NULL;
    basp2 = NULL;
    while (pbsfs)
     {
       if (pbsfs->id == feature_set_id)
        {
          basp2 = BiostrucAnnotSetNew();
          basp2->id = basp->id;
          basp->id = NULL; // unlink the id valnode from basp object
          basp2->descr = basp->descr;
          basp->descr = NULL;  // unlink the descr from basp object
          basp2->features = pbsfs;
          if (pbsfsLast) // relink next to prev 
            pbsfsLast->next = pbsfs->next;
          else
            basp->features = pbsfs->next;
          basp2->features->next = NULL;
          BiostrucAnnotSetFree(basp);
          return basp2;
        }
       pbsfsLast = pbsfs;
       pbsfs = pbsfs->next;
     }  
    BiostrucAnnotSetFree(basp);
    return basp2;

}	// LocalGetFeatureSet */



BiostrucAnnotSetPtr LIBCALL BiostrucAnnotSetGetByFid (
    BiostrucAnnotSetPtr basp, int feature_id, int feature_set_id)
{
  BiostrucAnnotSetPtr basp2 = NULL;
  BiostrucFeatureSetPtr pbsfs = NULL;
  BiostrucFeaturePtr pbsf = NULL;

  if (basp == NULL)
    return NULL;

  pbsfs = basp->features;
  while (pbsfs) {
    if (pbsfs->id == feature_set_id) {
      pbsf =  pbsfs->features;
      while(pbsf) {

        if (pbsf->id == feature_id) {  /* found it */
          basp2 = BiostrucAnnotSetNew();
          basp2->id = basp->id;
          basp->id = NULL; /* unlink the id valnode from basp object */
          basp2->descr = basp->descr;
          basp->descr = NULL;  /* unlink the descr from basp object */
          basp2->features = BiostrucFeatureSetNew();
          basp2->features->id = pbsfs->id;
          basp2->features->descr = pbsfs->descr;
          pbsfs->descr = NULL; /* unlink */
          basp2->features->features = BiostrucFeatureNew();
          basp2->features->features->id = pbsf->id;
          basp2->features->features->name = StringSave(pbsf->name);
          basp2->features->features->type = pbsf->type;
          basp2->features->features->Property_property=pbsf->Property_property;
          pbsf->Property_property = NULL; /* unlink */
          basp2->features->features->Location_location=pbsf->Location_location;
          pbsf->Location_location = NULL; /* unlink */
          BiostrucAnnotSetFree(basp);
          return basp2;
        }
        pbsf = pbsf->next;
      }
    }
    pbsfs = pbsfs->next;
  }
  BiostrucAnnotSetFree(basp);
  return basp2;
} /* BiostrucAnnotSetGetByFid() */



short LIBCALL Check_VastSearch_Password()
{
 
  char pcPassFile[24];
  char PassPath[PATH_MAX];
  FILE *passwdfile;
  char * pcPassOld;
  unsigned iPassLen;

  iPassLen = StringLen(Passwd.c_str());
  pcPassOld = StringSave(Passwd.c_str());
   
  sprintf(pcPassFile, "/%s.passwd", JobID.c_str());
  PassPath[0]='\0';
  StringCpy(PassPath, VSPATH);
  StringCat(PassPath, JobID.c_str());
  StringCat(PassPath, pcPassFile); 
  
  if ((passwdfile = FileOpen(PassPath, "r")) == NULL) {
        char str[100];

	sprintf(str, (char *)"%s", JobID.c_str());
      	PrtMesC(NULL, (char *)"VASTSRV", 
		(char *)"Password file doesn't exist for JobId = ", str, true);
	return 2;
  }
  
  fscanf(passwdfile, "%s", pcPassOld);
  FileClose(passwdfile);

  if (!StringNCmp(pcPassOld, Passwd.c_str(), iPassLen))
    return 1;
  else 
    return 0;
} /* end of Check_VastSearch_Password */



/* new Prune */
BiostrucAnnotSetPtr  PruneBiostrucAnnotHits(
    BiostrucAnnotSetPtr basp, unsigned FSID, ValNodePtr pvnFids)
{
  BiostrucAnnotSetPtr basp2 = NULL;
  BiostrucFeatureSetPtr pbsfs = NULL;
  BiostrucFeaturePtr pbsf = NULL, pre_pbsf = NULL;
  BiostrucFeaturePtr pbsfHold = NULL;
  BiostrucFeaturePtr pbsf2 = NULL;

  bool found = false;

  if ((basp == NULL) || (pvnFids == NULL) || (FSID == 0))
    return NULL;

  pbsfs = basp->features;
  while (pbsfs) {
    if (pbsfs->id == FSID) {
      basp2 = BiostrucAnnotSetNew();
      basp2->id = basp->id;
      basp->id = NULL; /* unlink the id valnode from basp object */
      basp2->descr = basp->descr;
      basp->descr = NULL;  /* unlink the descr from basp object */
      basp2->features = BiostrucFeatureSetNew();
      basp2->features->id = pbsfs->id;
      pbsfs->id = NULL;
      basp2->features->descr = pbsfs->descr;
      pbsfs->descr = NULL; /* unlink */
      pbsfHold = pbsfs->features;
      pbsfs->features = NULL;
      BiostrucAnnotSetFree(basp);
      pbsfs = NULL; 
      while (pvnFids) {

	found = false;
	pbsf = pbsfHold;
	pre_pbsf = NULL;

	while (pbsf) {
          if ((unsigned) pvnFids->data.intvalue == pbsf->id) {
            found = true;
            break;
          }
	  pre_pbsf = pbsf;
          pbsf = pbsf->next;
        }

        if (found) {
	  if (!pre_pbsf) pbsfHold = pbsf->next;
	  else pre_pbsf->next = pbsf->next;

          if (!basp2->features->features) {
            basp2->features->features = pbsf;
            pbsf2 = basp2->features->features;
          }
          else {
            pbsf2->next = pbsf;
            pbsf2 = pbsf;
          }
        }

	pvnFids = pvnFids->next;

      }  /* while pvbFids */

      pbsf2->next = NULL;
    }
    if(pbsfs) pbsfs = pbsfs->next;
  }
  return basp2;

} /* PruneBiostrucAnnotHits() */




void MakeNbrList(char * nbrString,  char * *(*SelNames), unsigned **SelSds,
	unsigned *iCount, unsigned ischar)
{

  char		SavedChar = NULLB;
  char  	*Ptr=NULL, *Ptr2=NULL;
  unsigned	i, j;

  Ptr2 = Ptr = SkipSpaces(nbrString);
  if (*Ptr == NULLB) return;

  *iCount = 0;

  while (*Ptr) {
	if (*Ptr != ',') (*iCount)++;
	Ptr = SkipToSet(Ptr, (char *)" ,");
	Ptr = SkipSet(Ptr, (char *)" ,");	
  }

  if (!(*iCount)) return;
	
  if (ischar) {
      *SelNames = new char* [*iCount];
      for (i=0; i< *iCount; i++) (*SelNames)[i] = new char [100];
  }
  else *SelSds = new unsigned [*iCount];

  Ptr = Ptr2;
  i = 0;
  while (*Ptr) {
	Ptr2 = SkipToSet(Ptr, (char *)" ,");
	SavedChar = *Ptr2;
	*Ptr2 = NULLB;
	if (ischar && Ptr[0] != '\0') {
	    for (j=0; j< StrLen(Ptr); j++) Ptr[j] = toupper(Ptr[j]);
	    if (Ptr[4] == '.') Ptr[4] = ' '; 
	    strcpy((*SelNames)[i++], Ptr);
	}
	else if (!ischar) {

	         if ( isInt(Ptr) ) (*SelSds)[i++] = atoi(Ptr);
		 else {
			
		      PrtMesC("", "VASTSRV", "Please use numbers as 3D-Domain ids, but not this: ", Ptr, false); 
		 }
	
	}

	*Ptr2 = SavedChar;
	Ptr = SkipSet(Ptr2, (char *)" ,");
  }

  *iCount = i;

} 	/* end of MakeNbrList */




void 
OrderCopyVpp(VastPageDataPtr vpptmp, VastPageDataPtr vpp, unsigned iKept, unsigned * KepId)
{
    unsigned 	i, j, k;

    for (i=0; i< iKept; i++) {
	j=0; 
	while (j<iKept &&
		KepId[i] != getfid(vpptmp[j].bSdi, vpptmp[j].bAlignId)) j++;
	vpp[i].id = vpptmp[j].id;
	vpp[i].acc = vpptmp[j].acc;
	vpp[i].aSdi = vpptmp[j].aSdi;
	vpp[i].bSdi = vpptmp[j].bSdi;
	vpp[i].bAlignId = vpptmp[j].bAlignId;
	sprintf(vpp[i].bDomName, vpptmp[j].bDomName);
	vpp[i].vScore = vpptmp[j].vScore;
	vpp[i].mlogp = vpptmp[j].mlogp;
	vpp[i].rmsd = vpptmp[j].rmsd;
	vpp[i].nres = vpptmp[j].nres;
	vpp[i].pcntId = vpptmp[j].pcntId;
	vpp[i].t11 = vpptmp[j].t11;
	vpp[i].t12 = vpptmp[j].t12;
	vpp[i].t13 = vpptmp[j].t13;
	vpp[i].r11 = vpptmp[j].r11;
	vpp[i].r12 = vpptmp[j].r12;
	vpp[i].r13 = vpptmp[j].r13;
	vpp[i].r21 = vpptmp[j].r21;
	vpp[i].r22 = vpptmp[j].r22;
	vpp[i].r23 = vpptmp[j].r23;
	vpp[i].r31 = vpptmp[j].r31;
	vpp[i].r32 = vpptmp[j].r32;
	vpp[i].r33 = vpptmp[j].r33;
	vpp[i].t21 = vpptmp[j].t21;
	vpp[i].t22 = vpptmp[j].t22;
	vpp[i].t23 = vpptmp[j].t23;
	sprintf(vpp[i].bDescr1, vpptmp[j].bDescr1);
	sprintf(vpp[i].bDescr2, vpptmp[j].bDescr2);
	sprintf(vpp[i].aDomName, vpptmp[j].aDomName);
	for (k=0; k< vpptmp[j].IpdpLen; k++) {
	     vpp[i].Ipdp[k].acc     = vpptmp[j].Ipdp[k].acc;
	     vpp[i].Ipdp[k].mstfrom = vpptmp[j].Ipdp[k].mstfrom;
	     vpp[i].Ipdp[k].mstto   = vpptmp[j].Ipdp[k].mstto;
	     vpp[i].Ipdp[k].slvfrom = vpptmp[j].Ipdp[k].slvfrom;
	     vpp[i].Ipdp[k].slvto   = vpptmp[j].Ipdp[k].slvto;
	}
	vpp[i].IpdpLen = vpptmp[j].IpdpLen;
	vpp[i].numHits = vpptmp[j].numHits;
    }

} /* OrderVpp */
