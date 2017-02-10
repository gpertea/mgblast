/*  lsqfetch.c
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
* File Name:  lsqfetch.c
*
* Author:  Jinghui Zhang
*
* Version Creation Date: 5/24/95
*
*
* File Description:  functions for fetching the local sequences
*
* Modifications:
* --------------------------------------------------------------------------
* Date     Name        Description of modification
*
* $Log: lsqfetch.c,v $
* Revision 6.30  2004/10/27 20:07:14  kans
* LsqFetch_AsnIoOpen to suppress missing file warning, similar to LsqFetch_FileOpen
*
* Revision 6.29  2004/10/26 14:45:30  kans
* LsqFetch_FileOpen suppresses FileOpen failure INFO message
*
* Revision 6.28  2004/10/05 19:11:23  kans
* separate internal CreateBinaryAsnIndex and CreateTextAsnIndex functions
*
* Revision 6.27  2004/10/05 18:56:48  kans
* AsnIndexedLibFetchEnable only works in text mode, backed out SeqEntryAsnRead change - will later implement use of catenated Seq-entry instead of se2bss processed input if text
*
* Revision 6.26  2004/10/05 17:34:16  kans
* protect all binary search functions against R out of range
*
* Revision 6.25  2004/10/05 17:25:04  kans
* AsnIndexedLibBioseqFetchFunc handles all Seq-id types, passes atp_se to SeqEntryAsnRead
*
* Revision 6.24  2004/10/05 16:21:37  kans
* LocalSeqFetchInit checks for INDEXED_TEXT_ASN and INDEXED_BIN_ASN, calls AsnIndexedLibFetchEnable
*
* Revision 6.23  2004/08/04 20:21:04  kans
* record alfp->binary at correct place in asn indexed fetch enable function
*
* Revision 6.22  2004/08/03 17:51:49  kans
* added AsnIndexedLibFetch enable and disable functions
*
* Revision 6.21  2004/08/02 19:10:14  kans
* added CreateAsnIndex for indexing Bioseq-set ftp release files
*
* Revision 6.20  2004/04/13 16:58:32  kans
* allow alt index to have gi numbers, also test .fsa if .fa fails
*
* Revision 6.19  2003/11/13 17:18:02  kans
* added SearchAltIndex, finished Alt fetch for chimp revision
*
* Revision 6.18  2003/11/12 23:49:11  kans
* SortIfpByID needed LIBCALLBACK for PC
*
* Revision 6.17  2003/11/12 23:38:48  kans
* changing AltIndexedFastaLibFetchEnable prototype, implementation not yet finished
*
* Revision 6.16  2003/08/27 21:24:05  kans
* enable alt indexed fasta looks up previously registered function, changes settings for new path
*
* Revision 6.15  2003/08/27 19:27:43  kans
* added AltIndexedFastaLibFetch functions for chimpanzee genome project
*
* Revision 6.14  2002/11/13 23:07:37  johnson
* Changed make_lib such that it looks to see if it matches the *whole* seq-id
* (defined by the next character being non-alphanumeric).
*
* Revision 6.13  2002/07/19 20:16:33  johnson
* bug fix in make_lib -- wasn't properly handling sequences >=1000 residues
*
* Revision 6.12  2002/01/22 19:50:32  kans
* IndexedFastaLibBioseqFetchFunc looks for prefix with upper case followed by lower case
*
* Revision 6.11  2001/09/21 20:02:04  kans
* allow U and u to be DNA in CheckDnaResidue for the BLAST guys, even though it is amino acid Selenocysteine in most places
*
* Revision 6.10  2001/03/13 16:48:58  kans
* fixes to saving path, binary search results
*
* Revision 6.9  2001/03/12 23:19:33  kans
* added IndexedFastaLib functions - currently uses genome contig naming conventions
*
* Revision 6.8  2001/01/09 00:12:39  kans
* now handles SEQID_GI
*
* Revision 6.7  1999/10/07 16:21:13  kans
* removed static AddSeqId and SeqIdDupList which were identical to public sequtil functions
*
* Revision 6.6  1999/07/20 21:18:39  sicotte
* add static AddSeqId for linker conflicts
*
* Revision 6.5  1999/07/20 21:16:49  sicotte
* add static SeqIdDupList for linker conflicts
*
* Revision 6.4  1999/04/01 22:26:13  sicotte
* Make lsqfetch Attempt to Parse the fasta defline, otherwise use the supplied SeqId
*
* Revision 6.3  1999/03/11 23:39:33  kans
* sprintf and sscanf casts
*
* Revision 6.2  1998/02/06 17:41:33  zjing
* make the function CheckDnaResidue external
*
* Revision 5.3  1997/07/18 21:31:36  zjing
* move some def to lsqfetch.h
*
* Revision 5.2  1997/06/19 18:38:16  vakatov
* [WIN32,MSVC++]  Adopted for the "NCBIOBJ.LIB" DLL'ization
*
* Revision 5.1  1997/02/12 22:42:43  zjing
* fix a memory leak
*
 * Revision 5.0  1996/05/28  13:23:23  ostell
 * Set to revision 5.0
 *
 * Revision 4.9  1996/03/06  18:28:53  zjing
 * .
 *
 * Revision 4.4  1995/10/11  19:29:28  zjing
 * add LIBCALL for find_big_bioseq
 *
 * Revision 4.3  1995/08/23  12:37:08  epstein
 * remove leading white space from conditional compilation
 *
 * Revision 4.2  1995/08/04  17:31:32  kans
 * JZ added LocalSeqFetchInit, LocalSeqFetchDisable
 *
 * Revision 4.1  1995/08/03  20:56:18  kans
 * paths can now be specified in a regular NCBI config file
 *
 * Revision 4.0  1995/07/26  13:49:01  ostell
 * force revision to 4.0
 *
 * Revision 1.4  1995/07/10  14:20:31  zjing
 * check in default search path for fasta files
 *
*
*
* ==========================================================================
*/


#ifndef _LSQFETCH_
#include <lsqfetch.h>
#endif

#include <sqnutils.h>

/* #include <accentr.h> */

/********************************************************************
*
*	names of config files in different platforms
*
*********************************************************************/
#ifdef WIN_MSWIN
	static CharPtr seqinfo_file = "seqinfo.dat";
#else
	static CharPtr seqinfo_file = ".seqinfo";
#endif

/********************************************************************
*
*	Local versions of FileOpen and AsnIoOpen suppress missing file error report
*
*********************************************************************/

static FILE* LIBCALL LsqFetch_FileOpen (const char *filename, const char *mode)

{
  FILE    *fp;
  ErrSev  sev;

  sev = ErrSetMessageLevel (SEV_ERROR);
  fp = FileOpen (filename, mode);
  ErrSetMessageLevel (sev);
  return fp;
}

static AsnIoPtr LIBCALL  LsqFetch_AsnIoOpen (CharPtr file_name, CharPtr mode)

{
  AsnIoPtr  aip;
  ErrSev    sev;

  sev = ErrSetMessageLevel (SEV_ERROR);
  aip = AsnIoOpen (file_name, mode);
  ErrSetMessageLevel (sev);
  return aip;
}


/***********************************************************************
***
*
*	fasta_sep(): making a Seq-entry from a FASTA formatted file
*
************************************************************************
***/
/***********************************************************************
*
*	Check if the sequence is a DNA or protein 
*	ck_len: the length for checking
*	pnon_DNA: store the number of non-DNA residue
*	return TRUE if it is a DNA sequence, FALSE for protein
*
***********************************************************************/
NLM_EXTERN Boolean CheckDnaResidue(CharPtr seq_ptr, Int4 ck_len, Int4Ptr pnon_DNA)
{
 
        Int4 i, non_DNA=0;
	Int4 len;
 
        for(i=0, len = 0; i<ck_len; ++i)
	{
		if(IS_ALPHA(seq_ptr[i]))
		{
             		if(StrChr("ACGTUNacgtun", seq_ptr[i])==NULL)
               			++non_DNA;
			++len;
		}
	}
 
	if(pnon_DNA != NULL)
		*pnon_DNA = non_DNA;
        if(non_DNA >= len/4)
          return FALSE;
        else return TRUE;
}


static SeqEntryPtr Sep_from_ByteStore(ByteStorePtr bsp, Int4 length, Boolean is_dna, SeqIdPtr sip) 
{
  	BioseqPtr biosp;
  	SeqEntryPtr sep;

	biosp = BioseqNew();
	biosp->id = SeqIdDupList(sip);
        biosp->seq_data = bsp;
        biosp->repr = Seq_repr_raw;
        biosp->mol = (is_dna) ? Seq_mol_dna : Seq_mol_aa;
        biosp->length = length;
        biosp->topology = 1;                    /** linear sequence**/
        biosp->seq_data_type = (is_dna) ? Seq_code_iupacna : Seq_code_ncbieaa;
	if(is_dna)
		BioseqRawConvert(biosp, Seq_code_ncbi4na);
        sep = SeqEntryNew();
        sep->choice = 1;
        sep->data.ptrvalue = biosp;
	SeqMgrSeqEntry(SM_BIOSEQ, (Pointer)biosp, sep);

        return sep;

}

static ByteStorePtr make_lib(FILE *ifp, CharPtr name, Int4Ptr length, BoolPtr is_DNA,SeqIdPtr PNTR sip)
{
   Int4 n_len, pos, seq_len;
   Char temp[1001];
   CharPtr p;
   ByteStorePtr bsp;
   Boolean is_found, is_end;
   Boolean check_DNA = FALSE;
   Int4 i;

   if(sip)
       *sip=NULL;
	if(name != NULL)
	{
		rewind(ifp);
   		n_len = StringLen(name);
	}
   	is_found = FALSE;
   	while(!is_found && FileGets(temp, 1000, ifp) != NULL)   /*find the right seq*/
   	{
		if(temp[0] == '>')
		{
			if(name!=NULL)
 			{
				i = 1;
				while(IS_WHITESP(temp[i]))
					++i;
				if(StringNCmp(temp+i, name, (size_t)n_len) ==0 && !isalnum(temp[i+n_len])) {
                                    is_found = TRUE;
                                    if(sip)
                                        *sip = SeqIdParse(temp+i);
				} else
					is_found = FALSE;
			}
			else {
                            i=1;
                            while(IS_WHITESP(temp[i]))
                                    ++i;
                            if(sip) {
                                CharPtr blank;
                                blank = StringChr(temp+i,' ');
                                if(blank)
                                    *blank='\0';
                                *sip = SeqIdParse(temp+i);
                                if(blank)
                                    *blank=' ';
                                
                            }                            
                            is_found = TRUE;
                        }
		}
		if(is_found)
		{
			pos = ftell (ifp);
			is_end = FALSE;
			seq_len =0;
			while(FileGets(temp, 1000, ifp) != NULL && !is_end)
			{
				if(temp[0] == '>')
					is_end = TRUE;
				else
					seq_len += StringLen(temp);
			}
                        if (seq_len > 0) /* don't count \n */
                            --seq_len;
		}
	}


    if(!is_found)
	return NULL;

    bsp = BSNew(seq_len+1);
    BSSeek(bsp, 0, SEEK_SET);

    fseek(ifp, pos, SEEK_SET);
    is_end = FALSE;
    seq_len = 0;
    while(FileGets(temp, 1000, ifp) != NULL && !is_end)
    {
	if(temp[0] == '>')
		is_end = TRUE;
	else
	{
	    if(!check_DNA)	/*check if it is a DNA sequence*/
	    {
		*is_DNA = CheckDnaResidue(temp, StringLen(temp), NULL);
		check_DNA = TRUE;
	    }
	    for(p=temp; *p!='\0' && *p != '\n'; ++p)
		if(IS_ALPHA(*p) || *p == '-' || *p == '*')
		{
			++seq_len;
			if(*p == '-')
				*p = '@';
			else
				*p = TO_UPPER(*p);
			BSPutByte(bsp, (Int2)(*p));
		}
        }
    }

    *length = seq_len;
    return bsp;
}


/*************************************************************************
***
*	fasta_lib_sep(): make a Seq-entry from the FASTA library
*
**************************************************************************
***/
NLM_EXTERN SeqEntryPtr fasta_lib_sep PROTO((FILE *fp, CharPtr seq_name, SeqIdPtr sip));

NLM_EXTERN SeqEntryPtr fasta_lib_sep(FILE *fp, CharPtr seq_name, SeqIdPtr sip)
{
	
	ByteStorePtr bsp;
	Int4 length;
	Boolean is_dna;
        SeqIdPtr sipnew;
	if((bsp = make_lib(fp, seq_name, &length, &is_dna,&sipnew)) != NULL) {
            if(sipnew) {
                sipnew=ValNodeLink(&sipnew,sip);
                sip=sipnew;
            }
            return Sep_from_ByteStore(bsp, length, is_dna, sip);
	} else
		return NULL;
}

	
static void FindBigCallback(SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
	ValNodePtr vnp;
	BioseqPtr bsp, curr;
	BioseqSetPtr bssp;

        if(data == NULL)
           return;
        vnp = (ValNodePtr)data;
        if(sep->choice != 1)
	{
		bssp = sep->data.ptrvalue;
		if(bssp->_class == 1)
			vnp->choice = 1;
		return;
	}

        bsp = sep->data.ptrvalue;
	if(bsp)
	{
		if(vnp->choice == 1)	/*needs the DNA sequence*/
		{
			if(bsp->mol != Seq_mol_aa)
			{
				if(vnp->data.ptrvalue == NULL)
					vnp->data.ptrvalue = bsp;
				else
				{
					curr = (BioseqPtr)(vnp->data.ptrvalue);
					if(bsp->length > curr->length)
						vnp->data.ptrvalue = bsp;
				}
				return;
			}
		}
		if(vnp->data.ptrvalue == NULL)
			vnp->data.ptrvalue = bsp;
		else
		{
			curr = (BioseqPtr)(vnp->data.ptrvalue);
			if(bsp->length > curr->length)
				vnp->data.ptrvalue = bsp;
		}
	}
 
}
 
 
NLM_EXTERN BioseqPtr LIBCALL find_big_bioseq(SeqEntryPtr sep)
{
	ValNode vn;

	vn.choice = 0;
	vn.data.ptrvalue = NULL;
	if (sep != NULL)
    		SeqEntryExplore (sep, (Pointer) (&vn), FindBigCallback);
 
  	return (BioseqPtr)(vn.data.ptrvalue);
}


 
/**********************************************************************
*
*	FastaLibBioseqFetchEnable(libs, now)
*	Initiate the function for fetch a Bioseq from a Fasta Library 
*	file. libs is a list of library file names. 
*	If now = TRUE, open the library files and set the state to 
*	FASTALIB_OPEN. return TRUE for success. 
*
***********************************************************************/
static CharPtr libproc = "FastaLibBioseqFetch";

static Pointer LIBCALLBACK FreeSeqName(Pointer data)
{
	MemFree(data);
	return NULL;
}
static Int2 LIBCALLBACK FastaLibBioseqFetchFunc (Pointer data)
{
	OMProcControlPtr ompcp;
	ObjMgrProcPtr ompp;
	FastaLibPtr flp;
	SeqIdPtr sip;
	OMUserDataPtr omdp;
	CharPtr seq_name = NULL;
	Char name[100];
	SeqEntryPtr sep = NULL;
	BioseqPtr bsp;

	ompcp = (OMProcControlPtr)data;
	ompp = ompcp->proc;
	flp = (FastaLibPtr)(ompp->procdata);
	sip = (SeqIdPtr)(ompcp->input_data);
	if(sip == NULL)
		return OM_MSG_RET_ERROR;
	/* if(sip->choice == SEQID_GI)
		return OM_MSG_RET_OK; */

	if(ompcp->input_entityID)
	{
		omdp = ObjMgrGetUserData(ompcp->input_entityID, ompp->procid, OMPROC_FETCH, 0);
		if(omdp != NULL)
			seq_name = omdp->userdata.ptrvalue;
	}

	if(seq_name == NULL)
	{
		seqid_to_string(sip, name, flp->use_locus);
		seq_name = name;
	}

	while(flp && sep==NULL)
	{
		if(flp->state == FASTALIB_CLOSE)
		{
			flp->fp = LsqFetch_FileOpen(flp->file_name, "r");
			if(flp->fp == NULL)
				flp->state = FASTALIB_ERROR;
			else
				flp->state = FASTALIB_OPEN;
		}
		if(flp->state == FASTALIB_OPEN)
			sep = fasta_lib_sep(flp->fp, seq_name, sip);
		if(sep == NULL)
			flp = flp->next;
	}
	
	if(sep == NULL)
		return OM_MSG_RET_OK;

	bsp = BioseqFindInSeqEntry(sip, sep);
	if(bsp == NULL)
		bsp = find_big_bioseq(sep);
	ompcp->output_data = (Pointer)bsp;
	ompcp->output_entityID = ObjMgrGetEntityIDForChoice(sep);
	omdp = ObjMgrAddUserData(ompcp->output_entityID, ompp->procid, OMPROC_FETCH, 0);
	omdp->userdata.ptrvalue = StringSave(seq_name);
	omdp->freefunc = FreeSeqName;
	return OM_MSG_RET_DONE;
}

NLM_EXTERN Boolean FastaLibBioseqFetchEnable(ValNodePtr libs, Boolean now)
{
	FastaLibPtr flp = NULL, new, curr;
	Boolean ok;
	FILE *fp;
	CharPtr file_name;
	
	while(libs)
	{
		ok = TRUE;
		file_name = libs->data.ptrvalue;
		if(now)
		{
			if((fp = LsqFetch_FileOpen(file_name, "r")) == NULL)
				ok = FALSE;
		}
		if(ok)
		{
			new = MemNew(sizeof(FastaLib));
			new->use_locus = FALSE;
			StringCpy(new->file_name, file_name);
			if(now)
			{
				new->state = FASTALIB_OPEN;
				new->fp = fp;
			}
			else
				new->state = FASTALIB_CLOSE;
			new->next = NULL;
			if(flp == NULL)
				flp = new;
			else
			{
				curr = flp;
				while(curr->next != NULL)
					curr = curr->next;
				curr->next = new;
			}
		}
		libs = libs->next;
	}
	
	if(flp == NULL)
		return FALSE;
	ObjMgrProcLoad(OMPROC_FETCH, libproc, libproc, OBJ_SEQID, 0, OBJ_BIOSEQ, 0, (Pointer)flp, FastaLibBioseqFetchFunc, PROC_PRIORITY_DEFAULT);
	return TRUE;
}


/***********************************************************************
*
*	FastaLibBioseqFetchDisable()
*	Free the data assoicated with the proc and Free the user data
*	as well.
*
***********************************************************************/
NLM_EXTERN void FastaLibBioseqFetchDisable(void)
{
	ObjMgrPtr omp;
	ObjMgrProcPtr ompp;
	FastaLibPtr flp, next;

	omp = ObjMgrGet();
	ompp = ObjMgrProcFind(omp, 0, libproc, OMPROC_FETCH);
	if(ompp == NULL)
		return;
	ObjMgrFreeUserData(0, ompp->procid, OMPROC_FETCH, 0);

	flp = (FastaLibPtr)(ompp->procdata);
	while(flp)
	{
		if(flp->state == FASTALIB_OPEN)
			FileClose(flp->fp);
		next = flp->next;
		MemFree(flp);
		flp = next;
	}

	return;
}

 

/*********************************************************************
*
*	seqid_to_string(sip, name, use_locus)
*	print the most important field in Seqid to a string stored in 
*	name. 
*
**********************************************************************/
NLM_EXTERN Boolean seqid_to_string(SeqIdPtr sip, CharPtr name, Boolean use_locus)
{
  DbtagPtr db_tag;
  ObjectIdPtr obj_id;
  TextSeqIdPtr tsip;
  PDBSeqIdPtr pip;
  GiimPtr gip;

        switch(sip->choice)
	{
          case 1:       /**local**/
            obj_id = sip->data.ptrvalue;
            if(obj_id->str)
                StringCpy(name, obj_id->str);
            else
                sprintf(name, "%ld", (long) obj_id->id);
            break;

          case 5:       /**genbank**/
          case 6:       /**EMBL**/
          case 7:       /**PIR**/
          case 8:       /**SwissProt**/
          case 10:      /**Other**/
          case 13:      /**DDBJ**/
          case 14:      /**PRF**/
            tsip = sip->data.ptrvalue;
            if(tsip->accession)
                StringCpy(name, tsip->accession);
            if((tsip->name && use_locus) || tsip->accession == NULL)
                StringCpy(name, tsip->name);
 
            break;
 
          case 11:      /**general**/
            db_tag = sip->data.ptrvalue;
            obj_id = db_tag->tag;
            if(obj_id->str)
              StringCpy(name, obj_id->str);
            else
                sprintf(name, "%ld", (long) obj_id->id);
            break;
 
          case 4:       /**giim**/
            gip = sip->data.ptrvalue;
            sprintf(name, "%ld", (long)(gip->id));
            break;

          case 2:     	/*gibbseq*/
	  case 3:	/*gibbmt*/
	  case 12:	/*gi*/
            sprintf(name, "%ld", (long)(sip->data.intvalue));
            break;

          case 15:      /*pdb*/
            pip = sip->data.ptrvalue;
            StringCpy(name, pip->mol);
            break;
	  default:
	    return FALSE;
	}

	return TRUE;
}

/*********************************************************************
*
*	FileBioseqFetchEnable(path, ext)
*	Initiate a BioseqFetch function by either reading an ASN.1
*	Seq-entry file or FASTA file. path->choice determines the
*	type of the file, such as text ASN, binary ASN and FASTA file
*	ext is the extension that is needed to add to the end of the
*	sequence name to make the sequence file
*
*********************************************************************/
 

static CharPtr fileproc = "FileBioseqFetch";
static Int2 LIBCALLBACK FileBioseqFetchFunc (Pointer data)
{
	OMProcControlPtr ompcp;
	ObjMgrProcPtr ompp;
	FileBspPtr fbp;
	SingleBspFilePtr sbfp;
	SeqIdPtr sip;
	OMUserDataPtr omdp;
	CharPtr file_name = NULL;
	Char name[100], f_name[100];
	CharPtr c_name;
	FILE *fp;
	AsnIoPtr aip;
	SeqEntryPtr sep = NULL;
	BioseqPtr bsp;
	Boolean bin;

	ompcp = (OMProcControlPtr)data;
	ompp = ompcp->proc;
	sbfp = (SingleBspFilePtr)(ompp->procdata);

	if(ompcp->input_entityID)
	{
		omdp = ObjMgrGetUserData(ompcp->input_entityID, ompp->procid, OMPROC_FETCH, 0);
		if(omdp != NULL)
			file_name = omdp->userdata.ptrvalue;
	}

	sip = (SeqIdPtr)(ompcp->input_data);
	if(sip == NULL)
		return OM_MSG_RET_ERROR;
	/* if(sip->choice == SEQID_GI)
		return OM_MSG_RET_OK; */
	while(sbfp && sep==NULL)
	{
		fbp = sbfp->data.ptrvalue;
		if(file_name == NULL)
		{
			seqid_to_string(sip, name, fbp->use_locus);
			if(fbp->path)
				sprintf(f_name, "%s%s", fbp->path, name);
			else
				StringCpy(f_name, name);
			if(fbp->ext)
				StringCat(f_name, fbp->ext);
			c_name= f_name;
		}
		else
			c_name = file_name;
		switch(sbfp->choice)
		{
			case FASTA_FILE:
				if((fp = LsqFetch_FileOpen(c_name, "r")) != NULL)
				{
					sep = fasta_lib_sep(fp, NULL, sip);
					FileClose(fp);
				}
				break;

			case TEXT_ASN:	
			case BIN_ASN:
				bin = (sbfp->choice == BIN_ASN);
				if((aip = LsqFetch_AsnIoOpen(c_name, bin?"rb":"r")) != NULL) 
				{
					sep = SeqEntryAsnRead(aip, NULL);
					AsnIoClose(aip);
				}
				break;

			default:
				break;
		}
		sbfp = sbfp->next;
	}

	if(sep == NULL)
		return OM_MSG_RET_OK;

	bsp = BioseqFindInSeqEntry(sip, sep);
	if(bsp == NULL)
		bsp = find_big_bioseq(sep);
	ompcp->output_data = (Pointer)bsp;
	ompcp->output_entityID = ObjMgrGetEntityIDForChoice(sep);
	omdp = ObjMgrAddUserData(ompcp->output_entityID, ompp->procid, OMPROC_FETCH, 0);
	omdp->userdata.ptrvalue = StringSave(file_name);
	omdp->freefunc = FreeSeqName;
	return OM_MSG_RET_DONE;
}

static Boolean path_is_loaded(SingleBspFilePtr head, Uint1 choice,  CharPtr path, CharPtr ext)
{
	FileBspPtr fbp;

	while(head)
	{
		if(head->choice == choice)
		{
			fbp = head->data.ptrvalue;
			if(StringCmp(path, fbp->path) ==0)
				if(StringCmp(ext, fbp->ext) ==0)
					return TRUE;
		}
		head = head->next;
	}
	return FALSE;
}
	

NLM_EXTERN Boolean FileBioseqFetchEnable(ValNodePtr path, ValNodePtr ext)
{
	SingleBspFilePtr sbfp = NULL; 
	FileBspPtr new;
	Char c_path[100], c_ext[20]; 
	CharPtr str;
	Int4 len;
	Char delimiter;
	
	if(path == NULL || ext == NULL)
		return FALSE;
	
	while(path && ext)
	{
		new = MemNew(sizeof(FileBsp));
		new->use_locus = FALSE;
		new->path = NULL;
		new->ext = NULL;
		c_path[0] = '\0';
		c_ext[0] = '\0';

		str = path->data.ptrvalue;
		if(str !=NULL)
		{
			delimiter = DIRDELIMCHR;
			len = StringLen(str);
			if(str[len-1] != delimiter)
				sprintf(c_path, "%s%c", str, delimiter);
			else
				StringCpy(c_path, str);
		}
			
		str = ext->data.ptrvalue;
		if(str !=NULL)
		{
			if(str[0] != '.')
				sprintf(c_ext, ".%s", str);
			else
				StringCpy(c_ext, str);
		}

		if(c_path[0] != '\0' &&  !path_is_loaded(sbfp, path->choice, c_path, c_ext))
		{
			new->path = StringSave(c_path);
			if(c_ext[0] != '\0')
				new->ext = StringSave(c_ext);
			ValNodeAddPointer(&sbfp, path->choice, new);
		}
		else
			MemFree(new);
		path = path->next;
		ext = ext->next;
	}

	ErrSetFatalLevel(SEV_MAX);
	
	ObjMgrProcLoad(OMPROC_FETCH, fileproc, fileproc, OBJ_SEQID, 0, OBJ_BIOSEQ, 0, (Pointer)sbfp, FileBioseqFetchFunc, PROC_PRIORITY_DEFAULT);
	return TRUE;
}


/**********************************************************************
*
*	FileBioseqFetchDisable()
*	Remove the proc associated with FileBioseqFetch and free all the 
*	sequence names in userdata
*
***********************************************************************/
NLM_EXTERN void FileBioseqFetchDisable(void)
{
	ObjMgrPtr omp;
	ObjMgrProcPtr ompp;
	SingleBspFilePtr sbfp, curr;
	FileBspPtr fbp;

	omp = ObjMgrGet();
	ompp = ObjMgrProcFind(omp, 0, fileproc, OMPROC_FETCH);
	if(ompp == NULL)
		return;
	ObjMgrFreeUserData(0, ompp->procid, OMPROC_FETCH, 0);

	sbfp= (SingleBspFilePtr)(ompp->procdata);
	for(curr = sbfp; curr !=NULL; curr = curr->next)
	{
		fbp = curr->data.ptrvalue;
		MemFree(fbp->path);
		MemFree(fbp->ext);
		MemFree(fbp);
	}
	ValNodeFree(sbfp);
	return;
}

static Boolean lib_is_loaded(ValNodePtr head, CharPtr lib_name)
{
	CharPtr str;

	while(head!=NULL)
	{
		str = (CharPtr)(head->data.ptrvalue);
		if(StringCmp(str, lib_name) == 0)
			return TRUE;
		head = head->next;
	}
	return FALSE;
}


static Boolean load_seq_info(CharPtr word, CharPtr val, ValNodePtr PNTR libs, ValNodePtr PNTR path, ValNodePtr PNTR ext)
{

	if(StringICmp(word, "FASTA_LIB") ==0)
	{
		if(!lib_is_loaded(*libs, val))
			ValNodeAddStr(libs, 0, StringSave(val));
		return TRUE;
	}

	if(StringICmp(word, "FASTA_FILE") ==0)
	{
		ValNodeAddStr(path, FASTA_FILE, StringSave(val));
		return TRUE;
	}

	if(StringICmp(word, "BIN_ASN") == 0)
	{
		ValNodeAddStr(path, BIN_ASN, StringSave(val));
		return TRUE;
	}

	if(StringICmp(word, "TEXT_ASN") == 0)
	{
		ValNodeAddStr(path, TEXT_ASN, StringSave(val));
		return TRUE;
	}

	if(StringICmp(word, "EXT") == 0)
	{
		ValNodeAddStr(ext, 0, StringSave(val));
		return TRUE;
	}

	return FALSE;
}



/*********************************************************************
*
*	ReadFetchConfigFile()
*	Reads TYPE, PATH and EXT fields from the LSQFETCH configuration file 
*	TYPE is TEXT_ASN, BIN_ASN, FASTA_FILE, or FASTA_LIB.
*
*   Now also looks for INDEXED_BIN_ASN for easy
*   access to AsnIndexedLib fetch mechanism.
*
*********************************************************************/
		
static void ReadFetchConfigFile (CharPtr sect, ValNodePtr PNTR libs, ValNodePtr PNTR path, ValNodePtr PNTR ext)

{
	Char str [PATH_MAX + 10];
	Char temp [20];

	if (GetAppParam ("LSQFETCH", sect, "PATH", "", str, sizeof (str) - 1)) {
		if (GetAppParam ("LSQFETCH", sect, "TYPE", "", temp, sizeof (temp) - 1)) {

			/* first check for indexed asn tags */

			if (StringICmp (temp, "INDEXED_BIN_ASN") == 0) {
				AsnIndexedLibFetchEnable (str, TRUE);
				return;
			}
			if (StringICmp (temp, "INDEXED_TEXT_ASN") == 0) {
				AsnIndexedLibFetchEnable (str, FALSE);
				return;
			}

			/* now look for regular lsqfetch files */

			load_seq_info (temp, str, libs, path, ext);
			if (GetAppParam ("LSQFETCH", sect, "EXT", "", str, sizeof (str) - 1)) {
				load_seq_info ("EXT", str, libs, path, ext);
			} else {
				load_seq_info ("EXT", NULL, libs, path, ext);
			}
		}
	}
}

/*********************************************************************
*
*	BioseqFetchInit()
*	Initiate BioseqFetch functions from local data and Entrez. 
*	Local data files are stored in a special file, but also in a 
*	config file.  If non is successful, return FALSE
*
*********************************************************************/
NLM_EXTERN Boolean LocalSeqFetchInit(Boolean now)
{
	CharPtr seq_file;
	Boolean success = FALSE;
	ValNodePtr path = NULL, ext = NULL, libs = NULL;
	FILE *fp;
	Char temp[100], str[PATH_MAX + 10];
	CharPtr word, val;
	Char delimiter;
	Int2 i;

	i = 1;
	sprintf (temp, "ORDER_%d", (int) i);
	while (GetAppParam ("LSQFETCH", "ORDER", temp, "", str, sizeof (str) - 1)) {
		ReadFetchConfigFile (str, &libs, &path, &ext);
		i++;
		sprintf (temp, "ORDER_%d", (int) i);
	}

	seq_file = seqinfo_file;	/*check the current search path*/
	if((fp = LsqFetch_FileOpen(seq_file, "r")) != NULL)
	{
   		while(FileGets(str, 100, fp) != NULL)   /*find the right seq*/
		{
			sscanf(str, "%s\n", temp);
			word = strtok(temp, "=");
			if(word !=NULL)
			{
				val = strtok(NULL, "=");
				load_seq_info(word, val, &libs, &path, &ext);
			}
		}
		FileClose(fp);
	}
	else	/*take the current path for searching path of FASTA*/
	{
		delimiter = DIRDELIMCHR;
#ifndef WIN_MAC
			sprintf(str, ".%c", delimiter);
#else
			sprintf(str, "%c", delimiter);
#endif
		ValNodeCopyStr(&path, FASTA_FILE, str);
		ValNodeCopyStr(&ext, 0, ".seq");
	}

	if(libs !=NULL)
		if(FastaLibBioseqFetchEnable(libs, now))
			success = TRUE;
	if(FileBioseqFetchEnable(path, ext))
		success = TRUE;

	ValNodeFreeData(libs);
	ValNodeFreeData(path);
	ValNodeFreeData(ext);

	return success;
}

NLM_EXTERN Boolean BioseqFetchInit(Boolean now)
{
	Boolean success = FALSE;

	/*if(EntrezBioseqFetchEnable ("testseq", now))
		success = TRUE;*/
	if(LocalSeqFetchInit(now))
		success = TRUE;

	return success;
}
/***********************************************************************
*
*	BioseqFetchDisable(): Remove all the functions associated with 
*	BioseqFetch
*
**********************************************************************/
NLM_EXTERN void LocalSeqFetchDisable(void)
{
	FastaLibBioseqFetchDisable();
	FileBioseqFetchDisable();
}

NLM_EXTERN void BioseqFetchDisable(void)
{
	LocalSeqFetchDisable();
	/*EntrezBioseqFetchDisable();*/
}


/**********************************************************************/

typedef struct fastaidx {
  CharPtr       image;
  CharPtr PNTR  seqids;
  CharPtr PNTR  offsets;
  Int4          numlines;
  CharPtr       path;
  CharPtr       file;
} FastaIndex, PNTR FastaIndexPtr;

static Int4 SearchFastaIndex (
  FastaIndexPtr fip,
  CharPtr seqid
)

{
  int       compare;
  Int4      L, R, mid;
  long int  val;

  if (fip == NULL || fip->seqids == NULL || fip->offsets == NULL) return -1;
  if (StringHasNoText (seqid)) return -1;

  L = 0;
  R = fip->numlines - 1;
  while (L < R) {
    mid = (L + R) / 2;
    compare = StringICmp (fip->seqids [mid], seqid);
    if (compare < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  if (R >= 0 && R < fip->numlines) {
    if (StringICmp (fip->seqids [R], seqid) == 0) {
      if (fip->offsets [R] != NULL &&
          sscanf (fip->offsets [R], "%ld", &val) == 1) {
        return (Int4) val;
      }
    }
  }

  return -1;
}

static FastaIndexPtr FreeFastaIndex (
  FastaIndexPtr fip
)

{
  if (fip == NULL) return NULL;

  MemFree (fip->seqids);
  MemFree (fip->offsets);
  MemFree (fip->image);
  MemFree (fip->path);
  MemFree (fip->file);
  MemFree (fip);

  return NULL;
}

static FastaIndexPtr ReadFastaIndex (
  CharPtr file
)

{
  Char           ch;
  FastaIndexPtr  fip;
  FILE           *fp;
  Int4           idx, len;
  CharPtr        last, ptr, tmp;
  Boolean        outOfOrder;
  Char           path [PATH_MAX];

  if (StringHasNoText (file)) return NULL;
  len = FileLength (file);
  if (len < 1) return NULL;

  fip = (FastaIndexPtr) MemNew (sizeof (FastaIndex));
  if (fip == NULL) return NULL;
  
  fp = LsqFetch_FileOpen (file, "r");
  if (fp == NULL) {
    MemFree (fip);
    return NULL;
  }

  fip->image = MemNew ((size_t) (len + 5));
  if (fip->image == NULL) {
    FileClose (fp);
    MemFree (fip);
    return NULL;
  }

  /* read already sorted file into allocated buffer */

  FileRead (fip->image, (size_t) len, 1, fp);

  FileClose (fp);

  /* count lines */

  fip->numlines = 0;
  tmp = fip->image;
  ch = *tmp;
  while (ch != '\0') {
    if (ch == '\n') {
      (fip->numlines++);
    }
    tmp++;
    ch = *tmp;
  }

  fip->seqids = MemNew (sizeof (CharPtr) * (size_t) (fip->numlines + 2));
  fip->offsets = MemNew (sizeof (CharPtr) * (size_t) (fip->numlines + 2));
  if (fip->seqids == NULL || fip->offsets == NULL) {
    return FreeFastaIndex (fip);
  }

  /* initialize seqids and offsets arrays */

  idx = 0;
  tmp = fip->image;
  ch = *tmp;
  fip->seqids [idx] = tmp;
  while (ch != '\0') {
    if (ch == '\n') {
      *tmp = '\0';
      tmp++;
      ch = *tmp;
      idx++;
      fip->seqids [idx] = tmp;
    } else if (ch == '\t') {
      *tmp = '\0';
      tmp++;
      ch = *tmp;
      fip->offsets [idx] = tmp;
    } else {
      tmp++;
      ch = *tmp;
    }
  }

  /* confirm sorted order and uniqueness */

  outOfOrder = FALSE;
  last = fip->seqids [0];
  for (idx = 1; idx < fip->numlines; idx++) {
    if (StringCmp (last, fip->seqids [idx]) >= 0) {
      outOfOrder = TRUE;
    }
    last = fip->seqids [idx];
  }

  if (outOfOrder) {
    ErrPostEx (SEV_ERROR, 0, 0, "FASTA index file out of order");
  }

  /* multi-FASTA file is in same directory as .idx file */

  StringNCpy_0 (path, file, sizeof (path));
  ptr = StringRChr (path, DIRDELIMCHR);
  if (ptr != NULL) {
    *ptr = '\0';
    if (! StringHasNoText (path)) {
      fip->path = StringSave (path);
    }
    ptr++;
    if (! StringHasNoText (ptr)) {
      fip->file = StringSave (ptr);
    }
  }

  return fip;
}

/* human genome object manager registerable fetch function */

static CharPtr fastalibfetchproc = "IndexedFastaLibBioseqFetch";

typedef struct flibftch {
  CharPtr        path;
  CharPtr        fastaname;
  FastaIndexPtr  currentfip;
} FastaLibFetchData, PNTR FastaLibFetchPtr;

static Int2 LIBCALLBACK IndexedFastaLibBioseqFetchFunc (Pointer data)

{
  BioseqPtr         bsp;
  Pointer           dataptr = NULL;
  Uint2             datatype, entityID = 0;
  Char              file [FILENAME_MAX], path [PATH_MAX], id [41], tmp [41];
  FastaLibFetchPtr  flfp;
  FILE              *fp;
  Int4              offset;
  OMProcControlPtr  ompcp;
  ObjMgrProcPtr     ompp;
  CharPtr           ptr;
  SeqEntryPtr       sep = NULL;
  SeqIdPtr          sip;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  ompp = ompcp->proc;
  if (ompp == NULL) return OM_MSG_RET_ERROR;
  flfp = (FastaLibFetchPtr) ompp->procdata;
  if (flfp == NULL) return OM_MSG_RET_ERROR;
  sip = (SeqIdPtr) ompcp->input_data;
  if (sip == NULL) return OM_MSG_RET_ERROR;

  if (sip->choice == SEQID_LOCAL) {

    SeqIdWrite (sip, id, PRINTID_REPORT, sizeof (id));
    if (StringLen (id) > 2 && IS_UPPER (id [0]) && IS_LOWER (id [1])) {
      StringCpy (tmp, id + 2);
      ptr = StringChr (tmp, '_');
      if (ptr != NULL) {
        *ptr = '\0';
        sprintf (file, "chr%s.idx", tmp);

        if (flfp->currentfip == NULL ||
            StringICmp (file, flfp->currentfip->file) != 0) {
          flfp->currentfip = FreeFastaIndex (flfp->currentfip);
          StringNCpy_0 (path, flfp->path, sizeof (path));
          FileBuildPath (path, NULL, file);
          flfp->currentfip = ReadFastaIndex (path);
        }

        if (flfp->currentfip != NULL) {
          offset = SearchFastaIndex (flfp->currentfip, id);
          if (offset < 0) return OM_MSG_RET_ERROR;
          sprintf (file, "chr%s.fa", tmp);
          StringNCpy_0 (path, flfp->path, sizeof (path));
          FileBuildPath (path, NULL, file);
          fp = LsqFetch_FileOpen (path, "r");
          if (fp == NULL) return OM_MSG_RET_ERROR;
          fseek (fp, offset, SEEK_SET);
          dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, &entityID,
                                            FALSE, FALSE, TRUE, FALSE);
          if (dataptr != NULL) {
            sep = GetTopSeqEntryForEntityID (entityID);
          }
          FileClose (fp);
        }
      }
    }
  }

  if (sep == NULL) return OM_MSG_RET_ERROR;
  bsp = BioseqFindInSeqEntry (sip, sep);
  ompcp->output_data = (Pointer) bsp;
  ompcp->output_entityID = ObjMgrGetEntityIDForChoice (sep);
  return OM_MSG_RET_DONE;
}

NLM_EXTERN Boolean IndexedFastaLibFetchEnable (CharPtr path)

{
  FastaLibFetchPtr  flfp;
  Char              str [PATH_MAX];

  StringNCpy_0 (str, path, sizeof (str));
  TrimSpacesAroundString (str);
  flfp = (FastaLibFetchPtr) MemNew (sizeof (FastaLibFetchData));
  if (flfp != NULL) {
    flfp->path = StringSave (str);
    flfp->currentfip = NULL;
  }
  ObjMgrProcLoad (OMPROC_FETCH, fastalibfetchproc, fastalibfetchproc,
                  OBJ_SEQID, 0, OBJ_BIOSEQ, 0, (Pointer) flfp,
                  IndexedFastaLibBioseqFetchFunc, PROC_PRIORITY_DEFAULT);
  return TRUE;
}

NLM_EXTERN void IndexedFastaLibFetchDisable (void)

{
  FastaLibFetchPtr  flfp;
  ObjMgrPtr         omp;
  ObjMgrProcPtr     ompp;

  omp = ObjMgrGet ();
  ompp = ObjMgrProcFind (omp, 0, fastalibfetchproc, OMPROC_FETCH);
  if (ompp == NULL) return;
  ObjMgrFreeUserData (0, ompp->procid, OMPROC_FETCH, 0);
  flfp = (FastaLibFetchPtr) ompp->procdata;
  if (flfp == NULL) return;
  MemFree (flfp->path);
  /* MemFree (flfp->fastaname); */
  FreeFastaIndex (flfp->currentfip);
  MemFree (flfp);
}

/* chimpanzee genome object manager registerable fetch function */

static CharPtr altfastalibfetchproc = "AltIndexedFastaLibBioseqFetch";

typedef struct idfip {
  CharPtr        seqid;
  FastaIndexPtr  fip;
} IdFip, PNTR IdFipPtr;

typedef struct alibftch {
  CharPtr     path;
  ValNodePtr  fiplist;
  IdFipPtr    index;
  Int4        numids;
} AltLibFetchData, PNTR AltLibFetchPtr;

static void ChangeLocalToGenbank (BioseqPtr bsp, Pointer userdata)

{
  Char      id [41], tmp [41];
  SeqIdPtr  sip;

  for (sip = bsp->id; sip != NULL && sip->choice != SEQID_LOCAL; sip = sip->next) continue;
  if (sip == NULL) return;
  SeqIdWrite (sip, id, PRINTID_REPORT, sizeof (id));
  sprintf (tmp, "gb|%s", id);
  sip = SeqIdParse (tmp);
  bsp->id = SeqIdSetFree (bsp->id);
  bsp->id = sip;
  SeqMgrReplaceInBioseqIndex (bsp);
}

static FastaIndexPtr SearchAltIndex (
  AltLibFetchPtr alfp,
  CharPtr seqid
)

{
  int       compare;
  IdFipPtr  ifp;
  Int4      L, R, mid;

  if (alfp == NULL || alfp->index == NULL) return NULL;
  ifp = alfp->index;
  if (StringHasNoText (seqid)) return NULL;

  L = 0;
  R = alfp->numids - 1;
  while (L < R) {
    mid = (L + R) / 2;
    compare = StringICmp (ifp [mid].seqid, seqid);
    if (compare < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  if (R >= 0 && R < alfp->numids) {
    if (StringICmp (ifp [R].seqid, seqid) == 0) {
      return ifp [R].fip;
    }
  }

  return NULL;
}

static Int2 LIBCALLBACK AltIndexedFastaLibBioseqFetchFunc (Pointer data)

{
  AltLibFetchPtr    alfp;
  BioseqPtr         bsp;
  Pointer           dataptr = NULL;
  Uint2             datatype, entityID = 0;
  Char              file [FILENAME_MAX], path [PATH_MAX], id [41];
  FastaIndexPtr     fip;
  FILE              *fp;
  Int4              offset;
  OMProcControlPtr  ompcp;
  ObjMgrProcPtr     ompp;
  SeqEntryPtr       sep = NULL;
  SeqIdPtr          sip;
  CharPtr           tmp;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  ompp = ompcp->proc;
  if (ompp == NULL) return OM_MSG_RET_ERROR;
  alfp = (AltLibFetchPtr) ompp->procdata;
  if (alfp == NULL) return OM_MSG_RET_ERROR;
  sip = (SeqIdPtr) ompcp->input_data;
  if (sip == NULL) return OM_MSG_RET_ERROR;

  if (sip->choice == SEQID_GENBANK || sip->choice == SEQID_GI) {

    SeqIdWrite (sip, id, PRINTID_REPORT, sizeof (id));
    fip = SearchAltIndex (alfp, id);
    if (fip != NULL) {
      offset = SearchFastaIndex (fip, id);
      if (offset < 0) return OM_MSG_RET_ERROR;
      StringCpy (file, fip->file);
      tmp = StringStr (file, ".idx");
      if (tmp != NULL) {
        *tmp = '\0';
      }
      StringCat (file, ".fa");
      StringNCpy_0 (path, fip->path, sizeof (path));
      FileBuildPath (path, NULL, file);
      fp = LsqFetch_FileOpen (path, "r");
      if (fp == NULL) {
        tmp = StringStr (file, ".fa");
        if (tmp != NULL) {
          *tmp = '\0';
          StringCat (file, ".fsa");
          StringNCpy_0 (path, fip->path, sizeof (path));
          FileBuildPath (path, NULL, file);
          fp = LsqFetch_FileOpen (path, "r");
        }
      }
      if (fp == NULL) return OM_MSG_RET_ERROR;
      fseek (fp, offset, SEEK_SET);
      dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, &entityID,
                                        FALSE, FALSE, TRUE, FALSE);
      if (dataptr != NULL) {
        sep = GetTopSeqEntryForEntityID (entityID);
      }
      FileClose (fp);
    }
  }

  if (sep == NULL) return OM_MSG_RET_ERROR;
  VisitBioseqsInSep (sep, NULL, ChangeLocalToGenbank);
  bsp = BioseqFindInSeqEntry (sip, sep);
  ompcp->output_data = (Pointer) bsp;
  ompcp->output_entityID = ObjMgrGetEntityIDForChoice (sep);
  return OM_MSG_RET_DONE;
}

static int LIBCALLBACK SortIfpByID (VoidPtr vp1, VoidPtr vp2)

{
  IdFipPtr ifp1, ifp2;

  if (vp1 == NULL || vp2 == NULL) return 0;
  ifp1 = (IdFipPtr) vp1;
  ifp2 = (IdFipPtr) vp2;
  if (ifp1 == NULL || ifp2 == NULL) return 0;
  return StringICmp (ifp1->seqid, ifp2->seqid);
}

NLM_EXTERN Boolean AltIndexedFastaLibFetchEnable (CharPtr path)

{
  AltLibFetchPtr  alfp = NULL;
  Char            file [FILENAME_MAX];
  FastaIndexPtr   fip;
  ValNodePtr      head;
  Int4            i;
  IdFipPtr        ifp;
  Boolean         is_new = FALSE;
  Int4            j;
  Int4            numids = 0;
  ObjMgrPtr       omp;
  ObjMgrProcPtr   ompp;
  Char            str [PATH_MAX];
  CharPtr         tmp;
  ValNodePtr      vnp;

  StringNCpy_0 (str, path, sizeof (str));
  TrimSpacesAroundString (str);
  omp = ObjMgrGet ();
  ompp = ObjMgrProcFind (omp, 0, altfastalibfetchproc, OMPROC_FETCH);
  if (ompp != NULL) {
    alfp = (AltLibFetchPtr) ompp->procdata;
    if (alfp != NULL) {
      alfp->path = MemFree (alfp->path);
      for (vnp = alfp->fiplist; vnp != NULL; vnp = vnp->next) {
        fip = (FastaIndexPtr) vnp->data.ptrvalue;
        FreeFastaIndex (fip);
      }
      alfp->fiplist = ValNodeFree (alfp->fiplist);
      alfp->index = MemFree (alfp->index);
    }
  } else {
    alfp = (AltLibFetchPtr) MemNew (sizeof (AltLibFetchData));
    is_new = TRUE;
  }
  if (alfp != NULL) {
    alfp->path = StringSave (str);
    head = DirCatalog (str);
    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      if (vnp->choice == 0) {
        tmp = (CharPtr) vnp->data.ptrvalue;
        if (StringStr (tmp, ".idx") != NULL) {
          StringCpy (str, alfp->path);
          sprintf (file, "%s", tmp);
          FileBuildPath (str, NULL, file);
          fip = ReadFastaIndex (str);
          if (fip != NULL) {
            ValNodeAddPointer (&(alfp->fiplist), 0, (Pointer) fip);
            numids += fip->numlines;
          }
        }
      }
    }
    ValNodeFreeData (head);
    ifp = (IdFipPtr) MemNew (sizeof (IdFip) * (numids + 2));
    alfp->index = ifp;
    alfp->numids = numids;
    if (ifp != NULL) {
      i = 0;
      for (vnp = alfp->fiplist; vnp != NULL; vnp = vnp->next) {
        fip = (FastaIndexPtr) vnp->data.ptrvalue;
        if (fip != NULL) {
          for (j = 0; j < fip->numlines; j++, i++) {
            ifp [i].seqid = fip->seqids [j];
            ifp [i].fip = fip;
          }
        }
      }
      HeapSort (ifp, (size_t) numids, sizeof (IdFip), SortIfpByID);
    }
  }
  if (is_new) {
    ObjMgrProcLoad (OMPROC_FETCH, altfastalibfetchproc, altfastalibfetchproc,
                    OBJ_SEQID, 0, OBJ_BIOSEQ, 0, (Pointer) alfp,
                    AltIndexedFastaLibBioseqFetchFunc, PROC_PRIORITY_DEFAULT);
  }
  return TRUE;
}

NLM_EXTERN void AltIndexedFastaLibFetchDisable (void)

{
  AltLibFetchPtr  alfp;
  FastaIndexPtr   fip;
  ObjMgrPtr       omp;
  ObjMgrProcPtr   ompp;
  ValNodePtr      vnp;

  omp = ObjMgrGet ();
  ompp = ObjMgrProcFind (omp, 0, altfastalibfetchproc, OMPROC_FETCH);
  if (ompp == NULL) return;
  ObjMgrFreeUserData (0, ompp->procid, OMPROC_FETCH, 0);
  alfp = (AltLibFetchPtr) ompp->procdata;
  if (alfp == NULL) return;
  alfp->path = MemFree (alfp->path);
  for (vnp = alfp->fiplist; vnp != NULL; vnp = vnp->next) {
    fip = (FastaIndexPtr) vnp->data.ptrvalue;
    FreeFastaIndex (fip);
  }
  alfp->fiplist = ValNodeFree (alfp->fiplist);
  alfp->index = MemFree (alfp->index);
  MemFree (alfp);
}

/* common function for creating indexes of fasta library files */

NLM_EXTERN void CreateFastaIndex (
  CharPtr file
)

{
  BioseqPtr    bsp;
  Pointer      dataptr = NULL;
  Uint2        datatype, entityID = 0;
  ValNodePtr   head = NULL, last = NULL, vnp;
  Char         id [41], path [PATH_MAX], tmp [64];
  FILE         *ifp, *ofp;
  Int4         offset;
  CharPtr      ptr;
  SeqEntryPtr  sep;

  if (StringHasNoText (file)) return;

  /* replace extension by .idx for index file */

  StringNCpy_0 (path, file, sizeof (path));
  ptr = StringRChr (path, '.');
  if (ptr != NULL) {
    *ptr = '\0';
  }
  StringCat (path, ".idx");

  ifp = LsqFetch_FileOpen (file, "r");
  if (ifp == NULL) return;

  ofp = LsqFetch_FileOpen (path, "w");
  if (ofp != NULL) {

    /* get initial file offset */

    offset = ftell (ifp);

    /* read next FASTA component */

    while ((dataptr = ReadAsnFastaOrFlatFile (ifp, &datatype, &entityID,
                                              FALSE, FALSE, TRUE, FALSE)) != NULL) {

      sep = GetTopSeqEntryForEntityID (entityID);
      if (sep != NULL) {
        bsp = FindNucBioseq (sep);
        if (bsp != NULL && bsp->length > 0) {
          SeqIdWrite (bsp->id, id, PRINTID_REPORT, sizeof (id));
          if (! StringHasNoText (id)) {

            /* save ID and offset separated by tab character */

            sprintf (tmp, "%s\t%ld", id, (long) offset);
            last = ValNodeNew (last);
            if (head == NULL) {
              head = last;
            }
            if (last != NULL) {
              last->data.ptrvalue = StringSave (tmp);
            }
          }
        }
      }

      ObjMgrFreeByEntityID (entityID);

      /* get file offset of next FASTA component */

      offset = ftell (ifp);
    }

    /* sort by ID */

    head = ValNodeSort (head, SortVnpByString);
    head = UniqueValNode (head);

    /* write ID and offset index */

    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      fprintf (ofp, "%s\n", (CharPtr) vnp->data.ptrvalue);
    }

    FileClose (ofp);
  }

  FileClose (ifp);

  ValNodeFreeData (head);
}

/* object manager registerable fetch function for local ASN.1 indexed files */

static CharPtr asnlibfetchproc = "AsnIndexedLibBioseqFetch";

typedef struct asnlibftch {
  CharPtr     path;
  ValNodePtr  fiplist;
  IdFipPtr    index;
  Int4        numids;
  Boolean     binary;
} AsnLibFetchData, PNTR AsnLibFetchPtr;

static FastaIndexPtr SearchAsnIndex (
  AsnLibFetchPtr alfp,
  CharPtr seqid
)

{
  int       compare;
  IdFipPtr  ifp;
  Int4      L, R, mid;

  if (alfp == NULL || alfp->index == NULL) return NULL;
  ifp = alfp->index;
  if (StringHasNoText (seqid)) return NULL;

  L = 0;
  R = alfp->numids - 1;
  while (L < R) {
    mid = (L + R) / 2;
    compare = StringICmp (ifp [mid].seqid, seqid);
    if (compare < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  if (R >= 0 && R < alfp->numids) {
    if (StringICmp (ifp [R].seqid, seqid) == 0) {
      return ifp [R].fip;
    }
  }

  return NULL;
}

static Int2 LIBCALLBACK AsnIndexedLibBioseqFetchFunc (Pointer data)

{
  AsnIoPtr          aip;
  AsnLibFetchPtr    alfp;
  BioseqPtr         bsp;
  Char              file [FILENAME_MAX], path [PATH_MAX], id [41];
  FastaIndexPtr     fip;
  Int4              offset;
  OMProcControlPtr  ompcp;
  ObjMgrProcPtr     ompp;
  SeqEntryPtr       sep = NULL;
  SeqIdPtr          sip;
  CharPtr           tmp;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  ompp = ompcp->proc;
  if (ompp == NULL) return OM_MSG_RET_ERROR;
  alfp = (AsnLibFetchPtr) ompp->procdata;
  if (alfp == NULL) return OM_MSG_RET_ERROR;
  sip = (SeqIdPtr) ompcp->input_data;
  if (sip == NULL) return OM_MSG_RET_ERROR;

  SeqIdWrite (sip, id, PRINTID_REPORT, sizeof (id));
  fip = SearchAsnIndex (alfp, id);
  if (fip != NULL) {
    offset = SearchFastaIndex (fip, id);
    if (offset < 0) return OM_MSG_RET_ERROR;
    StringCpy (file, fip->file);
    tmp = StringStr (file, ".idx");
    if (tmp != NULL) {
      *tmp = '\0';
    }
    StringCat (file, ".aso");
    StringNCpy_0 (path, fip->path, sizeof (path));
    FileBuildPath (path, NULL, file);
    aip = LsqFetch_AsnIoOpen (path, alfp->binary? "rb" : "r");
    if (aip == NULL) {
      tmp = StringStr (file, ".aso");
      if (tmp != NULL) {
        *tmp = '\0';
        StringCat (file, ".asn");
        StringNCpy_0 (path, fip->path, sizeof (path));
        FileBuildPath (path, NULL, file);
        aip = LsqFetch_AsnIoOpen (path, alfp->binary? "rb" : "r");
      }
    }
    if (aip == NULL) return OM_MSG_RET_ERROR;
    AsnIoSeek (aip, offset);
    sep = SeqEntryAsnRead (aip, NULL);
    AsnIoClose (aip);
  }

  if (sep == NULL) return OM_MSG_RET_ERROR;
  bsp = BioseqFindInSeqEntry (sip, sep);
  ompcp->output_data = (Pointer) bsp;
  ompcp->output_entityID = ObjMgrGetEntityIDForChoice (sep);
  return OM_MSG_RET_DONE;
}

NLM_EXTERN Boolean AsnIndexedLibFetchEnable (CharPtr path, Boolean binary)

{
  AsnLibFetchPtr  alfp = NULL;
  Char            file [FILENAME_MAX];
  FastaIndexPtr   fip;
  ValNodePtr      head;
  Int4            i;
  IdFipPtr        ifp;
  Boolean         is_new = FALSE;
  Int4            j;
  Int4            numids = 0;
  ObjMgrPtr       omp;
  ObjMgrProcPtr   ompp;
  Char            str [PATH_MAX];
  CharPtr         tmp;
  ValNodePtr      vnp;

  StringNCpy_0 (str, path, sizeof (str));
  TrimSpacesAroundString (str);
  omp = ObjMgrGet ();
  ompp = ObjMgrProcFind (omp, 0, asnlibfetchproc, OMPROC_FETCH);
  if (ompp != NULL) {
    alfp = (AsnLibFetchPtr) ompp->procdata;
    if (alfp != NULL) {
      alfp->path = MemFree (alfp->path);
      for (vnp = alfp->fiplist; vnp != NULL; vnp = vnp->next) {
        fip = (FastaIndexPtr) vnp->data.ptrvalue;
        FreeFastaIndex (fip);
      }
      alfp->fiplist = ValNodeFree (alfp->fiplist);
      alfp->index = MemFree (alfp->index);
    }
  } else {
    alfp = (AsnLibFetchPtr) MemNew (sizeof (AsnLibFetchData));
    is_new = TRUE;
    if (alfp != NULL) {
      alfp->binary = binary;
    }
  }
  if (alfp != NULL) {
    alfp->path = StringSave (str);
    head = DirCatalog (str);
    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      if (vnp->choice == 0) {
        tmp = (CharPtr) vnp->data.ptrvalue;
        if (StringStr (tmp, ".idx") != NULL) {
          StringCpy (str, alfp->path);
          sprintf (file, "%s", tmp);
          FileBuildPath (str, NULL, file);
          fip = ReadFastaIndex (str);
          if (fip != NULL) {
            ValNodeAddPointer (&(alfp->fiplist), 0, (Pointer) fip);
            numids += fip->numlines;
          }
        }
      }
    }
    ValNodeFreeData (head);
    ifp = (IdFipPtr) MemNew (sizeof (IdFip) * (numids + 2));
    alfp->index = ifp;
    alfp->numids = numids;
    if (ifp != NULL) {
      i = 0;
      for (vnp = alfp->fiplist; vnp != NULL; vnp = vnp->next) {
        fip = (FastaIndexPtr) vnp->data.ptrvalue;
        if (fip != NULL) {
          for (j = 0; j < fip->numlines; j++, i++) {
            ifp [i].seqid = fip->seqids [j];
            ifp [i].fip = fip;
          }
        }
      }
      HeapSort (ifp, (size_t) numids, sizeof (IdFip), SortIfpByID);
    }
  }
  if (is_new) {
    ObjMgrProcLoad (OMPROC_FETCH, asnlibfetchproc, asnlibfetchproc,
                    OBJ_SEQID, 0, OBJ_BIOSEQ, 0, (Pointer) alfp,
                    AsnIndexedLibBioseqFetchFunc, PROC_PRIORITY_DEFAULT);
  }
  return TRUE;
}

NLM_EXTERN void AsnIndexedLibFetchDisable (void)

{
  AsnLibFetchPtr  alfp;
  FastaIndexPtr   fip;
  ObjMgrPtr       omp;
  ObjMgrProcPtr   ompp;
  ValNodePtr      vnp;

  omp = ObjMgrGet ();
  ompp = ObjMgrProcFind (omp, 0, asnlibfetchproc, OMPROC_FETCH);
  if (ompp == NULL) return;
  ObjMgrFreeUserData (0, ompp->procid, OMPROC_FETCH, 0);
  alfp = (AsnLibFetchPtr) ompp->procdata;
  if (alfp == NULL) return;
  alfp->path = MemFree (alfp->path);
  for (vnp = alfp->fiplist; vnp != NULL; vnp = vnp->next) {
    fip = (FastaIndexPtr) vnp->data.ptrvalue;
    FreeFastaIndex (fip);
  }
  alfp->fiplist = ValNodeFree (alfp->fiplist);
  alfp->index = MemFree (alfp->index);
  MemFree (alfp);
}

/* common function for creating indexes of ASN.1 Bioseq-set ftp release files */

typedef struct asnidxdata {
  FILE        *ofp;
  Int4        offset;
  ValNodePtr  head;
  ValNodePtr  last;
} AsnIdxData, PNTR AsnIdxPtr;

static void SaveAsnIdxOffset (
  BioseqPtr bsp,
  Pointer userdata
)

{
  AsnIdxPtr  aip;
  Char       id [41], tmp [64];
  SeqIdPtr   sip;

  aip = (AsnIdxPtr) userdata;
  if (bsp == NULL || aip == NULL) return;

  sip = SeqIdFindBest (bsp->id, SEQID_GI);
  if (sip == NULL) {
    sip = SeqIdFindBest (bsp->id, 0);
  }

  SeqIdWrite (sip, id, PRINTID_REPORT, sizeof (id));
  if (! StringHasNoText (id)) {

    /* save ID and offset separated by tab character */

    sprintf (tmp, "%s\t%ld", id, (long) aip->offset);
    aip->last = ValNodeNew (aip->last);
    if (aip->head == NULL) {
      aip->head = aip->last;
    }
    if (aip->last != NULL) {
      aip->last->data.ptrvalue = StringSave (tmp);
    }
  }
}

static void CreateBinaryAsnIndex (
  CharPtr file
)

{
  AsnIdxData    aid;
  AsnIoPtr      aip;
  AsnModulePtr  amp;
  AsnTypePtr    atp, atp_bss, atp_se;
  FILE          *ofp;
  ObjMgrPtr     omp;
  Char          path [PATH_MAX];
  CharPtr       ptr;
  SeqEntryPtr   sep;
  ValNodePtr    vnp;

  if (StringHasNoText (file)) return;

  /* replace extension by .idx for index file */

  StringNCpy_0 (path, file, sizeof (path));
  ptr = StringRChr (path, '.');
  if (ptr != NULL) {
    *ptr = '\0';
  }
  StringCat (path, ".idx");

  aip = LsqFetch_AsnIoOpen (file, "rb");
  if (aip == NULL) return;

  ofp = LsqFetch_FileOpen (path, "w");
  if (ofp != NULL) {

    MemSet ((Pointer) &aid, 0, sizeof (AsnIdxData));
    aid.head = NULL;
    aid.last = NULL;
    aid.ofp = ofp;

    amp = AsnAllModPtr ();

    atp_bss = AsnFind ("Bioseq-set");
    atp_se = AsnFind ("Bioseq-set.seq-set.E");

    atp = atp_bss;

    /* get initial file offset */

    aid.offset = AsnIoTell (aip);

    /* read next ASN.1 component */

    while ((atp = AsnReadId (aip, amp, atp)) != NULL) {
      if (atp == atp_se) {

        sep = SeqEntryAsnRead (aip, atp);
        VisitBioseqsInSep (sep, (Pointer) &aid, SaveAsnIdxOffset);

        SeqEntryFree (sep);
        omp = ObjMgrGet ();
        ObjMgrReapOne (omp);
        ObjMgrFreeCache (0);

      } else {

        AsnReadVal (aip, atp, NULL);
      }

      /* get file offset of next ASN.1 component */

      aid.offset = AsnIoTell (aip);
    }

    /* sort by ID */

    aid.head = ValNodeSort (aid.head, SortVnpByString);
    aid.head = UniqueValNode (aid.head);

    /* write ID and offset index */

    for (vnp = aid.head; vnp != NULL; vnp = vnp->next) {
      fprintf (ofp, "%s\n", (CharPtr) vnp->data.ptrvalue);
    }

    FileClose (ofp);
  }

  AsnIoClose (aip);

  ValNodeFreeData (aid.head);
}

static void CreateTextAsnIndex (
  CharPtr file
)

{
  AsnIdxData   aid;
  Pointer      dataptr = NULL;
  Uint2        datatype, entityID = 0;
  FILE         *ifp, *ofp;
  Char         path [PATH_MAX];
  CharPtr      ptr;
  SeqEntryPtr  sep;
  ValNodePtr   vnp;

  if (StringHasNoText (file)) return;

  /* replace extension by .idx for index file */

  StringNCpy_0 (path, file, sizeof (path));
  ptr = StringRChr (path, '.');
  if (ptr != NULL) {
    *ptr = '\0';
  }
  StringCat (path, ".idx");

  ifp = LsqFetch_FileOpen (file, "r");
  if (ifp == NULL) return;

  ofp = LsqFetch_FileOpen (path, "w");
  if (ofp != NULL) {

    MemSet ((Pointer) &aid, 0, sizeof (AsnIdxData));
    aid.head = NULL;
    aid.last = NULL;
    aid.ofp = ofp;

    /* get initial file offset */

    aid.offset = ftell (ifp);

    /* read next ASN.1 component */

    while ((dataptr = ReadAsnFastaOrFlatFile (ifp, &datatype, &entityID,
                                              FALSE, FALSE, TRUE, FALSE)) != NULL) {

      sep = GetTopSeqEntryForEntityID (entityID);
      VisitBioseqsInSep (sep, (Pointer) &aid, SaveAsnIdxOffset);

      ObjMgrFreeByEntityID (entityID);

      /* get file offset of next ASN.1 component */

      aid.offset = ftell (ifp);
    }

    /* sort by ID */

    aid.head = ValNodeSort (aid.head, SortVnpByString);
    aid.head = UniqueValNode (aid.head);

    /* write ID and offset index */

    for (vnp = aid.head; vnp != NULL; vnp = vnp->next) {
      fprintf (ofp, "%s\n", (CharPtr) vnp->data.ptrvalue);
    }

    FileClose (ofp);
  }

  FileClose (ifp);

  ValNodeFreeData (aid.head);
}

NLM_EXTERN void CreateAsnIndex (
  CharPtr file,
  Boolean binary
)

{
  if (binary) {
    CreateBinaryAsnIndex (file);
  } else {
    CreateTextAsnIndex (file);
  }
}

