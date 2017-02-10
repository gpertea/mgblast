/*  valid.h
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
* File Name:  valid.h
*
* Author:  James Ostell
*   
* Version Creation Date: 1/1/94
*
* $Revision: 6.23 $
*
* File Description:  Sequence editing utilities
*
* Modifications:  
* --------------------------------------------------------------------------
* Date	   Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
* $Log: valid.h,v $
* Revision 6.23  2006/04/21 17:59:18  kans
* added ignoreExceptions flag to vsp - for MrnaTransCheck and CdTransCheck
*
* Revision 6.22  2006/02/16 19:34:47  kans
* use vsp->is_smupd_in_sep to suppress ERR_SEQ_FEAT_FeatureRefersToAccession
*
* Revision 6.21  2006/01/26 19:54:26  kans
* added ERR_SEQ_FEAT_FeatureRefersToAccession to look for inconsistent use of gi and accession (with or without version) for sfp->location or sfp->product references in a single blob
*
* Revision 6.20  2005/06/08 15:26:06  kans
* added is_htg_in_sep and is_refseq_in_sep flags to vsp to avoid repetitive checks
*
* Revision 6.19  2004/12/23 20:50:51  kans
* added context field to new callback
*
* Revision 6.18  2004/12/22 21:56:40  kans
* CustValErr supports ValidErrorFunc callback for finer error reporting
*
* Revision 6.17  2004/12/20 22:57:16  kans
* added verbosityLevel argument - to be used for finer control over error reporting by asn2val
*
* Revision 6.16  2004/10/04 15:50:22  kans
* added vsp->justShowAccession for extremely terse output
*
* Revision 6.15  2004/09/10 17:52:05  kans
* changed ValidateLimit enum to Int2 defines
*
* Revision 6.14  2004/09/10 15:31:43  kans
* added farFetchMRNAproducts, locusTagGeneralMatch, validationLimit flags to vsp
*
* Revision 6.13  2004/05/06 19:42:22  kans
* new function GetValidCountryList for access to country code list, which is now NULL terminated
*
* Revision 6.12  2003/12/02 15:37:37  kans
* added vsp->seqSubmitParent for use by tbl2asn, which usually has a Seq-submit wrapper that is added on-the-fly and not indexed
*
* Revision 6.11  2003/11/14 18:07:17  kans
* alignment parameters to find remote bsp, do seqhist assembly
*
* Revision 6.10  2003/03/05 18:47:45  ford
* Added prototype for IsNuclAcc().
*
* Revision 6.9  2003/02/18 20:19:29  kans
* added vsp->validateIDSet, initial work on validating update against ID set in database
*
* Revision 6.8  2002/10/28 19:30:34  kans
* added farFetchCDSproducts to vsp
*
* Revision 6.7  2002/07/16 17:13:09  kans
* added sourceQualTags to vsp, ERR_SEQ_DESCR_StructuredSourceNote by finite state machine text search
*
* Revision 6.6  2002/03/12 22:27:00  kans
* added alwaysRequireIsoJTA
*
* Revision 6.5  2000/02/14 15:00:19  kans
* added vsp->farIDsInAlignments for use by alignment validator
*
* Revision 6.4  1999/12/24 01:21:06  kans
* added validateAlignments flag controlling call to ValidateSeqAlignWithinValidator
*
* Revision 6.3  1999/07/22 22:04:36  kans
* added suppressContext flag
*
* Revision 6.2  1999/05/03 20:06:36  kans
* if no pubs or no biosource, report only once, not once per bioseq
*
* Revision 6.1  1998/07/02 17:53:47  kans
* useSeqMgrIndexes field added to ValidStructPtr, validator can use either old (nested gathers) or new (SeqMgr indexing) method
*
* Revision 6.0  1997/08/25 18:08:29  madden
* Revision changed to 6.0
*
* Revision 5.1  1997/06/19 18:39:54  vakatov
* [WIN32,MSVC++]  Adopted for the "NCBIOBJ.LIB" DLL'ization
*
* Revision 5.0  1996/05/28 13:23:23  ostell
* Set to revision 5.0
*
 * Revision 4.1  1996/01/23  23:10:10  kans
 * added onlyspell and justwarnonspell to ValidStructPtr
 *
 * Revision 4.0  1995/07/26  13:49:01  ostell
 * force revision to 4.0
 *
 * Revision 1.4  1995/06/03  13:45:47  ostell
 * changes made in valid to use gather functions and ErrPostItem instead
 * of previous custom functions
 *
 * Revision 1.3  1995/05/15  21:46:05  ostell
 * added Log line
 *
*
*
* ==========================================================================
*/

#ifndef _NCBI_Valid_
#define _NCBI_Valid_

/*****************************************************************************
*
*   valid.h
*   	values for cutoff are
*       0 INFO
*       1 WARN
*       2 ERROR
*       3 FATAL
*
*****************************************************************************/
#ifndef _NCBI_Seqport_
#include <seqport.h>
#endif

#ifndef _GATHER_
#include <gather.h>
#endif

#ifndef _SQNUTILS_
#include <sqnutils.h>
#endif

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef void (*SpellCallBackFunc) (char * str);
typedef int (* SpellCheckFunc) (char *String, SpellCallBackFunc);

/* callback type for finer error reporting */

typedef void (LIBCALLBACK *ValidErrorFunc) (
  ErrSev severity,
  int errcode,
  int subcode,
  CharPtr accession,
  CharPtr message,
  CharPtr objtype,
  CharPtr label,
  CharPtr context,
  CharPtr location,
  CharPtr product,
  Pointer userdata
);

#define SET_DEPTH 20

#define VALIDATE_ALL 0
#define VALIDATE_INST 1
#define VALIDATE_HIST 2
#define VALIDATE_CONTEXT 3
#define VALIDATE_GRAPH 4
#define VALIDATE_SET 5
#define VALIDATE_FEAT 6
#define VALIDATE_DESC 7

typedef struct validstruct {
	Int2 cutoff;                   /* lowest errmsg to show 0=default */
	Int2 errors[6];
	SeqEntryPtr sep;               /* top level SeqEntryPtr */
	BioseqSetPtr bssp;			   /* current bioseqset */
	BioseqPtr bsp;                 /* current bioseq */
	SeqFeatPtr sfp;                /* current feature */
	ValNodePtr descr;              /* current descriptor */
	Uint4 descrs [SET_DEPTH];      /* bit flags set by descriptor type */
	Int2 protcnt, nuccnt, segcnt;
	CharPtr errbuf;
	Boolean patch_seq;             /* repair invalid sequence residues? */
	Boolean non_ascii_chars;       /* non ascii chars found in read? */
	Boolean suppress_no_pubs;      /* internal use for no pub anywhere message */
	Boolean suppress_no_biosrc;    /* internal use for no biosource anywhere message */
	SpellCheckFunc spellfunc;
	SpellCallBackFunc spellcallback;
	GatherContextPtr gcp;          /* used for reporting the errors */
								   /* this section used for checking Bioseqs */
	Uint2 bsp_partial_val;         /* return from SeqLocPartial on segmented SeqLoc */
	Boolean onlyspell;             /* only do spell check */
	Boolean justwarnonspell;       /* severity WARNING instead of ERROR on spell */
	Boolean useSeqMgrIndexes;      /* new style indexing to speed up validation */
	Boolean suppressContext;       /* suppress context part of message */
	Boolean validateAlignments;    /* call alignval test suite */
	Boolean farIDsInAlignments;    /* fetch to get far IDs in alignments */
	Boolean alignFindRemoteBsp;    /* do remote fetching in alignment validation */
	Boolean doSeqHistAssembly;     /* do alignment validation in Seq-hist.assembly */
	Boolean alwaysRequireIsoJTA;   /* force check for iso_jta */
	Boolean farFetchCDSproducts;   /* lock CDS->products for CdTransCheck, if necessary */
	Boolean farFetchMRNAproducts;  /* lock MRNA->products for MrnaTransCheck, if necessary */
	Boolean locusTagGeneralMatch;  /* expect locus_tag to match Seq-id.general of CDS and mRNA product */
	Boolean validateIDSet;         /* look for gain or loss of general IDs on sequence update */
	Boolean seqSubmitParent;       /* flag from tbl2asn to suppress no pub message */
	Boolean justShowAccession;     /* extremely terse output with accession and error type */
	Boolean ignoreExceptions;      /* report translation and transcription problems even if exception set */
	Int2 validationLimit;          /* limit validation to major classes in Valid1GatherProc */
								   /* this section used for finer error reporting callback */
	ValidErrorFunc errfunc;
	Pointer userdata;
	Boolean convertGiToAccn;
								   /* this section used for internal flags */
	TextFsaPtr sourceQualTags;     /* for detecting structured qual tags in notes */
	Boolean is_htg_in_sep;         /* record has technique of htgs 0 through htgs 3 */
	Boolean is_refseq_in_sep;      /* record has seqid of type other (refseq) */
	Boolean is_smupd_in_sep;       /* record in INSD internal processing */
	Boolean feat_loc_has_gi;       /* at least one feature has a gi location reference */
	Boolean feat_prod_has_gi;      /* at least one feature has a gi product reference */
} ValidStruct, PNTR ValidStructPtr;

NLM_EXTERN Boolean ValidateSeqEntry PROTO((SeqEntryPtr sep, ValidStructPtr vsp));
NLM_EXTERN void ValidStructClear (ValidStructPtr vsp);  /* 0 out a ValidStruct */
NLM_EXTERN ValidStructPtr ValidStructNew (void);
NLM_EXTERN ValidStructPtr ValidStructFree (ValidStructPtr vsp);
NLM_EXTERN void SpellCallBack (char * str);
NLM_EXTERN Boolean IsNuclAcc (CharPtr name);

NLM_EXTERN CharPtr GetValidCategoryName (int errcode);
NLM_EXTERN CharPtr GetValidErrorName (int errcode, int subcode);

NLM_EXTERN CharPtr PNTR GetValidCountryList (void);

#ifdef __cplusplus
}
#endif

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif
