/*
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
* File Name: findrepl.c
*
* Author:  Jonathan Kans, Tim Ford
*
* Version Creation Date:   10/17/00
*
* File Description:
*   Complete redesign of find/replace from original of Yuri Sadykov
*
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
* ==========================================================================
*
*
* RCS Modification History:
* -------------------------
* $Log: findrepl.c,v $
* Revision 6.21  2006/01/17 17:50:01  bollin
* allow FindReplaceInEntity to search for a string made up of whitespace, as
* long as whole_word is not specified
*
* Revision 6.20  2006/01/10 18:13:56  kans
* FindReplAligns does not have case for SAS_DISC, since visit function recursively presents these components separately
*
* Revision 6.19  2006/01/09 21:15:03  bollin
* allow punctuation to terminate a word in find replace
*
* Revision 6.18  2006/01/04 21:26:57  kans
* FSA hit does not need code from validator unstructured source test, cleaned up variable names
*
* Revision 6.17  2006/01/04 20:39:41  kans
* added FindStringsInEntity using finite state machine, general cleanup of code
*
* Revision 6.16  2005/12/29 21:42:06  kans
* only call callback if text was found or replaced
*
* Revision 6.15  2005/12/29 20:54:41  kans
* FindReplaceInEntity takes callback and userdata
*
* Revision 6.14  2005/09/21 14:39:09  bollin
* fixed bug in FindReplace where if the whole-word flag was specified but
* the substring was found in a not-whole-word context earlier in the string
* the search terminated early with a null result
*
* Revision 6.13  2005/04/26 21:33:52  kans
* added SEQID_GPIPE
*
* Revision 6.12  2004/04/01 13:43:05  lavr
* Spell "occurred", "occurrence", and "occurring"
*
* Revision 6.11  2003/11/21 17:58:46  bollin
* remove tax ref and common name when changing taxonomy name via ASN Find/Replace
*
* Revision 6.10  2003/07/31 20:54:54  kans
* FindReplaceString does not need do_replace argument
*
* Revision 6.9  2003/07/31 18:18:03  kans
* added FindReplaceString
*
* Revision 6.8  2003/05/11 21:12:50  kans
* FindReplAligns loops through StdSegPtr chain, also does ssp->ids within
*
* Revision 6.7  2002/06/11 14:41:20  kans
* added support for locus_tag
*
* Revision 6.6  2002/03/05 21:11:04  kans
* do not replace ifp->key
*
* Revision 6.5  2001/12/12 17:38:38  kans
* added new subsource qualifiers, four now have empty name
*
* Revision 6.4  2001/12/07 13:49:34  kans
* workingBuffer needs to be large enough for terminal null byte
*
* Revision 6.3  2001/08/06 22:13:12  kans
* using NUM_SEQID, added TPA ids to arrays
*
* Revision 6.2  2000/11/03 20:36:00  kans
* FindReplaceInEntity replaces FindInEntity and FindInEntityX - complete redesign,
* no longer using AsnExpOptExplore because of the difficulty of replacing with a
* larger string (TF + JK)
*
* Revision 6.1  1999/03/05 23:31:07  kans
* FindInEntityX was not initializing flen, replen
*
* Revision 6.0  1997/08/25 18:05:38  madden
* Revision changed to 6.0
*
* Revision 5.3  1997/06/19 18:37:41  vakatov
* [WIN32,MSVC++]  Adopted for the "NCBIOBJ.LIB" DLL'ization
*
* Revision 5.2  1997/03/17 23:44:39  kans
* added whole_word parameter to FindInEntity and FindInEntityX, and protected
* against multiple ObjMgrAlsoSelects on a single itemID
*
* Revision 5.1  1996/09/06  20:20:41  kans
* keeps going even if ObjMgrTypeFind returns NULL (e.g., on OBJ_BIOSEQ_SEG),
* and adds a case_counts parameter for case sensitive/insensitive searches.
*
* Revision 5.0  1996/05/28  13:23:23  ostell
* Set to revision 5.0
*
* Revision 1.7  1996/02/28  04:53:06  ostell
* fix to prevernt recursion on substring replaces
*
* Revision 1.6  1996/02/26  20:24:05  kans
* replace needs MemCopy instead of StringMove (JO), and set dirty flag
*
* Revision 1.5  1996/01/03  23:06:32  ostell
* support for longer replaces, controlled updating
*
* Revision 1.3  1996/01/02  18:40:07  ostell
* simplified code.
*
* Revision 1.2  1996/01/01  00:05:14  kans
* replaced StringStr with StringISearch to ignore case
*
* Revision 1.1  1995/12/31  18:13:14  kans
* Initial revision
*
* Revision 1.1.1.1  1995/10/19 18:42:10  sad
* Initial version
*
*/

#include <ncbi.h>
#include <objall.h>
#include <objfdef.h>
#include <objsub.h>
#include <gather.h>
#include <sqnutils.h>
#include <subutil.h>
#include <findrepl.h>

/* callback type for search/replace functions */

typedef void (*FindReplFunc) (CharPtr PNTR strp, Pointer fspdata);

/* internal data structure */

typedef struct findstruct {
  Uint2         entityID;
  FindReplFunc  func;
  FindReplProc  callback;
  Pointer       userdata;

  CharPtr       find_string;
  CharPtr       replace_string;
  Boolean       case_counts;
  Boolean       whole_word;
  Int4          findLen;
  Int4          replaceLen;

  Boolean       select_item;
  Int2          send_update;
  Boolean       did_find;
  Boolean       did_replace;
  Boolean       dirty;

  Boolean       descFilter [SEQDESCR_MAX];
  Boolean       featFilter [FEATDEF_MAX];
  Boolean       seqidFilter [NUM_SEQID];

  int           d [256];
  TextFsaPtr    fsa;
} FindStruct, PNTR FindStructPtr;

#define FINDREPL_BUFFER_MAX  1000000

/* BOYER-MOORE SEARCH FUNCTIONS */

/* StringSearch and StringISearch use the Boyer-Moore algorithm, as described
   in Niklaus Wirth, Algorithms and Data Structures, Prentice- Hall, Inc.,
   Englewood Cliffs, NJ., 1986, p. 69.  The original had an error, where
   UNTIL (j < 0) OR (p[j] # s[i]) should be UNTIL (j < 0) OR (p[j] # s[k]). */

static CharPtr FindSubString (
  CharPtr str,
  CharPtr sub,
  Boolean case_counts,
  size_t strLen,
  size_t subLen,
  int *d
)

{
  int  ch;
  int  i;
  int  j;
  int  k;

  i = subLen;
  do {
    j = subLen;
    k = i;
    do {
      k--;
      j--;
    } while (j >= 0 &&
             (case_counts ? sub [j] : TO_UPPER (sub [j])) ==
             (case_counts ? str [k] : TO_UPPER (str [k])));
    if (j >= 0) {
      ch = (int) (case_counts ? str [i - 1] : TO_UPPER (str [i - 1]));
      if (ch >= 0 && ch <= 255) {
        i += d [ch];
      } else {
        i++;
      }
    }
  } while (j >= 0 && i <= (int) strLen);

  if (j < 0) {
    i -= subLen;
    return (CharPtr) (str + i);
  }

  return NULL;
}

/* passed subLen and d array to avoid repeated initialization
   of the Boyer-Moore displacement table */

static CharPtr SearchForString (
  CharPtr str,
  CharPtr sub,
  Boolean case_counts,
  Boolean whole_word,
  size_t subLen,
  int * d
)

{
  CharPtr  ptr = NULL;
  size_t   strLen;
  CharPtr  tmp;
  Boolean  keep_looking = TRUE;

  if (str == NULL || *str == '\0') return NULL;
  strLen = StringLen (str);
  if (subLen > strLen) return NULL;

  ptr = FindSubString (str, sub, case_counts, strLen, subLen, d);
  if (ptr == NULL) return NULL;

  if (! whole_word) return ptr;

  while (keep_looking && ptr != NULL) {
    keep_looking = FALSE;
    if (ptr > str) {
      tmp = ptr - 1;
      if (! IS_WHITESP (*tmp)) {
        keep_looking = TRUE;
      }
    }
    if (! keep_looking) {
      tmp = ptr + StringLen (sub);
      if (*tmp != '\0' && (! IS_WHITESP (*tmp)) && (! ispunct (*tmp))) {
        keep_looking = TRUE;
      }
    }
    if (keep_looking) {
      ptr = FindSubString (ptr + subLen, sub, case_counts, strLen, subLen, d);
    }
  }

  return ptr;
}

static void BoyerMooreFindString (
  CharPtr PNTR strp,
  Pointer userdata
)

{
  FindStructPtr  fsp;
  CharPtr        searchString;

  if (strp == NULL || userdata == NULL) return;
  fsp = (FindStructPtr) userdata;

  searchString = *strp;
  if (SearchForString (searchString, fsp->find_string, fsp->case_counts,
                       fsp->whole_word, fsp->findLen, fsp->d) != NULL) {
    fsp->did_find = TRUE;
  }
}

static void BoyerMooreReplaceString (
  CharPtr PNTR strp,
  Pointer userdata
)

{
  Int4           buffSize;
  FindStructPtr  fsp;
  Int4           searchLen;
  CharPtr        searchString;
  CharPtr        substringPtr;
  Boolean        wasChanged;
  CharPtr        workingBuffer;

  if (strp == NULL || userdata == NULL) return;
  fsp = (FindStructPtr) userdata;

  searchString = *strp;
  searchLen  = StringLen (searchString);

  wasChanged = FALSE;

  /*------------------------------------------------*/
  /* Make a guess of how big a working buffer we'll */
  /* need based on a worst case scenario.           */
  /*                                                */
  /*   A = Max occurrences of find string =         */
  /*               searchLen / findLen              */
  /*                                                */
  /*   B = Size increase for each replacement =     */
  /*               replaceLen - findLen             */
  /*                                                */
  /*   Maximum resultant string size =              */
  /*               searchLen + (A * B)              */
  /*                                                */
  /*------------------------------------------------*/

  if (fsp->replaceLen > fsp->findLen) {
      buffSize = searchLen + ((searchLen/fsp->findLen) * (fsp->replaceLen - fsp->findLen));
      if (buffSize > FINDREPL_BUFFER_MAX) {
        buffSize = FINDREPL_BUFFER_MAX;
      }
  } else {
    buffSize = searchLen;
  }

  workingBuffer = (CharPtr) MemNew (buffSize + 2);
  if (workingBuffer == NULL) return;

  workingBuffer[0] = '\0';

  /*----------------------------------------*/
  /* Create a new string with all instances */
  /* of the find string replaced by the     */
  /* replace string.                        */
  /*----------------------------------------*/

  while ((substringPtr = SearchForString (searchString, fsp->find_string,
          fsp->case_counts, fsp->whole_word, fsp->findLen, fsp->d)) != NULL) {
    wasChanged = TRUE;
    substringPtr [0] = '\0';

    if (StringLen (workingBuffer) + StringLen (searchString) > buffSize) return;

    StringCat (workingBuffer, searchString);
    StringCat (workingBuffer, fsp->replace_string);
    substringPtr [0] = 'x';
    searchString = substringPtr + fsp->findLen;
  }

  if (searchString != NULL) {
    StringCat (workingBuffer, searchString);
  }

  /*-------------------------------------*/
  /* If any replacements were made, then */
  /* swap in the new string for the old. */
  /*-------------------------------------*/

  if (wasChanged) {
    MemFree (*strp);
    (*strp) = workingBuffer;

    fsp->did_replace = TRUE;
    fsp->dirty = TRUE;
  } else {
    MemFree (workingBuffer);
  }
}

/* FINITE-STATE AUTOMATON SEARCH FUNCTION */

static void FSAFindStrings (
  CharPtr PNTR strp,
  Pointer userdata
)

{
  Char           ch;
  FindStructPtr  fsp;
  CharPtr        ptr;
  CharPtr        searchString;
  Int2           state;
  ValNodePtr     matches;

  if (strp == NULL || userdata == NULL) return;
  fsp = (FindStructPtr) userdata;

  searchString = *strp;
  if (searchString == NULL) return;

  state = 0;
  ptr = searchString;
  ch = *ptr;

  while (ch != '\0') {
    matches = NULL;
    state = TextFsaNext (fsp->fsa, state, ch, &matches);
    if (matches != NULL) {
      fsp->did_find = TRUE;
      return;
    }
    ptr++;
    ch = *ptr;
  }
}

/* MASTER SEARCH FUNCTION CALLS DESIGNATED FUNC CALLBACK */

/*=======================================================================*/
/*                                                                       */
/* FindReplString () - Does a search and replace in a given string.      */
/*                                                                       */
/*    Main Parameters:                                                   */
/*                                                                       */
/*         strp : The string to operate on. Passed as a pointer to       */
/*                a string so that it can be replaced by the             */
/*                resulting string.                                      */
/*                                                                       */
/*         fsp->find_string : The substring that is being replaced       */
/*                            in strp.                                   */
/*                                                                       */
/*         fsp->replace_string : The substring that is replacing         */
/*                               find_string in strp.                    */
/*                                                                       */
/*=======================================================================*/

static void FindReplString (
  CharPtr PNTR strp,
  FindStructPtr fsp
)

{
  if (strp == NULL || fsp == NULL || fsp->func == NULL) return;

  fsp->func (strp, (Pointer) fsp);
}

/*=======================================================================*/
/*                                                                       */
/*  FindReplStringList()                                                 */
/*                                                                       */
/*=======================================================================*/

static void FindReplStringList (
  ValNodePtr vnp,
  FindStructPtr fsp
)

{
  while (vnp != NULL) {
    FindReplString ((CharPtr PNTR) &(vnp->data.ptrvalue), fsp);
    vnp = vnp->next;
  }
}

/*=======================================================================*/
/*                                                                       */
/*  FindReplDbxrefs() -                                                  */
/*                                                                       */
/*=======================================================================*/

static void FindReplDbxrefs (
  ValNodePtr vnp,
  FindStructPtr fsp
)

{
  DbtagPtr     dbt;
  ObjectIdPtr  oip;

  while (vnp != NULL) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL) {
      FindReplString (&(dbt->db), fsp);
      oip = dbt->tag;
      if (oip != NULL && oip->str != NULL) {
        FindReplString (&(oip->str), fsp);
      }
    }
    vnp = vnp->next;
  }
}

/*=======================================================================*/
/*                                                                       */
/*  FindReplAffil() -                                                    */
/*                                                                       */
/*=======================================================================*/

static void FindReplAffil (
  AffilPtr pAffil,
  FindStructPtr fsp
)

{
  if (pAffil == NULL) return;

  if (pAffil->choice == 1) {
    FindReplString (&(pAffil->affil)      , fsp);
  } else {
    FindReplString (&(pAffil->affil)      , fsp);
    FindReplString (&(pAffil->div)        , fsp);
    FindReplString (&(pAffil->city)       , fsp);
    FindReplString (&(pAffil->sub)        , fsp);
    FindReplString (&(pAffil->country)    , fsp);
    FindReplString (&(pAffil->street)     , fsp);
    FindReplString (&(pAffil->email)      , fsp);
    FindReplString (&(pAffil->fax)        , fsp);
    FindReplString (&(pAffil->phone)      , fsp);
    FindReplString (&(pAffil->postal_code), fsp);
  }
}

/*=======================================================================*/
/*                                                                       */
/*  FindReplAuthor() -                                                   */
/*                                                                       */
/*=======================================================================*/

#define NAMESTD_LAST     0
#define NAMESTD_FIRST    1
#define NAMESTD_MIDDLE   2
#define NAMESTD_FULL     3
#define NAMESTD_INITIALS 4
#define NAMESTD_SUFFIX   5
#define NAMESTD_TITLE    6

#define PID_NOTSET 0
#define PID_DBTAG  1
#define PID_NAME   2
#define PID_ML     3
#define PID_STR    4

static void FindReplAuthor (
  AuthorPtr pAuthor,
  FindStructPtr fsp
)

{
  NameStdPtr pNameStandard;
  CharPtr    pNameStr;
  ValNodePtr pDbxref;

  if (pAuthor == NULL) return;

  FindReplAffil (pAuthor->affil, fsp);
  
  switch (pAuthor->name->choice) {
    case PID_NOTSET :
      break;
    case PID_DBTAG :
      pDbxref = pAuthor->name->data;
      FindReplDbxrefs (pDbxref, fsp);
      break;
    case PID_NAME :
      pNameStandard = pAuthor->name->data;
      if (pNameStandard != NULL) {
        FindReplString (&(pNameStandard->names [NAMESTD_LAST])    , fsp);
        FindReplString (&(pNameStandard->names [NAMESTD_FIRST])   , fsp);
        FindReplString (&(pNameStandard->names [NAMESTD_MIDDLE])  , fsp);
        FindReplString (&(pNameStandard->names [NAMESTD_FULL])    , fsp);
        FindReplString (&(pNameStandard->names [NAMESTD_INITIALS]), fsp);
        FindReplString (&(pNameStandard->names [NAMESTD_SUFFIX])  , fsp);
        FindReplString (&(pNameStandard->names [NAMESTD_TITLE])   , fsp);
      }
      break;
    case PID_ML :
    case PID_STR :
      pNameStr = pAuthor->name->data;
      FindReplString (&pNameStr, fsp);
      break;
    default:
      break;
    }
}

/*=======================================================================*/
/*                                                                       */
/*  FindReplAuthList() -                                                 */
/*                                                                       */
/*=======================================================================*/

#define AUTHLIST_STRUCTURED 1
#define AUTHLIST_ML         2
#define AUTHLIST_STRING     3
  
static void FindReplAuthlist (
  AuthListPtr alp,
  FindStructPtr fsp
)

{
  ValNodePtr vnpNames;
  CharPtr    szAuthor;
  AuthorPtr  pAuthor;

  if (alp == NULL) return;

  FindReplAffil (alp->affil, fsp);
  vnpNames = alp->names;
  while (vnpNames != NULL) {
    if (alp->choice == AUTHLIST_STRUCTURED) {
      pAuthor = (AuthorPtr) vnpNames->data.ptrvalue;
      if (pAuthor != NULL) {
        FindReplAuthor (pAuthor, fsp);
      }
    } else {
      szAuthor = (CharPtr) vnpNames->data.ptrvalue;
      if (szAuthor != NULL) {
        FindReplString (&szAuthor, fsp);
        vnpNames->data.ptrvalue = szAuthor;
      }
    }
    vnpNames = vnpNames->next;
  }
}

/*=======================================================================*/
/*                                                                       */
/*  FindReplCitRetract() -                                               */
/*                                                                       */
/*=======================================================================*/

static void FindReplCitRetract (
  CitRetractPtr pCitRetract,
  FindStructPtr fsp
)

{
  if (pCitRetract == NULL) return;

  FindReplString (&(pCitRetract->exp), fsp);
}

/*=======================================================================*/
/*                                                                       */
/*  FindReplImprint() -                                                  */
/*                                                                       */
/*=======================================================================*/

static void FindReplImprint (
  ImprintPtr pImprint,
  FindStructPtr fsp
)

{
  if (pImprint == NULL) return;

  FindReplString (&(pImprint->volume)   , fsp);
  FindReplString (&(pImprint->issue)    , fsp);
  FindReplString (&(pImprint->pages)    , fsp);
  FindReplString (&(pImprint->section)  , fsp);
  FindReplString (&(pImprint->part_sup) , fsp);
  FindReplString (&(pImprint->language) , fsp);
  FindReplString (&(pImprint->part_supi), fsp);

  FindReplAffil (pImprint->pub, fsp);
  FindReplCitRetract (pImprint->retract, fsp);
}

/*=======================================================================*/
/*                                                                       */
/*  FindReplCitBook() -                                                  */
/*                                                                       */
/*=======================================================================*/

static void FindReplCitBook (
  CitBookPtr pCitBook,
  FindStructPtr fsp
)

{
  AffilPtr    afp;
  CharPtr     tmpStr;
  ValNodePtr  vnp;

  if (pCitBook == NULL) return;

  FindReplStringList (pCitBook->title, fsp);
  FindReplImprint (pCitBook->imp, fsp);
  FindReplAuthlist (pCitBook->authors, fsp);
  FindReplStringList (pCitBook->title, fsp);
  FindReplStringList (pCitBook->coll, fsp);

  if (pCitBook->othertype == 1) {
    for (vnp = (ValNodePtr) pCitBook->otherdata; vnp != NULL; vnp = vnp->next) {
      switch (vnp->choice) {
        case 1 :
          FindReplString ((CharPtr PNTR) &(vnp->data.ptrvalue), fsp);
          break;
        case 3 :
          afp = (AffilPtr) vnp->data.ptrvalue;
          FindReplAffil (afp, fsp);
          break;
        default :
          break;
      }
    }
  } else if (pCitBook->othertype == 2) {
    tmpStr = (CharPtr) pCitBook->otherdata;
    FindReplString (&tmpStr, fsp);
    pCitBook->otherdata = tmpStr;
  }
}

static void FindReplCitArt (
  CitArtPtr pCitArt,
  FindStructPtr fsp
)

{
  CitBookPtr  pCitBook;
  CitJourPtr  pCitJournal;

  if (pCitArt == NULL) return;

  FindReplAuthlist (pCitArt->authors, fsp);
  if (pCitArt->fromptr != NULL) {
    switch (pCitArt->from) {
    case 1 :
      pCitJournal = (CitJourPtr) pCitArt->fromptr;
      FindReplStringList (pCitArt->title, fsp);
      FindReplImprint (pCitJournal->imp, fsp);
      break;
    case 2 :
    case 3 :
      pCitBook = (CitBookPtr) pCitArt->fromptr;
      FindReplCitBook (pCitBook, fsp);
      break;
    default :
      break;
    }
  }
}

/*=======================================================================*/
/*                                                                       */
/*  FindReplMedlineEntry() -                                             */
/*                                                                       */
/*=======================================================================*/

static void FindReplMedlineEntry (
  MedlineEntryPtr pMedlineEntry,
  FindStructPtr fsp
)

{
  MedlineFieldPtr pField;
  MedlineMeshPtr  pMesh;
  MedlineRnPtr    pRn;
  CharPtr         tmpStr;

  if (pMedlineEntry == NULL) return;

  FindReplCitArt(pMedlineEntry->cit, fsp);
  FindReplString (&(pMedlineEntry->abstract), fsp);

  pRn = pMedlineEntry->substance;
  while (pRn != NULL) {
    FindReplString (&(pRn->cit), fsp);
    FindReplString (&(pRn->name), fsp);
    pRn = pRn->next;
  }

  pMesh = pMedlineEntry->mesh;
  while (pMesh != NULL) {
    FindReplString (&(pMesh->term), fsp);
    pMesh = pMesh->next;
  }

  if (pMedlineEntry->xref != NULL) {
    tmpStr = (CharPtr) pMedlineEntry->xref->data.ptrvalue;
    FindReplString (&tmpStr, fsp);
    pMedlineEntry->xref->data.ptrvalue = tmpStr;
  }

  if (pMedlineEntry->idnum != NULL) {
    tmpStr = (CharPtr) pMedlineEntry->idnum->data.ptrvalue;
    FindReplString (&tmpStr, fsp);
    pMedlineEntry->idnum->data.ptrvalue = tmpStr;
  }

  if (pMedlineEntry->pub_type != NULL) {
    tmpStr = (CharPtr) pMedlineEntry->pub_type->data.ptrvalue;
    FindReplString (&tmpStr, fsp);
    pMedlineEntry->pub_type->data.ptrvalue = tmpStr;
  }

  if (pMedlineEntry->gene != NULL) {
    tmpStr = (CharPtr) pMedlineEntry->gene->data.ptrvalue;
    FindReplString (&tmpStr, fsp);
    pMedlineEntry->gene->data.ptrvalue = tmpStr;
  }

  pField = pMedlineEntry->mlfield;
  while (pField != NULL) {
    FindReplString (&(pField->str), fsp);
    pField = pField->next;
  }
}

/*=======================================================================*/
/*                                                                       */
/*  FindReplPub() -                                                      */
/*                                                                       */
/*=======================================================================*/

static void FindReplPub (
  ValNodePtr vnp,
  FindStructPtr fsp
)

{
  CitArtPtr       cap;
  CitBookPtr      cbp;
  CitGenPtr       cgp;
  CitJourPtr      cjp;
  CitPatPtr       cpp;
  ValNodePtr      cpvnp;
  CitSubPtr       csp;
  IdPatPtr        ipp;
  MedlineEntryPtr mep;
  CharPtr         tmpStr;
  ValNodePtr      pub;

  if (vnp == NULL) return;

  /* check for numerical pub types, NULL ptrvalue */

  switch (vnp->choice) {
    case PUB_PMid :
    case PUB_Muid :
      return;
    default :
      break;
  }
  if (vnp->data.ptrvalue == NULL) return;

  switch (vnp->choice) {
    case PUB_Gen :
      cgp = (CitGenPtr) vnp->data.ptrvalue;
      FindReplAuthlist (cgp->authors, fsp);
      FindReplString (&(cgp->cit), fsp);
      FindReplString (&(cgp->volume), fsp);
      FindReplString (&(cgp->issue), fsp);
      FindReplString (&(cgp->pages), fsp);
      FindReplString (&(cgp->title), fsp);
      if (cgp->journal != NULL) {
        tmpStr = (CharPtr) cgp->journal->data.ptrvalue;
        FindReplString (&tmpStr, fsp);
        cgp->journal->data.ptrvalue = tmpStr;
      }
      break;
    case PUB_Sub :
      csp = (CitSubPtr) vnp->data.ptrvalue;
      FindReplAuthlist (csp->authors, fsp);
      FindReplString (&(csp->descr), fsp);
      break;
    case PUB_Medline :
      mep = (MedlineEntryPtr) vnp->data.ptrvalue;
      FindReplMedlineEntry(mep, fsp);
      break;
    case PUB_Article :
      cap = (CitArtPtr) vnp->data.ptrvalue;
      FindReplCitArt(cap,fsp);
      break;
    case PUB_Journal :
      cjp = (CitJourPtr) vnp->data.ptrvalue;
      if (cjp->title != NULL) {
        tmpStr = (CharPtr) cjp->title->data.ptrvalue;
        FindReplString (&tmpStr, fsp);
        cjp->title->data.ptrvalue = tmpStr;
      }
      FindReplImprint (cjp->imp, fsp);
      break;
    case PUB_Book :
      cbp = (CitBookPtr) vnp->data.ptrvalue;
      FindReplCitBook (cbp, fsp);
      break;
    case PUB_Proc :
      cbp = (CitBookPtr) vnp->data.ptrvalue;
      cpvnp = cbp->otherdata;
      while (cpvnp != NULL) {
        if (cpvnp->choice == 1) {
          tmpStr = (CharPtr) cpvnp->data.ptrvalue;
          FindReplString (&tmpStr, fsp);
          cpvnp->data.ptrvalue = tmpStr;
        } else if (cpvnp->choice == 3) {
          FindReplAffil((AffilPtr) cpvnp->data.ptrvalue, fsp);
        }
        cpvnp = cpvnp->next;
      }
      break;
    case PUB_Patent :
      cpp = (CitPatPtr) vnp->data.ptrvalue;
      FindReplAuthlist (cpp->authors, fsp);
      FindReplAuthlist (cpp->applicants, fsp);
      FindReplAuthlist (cpp->assignees, fsp);
      FindReplString (&(cpp->country), fsp);
      FindReplString (&(cpp->doc_type), fsp);
      FindReplString (&(cpp->title), fsp);
      FindReplString (&(cpp->number), fsp);
      FindReplString (&(cpp->app_number), fsp);
      FindReplString (&(cpp->abstract), fsp);
      break;
    case PUB_Pat_id :
      ipp = (IdPatPtr) vnp->data.ptrvalue;
      FindReplString (&(ipp->country), fsp);
      FindReplString (&(ipp->number), fsp);
      FindReplString (&(ipp->app_number), fsp);
      FindReplString (&(ipp->doc_type), fsp);
      break;
    case PUB_Man :
      cbp = (CitBookPtr) vnp->data.ptrvalue;
      FindReplCitBook (cbp, fsp);
      break;
    case PUB_Equiv :
      /* recursive */
      for (pub = vnp->data.ptrvalue; pub != NULL; pub = pub->next) {
        FindReplPub (pub, fsp);
      }
    default:
      break;
    }
}

/*=======================================================================*/
/*                                                                       */
/*  FindReplPubDesc() -                                                  */
/*                                                                       */
/*=======================================================================*/

static void FindReplPubdesc (
  PubdescPtr pdp,
  FindStructPtr fsp
)

{
  ValNodePtr  vnp;

  if (pdp == NULL) return;

  FindReplString (&(pdp->comment), fsp);

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    FindReplPub (vnp, fsp);
  }
}

static void RemoveTaxRef (OrgRefPtr orp)
{
  ValNodePtr      vnp, next;
  ValNodePtr PNTR prev;
  DbtagPtr        dbt;

  vnp = orp->db;
  if (vnp == NULL) return;
  prev = (ValNodePtr PNTR) &(orp->db);
  while (vnp != NULL) {
    next = vnp->next;
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL && StringICmp ((CharPtr) dbt->db, "taxon") == 0) {
      *prev = vnp->next;
      vnp->next = NULL;
      DbtagFree (dbt);
      ValNodeFree (vnp);
    } else {
      prev = (ValNodePtr PNTR) &(vnp->next);
    }
    vnp = next;
  }
}

static void FindReplBioSource (
  BioSourcePtr biop,
  OrgRefPtr orp,
  FindStructPtr fsp
)

{
  OrgModPtr      omp;
  OrgNamePtr     onp;
  SubSourcePtr   ssp;
  CharPtr        old_taxname;

  if (biop != NULL) {
    orp = biop->org;
    for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
      if (ssp->subtype != SUBSRC_germline &&
          ssp->subtype != SUBSRC_rearranged &&
          ssp->subtype != SUBSRC_transgenic &&
          ssp->subtype != SUBSRC_environmental_sample) {
        FindReplString (&(ssp->name), fsp);
      }
      FindReplString (&(ssp->attrib), fsp);
    }
  }
  if (orp != NULL) {
    old_taxname = StringSave (orp->taxname);
    FindReplString (&(orp->taxname), fsp);
    if (StringCmp (old_taxname, orp->taxname) != 0)
    {
      RemoveTaxRef (orp);
      orp->common = MemFree (orp->common);
    }
    MemFree (old_taxname);
    FindReplString (&(orp->common), fsp);
    FindReplStringList (orp->mod, fsp);
    FindReplStringList (orp->syn, fsp);
    FindReplDbxrefs (orp->db, fsp);
    onp = orp->orgname;
    while (onp != NULL) {
      FindReplString (&(onp->attrib), fsp);
      FindReplString (&(onp->lineage), fsp);
      FindReplString (&(onp->div), fsp);
      for (omp = onp->mod; omp != NULL; omp = omp->next) {
        FindReplString (&(omp->subname), fsp);
        FindReplString (&(omp->attrib), fsp);
      }
      onp = onp->next;
    }
  }
}

/*=======================================================================*/
/*                                                                       */
/*  FindReplPatentSeqId() -                                              */
/*                                                                       */
/*=======================================================================*/

static void FindReplPatentSeqId (
  PatentSeqIdPtr pPatentSeqId,
  FindStructPtr fsp
)

{
  if (pPatentSeqId == NULL) return;
  if (pPatentSeqId->cit == NULL) return;

  FindReplString (&(pPatentSeqId->cit->country), fsp);
  FindReplString (&(pPatentSeqId->cit->number), fsp);
  FindReplString (&(pPatentSeqId->cit->app_number), fsp);
  FindReplString (&(pPatentSeqId->cit->doc_type), fsp);
}

/*=======================================================================*/
/*                                                                       */
/*  FindReplTextSeqId() -                                                */
/*                                                                       */
/*=======================================================================*/

static void FindReplTextSeqId (
  TextSeqIdPtr pTextSeqId, 
  FindStructPtr fsp
)

{
  if (pTextSeqId == NULL) return;

  FindReplString (&(pTextSeqId->name), fsp);
  FindReplString (&(pTextSeqId->accession), fsp);
  FindReplString (&(pTextSeqId->release), fsp);
} 

/*=======================================================================*/
/*                                                                       */
/*  FindReplGiim() -                                                     */
/*                                                                       */
/*=======================================================================*/

static void FindReplGiim (
  GiimPtr pGiim, 
  FindStructPtr fsp
)

{
  if (pGiim == NULL) return;

  FindReplString (&(pGiim->db), fsp);
  FindReplString (&(pGiim->release), fsp);
} 

/*=======================================================================*/
/*                                                                       */
/*  FindReplPDBSeqId() -                                                 */
/*                                                                       */
/*=======================================================================*/

static void FindReplPDBSeqId (
  PDBSeqIdPtr pPDBSeqId, 
  FindStructPtr fsp
)

{
  if (pPDBSeqId == NULL) return;

  FindReplString (&(pPDBSeqId->mol), fsp);
}

/*=======================================================================*/
/*                                                                       */
/*  FindReplObjectId() -                                                 */
/*                                                                       */
/*=======================================================================*/

static void FindReplObjectId (
  ObjectIdPtr pObjectId, 
  FindStructPtr fsp
)

{
  if (pObjectId == NULL) return;

  FindReplString (&(pObjectId->str), fsp);
}

/*=======================================================================*/
/*                                                                       */
/*  FindReplSeqId() -                                                    */
/*                                                                       */
/*=======================================================================*/

static void FindReplSeqId (
  SeqIdPtr sip,
  Pointer userdata
)

{
  FindStructPtr  fsp;
  Uint1          subtype;

  if (sip == NULL) return;
  fsp = (FindStructPtr) userdata;

  /* check subtype against filter */

  subtype = sip->choice;
  if (subtype >= NUM_SEQID) return;
  if (! fsp->seqidFilter [subtype]) return;

  switch (subtype) {
    case SEQID_NOT_SET :
      break;
    case SEQID_LOCAL :
      FindReplObjectId((ObjectIdPtr) sip->data.ptrvalue, fsp);
      break;
    case SEQID_GIBBSQ :
    case SEQID_GIBBMT :
      break;
    case SEQID_GIIM :
      FindReplGiim((GiimPtr) sip->data.ptrvalue, fsp);
      break;
    case SEQID_GENBANK :
    case SEQID_EMBL :
    case SEQID_PIR :
    case SEQID_SWISSPROT :
    case SEQID_OTHER :
    case SEQID_DDBJ :
    case SEQID_PRF :
    case SEQID_TPG :
    case SEQID_TPE :
    case SEQID_TPD :
    case SEQID_GPIPE :
      FindReplTextSeqId((TextSeqIdPtr) sip->data.ptrvalue, fsp);
      break;
    case SEQID_PATENT :
      FindReplPatentSeqId((PatentSeqIdPtr) sip->data.ptrvalue, fsp);
      break;
    case SEQID_GENERAL :
      FindReplDbxrefs((ValNodePtr) sip->data.ptrvalue, fsp);
      break;
    case SEQID_GI :
      break;
    case SEQID_PDB :
      FindReplPDBSeqId((PDBSeqIdPtr) sip->data.ptrvalue, fsp);
      break;
    default:
      break;
    }
}

/*=======================================================================*/
/*                                                                       */
/*  FindReplSendMessages() -                                             */
/*                                                                       */
/*=======================================================================*/

static void FindReplSendMessages (
  FindStructPtr fsp,
  Uint2 itemID,
  Uint2 itemtype
)

{
  if (fsp->send_update == UPDATE_EACH && fsp->did_replace) {
    ObjMgrSetDirtyFlag (fsp->entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, fsp->entityID, 0, 0);
  }
  if (fsp->select_item && (fsp->did_find || fsp->did_replace)) {
    ObjMgrAlsoSelect (fsp->entityID, itemID, itemtype, 0, NULL);
  }
  if (fsp->callback != NULL && (fsp->did_find || fsp->did_replace)) {
    fsp->callback (fsp->entityID, itemID, itemtype, fsp->userdata);
  }
}

/*=======================================================================*/
/*                                                                       */
/*  FindReplBioseqs() -                                                  */
/*                                                                       */
/*=======================================================================*/

static void FindReplBioseqs (
  BioseqPtr bsp,
  Pointer userdata
)

{
  FindStructPtr  fsp;
  SeqIdPtr       sip;

  if (bsp == NULL) return;

  fsp = (FindStructPtr) userdata;
  fsp->did_find = FALSE;
  fsp->did_replace = FALSE;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    FindReplSeqId (sip, userdata);
  }

  if (fsp->did_replace) {
    SeqMgrReplaceInBioseqIndex (bsp);
  }

  FindReplSendMessages (fsp, bsp->idx.itemID, bsp->idx.itemtype);
}

/*=======================================================================*/
/*                                                                       */
/*  FindReplAligns() -                                                   */
/*                                                                       */
/*=======================================================================*/

static void FindReplAligns (
  SeqAlignPtr sap,
  Pointer userdata
)

{
  DenseDiagPtr   ddp;
  DenseSegPtr    dsp;
  FindStructPtr  fsp;
  SeqIdPtr       sip;
  SeqLocPtr      slp;
  StdSegPtr      ssp;

  if (sap == NULL) return;

  fsp = (FindStructPtr) userdata;
  fsp->did_find = FALSE;
  fsp->did_replace = FALSE;

  VisitSeqIdsInSeqLoc (sap->bounds, userdata, FindReplSeqId);
  FindReplSeqId (sap->master, userdata);

  if (sap->segs == NULL) return;

  /* SAS_DISC recursively presented by visit function, so removed here */

  switch (sap->segtype) {
    case SAS_DENDIAG :
      ddp = (DenseDiagPtr) sap->segs;
      for (sip = ddp->id; sip != NULL; sip = sip->next) {
        FindReplSeqId (sip, userdata);
      }
      break;
    case SAS_DENSEG :
      dsp = (DenseSegPtr) sap->segs;
      for (sip = dsp->ids; sip != NULL; sip = sip->next) {
        FindReplSeqId (sip, userdata);
      }
      break;
    case SAS_STD :
      for (ssp = (StdSegPtr) sap->segs; ssp != NULL; ssp = ssp->next) {
        for (sip = ssp->ids; sip != NULL; sip = sip->next) {
          FindReplSeqId (sip, userdata);
        }
        for (slp = ssp->loc; slp != NULL; slp = slp->next) {
          VisitSeqIdsInSeqLoc (slp, userdata, FindReplSeqId);
        }
      }
      break;
    default :
      break;
  }

  FindReplSendMessages (fsp, sap->idx.itemID, sap->idx.itemtype);
}

/*=======================================================================*/
/*                                                                       */
/*  FindReplGraphs() -                                                   */
/*                                                                       */
/*=======================================================================*/

static void FindReplGraphs (
  SeqGraphPtr sgp,
  Pointer userdata
)

{
  FindStructPtr  fsp;

  if (sgp == NULL) return;

  fsp = (FindStructPtr) userdata;
  fsp->did_find = FALSE;
  fsp->did_replace = FALSE;

  VisitSeqIdsInSeqLoc (sgp->loc, userdata, FindReplSeqId);

  FindReplSendMessages (fsp, sgp->idx.itemID, sgp->idx.itemtype);
}

/*=======================================================================*/
/*                                                                       */
/*  FindReplFeats() -                                                    */
/*                                                                       */
/*=======================================================================*/

static void FindReplFeats (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  BioSourcePtr   biop;
  CodeBreakPtr   cbp;
  CdRegionPtr    crp;
  FindStructPtr  fsp;
  GBQualPtr      gbq;
  GeneRefPtr     grp;
  ImpFeatPtr     ifp;
  OrgRefPtr      orp = NULL;
  PubdescPtr     pdp;
  ProtRefPtr     prp;
  RnaRefPtr      rrp;
  Uint1          subtype;
  tRNAPtr        trp;

  if (sfp == NULL) return;

  fsp = (FindStructPtr) userdata;
  fsp->did_find = FALSE;
  fsp->did_replace = FALSE;

  /* change seqids on location and product */

  VisitSeqIdsInSeqLoc (sfp->location, userdata, FindReplSeqId);
  VisitSeqIdsInSeqLoc (sfp->product, userdata, FindReplSeqId);

  /* check subtype against filter */

  subtype = sfp->idx.subtype;
  if (subtype >= FEATDEF_MAX) return;
  if (! fsp->featFilter [subtype]) return;

  /* common fields */

  FindReplString (&(sfp->comment), fsp);
  FindReplString (&(sfp->title), fsp);
  FindReplString (&(sfp->except_text), fsp);

  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    FindReplString (&(gbq->qual), fsp);
    FindReplString (&(gbq->val), fsp);
  }

  FindReplDbxrefs (sfp->dbxref, fsp);

  /* check for numerical features, NULL ptrvalue */

  switch (sfp->data.choice) {
    case SEQFEAT_BOND :
    case SEQFEAT_SITE :
    case SEQFEAT_PSEC_STR :
    case SEQFEAT_COMMENT:
      return;
    default :
      break;
  }
  if (sfp->data.value.ptrvalue == NULL) return;

  /* feature-specific fields */

  switch (sfp->data.choice) {
    case SEQFEAT_GENE :
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      FindReplString (&(grp->locus), fsp);
      FindReplString (&(grp->allele), fsp);
      FindReplString (&(grp->desc), fsp);
      FindReplString (&(grp->maploc), fsp);
      FindReplString (&(grp->locus_tag), fsp);
      FindReplStringList (grp->syn, fsp);
      FindReplDbxrefs (grp->db, fsp);
      break;
    case SEQFEAT_ORG :
      orp = (OrgRefPtr) sfp->data.value.ptrvalue;
      FindReplBioSource (NULL, orp, fsp);
      break;
    case SEQFEAT_CDREGION :
      crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
        VisitSeqIdsInSeqLoc (cbp->loc, userdata, FindReplSeqId);
      }
      break;
    case SEQFEAT_PROT :
      prp = (ProtRefPtr) sfp->data.value.ptrvalue;
      FindReplString (&(prp->desc), fsp);
      FindReplStringList (prp->name, fsp);
      FindReplStringList (prp->ec, fsp);
      FindReplStringList (prp->activity, fsp);
      FindReplDbxrefs (prp->db, fsp);
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp->ext.choice == 1) {
        FindReplString ((CharPtr PNTR) &(rrp->ext.value.ptrvalue), fsp);
      } else if (rrp->ext.choice == 2) {
        trp = (tRNAPtr) rrp->ext.value.ptrvalue;
        VisitSeqIdsInSeqLoc (trp->anticodon, userdata, FindReplSeqId);
      }
      break;
    case SEQFEAT_PUB :
      pdp = (PubdescPtr) sfp->data.value.ptrvalue;
      FindReplPubdesc (pdp, fsp);
      break;
    case SEQFEAT_SEQ :
      break;
    case SEQFEAT_IMP :
      ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
      /* FindReplString (&(ifp->key), fsp); */
      FindReplString (&(ifp->loc), fsp);
      FindReplString (&(ifp->descr), fsp);
      break;
    case SEQFEAT_REGION :
      FindReplString ((CharPtr PNTR) &(sfp->data.value.ptrvalue), fsp);
      break;
    case SEQFEAT_RSITE :
      break;
    case SEQFEAT_USER :
      break;
    case SEQFEAT_TXINIT :
      break;
    case SEQFEAT_NUM :
      break;
    case SEQFEAT_NON_STD_RESIDUE :
      break;
    case SEQFEAT_HET :
      break;
    case SEQFEAT_BIOSRC :
      biop = (BioSourcePtr) sfp->data.value.ptrvalue;
      FindReplBioSource (biop, NULL, fsp);
      break;
    default :
      break;
  }

  FindReplSendMessages (fsp, sfp->idx.itemID, sfp->idx.itemtype);
}

/*=======================================================================*/
/*                                                                       */
/*  FindReplDescs() -                                                    */
/*                                                                       */
/*=======================================================================*/

static void FindReplDescs (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  BioSourcePtr   biop;
  FindStructPtr  fsp;
  GBBlockPtr     gbp;
  OrgRefPtr      orp = NULL;
  ObjValNodePtr  ovp;
  PubdescPtr     pdp;
  Uint1          subtype;

  if (sdp == NULL) return;

  fsp = (FindStructPtr) userdata;
  fsp->did_find = FALSE;
  fsp->did_replace = FALSE;

  /* check subtype against filter */

  subtype = sdp->choice;
  if (subtype >= SEQDESCR_MAX) return;
  if (! fsp->descFilter [subtype]) return;

  /* check for numerical descriptors, NULL ptrvalue */

  switch (sdp->choice) {
    case Seq_descr_mol_type :
    case Seq_descr_method :
      return;
    default :
      break;
  }
  if (sdp->data.ptrvalue == NULL) return;

  /* descriptor-specific fields */

  switch (sdp->choice) {
    case Seq_descr_modif :
      break;
    case Seq_descr_name :
      FindReplString ((CharPtr PNTR) &(sdp->data.ptrvalue), fsp);
      break;
    case Seq_descr_title :
      FindReplString ((CharPtr PNTR) &(sdp->data.ptrvalue), fsp);
      break;
    case Seq_descr_org :
      orp = (OrgRefPtr) sdp->data.ptrvalue;
      FindReplBioSource (NULL, orp, fsp);
      break;
    case Seq_descr_comment :
      FindReplString ((CharPtr PNTR) &(sdp->data.ptrvalue), fsp);
      break;
    case Seq_descr_num :
      break;
    case Seq_descr_maploc :
      break;
    case Seq_descr_pir :
      break;
    case Seq_descr_genbank :
      gbp = (GBBlockPtr) sdp->data.ptrvalue;
      FindReplStringList (gbp->extra_accessions, fsp);
      FindReplStringList (gbp->keywords, fsp);
      FindReplString (&(gbp->source), fsp);
      FindReplString (&(gbp->origin), fsp);
      FindReplString (&(gbp->date), fsp);
      FindReplString (&(gbp->div), fsp);
      FindReplString (&(gbp->taxonomy), fsp);
      break;
    case Seq_descr_pub :
      pdp = (PubdescPtr) sdp->data.ptrvalue;
      FindReplPubdesc (pdp, fsp);
      break;
    case Seq_descr_region :
      FindReplString ((CharPtr PNTR) &(sdp->data.ptrvalue), fsp);
      break;
    case Seq_descr_user :
      break;
    case Seq_descr_sp :
      break;
    case Seq_descr_dbxref :
      break;
    case Seq_descr_embl :
      break;
    case Seq_descr_create_date :
      break;
    case Seq_descr_update_date :
      break;
    case Seq_descr_prf :
      break;
    case Seq_descr_pdb :
      break;
    case Seq_descr_het :
      break;
    case Seq_descr_source :
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      FindReplBioSource (biop, NULL, fsp);
      break;
    case Seq_descr_molinfo :
      break;
    default :
      break;
  }

  if (sdp->extended != 0) {
    ovp = (ObjValNodePtr) sdp;
    FindReplSendMessages (fsp, ovp->idx.itemID, ovp->idx.itemtype);
  }
}

/*=======================================================================*/
/*                                                                       */
/*  FindReplSubmitBlock() -                                              */
/*                                                                       */
/*=======================================================================*/

static void FindReplSubmitBlock (
  SeqSubmitPtr ssp,
  FindStructPtr fsp
)

{
  ContactInfoPtr  cip;
  CitSubPtr       csp;
  SubmitBlockPtr  sub;

  if (ssp == NULL) return;
  sub = ssp->sub;
  if (sub == NULL) return;

  fsp->did_find = FALSE;
  fsp->did_replace = FALSE;

  FindReplString (&(sub->tool), fsp);
  FindReplString (&(sub->user_tag), fsp);
  FindReplString (&(sub->comment), fsp);

  cip = sub->contact;
  if (cip != NULL) {
    FindReplString (&(cip->name), fsp);
    FindReplStringList (cip->address, fsp);
    FindReplString (&(cip->phone), fsp);
    FindReplString (&(cip->fax), fsp);
    FindReplString (&(cip->email), fsp);
    FindReplString (&(cip->telex), fsp);
    FindReplObjectId (cip->owner_id, fsp);
    FindReplString (&(cip->last_name), fsp);
    FindReplString (&(cip->first_name), fsp);
    FindReplString (&(cip->middle_initial), fsp);
    FindReplAuthor (cip->contact, fsp);
  }

  csp = sub->cit;
  if (csp != NULL) {
    FindReplAuthlist (csp->authors, fsp);
    FindReplString (&(csp->descr), fsp);
  }

  FindReplSendMessages (fsp, ssp->idx.itemID, ssp->idx.itemtype);
}

/* EXTERNAL FIND-REPLACE FUNCTIONS */

/*=======================================================================*/
/*                                                                       */
/*  FindReplaceInEntity() - New find/replace function.                   */
/*                                                                       */
/*=======================================================================*/

NLM_EXTERN void FindReplaceInEntity (
  Uint2 entityID,
  CharPtr find_string,
  CharPtr replace_string,
  Boolean case_counts,
  Boolean whole_word,
  Boolean do_replace,
  Boolean select_item,
  Int2 send_update,
  BoolPtr descFilter,
  BoolPtr featFilter,
  BoolPtr seqidFilter,
  Boolean do_seqid_local,
  FindReplProc callback,
  Pointer userdata
)

{
  int            ch;
  FindStruct     fs;
  int            j;
  ObjMgrDataPtr  omdp;
  SeqEntryPtr    sep = NULL;
  SeqSubmitPtr   ssp = NULL;
  size_t         subLen;

  if (entityID == 0 || find_string == NULL 
      || (whole_word && StringHasNoText (find_string))) return;

  omdp = ObjMgrGetData (entityID);
  if (omdp != NULL) {
    switch (omdp->datatype) {
      case OBJ_SEQSUB :
        ssp = (SeqSubmitPtr) omdp->dataptr;
        if (ssp != NULL && ssp->datatype == 1) {
          sep = (SeqEntryPtr) ssp->data;
        }
        break;
      case OBJ_BIOSEQ :
        sep = (SeqEntryPtr) omdp->choice;
      case OBJ_BIOSEQSET :
        sep = (SeqEntryPtr) omdp->choice;
      default :
        break;
    }
  }
  /* sep = GetTopSeqEntryForEntityID (entityID); */
  if (sep == NULL) return;

  MemSet ((Pointer) &fs, 0, sizeof (FindStruct));

  fs.entityID = entityID;
  if (do_replace) {
    fs.func = BoyerMooreReplaceString;
  } else {
    fs.func = BoyerMooreFindString;
  }
  fs.callback = callback;
  fs.userdata = userdata;

  fs.find_string = find_string;
  fs.replace_string = replace_string;
  fs.case_counts = case_counts;
  fs.whole_word = whole_word;
  fs.findLen = StringLen (find_string);
  fs.replaceLen = StringLen (replace_string);

  fs.select_item = select_item;
  fs.send_update = send_update;
  fs.did_find = FALSE;
  fs.did_replace = FALSE;
  fs.dirty = FALSE;

  /* build Boyer-Moore displacement array in advance */

  subLen = StringLen (find_string);

  for (ch = 0; ch < 256; ch++) {
    fs.d [ch] = subLen;
  }
  for (j = 0; j < (int) (subLen - 1); j++) {
    ch = (int) (case_counts ? find_string [j] : TO_UPPER (find_string [j]));
    if (ch >= 0 && ch <= 255) {
      fs.d [ch] = subLen - j - 1;
    }
  }

  /* if desc or feat filter arrays not supplied, default to all TRUE */

  if (descFilter != NULL) {
    MemCopy ((Pointer) &fs.descFilter, (Pointer) descFilter, sizeof (fs.descFilter));
  } else {
    MemSet ((Pointer) &fs.descFilter, (int) TRUE, sizeof (fs.descFilter));
  }

  if (featFilter != NULL) {
    MemCopy ((Pointer) &fs.featFilter, (Pointer) featFilter, sizeof (fs.featFilter));
  } else {
    MemSet ((Pointer) &fs.featFilter, (int) TRUE, sizeof (fs.featFilter));
  }

  /* if seqid filter array not supplied, default to all FALSE */

  if (seqidFilter != NULL) {
    MemCopy ((Pointer) &fs.seqidFilter, (Pointer) seqidFilter, sizeof (fs.seqidFilter));
  } else if (do_seqid_local) {
    MemSet ((Pointer) &fs.seqidFilter, (int) FALSE, sizeof (fs.seqidFilter));
    fs.seqidFilter [SEQID_LOCAL] = TRUE;
  } else {
    MemSet ((Pointer) &fs.seqidFilter, (int) FALSE, sizeof (fs.seqidFilter));
  }

  /* ensure feature subtype is set in sfp->idx block */

  AssignIDsInEntity (entityID, 0, NULL);

  /* visit callbacks that find/replace specific fields */

  VisitBioseqsInSep (sep, (Pointer) &fs, FindReplBioseqs);

  VisitFeaturesInSep (sep, (Pointer) &fs, FindReplFeats);

  VisitAlignmentsInSep (sep, (Pointer) &fs, FindReplAligns);

  VisitGraphsInSep (sep, (Pointer) &fs, FindReplGraphs);

  VisitDescriptorsInSep (sep, (Pointer) &fs, FindReplDescs);

  if (ssp != NULL) {
    FindReplSubmitBlock (ssp, &fs);
  }

  /* send select message, if applicable */

  if (fs.send_update == UPDATE_ONCE && fs.dirty) {
    ObjMgrSetDirtyFlag (entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
  }
}

/*=======================================================================*/
/*                                                                       */
/*  FindStringsInEntity() - Multi-string find function.                  */
/*                                                                       */
/*=======================================================================*/

NLM_EXTERN void FindStringsInEntity (
  Uint2 entityID,
  CharPtr PNTR find_strings,
  Boolean case_counts,
  Boolean whole_word,
  Boolean select_item,
  Int2 send_update,
  BoolPtr descFilter,
  BoolPtr featFilter,
  BoolPtr seqidFilter,
  Boolean do_seqid_local,
  FindReplProc callback,
  Pointer userdata
)

{
  FindStruct     fs;
  int            j;
  ObjMgrDataPtr  omdp;
  SeqEntryPtr    sep = NULL;
  SeqSubmitPtr   ssp = NULL;

  if (entityID == 0 || find_strings == NULL) return;

  omdp = ObjMgrGetData (entityID);
  if (omdp != NULL) {
    switch (omdp->datatype) {
      case OBJ_SEQSUB :
        ssp = (SeqSubmitPtr) omdp->dataptr;
        if (ssp != NULL && ssp->datatype == 1) {
          sep = (SeqEntryPtr) ssp->data;
        }
        break;
      case OBJ_BIOSEQ :
        sep = (SeqEntryPtr) omdp->choice;
      case OBJ_BIOSEQSET :
        sep = (SeqEntryPtr) omdp->choice;
      default :
        break;
    }
  }
  /* sep = GetTopSeqEntryForEntityID (entityID); */
  if (sep == NULL) return;

  MemSet ((Pointer) &fs, 0, sizeof (FindStruct));

  fs.entityID = entityID;
  fs.func = FSAFindStrings;
  fs.callback = callback;
  fs.userdata = userdata;

  fs.find_string = NULL;
  fs.replace_string = NULL;
  fs.case_counts = case_counts;
  fs.whole_word = whole_word;
  fs.findLen = 0;
  fs.replaceLen = 0;

  fs.select_item = select_item;
  fs.send_update = send_update;
  fs.did_find = FALSE;
  fs.did_replace = FALSE;
  fs.dirty = FALSE;

  /* build finite state machine in advance */

  fs.fsa = TextFsaNew ();

  for (j = 0; find_strings [j] != NULL; j++) {
    TextFsaAdd (fs.fsa, find_strings [j]);
  }

  /* if desc or feat filter arrays not supplied, default to all TRUE */

  if (descFilter != NULL) {
    MemCopy ((Pointer) &fs.descFilter, (Pointer) descFilter, sizeof (fs.descFilter));
  } else {
    MemSet ((Pointer) &fs.descFilter, (int) TRUE, sizeof (fs.descFilter));
  }

  if (featFilter != NULL) {
    MemCopy ((Pointer) &fs.featFilter, (Pointer) featFilter, sizeof (fs.featFilter));
  } else {
    MemSet ((Pointer) &fs.featFilter, (int) TRUE, sizeof (fs.featFilter));
  }

  /* if seqid filter array not supplied, default to all FALSE */

  if (seqidFilter != NULL) {
    MemCopy ((Pointer) &fs.seqidFilter, (Pointer) seqidFilter, sizeof (fs.seqidFilter));
  } else if (do_seqid_local) {
    MemSet ((Pointer) &fs.seqidFilter, (int) FALSE, sizeof (fs.seqidFilter));
    fs.seqidFilter [SEQID_LOCAL] = TRUE;
  } else {
    MemSet ((Pointer) &fs.seqidFilter, (int) FALSE, sizeof (fs.seqidFilter));
  }

  /* ensure feature subtype is set in sfp->idx block */

  AssignIDsInEntity (entityID, 0, NULL);

  /* visit callbacks that find/replace specific fields */

  VisitBioseqsInSep (sep, (Pointer) &fs, FindReplBioseqs);

  VisitFeaturesInSep (sep, (Pointer) &fs, FindReplFeats);

  VisitAlignmentsInSep (sep, (Pointer) &fs, FindReplAligns);

  VisitGraphsInSep (sep, (Pointer) &fs, FindReplGraphs);

  VisitDescriptorsInSep (sep, (Pointer) &fs, FindReplDescs);

  if (ssp != NULL) {
    FindReplSubmitBlock (ssp, &fs);
  }

  /* clean up finite state machine */

  TextFsaFree (fs.fsa);

  /* send select message, if applicable */

  if (fs.send_update == UPDATE_ONCE && fs.dirty) {
    ObjMgrSetDirtyFlag (entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
  }
}

/*=======================================================================*/
/*                                                                       */
/*  FindReplaceString() - find/replace just one string.                  */
/*                                                                       */
/*=======================================================================*/

NLM_EXTERN void FindReplaceString (
  CharPtr PNTR strp,
  CharPtr find_string,
  CharPtr replace_string,
  Boolean case_counts,
  Boolean whole_word
)

{
  int         ch;
  FindStruct  fs;
  int         j;
  size_t      subLen;

  if (strp == NULL || StringHasNoText (find_string)) return;

  MemSet ((Pointer) &fs, 0, sizeof (FindStruct));

  fs.entityID = 0;
  fs.func = BoyerMooreReplaceString;
  fs.callback = NULL;
  fs.userdata = NULL;

  fs.find_string = find_string;
  fs.replace_string = replace_string;
  fs.case_counts = case_counts;
  fs.whole_word = whole_word;
  fs.findLen = StringLen (find_string);
  fs.replaceLen = StringLen (replace_string);

  fs.select_item = FALSE;
  fs.send_update = UPDATE_NEVER;
  fs.did_find = FALSE;
  fs.did_replace = FALSE;
  fs.dirty = FALSE;

  /* build Boyer-Moore displacement array in advance */

  subLen = StringLen (find_string);

  for (ch = 0; ch < 256; ch++) {
    fs.d [ch] = subLen;
  }
  for (j = 0; j < (int) (subLen - 1); j++) {
    ch = (int) (case_counts ? find_string [j] : TO_UPPER (find_string [j]));
    if (ch >= 0 && ch <= 255) {
      fs.d [ch] = subLen - j - 1;
    }
  }

  FindReplString (strp, &fs);
}

