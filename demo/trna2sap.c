/*   trna2sap.c
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
* File Name:  trna2sap.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   7/21/05
*
* $Revision: 1.5 $
*
* File Description:
*
* Modifications:
* --------------------------------------------------------------------------
* ==========================================================================
*/

#include <ncbi.h>
#include <objall.h>
#include <objsset.h>
#include <objsub.h>
#include <objfdef.h>
#include <sqnutils.h>

#define TRNA2SAPAPP_VER "1.2"

CharPtr TRNA2SAPAPPLICATION = TRNA2SAPAPP_VER;

typedef struct t2sdata {
  CharPtr  results;
  CharPtr  outfile;
  Boolean  failure;
  Boolean  allowPseudo;
  Boolean  allowUndet;
  Boolean  justFtable;
  Boolean  addCitation;
  CharPtr  name;
  CharPtr  title;
  CharPtr  comment;
  CharPtr  remark;
  Char     temp [PATH_MAX];
} T2SData, PNTR T2SPtr;

#define MAX_FIELDS  9

static void ProcessTrnaScanResult (
  FILE *ifp,
  FILE *ofp,
  Boolean allowPseudo,
  Boolean allowUndet,
  CharPtr name
)

{
  CharPtr   aa;
  CharPtr   beg;
  Char      buf [256];
  CharPtr   end;
  CharPtr   field [MAX_FIELDS];
  CharPtr   id;
  Int2      idNotSent = TRUE;
  Int2      inBody = FALSE;
  CharPtr   intronBeg;
  CharPtr   intronEnd;
  long int  intronStart;
  long int  intronStop;
  Int2      numFields = 0;
  CharPtr   ptr;
  long int  start;
  long int  stop;

  if (ifp == NULL || ofp == NULL) return;

/* line by line processing of tRNAscan-SE output table */

  while (fgets (buf, sizeof (buf), ifp) != NULL) {

    if (inBody) {
      memset (field, 0, sizeof (field));

/*
*  parse tab-delimited output line into array of fields, avoiding use of
*  strtok so that empty columns (adjacent tabs) are properly assigned to
*  field array
*/

      ptr = buf;
      for (numFields = 0; numFields < MAX_FIELDS && ptr != NULL; numFields++) {
        field [numFields] = ptr;
        ptr = strchr (ptr, '\t');
        if (ptr != NULL) {
          *ptr = '\0';
          ptr++;
        }
      }

/* interested in ID, start, stop, amino acid, and intron start and stop */

      id = field [0];
      beg = field [2];
      end = field [3];
      aa = field [4];
      intronBeg = field [6];
      intronEnd = field [7];

      if (numFields > 7 &&
          sscanf (beg, "%ld", &start) == 1 &&
          sscanf (end, "%ld", &stop) == 1 &&
          sscanf (intronBeg, "%ld", &intronStart) == 1 &&
          sscanf (intronEnd, "%ld", &intronStop) == 1 &&
          (strstr (aa, "Pseudo") == NULL || allowPseudo) &&
          (strstr (aa, "Undet") == NULL || allowUndet)) {

/* first line of output gives SeqId from FASTA definition line */

        if (idNotSent) {
          if (StringDoesHaveText (name)) {
            fprintf (ofp, ">Features %s %s\n", id, name);
          } else {
            fprintf (ofp, ">Features %s\n", id);
          }
          fflush (ofp);
          idNotSent = FALSE;
        }

/* first line of feature has start (tab) stop (tab) feature key */
/* multiple intervals would have lines of start (tab) stop */

        if (intronStart == 0 && intronStop == 0) {
          fprintf (ofp, "%ld\t%ld\ttRNA\n", (long) start, (long) stop);
          fflush (ofp);
        } else {
          fprintf (ofp, "%ld\t%ld\ttRNA\n", (long) start, (long) (intronStart - 1));
          fprintf (ofp, "%ld\t%ld\t\n", (long) (intronStop + 1), (long) stop);
          fflush (ofp);
        }

/* qualifier lines are (tab) (tab) (tab) qualifier key (tab) value */

        if (strstr (aa, "Pseudo") != NULL) {
          fprintf (ofp, "\t\t\tnote\ttRNA-Pseudo\n");
          fflush (ofp);
        } else if (strstr (aa, "Undet") != NULL) {
          fprintf (ofp, "\t\t\tproduct\tXxx\n");
          fflush (ofp);
        } else {
          fprintf (ofp, "\t\t\tproduct\t%s\n", aa);
          fflush (ofp);
        }

/* dash (formerly empty) gene qualifier to suppress /gene (e.g., if tRNA is in an intron) */

        fprintf (ofp, "\t\t\tgene\t-\n");
        fflush (ofp);
      }
    }

/* detect last line of table header, ignoring everything before data section */

    if (strstr (buf, "-----") != NULL) {
      inBody = TRUE;
    }
  }
}

static void PrepareOutputFile (
  T2SPtr tsp,
  CharPtr filename,
  CharPtr path
)

{
  Char     file [FILENAME_MAX];
  CharPtr  ptr, str;

  *path = '\0';
  str = StringRChr (filename, DIRDELIMCHR);
  if (str != NULL) {
    str++;
  } else {
    str = filename;
  }
  StringCpy (file, str);

  if (StringDoesHaveText (tsp->outfile)) {
    StringCpy (path, tsp->outfile);
  } else {
    str = StringRChr (filename, DIRDELIMCHR);
    if (str != NULL) {
      *str = '\0';
      str++;
    } else {
      str = filename;
    }
    if (StringDoesHaveText (tsp->results)) {
      StringCpy (path, tsp->results);
    } else {
      StringCpy (path, filename);
    }
    if (str != NULL) {
      ptr = StringRChr (str, '.');
      if (ptr != NULL) {
        *ptr = '\0';
      }
      StringCpy (file, str);
      if (tsp->justFtable) {
        StringCat (file, ".tbl");
      } else {
        StringCat (file, ".sap");
      }
      FileBuildPath (path, NULL, file);
    }
  }
}

static CharPtr tsePubStr = "Pubdesc ::= { \n" \
"pub { \n" \
"article { \n" \
"title { \n" \
"name \"tRNAscan-SE: a program for improved detection of transfer RNA \n" \
"genes in genomic sequence.\" } , \n" \
"authors { \n" \
"names \n" \
"std { \n" \
"{ \n" \
"name \n" \
"name { \n" \
"last \"Lowe\" , \n" \
"first \"T\" , \n" \
"initials \"T.M.\" } } , \n" \
"{ \n" \
"name \n" \
"name { \n" \
"last \"Eddy\" , \n" \
"first \"S\" , \n" \
"initials \"S.R.\" } } } , \n" \
"affil \n" \
"str \"Department of Genetics, Washington University School of \n" \
"Medicine, 660 South Euclid, Box 8232, St Louis, MO 63110, USA.\" } , \n" \
"from \n" \
"journal { \n" \
"title { \n" \
"iso-jta \"Nucleic Acids Res.\" } , \n" \
"imp { \n" \
"date \n" \
"std { \n" \
"year 1997 , \n" \
"month 3 , \n" \
"day 1 } , \n" \
"volume \"25\" , \n" \
"issue \"5\" , \n" \
"pages \"955-964\" , \n" \
"pubstatus ppublish } } } , \n" \
"pmid 9023104 } , \n" \
"reftype feats }";

static void ProcessOneRecord (
  CharPtr filename,
  Pointer userdata
)

{
  AnnotDescrPtr  adp;
  AsnIoMemPtr    aimp;
  AsnIoPtr       aip = NULL;
  Pointer        dataptr = NULL;
  Uint2          datatype, entityID = 0;
  Char           path [PATH_MAX];
  PubdescPtr     pdp;
  FILE           *ifp, *ofp = NULL;
  SeqAnnotPtr    sap;
  T2SPtr         tsp;

  if (StringHasNoText (filename)) return;
  tsp = (T2SPtr) userdata;
  if (tsp == NULL) return;

  ifp = FileOpen (filename, "r");
  if (ifp == NULL) {
    Message (MSG_POSTERR, "Failed to open input file '%s'", filename);
    tsp->failure = TRUE;
    return;
  }

  if (tsp->justFtable) {
    PrepareOutputFile (tsp, filename, path);
    ofp = FileOpen (path, "w");
    if (ofp == NULL) {
      Message (MSG_POSTERR, "Failed to open '%s'", path);
      FileClose (ifp);
      tsp->failure = TRUE;
      return;
    }

    ProcessTrnaScanResult (ifp, ofp, tsp->allowPseudo, tsp->allowUndet, tsp->name);

    FileClose (ofp);
    FileClose (ifp);

    return;
  }

  ofp = FileOpen (tsp->temp, "w");
  if (ofp == NULL) {
    Message (MSG_POSTERR, "Failed to open temporary file '%s' for output", tsp->temp);
    FileClose (ifp);
    tsp->failure = TRUE;
    return;
  }

  ProcessTrnaScanResult (ifp, ofp, tsp->allowPseudo, tsp->allowUndet, tsp->name);

  FileClose (ofp);
  FileClose (ifp);

  ifp = FileOpen (tsp->temp, "r");
  if (ifp == NULL) {
    Message (MSG_POSTERR, "Failed to open temporary file '%s' for input", tsp->temp);
    tsp->failure = TRUE;
    return;
  }

  dataptr = ReadAsnFastaOrFlatFile (ifp, &datatype, NULL, FALSE, FALSE, FALSE, FALSE);

  FileClose (ifp);

  if (datatype != OBJ_SEQANNOT) {
    Message (MSG_POSTERR, "Expected Seq-annot, got datatype %d instead", (int) datatype);
    tsp->failure = TRUE;
    return;
  }
  sap = (SeqAnnotPtr) dataptr;
  if (sap == NULL) {
    Message (MSG_POSTERR, "Seq-annot was not read");
    tsp->failure = TRUE;
    return;
  }
  if (StringDoesHaveText (tsp->title)) {
    adp = AnnotDescrNew (sap->desc);
    if (adp != NULL) {
      adp->choice = Annot_descr_title;
      adp->data.ptrvalue = StringSave (tsp->title);
    }
  }
  if (StringDoesHaveText (tsp->comment)) {
    adp = AnnotDescrNew (sap->desc);
    if (adp != NULL) {
      adp->choice = Annot_descr_comment;
      adp->data.ptrvalue = StringSave (tsp->comment);
    }
  }
  if (tsp->addCitation) {
    aimp = AsnIoMemOpen ("r", (BytePtr) tsePubStr, (Int4) StringLen (tsePubStr));
    if (aimp != NULL && aimp->aip != NULL) {
      pdp = PubdescAsnRead (aimp->aip, NULL);
      AsnIoMemClose (aimp);
      if (pdp != NULL) {
        adp = AnnotDescrNew (sap->desc);
        if (adp != NULL) {
          adp->choice = Annot_descr_pub;
          adp->data.ptrvalue = pdp;
          if (StringDoesHaveText (tsp->remark)) {
            pdp->comment = StringSave (tsp->remark);
          }
        }
      }
    }
  }

  PrepareOutputFile (tsp, filename, path);
  aip = AsnIoOpen (path, "w");

  if (aip != NULL) {
    SeqAnnotAsnWrite (sap, aip, NULL);
    AsnIoClose (aip);
  }

  ObjMgrFree (datatype, dataptr);
}

/* Args structure contains command-line arguments */

#define p_argInputPath     0
#define r_argOutputPath    1
#define i_argInputFile     2
#define o_argOutputFile    3
#define f_argFilter        4
#define x_argSuffix        5
#define s_argIgnorePseudo  6
#define u_argIgnoreUndet   7
#define j_argJustFtable    8
#define a_argAddCitation   9
#define n_argName         10
#define t_argTitle        11
#define c_argComment      12
#define m_argRemark       13


Args myargs [] = {
  {"Path to Files", NULL, NULL, NULL,
    TRUE, 'p', ARG_STRING, 0.0, 0, NULL},
  {"Path for Results", NULL, NULL, NULL,
    TRUE, 'r', ARG_STRING, 0.0, 0, NULL},
  {"Single Input File", "stdin", NULL, NULL,
    TRUE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"Single Output File", "stdout", NULL, NULL,
    TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Substring Filter", NULL, NULL, NULL,
    TRUE, 'f', ARG_STRING, 0.0, 0, NULL},
  {"File Selection Suffix", ".trna", NULL, NULL,
    TRUE, 'x', ARG_STRING, 0.0, 0, NULL},
  {"Ignore Pseudo tRNAs", "F", NULL, NULL,
    TRUE, 's', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Ignore Undetermined tRNAs", "F", NULL, NULL,
    TRUE, 'u', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Just Make Ftable", "F", NULL, NULL,
    TRUE, 'j', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Add tRNAscan-SE Citation", "F", NULL, NULL,
    TRUE, 'a', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Name", "tRNA", NULL, NULL,
    TRUE, 'n', ARG_STRING, 0.0, 0, NULL},
  {"Title", "tRNAscan-SE", NULL, NULL,
    TRUE, 't', ARG_STRING, 0.0, 0, NULL},
  {"Comment", NULL, NULL, NULL,
    TRUE, 'c', ARG_STRING, 0.0, 0, NULL},
  {"Remark", NULL, NULL, NULL,
    TRUE, 'm', ARG_STRING, 0.0, 0, NULL},
};

Int2 Main (void)

{
  Char     app [64];
  CharPtr  directory, filter, infile, outfile, results, suffix;
  Boolean  ignorePseudo, ignoreUndet;
  T2SData  tsd;
  T2SPtr   tsp;

  /* standard setup */

  ErrSetFatalLevel (SEV_MAX);
  ErrClearOptFlags (EO_SHOW_USERSTR);
  ErrSetLogfile ("stderr", ELOG_APPEND);
  UseLocalAsnloadDataAndErrMsg ();
  ErrPathReset ();

  if (! AllObjLoad ()) {
    Message (MSG_FATAL, "AllObjLoad failed");
    return 1;
  }
  if (! SubmitAsnLoad ()) {
    Message (MSG_FATAL, "SubmitAsnLoad failed");
    return 1;
  }
  if (! FeatDefSetLoad ()) {
    Message (MSG_FATAL, "FeatDefSetLoad failed");
    return 1;
  }
  if (! SeqCodeSetLoad ()) {
    Message (MSG_FATAL, "SeqCodeSetLoad failed");
    return 1;
  }
  if (! GeneticCodeTableLoad ()) {
    Message (MSG_FATAL, "GeneticCodeTableLoad failed");
    return 1;
  }

  /* process command line arguments */

  sprintf (app, "trna2sap %s", TRNA2SAPAPPLICATION);
  if (! GetArgs (app, sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  MemSet ((Pointer) &tsd, 0, sizeof (T2SData));
  tsp = &tsd;
  tsd.failure = FALSE;
  ignorePseudo = (Boolean) myargs [s_argIgnorePseudo].intvalue;
  tsd.allowPseudo = (Boolean) (! ignorePseudo);
  ignoreUndet = (Boolean) myargs [u_argIgnoreUndet].intvalue;
  tsd.allowUndet = (Boolean) (! ignoreUndet);
  tsd.justFtable = (Boolean) myargs [j_argJustFtable].intvalue;
  tsd.addCitation = (Boolean) myargs [a_argAddCitation].intvalue;
  tsd.name = (CharPtr) myargs [n_argName].strvalue;
  tsd.title = (CharPtr) myargs [t_argTitle].strvalue;
  tsd.comment = (CharPtr) myargs [c_argComment].strvalue;
  tsd.remark = (CharPtr) myargs [m_argRemark].strvalue;

  directory = (CharPtr) myargs [p_argInputPath].strvalue;
  results = (CharPtr) myargs [r_argOutputPath].strvalue;
  if (StringHasNoText (results)) {
    results = NULL;
  }
  infile = (CharPtr) myargs [i_argInputFile].strvalue;
  outfile = (CharPtr) myargs [o_argOutputFile].strvalue;
  filter = (CharPtr) myargs [f_argFilter].strvalue;
  suffix = (CharPtr) myargs [x_argSuffix].strvalue;

  TmpNam (tsd.temp);

  /* process input file */

  if (StringDoesHaveText (directory)) {

    tsd.results = results;

    DirExplore (directory, filter, suffix, TRUE, ProcessOneRecord, (Pointer) &tsd);

  } else if (StringDoesHaveText (infile) && StringDoesHaveText (outfile)) {

    tsd.outfile = outfile;

    ProcessOneRecord (infile, (Pointer) &tsd);
  }

  FileRemove (tsd.temp);

  if (tsd.failure) return 1;

  return 0;
}

