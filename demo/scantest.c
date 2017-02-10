/*   scantest.c
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
* File Name:  scantest.c
*
* Author:  Kans
*
* Version Creation Date:   1/20/95
*
* $Revision: 6.3 $
*
* File Description: 
*       template for custom scans of ASN.1 release files
*       (was - scans through sequence records on the Entrez discs)
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* ==========================================================================
*/

#include <ncbi.h>
#include <objall.h>
#include <objsset.h>
#include <objsub.h>
#include <objfdef.h>
#include <sequtil.h>
#include <sqnutils.h>
#include <explore.h>

typedef struct appflags {
  Boolean  binary;
  Boolean  compressed;
  Boolean  verbose;
  FILE     *fp;
  Char     id [64];
} AppFlagData, PNTR AppFlagPtr;

static void DoOneUser (UserObjectPtr uop, Pointer userdata)

{
  AppFlagPtr   afp;
  Char         buf [128];
  ObjectIdPtr  oip;

  if (uop == NULL) return;
  afp = (AppFlagPtr) userdata;
  if (afp == NULL) return;

  buf [0] = '\0';
  if (StringDoesHaveText (uop->_class)) {
    StringCat (buf, uop->_class);
  }
  StringCat (buf, "                                ");
  buf [30] = '\0';
  fprintf (afp->fp, "%s", buf);

  buf [0] = '\0';
  oip = uop->type;
  if (oip != NULL) {
    if (StringDoesHaveText (oip->str)) {
      StringCat (buf, oip->str);
    } else if (oip->id > 0) {
      sprintf (buf, "%ld", (long) oip->id);
    }
  }
  StringCat (buf, "                                ");
  buf [30] = '\0';
  fprintf (afp->fp, "%s", buf);

  if (afp->verbose) {
    fprintf (afp->fp, "     %s", afp->id);
  }

  fprintf (afp->fp, "\n");
  fflush (afp->fp);
}

static void DoOneDescriptor (SeqDescrPtr sdp, Pointer userdata)

{
  AppFlagPtr     afp;
  UserObjectPtr  uop;

  if (sdp == NULL || sdp->choice != Seq_descr_user) return;
  afp = (AppFlagPtr) userdata;
  if (afp == NULL) return;

  uop = (UserObjectPtr) sdp->data.ptrvalue;
  if (uop == NULL) return;

  VisitUserObjectsInUop (uop, (Pointer) afp, DoOneUser);
}

static void DoOneFeature (SeqFeatPtr sfp, Pointer userdata)

{
  AppFlagPtr     afp;
  UserObjectPtr  uop;

  if (sfp == NULL) return;
  afp = (AppFlagPtr) userdata;
  if (afp == NULL) return;

  uop = sfp->ext;
  if (uop != NULL) {
    VisitUserObjectsInUop (uop, (Pointer) afp, DoOneUser);
  }

  for (uop = sfp->exts; uop != NULL; uop = uop->next) {
    VisitUserObjectsInUop (uop, (Pointer) afp, DoOneUser);
  }
}

static void DoRecord (SeqEntryPtr sep, Pointer userdata)

{
  AppFlagPtr   afp;
  BioseqPtr    fbsp;
  SeqEntryPtr  fsep;

  if (sep == NULL) return;
  afp = (AppFlagPtr) userdata;
  if (afp == NULL) return;

  fsep = FindNthBioseq (sep, 1);
  if (fsep == NULL) return;
  fbsp = (BioseqPtr) fsep->data.ptrvalue;
  if (fbsp == NULL) return;

  SeqIdWrite (fbsp->id, afp->id, PRINTID_FASTA_LONG, 64);

  VisitDescriptorsInSep (sep, (Pointer) afp, DoOneDescriptor);
  VisitFeaturesInSep (sep, (Pointer) afp, DoOneFeature);
}

static void ProcessOneRecord (
  CharPtr filename,
  Pointer userdata
)

{
  AppFlagPtr  afp;

  if (StringHasNoText (filename)) return;
  afp = (AppFlagPtr) userdata;
  if (afp == NULL) return;

  if (StringStr (filename, "gbest") != NULL ||
      StringStr (filename, "gbgss") != NULL ||
      StringStr (filename, "gbhtg") != NULL) {
    printf ("Skipping %s\n", filename);
    return;
  }

  printf ("%s\n", filename);
  fflush (stdout);

  fprintf (afp->fp, "%s\n", filename);
  fflush (afp->fp);

  ScanBioseqSetRelease (filename, afp->binary, afp->compressed, (Pointer) afp, DoRecord);

  fprintf (afp->fp, "\n");
  fflush (afp->fp);
}

#define p_argInputPath    0
#define i_argInputFile    1
#define o_argOutputFile   2
#define f_argFilter       3
#define x_argSuffix       4
#define u_argRecurse      5
#define b_argBinary       6
#define c_argCompressed   7
#define v_argVerbose      8

Args myargs [] = {
  {"Path to Files", NULL, NULL, NULL,
    TRUE, 'p', ARG_STRING, 0.0, 0, NULL},
  {"Input File Name", NULL, NULL, NULL,
    TRUE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"Output File Name", NULL, NULL, NULL,
    TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Substring Filter", NULL, NULL, NULL,
    TRUE, 'f', ARG_STRING, 0.0, 0, NULL},
  {"File Selection Suffix", ".aso", NULL, NULL,
    TRUE, 'x', ARG_STRING, 0.0, 0, NULL},
  {"Recurse", "F", NULL, NULL,
    TRUE, 'u', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Bioseq-set is Binary", "F", NULL, NULL,
    TRUE, 'b', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Bioseq-set is Compressed", "F", NULL, NULL,
    TRUE, 'c', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Verbose", "F", NULL, NULL,
    TRUE, 'v', ARG_BOOLEAN, 0.0, 0, NULL},
};

extern Int2 Main (void)

{
  AppFlagData  afd;
  Boolean      dorecurse;
  CharPtr      filter, infile, outfile, directory, suffix;

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

  if (! GetArgs ("scantest", sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  MemSet ((Pointer) &afd, 0, sizeof (AppFlagData));

  directory = (CharPtr) myargs [p_argInputPath].strvalue;
  infile = (CharPtr) myargs [i_argInputFile].strvalue;
  outfile = (CharPtr) myargs [o_argOutputFile].strvalue;
  filter = (CharPtr) myargs [f_argFilter].strvalue;
  suffix = (CharPtr) myargs [x_argSuffix].strvalue;
  dorecurse = (Boolean) myargs [u_argRecurse].intvalue;
  afd.binary = (Boolean) myargs [b_argBinary].intvalue;
  afd.compressed = (Boolean) myargs [c_argCompressed].intvalue;
  afd.verbose = (Boolean) myargs[v_argVerbose].intvalue;

  afd.fp = FileOpen (outfile, "w");
  if (afd.fp == NULL) {
    return 0;
  }

  if (StringDoesHaveText (directory)) {

    DirExplore (directory, NULL, suffix, dorecurse, ProcessOneRecord, (Pointer) &afd);

  } else if (StringDoesHaveText (infile)) {

    ProcessOneRecord (infile, &afd);
  }

  FileClose (afd.fp);

  return 0;
}
