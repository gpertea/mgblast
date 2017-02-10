/*   sequin9.c
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
* File Name:  sequin9.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   4/20/99
*
* $Revision: 6.427 $
*
* File Description: 
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* ==========================================================================
*/

#include <subutil.h>
#include <explore.h>
#include <alignmgr.h>
#include <urkptpf.h>
#include <entrez.h>
#include <accentr.h>
#include <urlquery.h>
#include <vecscrn.h>
#include <vecscnapi.h>
#include <qblastapi.h>
#include <edutil.h>
#include <actutils.h>
#include <findrepl.h>
#include <rpsutil.h>
#include "sequin.h"
#include <seqpanel.h>
#include <salpanel.h>
#include <assert.h>
#include <pmfapi.h>

/*-------------------*/
/* Defined Constants */
/*-------------------*/

#define MAX_ID_LEN 41
#define FASTA_READ_OK     0
#define FASTA_READ_ERROR -1
#define FASTA_READ_DONE   1

#define CONVERTPUBS_NOT_SET 0
#define CONVERTPUBS_YES     1
#define CONVERTPUBS_NO      2

/* constants for update sequence */
#define UPDATE_CHOICE_NOT_SET        0
#define UPDATE_SEQUENCE_ONLY         1
#define UPDATE_FEATURES_ONLY         2
#define UPDATE_SEQUENCE_AND_FEATURES 3

#define UPDATE_REPLACE               1
#define UPDATE_EXTEND5               2
#define UPDATE_EXTEND3               3
#define UPDATE_PATCH                 4

enum update_dup_feat_type
{
  UPDATE_FEAT_DUP_NOT_SET = 0,
  UPDATE_FEAT_DUP_USE_NEW,
  UPDATE_FEAT_DUP_USE_OLD,
  UPDATE_FEAT_DUP_USE_BOTH,
  UPDATE_FEAT_DUP_MERGE,
  UPDATE_FEAT_DUP_REPLACE
};


/*-----------------*/
/* Data Structures */
/*-----------------*/

typedef struct {
  Char       newId[MAX_ID_LEN];
  BioseqPtr  matchingBsp;
} UpdateData, PNTR UpdateDataPtr;

typedef struct upsdata {
  FORM_MESSAGE_BLOCK
  ButtoN              accept;
  ButtoN              acceptAll;
  CharPtr             aln1;
  CharPtr             aln2;
  Int4                aln_length;
  Int2                charwidth;
  Int2                convertPubs;
  VieweR              details;
  Boolean             diffOrgs;
  SegmenT             dtpict;
  FILE                *fp;
  ValNodePtr          indels;
  Boolean             isSet;
  ButtoN              keepProteinIDs;
  PaneL               letters;
  Int2                lineheight;
  Int4                log10_aln_length;
  Int2                maxchars;
  ValNodePtr          mismatches;
  Int4                new3;
  Int4                new5;
  Int4                newa;
  BioseqPtr           newbsp;
  ButtoN              replace_all;
  GrouP               nobm;
  Int4                old3;
  Int4                old5;
  Int4                olda;
  BioseqPtr           oldbsp;
  VieweR              overview;
  SegmenT             ovpict;
  Int4                recomb1;
  Int4                recomb2;
  Boolean             revcomp;
  GrouP               rmc;
  SeqAlignPtr         salp;
  Int4                scaleX;
  CharPtr             seq1;
  CharPtr             seq2;
  GrouP               sfb;
  Int4                startmax;
  Int4                stopmax;
  Uint1               strandnew;
  Uint1               strandold;
  Boolean             useGUI;
  Boolean             do_update;
  ButtoN              add_cit_subs;
  ButtoN              update_proteins;
  Boolean             suppress_continue_msg;
  Boolean             suppress_instant_refresh;
  ButtoN              update_quality_scores_btn;
  FILE                *log_fp;
  Char                log_path [PATH_MAX];
  Boolean             data_in_log;
  ButtoN              truncate_proteins_btn;
  ButtoN              extend_proteins5_btn;
  ButtoN              extend_proteins3_btn;
  ButtoN              correct_cds_genes_btn;
  Boolean             truncate_proteins;
  Boolean             extend_proteins5;
  Boolean             extend_proteins3;
  Boolean             correct_cds_genes;
  ValNodePtr          transl_except_list;
  Int2                rmcval;      /* This is the choice selected from the
                                    * rmc radio button group.
                                    * We store it so that we can use it in
                                    * Accept All.
                                    */
  
  SeqEntryPtr         seqsubsep;   /* when we are updating from a Sequin record
                                    * that contains a set, store the set here.
                                    */
  Int4                seqsubpos;   /* when we are updating from a Sequin record
                                    * that contains a set, store the number of
                                    * Bioseqs already processed here.
                                    */
                                    
  Int2                no_aln_choice; /* what to do when updating a set and no
                                      * alignment is found.
                                      */
  
  ValNodePtr          affected_variation_features;   /* list of variation features
                                                      * affected by this update.
                                                      * by "affected" we mean that
                                                      * an insertion, deletion,
                                                      * or replacement takes place
                                                      * either inside the variation
                                                      * location, or immediately
                                                      * to the left or right of the
                                                      * location.
                                                      */
} UpsData, PNTR UpsDataPtr;

/*---------------------*/
/* Function prototypes */
/*---------------------*/

static Int2 UpdateNextBioseqInFastaSet (UpsDataPtr udp);
static void FreeUdpFields (UpsDataPtr udp);
static void TruncateCDS (SeqFeatPtr sfp, Uint1 frame, BioseqPtr pbsp);

extern void SubmitToNCBI (IteM i);
extern void QuitProc (void);
extern void NewFeaturePropagate (IteM i);
extern void FixCdsAfterPropagate (IteM i);
extern void MakeRedundantGPSwithProp (IteM i);
extern void MakeRedundantGPSnoProp (IteM i);
extern void MakeRedundantGPSjustXref (IteM i);
extern void FuseSlpJoins (IteM i);

/*-----------*/
/* Functions */
/*-----------*/

/*
 * This function is called by HandleOneNewAsnProc, CommonHandleMultBioseqs,
 * and DoReadAnythingLoop in sequin1.c
 */
extern void HandleProjectAsn (ProjectPtr proj, Uint2 entityID)

{
  Int2              db = -1;
  EntrezGlobalsPtr  egp;
  Int4              i;
  ValNodePtr        list;
  Int4              num;
  ValNodePtr        pip;
  Int4Ptr           uids;
  ValNodePtr        vnp;

  if (proj == NULL) return;
  if (! useEntrez) return;
  egp = (EntrezGlobalsPtr) GetAppProperty ("EntrezGlobals");
  if (egp == NULL) return;
  pip = proj->data;
  if (pip == NULL) return;
  list = (ValNodePtr) pip->data.ptrvalue;
  if (list == NULL) return;
  if (pip->choice >= ProjectItem_protent && pip->choice <= ProjectItem_genomeent) {
    if (egp->retrieveProjectProc == NULL) return;
    if (! EntrezIsInited ()) {
      SequinEntrezInit ("Sequin", FALSE, NULL);
    }
    egp->retrieveProjectProc (NULL, (Pointer) proj);
    Update ();
    return;
  }
  if (egp->retrieveDocsProc == NULL) return;
  switch (pip->choice) {
    case ProjectItem_pmuid :
      db = 0;
      break;
    case ProjectItem_protuid :
      db = 1;
      break;
    case ProjectItem_nucuid :
      db = 2;
      break;
    case ProjectItem_genomeuid :
      db = 4;
      break;
    case ProjectItem_structuid :
      db = 3;
      break;
    default :
      break;
  }
  if (db == -1) return;
  if (! EntrezIsInited ()) {
    SequinEntrezInit ("Sequin", FALSE, NULL);
  }
  num = ValNodeLen (list);
  uids = MemNew ((size_t) (num * sizeof (Int4)));
  if (uids == NULL) return;
  for (i = 0, vnp = list; i < 32766 && vnp != NULL; i++, vnp = vnp->next) {
    uids [i] = vnp->data.intvalue;
  }
  if (egp->retrieveDocsProc != NULL) {
    egp->retrieveDocsProc (NULL, (Int2) num, 0, uids, db);
  }
  MemFree (uids);
  Update ();
}

/* BioseqViewOrDocSumChoice allows a single callback for each analysis item */
static Int2 BioseqViewOrDocSumChoice (NewObjectPtr nop)

{
  Int2  which = 0;   /* 1 = bioseq viewer, 2 = docsum window */

  if (nop == NULL) return 0;
#ifdef WIN_MAC
  if (bioseqViewUp) {
    which = 1;
  } else if (docSumUp) {
    which = 2;
  }
#else
  if (nop->bspOK) {
    which = 1;
  } else if (nop->dsmOK) {
    which = 2;
  }
#endif
  return which;
}

/*
#define TEST_APPLE_EVENT_MESSAGING
*/

#ifndef TEST_APPLE_EVENT_MESSAGING
static void AddRestrictionSite (SeqAnnotPtr annot, PackSeqPntPtr pspp, CharPtr name)

{
  DbtagPtr     dbt;
  ObjectIdPtr  oip;
  RsiteRefPtr  rrp;
  SeqFeatPtr   sfp, lastsfp;
  SeqLocPtr    slp;

  if (annot == NULL || pspp == NULL || name == NULL) return;
  slp = ValNodeNew (NULL);
  if (slp == NULL) return;
  slp->choice = SEQLOC_PACKED_PNT;
  slp->data.ptrvalue = (Pointer) pspp;
  sfp = SeqFeatNew ();
  if (sfp == NULL) return;

  sfp->data.choice = SEQFEAT_RSITE;
  sfp->location = slp;
  dbt = DbtagNew ();
  if (dbt != NULL) {
    dbt->db = StringSave ("REBASE");
    oip = ObjectIdNew ();
    if (oip != NULL) {
      oip->str = StringSave (name);
    }
    dbt->tag = oip;
  }
  rrp = ValNodeNew (NULL);
  if (rrp != NULL) {
    rrp->choice = 2;
    rrp->data.ptrvalue = dbt;
  }
  sfp->data.value.ptrvalue = (Pointer) rrp;

  if (annot->data == NULL) {
    annot->data = (Pointer) sfp;
  } else {
    lastsfp = (SeqFeatPtr) annot->data;
    while (lastsfp->next != NULL) {
      lastsfp = lastsfp->next;
    }
    lastsfp->next = sfp;
  }
}

static void RestrictionSearchProc (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  SeqAnnotPtr    annot;
  BioseqPtr      bsp;
  ComPatPtr      cpp, cpph;
  ValNodePtr     desc;
  SeqAnnotPtr    lastannot;
  PackSeqPntPtr  pspp;
  Int4           pt;
  SeqAlignPtr    sap;
  SeqLocPtr      slp, slpn;

  if (sep == NULL) return;
  if (! IS_Bioseq (sep)) return;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL) return;
  if (! ISA_na (bsp->mol)) return;

  desc = SeqDescrNew (NULL);
  desc->choice = Annot_descr_name;
  desc->data.ptrvalue = StringSave ("cutsites");

  annot = SeqAnnotNew ();
  annot->type = 1;
  annot->desc = desc;
  annot->data = NULL;

  cpph = (ComPatPtr) mydata;
  cpp = cpph;
  while (cpp != NULL) {
    sap = PatternMatchBioseq (bsp, cpp, 0);
    slp = MatchSa2Sl (&sap);
    if (slp != NULL) {
      pspp = PackSeqPntNew ();
      pspp->id = SeqIdDup (SeqIdFindBest (bsp->id, 0));
      while (slp != NULL) {
        pt = SeqLocStart (slp);
        PackSeqPntPut (pspp, pt);
        slpn = slp->next;
        SeqLocFree (slp);
        slp = slpn;
      }
      AddRestrictionSite (annot, pspp, cpp->name);
    }
    cpp = cpp->nextpattern;
  }

  if (annot->data == NULL) {
    SeqAnnotFree (annot);
    return;
  }
  if (bsp->annot == NULL) {
    bsp->annot = annot;
  } else {
    lastannot = bsp->annot;
    while (lastannot->next != NULL) {
      lastannot = lastannot->next;
    }
    lastannot->next = annot;
  }
}
#endif

static void SimpleRsiteProc (IteM i)

{
  BaseFormPtr   bfp;
  BioseqPtr     bsp;
  NewObjectPtr  nop;
  SeqEntryPtr   sep = NULL;
  Int2          which;
#ifdef TEST_APPLE_EVENT_MESSAGING
  AsnIoPtr      aip;
  Char          tmp [PATH_MAX];
#else
  ComPatPtr     cpph;
  Char          enzyme [64];
  Int2          j;
  Char          temp [32];
  ValNodePtr    enzymes;
#endif

  nop = (NewObjectPtr) GetObjectExtra (i);
  if (nop == NULL) return;
#ifdef WIN_MAC
  bfp = (BaseFormPtr) currentFormDataPtr;
#else
  bfp = nop->bfp;
#endif
  if (bfp == NULL) return;
  which = BioseqViewOrDocSumChoice (nop);
  if (which != 1) return;
  bsp =  GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
  if (bsp != NULL) {
    sep = SeqMgrGetSeqEntryForData (bsp);
  } else {
    sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  }
  if (sep == NULL) return;

#ifdef TEST_APPLE_EVENT_MESSAGING
  TmpNam (tmp);
  aip = AsnIoOpen (tmp, "w");
  if (aip != NULL) {
    SeqEntryAsnWrite (sep, aip, NULL);
    AsnIoClose (aip);
    /* Nlm_SendOpenDocAppleEventEx (tmp, "REST", NULL, TRUE); */
    Nlm_SendOpenDocAppleEventEx (tmp, NULL, "RsiteFind", TRUE);
  }
#else
  enzymes = NULL;
  j = 1;
  sprintf (temp, "ENZ_%d", (int) j);
  while (GetAppParam ("SEQNCGIS", "ENZYMES", temp, NULL, enzyme, sizeof (enzyme) - 1)) {
    ValNodeCopyStr (&enzymes, 0, enzyme);
    j++;
    sprintf (temp, "ENZ_%d", (int) j);
  }
  if (enzymes == NULL) {
    ValNodeCopyStr (&enzymes, 0, "BamHI");
    ValNodeCopyStr (&enzymes, 0, "EcoRI");
    ValNodeCopyStr (&enzymes, 0, "HindIII");
  }
  cpph = ReadRENPattern ("ncbiren.dat", FALSE, enzymes);
  /* PalindromeCheck (cpph); */
  SeqEntryExplore (sep, (Pointer) cpph, RestrictionSearchProc);
  ComPatFree (cpph);
  ValNodeFreeData (enzymes);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();
#endif
}

static VQUEUE  vsquerylist = NULL;

static Int2 vsquerynum = 0;

/*
static void LIBCALLBACK AnnounceCallback (CharPtr requestID, CharPtr seqID, Int2 estimatedSeconds)

{
  if (StringHasNoText (requestID)) {
    requestID = "?";
  }
  if (StringHasNoText (seqID)) {
    seqID = "?";
  }
  Message (MSG_POST, "Queued rID %s, seqID %s, estimated seconds = %d",
           requestID, seqID, (int) estimatedSeconds);

  vsquerynum++;
}

static Boolean LIBCALLBACK VecScreenCallback (
  CharPtr filename,
  VoidPtr userdata,
  CharPtr requestID,
  CharPtr seqID,
  Boolean success
)

{
  if (StringHasNoText (requestID)) {
    requestID = "?";
  }
  if (StringHasNoText (seqID)) {
    seqID = "?";
  }
  if (success) {
    if (! SequinHandleNetResults (filename)) {
    }
  } else {
    Message (MSG_POST, "Failure of rID '%s', seqID %s", requestID, seqID);
  }
  return TRUE;
}

static void DoVecScreens (BioseqPtr bsp, Pointer userdata)

{
  CharPtr  service;

  if (bsp == NULL || ISA_aa (bsp->mol)) return;
  service = (CharPtr) userdata;
  VecScreenAsynchronousRequest (service, bsp, &vsquerylist, VecScreenCallback, AnnounceCallback, NULL);
}

static void SimpleVecScreenCommon (IteM i, CharPtr service)

{
  BaseFormPtr   bfp;
  BioseqPtr     bsp;
  NewObjectPtr  nop;
  SeqEntryPtr   sep = NULL;
  Int2          which;

  nop = (NewObjectPtr) GetObjectExtra (i);
  if (nop == NULL) return;
#ifdef WIN_MAC
  bfp = (BaseFormPtr) currentFormDataPtr;
#else
  bfp = nop->bfp;
#endif
  if (bfp == NULL) return;
  which = BioseqViewOrDocSumChoice (nop);
  if (which != 1) return;
  bsp =  GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
  if (bsp != NULL) {
    sep = SeqMgrGetSeqEntryForData (bsp);
    if (sep == NULL) return;
    VecScreenAsynchronousRequest (service, bsp, &vsquerylist, VecScreenCallback, AnnounceCallback, NULL);
  } else {
    sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
    if (sep == NULL) return;
    VisitBioseqsInSep (sep, (Pointer) service, DoVecScreens);
  }
}
*/

typedef struct vsdata {
  CharPtr     date;
  CharPtr     path;
  CharPtr     database;
  Boolean     changed;
  Int4        count;
  MonitorPtr  mon;
} VsData, PNTR VsDataPtr;

static CharPtr vectorStrengths [6] = {
  NULL,
  "Strong match",
  "Moderate match",
  "Weak match",
  "Suspect origin",
  "Unknown vector match"
};

static void CountVecScreens (BioseqPtr bsp, Pointer userdata)

{
  Int4Ptr  maxp;

  if (bsp == NULL || ISA_aa (bsp->mol)) return;
  maxp = (Int4Ptr) userdata;
  (*maxp)++;
}

static void DoVecScreens (BioseqPtr bsp, Pointer userdata)

{
  AnnotDescrPtr   desc;
  GeneRefPtr      grp;
  Int2            hits;
  ImpFeatPtr      ifp;
  ValNodePtr      locs = NULL;
  Char            note [128];
  SeqAnnotPtr     prevsap;
  SeqFeatPtr      prevsfp;
  SeqAnnotPtr     sap = NULL;
  SeqFeatPtr      sfp = NULL;
  SeqIdPtr        sip;
  SeqLocPtr       slp;
  Int4            start;
  Int4            stop;
  Uint1           strength;
  SeqLocPtr       tmp;
  VsDataPtr       vsp;
  ValNodePtr      vnp;
  SeqFeatXrefPtr  xref;

  if (bsp == NULL || ISA_aa (bsp->mol)) return;
  vsp = (VsDataPtr) userdata;
  if (vsp == NULL) return;
  sip = SeqIdFindBest (bsp->id, 0);
  if (sip == NULL) return;

  if (vsp->mon != NULL) {
    MonitorIntValue (vsp->mon, vsp->count);
  }
  (vsp->count)++;

  hits = VSScreenSequence (bsp, NULL, vsp->path, NULL, &locs, NULL, NULL);
  for (vnp = locs; vnp != NULL; vnp = vnp->next) {
    strength = vnp->choice;
    if (strength < 1 || strength > 4) {
      strength = 5;
    }
    slp = (SeqLocPtr) vnp->data.ptrvalue;
    if (slp == NULL) continue;
    if (slp->choice == SEQLOC_PACKED_INT) {
      for (tmp = (SeqLocPtr) slp->data.ptrvalue; tmp != NULL; tmp = tmp->next) {
        start = SeqLocStart (tmp);
        stop = SeqLocStop (tmp);
        if (start < 0 || stop < 0) continue;
        vsp->changed = TRUE;

        if (sap == NULL) {
          sap = SeqAnnotNew ();
          if (sap != NULL) {
            sap->type = 1;
            desc = AnnotDescrNew (NULL);
            if (desc != NULL) {
              desc->choice = Annot_descr_name;
              desc->data.ptrvalue = StringSave ("VecScreen");
              sap->desc = desc;
            }
          }
        }

        sfp = SeqFeatNew ();
        if (sfp != NULL) {

          /* make misc_feature for now */

          sfp->data.choice = SEQFEAT_IMP;
          ifp = ImpFeatNew ();
          if (ifp != NULL) {
            ifp->key = StringSave ("misc_feature");
          }
          AddQualifierToFeature (sfp, "standard_name", "Vector Contamination");
          AddQualifierToFeature (sfp, "phenotype", vectorStrengths [(int) strength]);

          sprintf (note, "Screened against %s using VecScreen on %s", vsp->database, vsp->date);
          sfp->comment = StringSave (note);

          /* suppress /gene */

          grp = GeneRefNew ();
          if (grp != NULL) {
            xref = SeqFeatXrefNew ();
            sfp->xref = xref;
            if (xref != NULL) {
              xref->data.choice = SEQFEAT_GENE;
              xref->data.value.ptrvalue = (Pointer) grp;
            }
          }

          sfp->data.value.ptrvalue = (Pointer) ifp;

          if (sap != NULL) {
            if (sap->data != NULL) {
              prevsfp = sap->data;
              while (prevsfp->next != NULL) {
                prevsfp = prevsfp->next;
              }
              prevsfp->next = sfp;
            } else {
              sap->data = (Pointer) sfp;
            }
          }

          sfp->location = AddIntervalToLocation (NULL, sip, (Int4) start, (Int4) stop, FALSE, FALSE);
        }
      }
    }
    SeqLocFree (slp);
  }
  locs = ValNodeFree (locs);

  if (sap != NULL) {
    if (bsp->annot != NULL) {
      prevsap = bsp->annot;
      while (prevsap->next != NULL) {
        prevsap = prevsap->next;
      }
      prevsap->next = sap;
    } else {
      bsp->annot = sap;
    }
  }
}

static void SimpleVecScreenCommon (IteM i, CharPtr database)

{
  BaseFormPtr   bfp;
  BioseqPtr     bsp;
  Char          date [32];
  DatePtr       dp;
  Int4          max;
  Char          path [PATH_MAX];
  SeqEntryPtr   sep;
  VsData        vsd;

#ifdef WIN_MAC
  bfp = (BaseFormPtr) currentFormDataPtr;
#else
  bfp = (BaseFormPtr) GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  path [0] = '\0';
  GetAppParam ("NCBI", "NCBI", "DATA", "", path, sizeof (path));
  FileBuildPath (path, NULL, database);

  date [0] = '\0';
  dp = DateCurr ();
  DatePrint (dp, date);
  DateFree (dp);

  vsd.date = date;
  vsd.path = path;
  vsd.database = database;
  vsd.changed = FALSE;
  vsd.count = 0;
  vsd.mon = NULL;

  bsp =  GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
  if (bsp != NULL) {
    DoVecScreens (bsp, (Pointer) &vsd);
  } else {
    sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
    if (sep == NULL) return;
    max = 0;
    VisitBioseqsInSep (sep, (Pointer) &max, CountVecScreens);
    if (max > 2) {
      vsd.mon = MonitorIntNewEx ("VecScreen Progress", 0, max - 1, FALSE);
    }
    VisitBioseqsInSep (sep, (Pointer) &vsd, DoVecScreens);
    if (vsd.mon != NULL) {
      vsd.mon = MonitorFree (vsd.mon);
      Update ();
    }
  }

  if (vsd.changed) {
    ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
    Update ();
  } else {
    Message (MSG_POST, "No vector contamination found");
  }
}

extern void SimpleUniVecScreenProc (IteM i)

{
  SimpleVecScreenCommon (i, "UniVec");
}

extern void SimpleUniVecCoreScreenProc (IteM i)

{
  SimpleVecScreenCommon (i, "UniVec_Core");
}

static QBQUEUE  qbquerylist = NULL;

static void LIBCALLBACK QBAnnounceCallback (CharPtr requestID, CharPtr seqID, Int2 estimatedSeconds)

{
  if (StringHasNoText (requestID)) {
    requestID = "?";
  }
  if (StringHasNoText (seqID)) {
    seqID = "?";
  }
  Message (MSG_POST, "Queued rID %s, seqID %s, estimated seconds = %d",
           requestID, seqID, (int) estimatedSeconds);
}

static Boolean LIBCALLBACK QBCallback (
  CharPtr filename,
  VoidPtr userdata,
  CharPtr requestID,
  CharPtr seqID,
  Boolean success
)

{
  if (StringHasNoText (requestID)) {
    requestID = "?";
  }
  if (StringHasNoText (seqID)) {
    seqID = "?";
  }
  if (success) {
    if (! SequinHandleNetResults (filename)) {
      /* LaunchGeneralTextViewer (path, "QueueFastaQueryToURL failed"); */
    }
  } else {
    Message (MSG_POST, "Failure of rID '%s', seqID %s", requestID, seqID);
  }
  return TRUE;
}

static void SimpleQBlastProc (IteM i)

{
  BaseFormPtr   bfp;
  BioseqPtr     bsp;
  NewObjectPtr  nop;
  SeqEntryPtr   sep = NULL;
  Int2          which;

  nop = (NewObjectPtr) GetObjectExtra (i);
  if (nop == NULL) return;
#ifdef WIN_MAC
  bfp = (BaseFormPtr) currentFormDataPtr;
#else
  bfp = nop->bfp;
#endif
  if (bfp == NULL) return;
  which = BioseqViewOrDocSumChoice (nop);
  if (which != 1) return;
  bsp =  GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
  if (bsp != NULL) {
    sep = SeqMgrGetSeqEntryForData (bsp);
  } else {
    sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  }
  if (sep == NULL) return;

  QBlastAsynchronousRequest ("nr", "blastn", bsp, &qbquerylist, QBCallback, QBAnnounceCallback, NULL);
}

/* Analysis menu can launch external programs or use Web services */

static QUEUE  urlquerylist = NULL;

static Int4 pendingqueries = 0;

static Boolean LIBCALLBACK SubmitToNCBIResultProc (CONN conn, VoidPtr userdata, EIO_Status status)

{
  AsnIoPtr     aop;
  FILE         *fp;
  Char         path [PATH_MAX];
  SeqEntryPtr  sep;

  TmpNam (path);
  fp = FileOpen (path, "wb");
  QUERY_CopyResultsToFile (conn, fp);
  FileClose (fp);
  aop = AsnIoOpen (path, "rb");
  sep = SeqEntryAsnRead (aop, NULL);
  AsnIoClose (aop);
  aop = AsnIoOpen (path, "w");
  SeqEntryAsnWrite (sep, aop, NULL);
  AsnIoFlush (aop);
  AsnIoClose (aop);
  LaunchGeneralTextViewer (path, "Echo binary transformation of Seq-entry");
  FileRemove (path);
  return TRUE;
}

extern void SubmitToNCBI (IteM i)

{
  AsnIoPtr     aop;
  BaseFormPtr  bfp;
  CONN         conn;
  FILE         *fp;
  Char         path [PATH_MAX];
  Char         progname [64];
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  sprintf (progname, "Sequin/%s", SEQUIN_APPLICATION);

  TmpNam (path);

  aop = AsnIoOpen (path, "wb");
  SeqEntryAsnWrite (sep, aop, NULL);
  AsnIoFlush (aop);
  AsnIoClose (aop);

  conn = QUERY_OpenUrlQuery ("cruncher.nlm.nih.gov", 80,
                             "/cgi-bin/Sequin/testcgi.cgi", "request=echo",
                             progname, 30, eMIME_T_NcbiData, eMIME_AsnBinary,
                             eENCOD_Url,
                             fHCC_UrlDecodeInput | fHCC_UrlEncodeOutput);

  fp = FileOpen (path, "rb");
  QUERY_CopyFileToQuery (conn, fp);
  FileClose (fp);

  QUERY_SendQuery (conn);

  QUERY_AddToQueue (&urlquerylist, conn, SubmitToNCBIResultProc, NULL, TRUE);

  pendingqueries++;

  FileRemove (path);
}

static QUEUE  cddquerylist = NULL;

static Int2  cddquerynum = 0;

#include <cddapi.h>

#define CDD_EXPECT_VALUE 0.01

static Boolean LIBCALLBACK CddProc (
  CONN conn,
  VoidPtr userdata,
  EIO_Status status
)

{
  BioseqPtr    bsp;
  SeqAnnotPtr  prev;
  SeqAnnotPtr  sap;

  if (conn == NULL || userdata == NULL) return TRUE;
  if (status != eIO_Success) return TRUE;
  sap = CddReadReply (conn, status);
  if (sap == NULL) return FALSE;

  bsp = (BioseqPtr) userdata;
  CddCorrectIDs (bsp, sap);

  prev = bsp->annot;
  if (prev == NULL) {
    bsp->annot = sap;
  } else {
    while (prev->next != NULL) {
      prev = prev->next;
    }
    prev->next = sap;
  }

  ObjMgrSetDirtyFlag (bsp->idx.entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bsp->idx.entityID, 0, 0);

  return TRUE;
}

static void SearchCDD (BioseqPtr bsp, Pointer userdata)

{
  BoolPtr  dofeats;

  if (bsp == NULL) return;
  if (! ISA_aa (bsp->mol)) return;

  dofeats = (BoolPtr) userdata;
  if (! CddAsynchronousQuery (bsp, CDD_EXPECT_VALUE, TRUE, TRUE, *dofeats, "cdd", FALSE, &cddquerylist, CddProc, (Pointer) bsp)) {
    ErrPostEx (SEV_ERROR, 0, 0, "Unable to run CDD search");
  } else {
    cddquerynum++;
  }
}

/*
static void SearchCDD (BioseqPtr bsp, Pointer userdata)

{
  BoolPtr      dofeats;
  SeqAnnotPtr  prev;
  SeqAnnotPtr  sap;

  if (bsp == NULL) return;
  if (! ISA_aa (bsp->mol)) return;

  dofeats = (BoolPtr) userdata;
  sap = CddSynchronousQuery (bsp, CDD_EXPECT_VALUE, TRUE, TRUE, *dofeats, "cdd", FALSE);
  if (sap == NULL) return;
  CddCorrectIDs (bsp, sap);

  prev = bsp->annot;
  if (prev == NULL) {
    bsp->annot = sap;
  } else {
    while (prev->next != NULL) {
      prev = prev->next;
    }
    prev->next = sap;
  }
}
*/

static void SimpleCDDSearchCommonProc (IteM i, Boolean makeFeats)

{
  BaseFormPtr  bfp;
  Boolean      dofeats;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  WatchCursor ();
  Update ();

  if (makeFeats) {
    FreeCDDRegions (sep);
  } else {
    FreeCDDAligns (sep);
  }

  dofeats = makeFeats;
  VisitBioseqsInSep (sep, (Pointer) &dofeats, SearchCDD);

  /*
  RemoveDuplicateCDDs (sep);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  */
  ArrowCursor ();
  Update ();
 }

void SimpleCDDSearchFeatProc (IteM i)

{
  SimpleCDDSearchCommonProc (i, TRUE);
}

void SimpleCDDSearchAlignProc (IteM i)

{
  SimpleCDDSearchCommonProc (i, FALSE);
}

extern void SequinCheckSocketsProc (void)

{
  Int4  remaining;

  remaining = QUERY_CheckQueue (&urlquerylist);
  if (remaining < pendingqueries) {
    Beep ();
    pendingqueries--;
  }
  remaining = VecScreenCheckQueue (&vsquerylist);
  if (remaining < vsquerynum) {
    vsquerynum = remaining;
    Message (MSG_POST, "There are %d vector screens remaining", (int) vsquerynum);
  }
  remaining = QBlastCheckQueue (&qbquerylist);
  remaining = CddCheckQueue (&cddquerylist);
  if (remaining < cddquerynum) {
    cddquerynum = remaining;
    Message (MSG_POST, "There are %d cdd searches remaining", (int) cddquerynum);
  }
}

static Boolean LIBCALLBACK DemoModeResultProc (CONN conn, VoidPtr userdata, EIO_Status status)

{
  FILE  *fp;
  Char  path [PATH_MAX];

  TmpNam (path);
  fp = FileOpen (path, "w");
  QUERY_CopyResultsToFile (conn, fp);
  FileClose (fp);
  LaunchGeneralTextViewer (path, "QueueFastaQueryToURL results");
  FileRemove (path);
  return TRUE;
}

static Boolean LIBCALLBACK SequinHandleURLResults (CONN conn, VoidPtr userdata, EIO_Status status)

{
  FILE  *fp;
  Char  path [PATH_MAX];

  TmpNam (path);
  fp = FileOpen (path, "w");
  QUERY_CopyResultsToFile (conn, fp);
  FileClose (fp);
  if (! SequinHandleNetResults (path)) {
    /* LaunchGeneralTextViewer (path, "QueueFastaQueryToURL failed"); */
  }
  FileRemove (path);
  return TRUE;
}

static void FinishURLProc (NewObjectPtr nop, CharPtr arguments, CharPtr path)

{
  CONN             conn;
  FILE             *fp;
  Char             progname [64];
  QueryResultProc  resultproc;
  EMIME_SubType    subtype;

  sprintf (progname, "Sequin/%s", SEQUIN_APPLICATION);

  if (nop->demomode) {
    resultproc = DemoModeResultProc;
  } else {
    resultproc = nop->resultproc;
  }
  if (nop->format == 1) {
    subtype = eMIME_Fasta;
  } else if (nop->format == 2) {
    subtype = eMIME_AsnText;
  } else {
    subtype = eMIME_Unknown;
  }

  conn = QUERY_OpenUrlQuery (nop->host_machine, nop->host_port,
                             nop->host_path, arguments,
                             progname, nop->timeoutsec,
                             eMIME_T_NcbiData, subtype, eENCOD_Url,
                             fHCC_UrlDecodeInput | fHCC_UrlEncodeOutput);
  if (conn == NULL) return;

  fp = FileOpen (path, "r");
  QUERY_CopyFileToQuery (conn, fp);
  FileClose (fp);

  QUERY_SendQuery (conn);

  QUERY_AddToQueue (&urlquerylist, conn, resultproc, NULL, TRUE);

  pendingqueries++;
}

static void DoAnalysisProc (NewObjectPtr nop, BaseFormPtr bfp, Int2 which, CharPtr arguments, ResultProc dotheanalysis)

{
  AsnIoPtr     aop;
  BioseqPtr    bsp;
  Char         path1 [PATH_MAX];
  SeqEntryPtr  sep;

  if (nop == NULL || bfp == NULL) return;
  switch (which) {
    case 1 :
      if (BioseqViewCanSaveFasta (bfp->form, nop->fastaNucOK, nop->fastaProtOK, nop->onlyBspTarget)) {
        TmpNam (path1);
        switch (nop->format) {
          case 1 :
            ExportBioseqViewFasta (bfp->form, path1, nop->fastaNucOK, nop->fastaProtOK, nop->onlyBspTarget);
            break;
          case 2 :
            sep = NULL;
            if (nop->onlyBspTarget) {
              bsp =  GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
              sep = SeqMgrGetSeqEntryForData (bsp);
            } else {
              sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
            }
            if (sep != NULL) {
              aop = AsnIoOpen (path1, "w");
              SeqEntryAsnWrite (sep, aop, NULL);
              AsnIoFlush (aop);
              AsnIoClose (aop);
            }
            break;
          default :
            break;
        }
        if (dotheanalysis != NULL) {
          dotheanalysis (path1);
        } else {
          FinishURLProc (nop, arguments, path1);
        }
        FileRemove (path1);
      } else {
        ErrPostEx (SEV_ERROR, 0, 0, "BioseqView cannot save fasta format");
      }
      break;
    case 2 :
      if (DocSumCanSaveFasta (bfp->form, nop->fastaNucOK, nop->fastaProtOK)) {
        TmpNam (path1);
        ExportDocSumFasta (bfp->form, path1, nop->fastaNucOK, nop->fastaProtOK);
        if (dotheanalysis != NULL) {
          dotheanalysis (path1);
        } else {
          FinishURLProc (nop, arguments, path1);
        }
        FileRemove (path1);
      } else {
        ErrPostEx (SEV_ERROR, 0, 0, "DocSum cannot save fasta format");
      }
      break;
    default :
      break;
  }
}

/* encodes spaces as %20 in URLs */
static CharPtr StrSaveNoNullEncodeSpaces (CharPtr from)

{
  Char     ch;
  size_t   len = 0;
  CharPtr  p;
  CharPtr  q;
  CharPtr  to;

  if (StringHasNoText (from)) return NULL;
  p = from;
  ch = *p;
  while (ch != '\0') {
    if (ch == ' ') {
      len += 3;
    } else {
      len++;
    }
    p++;
    ch = *p;
  }
  to = MemNew (len + 1);
  if (to == NULL) return NULL;

  q = to;
  p = from;
  ch = *p;
  while (ch != '\0') {
    if (ch == ' ') {
      *q = '%';
      q++;
      *q = '2';
      q++;
      *q = '0';
      q++;
    } else {
      *q = ch;
      q++;
    }
    p++;
    ch = *p;
  }
  *q = '\0';
  return to;
}

typedef struct urlargform {
  FORM_MESSAGE_BLOCK

  NewObjectPtr nop;
  BaseFormPtr  bfp;
  ValNodePtr   controls;
  ValNodePtr   helps;
  Int2         which;
} UrlArgForm, PNTR UrlArgFormPtr;

static void AcceptArgumentFormProc (ButtoN b)

{
  CharPtr        args = NULL;
  CharPtr        arguments = NULL;
  ButtoN         btn;
  Char           ch;
  Int2           choice;
  Char           cpy [256];
  GrouP          grp;
  ValNodePtr     head = NULL;
  Int2           i;
  CharPtr        itms;
  CharPtr        last;
  size_t         len;
  LisT           lst;
  NewObjectPtr   nop;
  Boolean        notFirst = FALSE;
  PopuP          pop;
  ValNodePtr     ppt;
  CharPtr        ptr;
  CharPtr        str;
  Char           tmp [256];
  TexT           txt;
  UrlArgFormPtr  ufp;
  UrlParamPtr    upp;
  Int2           val;
  ValNodePtr     vnp;

  ufp = (UrlArgFormPtr) GetObjectExtra (b);
  if (ufp == NULL) return;
  Hide (ufp->form);
  Update ();
  nop = ufp->nop;
  if (nop != NULL) {
    if (! StringHasNoText (nop->prefix)) {
      ValNodeCopyStr (&head, 0, nop->prefix);
    }
    for (vnp = ufp->controls, ppt = nop->paramlist;
         vnp != NULL && ppt != NULL;
         vnp = vnp->next, ppt = ppt->next) {
      upp = (UrlParamPtr) ppt->data.ptrvalue;
      if (upp == NULL) continue;
      choice = vnp->choice;
      switch (upp->type) {
        case 1 :
          txt = (TexT) vnp->data.ptrvalue;
          str = SaveStringFromText (txt);
          if (str != NULL) {
            sprintf (tmp, "%s=%s", upp->param, str);
            ValNodeCopyStr (&head, ppt->choice, tmp);
            MemFree (str);
          }
          break;
        case 2 :
          btn = (ButtoN) vnp->data.ptrvalue;
          if (GetStatus (btn)) {
            sprintf (tmp, "%s=TRUE", upp->param);
          } else {
            sprintf (tmp, "%s=FALSE", upp->param);
          }
          ValNodeCopyStr (&head, ppt->choice, tmp);
          break;
        case 3 :
          pop = (PopuP) vnp->data.ptrvalue;
          val = GetValue (pop);
          if (val > 0) {
            i = 0;
            itms = upp->choices;
            StringNCpy_0 (tmp, itms, sizeof (tmp));
            last = tmp;
            ptr = last;
            ch = *ptr;
            while (ch != '\0') {
              if (ch == ',') {
                *ptr = '\0';
                i++;
                if (val == i) {
                  sprintf (cpy, "%s=%s", upp->param, last);
                  ValNodeCopyStr (&head, ppt->choice, cpy);
                }
                ptr++;
                last = ptr;
                ch = *ptr;
              } else {
                ptr++;
                ch = *ptr;
              }
            }
            if (! StringHasNoText (last)) {
              i++;
              if (val == i) {
                sprintf (cpy, "%s=%s", upp->param, last);
                ValNodeCopyStr (&head, ppt->choice, cpy);
              }
            }
          }
          break;
        case 4 :
          grp = (GrouP) vnp->data.ptrvalue;
          val = GetValue (grp);
          if (val > 0) {
            i = 0;
            itms = upp->choices;
            StringNCpy_0 (tmp, itms, sizeof (tmp));
            last = tmp;
            ptr = last;
            ch = *ptr;
            while (ch != '\0') {
              if (ch == ',') {
                *ptr = '\0';
                i++;
                if (val == i) {
                  sprintf (cpy, "%s=%s", upp->param, last);
                  ValNodeCopyStr (&head, ppt->choice, cpy);
                }
                ptr++;
                last = ptr;
                ch = *ptr;
              } else {
                ptr++;
                ch = *ptr;
              }
            }
            if (! StringHasNoText (last)) {
              i++;
              if (val == i) {
                sprintf (cpy, "%s=%s", upp->param, last);
                ValNodeCopyStr (&head, ppt->choice, cpy);
              }
            }
          }
          break;
        case 5 :
          lst = (LisT) vnp->data.ptrvalue;
          val = GetValue (lst);
          if (val > 0) {
            i = 0;
            itms = upp->choices;
            StringNCpy_0 (tmp, itms, sizeof (tmp));
            last = tmp;
            ptr = last;
            ch = *ptr;
            while (ch != '\0') {
              if (ch == ',') {
                *ptr = '\0';
                i++;
                if (val == i) {
                  sprintf (cpy, "%s=%s", upp->param, last);
                  ValNodeCopyStr (&head, ppt->choice, cpy);
                }
                ptr++;
                last = ptr;
                ch = *ptr;
              } else {
                ptr++;
                ch = *ptr;
              }
            }
            if (! StringHasNoText (last)) {
              i++;
              if (val == i) {
                sprintf (cpy, "%s=%s", upp->param, last);
                ValNodeCopyStr (&head, ppt->choice, cpy);
              }
            }
          }
          break;
        default :
          break;
      }
    }
    head = SortValNode (head, SortByVnpChoice);
    if (! StringHasNoText (nop->suffix)) {
      ValNodeCopyStr (&head, 0, nop->suffix);
    }
    for (len = 0, vnp = head; vnp != NULL; vnp = vnp->next) {
      len += StringLen ((CharPtr) vnp->data.ptrvalue) + 1;
    }
    if (len > 0) {
      arguments = MemNew (len + 5);
      if (arguments != NULL) {
        vnp = head;
        notFirst = FALSE;
        while (vnp != NULL) {
          if (notFirst) {
            StringCat (arguments, "&");
          }
          StringCat (arguments, (CharPtr) vnp->data.ptrvalue);
          notFirst = TRUE;
          vnp = vnp->next;
        }
      }
    }
    args = /* StrSaveNoNullEncodeSpaces */ StringSave (arguments);
    MemFree (arguments);
    DoAnalysisProc (nop, ufp->bfp, ufp->which, args, NULL);
    MemFree (args);
  }
  Remove (ufp->form);
}

static void ShowArgumentHelp (ButtoN b)

{
  NewObjectPtr   nop;
  ValNodePtr     ppt;
  CharPtr        str;
  UrlArgFormPtr  ufp;
  UrlParamPtr    upp;
  ValNodePtr     vnp;

  ufp = (UrlArgFormPtr) GetObjectExtra (b);
  if (ufp == NULL) return;
  nop = ufp->nop;
  if (nop == NULL) return;
  for (vnp = ufp->helps, ppt = nop->paramlist;
         vnp != NULL && ppt != NULL;
         vnp = vnp->next, ppt = ppt->next) {
    upp = (UrlParamPtr) ppt->data.ptrvalue;
    if (upp == NULL) continue;
    if ((Pointer) b == (Pointer) vnp->data.ptrvalue) {
      str = upp->help;
      Message (MSG_OK, "%s", str);
      return;
    }
  }
  Beep ();
}

static void ArgumentFormMessage (ForM f, Int2 mssg)

{
  UrlArgFormPtr  ufp;

  ufp = (UrlArgFormPtr) GetObjectExtra (f);
  if (ufp != NULL) {
    switch (mssg) {
      case VIB_MSG_CLOSE :
        Remove (f);
        break;
      default :
        if (ufp->appmessage != NULL) {
          ufp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static void CleanupArgumentForm (GraphiC g, VoidPtr data)

{
  UrlArgFormPtr  ufp;

  ufp = (UrlArgFormPtr) data;
  if (ufp != NULL) {
    ValNodeFree (ufp->controls);
    ValNodeFree (ufp->helps);
  }
  StdCleanupFormProc (g, data);
}

static ValNodePtr RearrangeParamList (ValNodePtr paramlist)

{
  ValNodePtr       curr;
  CharPtr          group;
  ValNodePtr       head = NULL;
  ValNodePtr       list;
  ValNodePtr       next;
  ValNodePtr PNTR  prev;
  ValNodePtr       ppt;
  UrlParamPtr      upp;

  ppt = paramlist;
  while (ppt != NULL) {
    list = ppt->next;
    ppt->next = NULL;
    ValNodeLink (&head, ppt);
    upp = (UrlParamPtr) ppt->data.ptrvalue;
    if (upp == NULL) {
      ppt = list;
      continue;
    }
    group = upp->group;
    curr = list;
    prev = &list;
    while (curr != NULL) {
      next = curr->next;
      upp = (UrlParamPtr) curr->data.ptrvalue;
      if (upp == NULL) {
        prev = &(curr->next);
        curr = next;
        continue;
      }
      if (StringICmp (upp->group, group) == 0) {
        *prev = next;
        curr->next = NULL;
        ValNodeLink (&head, curr);
      } else {
        prev = &(curr->next);
      }
      curr = next;
    }
    ppt = list;
  }
  return head;
}

static void BuildArgumentForm (NewObjectPtr nop, BaseFormPtr bfp, Int2 which)

{
  ButtoN             b;
  ButtoN             btn;
  GrouP              c;
  Char               ch;
  CharPtr            def;
  Int2               delta;
  TexT               first = NULL;
  GrouP              g;
  GrouP              grp;
  GrouP              h;
  ValNodePtr         hlp;
  Int2               i;
  CharPtr            itms;
  CharPtr            last;
  CharPtr            lastGroup = " ";
  LisT               lst;
  GrouP              m;
  Int2               max;
  Int2               min;
  ValNodePtr         moveMe = NULL;
  Nlm_Handle         obj1, obj2;
  PopuP              pop;
  PrompT             prmpt;
  ValNodePtr         ppt;
  CharPtr            ptr;
  RecT               r;
  StdEditorProcsPtr  sepp;
  CharPtr            str;
  Char               tmp [128];
  TexT               txt;
  UrlArgFormPtr      ufp;
  UrlParamPtr        upp;
  Int2               val;
  ValNodePtr         vnp;
  WindoW             w;

  if (nop == NULL || bfp == NULL) return;
  ufp = (UrlArgFormPtr) MemNew (sizeof (UrlArgForm));
  if (ufp == NULL) return;

  nop->paramlist = RearrangeParamList (nop->paramlist);

  w = FixedWindow (-50, -33, -10, -10, "Arguments", NULL);
  SetObjectExtra (w, ufp, CleanupArgumentForm);
  ufp->form = (ForM) w;
  ufp->formmessage = ArgumentFormMessage;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    ufp->appmessage = sepp->handleMessages;
  }

  ufp->bfp = bfp;
  ufp->nop = nop;
  ufp->which = which;

  m = HiddenGroup (w, 1, 0, NULL);

  g = NULL;
  for (ppt = nop->paramlist;
       ppt != NULL;
       ppt = ppt->next) {
    upp = (UrlParamPtr) ppt->data.ptrvalue;
    if (upp == NULL) continue;
    if (StringICmp (upp->group, lastGroup) != 0) {
      if (StringHasNoText (upp->group)) {
        if (StringHasNoText (lastGroup)) {
          g = HiddenGroup (m, 3, 0, NULL);
        } else {
          g = NormalGroup (m, 3, 0, "", programFont, NULL);
        }
      } else {
        g = NormalGroup (m, 3, 0, upp->group, programFont, NULL);
      }
      lastGroup = upp->group;
    }
    if (g == NULL) {
      g = HiddenGroup (m, 3, 0, NULL);
    }
    switch (upp->type) {
      case 1 :
        str = upp->prompt;
        StaticPrompt (g, str, 0, dialogTextHeight, programFont, 'l');
        def = upp->dfault;
        if (StringHasNoText (def)) {
          def = "";
        }
        txt = DialogText (g, def, 10, NULL);
        if (first == NULL) {
          first = txt;
        }
        ValNodeAddPointer (&(ufp->controls), 1, (Pointer) txt);
        ValNodeAddPointer (&moveMe, 0, (Pointer) txt);
        b = PushButton (g, "?", ShowArgumentHelp);
        SetObjectExtra (b, ufp, NULL);
        ValNodeAddPointer (&(ufp->helps), 0, (Pointer) b);
        break;
      case 2 :
        str = upp->prompt;
        btn = CheckBox (g, str, NULL);
        def = upp->dfault;
        if (StringICmp (def, "TRUE") == 0) {
          SetStatus (btn, TRUE);
        }
        prmpt = StaticPrompt (g, "", 0, 0, programFont, 'l');
        ValNodeAddPointer (&moveMe, 0, (Pointer) prmpt);
        ValNodeAddPointer (&(ufp->controls), 2, (Pointer) btn);
        b = PushButton (g, "?", ShowArgumentHelp);
        SetObjectExtra (b, ufp, NULL);
        ValNodeAddPointer (&(ufp->helps), 0, (Pointer) b);
        break;
      case 3 :
        str = upp->prompt;
        StaticPrompt (g, str, 0, dialogTextHeight, programFont, 'l');
        h = HiddenGroup (g, 1, 0, NULL);
        pop = PopupList (h, TRUE, NULL);
        def = upp->dfault;
        val = 0;
        i = 0;
        itms = upp->choices;
        StringNCpy_0 (tmp, itms, sizeof (tmp));
        last = tmp;
        ptr = last;
        ch = *ptr;
        while (ch != '\0') {
          if (ch == ',') {
            *ptr = '\0';
            PopupItem (pop, last);
            i++;
            if (StringICmp (def, last) == 0) {
              val = i;
            }
            ptr++;
            last = ptr;
            ch = *ptr;
          } else {
            ptr++;
            ch = *ptr;
          }
        }
        if (! StringHasNoText (last)) {
          PopupItem (pop, last);
          i++;
          if (StringICmp (def, last) == 0) {
            val = i;
          }
        }
        if (val > 0) {
          SetValue (pop, val);
        }
        ValNodeAddPointer (&(ufp->controls), 3, (Pointer) pop);
        ValNodeAddPointer (&moveMe, 0, (Pointer) pop);
        b = PushButton (g, "?", ShowArgumentHelp);
        SetObjectExtra (b, ufp, NULL);
        ValNodeAddPointer (&(ufp->helps), 0, (Pointer) b);
        break;
      case 4 :
        str = upp->prompt;
        StaticPrompt (g, str, 0, dialogTextHeight, programFont, 'l');
        h = HiddenGroup (g, 1, 0, NULL);
        grp = HiddenGroup (h, -3, 0, NULL);
        def = upp->dfault;
        val = 0;
        i = 0;
        itms = upp->choices;
        StringNCpy_0 (tmp, itms, sizeof (tmp));
        last = tmp;
        ptr = last;
        ch = *ptr;
        while (ch != '\0') {
          if (ch == ',') {
            *ptr = '\0';
            RadioButton (grp, last);
            i++;
            if (StringICmp (def, last) == 0) {
              val = i;
            }
            ptr++;
            last = ptr;
            ch = *ptr;
          } else {
            ptr++;
            ch = *ptr;
          }
        }
        if (! StringHasNoText (last)) {
          RadioButton (grp, last);
          i++;
          if (StringICmp (def, last) == 0) {
            val = i;
          }
        }
        if (val > 0) {
          SetValue (grp, val);
        }
        ValNodeAddPointer (&(ufp->controls), 4, (Pointer) grp);
        ValNodeAddPointer (&moveMe, 0, (Pointer) grp);
        b = PushButton (g, "?", ShowArgumentHelp);
        SetObjectExtra (b, ufp, NULL);
        ValNodeAddPointer (&(ufp->helps), 0, (Pointer) b);
        break;
      case 5 :
        str = upp->prompt;
        StaticPrompt (g, str, 0, dialogTextHeight, programFont, 'l');
        h = HiddenGroup (g, 1, 0, NULL);
        lst = SingleList (h, 10, 3, NULL);
        def = upp->dfault;
        val = 0;
        i = 0;
        itms = upp->choices;
        StringNCpy_0 (tmp, itms, sizeof (tmp));
        last = tmp;
        ptr = last;
        ch = *ptr;
        while (ch != '\0') {
          if (ch == ',') {
            *ptr = '\0';
            ListItem (lst, last);
            i++;
            if (StringICmp (def, last) == 0) {
              val = i;
            }
            ptr++;
            last = ptr;
            ch = *ptr;
          } else {
            ptr++;
            ch = *ptr;
          }
        }
        if (! StringHasNoText (last)) {
          ListItem (lst, last);
          i++;
          if (StringICmp (def, last) == 0) {
            val = i;
          }
        }
        if (val > 0) {
          SetValue (lst, val);
        }
        ValNodeAddPointer (&(ufp->controls), 5, (Pointer) lst);
        ValNodeAddPointer (&moveMe, 0, (Pointer) lst);
        b = PushButton (g, "?", ShowArgumentHelp);
        SetObjectExtra (b, ufp, NULL);
        ValNodeAddPointer (&(ufp->helps), 0, (Pointer) b);
        break;
      default :
        break;
    }
  }

  min = 0;
  max = 0;
  for (vnp = moveMe; vnp != NULL; vnp = vnp->next) {
    obj1 = (Nlm_Handle) vnp->data.ptrvalue;
    GetPosition (obj1, &r);
    min = MAX (min, r.left);
  }
  for (vnp = moveMe; vnp != NULL; vnp = vnp->next) {
    obj1 = (Nlm_Handle) vnp->data.ptrvalue;
    GetPosition (obj1, &r);
    delta = min - r.left;
    OffsetRect (&r, delta, 0);
    SetPosition (obj1, &r);
    AdjustPrnt (obj1, &r, FALSE);
    max = MAX (max, r.right);
  }
  max += 3;
  for (vnp = moveMe, hlp = ufp->helps;
       vnp != NULL && hlp != NULL;
       vnp = vnp->next, hlp = hlp->next) {
    obj2 = (Nlm_Handle) hlp->data.ptrvalue;
    GetPosition (obj2, &r);
    delta = max - r.left;
    OffsetRect (&r, delta, 0);
    SetPosition (obj2, &r);
    AdjustPrnt (obj2, &r, TRUE);
  }

  c = HiddenGroup (w, 2, 0, NULL);
  SetGroupSpacing (c, 10, 3);
  b = DefaultButton (c, "Accept", AcceptArgumentFormProc);
  SetObjectExtra (b, ufp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) m, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Select (w);
  if (first != NULL) {
    Select (first);
  }
}

static void DoURLProc (IteM i)

{
  CharPtr       args = NULL;
  BaseFormPtr   bfp;
  size_t        len;
  NewObjectPtr  nop;
  Int2          which;

  nop = (NewObjectPtr) GetObjectExtra (i);
  if (nop == NULL) return;
#ifdef WIN_MAC
  bfp = (BaseFormPtr) currentFormDataPtr;
#else
  bfp = nop->bfp;
#endif
  if (bfp == NULL) return;
  which = BioseqViewOrDocSumChoice (nop);
  if (nop->paramlist == NULL) {
    len = StringLen (nop->prefix) + StringLen (nop->suffix);
    if (len > 0) {
      args = MemNew (sizeof (Char) * (len + 2));
      StringCpy (args, nop->prefix);
      if (! StringHasNoText (nop->suffix)) {
        StringCat (args, "&");
        StringCat (args, nop->suffix);
      }
    }
    DoAnalysisProc (nop, bfp, which, args, NULL);
  } else {
    BuildArgumentForm (nop, bfp, which);
  }
}

extern void EnableAnalysisItems (BaseFormPtr bfp, Boolean isDocSum)

{
  Boolean       hasFastaNuc;
  Boolean       hasFastaProt;
  NewObjectPtr  nop;

  if (bfp == NULL) return;
#ifdef WIN_MAC
  nop = (NewObjectPtr) macUserDataPtr;
#else
  nop = (NewObjectPtr) bfp->userDataPtr;
#endif
  if (isDocSum) {
  } else {
  }
  while (nop != NULL) {
    if (nop->kind == 1) {
      /* annotate menu item, ignore it */
    } else if (isDocSum) {
      if (nop->dsmOK) {
        hasFastaNuc = DocSumCanSaveFasta (bfp->form, TRUE, FALSE);
        hasFastaProt = DocSumCanSaveFasta (bfp->form, FALSE, TRUE);
        if (nop->fastaNucOK && hasFastaNuc) {
          SafeEnable (nop->item);
        } else if (nop->fastaProtOK && hasFastaProt) {
          SafeEnable (nop->item);
        } else {
          SafeDisable (nop->item);
        }
      } else {
        SafeDisable (nop->item);
      }
    } else {
      if (nop->bspOK) {
        hasFastaNuc = BioseqViewCanSaveFasta (bfp->form, TRUE, FALSE, nop->onlyBspTarget);
        hasFastaProt = BioseqViewCanSaveFasta (bfp->form, FALSE, TRUE, nop->onlyBspTarget);
        if (nop->fastaNucOK && hasFastaNuc) {
          SafeEnable (nop->item);
        } else if (nop->fastaProtOK && hasFastaProt) {
          SafeEnable (nop->item);
        } else {
          SafeDisable (nop->item);
        }
      } else {
        SafeDisable (nop->item);
      }
    }
    nop = nop->next;
  }
}

static VoidPtr LinkNewObjectLists (NewObjectPtr list1, NewObjectPtr list2)

{
  NewObjectPtr  nop;

  if (list1 == NULL) return list2;
  nop = list1;
  while (nop->next != NULL) {
    nop = nop->next;
  }
  nop->next = list2;
  return list1;
}

static void CleanupAnalysisExtraProc (GraphiC g, VoidPtr data)

{
  NewObjectPtr  nop;
  ValNodePtr    ppt;
  UrlParamPtr   upp;

  nop = (NewObjectPtr) data;
  if (nop != NULL) {
    MemFree (nop->host_machine);
    MemFree (nop->host_path);
    for (ppt = nop->paramlist; ppt != NULL; ppt = ppt->next) {
      upp = (UrlParamPtr) ppt->data.ptrvalue;
      if (upp == NULL) continue;
      MemFree (upp->param);
      MemFree (upp->prompt);
      MemFree (upp->dfault);
      MemFree (upp->choices);
      MemFree (upp->group);
      MemFree (upp->help);
    }
    ValNodeFreeData (nop->paramlist);
    MemFree (nop->prefix);
    MemFree (nop->suffix);
  }
  MemFree (data);
}

typedef struct sbstruc {
  CharPtr    name;
  MenU       menu;
} Sbstruc, PNTR SbstrucPtr;

static ValNodePtr  analysissubmenulist = NULL;

static void AddAnalysisItem (MenU m, BaseFormPtr bfp,
                             Boolean bspviewOK, Boolean docsumOK,
                             Boolean nucOK, Boolean protOK, Boolean onlyBspTarget,
                             CharPtr host_machine, Uint2 host_port,
                             CharPtr host_path, CharPtr program,
                             Uint2 timeoutsec, Int2 format, Boolean demomode,
                             QueryResultProc resultproc, ValNodePtr paramlist,
                             CharPtr prefix, CharPtr suffix,
                             CharPtr title, CharPtr submenu,
                             ItmActnProc actn, NewObjectPtr PNTR head)

{
  IteM          i;
  NewObjectPtr  last;
  size_t        len;
  NewObjectPtr  nop;
  SbstrucPtr    sbp;
  CharPtr       tmp;
  ValNodePtr    vnp;
  MenU          x;

  if (m == NULL || actn == NULL) return;
  x = NULL;
  if (! StringHasNoText (submenu)) {
    vnp = analysissubmenulist;
    while (vnp != NULL && x == NULL) {
      sbp = (SbstrucPtr) vnp->data.ptrvalue;
      if (sbp != NULL && StringICmp (sbp->name, submenu) == 0) {
        x = sbp->menu;
      }
      vnp = vnp->next;
    }
    if (x == NULL) {
      sbp = (SbstrucPtr) MemNew (sizeof (Sbstruc));
      if (sbp != NULL) {
        sbp->name = StringSave (submenu);
        sbp->menu = SubMenu (m, sbp->name);
        x = sbp->menu;
        ValNodeAddPointer (&analysissubmenulist, 0, (VoidPtr) sbp);
      }
    }
  }
  if (x == NULL) {
    x = m;
  }
  i = CommandItem (x, title, actn);
  nop = (NewObjectPtr) MemNew (sizeof (NewObjectData));
  if (nop != NULL) {
    nop->kind = 2; /* analysis menu item */
    nop->bfp = bfp;
    nop->item = i;
    nop->bspOK = bspviewOK;
    nop->dsmOK = docsumOK;
    nop->fastaNucOK = nucOK;
    nop->fastaProtOK = protOK;
    nop->onlyBspTarget = onlyBspTarget;
    nop->host_machine = /* StrSaveNoNullEncodeSpaces */ StringSave (host_machine);
    nop->host_port = host_port;
    len = StringLen (host_path);
    tmp = MemNew (len + StringLen (program) + 5);
    if (tmp != NULL) {
      StringCpy (tmp, host_path);
      if (len > 1 && tmp [len - 1] != '/') {
        StringCat (tmp, "/");
      }
      StringCat (tmp, program);
    }
    nop->host_path = /* StrSaveNoNullEncodeSpaces */ StringSave (tmp);
    MemFree (tmp);
    nop->query = NULL;
    /*
    nop->host_path = StrSaveNoNullEncodeSpaces (host_path);
    nop->query = StrSaveNoNullEncodeSpaces (program);
    */
    nop->timeoutsec = timeoutsec;
    nop->format = format;
    nop->demomode = demomode;
    nop->resultproc = resultproc;
    nop->paramlist = paramlist;
    nop->prefix = StringSaveNoNull (prefix);
    nop->suffix = StringSaveNoNull (suffix);
  }
  SetObjectExtra (i, (Pointer) nop, CleanupAnalysisExtraProc);
  if (head == NULL) return;
  last = *head;
  if (last != NULL) {
    while (last->next != NULL) {
     last = last->next;
    }
    last->next = nop;
  } else {
    *head = nop;
  }
}

/* Sample seqncgis.cnf/seqncgis.ini/.seqncgisrc/sequincgi.cfg config file.
   PATH can contain query (separated by ? symbol), or separate QUERY item can
   be used, or multiple QUERY and TITLE items can also be used.

[SERVICES]
PATH=mydisk:Common Files:services:

[ORDER]
ORDER_1=tRNAscan
ORDER_2=Seg

[tRNAscan]
PROGRAM=testcgi.cgi?request=trnascan
HOST=www.myserver.myschool.edu
PORT=80
PATH=/MyServices/cgi-bin/testcgi.cgi
SUBMENU=Search
FORMATIN=FASTA
FLAGS=SEQ,NUC,TRG
TIMEOUT=30

[Seg]
PROGRAM=segify
HOST=www.myserver.myschool.edu
PORT=80
PATH=/MyServices/cgi-bin/testcgi.cgi
FORMATIN=fasta
FLAGS=SEQ,DOC,PRT,TRG
SUBMENU=Secondary structure prediction
PROMPT_1=Window Size
PARAM_1=window
DESCRIPTION_1=window size for determining low-complexity segments
TYPE_1=text
DEFAULT_1=12
REQUIRED_1=FALSE
IMPORTANCE_1=
GROUP_1=Algorithm
HELP_1=window size for determining low-complexity segments
PROMPT_2=Trigger Complexity
PARAM_2=trigger
DESCRIPTION_2=trigger complexity for determining low-complexity segments
TYPE_2=text
DEFAULT_2=2.2
REQUIRED_2=FALSE
IMPORTANCE_2=
GROUP_2=Algorithm
HELP_2=trigger complexity for determining low-complexity segments
...

[ENZYMES]
ENZ_1=BamHI
ENZ_2=EcoRI
ENZ_3=HindIII

*/

static Int2 GetServiceParam (ValNodePtr head, CharPtr type, CharPtr buf, Int2 buflen)

{
  size_t      len;
  Boolean     seenBracket = FALSE;
  CharPtr     str;
  ValNodePtr  vnp;

  if (buf == NULL || buflen <= 0) return 0;
  *buf = '\0';
  len = StringLen (type);
  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (str != NULL) {
      if (str [0] == '[') {
        if (seenBracket) return 0;
        seenBracket = TRUE;
      } else if (StringNICmp (str, type, len) == 0) {
        str += len;
        StringNCpy_0 (buf, str, buflen);
        return (Int2) StringLen (buf);
      }
    }
  }
  return 0;
}

static ValNodePtr GetConfigParamAndPromptLists (CharPtr sect)

{
  Int2         i;
  ValNodePtr   paramlist = NULL;
  Uint1        paramtype;
  Char         title [512];
  Char         tmp [32];
  UrlParamPtr  upp;

  if (sect == NULL) return NULL;
  i = 1;
  sprintf (tmp, "PARAM_%d", (int) i);
  while (GetAppParam ("SEQNCGIS", sect, tmp, NULL, title, sizeof (title) - 1)) {
    upp = (UrlParamPtr) MemNew (sizeof (UrlParamData));
    if (upp == NULL) continue;
    upp->param = StringSave (title);
    sprintf (tmp, "TYPE_%d", (int) i);
    paramtype = 1;
    if (GetAppParam ("SEQNCGIS", sect, tmp, NULL, title, sizeof (title) - 1)) {
      if (StringICmp (title, "text") == 0) {
        paramtype = 1;
      } else if (StringICmp (title, "checkbox") == 0) {
        paramtype = 2;
      } else if (StringICmp (title, "popup") == 0) {
        paramtype = 3;
      } else if (StringICmp (title, "radio") == 0) {
        paramtype = 4;
      } else if (StringICmp (title, "list") == 0) {
        paramtype = 5;
      }
    }
    upp->type = paramtype;
    sprintf (tmp, "PROMPT_%d", (int) i);
    if (GetAppParam ("SEQNCGIS", sect, tmp, NULL, title, sizeof (title) - 1)) {
      upp->prompt = StringSave (title);
    } else {
      upp->prompt = StringSave (upp->param);
    }
    sprintf (tmp, "DEFAULT_%d", (int) i);
    if (GetAppParam ("SEQNCGIS", sect, tmp, NULL, title, sizeof (title) - 1)) {
      upp->dfault = StringSave (title);
    } else {
      upp->dfault = StringSave (" ");
    }
    sprintf (tmp, "CHOICES_%d", (int) i);
    if (GetAppParam ("SEQNCGIS", sect, tmp, NULL, title, sizeof (title) - 1)) {
      upp->choices = StringSave (title);
    } else {
      upp->choices = StringSave (" ");
    }
    sprintf (tmp, "GROUP_%d", (int) i);
    if (GetAppParam ("SEQNCGIS", sect, tmp, NULL, title, sizeof (title) - 1)) {
      upp->group = StringSave (title);
    } else {
      upp->group = StringSave (" ");
    }
    sprintf (tmp, "HELP_%d", (int) i);
    if (GetAppParam ("SEQNCGIS", sect, tmp, NULL, title, sizeof (title) - 1)) {
      upp->help = StringSave (title);
    } else {
      upp->help = StringSave (" ");
    }
    ValNodeAddPointer (&paramlist, i, (Pointer) upp);
    i++;
    sprintf (tmp, "PARAM_%d", (int) i);
  }
  return paramlist;
}

static ValNodePtr GetServiceParamAndPromptLists (ValNodePtr list)

{
  Int2         i;
  ValNodePtr   paramlist = NULL;
  Uint1        paramtype;
  Char         title [512];
  Char         tmp [32];
  UrlParamPtr  upp;

  if (list == NULL) return NULL;
  i = 1;
  sprintf (tmp, "PARAM_%d=", (int) i);
  while (GetServiceParam (list, tmp, title, sizeof (title) - 1)) {
    upp = (UrlParamPtr) MemNew (sizeof (UrlParamData));
    if (upp == NULL) continue;
    upp->param = StringSave (title);
    sprintf (tmp, "TYPE_%d", (int) i);
    paramtype = 1;
    if (GetServiceParam (list, tmp, title, sizeof (title) - 1)) {
      if (StringICmp (title, "text") == 0) {
        paramtype = 1;
      } else if (StringICmp (title, "checkbox") == 0) {
        paramtype = 2;
      } else if (StringICmp (title, "popup") == 0) {
        paramtype = 3;
      } else if (StringICmp (title, "radio") == 0) {
        paramtype = 4;
      } else if (StringICmp (title, "list") == 0) {
        paramtype = 5;
      }
    }
    upp->type = paramtype;
    sprintf (tmp, "PROMPT_%d=", (int) i);
    if (GetServiceParam (list, tmp, title, sizeof (title) - 1)) {
      upp->prompt = StringSave (title);
    } else {
      upp->prompt = StringSave (upp->param);
    }
    sprintf (tmp, "DEFAULT_%d=", (int) i);
    if (GetServiceParam (list, tmp, title, sizeof (title) - 1)) {
      upp->dfault = StringSave (title);
    } else {
      upp->dfault = StringSave (" ");
    }
    sprintf (tmp, "CHOICES_%d=", (int) i);
    if (GetServiceParam (list, tmp, title, sizeof (title) - 1)) {
      upp->choices = StringSave (title);
    } else {
      upp->choices = StringSave (" ");
    }
    sprintf (tmp, "GROUP_%d=", (int) i);
    if (GetServiceParam (list, tmp, title, sizeof (title) - 1)) {
      upp->group = StringSave (title);
    } else {
      upp->group = StringSave (" ");
    }
    sprintf (tmp, "HELP_%d=", (int) i);
    if (GetServiceParam (list, tmp, title, sizeof (title) - 1)) {
      upp->help = StringSave (title);
    } else {
      upp->help = StringSave (" ");
    }
    ValNodeAddPointer (&paramlist, i, (Pointer) upp);
    i++;
    sprintf (tmp, "PARAM_%d=", (int) i);
  }
  return paramlist;
}

static void ReadAnalysisConfigFile (CharPtr sect, MenU m, BaseFormPtr bfp,
                                    Boolean bspviewOK, Boolean docsumOK,
                                    NewObjectPtr PNTR head)

{
  Boolean     demomode = FALSE;
  Int2        format = 1;
  Char        host [128];
  Boolean     nucOK = FALSE;
  Boolean     onlyBspTarget = FALSE;
  ValNodePtr  paramlist = NULL;
  Char        program [128];
  Char        path [256];
  Uint2       port = 80;
  Char        prefix [128];
  Boolean     protOK = FALSE;
  Char        submenu [128];
  Char        suffix [128];
  Uint2       timeoutsec = 30;
  Char        title [128];
  Char        tmp [32];
  unsigned    int  val;

  if (! GetAppParam ("SEQNCGIS", sect, "TITLE", NULL, title, sizeof (title) - 1)) {
    StringNCpy_0 (title, sect, sizeof (title));
  }
  if (GetAppParam ("SEQNCGIS", sect, "HOST", NULL, host, sizeof (host) - 1)) {
    if (GetAppParam ("SEQNCGIS", sect, "FLAGS", NULL, tmp, sizeof (tmp) - 1)) {
      if (StringStr (tmp, "SEQ") == NULL) {
        bspviewOK = FALSE;
      }
      if (StringStr (tmp, "DOC") == NULL) {
        docsumOK = FALSE;
      }
      if (StringStr (tmp, "NUC") != NULL) {
        nucOK = TRUE;
      }
      if (StringStr (tmp, "PRT") != NULL) {
        protOK = TRUE;
      }
      if (StringStr (tmp, "TRG") != NULL) {
        onlyBspTarget = TRUE;
      }
    }

    if ((! bspviewOK) && (! docsumOK)) return;

    if (GetAppParam ("SEQNCGIS", sect, "PORT", NULL, tmp, sizeof (tmp) - 1) && 
        sscanf (tmp, "%u", &val) == 1) {
      port = (Uint2) val;
    } else {
      port = 80;
    }
    if (GetAppParam ("SEQNCGIS", sect, "FORMATIN", NULL, tmp, sizeof (tmp) - 1)) {
      if (StringICmp (tmp, "FASTA") == 0) {
        format = 1;
      } else if (StringICmp (tmp, "ASN.1") == 0) {
        format = 2;
      }
    }
    if (GetAppParam ("SEQNCGIS", sect, "TIMEOUT", NULL, tmp, sizeof (tmp) - 1) && 
        sscanf (tmp, "%u", &val) == 1) {
      timeoutsec = (Uint2) val;
    } else {
      timeoutsec = 30;
    }
    submenu [0] = '\0';
    GetAppParam ("SEQNCGIS", sect, "SUBMENU", NULL, submenu, sizeof (submenu) - 1);
    if (GetAppParam ("SEQNCGIS", sect, "DEMO", NULL, tmp, sizeof (tmp) - 1)) {
      if (StringICmp (tmp, "TRUE") == 0) {
        demomode = TRUE;
      }
    }

    if (GetAppParam ("SEQNCGIS", sect, "PATH", NULL, path, sizeof (path) - 1)) {
      if (GetAppParam ("SEQNCGIS", sect, "PROGRAM", NULL, program, sizeof (program) - 1)) {
        paramlist = GetConfigParamAndPromptLists (sect);
        prefix [0] = '\0';
        GetAppParam ("SEQNCGIS", sect, "PREFIX", NULL, prefix, sizeof (prefix) - 1);
        suffix [0] = '\0';
        GetAppParam ("SEQNCGIS", sect, "SUFFIX", NULL, suffix, sizeof (suffix) - 1);
        AddAnalysisItem (m, bfp, bspviewOK, docsumOK,
                         nucOK, protOK, onlyBspTarget,
                         host, port, path, program, timeoutsec, format, demomode,
                         SequinHandleURLResults, paramlist, prefix, suffix,
                         title, submenu, DoURLProc, head);
      }
    }
  }
}

static void ReadServiceConfigFile (CharPtr pathbase, ValNodePtr config,
                                   MenU m, BaseFormPtr bfp,
                                   Boolean bspviewOK, Boolean docsumOK,
                                   NewObjectPtr PNTR head)

{
  Char          ch;
  Boolean       demomode = FALSE;
  Int2          format = 1;
  FILE          *fp;
  Boolean       goOn = TRUE;
  Char          host [128];
  Boolean       keepGoing;
  ValNodePtr    list = NULL;
  Boolean       nucOK = FALSE;
  Boolean       onlyBspTarget = FALSE;
  ValNodePtr    paramlist = NULL;
  Char          program [128];
  Char          path [PATH_MAX];
  Uint2         port = 80;
  Char          prefix [128];
  Boolean       protOK = FALSE;
  CharPtr       ptr;
  Boolean       seenBracket;
  Char          str [256];
  Char          submenu [128];
  Char          suffix [128];
  Uint2         timeoutsec = 30;
  Char          title [128];
  Char          tmp [32];
  unsigned int  val;
  ValNodePtr    vnp;

  if (path == NULL || config == NULL || config->data.ptrvalue == NULL) return;
  StringNCpy_0 (path, pathbase, sizeof (path));
  FileBuildPath (path, NULL, (CharPtr) config->data.ptrvalue);
  fp = FileOpen (path, "r");
  if (fp == NULL) return;
  while (FileGets (str, sizeof (str), fp) != NULL) {
    ptr = str;
    ch = *ptr;
    while (ch != '\0' && ch != '\n' && ch != '\r') {
      ptr++;
      ch = *ptr;
    }
    *ptr = '\0';
    ValNodeCopyStr (&list, 1, str);
  }
  FileClose (fp);
  while (goOn) {
    goOn = FALSE;
    title [0] = '\0';
    if (GetServiceParam (list, "TITLE=", tmp, sizeof (tmp) - 1)) {
      StringNCpy_0 (title, tmp, sizeof (title));
    }
    if (StringHasNoText (title)) {
      if (GetServiceParam (list, "[", title, sizeof (title) - 1)) {
        ptr = StringChr (title, ']');
        if (ptr != NULL) {
          *ptr = '\0';
        }
      }
    }
    if (title [0] != '\0' && GetServiceParam (list, "HOST=", host, sizeof (host) - 1)) {
      if (GetServiceParam (list, "FLAGS=", tmp, sizeof (tmp) - 1)) {
        if (StringStr (tmp, "SEQ") == NULL) {
          bspviewOK= FALSE;
        }
        if (StringStr (tmp, "DOC") == NULL) {
          docsumOK= FALSE;
        }
        if (StringStr (tmp, "NUC") != NULL) {
          nucOK= TRUE;
        }
        if (StringStr (tmp, "PRT") != NULL) {
          protOK= TRUE;
        }
        if (StringStr (tmp, "TRG") != NULL) {
          onlyBspTarget= TRUE;
        }
      }

      if (bspviewOK || docsumOK) {

        if (GetServiceParam (list, "PORT=", tmp, sizeof (tmp) - 1) && 
            sscanf (tmp, "%u", &val) == 1) {
          port = (Uint2) val;
        } else {
          port = 80;
        }
        if (GetServiceParam (list, "FORMATIN=", tmp, sizeof (tmp) - 1)) {
          if (StringICmp (tmp, "FASTA") == 0) {
            format = 1;
          } else if (StringICmp (tmp, "ASN.1") == 0) {
            format = 2;
          }
        }
       if (GetServiceParam (list, "TIMEOUT=", tmp, sizeof (tmp) - 1) && 
            sscanf (tmp, "%u", &val) == 1) {
          timeoutsec = (Uint2) val;
        } else {
          timeoutsec = 30;
        }
        submenu [0] = '\0';
        GetServiceParam (list, "SUBMENU=", submenu, sizeof (submenu) - 1);
        if (GetServiceParam (list, "DEMO=", tmp, sizeof (tmp) - 1)) {
          if (StringICmp (tmp, "TRUE") == 0) {
            demomode = TRUE;
          }
        }

        if (GetServiceParam (list, "PATH=", path, sizeof (path) - 1)) {
          if (GetServiceParam (list, "PROGRAM=", program, sizeof (program) - 1)) {
            paramlist = GetServiceParamAndPromptLists (list);
            prefix [0] = '\0';
            GetServiceParam (list, "PREFIX=", prefix, sizeof (prefix) - 1);
            suffix [0] = '\0';
            GetServiceParam (list, "SUFFIX=", suffix, sizeof (suffix) - 1);
            AddAnalysisItem (m, bfp, bspviewOK, docsumOK,
                             nucOK, protOK, onlyBspTarget,
                             host, port, path, program, timeoutsec, format, demomode,
                             SequinHandleURLResults, paramlist, prefix, suffix,
                             title, submenu, DoURLProc, head);
          }
        }

      }
    }

    seenBracket = FALSE;
    keepGoing = TRUE;
    for (vnp = list; vnp != NULL && keepGoing; vnp = vnp->next) {
      ptr = (CharPtr) vnp->data.ptrvalue;
      if (ptr != NULL) {
        if (ptr [0] == '[') {
          if (seenBracket) {
            keepGoing = FALSE;
          } else {
            seenBracket = TRUE;
          }
        }
        if (keepGoing) {
          vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
        }
      }
    }

  }

  ValNodeFreeData (list);
}

extern MenU CreateAnalysisMenu (WindoW w, BaseFormPtr bfp, Boolean bspviewOK, Boolean docsumOK)

{
  NewObjectPtr  first;
  ValNodePtr    head1 = NULL, head2 = NULL;
  Int2          i;
  size_t        len;
  MenU          m;
  Char          path1 [PATH_MAX];
  Char          path2 [PATH_MAX];
  CharPtr       ptr;
  SbstrucPtr    sbp;
  Char          sect [256];
  Char          temp [32];
  ValNodePtr    vnp;

  ProgramPath (path1, sizeof (path1));
  ptr = StringRChr (path1, DIRDELIMCHR);
  if (ptr != NULL) {
    ptr++;
    *ptr = '\0';
  }
  FileBuildPath (path1, "services", NULL);
  head1 = DirCatalog (path1);

  if (GetAppParam ("SEQNCGIS", "SERVICES", "PATH", NULL, path2, sizeof (path2) - 1)) {
    len = StringLen (path2);
    if (path2 [len - 1] != DIRDELIMCHR) {
      StringCat (path2, DIRDELIMSTR);
    }
    if (StringCmp (path1, path2) != 0) {
      head2 = DirCatalog (path2);
    }
  }

  if ((! extraServices) && (! indexerVersion) && (! genomeCenter) &&
      head1 == NULL && head2 == NULL) {
    if (! GetAppParam ("SEQNCGIS", "ORDER", NULL, NULL, sect, sizeof (sect) - 1)) {
      return NULL;
    }
  }
  m = PulldownMenu (w, "Analysis");
  if (m == NULL) return NULL;
  analysissubmenulist = NULL;
  first = NULL;
  if (bspviewOK) {
    AddAnalysisItem (m, bfp, bspviewOK, FALSE, TRUE, FALSE, TRUE,
                     NULL, 0, NULL, NULL, 0, 0, FALSE, NULL, NULL, NULL, NULL,
                     "Restriction Search", "Search",
                     SimpleRsiteProc, &first);
    if (indexerVersion) {
      AddAnalysisItem (m, bfp, bspviewOK, FALSE, TRUE, FALSE, TRUE,
                       NULL, 0, NULL, NULL, 0, 0, FALSE, NULL, NULL, NULL, NULL,
                       "QBlast Test", "Search",
                       SimpleQBlastProc, &first);
    }
  }
  if (bspviewOK || docsumOK) {
    if (useEntrez) {
      i = 1;
      sprintf (temp, "ORDER_%d", (int) i);
      while (GetAppParam ("SEQNCGIS", "ORDER", temp, NULL, sect, sizeof (sect) - 1)) {
        ReadAnalysisConfigFile (sect, m, bfp, bspviewOK, docsumOK, &first);
        i++;
        sprintf (temp, "ORDER_%d", (int) i);
      }
      for (vnp = head1; vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == 0) {
          ReadServiceConfigFile (path1, vnp, m, bfp, bspviewOK, docsumOK, &first);
        }
      }
      for (vnp = head2; vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == 0) {
          ReadServiceConfigFile (path2, vnp, m, bfp, bspviewOK, docsumOK, &first);
        }
      }
    }
  }
  if (bspviewOK) {
  }
  if (docsumOK) {
  }
#ifdef WIN_MAC
  macUserDataPtr = LinkNewObjectLists (macUserDataPtr, first);
#else
  bfp->userDataPtr = LinkNewObjectLists (bfp->userDataPtr, first);
#endif
  for (vnp = analysissubmenulist; vnp != NULL; vnp = vnp->next) {
    sbp = (SbstrucPtr) vnp->data.ptrvalue;
    if (sbp != NULL) {
      sbp->name = MemFree (sbp->name);
    }
  }
  analysissubmenulist = ValNodeFreeData (analysissubmenulist);
  ValNodeFreeData (head1);
  ValNodeFreeData (head2);
  return m;
}

/* NEW UPDATE SEQUENCE SECTION */


#define SQN_LEFT    1
#define SQN_RIGHT   2
#define SQN_MIDDLE  3

static Uint4 sqn_binary_search_on_uint4_list(Uint4Ptr list, Uint4 pos, Uint4 listlen)
{
   Uint4  L;
   Uint4  mid;
   Uint4  R;

   if (list == NULL || listlen == 0)
      return 0;
   L = 0;
   R = listlen - 1;
   while (L < R)
   {
      mid = (L+R)/2;
      if (list[mid + 1] <= pos)
      {
         L = mid + 1;
      } else
      {
         R = mid;
      }
   }
   return R;
}

static Int4 MapRowCoordsSpecial(SeqAlignPtr sap, Uint4 pos, Int4 row, Int4 which_end)
{
   DenseSegPtr  dsp;
   Int4         idx;
   Int4         offset;
   SAIndexPtr   saip;
   Int4         start;

   if (sap == NULL || row < 0)
      return -1;
   if (sap->saip == NULL)
      return -1;
   saip = (SAIndexPtr)sap->saip;
   dsp = (DenseSegPtr)sap->segs;
   start = sqn_binary_search_on_uint4_list(saip->aligncoords, pos, dsp->numseg);
   offset = pos - saip->aligncoords[start];
   idx = (dsp->dim*start) + row - 1;
   if (dsp->starts[idx] == -1)
   {
      if (which_end == SQN_RIGHT)
      {
         /* round down */
         while (start >= 0) {
            idx = (dsp->dim*start) + row - 1;
            if (dsp->starts[idx] != -1)
               return (dsp->starts[idx] + dsp->lens[start] - 1);
            start--;
         }
         return -2;
      } else if (which_end == SQN_LEFT)
      {
         /* round up */
         while (start < dsp->numseg) {
            idx = (dsp->dim*start) + row - 1;
            if (dsp->starts[idx] != -1)
               return (dsp->starts[idx]);
            start++;
         }
         return -2;
      }
   } else
   {
      idx = (dsp->dim*start) + row - 1;
      if (dsp->strands[idx] != Seq_strand_minus)
         return (dsp->starts[idx] + offset);
      else
         return (dsp->starts[idx] + dsp->lens[start] - 1 - offset);
   }
   return -1;
}

static Int4 MapBioseqToBioseqSpecial(SeqAlignPtr sap, Int4 begin, Int4 fin, Int4 pos, Int4 which_end)
{
   Int4  bspos;
   Int4  sapos;
   Int4  start1;
   Int4  start2;
   Int4  stop1;
   Int4  stop2;

   if (sap == NULL || sap->saip == NULL)
      return -2;
   AlnMgr2GetNthSeqRangeInSA(sap, begin, &start1, &stop1);
   AlnMgr2GetNthSeqRangeInSA(sap, fin, &start2, &stop2);
   /* check to see whether the position is outside the alignment */
   if (pos < start1)
      return (start2 - (start1 - pos));
   else if (pos > stop1)
      return (stop2 + (pos-stop1));
   sapos = AlnMgr2MapBioseqToSeqAlign(sap, pos, begin);
   bspos = MapRowCoordsSpecial(sap, sapos, fin, which_end);
   if (bspos >= 0)
      return bspos;
   else if (which_end == SQN_LEFT)
      return (start2-1);
   else if (which_end == SQN_RIGHT)
      return (stop2+1);
   else
      return 0;
}

static void ListPhrapGraphsCallback (SeqGraphPtr sgp, Pointer userdata)
{
  ValNodePtr PNTR vnpp;
  
  if (sgp == NULL || userdata == NULL) return;
  if (StringICmp (sgp->title, "Phrap Quality") == 0)
  {
    vnpp = (ValNodePtr PNTR) userdata;
    ValNodeAddPointer (vnpp, 0, sgp);
  }
}

/* THOUGHTS:
 * Can we/must we update quality scores before/after the old Bioseq has been replaced?
 * If we replace quality scores after the bioseq has been replaced, the oldbsp->length
 * is the length of the buffer we need to hold the quality scores,
 * otherwise use the newbsp->length.
 * Useful functions:
  aln_len = AlnMgr2GetAlnLength(salp, FALSE);

NLM_EXTERN Int4 AlnMgr2GetNumAlnBlocks(SeqAlignPtr sap)
NLM_EXTERN Boolean AlnMgr2GetNthBlockRange(SeqAlignPtr sap, Int4 n, Int4Ptr start, Int4Ptr stop)
  
  
 * Assumptions: data replacement has already taken place, oldbsp is in row 1 of salp,
 *              newbsp is in row 2 of salp.
 
 */
static Boolean 
ReplaceQualityScores
(BioseqPtr   oldbsp,
 BioseqPtr   newbsp,
 SeqAlignPtr salp,
 FILE        *log_fp,
 BoolPtr     data_in_log)
{
  ValNodePtr    oldhead = NULL, newhead = NULL, vnp;
  BytePtr       new_store = NULL;  
  Int4          cur_pos;
  SeqGraphPtr   sgp, sgp_list = NULL, last_sgp = NULL;
  ByteStorePtr  bs;
  Int4          len;
  Int4          graph_left, i, graph_len;
  SeqIntPtr     sintp;
  SeqAnnotPtr   sap, last_sap;
  Int2          min, max;
  Char          acc_str [256];
  
  if (oldbsp == NULL || newbsp == NULL || salp == NULL
      || !ISA_na (oldbsp->mol) || !ISA_na (newbsp->mol) || oldbsp->length >= MAXALLOC)
  {
    return FALSE;
  }
  
  VisitGraphsOnBsp (oldbsp, &oldhead, ListPhrapGraphsCallback);
  VisitGraphsOnBsp (newbsp, &newhead, ListPhrapGraphsCallback);
  
  /* prepare buffer to hold scores from both */
  len = oldbsp->length;
  new_store = MemNew (sizeof (Byte) * (len + 2));
  if (new_store == NULL) 
  {
    oldhead = ValNodeFreeData (oldhead);
    newhead = ValNodeFreeData (newhead);
    return FALSE;
  }
  /* init byte store */
  for (cur_pos = 0; cur_pos < len; cur_pos ++)
  {
    new_store [cur_pos] = 255;
  }

  /* copy in old scores */
  for (vnp = oldhead; vnp != NULL; vnp = vnp->next) 
  {
    sgp = vnp->data.ptrvalue;
    if (sgp != NULL)
    {
      sgp->idx.deleteme = TRUE;
    }
  }

  if (oldhead != NULL && newhead == NULL && log_fp != NULL && data_in_log != NULL)
  {
    SeqIdWrite (newbsp->id, acc_str, PRINTID_REPORT, sizeof (acc_str));
    fprintf (log_fp, "Quality scores cleared for %s\n", acc_str);
    *data_in_log = TRUE;
  }
  oldhead = ValNodeFree (oldhead);
  
  /* now copy in new scores */
  for (vnp = newhead; vnp != NULL; vnp = vnp->next) {
    sgp = vnp->data.ptrvalue;
    bs = (ByteStorePtr) sgp->values;
    BSSeek (bs, 0, SEEK_SET);
    cur_pos = GetOffsetInBioseq (sgp->loc, newbsp, SEQLOC_LEFT_END);
    for (i = 0; i < sgp->numval && cur_pos < len; i++, cur_pos++)
    {
      new_store [cur_pos] = (Byte) BSGetByte (bs);        
    }
  }
  
  newhead = ValNodeFree (newhead);

  /* Now we have a Byte array that contains the quality scores.
   * Time to create graphs from the Byte array.
   */
  i = 0;
  graph_left = -1;
  while (i < len)
  {
    if (new_store [i] == (Byte)-1)
    {
      if (graph_left > -1)
      {
        /* add new SeqGraph to list */
        graph_len = i - graph_left;
        bs = BSNew (1000);
        if (bs != NULL)
        {
          BSSeek (bs, 0, SEEK_SET);
          BSWrite (bs, new_store + graph_left, graph_len);
          sgp = SeqGraphNew ();
          sgp->numval = BSLen (bs);
          BSPutByte (bs, EOF);
          sgp->title = StringSave ("Phrap Quality");
          if (len != sgp->numval) {
            sgp->flags [0] = 1;
            sgp->compr = (len) / sgp->numval;
          } else {
            sgp->flags [0] = 0;
            sgp->compr = 1;
          }
          sgp->flags [1] = 0;
          sgp->flags [2] = 3;
          sgp->axis.intvalue = 0;
          sgp->min.intvalue = min;
          sgp->max.intvalue = max;
          sgp->a = 1.0;
          sgp->b = 0;
          sgp->values = (Pointer) bs;

          sintp = SeqIntNew ();
          sintp->from = 0;
          sintp->to = len - 1;
          sintp->id = SeqIdDup (oldbsp->id);
          ValNodeAddPointer (&(sgp->loc), SEQLOC_INT, (Pointer) sintp);
          if (last_sgp == NULL)
          {
            sgp_list = sgp;
          }
          else
          {
            last_sgp->next = sgp;
          }
          last_sgp = sgp;
        }
        graph_left = -1;
      }
    }
    else
    {
      if (graph_left == -1)
      {
        graph_left = i;
        min = new_store [i];
        max = min;
      }
      else
      {
        min = MIN (min, new_store[i]);
        max = MAX (max, new_store[i]);
      }
    }
    i++;
  }
  if (graph_left > -1)
  {
    /* add new SeqGraph to list */
    graph_len = i - graph_left;
    bs = BSNew (1000);
    if (bs != NULL)
    {
      BSSeek (bs, 0, SEEK_SET);
      BSWrite (bs, new_store + graph_left, graph_len);
      sgp = SeqGraphNew ();
      sgp->numval = BSLen (bs);
      BSPutByte (bs, EOF);
      sgp->title = StringSave ("Phrap Quality");
      if (len != sgp->numval) {
        sgp->flags [0] = 1;
        sgp->compr = (len) / sgp->numval;
      } else {
        sgp->flags [0] = 0;
        sgp->compr = 1;
      }
      sgp->flags [1] = 0;
      sgp->flags [2] = 3;
      sgp->axis.intvalue = 0;
      sgp->min.intvalue = min;
      sgp->max.intvalue = max;
      sgp->a = 1.0;
      sgp->b = 0;
      sgp->values = (Pointer) bs;

      sintp = SeqIntNew ();
      sintp->from = 0;
      sintp->to = len - 1;
      sintp->id = SeqIdDup (oldbsp->id);
      ValNodeAddPointer (&(sgp->loc), SEQLOC_INT, (Pointer) sintp);
      if (last_sgp == NULL)
      {
        sgp_list = sgp;
      }
      else
      {
        last_sgp->next = sgp;
      }
      last_sgp = sgp;
    }
  }
  
  if (sgp_list != NULL)
  {
    /* now add in new phrap annotations */
    last_sap = NULL;
    last_sgp = NULL;
    for (sap = oldbsp->annot; sap != NULL && sap->type !=3; sap = sap->next) 
    {
      last_sap = sap;
    }
    
    if (sap == NULL)
    {
      sap = SeqAnnotNew ();
      sap->type = 3;
      sap->data = sgp_list;
      if (last_sap == NULL)
      {
        oldbsp->annot = sap;
      }
      else
      {
        last_sap->next = sap;
      }
    }
    else
    {
      for (sgp = (SeqGraphPtr) sap->data; sgp != NULL; sgp = sgp->next) 
      {
        last_sgp = sgp;  
      }
      if (last_sgp == NULL)
      {
        sap->data = sgp_list;
      }
      else
      {
        last_sgp->next = sgp_list;
      }
    }
    
  }

  /* now remove old phrap annotations */
  DeleteMarkedObjects (0, OBJ_BIOSEQ, (Pointer) oldbsp);
  return TRUE;
}


static Boolean AdjustAlignment (
  UpsDataPtr udp,
  Int2 choice
)

{
  DenseSegPtr  dsp;
  Int2         j;
  SeqAlignPtr  sap;

  if (udp == NULL) return FALSE;

  sap = udp->salp;
  if (sap == NULL) return FALSE;
  AMFreeAllIndexes (sap);

  if (sap->segtype == SAS_DENSEG) {
    dsp = (DenseSegPtr) sap->segs;

    switch (choice) {
      case 2 :
        /* adjust alignment 5' */
        if (dsp != NULL && dsp->lens != NULL && dsp->numseg > 0) {
          dsp->lens [dsp->numseg - 1] += udp->old3;
        }
        break;
      case 3 :
        /* adjust alignment 3' */
        if (dsp != NULL && dsp->lens != NULL && dsp->starts != NULL && dsp->numseg > 0) {
          dsp->lens [0] += udp->old5;
          dsp->starts [0] = 0;
          dsp->starts [1] = 0;
          for (j = 1; j < dsp->numseg; j++) {
            if (dsp->starts [1 + j * 2] != -1) {
              dsp->starts [1 + j * 2] += udp->old5 - udp->new5;
            }
          }
        }
        break;
      case 4 :
        /* adjust alignment patch */
        if (dsp != NULL && dsp->lens != NULL && dsp->starts != NULL && dsp->numseg > 0) {
          dsp->lens [dsp->numseg - 1] += udp->old3;
          dsp->lens [0] += udp->old5;
          dsp->starts [0] = 0;
          dsp->starts [1] = 0;
          for (j = 1; j < dsp->numseg; j++) {
            if (dsp->starts [1 + j * 2] != -1) {
              dsp->starts [1 + j * 2] += udp->old5 - udp->new5;
            }
          }
        }
        break;
      default :
        break;
    }
  }

  AlnMgr2IndexSingleChildSeqAlign (sap);

  return TRUE;
}

static void OffsetLoc (SeqLocPtr slp, Int4 offset, SeqIdPtr sip)

{
  PackSeqPntPtr  psp;
  SeqIntPtr      sinp;
  SeqPntPtr      spp;
  Uint1          used;

  if (slp == NULL) return;
  switch (slp->choice) {
    case SEQLOC_INT :
      sinp = (SeqIntPtr) slp->data.ptrvalue;
      if (sinp != NULL) {
        sinp->from += offset;
        sinp->to += offset;
        if (sip != NULL) {
          sinp->id = SeqIdFree (sinp->id);
          sinp->id = SeqIdDup (sip);
        }
      }
      break;
    case SEQLOC_PNT :
      spp = (SeqPntPtr) slp->data.ptrvalue;
      if (spp != NULL) {
        spp->point += offset;
        if (sip != NULL) {
          spp->id = SeqIdFree (spp->id);
          spp->id = SeqIdDup (sip);
        }
      }
      break;
    case SEQLOC_PACKED_PNT :
      psp = (PackSeqPntPtr) slp->data.ptrvalue;
      if (psp != NULL) {
        for (used = 0; used < psp->used; used++) {
          psp->pnts [used] += offset;
        }
        if (sip != NULL) {
          psp->id = SeqIdFree (psp->id);
          psp->id = SeqIdDup (sip);
        }
      }
      break;
    default :
      break;
  }
}

extern void OffsetLocation (SeqLocPtr loc, Int4 offset, SeqIdPtr sip)

{
  SeqLocPtr  slp;

  slp = SeqLocFindNext (loc, NULL);
  while (slp != NULL) {
    OffsetLoc (slp, offset, sip);
    slp = SeqLocFindNext (loc, slp);
  }
}

static void PromoteSeqId (SeqIdPtr sip, Pointer userdata)

{
  SeqIdPtr  bestid, newid, oldid;

  bestid = (SeqIdPtr) userdata;

  newid = SeqIdDup (bestid);
  if (newid == NULL) return;

  oldid = ValNodeNew (NULL);
  if (oldid == NULL) return;

  MemCopy (oldid, sip, sizeof (ValNode));
  oldid->next = NULL;

  sip->choice = newid->choice;
  sip->data.ptrvalue = newid->data.ptrvalue;

  SeqIdFree (oldid);
  ValNodeFree (newid);

  SeqIdStripLocus (sip);
}

static void CorrectFeatureSeqIds (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  VisitSeqIdsInSeqLoc (sfp->location, userdata, PromoteSeqId);
}

static Boolean DoFeaturePropWithOffset (
  UpsDataPtr udp,
  Int4 offset,
  SeqAnnotPtr PNTR sapp,
  Boolean patch
)

{
  BioseqPtr          bsp, newbsp, oldbsp;
  CodeBreakPtr       cbp;
  SeqMgrFeatContext  context;
  CdRegionPtr        crp;
  SeqFeatPtr         dup, sfp, last = NULL;
  Uint2              entityID;
  Boolean            keepProteinIDs;
  SeqEntryPtr        newsep, prdsep, top;
  RnaRefPtr          rrp;
  SeqAnnotPtr        sap = NULL, saptmp;
  SeqDescrPtr        sdp;
  SeqIdPtr           sip;
  tRNAPtr            trp;

  if (udp == NULL) return FALSE;

  SeqEntrySetScope (NULL);

  sfp = SeqMgrGetNextFeature (udp->newbsp, NULL, 0, 0, &context);
  if (sfp == NULL) return FALSE;

  if (udp->diffOrgs) {
    keepProteinIDs = FALSE;
  } else {
    keepProteinIDs = GetStatus (udp->keepProteinIDs);
  }

  oldbsp = udp->oldbsp;

  entityID = ObjMgrGetEntityIDForPointer (udp->oldbsp);
  top = GetBestTopParentForData (entityID, udp->oldbsp);

  sdp = ExtractBioSourceAndPubs (top);

  sip = SeqIdFindBest (oldbsp->id, 0);

  while (sfp != NULL) {

    if ((! patch) || (context.right >= udp->new5 && context.left <= udp->new5 + udp->newa)) {

    dup = AsnIoMemCopy ((Pointer) sfp,
                        (AsnReadFunc) SeqFeatAsnRead,
                        (AsnWriteFunc) SeqFeatAsnWrite);

    if (last == NULL) {
      sap = SeqAnnotNew ();
      if (oldbsp->annot == NULL) {
        oldbsp->annot = sap;
      } else {
        for (saptmp = oldbsp->annot; saptmp->next != NULL; saptmp = saptmp->next) continue;
        saptmp->next = sap;
      }
      sap->type = 1;
      sap->data = (Pointer) dup;
    } else {
      last->next = dup;
    }
    last = dup;

    /*
    sep = SeqMgrGetSeqEntryForData (oldbsp);
    CreateNewFeature (sep, NULL, dup->data.choice, dup);
    */

    OffsetLocation (dup->location, offset, sip);
    switch (dup->data.choice) {
      case SEQFEAT_CDREGION :
        crp = (CdRegionPtr) dup->data.value.ptrvalue;
        if (crp != NULL) {
          for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
            OffsetLocation (cbp->loc, offset, sip);
          }
        }
        break;
      case SEQFEAT_RNA :
        rrp = (RnaRefPtr) dup->data.value.ptrvalue;
        if (rrp != NULL && rrp->ext.choice == 2) {
          trp = (tRNAPtr) rrp->ext.value.ptrvalue;
          if (trp != NULL && trp->anticodon != NULL) {
            OffsetLocation (trp->anticodon, offset, sip);
          }
        }
        break;
      default :
        break;
    }
    if (dup->product != NULL) {
      SeqEntrySetScope (NULL);
      bsp = BioseqFindFromSeqLoc (dup->product);
      if (bsp != NULL) {
        prdsep = SeqMgrGetSeqEntryForData (bsp);
        if (prdsep != NULL) {
          newsep = AsnIoMemCopy ((Pointer) prdsep,
                                 (AsnReadFunc) SeqEntryAsnRead,
                                 (AsnWriteFunc) SeqEntryAsnWrite);
          if (newsep != NULL) {
            if (IS_Bioseq (newsep)) {
              newbsp = (BioseqPtr) newsep->data.ptrvalue;
              if (newbsp != NULL) {
                if (! keepProteinIDs) {
                  newbsp->id = SeqIdSetFree (newbsp->id);
                  newbsp->id = MakeNewProteinSeqId (NULL, sip);
                  newbsp->hist = SeqHistFree (newbsp->hist);
                  VisitFeaturesOnBsp (newbsp, (Pointer) newbsp->id, CorrectFeatureSeqIds);
                  SetSeqFeatProduct (dup, newbsp);
                  /*
                  dup->product = SeqLocFree (dup->product);
                  dup->product = CreateWholeInterval (newsep);
                  */
                }
                SeqMgrReplaceInBioseqIndex (newbsp);
              }
            }
            AddSeqEntryToSeqEntry (top, newsep, TRUE);
          }
        }
      }
    }
    }

    sfp = SeqMgrGetNextFeature (udp->newbsp, sfp, 0, 0, &context);
  }

  ReplaceBioSourceAndPubs (top, sdp);

  if (sapp != NULL) {
    *sapp = sap;
  }

  return TRUE;
}

/* This function adjusts the endpoints of a location, as long as the
 * endpoints are in the area represented by the alignment.
 * When we are adjusting locations for an alignment of a part, we will 
 * be looking at all features indexed on the main segment, but we only
 * want to adjust feature endpoints located on the segment that we are
 * updating.
 */
static Int4 AdjustEndpoint 
(SeqAlignPtr salp,
 SeqLocPtr   slp, 
 Int4        max_length, 
 Int4        begin, 
 Int4        fin,
 Int4        endpoint,
 Int4        end)
{
  BioseqPtr            slp_bsp, parent_bsp, old_bsp;
  SeqMgrSegmentContext segcontext; 
  SeqIdPtr             old_sip; 
  Int4                 pt;
  
  if (slp == NULL || salp == NULL)
  {
    return endpoint;
  }
  
  old_sip = AlnMgr2GetNthSeqIdPtr (salp, begin);
  old_bsp = BioseqFind (old_sip);

  parent_bsp = SeqMgrGetParentOfPart (old_bsp, &segcontext);
  
  slp_bsp = BioseqFind (SeqLocId (slp));
  if (slp_bsp == old_bsp
      || (slp_bsp == parent_bsp 
          && endpoint >= segcontext.cumOffset + segcontext.from
          && endpoint < segcontext.cumOffset + segcontext.to))
  {
    if (slp_bsp == parent_bsp)
    {
      endpoint -= segcontext.cumOffset + segcontext.from;
    }
    pt = MapBioseqToBioseqSpecial (salp, begin, fin, endpoint, end);
    if (pt < 0) {
      pt = 0;
    } else if (pt >= max_length) {
      pt = max_length - 1;
    }
    if (slp_bsp == parent_bsp)
    {
      pt += segcontext.cumOffset + segcontext.from;
    }
  }
  else
  {
    pt = endpoint;
  }
  
  return pt;  
}

static void ReplaceLocation (SeqAlignPtr salp, SeqLocPtr slp, Int4 length, Int4 begin, Int4 fin)

{
  PackSeqPntPtr  psp;
  SeqIntPtr      sinp;
  SeqPntPtr      spp;
  Uint1          used;

  if (slp == NULL) return;
  switch (slp->choice) {
    case SEQLOC_INT :
      sinp = (SeqIntPtr) slp->data.ptrvalue;
      if (sinp != NULL) {
        sinp->from = AdjustEndpoint (salp, slp, length, begin, fin,
                                     sinp->from, SQN_LEFT);
        sinp->to = AdjustEndpoint (salp, slp, length, begin, fin,
                                   sinp->to, SQN_RIGHT);
      }
      break;
    case SEQLOC_PNT :
      spp = (SeqPntPtr) slp->data.ptrvalue;
      if (spp != NULL) {
        spp->point = AdjustEndpoint (salp, slp, length, begin, fin,
                                     spp->point, SQN_LEFT);
      }
      break;
    case SEQLOC_PACKED_PNT :
      psp = (PackSeqPntPtr) slp->data.ptrvalue;
      if (psp != NULL) {
        for (used = 0; used < psp->used; used++) {
          psp->pnts [used] = AdjustEndpoint (salp, slp, length, begin, fin,
                                             psp->pnts [used], SQN_LEFT);
        }
      }
      break;
    default :
      break;
  }
}

/* this function iterates through the pieces of a complex location
 * and calls ReplaceLocation for each one.  ReplaceLocation will only
 * act on SEQLOC_INT, SEQLOC_PNT, and SEQLOC_PACKED_PNT and will ignore
 * other types.
 */
static void 
ReplaceComplexLocation 
(SeqLocPtr   slp,
 SeqAlignPtr salp,
 Int4        new_len,
 Int4        begin,
 Int4        fin)
{
  SeqLocPtr subslp;
  
  if (slp == NULL || salp == NULL)
  {
    return;
  }
  
  subslp = SeqLocFindNext (slp, NULL);
  while (subslp != NULL) {
    ReplaceLocation (salp, subslp, new_len, begin, fin);
    subslp = SeqLocFindNext (slp, subslp);
  }
}

static Int4 LengthForNewSequence (UpsDataPtr udp)
{
  Int4 new_len = 0;
  
  if (udp == NULL)
  {
    return 0;
  }
  else if (GetValue (udp->sfb) == UPDATE_FEATURES_ONLY)
  {
    new_len = udp->oldbsp->length;
  }
  else if (udp->rmcval == UPDATE_PATCH)
  {
    new_len = udp->old5 + udp->newa + udp->old3;
  }
  else if (udp->rmcval == UPDATE_REPLACE)
  {
    new_len = udp->new5 + udp->newa + udp->new3;
  }
  else if (udp->rmcval == UPDATE_EXTEND5)
  {
    new_len = udp->new5 + udp->olda + udp->old3;
  }
  else if (udp->rmcval == UPDATE_EXTEND3)
  {
    new_len = udp->old5 + udp->olda + udp->new3;
  }
  return new_len;
}

static SeqLocPtr 
GetPropagatedLocation 
(SeqLocPtr   orig_loc, 
 BioseqPtr   newbsp,
 BioseqPtr   oldbsp,
 Int4        new_len,
 SeqAlignPtr salp)
{
  SeqLocPtr tmp_loc, new_loc;
  Boolean   split;
  
  tmp_loc = SeqLocCopy (orig_loc);
  ReplaceComplexLocation (tmp_loc, salp, new_len, 2, 1);   

  new_loc = SeqLocCopyRegion (oldbsp->id, tmp_loc, newbsp, 0, oldbsp->length - 1, Seq_strand_plus, &split);
  
  tmp_loc = SeqLocFree (tmp_loc);
  
  return new_loc;
}

static Boolean DoFeaturePropThruAlign (
  UpsDataPtr udp,
  SeqAnnotPtr PNTR sapp
)

{
  BioseqPtr          bsp, newbsp, oldbsp;
  CodeBreakPtr       cbp, prevcbp, nextcbp;
  SeqMgrFeatContext  context;
  CdRegionPtr        crp;
  SeqFeatPtr         dup, sfp, last = NULL;
  Uint2              entityID;
  Int4               from, to;
  Boolean            keepProteinIDs;
  SeqLocPtr          newloc;
  SeqEntryPtr        newsep, prdsep, top;
  RnaRefPtr          rrp;
  SeqAnnotPtr        sap = NULL, saptmp;
  SeqDescrPtr        sdp;
  SeqIdPtr           sip;
  Boolean            split;
  tRNAPtr            trp;
  Boolean            partial5, partial3;

  if (udp == NULL) return FALSE;

  SeqEntrySetScope (NULL);

  sfp = SeqMgrGetNextFeature (udp->newbsp, NULL, 0, 0, &context);
  if (sfp == NULL) return FALSE;

  keepProteinIDs = GetStatus (udp->keepProteinIDs);

  oldbsp = udp->oldbsp;

  entityID = ObjMgrGetEntityIDForPointer (oldbsp);
  top = GetBestTopParentForData (entityID, oldbsp);

  sdp = ExtractBioSourceAndPubs (top);

  sip = SeqIdFindBest (oldbsp->id, 0);

  from = udp->new5;
  to = udp->new5 + udp->newa;

  while (sfp != NULL) {

    if (context.right >= from && context.left <= to) {
      split = FALSE;
      newloc = GetPropagatedLocation (sfp->location, udp->newbsp, udp->oldbsp, 
                                      LengthForNewSequence (udp), udp->salp);
      if (newloc != NULL) {
        CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
        SetSeqLocPartial (newloc, partial5, partial3);
        dup = AsnIoMemCopy ((Pointer) sfp,
                            (AsnReadFunc) SeqFeatAsnRead,
                            (AsnWriteFunc) SeqFeatAsnWrite);
                            
        SeqLocFree (dup->location);
        dup->location = newloc;
        if (split) {
          dup->partial = TRUE;
        }
        dup->partial |= partial5;
        dup->partial |= partial3;

        if (last == NULL) {
          sap = SeqAnnotNew ();
          if (oldbsp->annot == NULL) {
            oldbsp->annot = sap;
          } else {
            for (saptmp = oldbsp->annot; saptmp->next != NULL; saptmp = saptmp->next) continue;
            saptmp->next = sap;
          }
          sap->type = 1;
          sap->data = (Pointer) dup;
        } else {
          last->next = dup;
        }
        last = dup;

        switch (dup->data.choice) {
          case SEQFEAT_CDREGION :
            crp = (CdRegionPtr) dup->data.value.ptrvalue;
            if (crp != NULL) {
              prevcbp = NULL;
              for (cbp = crp->code_break; cbp != NULL; cbp = nextcbp) {
                nextcbp = cbp->next;
                newloc = GetPropagatedLocation (cbp->loc, udp->newbsp, udp->oldbsp, 
                                                LengthForNewSequence (udp), udp->salp);
                SeqLocFree (cbp->loc);
                cbp->loc = newloc;
                if (cbp->loc == NULL) {
                  if (prevcbp != NULL) {
                    prevcbp->next = nextcbp;
                  } else {
                    crp->code_break = nextcbp;
                  }
                  cbp->next = NULL;
                  CodeBreakFree (cbp);
                } else {
                  prevcbp = cbp;
                }
              }
            }
            break;
          case SEQFEAT_RNA :
            rrp = (RnaRefPtr) dup->data.value.ptrvalue;
            if (rrp != NULL && rrp->ext.choice == 2) {
              trp = (tRNAPtr) rrp->ext.value.ptrvalue;
              if (trp != NULL && trp->anticodon != NULL) {
                newloc = GetPropagatedLocation (trp->anticodon, udp->newbsp, udp->oldbsp, 
                                                LengthForNewSequence (udp), udp->salp);
                SeqLocFree (trp->anticodon);
                trp->anticodon = newloc;
              }
            }
            break;
          default :
            break;
        }
        if (dup->product != NULL) {
          SeqEntrySetScope (NULL);
          bsp = BioseqFindFromSeqLoc (dup->product);
          if (bsp != NULL) {
            prdsep = SeqMgrGetSeqEntryForData (bsp);
            if (prdsep != NULL) {
              newsep = AsnIoMemCopy ((Pointer) prdsep,
                                     (AsnReadFunc) SeqEntryAsnRead,
                                     (AsnWriteFunc) SeqEntryAsnWrite);
              if (newsep != NULL) {
                if (IS_Bioseq (newsep)) {
                  newbsp = (BioseqPtr) newsep->data.ptrvalue;
                  if (newbsp != NULL) {
                    if (! keepProteinIDs) {
                      newbsp->id = SeqIdSetFree (newbsp->id);
                      newbsp->id = MakeNewProteinSeqId (NULL, sip);
                      VisitFeaturesOnBsp (newbsp, (Pointer) newbsp->id, CorrectFeatureSeqIds);
                      SetSeqFeatProduct (dup, newbsp);
                      /*
                      dup->product = SeqLocFree (dup->product);
                      dup->product = CreateWholeInterval (newsep);
                      */
                    }
                    SeqMgrReplaceInBioseqIndex (newbsp);
                  }
                }
                AddSeqEntryToSeqEntry (top, newsep, TRUE);
              }
            }
          }
        }
      }
    }

    sfp = SeqMgrGetNextFeature (udp->newbsp, sfp, 0, 0, &context);
  }

  ReplaceBioSourceAndPubs (top, sdp);

  if (sapp != NULL) {
    *sapp = sap;
  }

  return TRUE;
}

static void UpdateOneFeatureForSequenceReplace 
(SeqFeatPtr  sfp, 
 SeqAlignPtr salp,
 BioseqPtr   oldbsp,
 Int4        new_len)
{
  CodeBreakPtr cbp;
  CdRegionPtr  crp;
  RnaRefPtr    rrp;
  tRNAPtr      trp;
  
  if (sfp == NULL || salp == NULL)
  {
    return;
  }
  
  ReplaceComplexLocation (sfp->location, salp, new_len, 1, 2);

  switch (sfp->data.choice) {
    case SEQFEAT_CDREGION :
      crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      if (crp != NULL) 
      {
        for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) 
        {
          ReplaceComplexLocation (cbp->loc, salp, new_len, 1, 2);
        }
      }
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp != NULL && rrp->ext.choice == 2) {
        trp = (tRNAPtr) rrp->ext.value.ptrvalue;
        if (trp != NULL && trp->anticodon != NULL) {
          ReplaceComplexLocation (trp->anticodon, salp, new_len, 1, 2);
        }
      }
      break;
    default :
      break;
  }
}

static void UpdateLocationsForSequenceReplace 
(SeqAlignPtr salp,
 BioseqPtr oldbsp,
 BioseqPtr newbsp)
{
  BioseqPtr    parentbsp;
  SeqMgrFeatContext  context;
  SeqMgrSegmentContext segcontext;
  SeqFeatPtr           sfp;
  
  if (salp == NULL || oldbsp == NULL || newbsp == NULL)
  {
    return;
  }
    
  /* if this sequence is a part, the features will be indexed on
   * the parent.
   */
  parentbsp = SeqMgrGetParentOfPart (oldbsp, &segcontext);
  if (parentbsp == NULL)
  {
    sfp = SeqMgrGetNextFeature (oldbsp, NULL, 0, 0, &context);
    while (sfp != NULL)
    {
      UpdateOneFeatureForSequenceReplace (sfp, salp, oldbsp, newbsp->length);
      sfp = SeqMgrGetNextFeature (oldbsp, sfp, 0, 0, &context);
                                  
    }
  }
  else
  {
    sfp = SeqMgrGetNextFeature (parentbsp, NULL, 0, 0, &context);
    while (sfp != NULL)
    {
      UpdateOneFeatureForSequenceReplace (sfp, salp, oldbsp, newbsp->length);
      sfp = SeqMgrGetNextFeature (parentbsp, sfp, 0, 0, &context);
    }
  }  
}

static void 
ReplaceOneSequence 
(SeqAlignPtr salp,
 BioseqPtr oldbsp,
 BioseqPtr newbsp)
{
  ByteStorePtr       bs;
  Int4               len, len_change;
  Uint1              seq_data_type, seq_ext_type;
  Pointer            seq_ext;
  Uint1              repr;
  BioseqPtr          parent_bsp;
  SeqMgrSegmentContext context;
  
  if (oldbsp == NULL || newbsp == NULL)
  {
    return;
  }

  UpdateLocationsForSequenceReplace (salp, oldbsp, newbsp);
  len_change = newbsp->length - oldbsp->length;

  /* switch bioseqs to finish update */

  bs = oldbsp->seq_data;
  oldbsp->seq_data = newbsp->seq_data;
  newbsp->seq_data = bs;
  len = oldbsp->length;
  oldbsp->length = newbsp->length;
  newbsp->length = len;
  seq_data_type = oldbsp->seq_data_type;
  oldbsp->seq_data_type = newbsp->seq_data_type;
  newbsp->seq_data_type = seq_data_type;  
  /* also move seq_ext, for delta sequences */
  seq_ext_type = oldbsp->seq_ext_type;
  seq_ext = oldbsp->seq_ext;
  oldbsp->seq_ext_type = newbsp->seq_ext_type;
  oldbsp->seq_ext = newbsp->seq_ext;
  newbsp->seq_ext_type = seq_ext_type;
  newbsp->seq_ext = seq_ext;
  
  /* swap repr */
  repr = oldbsp->repr;
  oldbsp->repr = newbsp->repr;
  newbsp->repr = repr;
  
  /* if this was part of a segmented set, update the parent length */
  parent_bsp = SeqMgrGetParentOfPart (oldbsp, &context);
  if (parent_bsp != NULL)
  {
    parent_bsp->length += len_change;
  }
}

static Boolean ReplaceSequence (UpsDataPtr udp)

{
  MsgAnswer          ans;

  if (udp == NULL)
  {
    return TRUE;
  }
  
  if (FALSE == udp->isSet)
  {
    if ((udp->seq1 != NULL || udp->seq2 != NULL)
        && StringICmp (udp->seq1, udp->seq2) == 0
        && ! udp->revcomp) 
    {
	    ans = Message (MSG_OKC, "Replacement sequence is identical to"
	                           " original - possible error");
	    if (ans == ANS_CANCEL) return FALSE;
    }    
  }

  ReplaceOneSequence (udp->salp, udp->oldbsp, udp->newbsp);
  return TRUE;
}

static Boolean Merge5Prime (UpsDataPtr udp)

{
  ByteStorePtr  bs;
  Char          ch;
  Int4          i, newlen;
  BioseqPtr     newbsp;
  CharPtr       ptr, str, tmp;

  /* construct replacement sequence by recombining between old and overlap */

  tmp = udp->seq2;

  newlen = udp->new5 + udp->newa + udp->old3;
  str = (CharPtr) MemNew (sizeof (Char) * (size_t) (newlen + 5));
  if (str == NULL)
    return FALSE;
  ptr = str;

  for (i = 0; i < udp->new5 + udp->newa; i++) {
    ch = *tmp;
    *ptr = ch;
    tmp++;
    ptr++;
  }
    
  tmp = udp->seq1 + udp->old5 + udp->olda;
  for (i = 0; i < udp->old3; i++) {
    ch = *tmp;
    *ptr = ch;
    tmp++;
    ptr++;
  }

  *ptr = '\0';
  bs = BSNew (newlen);
  BSWrite (bs, (VoidPtr) str, newlen);

  udp->seq2 = MemFree (udp->seq2);
  udp->seq2 = str;

  if (bs != NULL && BSLen (bs) < 1) {
    bs = BSFree (bs);
  }
  if (bs == NULL) return FALSE;

  /* overlap turned into replacement sequence */

  newbsp = udp->newbsp;
  newbsp->seq_data = BSFree (newbsp->seq_data);
  newbsp->seq_data = bs;
  newbsp->seq_data_type = Seq_code_iupacna;
  newbsp->length = newlen;

  /* adjust alignment and reindex */

  if (! AdjustAlignment (udp, 2)) return FALSE;

  /* then finish by replacing with new sequence */

  return ReplaceSequence (udp);
}

static Boolean Merge3Prime (UpsDataPtr udp)

{
  ByteStorePtr  bs;
  Char          ch;
  Int4          i, newlen;
  BioseqPtr     newbsp;
  CharPtr       ptr, str, tmp;

  /* construct replacement sequence by recombining between old and overlap */

  newlen = udp->old5 + udp->newa + udp->new3;
  str = (CharPtr) MemNew (sizeof (Char) * (size_t) (newlen + 5));
  if (str == NULL)
    return FALSE;
  ptr = str;

  tmp = udp->seq1;
  for (i = 0; i < udp->old5; i++) {
    ch = *tmp;
    *ptr = ch;
    tmp++;
    ptr++;
  }
    
  tmp = udp->seq2 + udp->new5;
  for (i = 0; i < udp->newa + udp->new3; i++) {
    ch = *tmp;
    *ptr = ch;
    tmp++;
    ptr++;
  }

  *ptr = '\0';
  bs = BSNew (newlen);
  BSWrite (bs, (VoidPtr) str, newlen);

  udp->seq2 = MemFree (udp->seq2);
  udp->seq2 = str;

  if (bs != NULL && BSLen (bs) < 1) {
    bs = BSFree (bs);
  }
  if (bs == NULL) return FALSE;

  /* overlap turned into replacement sequence */

  newbsp = udp->newbsp;
  newbsp->seq_data = BSFree (newbsp->seq_data);
  newbsp->seq_data = bs;
  newbsp->seq_data_type = Seq_code_iupacna;
  newbsp->length = newlen;

  /* adjust alignment and reindex */

  if (! AdjustAlignment (udp, 3)) return FALSE;

  /* then finish by replacing with new sequence */

  return ReplaceSequence (udp);
}

static Boolean ExtendFeatures (UpsDataPtr udp, Int4 offset);

/*------------------------------------------------------------------*/
/*                                                                  */
/*  Merge5PrimeNoOverlap () -- Merge a new sequence onto the 5' end */
/*                             of an existing sequence.             */
/*                                                                  */
/*                             Performs a similar function to       */
/*                             Merge5Prime() except works when      */
/*                             there is no alignment between the    */
/*                             two sequences.                       */
/*                                                                  */
/*------------------------------------------------------------------*/

static Boolean Merge5PrimeNoOverlap (UpsDataPtr udp)

{
  CharPtr       origSeqStr;
  CharPtr       newSeqStr;
  CharPtr       mergedSeqStr;
  Int4          mergedLen;
  ByteStorePtr  mergedBS;

  /* Get original and new sequences */

  origSeqStr = GetSequenceByBsp (udp->oldbsp);
  newSeqStr = GetSequenceByBsp (udp->newbsp);

  /* Concatenate the new sequence onto the beginning */
  /* (i.e. the 5' end) of the original sequence.     */

  mergedLen =  StringLen (newSeqStr) + StringLen (origSeqStr);
  mergedSeqStr = (CharPtr) MemNew (mergedLen + 1);
  sprintf (mergedSeqStr, "%s%s", newSeqStr, origSeqStr);

  /* Convert the new sequence into a ByteStore */

  mergedBS = BSNew (mergedLen);
  BSWrite (mergedBS, (VoidPtr) mergedSeqStr, mergedLen);

  /* Replace the original sequence with the */
  /* new concatenated sequence.             */

  udp->newbsp->seq_data      = BSFree (udp->newbsp->seq_data);
  udp->newbsp->seq_data      = mergedBS;
  udp->newbsp->seq_data_type = Seq_code_iupacna;
  udp->newbsp->length        = mergedLen;

  /* Replace the merged sequence and return */

  return ExtendFeatures (udp, StringLen (newSeqStr));
}

/*------------------------------------------------------------------*/
/*                                                                  */
/*  Merge3PrimeNoOverlap () -- Merge a new sequence onto the 3' end */
/*                             of an existing sequence.             */
/*                                                                  */
/*                             Performs a similar function to       */
/*                             Merge3Prime() except works when      */
/*                             there is no alignment between the    */
/*                             two sequences.                       */
/*                                                                  */
/*------------------------------------------------------------------*/

static Boolean Merge3PrimeNoOverlap (UpsDataPtr udp)

{
  CharPtr       origSeqStr;
  CharPtr       newSeqStr;
  CharPtr       mergedSeqStr;
  Int4          mergedLen;
  ByteStorePtr  mergedBS;

  /* Get original and new sequences */

  origSeqStr = GetSequenceByBsp (udp->oldbsp);
  newSeqStr = GetSequenceByBsp (udp->newbsp);

  /* Concatenate the new sequence onto the end   */
  /* (i.e. the 3' end) of the original sequence. */

  mergedLen =  StringLen (newSeqStr) + StringLen (origSeqStr);
  mergedSeqStr = (CharPtr) MemNew (mergedLen + 1);
  sprintf (mergedSeqStr, "%s%s", origSeqStr, newSeqStr);

  /* Convert the new sequence into a ByteStore */

  mergedBS = BSNew (mergedLen);
  BSWrite (mergedBS, (VoidPtr) mergedSeqStr, mergedLen);

  /* Replace the original sequence with the */
  /* new concatenated sequence.             */

  udp->newbsp->seq_data      = BSFree (udp->newbsp->seq_data);
  udp->newbsp->seq_data      = mergedBS;
  udp->newbsp->seq_data_type = Seq_code_iupacna;
  udp->newbsp->length        = mergedLen;

  /* Replace the merged sequence and return */

  return ExtendFeatures (udp, 0);
}

static Boolean OkToPatchDelta (UpsDataPtr udp)
{
  Boolean rval = TRUE;
  
  if (udp == NULL || udp->oldbsp == NULL || udp->newbsp == NULL
      || udp->oldbsp->repr != Seq_repr_delta || udp->newbsp->repr != Seq_repr_delta
      || udp->oldbsp->seq_ext_type != 4 || udp->newbsp->seq_ext_type != 4)
  {
    rval = FALSE;
  }

  return rval;
}

static void SplitDeltaSeq (DeltaSeqPtr dsp, Int4 offset)
{
  SeqLocPtr   slp1, slp2;
  SeqLitPtr   slip1, slip2;
  Int4        len;
  Boolean     changed;
  DeltaSeqPtr dsp_new;
  ByteStorePtr bs_1, bs_2;
  Int2         residue;
  Int4         pos;
  
  if (dsp == NULL || dsp->data.ptrvalue == NULL || offset == 0)
  {
    return;
  }
  
  if (dsp->choice == 1)
  {
    slp1 = (SeqLocPtr)(dsp->data.ptrvalue);
    len = SeqLocLen (slp1);
    if (offset > len)
    {
      return;
    }
    slp2 = (SeqLocPtr) AsnIoMemCopy (slp1, (AsnReadFunc) SeqLocAsnRead,
                                     (AsnWriteFunc) SeqLocAsnWrite);
    slp1 = SeqLocDelete (slp1, SeqLocId (slp1), 
                          offset, len - 1, FALSE, &changed);
    slp2 = SeqLocDelete (slp2, SeqLocId (slp2),
                         0, offset, FALSE, &changed);
    dsp_new = ValNodeNew (NULL);
    dsp_new->choice = 1;
    dsp_new->data.ptrvalue = slp2;
    dsp_new->next = dsp->next;
    dsp->next = dsp_new;
  }
  else if (dsp->choice == 2)
  {
    slip1 = (SeqLitPtr) dsp->data.ptrvalue;
    if (offset > slip1->length)
    {
      return;
    }
    slip2 = SeqLitNew ();
    if (slip1->seq_data == NULL)
    {
      slip2->length = slip1->length - offset;
      slip1->length = offset;
    }
    else
    {
      if (slip1->seq_data_type == Seq_code_iupacna)
      {
        bs_1 = slip1->seq_data;
      }
      else
      {
        bs_1 = BSConvertSeq(slip1->seq_data, Seq_code_iupacna, 
                              slip1->seq_data_type, 
                              slip1->length);
        slip1->seq_data_type = Seq_code_iupacna;
        slip1->seq_data = bs_1;
      }
      bs_2 = BSNew (slip1->length - offset);
      pos = offset;
      BSSeek(bs_1, pos, SEEK_SET);
      BSSeek (bs_2, 0L, SEEK_SET);
      while (pos < slip1->length)
      {
        residue = BSGetByte (bs_1);
        BSPutByte (bs_2, residue);
        pos++;
      }
      BSSeek(bs_1, offset, SEEK_SET);
      BSDelete(bs_1, slip1->length - offset);
      
      slip2->seq_data = bs_2;
      slip2->seq_data_type = slip1->seq_data_type;
      slip2->length = slip1->length - offset;
      slip1->length = offset;
    }
    dsp_new = ValNodeNew (NULL);
    dsp_new->choice = 2;
    dsp_new->data.ptrvalue = slip2;
    dsp_new->next = dsp->next;
    dsp->next = dsp_new;
  }
}

/* This function will patch a delta sequence with another delta sequence.
 * The pieces in the overlap from the old sequence will be replaced by pieces
 * in the overlap from the new sequence.
 */
static Boolean PatchDeltaSequence (UpsDataPtr udp)

{
  Int4        currnew_pos = 0, currold_pos;
  SeqLitPtr   slip, slip_new;
  DeltaSeqPtr dspold, dspnew;
  Int4        seqstart;
  DeltaSeqPtr new_list = NULL;
  
  if (! OkToPatchDelta (udp))
  {
    return FALSE;
  }

  /* keep old 5' end intact */
  currold_pos = 0;
  seqstart = 0;
  dspold = (DeltaSeqPtr) udp->oldbsp->seq_ext;
  while (dspold != NULL && currold_pos < udp->old5)
  {
    seqstart = currold_pos;
    if (dspold->data.ptrvalue == NULL || dspold->choice != 2)
    {
      return FALSE;
    }
    slip = (SeqLitPtr) (dspold->data.ptrvalue);
	  currold_pos += slip->length;
		if (currold_pos > udp->old5)
		{
      SplitDeltaSeq (dspold, udp->old5 - seqstart);
      slip = (SeqLitPtr) (dspold->data.ptrvalue);
      currold_pos = udp->old5;		  
		}
		slip_new = (SeqLitPtr) AsnIoMemCopy (slip, (AsnReadFunc) SeqLitAsnRead,
		                                     (AsnWriteFunc) SeqLitAsnWrite);
		ValNodeAddPointer (&new_list, 2, slip_new);
	  dspold = dspold->next;
  }
  
  /* skip over new 5' end */
  currnew_pos = 0;
  seqstart = 0;
  dspnew = (DeltaSeqPtr) udp->newbsp->seq_ext;
  while (dspnew != NULL && currnew_pos < udp->new5)
  {
    seqstart = currold_pos;
    if (dspnew->data.ptrvalue == NULL || dspnew->choice != 2)
    {
      return FALSE;
    }
	  slip = (SeqLitPtr) (dspnew->data.ptrvalue);
	  currnew_pos += slip->length;
	  if (currnew_pos > udp->new5)
	  {
      SplitDeltaSeq (dspnew, udp->new5 - seqstart);
      currnew_pos = udp->new5;
	  }
	  dspnew = dspnew->next;
  }
  
  /* copy in new overlap */
  while (dspnew != NULL && currnew_pos < udp->new5 + udp->newa)
  {
    seqstart = currold_pos;
    if (dspnew->data.ptrvalue == NULL || dspnew->choice != 2)
    {
      return FALSE;
    }
	  slip = (SeqLitPtr) (dspnew->data.ptrvalue);
	  currnew_pos += slip->length;
		if (currnew_pos > udp->new5 + udp->newa)
		{
      SplitDeltaSeq (dspnew, udp->new5 + udp->newa - seqstart);
      slip = (SeqLitPtr) (dspnew->data.ptrvalue);
      currnew_pos = udp->new5 + udp->newa;		  
		}
		slip_new = (SeqLitPtr) AsnIoMemCopy (slip, (AsnReadFunc) SeqLitAsnRead,
		                                     (AsnWriteFunc) SeqLitAsnWrite);
		ValNodeAddPointer (&new_list, 2, slip_new);
		dspnew = dspnew->next;
  }
  
  /* skip over old overlap */
  
  while (dspold != NULL && currold_pos < udp->old5 + udp->olda)
  {
    seqstart = currold_pos;
    if (dspold->data.ptrvalue == NULL || dspold->choice != 2)
    {
      return FALSE;
    }
    slip = (SeqLitPtr) (dspold->data.ptrvalue);
    currold_pos += slip->length;
    if (currold_pos > udp->old5 + udp->olda)
    {
      SplitDeltaSeq (dspold, udp->new5 + udp->newa - seqstart);
      currold_pos = udp->old5 + udp->olda;		        
    }
    dspold = dspold->next;
  }
  
  /* copy in old 3' */
  
  while (dspold != NULL)
  {
    if (dspold->data.ptrvalue == NULL || dspold->choice != 2)
    {
      return FALSE;
    }
    slip = (SeqLitPtr) (dspold->data.ptrvalue);
		slip_new = (SeqLitPtr) AsnIoMemCopy (slip, (AsnReadFunc) SeqLitAsnRead,
		                                     (AsnWriteFunc) SeqLitAsnWrite);
		ValNodeAddPointer (&new_list, 2, slip_new);
		dspold = dspold->next;
  }
  
  /* free newbsp's old SeqLit List */
  for (dspnew = (DeltaSeqPtr) udp->newbsp->seq_ext;
       dspnew != NULL; 
       dspnew = dspnew->next)
  {
    slip = (SeqLitPtr) (dspnew->data.ptrvalue);
    SeqLitFree (slip);
  }
  udp->newbsp->seq_ext = ValNodeFree (udp->newbsp->seq_ext);
  udp->newbsp->seq_ext = new_list;
  udp->newbsp->length = udp->old5 + udp->newa + udp->old3;
  return TRUE;  
}

static Boolean OkToPatchRaw (UpsDataPtr udp)
{
  Boolean rval = TRUE;
  
  if (udp == NULL || udp->oldbsp == NULL || udp->newbsp == NULL
      || udp->oldbsp->repr != Seq_repr_raw || udp->newbsp->repr != Seq_repr_raw)
  {
    rval = FALSE;
  }

  return rval;
}

static Boolean PatchRawSequence (UpsDataPtr udp)

{
  ByteStorePtr  bs;
  Char          ch;
  Int4          i, newlen;
  BioseqPtr     newbsp;
  CharPtr       ptr, str, tmp;

  newlen = udp->old5 + udp->newa + udp->old3;
  str = (CharPtr) MemNew (sizeof (Char) * (size_t) (newlen + 5));
  if (str == NULL) return FALSE;

  /* construct replacement sequence by double recombination */
  ptr = str;

  tmp = udp->seq1;
  for (i = 0; i < udp->old5; i++) {
    ch = *tmp;
    *ptr = ch;
    tmp++;
    ptr++;
  }

  tmp = udp->seq2 + udp->new5;
  for (i = 0; i < udp->newa; i++) {
    ch = *tmp;
    *ptr = ch;
    tmp++;
    ptr++;
  }

  tmp = udp->seq1 + udp->old5 + udp->olda;
  for (i = 0; i < udp->old3; i++) {
    ch = *tmp;
    *ptr = ch;
    tmp++;
    ptr++;
  }

  *ptr = '\0';
  bs = BSNew (newlen);
  BSWrite (bs, (VoidPtr) str, newlen);

  udp->seq2 = MemFree (udp->seq2);
  udp->seq2 = str;

  if (bs != NULL && BSLen (bs) < 1) {
    bs = BSFree (bs);
  }
  if (bs == NULL) return FALSE;

  /* overlap turned into replacement sequence */

  newbsp = udp->newbsp;
  newbsp->seq_data = BSFree (newbsp->seq_data);
  newbsp->seq_data = bs;
  newbsp->seq_data_type = Seq_code_iupacna;
  newbsp->length = newlen;
  return TRUE;  
}

static Boolean PatchSequence (UpsDataPtr udp)

{
  Boolean rval = FALSE;

  if (OkToPatchRaw (udp))
  {
    rval = PatchRawSequence (udp);
  }
  else if (OkToPatchDelta (udp))
  {
    rval = PatchDeltaSequence (udp);
  }
  
  if (!rval)
  {
    return rval;
  }

  /* adjust alignment and reindex */

  if (! AdjustAlignment (udp, 4)) return FALSE;

  /* then finish by replacing with new sequence */

  return ReplaceSequence (udp);
}

static void MarkProductForDeletion (
  SeqLocPtr product
)

{
  BioseqPtr  bsp;
  SeqIdPtr   sip;

  if (product == NULL) return;
  sip = SeqLocId (product);
  if (sip == NULL) return;
  bsp = BioseqFind (sip);
  if (bsp == NULL) return;
  bsp->idx.deleteme = TRUE;
}

static void CombineTexts (
  CharPtr PNTR txtptr,
  CharPtr PNTR oldtxtptr
)

{
  size_t     len;
  CharPtr    str;

  if (txtptr == NULL || oldtxtptr == NULL) return;

  if (*txtptr == NULL) {

    *txtptr = *oldtxtptr;
    *oldtxtptr = NULL;

  } else if (*oldtxtptr != NULL && StringICmp (*txtptr, *oldtxtptr) != 0) {

    len = StringLen (*txtptr) + StringLen (*oldtxtptr) + 5;
    str = MemNew (sizeof (Char) * len);
    StringCpy (str, *txtptr);
    StringCat (str, "; ");
    StringCat (str, *oldtxtptr);
    *txtptr = MemFree (*txtptr);
    *txtptr = str;
  }
}

static void FuseCommonFeatureFields (
  SeqFeatPtr sfp,
  SeqFeatPtr oldsfp
)

{
  GBQualPtr       lastgbq;
  SeqFeatXrefPtr  lastxref;

  if (sfp == NULL || oldsfp == NULL) return;

  CombineTexts (&(sfp->comment), &(oldsfp->comment));
  CombineTexts (&(sfp->title), &(oldsfp->title));
  CombineTexts (&(sfp->except_text), &(oldsfp->except_text));

  if (sfp->qual == NULL) {
    sfp->qual = oldsfp->qual;
    oldsfp->qual = NULL;
  } else if (oldsfp->qual != NULL) {
    for (lastgbq = sfp->qual; lastgbq->next != NULL; lastgbq = lastgbq->next) continue;
    lastgbq->next = oldsfp->qual;
    oldsfp->qual = NULL;
  }

  ValNodeLink (&(sfp->dbxref), oldsfp->dbxref);
  oldsfp->dbxref = NULL;

  ValNodeLink (&(sfp->cit), oldsfp->cit);
  oldsfp->cit = NULL;

  if (sfp->xref == NULL) {
    sfp->xref = oldsfp->xref;
    oldsfp->xref = NULL;
  } else if (oldsfp->xref != NULL) {
    for (lastxref = sfp->xref; lastxref->next != NULL; lastxref = lastxref->next) continue;
    lastxref->next = oldsfp->xref;
    oldsfp->xref = NULL;
  }

  if (sfp->ext == NULL) {
    sfp->ext = oldsfp->ext;
    oldsfp->ext = NULL;
  } else if (oldsfp->ext != NULL) {
    sfp->ext = CombineUserObjects (sfp->ext, oldsfp->ext);
    oldsfp->ext = NULL;
  }

  sfp->partial |= oldsfp->partial;
  sfp->excpt |= oldsfp->excpt;
  sfp->pseudo |= oldsfp->pseudo;
}

static void FuseFeatures (
  SeqFeatPtr sfp,
  SeqFeatPtr oldsfp
)

{
  GeneRefPtr      grp, oldgrp;
  BioseqPtr       prod, oldprod;
  SeqFeatPtr      prot, oldprot;
  ProtRefPtr      prp, oldprp;
  RnaRefPtr       rrp, oldrrp;
  SeqIdPtr        sip;

  if (sfp == NULL || oldsfp == NULL) return;

  /* merge common fields */

  FuseCommonFeatureFields (sfp, oldsfp);

  /* now deal with type-specific data */

  switch (sfp->data.choice) {
    case SEQFEAT_GENE :
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      oldgrp = (GeneRefPtr) oldsfp->data.value.ptrvalue;
      if (grp == NULL || oldgrp == NULL) return;
      CombineTexts (&(grp->locus), &(oldgrp->locus));
      CombineTexts (&(grp->allele), &(oldgrp->allele));
      CombineTexts (&(grp->desc), &(oldgrp->desc));
      CombineTexts (&(grp->maploc), &(oldgrp->maploc));
      CombineTexts (&(grp->locus_tag), &(oldgrp->locus_tag));
      grp->pseudo |= oldgrp->pseudo;
      ValNodeLink (&(grp->db), oldgrp->db);
      oldgrp->db = NULL;
      ValNodeLink (&(grp->syn), oldgrp->syn);
      oldgrp->syn = NULL;
      break;
    case SEQFEAT_CDREGION :
      sip = SeqLocId (sfp->product);
      prod = BioseqFind (sip);
      sip = SeqLocId (oldsfp->product);
      oldprod = BioseqFind (sip);
      if (prod == NULL || oldprod == NULL) return;
      prot = SeqMgrGetBestProteinFeature (prod, NULL);
      oldprot = SeqMgrGetBestProteinFeature (oldprod, NULL);
      if (prot == NULL || oldprot == NULL) return;
      FuseCommonFeatureFields (prot, oldprot);
      prp = (ProtRefPtr) prot->data.value.ptrvalue;
      oldprp = (ProtRefPtr) oldprot->data.value.ptrvalue;
      if (prp == NULL || oldprp == NULL) return;
      ValNodeLink (&(prp->name), oldprp->name);
      oldprp->name = NULL;
      ValNodeLink (&(prp->ec), oldprp->ec);
      oldprp->ec = NULL;
      ValNodeLink (&(prp->activity), oldprp->activity);
      oldprp->activity = NULL;
      ValNodeLink (&(prp->db), oldprp->db);
      oldprp->db = NULL;
      CombineTexts (&(prp->desc), &(oldprp->desc));
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      oldrrp = (RnaRefPtr) oldsfp->data.value.ptrvalue;
      if (rrp == NULL || oldrrp == NULL) return;
      if (rrp->ext.choice == 1 && oldrrp->ext.choice == 1) {
        CombineTexts ((CharPtr PNTR) &(rrp->ext.value.ptrvalue), (CharPtr PNTR) &(oldrrp->ext.value.ptrvalue));
      }
      break;
    case SEQFEAT_REGION :
    case SEQFEAT_COMMENT :
      if (sfp->data.value.ptrvalue == NULL || oldsfp->data.value.ptrvalue == NULL) return;
      CombineTexts ((CharPtr PNTR) &(sfp->data.value.ptrvalue), (CharPtr PNTR) &(oldsfp->data.value.ptrvalue));
      break;
    default :
      break;
  }
}

static void RemoveOldFeatsInRegion (
  UpsDataPtr udp,
  BioseqPtr bsp,
  SeqAnnotPtr sap
)

{
  SeqMgrFeatContext  context;
  Int4               left, right;
  SeqFeatPtr         sfp;

  if (udp == NULL || bsp == NULL || sap == NULL) return;
  if (sap->type != 1) return;

  left = INT4_MAX;
  right = INT4_MIN;

  for (sfp = (SeqFeatPtr) sap->data; sfp != NULL; sfp = sfp->next) {
    if (sfp != SeqMgrGetDesiredFeature (0, bsp, 0, 0, sfp, &context)) continue;
    left = MIN (left, context.left);
    right = MAX (right, context.right);
  }

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &context);

  while (sfp != NULL) {

    if (context.sap != sap && context.right >= left && context.left <= right) {
      sfp->idx.deleteme = TRUE;
      MarkProductForDeletion (sfp->product);
    }

    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &context);
  }
}

static void RemoveOldFeats (BioseqPtr bsp)

{
  SeqMgrFeatContext  context;
  SeqFeatPtr         sfp;

  if (bsp == NULL) return;

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &context);

  while (sfp != NULL) 
  {
    sfp->idx.deleteme = TRUE;
    MarkProductForDeletion (sfp->product);
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &context);
  }
}


static void ResolveDuplicateFeats (
  UpsDataPtr udp,
  BioseqPtr bsp,
  SeqAnnotPtr sap
)

{
  SeqMgrFeatContext  context, lastcontext;
  Int2               i, j;
  Boolean            ivalssame;
  SeqFeatPtr         lastsfp = NULL, sfp;
  Int2               nobmval;

  if (udp == NULL || bsp == NULL || sap == NULL) return;

  nobmval = GetValue (udp->nobm);
  if (nobmval == UPDATE_FEAT_DUP_USE_BOTH) return; /* keep both */

  SeqMgrIndexFeatures (0, (Pointer) bsp);

  if (nobmval == UPDATE_FEAT_DUP_REPLACE) {
    RemoveOldFeatsInRegion (udp, bsp, sap);
    return;
  }

  lastsfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &context);
  if (lastsfp == NULL) return;

  MemCopy ((Pointer) &lastcontext, (Pointer) &context, sizeof (SeqMgrFeatContext));

  sfp = SeqMgrGetNextFeature (bsp, lastsfp, 0, 0, &context);
  if (sfp == NULL) return;

  while (sfp != NULL) {

    if (context.left == lastcontext.left &&
        context.right == lastcontext.right &&
        context.featdeftype == lastcontext.featdeftype) {

      if (context.strand == lastcontext.strand ||
          lastcontext.strand == Seq_strand_unknown ||
          context.strand == Seq_strand_unknown) {

        ivalssame = TRUE;
        if (context.numivals != lastcontext.numivals ||
            context.ivals == NULL ||
            lastcontext.ivals == NULL) {

          ivalssame = FALSE;

        } else {

          for (i = 0, j = 0; i < lastcontext.numivals; i++, j += 2) {
            if (context.ivals [j] != lastcontext.ivals [j]) {
              ivalssame = FALSE;
            }
            if (context.ivals [j + 1] != lastcontext.ivals [j + 1]) {
              ivalssame = FALSE;
            }
          }
        }

        if (ivalssame &&
            context.sap != lastcontext.sap &&
            (context.sap == sap || lastcontext.sap == sap)) {

          if (nobmval == UPDATE_FEAT_DUP_USE_NEW) { /* keep new */
            if (context.sap == sap) {
              lastsfp->idx.deleteme = TRUE;
              MarkProductForDeletion (lastsfp->product);
            } else if (lastcontext.sap == sap) {
              sfp->idx.deleteme = TRUE;
              MarkProductForDeletion (sfp->product);
            }

          } else if (nobmval == UPDATE_FEAT_DUP_USE_OLD) { /* keep old */
            if (context.sap == sap) {
              sfp->idx.deleteme = TRUE;
              MarkProductForDeletion (sfp->product);
            } else if (lastcontext.sap == sap) {
              lastsfp->idx.deleteme = TRUE;
              MarkProductForDeletion (lastsfp->product);
            }

          } else if (nobmval == UPDATE_FEAT_DUP_MERGE) { /* merge */
            if (context.sap == sap) {
              FuseFeatures (sfp, lastsfp);
              lastsfp->idx.deleteme = TRUE;
              MarkProductForDeletion (lastsfp->product);
            } else if (lastcontext.sap == sap) {
              FuseFeatures (lastsfp, sfp);
              sfp->idx.deleteme = TRUE;
              MarkProductForDeletion (sfp->product);
            }
          }
        }
      }
    }

    lastsfp = sfp;
    MemCopy ((Pointer) &lastcontext, (Pointer) &context, sizeof (SeqMgrFeatContext));

    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &context);
  }
}

extern void AddCitSubToUpdatedSequence (BioseqPtr upd_bsp, Uint2 input_entityID)
{
  SeqEntryPtr top_sep, upd_sep;

  upd_sep = SeqMgrGetSeqEntryForData (upd_bsp);
  if (upd_sep == NULL) return;
  top_sep = GetTopSeqEntryForEntityID ( input_entityID);
  if (top_sep == NULL) return;
  CreateUpdateCitSubFromBestTemplate (top_sep, upd_sep);
}

static Boolean ExtendFeatures (UpsDataPtr udp, Int4 offset)

{
  MsgAnswer          ans;
  ByteStorePtr       bs;
  CodeBreakPtr       cbp;
  SeqMgrFeatContext  context;
  CdRegionPtr        crp;
  Int4               len;
  BioseqPtr          newbsp;
  BioseqPtr          oldbsp;
  RnaRefPtr          rrp;
  Uint1              seq_data_type;
  SeqFeatPtr         sfp;
  SeqIdPtr           sip;
  tRNAPtr            trp;

  if (udp->salp != NULL)
    if (FALSE == udp->isSet)
      if (StringICmp (udp->seq1, udp->seq2) == 0) {
	ans = Message (MSG_OKC, "Replacement sequence is identical to"
		       " original - possible error");
	if (ans == ANS_CANCEL) return FALSE;
      }

  oldbsp = udp->oldbsp;
  newbsp = udp->newbsp;

  sip = SeqIdFindBest (oldbsp->id, 0);
  if (sip == NULL) return FALSE;

  if (offset > 0) {
    sfp = SeqMgrGetNextFeature (oldbsp, NULL, 0, 0, &context);
    while (sfp != NULL) {
      OffsetLocation (sfp->location, offset, sip);
      switch (sfp->data.choice) {
        case SEQFEAT_CDREGION :
          crp = (CdRegionPtr) sfp->data.value.ptrvalue;
          if (crp != NULL) {
            for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
              OffsetLocation (cbp->loc, offset, sip);
            }
          }
          break;
        case SEQFEAT_RNA :
          rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
          if (rrp != NULL && rrp->ext.choice == 2) {
            trp = (tRNAPtr) rrp->ext.value.ptrvalue;
            if (trp != NULL && trp->anticodon != NULL) {
              OffsetLocation (trp->anticodon, offset, sip);
            }
          }
          break;
        default :
          break;
      }
      sfp = SeqMgrGetNextFeature (oldbsp, sfp, 0, 0, &context);
    }
  }

  /* switch bioseqs to finish extension */

  bs = oldbsp->seq_data;
  oldbsp->seq_data = newbsp->seq_data;
  newbsp->seq_data = bs;
  len = oldbsp->length;
  oldbsp->length = newbsp->length;
  newbsp->length = len;
  seq_data_type = oldbsp->seq_data_type;
  oldbsp->seq_data_type = newbsp->seq_data_type;
  newbsp->seq_data_type = seq_data_type;

  return TRUE;
}

static Boolean ExtendBothEnds (UpsDataPtr udp)

{
  ByteStorePtr  bs;
  Char          ch;
  Int4          i, newlen;
  BioseqPtr     newbsp;
  BioseqPtr     oldbsp;
  CharPtr       ptr, str, tmp;

  /* construct replacement sequence by extending old sequence */

  newlen = udp->new5 + udp->olda + udp->new3;
  str = (CharPtr) MemNew (sizeof (Char) * (size_t) (newlen + 5));
  if (str == NULL)
    return FALSE;
  ptr = str;

  tmp = udp->seq2;
  for (i = 0; i < udp->new5; i++) {
    ch = *tmp;
    *ptr = ch;
    tmp++;
    ptr++;
  }
    
  tmp = udp->seq1 + udp->old5;
  for (i = 0; i < udp->olda; i++) {
    ch = *tmp;
    *ptr = ch;
    tmp++;
    ptr++;
  }

  tmp = udp->seq2 + udp->new5 + udp->newa;
  for (i = 0; i < udp->new3; i++) {
    ch = *tmp;
    *ptr = ch;
    tmp++;
    ptr++;
  }

  *ptr = '\0';
  bs = BSNew (newlen);
  BSWrite (bs, (VoidPtr) str, newlen);

  udp->seq2 = MemFree (udp->seq2);
  udp->seq2 = str;

  if (bs != NULL && BSLen (bs) < 1) {
    bs = BSFree (bs);
  }
  if (bs == NULL) return FALSE;

  /* overlap turned into replacement sequence */

  oldbsp = udp->oldbsp;
  newbsp = udp->newbsp;
  newbsp->seq_data = BSFree (newbsp->seq_data);
  newbsp->seq_data = bs;
  newbsp->seq_data_type = Seq_code_iupacna;
  newbsp->length = newlen;

  /* then finish by offsetting features */

  return ExtendFeatures (udp, udp->new5);
}

static Boolean Extend5Prime (UpsDataPtr udp)

{
  ByteStorePtr  bs;
  Char          ch;
  Int4          i, newlen;
  BioseqPtr     newbsp;
  BioseqPtr     oldbsp;
  CharPtr       ptr, str, tmp;

  /* construct replacement sequence by extending old sequence */

  newlen = udp->new5 + udp->olda + udp->old3;
  str = (CharPtr) MemNew (sizeof (Char) * (size_t) (newlen + 5));
  if (str == NULL)
    return FALSE;
  ptr = str;

  tmp = udp->seq2;
  for (i = 0; i < udp->new5; i++) {
    ch = *tmp;
    *ptr = ch;
    tmp++;
    ptr++;
  }
    
  tmp = udp->seq1 + udp->old5;
  for (i = 0; i < udp->olda + udp->old3; i++) {
    ch = *tmp;
    *ptr = ch;
    tmp++;
    ptr++;
  }

  *ptr = '\0';
  bs = BSNew (newlen);
  BSWrite (bs, (VoidPtr) str, newlen);

  udp->seq2 = MemFree (udp->seq2);
  udp->seq2 = str;

  if (bs != NULL && BSLen (bs) < 1) {
    bs = BSFree (bs);
  }
  if (bs == NULL) return FALSE;

  /* overlap turned into replacement sequence */

  oldbsp = udp->oldbsp;
  newbsp = udp->newbsp;
  newbsp->seq_data = BSFree (newbsp->seq_data);
  newbsp->seq_data = bs;
  newbsp->seq_data_type = Seq_code_iupacna;
  newbsp->length = newlen;

  /* then finish by offsetting features */

  return ExtendFeatures (udp, udp->new5);
}

static Boolean Extend3Prime (UpsDataPtr udp)

{
  ByteStorePtr  bs;
  Char          ch;
  Int4          i, newlen;
  BioseqPtr     newbsp;
  BioseqPtr     oldbsp;
  CharPtr       ptr, str, tmp;

  /* construct replacement sequence by extending old sequence */

  newlen = udp->old5 + udp->olda + udp->new3;
  str = (CharPtr) MemNew (sizeof (Char) * (size_t) (newlen + 5));
  if (str == NULL)
    return FALSE;
  ptr = str;

  tmp = udp->seq1;
  for (i = 0; i < udp->old5 + udp->olda; i++) {
    ch = *tmp;
    *ptr = ch;
    tmp++;
    ptr++;
  }
    
  tmp = udp->seq2 + udp->new5 + udp->newa;
  for (i = 0; i < udp->new3; i++) {
    ch = *tmp;
    *ptr = ch;
    tmp++;
    ptr++;
  }

  *ptr = '\0';
  bs = BSNew (newlen);
  BSWrite (bs, (VoidPtr) str, newlen);

  udp->seq2 = MemFree (udp->seq2);
  udp->seq2 = str;

  if (bs != NULL && BSLen (bs) < 1) {
    bs = BSFree (bs);
  }
  if (bs == NULL) return FALSE;

  /* overlap turned into replacement sequence */

  oldbsp = udp->oldbsp;
  newbsp = udp->newbsp;
  newbsp->seq_data = BSFree (newbsp->seq_data);
  newbsp->seq_data = bs;
  newbsp->seq_data_type = Seq_code_iupacna;
  newbsp->length = newlen;

  /* no offset, but ExtendFeatures finishes switch */

  return ExtendFeatures (udp, 0);
}

static Boolean ExtendOneSequence (UpsDataPtr udp)
{
  Boolean      update = FALSE;
  MsgAnswer    ans;
  Uint2        entityID;
  CharPtr      errmsg = NULL;
  
  
  if (udp == NULL) return FALSE;
  
  switch (udp->rmcval) {
    case UPDATE_REPLACE :
      if (udp->old5 > 0 && udp->old3 > 0) {
        errmsg = "Unaligned sequence at 5' and 3' ends.  Do you wish to proceed?";
      } else if (udp->old5 > 0) {
        errmsg = "Unaligned sequence at 5' end.  Do you wish to proceed?";
      } else if (udp->old3 > 0) {
        errmsg = "Unaligned sequence at 3' end.  Do you wish to proceed?";
      }
      break;
    case UPDATE_EXTEND5 :
      if (udp->old5 > 0) {
        errmsg = "Unaligned sequence at 5' end.  Do you wish to proceed?";
      }
      break;
    case UPDATE_EXTEND3 :
      if (udp->old3 > 0) {
        errmsg = "Unaligned sequence at 3' end.  Do you wish to proceed?";
      }
      break;
    default :
      break;
  }
  if (errmsg != NULL) {
    ans = Message (MSG_YN, "%s", errmsg);
    if (ans == ANS_NO) {
      return FALSE;
    }
  }

  switch (udp->rmcval) {
    case UPDATE_REPLACE :
      if (ExtendBothEnds (udp)) {
        update = TRUE;
      }
      break;
    case UPDATE_EXTEND5:
      if (udp->salp == NULL) {
        if (Merge5PrimeNoOverlap (udp)) {
          update = TRUE;
        }
      } else if (Extend5Prime (udp)) {
        update = TRUE;
      }
      break;
    case UPDATE_EXTEND3 :
      if (udp->salp == NULL) {
        if (Merge3PrimeNoOverlap (udp)) {
          update = TRUE;
        }
      } else if (Extend3Prime (udp)) {
        update = TRUE;
      }
      break;
    default :
      break;
  }

  if (update) {
      
    entityID = ObjMgrGetEntityIDForPointer (udp->oldbsp);
    if (GetStatus (udp->add_cit_subs))
    {
      AddCitSubToUpdatedSequence ( udp->oldbsp, entityID);
    }
  }
  return update;
}

static void OpenSequenceUpdateLog (UpsDataPtr udp)
{
  if (udp == NULL || udp->log_fp != NULL)
  {
  	return;
  }
  TmpNam (udp->log_path);
  udp->log_fp = FileOpen (udp->log_path, "wb");
}

static void CloseOutSequenceUpdateLog (UpsDataPtr udp)
{
  if (udp == NULL || udp->log_fp == NULL) 
  {
    return;
  }
  FileClose (udp->log_fp);
  udp->log_fp = NULL;
  if (udp->data_in_log) {
    LaunchGeneralTextViewer (udp->log_path, "Protein changes");
    udp->data_in_log = FALSE;
  }
  FileRemove (udp->log_path);  	
}


static Boolean PrepareToUpdateSequences (UpsDataPtr udp);
static Boolean PrepareUpdatePtr (UpsDataPtr    udp);
static ForM UpdateSequenceForm (UpsDataPtr udp);
static void UpdateOneSequence (
  UpsDataPtr   udp,
  Int2         sfbval,
  Boolean      add_cit_subs,
  Boolean      update_proteins);


static Int4 FindNewCDSStop (SeqLocPtr slp, BioseqPtr bsp, Int4 prot_len)
{
  Int4      loc_len, tot_len = 0;
  SeqLocPtr this_slp;
  Int4      curr_start;
  
  for (this_slp = SeqLocFindNext (slp, NULL);
       this_slp != NULL;
       this_slp = SeqLocFindNext (slp, this_slp))
  {
    loc_len = SeqLocLen (this_slp);
    curr_start = GetOffsetInBioseq (this_slp, bsp, SEQLOC_START);
    if (loc_len + tot_len > prot_len)
    {
      curr_start = GetOffsetInBioseq (this_slp, bsp, SEQLOC_START);
      if (SeqLocStrand (this_slp) == Seq_strand_minus)
      {
        return curr_start - (prot_len - tot_len) + 1;
      }
      else
      {
        return curr_start + prot_len - tot_len - 1;
      }
    }
    tot_len += loc_len;
  }
  return tot_len;
}


static CharPtr FixProteinString
(SeqFeatPtr    sfp,
 Uint1         strand,
 Boolean       truncate_proteins,
 Boolean PNTR  truncated,
 Boolean PNTR  contains_start,
 Boolean PNTR  contains_stop)
{
  ByteStorePtr bs;
  BioseqPtr    nucBsp;
  CharPtr      newprot;
  CharPtr      ptr;
  Char         ch;
  Int4         start, stop, new_stop;
  Boolean      changed;
  
  if (sfp == NULL || truncated == NULL 
      || contains_start == NULL
      || contains_stop == NULL) {
    return NULL;
  }
  
  bs = ProteinFromCdRegionEx (sfp, TRUE, FALSE);
  if (bs == NULL) return NULL;
  newprot = BSMerge (bs, NULL);
  bs = BSFree (bs);
  if (newprot == NULL) return NULL;

  ptr = newprot;
  ch = *ptr;
  if (ch == 'M') {
    *contains_start = TRUE;
  } else {
    *contains_start = FALSE;
  }
  *contains_stop = FALSE;
  *truncated = FALSE;
  while (ch != '\0')
  {
    *ptr = TO_UPPER (ch);
    if (ch == '*') {
      *contains_stop = TRUE;
      if (*(ptr + 1) == 0 || truncate_proteins) {
        *ptr = 0;
        if (truncate_proteins && *(ptr + 1) != 0) {
          *truncated = TRUE;
          nucBsp = BioseqFindFromSeqLoc (sfp->location);
          if (nucBsp == NULL) return newprot;
          start = GetOffsetInBioseq (sfp->location, nucBsp, SEQLOC_START);
          stop = GetOffsetInBioseq (sfp->location, nucBsp, SEQLOC_STOP);
          new_stop = FindNewCDSStop (sfp->location, nucBsp, (1 + ptr - newprot) * 3);
          if (strand == Seq_strand_minus) {
            sfp->location = SeqLocDelete (sfp->location, SeqLocId (sfp->location),
                          stop, new_stop - 1, FALSE, &changed);
          } else {
            sfp->location = SeqLocDelete (sfp->location, SeqLocId (sfp->location),
                          new_stop + 1, stop, FALSE, &changed);
          }
        }
        return newprot;
      }
    }
    ptr++;
    ch = *ptr;
  }
  return newprot;
}

/* This function will truncate base pairs to make the length of the CDS
 * a multiple of three.  If truncation occurs, the function returns TRUE,
 * otherwise the function returns FALSE.
 */
static SeqLocPtr ShiftStopForLengthChange
(SeqLocPtr slp,
 BioseqPtr nucBsp,
 Int4      transl_except_len,
 Boolean PNTR changed)
{
  Int4      start, stop, new_stop;
  Uint1     strand;
  SeqIdPtr  sip;
  Int4      len;

  if (slp == NULL || nucBsp == NULL || changed == NULL) {
    return NULL;
  }

  *changed = FALSE;

  start = GetOffsetInBioseq (slp, nucBsp, SEQLOC_START);
  stop = GetOffsetInBioseq (slp, nucBsp, SEQLOC_STOP);
  new_stop = stop;
  strand = SeqLocStrand (slp);
  len = SeqLocLen (slp);

  if (strand == Seq_strand_minus) {
    if (len % 3 != transl_except_len) {
      new_stop += len % 3 + transl_except_len;
      sip = SeqLocId (slp);
      slp = SeqLocDelete (slp, sip, stop, new_stop - 1, FALSE, changed);
      if (slp == NULL) {
        return NULL;
      }
    }
  } else {
    if (len % 3 != transl_except_len) {
      new_stop -= len % 3 - transl_except_len;
      sip = SeqLocId (slp);
      slp = SeqLocDelete (slp, sip, new_stop + 1, stop, FALSE, changed);
      if (slp == NULL) {
        return NULL;
      }
    }
  }
  return slp;
}

extern SeqLocPtr 
ExpandSeqLoc 
(Int4 start,
 Int4 stop,
 Uint1 strand,
 BioseqPtr bsp,
 SeqLocPtr slp)
{
  Int4      curr_start, curr_stop, tmp_start, tmp_stop;
  SeqLocPtr this_slp;

  if (slp == NULL) return NULL;

  if (slp->choice == SEQLOC_INT || slp->choice == SEQLOC_PNT) {
    slp = expand_seq_loc (start, stop, strand, slp);
  } else {
    curr_start = GetOffsetInBioseq (slp, bsp, SEQLOC_START);
    curr_stop = GetOffsetInBioseq (slp, bsp, SEQLOC_STOP);
    for (this_slp = SeqLocFindNext (slp, NULL);
         this_slp != NULL;
         this_slp = SeqLocFindNext (slp, this_slp))
    {
      tmp_start = GetOffsetInBioseq (this_slp, bsp, SEQLOC_START);
      tmp_stop = GetOffsetInBioseq (this_slp, bsp, SEQLOC_STOP);
      if (strand == Seq_strand_minus) {
        if (tmp_start == curr_start && tmp_start < stop) {
          this_slp = expand_seq_loc (tmp_stop, stop, strand, this_slp);
          tmp_start = stop;
        }
        if (tmp_stop == curr_stop && tmp_stop > start) {
          this_slp = expand_seq_loc (start, tmp_start, strand, this_slp);
        }
      } else {
        if (tmp_start == curr_start && tmp_start > start) {
          this_slp = expand_seq_loc (start, tmp_stop, strand, this_slp);
          tmp_start = start;
        }
        if (tmp_stop == curr_stop && tmp_stop < stop) {
          this_slp = expand_seq_loc (tmp_start, stop, strand, this_slp);
        }
      }
    }
  }
  return slp;
}

static SeqLocPtr ShiftLocationForFrame
(SeqLocPtr slp,
 Uint1 frame,
 BioseqPtr nucBsp)
{
  Int4      max_stop, start, stop;
  Uint1     strand;
  SeqIdPtr  sip;
  Int4      offset, new_start;
  Boolean   changed;
  

  if (slp == NULL || nucBsp == NULL) return NULL;
  if (frame == 0 || frame == 1) return slp;

  start = GetOffsetInBioseq (slp, nucBsp, SEQLOC_START);
  stop = GetOffsetInBioseq (slp, nucBsp, SEQLOC_STOP);
  strand = SeqLocStrand (slp);
  max_stop = nucBsp->length - 1;
  new_start = start;
 
  offset = frame - 1;
  sip = SeqLocId (slp);
  if (strand == Seq_strand_minus) {
    stop -= offset;
    new_start = start - offset;
    if ((1 + new_start - stop ) % 3 != 0) {
      new_start -= (1 + new_start - stop) % 3;
    }
    if (stop < 0) stop = 0;
    slp = ExpandSeqLoc (stop, start, strand, nucBsp, slp);
    if (new_start < 0) new_start = 0;
    if (new_start < start) {
      slp = SeqLocDelete (slp, sip, new_start + 1, start, FALSE, &changed);
    }
  } else {
    stop += offset;
    new_start = start + offset;
    if ((1 + stop - new_start) % 3 != 0) {
      new_start += (1 + stop - new_start) % 3;
    }
    if (stop > max_stop) stop = max_stop;
    slp = ExpandSeqLoc (start, stop, strand, nucBsp, slp);
    if (start > max_stop) start = max_stop;
    if (new_start > start) {
      slp = SeqLocDelete (slp, sip, start, new_start - 1, FALSE, &changed);
    }
  }
  return slp;
}

static TransTablePtr GetTranslationTable (CdRegionPtr crp, Boolean PNTR table_is_local)
{
  TransTablePtr tbl = NULL;
  ValNodePtr    vnp;
  Int2          genCode = 0;
  Char          str [32];

  if (crp == NULL || table_is_local == NULL) return NULL;

  *table_is_local = FALSE;
  /* find genetic code */

  if (crp->genetic_code != NULL) {
    vnp = (ValNodePtr) crp->genetic_code->data.ptrvalue;
    while (vnp != NULL) {
      if (vnp->choice == 2) {
        genCode = (Int2) vnp->data.intvalue;
      }
      vnp = vnp->next;
    }
  }

  if (genCode == 7) {
    genCode = 4;
  } else if (genCode == 8) {
    genCode = 1;
  } else if (genCode == 0) {
    genCode = 1;
  }

  /* set up translation table */
  /* set app property name for storing desired FSA */

  sprintf (str, "TransTableFSAforGenCode%d", (int) genCode);

  /* get FSA for desired genetic code if it already exists */

  tbl = (TransTablePtr) GetAppProperty (str);
  if (tbl == NULL) {
    tbl = TransTableNew (genCode);
    *table_is_local = TRUE;
  }
  return tbl;
}

static CharPtr ExtendProtein5
(SeqFeatPtr sfp,
 Uint2      input_entityID,
 Boolean    force_partial)
{
  CdRegionPtr   crp;
  TransTablePtr tbl = NULL;
  Boolean       table_is_local = FALSE;
  SeqLocPtr     test_slp;
  SeqIdPtr      sip;
  BioseqPtr     nucBsp;
  Int4          start, new_start;
  Int4          stop, new_stop;
  Int4          increment = 3000;
  Int4          offset;
  Uint1         strand;
  Boolean       found_start = FALSE;
  Boolean       found_stop = FALSE;
  Boolean       stop_looking = FALSE;
  Int4          dna_len;
  CharPtr       bases;
  Int4          total;
  Int2          state;
  CharPtr       codon_start;
  Boolean       partial3, partial5;
  Boolean       contains_start, contains_stop;
  Boolean       changed;
  CharPtr       newprot;
  Boolean       truncated;
 
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION 
      || (crp = (CdRegionPtr)sfp->data.value.ptrvalue) == NULL) {
    return NULL;
  }
  nucBsp = GetBioseqGivenSeqLoc (sfp->location, input_entityID);
  if (nucBsp == NULL) return NULL;

  tbl = GetTranslationTable (crp, &table_is_local);
  if (tbl == NULL) return NULL;

  test_slp = SeqLocMerge (nucBsp, sfp->location, NULL, FALSE, FALSE, FALSE);
  strand = SeqLocStrand (sfp->location);
  sip = SeqLocId (sfp->location);
  offset = -1;

  start = GetOffsetInBioseq (test_slp, nucBsp, SEQLOC_START);
  if (start == 0)
  {
  	stop_looking = TRUE;
  }
 
  while (((! found_start && ! found_stop) || force_partial) && ! stop_looking) {
    start = GetOffsetInBioseq (test_slp, nucBsp, SEQLOC_START);
    stop = GetOffsetInBioseq (test_slp, nucBsp, SEQLOC_STOP);
    if (strand == Seq_strand_minus) {
      new_start = start + increment;
      new_stop = start + 1;
      if (new_start > nucBsp->length - 1) {
        new_start = start + ((Int4)((nucBsp->length - 1 - start) / 3)) * 3;
        stop_looking = TRUE;
      }
      test_slp = ExpandSeqLoc (stop, new_start, strand, nucBsp, test_slp);
      test_slp = SeqLocDelete (test_slp, sip, stop, start, FALSE, &changed);
    } else {
      new_start = start - increment;
      new_stop = start - 1;
      if (new_start < 0) {
        new_start = start % 3;
        stop_looking = TRUE;
      }
      test_slp = ExpandSeqLoc (new_start, stop, strand, nucBsp, test_slp);
      test_slp = SeqLocDelete (test_slp, sip, start, stop, FALSE, &changed);
    }
    dna_len = SeqLocLen (test_slp);
    bases = ReadCodingRegionBases (test_slp, dna_len, crp->frame, &total);
    if (bases == NULL) {
      stop_looking = TRUE;
    } else {
      state = 0;
      codon_start = bases + StringLen (bases) - 3;
      while (codon_start >= bases && ! found_start && ! found_stop) {
        state = 0;
        state = NextCodonState (tbl, state, (Uint1)*codon_start);
        state = NextCodonState (tbl, state, (Uint1)*(codon_start + 1));
        state = NextCodonState (tbl, state, (Uint1)*(codon_start + 2));
        if (IsOrfStart (tbl, state, TTBL_TOP_STRAND)) {
          found_start = TRUE;
          if (strand == Seq_strand_minus) {
            offset = new_start - (codon_start - bases);
          } else {
            offset = new_start + codon_start - bases;
          }
        } else if (IsOrfStop (tbl, state, TTBL_TOP_STRAND)) {
          found_stop = TRUE;
        } else {
          codon_start -= 3;
        }
      }
      MemFree (bases);
    }
  }
  
  SeqLocFree (test_slp);
  if (! found_stop) { 
    start = GetOffsetInBioseq (sfp->location, nucBsp, SEQLOC_START);
    stop = GetOffsetInBioseq (sfp->location, nucBsp, SEQLOC_STOP);
    if (found_start) {
      if (strand == Seq_strand_minus) {
        sfp->location = ExpandSeqLoc (stop, offset, strand, nucBsp, sfp->location);
      } else {
        sfp->location = ExpandSeqLoc (offset, stop, strand, nucBsp, sfp->location);
      }
    } else {
      CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
      SetSeqLocPartial (sfp->location, TRUE, partial3);
      sfp->partial = TRUE;
      if (crp->frame == 0)
      {
        crp->frame = 1;
      }
      if (strand == Seq_strand_minus) {
        sfp->location = ExpandSeqLoc (stop, nucBsp->length - 1, strand, nucBsp, sfp->location);
        crp->frame = (nucBsp->length - 1 - start + crp->frame - 1) % 3 + 1;
      } else {
        sfp->location = ExpandSeqLoc (0, stop, strand, nucBsp, sfp->location);
        crp->frame = (start + crp->frame - 1) % 3 + 1;
      }
    }
  }
  newprot = FixProteinString (sfp, strand, FALSE, &truncated,
                              &contains_start, &contains_stop);
  if (table_is_local) {
    TransTableFree (tbl);
  }

  return newprot;
}


extern CharPtr ExtendProtein3 
(SeqFeatPtr sfp,
 Uint2      input_entityID,
 Boolean    force_partial)
{
  BioseqPtr nucBsp;
  Int4      max_stop, min_start, start, stop;
  Uint1     strand;
  Boolean   contains_start, contains_stop;
  Int4      increment = 3000;
  CharPtr   newprot;
  Boolean   partial5, partial3;
  Boolean   truncated;

  if (sfp == NULL) return NULL;

  nucBsp = GetBioseqGivenSeqLoc (sfp->location, input_entityID);
  if (nucBsp == NULL) return NULL;

  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  start = GetOffsetInBioseq (sfp->location, nucBsp, SEQLOC_START);
  stop = GetOffsetInBioseq (sfp->location, nucBsp, SEQLOC_STOP);
  strand = SeqLocStrand (sfp->location);
  max_stop = nucBsp->length - 1;
  if (stop > start) {
    while (max_stop % 3 != stop % 3) {
      max_stop --;
    }
    min_start = start %3;
  } else {
    while (max_stop % 3 != start % 3) {
      max_stop --;
    }
    min_start = stop % 3;
  } 
  
  contains_stop = FALSE;
  contains_start = FALSE;
  newprot = NULL;
  /* need to initialize newprot in case we're already at the edge */
  if ((strand != Seq_strand_minus && stop == max_stop)
      || (strand == Seq_strand_minus && stop == min_start))
  {
  	newprot = FixProteinString (sfp, strand, FALSE, &truncated,
                              &contains_start, &contains_stop);
  }
  while ((! contains_stop || force_partial) &&
         (   (strand == Seq_strand_minus && stop > min_start) 
          || (strand != Seq_strand_minus && stop < max_stop)))
  {
    if (newprot != NULL) {
      MemFree (newprot);
      newprot = NULL;
    }
    if (strand == Seq_strand_minus) {
      stop -= increment;
      if (stop < min_start) {
        stop = min_start;
      }
      sfp->location = ExpandSeqLoc (stop, start, strand, nucBsp, sfp->location);
    } else {
      stop += increment;
      if (stop > max_stop) {
        stop = max_stop;
      }
      sfp->location = ExpandSeqLoc (start, stop, strand, nucBsp, sfp->location);
    }
    newprot = FixProteinString (sfp, strand, TRUE, &truncated,
                              &contains_start, &contains_stop);
  }

  if (! contains_stop || force_partial) {
    start = GetOffsetInBioseq (sfp->location, nucBsp, SEQLOC_START);
    stop = GetOffsetInBioseq (sfp->location, nucBsp, SEQLOC_STOP);
    if (strand == Seq_strand_minus) {
      sfp->location = ExpandSeqLoc (0, start, strand, nucBsp, sfp->location);
    } else {
      sfp->location = ExpandSeqLoc (start, nucBsp->length - 1, strand, nucBsp, sfp->location);
    }
    partial3 = TRUE;
    SetSeqLocPartial (sfp->location, partial5, TRUE);
  }
  return newprot;
}

static Boolean MergeOverlapsForThisFeature (SeqFeatPtr sfp)
{  
  if (sfp == NULL) return FALSE;
  return ! sfp->excpt;
}

static Boolean 
AdjustCDSForUpdate
(SeqFeatPtr   sfp,
 Uint2        input_entityID,
 Uint1        frame,
 Int4         transl_except_len,
 Boolean      truncate_proteins,
 Boolean      extend_proteins5,
 Boolean      extend_proteins3,
 Boolean PNTR truncated,
 Boolean PNTR stop_shifted,
 Boolean PNTR extended5,
 Boolean PNTR extended3,
 BioseqPtr    protBsp,
 BioseqPtr PNTR newbsp)
{
  ByteStorePtr  bs;
  CharPtr       newprot, ptr;
  SeqEntryPtr   nwsep;
  BioseqPtr     newBsp;
  Boolean       contains_start, contains_stop;
  Uint1         strand;
  SeqLocPtr     newloc;
  BioseqPtr     nucBsp;
  Boolean       partial5, partial3;
  CharPtr       seqnew = NULL, seqold = NULL;
  Boolean       rval = FALSE;
 
  if (sfp == NULL
    || sfp->idx.subtype != FEATDEF_CDS
    || sfp->product == NULL
    || protBsp == NULL
    || truncated == NULL || stop_shifted == NULL
    || extended5 == NULL || extended3 == NULL
    || newbsp == NULL)
  {
    return FALSE;
  }
  
  CheckSeqLocForPartial (sfp->location, &partial3, &partial5);
  
  nucBsp = GetBioseqGivenSeqLoc (sfp->location, input_entityID);
  if (nucBsp == NULL) return FALSE;
  
  newloc = SeqLocMergeExEx (nucBsp, sfp->location, NULL, FALSE, FALSE,
                            MergeOverlapsForThisFeature (sfp), FALSE, FALSE);

  if (newloc == NULL) return FALSE;
  sfp->location = newloc;
  strand = SeqLocStrand (sfp->location);
  sfp->location = ShiftStopForLengthChange (sfp->location, nucBsp, transl_except_len, stop_shifted);
  if (sfp->location == NULL) {
    return FALSE;
  }

  sfp->location = ShiftLocationForFrame (sfp->location, frame, nucBsp);

  newprot = FixProteinString (sfp, strand, truncate_proteins, truncated,
                              &contains_start, &contains_stop);

  /* Must do 3' end first, otherwise may truncate at stops introduced by expanding 5' end for partiality */
  if ((! contains_stop && extend_proteins3 && transl_except_len == 0)
      || ((extend_proteins3 || partial3) && !truncate_proteins)) {
    MemFree (newprot);
    newprot = ExtendProtein3 (sfp, input_entityID, partial3 && !truncate_proteins);
    if (newprot == NULL) return FALSE;
    *extended3 = TRUE;
  } else {
    *extended3 = FALSE;
  }
  if (! contains_start && (extend_proteins5 || partial5)) {
    MemFree (newprot);
    newprot = ExtendProtein5 (sfp, input_entityID, partial5);
    if (newprot == NULL) return FALSE;
    *extended5 = TRUE;
  } else {
    *extended5 = FALSE;
  }

  sfp->partial = CheckSeqLocForPartial (sfp->location, &partial5, &partial3);

  bs = BSNew (1000);
  if (bs != NULL)
  {
    ptr = newprot;
    BSWrite (bs, (VoidPtr) ptr, (Int4) StringLen (ptr));
  } 
  MemFree (newprot);
  
  newBsp = BioseqNew ();
  if (newBsp == NULL) {
    return FALSE;
  }

  newBsp->id = SeqIdParse ("lcl|ProtAlign");
  newBsp->repr = Seq_repr_raw;
  newBsp->mol = Seq_mol_aa;
  newBsp->seq_data_type = Seq_code_ncbieaa;
  newBsp->seq_data = bs;
  newBsp->length = BSLen (bs);

  /* create SeqEntry for temporary protein bioseq to live in */
  nwsep = SeqEntryNew ();
  nwsep->choice = 1;
  nwsep->data.ptrvalue = newBsp;
  SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) newBsp, nwsep);

  seqold = GetSequenceByBsp (protBsp);
  seqnew = GetSequenceByBsp (newBsp);

  if (StringCmp (seqold, seqnew) == 0
    || (contains_stop && ! *truncated
      && StringNCmp (seqold,
                     seqnew,
                     StringLen (seqold)) == 0))
  {
    SeqEntryFree (nwsep);
    rval = FALSE;
  }
  else
  {
    *newbsp = newBsp;
    rval = TRUE;
  }

  return rval;   
}

  
extern SeqLocPtr SeqLocWholeNew (BioseqPtr bsp)
{
  ValNodePtr vnp;

  if (bsp == NULL) return NULL;

  vnp = ValNodeNew (NULL);

  if (vnp == NULL) return NULL;

  vnp->choice = SEQLOC_WHOLE;
  vnp->data.ptrvalue = (Pointer) SeqIdDup (SeqIdFindBest (bsp->id, 0));
  return (SeqLocPtr)vnp;
}


static void WarnNoProteinUpdate (BioseqPtr bsp)
{
  Char     acc_str [256];
  CharPtr  warn_format = "Protein %s has not been updated - you must manually re-translate the coding regions and move any protein features";
  CharPtr  warn_msg;

  if (bsp == NULL || bsp->id == NULL) return;
  SeqIdWrite (bsp->id, acc_str, PRINTID_REPORT, sizeof (acc_str));
  warn_msg = (CharPtr) MemNew (StringLen (warn_format) + StringLen (acc_str));
  if (warn_msg == NULL) return;
  sprintf (warn_msg, warn_format, acc_str);
  ErrPostEx (SEV_ERROR, 0, 0, warn_msg);
  MemFree (warn_msg);
}


static void LogFrameChange (FILE *fp, BioseqPtr bsp, Uint1 frame)
{
  Char     acc_str [256];

  if (fp == NULL || bsp == NULL || bsp->id == NULL) return;
  SeqIdWrite (bsp->id, acc_str, PRINTID_REPORT, sizeof (acc_str));
  fprintf (fp, "Changed frame for %s to %d\n", acc_str, frame);
}


static void LogProteinTruncate (FILE *fp, BioseqPtr bsp)
{
  Char     acc_str [256];

  if (fp == NULL || bsp == NULL || bsp->id == NULL) return;
  SeqIdWrite (bsp->id, acc_str, PRINTID_REPORT, sizeof (acc_str));
  fprintf (fp, "Truncated protein %s at stop\n", acc_str);
}


static void LogProteinStopShift (FILE *fp, BioseqPtr bsp)
{
  Char     acc_str [256];

  if (fp == NULL || bsp == NULL || bsp->id == NULL) return;
  SeqIdWrite (bsp->id, acc_str, PRINTID_REPORT, sizeof (acc_str));
  fprintf (fp, "Adjusted length of CDS for protein %s to be multiple of 3\n", acc_str);
}


static void LogProteinExtend5 (FILE *fp, BioseqPtr bsp)
{
  Char     acc_str [256];

  if (fp == NULL || bsp == NULL || bsp->id == NULL) return;
  SeqIdWrite (bsp->id, acc_str, PRINTID_REPORT, sizeof (acc_str));
  fprintf (fp, "Extended protein %s on 5' end\n", acc_str);
}


static void LogProteinExtend3 (FILE *fp, BioseqPtr bsp)
{
  Char     acc_str [256];

  if (fp == NULL || bsp == NULL || bsp->id == NULL) return;
  SeqIdWrite (bsp->id, acc_str, PRINTID_REPORT, sizeof (acc_str));
  fprintf (fp, "Extended protein %s on 3' end\n", acc_str);
}

static void LogGeneCorrection (FILE *fp, SeqFeatPtr gene, SeqLocPtr oldloc)
{
  CharPtr           loc1, loc2;
  SeqMgrFeatContext gcontext;
  SeqFeatPtr        found_gene;

  if (fp == NULL || gene == NULL || gene->location == NULL || oldloc == NULL) return;
  found_gene = SeqMgrGetDesiredFeature (gene->idx.entityID, NULL, 0, 0, gene, &gcontext);
  if (found_gene == NULL) return;
  loc1 = SeqLocPrint (gene->location);
  loc2 = SeqLocPrint (oldloc);
  if (loc1 != NULL && loc2 != NULL && StringCmp (loc1, loc2) != 0) {
    fprintf (fp, "Moved gene %s from %s to %s\n", gcontext.label, loc2, loc1);
  }
  MemFree (loc1);
  MemFree (loc2);
}


static void LogGeneBadCorrection (FILE *fp, SeqFeatPtr gene)
{
  SeqMgrFeatContext gcontext;
  SeqFeatPtr        found_gene;

  if (fp == NULL || gene == NULL || gene->location == NULL) return;
  found_gene = SeqMgrGetDesiredFeature (gene->idx.entityID, NULL, 0, 0, gene, &gcontext);
  if (found_gene == NULL) return;
  fprintf (fp, "Please examine gene %s, new gene may be too large\n", gcontext.label);
}


static void CorrectCDSGene
(SeqFeatPtr sfp,
 SeqFeatPtr gene_sfp,
 FILE *fp,
 Boolean PNTR data_in_log)
{
  BioseqPtr bsp;
  Boolean   partial5, partial3;
  SeqLocPtr log_slp, new_slp;
  Int2      res;

  if (sfp == NULL || gene_sfp == NULL || fp == NULL 
      || data_in_log == NULL){
    return;
  }
  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return;
  new_slp = SeqLocMerge (bsp, sfp->location, NULL, TRUE, FALSE, FALSE);
  res = SeqLocCompare (gene_sfp->location, new_slp);
  if (res == SLC_A_EQ_B) {
    SeqLocFree (new_slp);
    return;
  }
  if (SeqLocLen (gene_sfp->location) > SeqLocLen (new_slp)) {
    SeqLocFree (new_slp);
    new_slp = SeqLocMerge (bsp, sfp->location, gene_sfp->location, TRUE, FALSE, FALSE);
    LogGeneBadCorrection (fp, gene_sfp);
  }
  log_slp = SeqLocMerge (bsp, sfp->location, NULL, FALSE, FALSE, FALSE);

  *data_in_log = TRUE;
  SeqLocFree (gene_sfp->location);
  gene_sfp->location = new_slp;
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  SetSeqLocPartial (gene_sfp->location, partial5, partial3);
  gene_sfp->partial = sfp->partial;
  LogGeneCorrection (fp, gene_sfp, log_slp);
  SeqLocFree (log_slp);
}

static void FixProtRefPtrs (ValNodePtr feat_list)
{
  ValNodePtr vnp;
  SeqFeatPtr sfp;
  Uint1      strand;
  BioseqPtr  bsp;

  if (feat_list == NULL) return;
  for (vnp = feat_list; vnp != NULL; vnp = vnp->next) {
    sfp = (SeqFeatPtr)vnp->data.ptrvalue;
    if (sfp == NULL) continue;
    bsp = BioseqFindFromSeqLoc (sfp->location);
    if (bsp == NULL) continue;
    strand = SeqLocStrand (sfp->location);
    if (strand == Seq_strand_minus) {
      sfp->location = ExpandSeqLoc (bsp->length - 1, 0, strand, bsp, sfp->location);
    } else {
      sfp->location = ExpandSeqLoc (0, bsp->length - 1, strand, bsp, sfp->location);
    }
  }
}


static Boolean PrepareUpdateAlignmentForProtein
(SeqFeatPtr sfp,
 BioseqPtr  protBsp,
 Uint2      input_entityID,
 FILE *     fp,
 Boolean    truncate_proteins,
 Boolean    extend_proteins5,
 Boolean    extend_proteins3,
 Boolean    correct_cds_genes,
 Int4       transl_except_len,
 Boolean    *data_in_log,
 SeqAlignPtr PNTR salpptr,
 BioseqPtr PNTR newbspptr)
{
  Boolean       align_succeeded;
  Boolean       changed_frame;
  Uint1         frame_attempts;
  ErrSev        level;
  CdRegionPtr   crp;
  SeqLocPtr     orig_slp;
  Boolean       orig_partial;
  SeqFeatPtr    gene_sfp;
  SeqMgrFeatContext gcontext;
  Boolean       truncated, stop_shifted;
  Boolean       extended5, extended3;
  Uint1         old_frame;
  BioseqPtr     newBsp = NULL;
  Boolean       adjust_succeeded;
  SeqAlignPtr   salp = NULL;
  Boolean       revcomp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION
      || (crp = (CdRegionPtr)sfp->data.value.ptrvalue) == NULL
      || protBsp == NULL
      || data_in_log == NULL) {
    return FALSE;
  }

  old_frame = crp->frame;

  orig_slp = sfp->location;
  orig_partial = sfp->partial;
  if (correct_cds_genes) {
    gene_sfp = SeqMgrGetOverlappingGene (sfp->location, &gcontext);
  }
  align_succeeded = FALSE;
  changed_frame = FALSE;

  level = ErrSetMessageLevel (SEV_MAX);
  for (frame_attempts = 1;
       frame_attempts < 4 && ! align_succeeded;
       frame_attempts ++) {
    if (sfp->location != orig_slp) {
      if (sfp->location->choice == 0) {
        SeqLocFree (sfp->location);
      }
      sfp->location = orig_slp;
    }
    sfp->partial = orig_partial;    
    crp->frame = old_frame;
    if (newBsp != NULL) {
      newBsp = BioseqFree (newBsp);
    }
    truncated= FALSE;
    stop_shifted = FALSE;
    extended5 = FALSE;
    extended3 = FALSE;
    adjust_succeeded = AdjustCDSForUpdate (sfp,
                                           input_entityID,
                                           frame_attempts,
                                           transl_except_len,
                                           truncate_proteins,
                                           extend_proteins5,
                                           extend_proteins3,
                                           &truncated, &stop_shifted,
                                           &extended5, &extended3,
                                           protBsp, &newBsp);
    if (sfp->location == NULL || sfp->location->choice == 0) 
    {
      sfp->location = orig_slp;
    }
    if (adjust_succeeded)
    {
      salp = Sequin_GlobalAlign2Seq (protBsp, newBsp, &revcomp);
    }
    if (!adjust_succeeded || salp != NULL)
    {
      if (frame_attempts > 1) {
        changed_frame = TRUE;
        LogFrameChange (fp, protBsp, frame_attempts);
        *data_in_log = TRUE;
      }
      break;
    }    
  }
  ErrSetMessageLevel (level);
  if (adjust_succeeded && salp == NULL) 
  {
    /* put CD Region back to original state */
    crp->frame = old_frame;
    if (sfp->location != orig_slp) {
      SeqLocFree (sfp->location);
      sfp->location = orig_slp;
    }
    sfp->partial = orig_partial;
    WarnNoProteinUpdate (protBsp);
    newBsp = BioseqFree (newBsp);
    return FALSE;
  }

  if (truncated) {
    LogProteinTruncate (fp, protBsp);
    *data_in_log = TRUE;
  }
  if (stop_shifted) {
    LogProteinStopShift (fp, protBsp);
    *data_in_log = TRUE;
  }
  if (extended5) {
    LogProteinExtend5 (fp, protBsp);
    *data_in_log = TRUE;
  }
  if (extended3) {
    LogProteinExtend3 (fp, protBsp);
    *data_in_log = TRUE;
  }
  
  if (correct_cds_genes && gene_sfp != NULL) {
    CorrectCDSGene (sfp, gene_sfp, fp, data_in_log);
  }
  if (sfp->location != orig_slp) {
    SeqLocFree (orig_slp);
  }
  
  *salpptr = salp;
  *newbspptr = newBsp;
  return TRUE;
}

static Boolean CodingRegionHasTranslExcept (SeqFeatPtr sfp)
{
  CodeBreakPtr cbp;
  Int4         len;
  CdRegionPtr  crp;
  SeqLocPtr    slp;
  Int4         codon_start, codon_stop, pos, codon_length;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION
      || (crp = (CdRegionPtr)sfp->data.value.ptrvalue) == NULL
      || crp->code_break == NULL)
  {
  	return FALSE;
  }

  len = SeqLocLen (sfp->location);
  if (crp->frame > 1) 
  {
  	len -= crp->frame - 1;
  }
  if (len % 3 == 0) 
  {
  	return FALSE;
  }
  for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next)
  {
    codon_start = INT4_MAX;
    codon_stop = -10;
    slp = NULL;
    while ((slp = SeqLocFindNext (cbp->loc, slp)) != NULL) {
      pos = GetOffsetInLoc (slp, sfp->location, SEQLOC_START);
      if (pos < codon_start)
      {
        codon_start = pos;
        pos = GetOffsetInLoc (slp, sfp->location, SEQLOC_STOP);
        if (pos > codon_stop)
        {
          codon_stop = pos;
        }
        codon_length = codon_stop - codon_start;      /* codon length */
        if (codon_length >= 0 && codon_length <= 1 && codon_stop == len - 1)
        {                       /*  a codon */
          /* allowing a partial codon at the end */
          return TRUE;
        }
      }
    }
  }
  return FALSE;
}

typedef struct proteinfromcdsdata {
 Uint2      input_entityID;
 Uint4      input_itemID;
 Uint4      input_itemtype;
 FILE *     fp;
 Boolean    data_to_report;
} ProteinFromCDSData, PNTR ProteinFromCDSPtr;

static void UpdateOneProteinFromCDS (SeqFeatPtr sfp, Pointer userdata)
{
  ProteinFromCDSPtr pfcp;
  SeqLocPtr         new_product;
  Boolean           data_in_log;
  Int4              transl_except_len = 0;
  BioseqPtr         protBsp;
  SeqAlignPtr       salp = NULL;
  BioseqPtr         newbsp = NULL;
  Boolean           rval;
  
  pfcp = (ProteinFromCDSPtr) userdata;
  if ( pfcp == NULL
    || sfp == NULL
    || sfp->idx.subtype != FEATDEF_CDS
    || sfp->product == NULL)
  {
    return;
  }
  
  protBsp = BioseqFindFromSeqLoc (sfp->product);
  if (protBsp == NULL)
  {
    return;
  }
 
  data_in_log = FALSE;
  
  if (CodingRegionHasTranslExcept(sfp))
  {
  	transl_except_len = SeqLocLen (sfp->location) % 3;
  }
  
  rval = PrepareUpdateAlignmentForProtein (sfp,
                                           protBsp,
                                           pfcp->input_entityID,
                                           pfcp->fp,
                                           TRUE, TRUE, TRUE, TRUE,
                                           transl_except_len,
                                           &data_in_log,
                                           &salp,
                                           &newbsp);
  if (data_in_log) {
    pfcp->data_to_report = TRUE;
  }
  if (!rval) return;

  ReplaceOneSequence (salp, protBsp, newbsp);
  if (sfp->product->choice != SEQLOC_WHOLE)
  {
    new_product = SeqLocWholeNew (protBsp);
    if (new_product == NULL) return;
    SeqLocFree (sfp->product);
    sfp->product = new_product;
  }
}

static Boolean CheckForIDCollision (
  BioseqPtr oldbsp,
  BioseqPtr newbsp,
  BoolPtr islocal
)

{
  SeqIdPtr  sip;

  if (oldbsp == NULL || newbsp == NULL) return FALSE;
  for (sip = newbsp->id; sip != NULL; sip = sip->next) {
    if (SeqIdIn (sip, oldbsp->id)) {
      if (sip->choice == SEQID_LOCAL) {
        *islocal = TRUE;
      }
      return TRUE;
    }
  }
  return FALSE;
}

static CharPtr convPubDescMssg =
"Do you wish to convert publications to apply only to the appropriate ranges?";

#define UPDATE_SKIP_THIS    1
#define UPDATE_SKIP_ALL     2
#define UPDATE_REPLACE_THIS 3
#define UPDATE_REPLACE_ALL  4
#define UPDATE_CANCEL       5

typedef struct noalignmentchoice
{
  Boolean done;
  Boolean cancelled;
  GrouP   action_choice;
  Boolean skip_this;
  Boolean skip_all;
  Boolean replace_this;
  Boolean replace_all;
} NoAlignmentChoiceData, PNTR NoAlignmentChoicePtr;

static void NoAlignmentChoiceOk (ButtoN b)
{
  NoAlignmentChoicePtr nacp;
  
  nacp = (NoAlignmentChoicePtr) GetObjectExtra (b);
  if (nacp == NULL) return;
  nacp->cancelled = FALSE;
  nacp->done = TRUE;
}

static void NoAlignmentChoiceCancel (ButtoN b)
{
  NoAlignmentChoicePtr nacp;
  
  nacp = (NoAlignmentChoicePtr) GetObjectExtra (b);
  if (nacp == NULL) return;
  nacp->cancelled = TRUE;
  nacp->done = TRUE;
}

static Int2 GetNoAlignmentChoice (SeqIdPtr id, Int2 previous_choice)
{
  Char                  buf [64];
  WindoW                w;
  GrouP                 h, g, c;
  ButtoN                b;
  NoAlignmentChoiceData nacd;
  Char                  promptstr[115];
  Int2                  no_aln_choice;
  SeqIdPtr              sip, sip_next;

  sip = SeqIdFindBest (id, 0);
  if (sip == NULL) return UPDATE_CANCEL;
  w = ModalWindow(-20, -13, -10, -10, NULL);
  if (w == NULL) return 0;

  h = HiddenGroup (w, -1, 0, NULL);
  g = HiddenGroup(h, 0, 4, NULL);
  sip_next = sip->next;
  sip->next = NULL;
  SeqIdWrite (sip, buf, PRINTID_REPORT, sizeof (buf) - 1);
  sip->next = sip_next;
  sprintf (promptstr, "There is no alignment between the sequences for %s.", buf);
  StaticPrompt (g, promptstr, 0, popupMenuHeight, programFont, 'l');
  StaticPrompt (g, "You may choose to :", 0, popupMenuHeight, programFont, 'l');
  nacd.action_choice = HiddenGroup (g, 0, 4, NULL);
  RadioButton (nacd.action_choice, "Skip This Sequence");
  RadioButton (nacd.action_choice, "Skip All Sequences without Alignments");
  RadioButton (nacd.action_choice, "Replace This Sequence");
  RadioButton (nacd.action_choice, "Replace All Sequences without Alignments");
  if (previous_choice > 0 && previous_choice < UPDATE_CANCEL)
  {
    SetValue (nacd.action_choice, previous_choice);
  }
  else
  {
    SetValue (nacd.action_choice, 1);
  }
  c = HiddenGroup (h, 2, 0, NULL);
  b = PushButton (c, "Ok", NoAlignmentChoiceOk);
  SetObjectExtra (b, &nacd, NULL);
  b = PushButton (c, "Cancel", NoAlignmentChoiceCancel);
  SetObjectExtra (b, &nacd, NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);

  nacd.cancelled = FALSE;
  nacd.done = FALSE;
  Show(w); 
  Select (w);
  while (!nacd.done)
  {
    ProcessExternalEvent ();
    Update ();
  }
  ProcessAnEvent ();
  no_aln_choice = GetValue (nacd.action_choice);
  Remove (w);
  if (nacd.cancelled)
  {
    return UPDATE_CANCEL;
  }
  else
  {
    return no_aln_choice;
  }
}

static Boolean PrepareUpdatePtr (UpsDataPtr    udp)
{
  Uint2        entityID;
  SeqEntryPtr  oldsep, newsep;
  SeqIdPtr     sip;
  Boolean      islocal = FALSE;
  Char         buf [64];
  SeqAlignPtr  salp = NULL;
  Boolean      revcomp = FALSE;
  Boolean      asked_about_desc_prop = FALSE;
  Boolean      propagate_descriptors = FALSE;
  SeqIdPtr     id, id_next;
  MsgAnswer    ans;
  Char         collision_id [100];

  if (udp->oldbsp == NULL || udp->newbsp == NULL) return FALSE;
  if (ISA_na (udp->oldbsp->mol) != ISA_na (udp->newbsp->mol)) {
    Message (MSG_OK, "Both sequences must be either nucleotides or proteins");
    return FALSE;
  }

  entityID = ObjMgrGetEntityIDForPointer (udp->oldbsp);
  oldsep = GetBestTopParentForData (entityID, udp->oldbsp);
  entityID = ObjMgrGetEntityIDForPointer (udp->newbsp);
  newsep = GetBestTopParentForData (entityID, udp->newbsp);
  if (oldsep == NULL || newsep == NULL)
    return FALSE;

  if (CONVERTPUBS_NOT_SET == udp->convertPubs) 
  {
    if (Message (MSG_YN, convPubDescMssg) == ANS_YES) {
      ConvertPubSrcComDescsToFeats (oldsep, TRUE, FALSE, FALSE, FALSE, &asked_about_desc_prop, &propagate_descriptors, NULL);
      ConvertPubSrcComDescsToFeats (newsep, TRUE, FALSE, FALSE, FALSE, &asked_about_desc_prop, &propagate_descriptors, NULL);
      udp->convertPubs = CONVERTPUBS_YES;
    }
    else
      udp->convertPubs = CONVERTPUBS_NO;
  }

  if (CheckForIDCollision (udp->oldbsp, udp->newbsp, &islocal)) {
    sprintf (collision_id, "lcl|SequinUpdateSequence%d", udp->seqsubpos);
    sip = SeqIdParse (collision_id);
    if (sip != NULL) {
      BioseqReplaceID (udp->newbsp, sip);
      sip = SeqIdFree (sip);
    }
  }
  
  salp = Sequin_GlobalAlign2Seq (udp->oldbsp, udp->newbsp, &revcomp);

  if (salp == NULL) {
    if (udp->log_fp != NULL) {
      id = SeqIdFindBest (udp->oldbsp->id, 0);
      if (id != NULL)
      {
        id_next = id->next;
        id->next = NULL;
        SeqIdWrite (id, buf, PRINTID_REPORT, sizeof (buf) - 1);
        id->next = id_next;
        fprintf (udp->log_fp, "No sequence similarity for %s\n", buf);
        udp->data_in_log = TRUE;  
      }
    }
    if (udp->suppress_continue_msg) {
      return FALSE;
    } else {
      if (udp->isSet && udp->do_update)
      {
        if (udp->no_aln_choice != UPDATE_SKIP_ALL && udp->no_aln_choice != UPDATE_REPLACE_ALL)
        {
          udp->no_aln_choice = GetNoAlignmentChoice (udp->oldbsp->id, udp->no_aln_choice);
        }
        if (udp->no_aln_choice == UPDATE_CANCEL)
        {
          return FALSE;
        }
      }
      else if (udp->do_update)
      {
        ans = Message (MSG_YN, "There is no alignment between the sequences.  Do you wish to continue?");
        if (ans == ANS_YES)
        {
          udp->useGUI = TRUE;
        }
        else
        {
          return FALSE;
        }
      }
    }
  }
  udp->salp     = salp;
  udp->revcomp  = revcomp;
  udp->diffOrgs = FALSE;
  udp->recomb1  = -1;
  udp->recomb2  = -1;
  
  return TRUE;
}

static void UpdateProteinsOnNewBsp (SeqFeatPtr sfp, Pointer userdata);

static void
FindProtRefPtrsOnBsp
(BioseqPtr bsp,
 ValNodePtr PNTR feat_list)
{
  SeqMgrFeatContext context;
  SeqFeatPtr        sfp;
  ValNodePtr        vnp;
 
  if (bsp == NULL || feat_list == NULL) return;
  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_PROT, 0, &context);
  while (sfp != NULL) {
    if (context.left == 0 && context.right == bsp->length - 1) {
      vnp = ValNodeNew (*feat_list);
      if (vnp != NULL) {
        vnp->data.ptrvalue = sfp;
      }
      if (*feat_list == NULL) {
        *feat_list = vnp;
      }
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_PROT, 0, &context);
  }
}

static ValNodePtr FindProductProtRefs (BioseqPtr bsp)
{
  SeqFeatPtr        sfp;
  SeqMgrFeatContext context;
  ValNodePtr        feat_list = NULL;

  if (bsp == NULL) return NULL;
  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &context);
  while (sfp != NULL) {
    FindProtRefPtrsOnBsp (BioseqFindFromSeqLoc (sfp->product), &feat_list);
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &context);
  }
  return feat_list; 
}

static ValNodePtr FindTranslExceptCDSs (BioseqPtr bsp)
{
  SeqFeatPtr        sfp;
  SeqMgrFeatContext context;
  ValNodePtr        feat_list = NULL;
  ValNodePtr        vnp;
  Int4              cd_len;

  if (bsp == NULL) return NULL;
  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &context);
  while (sfp != NULL) {
    if (CodingRegionHasTranslExcept (sfp)) 
    {
      vnp = ValNodeNew(feat_list);
      if (feat_list == NULL) 
      {
      	feat_list = vnp;
      }
      if (vnp != NULL)
      {
      	cd_len = SeqLocLen (sfp->location);
      	vnp->choice = cd_len % 3;
      	vnp->data.ptrvalue = sfp;
      }
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &context);
  }
  return feat_list; 
}

static Int4 GetOriginalTranslExceptLen (SeqFeatPtr sfp, ValNodePtr list)
{
  ValNodePtr vnp;
  
  if (sfp == NULL || list == NULL) return 0;
  for (vnp = list; vnp != NULL; vnp = vnp->next)
  {
  	if (vnp->data.ptrvalue == sfp) 
  	{
  		return vnp->choice;
  	}
  }
  if (CodingRegionHasTranslExcept (sfp)) 
  {
  	return SeqLocLen (sfp->location) % 3;
  }
  return 0;
}

static void UpdateOneSequence (
  UpsDataPtr   udp,
  Int2         sfbval,
  Boolean      add_cit_subs,
  Boolean      update_proteins
)
{
  SeqAnnotPtr  sap = NULL;
  Uint2        entityID;
  SeqEntryPtr  sep;
  Boolean      update = FALSE;
  Boolean      feature_update = FALSE;
  ValNodePtr   prot_feat_list = NULL;

  if (udp != NULL)
  {
    if (GetStatus (udp->replace_all))
    {
  	  RemoveOldFeats (udp->oldbsp);
    }
  }

  if (update_proteins && udp != NULL) {
    prot_feat_list = FindProductProtRefs (udp->oldbsp);
    if (udp->transl_except_list != NULL)
    {
      ValNodeFree (udp->transl_except_list);
      udp->transl_except_list = NULL;
    }
    udp->transl_except_list = FindTranslExceptCDSs (udp->oldbsp);
  }

  if (sfbval == UPDATE_SEQUENCE_ONLY || sfbval == UPDATE_SEQUENCE_AND_FEATURES) {
    switch (udp->rmcval) {
      case UPDATE_REPLACE :
        if (ReplaceSequence (udp)) {
          update = TRUE;
        }
        break;
      case UPDATE_EXTEND5 :
	      if (NULL == udp->salp) {
	        if (Merge5PrimeNoOverlap (udp))
	          update = TRUE;
	      }
	      else 
	      {
	        if (Merge5Prime (udp))
	          update = TRUE;
        }
        break;
      case UPDATE_EXTEND3 :
	      if (NULL == udp->salp) {
	        if (Merge3PrimeNoOverlap (udp))
	          update = TRUE;
	      }
	      else 
	      {
	        if (Merge3Prime (udp))
	          update = TRUE;
        }
        break;
      case UPDATE_PATCH :
        if (PatchSequence (udp)) {
          update = TRUE;
        }
        break;
      default :
        break;
    }
    if ( sfbval == UPDATE_SEQUENCE_AND_FEATURES) {
      switch (udp->rmcval) {
        case UPDATE_REPLACE :
          if (DoFeaturePropWithOffset (udp, 0, &sap, FALSE)) {
            update = TRUE;
            feature_update = TRUE;
          }
          break;
        case UPDATE_EXTEND5 :
          if (DoFeaturePropWithOffset (udp, 0, &sap, FALSE)) {
            update = TRUE;
            feature_update = TRUE;
          }
          break;
        case UPDATE_EXTEND3 :
          if (DoFeaturePropWithOffset (udp, udp->old5 - udp->new5, &sap, FALSE)) {
            update = TRUE;
            feature_update = TRUE;
          }
          break;
        case UPDATE_PATCH :
          if (DoFeaturePropWithOffset (udp, udp->old5 - udp->new5, &sap, TRUE)) {
            update = TRUE;
            feature_update = TRUE;
          }
          break;
        default :
          break;
      }
    }
    if (udp->rmcval == UPDATE_REPLACE && udp->revcomp && update)
    {
      ReverseComplementBioseqAndFeats (udp->oldbsp, udp->input_entityID);
    }

  } else if (sfbval == UPDATE_FEATURES_ONLY) {
    switch (udp->rmcval) {
      case UPDATE_REPLACE :
        if (DoFeaturePropThruAlign (udp, &sap)) {
          update = TRUE;
          feature_update = TRUE;
        }
        break;
      case UPDATE_EXTEND5 :
        if (DoFeaturePropThruAlign (udp, &sap)) {
          update = TRUE;
          feature_update = TRUE;
        }
        break;
      case UPDATE_EXTEND3 :
        if (DoFeaturePropThruAlign (udp, &sap)) {
          update = TRUE;
          feature_update = TRUE;
        }
        break;
      case UPDATE_PATCH :
        if (DoFeaturePropThruAlign (udp, &sap)) {
          update = TRUE;
          feature_update = TRUE;
        }
        break;
      default :
        break;
    }
  }
  if (update) {
    entityID = ObjMgrGetEntityIDForPointer (udp->oldbsp);
    if (add_cit_subs
      && (feature_update || StringCmp (udp->seq1, udp->seq2) != 0))
    {
      AddCitSubToUpdatedSequence ( udp->oldbsp, entityID);
    }
    if (update_proteins)
    {
      SeqMgrClearFeatureIndexes (entityID, udp->oldbsp);
      SeqMgrIndexFeatures (entityID, NULL);
      sep = GetBestTopParentForData (entityID, udp->oldbsp);
      udp->truncate_proteins = GetStatus (udp->truncate_proteins_btn);
      udp->extend_proteins5 = GetStatus (udp->extend_proteins5_btn);
      udp->extend_proteins3 = GetStatus (udp->extend_proteins3_btn);
      udp->correct_cds_genes = GetStatus (udp->correct_cds_genes_btn);
      VisitFeaturesInSep (sep, udp, UpdateProteinsOnNewBsp);
      FixProtRefPtrs (prot_feat_list);
    }
    if (! udp->suppress_instant_refresh) {
      ObjMgrSetDirtyFlag (entityID, TRUE);
      ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
    }
  }
  
  if (GetStatus (udp->update_quality_scores_btn) && udp->rmcval == UPDATE_REPLACE 
      && (sfbval == UPDATE_SEQUENCE_ONLY || sfbval == UPDATE_SEQUENCE_AND_FEATURES))
  {
    ReplaceQualityScores (udp->oldbsp, udp->newbsp, udp->salp, udp->log_fp, &(udp->data_in_log));
  }
  
  ValNodeFree (prot_feat_list);
  ValNodeFree (udp->transl_except_list);
  udp->transl_except_list = NULL;
  if (sfbval == UPDATE_FEATURES_ONLY || sfbval == UPDATE_SEQUENCE_AND_FEATURES) {
    if (update) {
      entityID = ObjMgrGetEntityIDForPointer (udp->oldbsp);
      sep = GetTopSeqEntryForEntityID (entityID);
      /* need to set scope to make sure we mark the right bioseq for deletion */
      SeqEntrySetScope (sep);
      /* resolve features unless the policy was to remove all the old ones */
      if (!GetStatus (udp->replace_all))
      {
        ResolveDuplicateFeats (udp, udp->oldbsp, sap);      	
      }
      SeqEntrySetScope (NULL);
      DeleteMarkedObjects (entityID, 0, NULL);
      SeqMgrClearFeatureIndexes (entityID, NULL);
      ObjMgrSetDirtyFlag (entityID, TRUE);
      ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
    }
  }
}

static void UpdateProteinsOnNewBsp (SeqFeatPtr sfp, Pointer userdata)
{
  UpsDataPtr    udp_orig;
  SeqLocPtr     new_product;
  Boolean       data_in_log;
  Int4          transl_except_len = 0;
  Boolean       fix_products = TRUE;
  BioseqPtr     protBsp = NULL;
  SeqAlignPtr   salp = NULL;
  BioseqPtr     newbsp = NULL;
  Boolean       rval;

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_CDS || userdata == NULL)
  {
    return;
  }
  
  if (sfp->idx.deleteme)
  {
    return;
  }
  
  protBsp = BioseqFindFromSeqLoc (sfp->product);
  if (protBsp == NULL)
  {
    return;
  }
  udp_orig = (UpsDataPtr) userdata;

  data_in_log = FALSE;
  transl_except_len = GetOriginalTranslExceptLen (sfp, udp_orig->transl_except_list);
  rval = PrepareUpdateAlignmentForProtein (sfp,
                                           protBsp,
                                           udp_orig->input_entityID,
                                           udp_orig->log_fp,
                                           udp_orig->truncate_proteins,
                                           udp_orig->extend_proteins5,
                                           udp_orig->extend_proteins3,
                                           udp_orig->correct_cds_genes,
                                           transl_except_len,
                                           &data_in_log,
                                           &salp,
                                           &newbsp);
  if (data_in_log) {
    udp_orig->data_in_log = TRUE;
  }
  if (!rval) return;
  
  if (protBsp->idx.deleteme)
  {
    fix_products = FALSE;
  }

  ReplaceOneSequence (salp, protBsp, newbsp);
  
  if (fix_products)
  {
    if (sfp->product->choice != SEQLOC_WHOLE) {
      new_product = SeqLocWholeNew (protBsp);
      if (new_product == NULL) return;
      SeqLocFree (sfp->product);
      sfp->product = new_product;
    }
    newbsp = BioseqFree (newbsp);
  }
}

extern void UpdateProteinsFromCDS ( IteM i)
{
  BaseFormPtr        bfp;
  SeqEntryPtr        sep;
  ProteinFromCDSData pfcd;
  Char               path [PATH_MAX];

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif

  if (bfp == NULL)
    return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  pfcd.input_entityID = bfp->input_entityID;
  pfcd.input_itemID = bfp->input_itemID;
  pfcd.input_itemtype = bfp->input_itemtype;
  pfcd.data_to_report = FALSE;
  TmpNam (path);
  pfcd.fp = FileOpen (path, "wb");
  if (pfcd.fp == NULL) return;

  WatchCursor ();
  VisitFeaturesInSep (sep, (Pointer) &pfcd, UpdateOneProteinFromCDS);
  FileClose (pfcd.fp);
  if (pfcd.data_to_report) {
    LaunchGeneralTextViewer (path, "Update Log");
  }
  FileRemove (path);
  ArrowCursor ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

/* This section of code is used to warn the user when an alignment update
 * will affect the area inside a variation or the area immediately to the
 * left or right of a variation.
 */
static void VariationAlignmentCallback (SeqFeatPtr sfp, Pointer userdata)
{
  CharPtr    buf, cp;
  UpsDataPtr udp;
  Boolean    change_found = FALSE;
  
  if (sfp == NULL || userdata == NULL || sfp->idx.subtype != FEATDEF_variation)
  {
    return;
  }
  udp = (UpsDataPtr) userdata;
  if (udp->salp == NULL)
  {
    return;
  }
  
  buf = FeatureLocationAlignment (sfp, udp->salp, 1, 2);
  if (buf == NULL)
  {
    return;
  }
  cp = buf;
  /* skip over "Old    :" */
  if (StringLen (cp) > 8)
  {
    cp += 8;
  }
  while (*cp != 0 && *cp != '\n' && !change_found)
  {
    if (*cp == '-')
    {
      change_found = TRUE;
    }
    cp++;
  }
  if (*cp == '\n')
  {
    cp++;
  }
  /* skip over "New    :" */
  if (StringLen (cp) > 8)
  {
    cp += 8;
  }
  while (*cp != 0 && *cp != '\n' && !change_found)
  {
    if (*cp != '.')
    {
      change_found = TRUE;
    }
    cp++;
  }
  if (change_found)
  {
    ValNodeAddPointer (&udp->affected_variation_features, 0, sfp);
  }
  buf = MemFree (buf);
}

static CharPtr PrepareVariationFeatureLogEntry (SeqFeatPtr sfp, SeqAlignPtr salp)
{
  CharPtr loc_buf, feat_buf, total_buf;
  Int4    buf_len;
  
  if (sfp == NULL || salp == NULL)
  {
    return NULL;
  }
  
  loc_buf = SeqLocPrint (sfp->location);
  feat_buf = FeatureLocationAlignment (sfp, salp, 1, 2);
    
  buf_len = StringLen (sfp->comment) + StringLen (loc_buf) + StringLen (feat_buf) + 3;
  total_buf = (CharPtr) MemNew (buf_len * sizeof (Char));
  if (total_buf != NULL)
  {
    total_buf [0] = 0;
    /* add feature name */
    if (!StringHasNoText (sfp->comment))
    {
      StringCpy (total_buf, sfp->comment);
      StringCat (total_buf, "\n");
    }
    /* add feature location */
    if (!StringHasNoText (loc_buf))
    {
      StringCat (total_buf, loc_buf);
      StringCat (total_buf, "\n");
    }
    /* add alignment picture */
    if (!StringHasNoText (feat_buf))
    {
      StringCat (total_buf, feat_buf);
    }
    loc_buf = MemFree (loc_buf);
    feat_buf = MemFree (feat_buf);
  }
  return total_buf;
}

static Int4 
VariationFeatureChangesOk 
(ValNodePtr  variation_feature_list,
 SeqAlignPtr salp,
 Boolean     allow_skip)
{
  WindoW w;
  GrouP  h, c;
  ButtoN b;
  DoC    doc;
  ValNodePtr vnp;
  SeqFeatPtr sfp;
  CharPtr               total_buf;
  ModalAcceptCancelData acd;
  CharPtr               str_format = "%d variation features are affected by this update.";
  
  if (variation_feature_list == NULL || salp == NULL)
  {
    return TRUE;
  }
  
  w = MovableModalWindow (-20, -13, -10, -10, "Affected Variation Sequences", NULL);
  h = HiddenGroup(w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  doc = DocumentPanel (h, stdCharWidth * 40, stdLineHeight * 12);
  SetDocAutoAdjust (doc, TRUE);
  
  total_buf = MemNew (StringLen (str_format) + 15);
  if (total_buf != NULL)
  {
    sprintf (total_buf, str_format, ValNodeLen (variation_feature_list));
    AppendText (doc, total_buf, NULL, NULL, programFont);
    total_buf = MemFree (total_buf);
  }
  
  for (vnp = variation_feature_list; vnp != NULL; vnp = vnp->next)
  {
    sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    total_buf = PrepareVariationFeatureLogEntry (sfp, salp);
    if (!StringHasNoText (total_buf))
    {
      AppendText (doc, total_buf, NULL, NULL, programFont);
    }
    MemFree (total_buf);
  }
  UpdateDocument (doc, 0, 0);
  
  c = HiddenGroup (h, 3, 0, NULL);
  b = PushButton (c, "Proceed", ModalAcceptButton);
  SetObjectExtra (b, &acd, NULL);
  if (allow_skip)
  {
    b = PushButton (c, "Skip Update", ModalThirdOptionButton);
    SetObjectExtra (b, &acd, NULL);
  }
  
  b = PushButton (c, "Cancel Update", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
 
  
  AlignObjects (ALIGN_CENTER, (HANDLE) doc,
                              (HANDLE) c, 
                              NULL);

  Show (w);
  Select (w);
  
  acd.cancelled = FALSE;
  acd.third_option = FALSE;
  acd.accepted = FALSE;
  while (!acd.accepted && ! acd.cancelled && !acd.third_option)
  {
    ProcessExternalEvent ();
    Update ();
  }
  ProcessAnEvent ();
  Remove (w);
  
  if (acd.accepted)
  {
    return 1;
  }
  else if (acd.cancelled)
  {
    return 2;
  }
  else
  {
    return 3;
  }
}

static void 
AddVariationFeaturesToLog 
(FILE       *fp, 
 BoolPtr    data_in_log,
 ValNodePtr variation_feature_list,
 SeqAlignPtr salp)
{
  ValNodePtr vnp;
  CharPtr    total_buf;
  SeqFeatPtr sfp;
  
  if (fp == NULL || variation_feature_list == NULL)
  {
    return;
  }
  
  for (vnp = variation_feature_list; vnp != NULL; vnp = vnp->next)
  {
    sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    total_buf = PrepareVariationFeatureLogEntry (sfp, salp);
    if (!StringHasNoText (total_buf))
    {
      fprintf (fp, total_buf);
      if (data_in_log != NULL)
      {
        *data_in_log = TRUE;
      }
    }
    MemFree (total_buf);
  }
  
}

static void AcceptRMCOrExtend (ButtoN b)
{
  UpsDataPtr   udp;
  Int2         sfbval;
  Boolean      log_is_local;
  Boolean      update_proteins = FALSE;
  Uint2        entityID;

  udp = (UpsDataPtr) GetObjectExtra (b);
  if (udp == NULL) return;
  SafeHide (udp->form);

  sfbval = GetValue (udp->sfb);
  udp->rmcval = GetValue (udp->rmc);

  if (udp->log_fp == NULL) {
    OpenSequenceUpdateLog (udp);
    if (udp->log_fp == NULL) return;
    log_is_local = TRUE;
  } else {
    log_is_local = FALSE;
  }
  WatchCursor ();
  Update ();
  
  if (udp->do_update)
  {
    udp->affected_variation_features = ValNodeFree (udp->affected_variation_features);
    VisitFeaturesOnBsp (udp->oldbsp, udp, VariationAlignmentCallback);
    if (1 == VariationFeatureChangesOk (udp->affected_variation_features, udp->salp, FALSE))
    {
      if (udp->update_proteins != NULL 
        && Enabled (udp->update_proteins)
        && GetStatus (udp->update_proteins))
      {
        update_proteins = TRUE;
      }
      AddVariationFeaturesToLog (udp->log_fp, &(udp->data_in_log), 
                                 udp->affected_variation_features, udp->salp);
      UpdateOneSequence (udp, sfbval,
                         GetStatus (udp->add_cit_subs),
                         update_proteins);
    }
    udp->affected_variation_features = ValNodeFree (udp->affected_variation_features);
  }
  else
  {
    if (ExtendOneSequence (udp))
    {
      entityID = ObjMgrGetEntityIDForPointer (udp->oldbsp);
      ObjMgrSetDirtyFlag (entityID, TRUE);
      ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
    }
  }
  
  Remove (udp->form);
  if (log_is_local) {
    CloseOutSequenceUpdateLog (udp);
  }
  ArrowCursor ();
  Update();
}

static void DoAcceptRMCOrExtendSet (UpsDataPtr udp)
{
  Int2         sfbval;
  Boolean      log_is_local = FALSE;
  Boolean      update_proteins = FALSE;
  Boolean      do_update = TRUE;
  Boolean      old_useGUI;
  Char         acc_str [256];
  SeqIdPtr     sip, sip_next;
  Int2         prior_rmcval;
  Int4         variation_action = 1;

  if (udp == NULL) return;
  
  SafeHide (udp->form);
  old_useGUI = udp->useGUI;
  
  prior_rmcval = udp->rmcval;
  if (udp->do_update)
  {
    if (udp->salp == NULL 
        && (udp->no_aln_choice == UPDATE_REPLACE_THIS
         || udp->no_aln_choice == UPDATE_REPLACE_ALL))
    {
      sfbval = UPDATE_SEQUENCE_ONLY;
      udp->rmcval = UPDATE_REPLACE;
      udp->useGUI = FALSE;
    }
    else if (udp->salp == NULL
             && (udp->no_aln_choice == UPDATE_SKIP_THIS
              || udp->no_aln_choice == UPDATE_SKIP_ALL))
    {
      do_update = FALSE;
    }
    else
    {
      sfbval = GetValue (udp->sfb);
      udp->rmcval = GetValue (udp->rmc);
      prior_rmcval = udp->rmcval;
    }
  }
  else
  {
    udp->rmcval = GetValue (udp->rmc);
    prior_rmcval = udp->rmcval;
  }

  if (udp->log_fp == NULL) 
  {
    OpenSequenceUpdateLog (udp);
    if (udp->log_fp == NULL) return;
    log_is_local = TRUE;
  }
  if (do_update && udp->do_update)
  {
    udp->affected_variation_features = ValNodeFree (udp->affected_variation_features);
    VisitFeaturesOnBsp (udp->oldbsp, udp, VariationAlignmentCallback);
    variation_action = VariationFeatureChangesOk (udp->affected_variation_features, udp->salp, TRUE);
    if (1 != variation_action)
    {
      do_update = FALSE;
    }
    
    
  }
    
  if (do_update)
  {
    if (udp->do_update)
    {
      if (udp->update_proteins != NULL 
          && Enabled (udp->update_proteins)
          && GetStatus (udp->update_proteins))
      {
        update_proteins = TRUE;
      }

      AddVariationFeaturesToLog (udp->log_fp, &(udp->data_in_log), 
                                 udp->affected_variation_features, udp->salp);
      UpdateOneSequence (udp, sfbval, GetStatus (udp->add_cit_subs),
                         update_proteins);
    }
    else
    {
      ExtendOneSequence (udp);
    }

    Remove (udp->form);
  }
  else
  {
    sip = SeqIdFindBest (udp->oldbsp->id, 0);
    if (sip != NULL)
    {
      sip_next = sip->next;
      SeqIdWrite (sip, acc_str, PRINTID_REPORT, sizeof (acc_str));
      fprintf (udp->log_fp, "Skipped update for %s\n", acc_str);
	    udp->data_in_log = TRUE;
      sip->next = sip_next;
    }
  }
  udp->affected_variation_features = ValNodeFree (udp->affected_variation_features);
  
  if (log_is_local) 
  {
    CloseOutSequenceUpdateLog (udp);
  }
  
  if (2 == variation_action)
  {
    FileClose (udp->fp);
    CloseOutSequenceUpdateLog (udp);
  }
  
  udp->useGUI = old_useGUI;
  udp->rmcval = prior_rmcval;
  /* if we are updating a set from a SeqSub, we don't want to free the SeqSub yet */
  if (udp->seqsubsep != NULL)
  {
    udp->newbsp = NULL;
  }
  FreeUdpFields (udp);
}

/*------------------------------------------------------------------*/
/*                                                                  */
/* AcceptRMCAll () -- Breaks out of the GUI interface and updates   */
/*                    all remaining sequences in the file without   */
/*                    user intervention.                            */
/*                                                                  */
/*------------------------------------------------------------------*/

static void AcceptRMCOrExtendAll (ButtoN b)

{
  UpsDataPtr   udp;
  Int2         state;
  Int4         count;
  Char         msgStr[256];

  /* Get current data */

  udp = (UpsDataPtr) GetObjectExtra (b);
  if (udp == NULL)
    return;

  /* Process the current sequence */

  WatchCursor ();
  udp->useGUI = FALSE;

  if (udp->log_fp == NULL) {
    OpenSequenceUpdateLog (udp);
    if (udp->log_fp == NULL) return;
  }
  
  
  udp->rmcval = GetValue (udp->rmc);

  DoAcceptRMCOrExtendSet (udp);

  /* Loop through the file, processing all others */

  state = FASTA_READ_OK;
  count = 0;

  while (FASTA_READ_OK == state) {
    count++;
    WatchCursor ();
    state = UpdateNextBioseqInFastaSet (udp);
    if (udp->useGUI)
      return;
  }

  CloseOutSequenceUpdateLog (udp);
  ArrowCursor ();

  /* If there was an error, report it, otherwise */
  /* print a status message.                     */

  if (FASTA_READ_ERROR == state) {
    sprintf (msgStr, "Encountered error while updating.  Only %d sequences "
	     "were updated.", count);
    Message (MSG_OK, msgStr);
  }
  else {
    sprintf (msgStr, "Successfully processed %d sequences from file (not"
	     " counting any updated before hitting 'Accept All').", count);
    Message (MSG_OK, msgStr);
  }
}


static void AcceptRMCOrExtendSet (ButtoN b)

{
  UpsDataPtr   udp;

  udp = (UpsDataPtr) GetObjectExtra (b);
  if (udp == NULL)
    return;
  
  if (udp->log_fp == NULL) {
    OpenSequenceUpdateLog (udp);
    if (udp->log_fp == NULL) return;
  }

  udp->rmcval = GetValue (udp->rmc);
  
  DoAcceptRMCOrExtendSet (udp);

  UpdateNextBioseqInFastaSet (udp);
}

static void SetProteinOptionsEnable (UpsDataPtr udp)
{
  if (udp == NULL || udp->update_proteins == NULL) return;
  
  if (Enabled (udp->update_proteins) && GetStatus (udp->update_proteins))
  {
    Enable (udp->truncate_proteins_btn);
    Enable (udp->extend_proteins3_btn);
    Enable (udp->extend_proteins5_btn);
    Enable (udp->correct_cds_genes_btn);
  }
  else
  {
    Disable (udp->truncate_proteins_btn);
    Disable (udp->extend_proteins3_btn);
    Disable (udp->extend_proteins5_btn);
    Disable (udp->correct_cds_genes_btn);
  }
}

static void SetStatusUpdateAcceptBtns (UpsDataPtr udp, Boolean status)
{
  if (udp == NULL) return;
  
  if (status)
  {
    SafeEnable (udp->accept);
    SafeEnable (udp->acceptAll);
  }
  else
  {
    SafeDisable (udp->accept);
    SafeDisable (udp->acceptAll);
  }    
}

static void UpdateAccept (GrouP g)

{
  Int2        rmcval, sfbval;
  UpsDataPtr  udp;
  Int2        nobmval = UPDATE_FEAT_DUP_NOT_SET;

  udp = (UpsDataPtr) GetObjectExtra (g);
  if (udp == NULL) return;

  rmcval = GetValue (udp->rmc);
  sfbval = GetValue (udp->sfb);
  
  if (rmcval == UPDATE_REPLACE)
  {
    sfbval = GetValue (udp->sfb);
    if (sfbval == UPDATE_SEQUENCE_ONLY || sfbval == UPDATE_SEQUENCE_AND_FEATURES)
    {
      Enable (udp->update_quality_scores_btn);
    }
    else
    {
      Disable (udp->update_quality_scores_btn);
    }
  }
  else
  {
    Disable (udp->update_quality_scores_btn);
  }
  
  if (sfbval == UPDATE_FEATURES_ONLY || sfbval == UPDATE_SEQUENCE_AND_FEATURES)
  {
    if (!GetStatus (udp->replace_all))
    {
      nobmval = GetValue (udp->nobm);
    }
  }
  
  /* set enables for protein updates */
  if (sfbval == UPDATE_SEQUENCE_ONLY
     || (sfbval == UPDATE_SEQUENCE_AND_FEATURES
        && (nobmval == UPDATE_FEAT_DUP_USE_OLD || nobmval == UPDATE_FEAT_DUP_USE_BOTH)))
  {
    SafeEnable (udp->update_proteins);
  }
  else
  {
    SafeDisable (udp->update_proteins);
  }
  SetProteinOptionsEnable (udp);    

  if (rmcval == UPDATE_CHOICE_NOT_SET) {
    SetStatusUpdateAcceptBtns (udp, FALSE);
    return;
  }
  if (! udp->do_update) {
    SetStatusUpdateAcceptBtns (udp, TRUE);
    return;
  }
    
  if (sfbval == UPDATE_CHOICE_NOT_SET || sfbval == UPDATE_SEQUENCE_ONLY || udp->diffOrgs) {
    SafeDisable (udp->keepProteinIDs);
    if (sfbval == UPDATE_CHOICE_NOT_SET || sfbval == UPDATE_SEQUENCE_ONLY) {
      SafeDisable (udp->nobm);
      SafeDisable (udp->replace_all);
    } else {
      SafeEnable (udp->replace_all);
      if (GetStatus (udp->replace_all))
      {
        SafeDisable (udp->nobm);     	
      }
      else
      {
      	SafeEnable (udp->nobm);
      }
    }
  } else if (sfbval == UPDATE_FEATURES_ONLY || sfbval == UPDATE_SEQUENCE_AND_FEATURES) {
    SafeEnable (udp->keepProteinIDs);
    SafeEnable (udp->replace_all);
    if (GetStatus (udp->replace_all))
    {
      SafeDisable (udp->nobm);      	
    }
    else
    {
      SafeEnable (udp->nobm);
    }
  }
  if (sfbval == UPDATE_CHOICE_NOT_SET) {
    SetStatusUpdateAcceptBtns (udp, FALSE);
    return;
  }
  if (sfbval == UPDATE_FEATURES_ONLY || sfbval == UPDATE_SEQUENCE_AND_FEATURES) {
    if (!GetStatus (udp->replace_all) && nobmval == UPDATE_FEAT_DUP_NOT_SET)
    {
      SetStatusUpdateAcceptBtns (udp, FALSE);
      return;
    }
  }

  SetStatusUpdateAcceptBtns (udp, TRUE);
}

static void UpdateButtons (GrouP g)
{
  UpsDataPtr        udp;
  Int2              rmcval;
  Uint2             entityID;
  SeqFeatPtr        sfp;
  SeqMgrFeatContext fcontext;

  udp = (UpsDataPtr) GetObjectExtra (g);
  if (udp == NULL) return;

  rmcval = GetValue (udp->rmc);

  if (udp->new5 <= udp->old5 && udp->new3 <= udp->old3) 
  {
    if (rmcval == UPDATE_PATCH)
    {
      /* If patch sequence matches, must be feature propagation only */

      if (StringNICmp (udp->seq1 + udp->old5 - udp->new5,
                       udp->seq2,
                       StringLen (udp->seq2)) == 0) {
        SetValue (udp->sfb, UPDATE_FEATURES_ONLY);
        Disable (udp->sfb);
      }
    }
    else if (rmcval == UPDATE_REPLACE)
    {
      /* If no features, must be sequence update only */

      entityID = ObjMgrGetEntityIDForPointer (udp->newbsp);
      if (! SeqMgrFeaturesAreIndexed (entityID))
        SeqMgrIndexFeatures (entityID, NULL);

      sfp = SeqMgrGetNextFeature (udp->newbsp, NULL, 0, 0, &fcontext);
      if (sfp == NULL) 
      {
        SetValue (udp->sfb, UPDATE_SEQUENCE_ONLY);
        Disable (udp->sfb);
        Disable (udp->replace_all);
        Disable (udp->nobm);
      }
      else if (!indexerVersion &&
               (udp->newbsp->repr != Seq_repr_raw || udp->oldbsp->repr != Seq_repr_raw))
      {
        SetValue (udp->sfb, UPDATE_FEATURES_ONLY);
        Disable (udp->sfb);
      }
      else
      {
        Enable (udp->sfb);
      }
    }
  }
  UpdateAccept (g);
}

static void DrawAlignBlock (
  SegmenT pict,
  Int4 top,
  Int4 bottom,
  Int4 labelpt,
  Int2 labelaln,
  Int4 len5,
  Int4 lena,
  Int4 len3,
  Int4 aln_length
)

{
  Char  str [96];

  if (len5 > 0) {
    AddRectangle (pict, -len5, top, 0, bottom, NO_ARROW, FALSE, 0);
  }
  sprintf (str, "%ld", (long) len5);
  AddLabel (pict, -len5, (top + bottom) / 2, str, SMALL_TEXT, 5, MIDDLE_LEFT, 0);

  if (len3 > 0) {
    AddRectangle (pict, aln_length, top, aln_length + len3, bottom, NO_ARROW, FALSE, 0);
  }
  sprintf (str, "%ld", (long) len3);
  AddLabel (pict, aln_length + len3, (top + bottom) / 2, str, SMALL_TEXT, 5, MIDDLE_RIGHT, 0);

  AddRectangle (pict, 0, top, aln_length, bottom, NO_ARROW, TRUE, 0);
  sprintf (str, "%ld", (long) lena);
  AddLabel (pict, aln_length / 2, labelpt, str, SMALL_TEXT, 5, labelaln, 0);
}

static SegmenT MakeAlignPicture (
  UpsDataPtr udp,
  CharPtr strid1,
  CharPtr strid2,
  SeqAlignPtr sap
)

{
  SegmenT  pict;
  Char     str [96];
  Int4     top, bottom;

  pict = CreatePicture ();
  if (sap == NULL) return pict;

  top = 0;
  bottom = top - 10;

  DrawAlignBlock (pict, top, bottom, bottom, LOWER_CENTER, udp->old5, udp->olda, udp->old3, udp->aln_length);

  /*
  AddLabel (pict, (udp->stopmax - udp->startmax) / 2, bottom - 20, strid1, SMALL_TEXT, 5, LOWER_CENTER, 0);
  */


  sprintf (str, "%ld", (long) udp->aln_length);
  AddLabel (pict, udp->aln_length / 2, 10, str, SMALL_TEXT, 5, MIDDLE_CENTER, 0);


  top = 30;
  bottom = top - 10;

  DrawAlignBlock (pict, top, bottom, top, UPPER_CENTER, udp->new5, udp->newa, udp->new3, udp->aln_length);

  /*
  AddLabel (pict, (udp->stopmax - udp->startmax) / 2, top + 20, strid2, SMALL_TEXT, 5, UPPER_CENTER, 0);
  */

  return pict;
}

static void DrawAlignDiffs (
  UpsDataPtr udp,
  SegmenT pict,
  Int4 top,
  Int4 bottom,
  SeqAlignPtr sap
)

{
  AlnMsg2Ptr  amp1, amp2;
  SegmenT    seg;
  Int4       len1, len2, i;
  Int4       seg_i, seg_n, seg_start, seg_stop;

  if (udp->seq1 == NULL || udp->seq2 == NULL) return;
  len1 = StringLen (udp->seq1);
  len2 = StringLen (udp->seq2);

  seg = CreateSegment (pict, 0, 0);
  AddAttribute (seg, COLOR_ATT, RED_COLOR, 0, 0, 0, 0);

  seg_n = AlnMgr2GetNumSegs(sap);
  for (seg_i = 1; seg_i<=seg_n; seg_i++) {
    AlnMgr2GetNthSegmentRange(sap, seg_i, &seg_start, &seg_stop);

    amp1 = AlnMsgNew2 ();
    amp2 = AlnMsgNew2 ();
    if (amp1 == NULL || amp2 == NULL) return;

    amp1->from_aln = seg_start;
    amp1->to_aln = seg_stop;
    amp1->row_num = 1;

    amp2->from_aln = seg_start;
    amp2->to_aln = seg_stop;
    amp2->row_num = 2;

    AlnMgr2GetNextAlnBit (sap, amp1);
    AlnMgr2GetNextAlnBit (sap, amp2);

    if (amp1->to_row - amp1->from_row == amp2->to_row - amp2->from_row &&
        amp1->type == AM_SEQ && amp2->type == AM_SEQ) {
      for (i=0; i<seg_stop-seg_start+1; i++) {
        if (udp->seq1[amp1->from_row+i] != udp->seq2[amp2->from_row+i]) {

          /* record for accurate scrolling to text view */
          ValNodeAddInt (&(udp->mismatches), 0, i);

          AddLine (seg, seg_start+i, top, seg_start+i, bottom, FALSE, 0);
        }
      }
    }

    AlnMsgFree2 (amp1);
    AlnMsgFree2 (amp2);
  }
}

static void DrawAlignBits (
  UpsDataPtr udp,
  SegmenT pict,
  Int4 top,
  Int4 bottom,
  Int4 row,
  Int4 pos1,
  Int4 pos2,
  SeqAlignPtr sap
)

{
  AlnMsg2Ptr  amp;
  Int4       len, start, stop, from, to;
  Char       str [96];
  Boolean    wasgap;

  amp = AlnMsgNew2 ();
  if (amp == NULL) return;

  amp->from_aln = 0;
  amp->to_aln = -1;
  amp->row_num = row;

  start = 0;
  stop = 0;
  from = 0;
  to = 0;
  wasgap = FALSE;

  while (AlnMgr2GetNextAlnBit (sap, amp)) {
    len = amp->to_row - amp->from_row + 1;
    stop = start + len;
    if (amp->type == AM_GAP) {
      if (wasgap) {
        to = stop;
      } else {
        AddRectangle (pict, from, top, to, bottom, NO_ARROW, FALSE, 0);
        wasgap = TRUE;
        from = start;
        to = stop;
      }
    } else {
      if (wasgap) {

        /* record for accurate scrolling to text view */
        ValNodeAddInt (&(udp->indels), 0, from);

        AddLine (pict, from, (top + bottom) / 2, to, (top + bottom) / 2, FALSE, 0);
        wasgap = FALSE;
        from = start;
        to = stop;
      } else {
        to = stop;
      }
    }
    start += len;
  }

  if (to > from) {
    if (wasgap) {

      /* record for accurate scrolling to text view */
      ValNodeAddInt (&(udp->indels), 0, from);

      AddLine (pict, from, (top + bottom) / 2, to, (top + bottom) / 2, FALSE, 0);
    } else {
      AddRectangle (pict, from, top, to, bottom, NO_ARROW, FALSE, 0);
    }
  }

  AlnMsgFree2 (amp);

  sprintf (str, "%ld", (long) pos1);
  AddLabel (pict, 0, (top + bottom) / 2, str, SMALL_TEXT, 5, MIDDLE_LEFT, 0);

  sprintf (str, "%ld", (long) pos2);
  AddLabel (pict, to, (top + bottom) / 2, str, SMALL_TEXT, 5, MIDDLE_RIGHT, 0);
}

static SegmenT MakeAlignDetails (
  UpsDataPtr udp,
  CharPtr strid1,
  CharPtr strid2,
  SeqAlignPtr sap
)

{
  Int4     aln_length;
  SegmenT  pict;
  Int4     top, bottom;

  pict = CreatePicture ();
  if (sap == NULL) return pict;

  aln_length = udp->aln_length;

  top = 0;
  bottom = top - 10;

  DrawAlignBits (udp, pict, top, bottom, 1, udp->old5 + 1, udp->old5 + udp->olda, sap);

  /*
  AddLabel (pict, aln_length / 2, bottom, strid1, SMALL_TEXT, 5, LOWER_CENTER, 0);
  */

  top = 30;
  bottom = top - 10;

  if (udp->revcomp) {
    DrawAlignBits (udp, pict, top, bottom, 2, udp->new3 + udp->newa, udp->new3 + 1, sap);
  } else {
    DrawAlignBits (udp, pict, top, bottom, 2, udp->new5 + 1, udp->new5 + udp->newa, sap);
  }

  /*
  AddLabel (pict, aln_length / 2, top, strid2, SMALL_TEXT, 5, UPPER_CENTER, 0);
  */

  top = 15;
  bottom = top - 10;

  DrawAlignDiffs (udp, pict, top, bottom, sap);

  return pict;
}

static CharPtr MakeAlignSequence (
  UpsDataPtr udp,
  SeqAlignPtr sap,
  Int4 row,
  CharPtr seq
)

{
  CharPtr    aln;
  AlnMsg2Ptr  amp;
  Int4       aln_length, len, lens, start, stop, from, to, i, j;

  if (udp == NULL || sap == NULL || seq == NULL || udp->aln_length < 1) return NULL;
  lens = StringLen (seq);

  aln = (CharPtr) MemNew (sizeof (Char) * (udp->aln_length + 2));
  if (aln == NULL) return NULL;
  aln_length = udp->aln_length;
  MemSet ((Pointer) aln, '-', aln_length);

  amp = AlnMsgNew2 ();
  if (amp == NULL) return aln;

  amp->from_aln = 0;
  amp->to_aln = -1;
  amp->row_num = row;

  start = 0;
  stop = 0;
  from = 0;
  to = 0;

  while (AlnMgr2GetNextAlnBit (sap, amp)) {
    len = amp->to_row - amp->from_row + 1;
    stop = start + len;

    if (amp->type == AM_SEQ) {
      for (i = start, j = amp->from_row; i < stop && j < lens; i++, j++) {
        aln [i] = seq [j];
      }
    }
    start += len;
  }

  AlnMsgFree2 (amp);

  return aln;
}

static void PrintTAln (ButtoN b)

{
  AsnIoPtr    aip;
  Char        path [PATH_MAX];
  UpsDataPtr  udp;

  udp = (UpsDataPtr) GetObjectExtra (b);
  if (udp == NULL) return;
  TmpNam (path);
  aip = AsnIoOpen (path, "w");
  if (aip != NULL) {
    SeqAlignAsnWrite (udp->salp, aip, NULL);
    AsnIoClose (aip);
    LaunchGeneralTextViewer (path, "Update Sequence Alignment");
  }
  FileRemove (path);
}

static void PrintGAln (ButtoN b)

{
  UpsDataPtr  udp;

  udp = (UpsDataPtr) GetObjectExtra (b);
  if (udp == NULL) return;
  PrintViewer (udp->overview);
  PrintViewer (udp->details);
}

static void CalculateOverhangs (
  UpsDataPtr udp
)

{
  Int4         aln_length;
  Uint2        entityID;
  SeqAlignPtr  sap;
  SeqEntryPtr  sep;
  Int4         stopold, startold, lenold, stopnew, startnew, lennew;

  if (udp == NULL) return;
  sap = udp->salp;
  if (sap == NULL) return;
  aln_length = AlnMgr2GetAlnLength (sap, FALSE);
  AlnMgr2GetNthSeqRangeInSA (sap, 1, &startold, &stopold);
  AlnMgr2GetNthSeqRangeInSA (sap, 2, &startnew, &stopnew);
  lenold = udp->oldbsp->length;
  lennew = udp->newbsp->length;

  udp->old5 = startold;
  udp->old3 = lenold - stopold - 1;
  udp->olda = stopold - startold + 1;

  udp->new5 = startnew;
  udp->new3 = lennew - stopnew - 1;
  udp->newa = stopnew - startnew + 1;

  udp->aln_length = aln_length;
  udp->startmax = MAX (startold, startnew);
  udp->stopmax = MAX (aln_length + lenold - stopold, aln_length + lennew - stopnew);

  udp->strandold = AlnMgr2GetNthStrand (sap, 1);
  udp->strandnew = AlnMgr2GetNthStrand (sap, 2);

  entityID = ObjMgrGetEntityIDForPointer (udp->oldbsp);
  sep = GetTopSeqEntryForEntityID (entityID);
  SeqEntrySetScope (sep);
  udp->seq1 = GetSequenceByBsp (udp->oldbsp);
  SeqEntrySetScope (NULL);

  entityID = ObjMgrGetEntityIDForPointer (udp->oldbsp);
  sep = GetTopSeqEntryForEntityID (entityID);
  SeqEntrySetScope (sep);
  udp->seq2 = GetSequenceByBsp (udp->newbsp);
  SeqEntrySetScope (NULL);

  udp->aln1 = MakeAlignSequence (udp, sap, 1, udp->seq1);
  udp->aln2 = MakeAlignSequence (udp, sap, 2, udp->seq2);

  udp->log10_aln_length = 1;
  while (aln_length >= 10) {
    aln_length /= 10;
    (udp->log10_aln_length)++;
  }
}

static Int4 CalculateBestScale (
  UpsDataPtr udp,
  VieweR vwr,
  SegmenT pict
)

{
  BoxInfo  box;
  Int2     i;
  Int4     max, worldwid, portwid;
  RecT     r;
  Int4     scaleX, oldscaleX;
  Int4     wid;

  ObjectRect (vwr, &r);
  InsetRect (&r, 4, 4);
  wid = (Int4) (r.right - r.left + 1);

  SegmentBox (pict, &box);
  oldscaleX = (box.right - box.left + wid - 1) / wid;
  RecalculateSegment (pict, oldscaleX, 1);
  SegmentBox (pict, &box);
  portwid = wid * oldscaleX;
  worldwid = box.right - box.left + 20 * oldscaleX + 1;
  max = MAX (worldwid, portwid);
  scaleX = (max + wid - 1) / wid;
  i = 0;
  while (i < 10 && (scaleX > oldscaleX || portwid < worldwid)) {
    oldscaleX = scaleX;
    RecalculateSegment (pict, oldscaleX, 1);
    SegmentBox (pict, &box);
    portwid = wid * oldscaleX;
    worldwid = box.right - box.left + 20 * oldscaleX + 1;
    max = MAX (worldwid, portwid);
    scaleX = (max + wid - 1) / wid;
    i++;
  }

  return scaleX;
}

static Uint1 leftTriFillSym [] = {
  0x0C, 0x3C, 0xFC, 0x3C, 0x0C, 0x00, 0x00, 0x00
};
static Uint1 rightTriFillSym [] = {
  0xC0, 0xF0, 0xFC, 0xF0, 0xC0, 0x00, 0x00, 0x00
};

static void LetDraw (
  PaneL pnl
)

{
  Char        ch1, ch2;
  Int2        i, k, q, left, top, bottom, arrowwidth;
  size_t      len;
  Int4        offset, j, pos, realpos;
  RecT        r, x;
  BaR         sb;
  Char        str [32];
  UpsDataPtr  udp;

  udp = (UpsDataPtr) GetObjectExtra (pnl);
  if (udp == NULL) return;

  ObjectRect (pnl, &r);
  InsetRect (&r, 4, 4);

  sb = GetSlateHScrollBar ((SlatE) pnl);
  offset = GetBarValue (sb);

  SelectFont (SetSmallFont ());

  /* draw top (new) letters */

  if (udp->aln2 != NULL)
  {
    MoveTo (r.left, r.top + 8 + 3 * udp->lineheight);
    for (i = 0, j = offset; i < udp->maxchars && j < udp->aln_length; i++, j++) {
      PaintChar (udp->aln2 [j]);
    }
  }

  /* draw bottom (old) letters */

  if (udp->aln1 != NULL) 
  {
    MoveTo (r.left, r.top + 8 + 5 * udp->lineheight);
    for (i = 0, j = offset; i < udp->maxchars && j < udp->aln_length; i++, j++) {
      PaintChar (udp->aln1 [j]);
    }
  }

  /* draw recombination arrows */

  arrowwidth = MIN (6, udp->charwidth);
  if (udp->recomb1 >= offset && udp->recomb1 <= offset + udp->maxchars) {
    left = r.left + udp->charwidth * (udp->recomb1 - offset);
    LoadRect (&x, left, r.top, left + arrowwidth, r.top + 6);
    CopyBits (&x, leftTriFillSym);
  }

  if (udp->recomb2 >= offset && udp->recomb2 <= offset + udp->maxchars) {
    left = r.left + udp->charwidth * (udp->recomb2 - offset - 1);
    LoadRect (&x, left, r.top, left + arrowwidth, r.top + 6);
    CopyBits (&x, rightTriFillSym);
  }

  if (udp->aln1 == NULL || udp->aln2 == NULL) 
  {
  	return;
  }
  /* draw red mismatch lines */

  Red ();
  top = r.top + 8 + 4 * udp->lineheight - Ascent ();
  bottom = top + udp->lineheight - 2;

  for (i = 0, j = offset; i < udp->maxchars && j < udp->aln_length; i++, j++) {
    ch1 = udp->aln1 [j];
    ch2 = udp->aln2 [j];
    if (ch1 == ch2) {
    } else if (ch1 == '-' || ch2 == '-') {
    } else {
      left = r.left + i * udp->charwidth + udp->charwidth / 2 - 1;
      MoveTo (left, top);
      LineTo (left, bottom);
    }
  }
  Black ();

  /* draw top (new) tick marks and coordinates */

  bottom = r.top + 8 + 3 * udp->lineheight - Ascent () - 2;
  top = bottom - 5;
  i = 0;
  j = offset;
  pos = AlnMgr2MapSeqAlignToBioseq (udp->salp, j, 2);
  while (pos < 1 && i < udp->maxchars && j < udp->aln_length) {
    i++;
    j++;
    pos = AlnMgr2MapSeqAlignToBioseq (udp->salp, j, 2);
  }
  for (; i < udp->maxchars + udp->log10_aln_length && j < udp->aln_length; i++, j++) {
    ch1 = udp->aln2 [j];
    if (ch1 != '-') {
      if (udp->revcomp) {
        realpos = (udp->newbsp->length - pos - 1);
      } else {
        realpos = pos;
      }
      if (((realpos + 1) % 10) == 0) {
        left = r.left + i * udp->charwidth + udp->charwidth / 2 - 1;
        if (i < udp->maxchars) {
          MoveTo (left, top);
          LineTo (left, bottom);
        }
        sprintf (str, "%ld", (long) (realpos + 1));
        len = StringLen (str);
        if (len <= j + 1) {
          k = i - len + 1;
          q = 0;
          if (k < 0) {
            q -= k;
            k = 0;
          }
          if (q < len) {
            left = r.left + k * udp->charwidth;
            MoveTo (left, r.top + 8 + udp->lineheight);
            while (k < udp->maxchars && q < len) {
              PaintChar (str [q]);
              k++;
              q++;
            }
          }
        }
      } else if (((realpos + 1) % 5) == 0) {
        left = r.left + i * udp->charwidth + udp->charwidth / 2 - 1;
        if (i < udp->maxchars) {
          MoveTo (left, top + 3);
          LineTo (left, bottom);
        }
      }
      pos++;
    }
  }

  /* draw bottom (old) tick marks and coordinates */

  top = r.top + 8 + 6 * udp->lineheight - Ascent () + 2;
  bottom = top + 5;
  i = 0;
  j = offset;
  pos = AlnMgr2MapSeqAlignToBioseq (udp->salp, j, 1);
  while (pos < 1 && i < udp->maxchars && j < udp->aln_length) {
    i++;
    j++;
    pos = AlnMgr2MapSeqAlignToBioseq (udp->salp, j, 1);
  }
  for (; i < udp->maxchars + udp->log10_aln_length && j < udp->aln_length; i++, j++) {
    ch1 = udp->aln1 [j];
    if (ch1 != '-') {
      if (((pos + 1) % 10) == 0) {
        left = r.left + i * udp->charwidth + udp->charwidth / 2 - 1;
        if (i < udp->maxchars) {
          MoveTo (left, top);
          LineTo (left, bottom);
        }
        sprintf (str, "%ld", (long) (pos + 1));
        len = StringLen (str);
        if (len <= j + 1) {
          k = i - len + 1;
          q = 0;
          if (k < 0) {
            q -= k;
            k = 0;
          }
          if (q < len) {
            left = r.left + k * udp->charwidth;
            MoveTo (left, r.top + 8 + 7 * udp->lineheight);
            while (k < udp->maxchars && q < len) {
              PaintChar (str [q]);
              k++;
              q++;
            }
          }
        }
      } else if (((pos + 1) % 5) == 0) {
        left = r.left + i * udp->charwidth + udp->charwidth / 2 - 1;
        if (i < udp->maxchars) {
          MoveTo (left, top);
          LineTo (left, bottom - 3);
        }
      }
      pos++;
    }
  }
  SelectFont (systemFont);
}

static void LetScrl (
  BaR sb,
  SlatE slt,
  Int4 newval,
  Int4 oldval
)

{
  RecT        r;
  UpsDataPtr  udp;

  udp = (UpsDataPtr) GetObjectExtra (slt);
  if (udp == NULL) return;

  ObjectRect (udp->letters, &r);
  InsetRect (&r, 4, 4);
  Select (udp->letters);
  if (ABS (oldval - newval) < udp->maxchars) {
    ScrollRect (&r, (oldval - newval) * udp->charwidth, 0);
  } else {
    InsetRect (&r, -2, -2);
    InvalRect (&r);
  }
  Update ();
}

static void DtlClck (
  VieweR vwr,
  SegmenT pict,
  PoinT pt
)

{
  Int4        goHere;
  Int4        offset;
  Int4        maxover2;
  PntInfo     pnt;
  BaR         sb;
  UpsDataPtr  udp;
  ValNodePtr  vnp;

  udp = (UpsDataPtr) GetViewerData (vwr);
  if (udp == NULL) return;

  sb = GetSlateHScrollBar ((SlatE) udp->letters);

  MapViewerToWorld (vwr, pt, &pnt);
  maxover2 = udp->maxchars / 2;
  if (pnt.x <= 0) {
    pnt.x = 0;
  } else if (pnt.x >= udp->aln_length) {
    pnt.x = udp->aln_length  - udp->maxchars;
  } else if (pnt.x >= maxover2) {

    offset = GetBarValue (sb);

    /* look for clicks within 5 pixels of an indel start or a mismatch */

    goHere = -1;
    for (vnp = udp->indels; vnp != NULL && goHere < 0; vnp = vnp->next) {
      if (ABS (pnt.x - vnp->data.intvalue) < udp->scaleX * 5) {
        goHere = vnp->data.intvalue;
      }
    }
    for (vnp = udp->mismatches; vnp != NULL && goHere < 0; vnp = vnp->next) {
      if (ABS (pnt.x - vnp->data.intvalue) < udp->scaleX * 5) {
        goHere = vnp->data.intvalue;
      }
    }

    if (goHere >= 0) {
      pnt.x = goHere;
    } else {
      /* if already visible, no need to scroll */
      if (pnt.x - maxover2 > offset && pnt.x - maxover2 < offset + maxover2 - 5) return;
      if (pnt.x - maxover2 < offset && pnt.x - maxover2 > offset - maxover2 + 5) return;
    }

    /* go left 1/2 screen so desired point is in the middle */

    pnt.x -= maxover2;
  }

  ResetClip ();
  SetBarValue (sb, pnt.x);
  Update ();
}

static void FrameVwr (
  VieweR vwr,
  SegmenT pict
)

{
  RecT  r;

  ResetClip ();
  ObjectRect (vwr, &r);
  FrameRect (&r);
}

static int LIBCALLBACK SortVnpByInt (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;

  if (vnp1->data.intvalue > vnp2->data.intvalue) {
    return 1;
  } else if (vnp1->data.intvalue < vnp2->data.intvalue) {
    return -1;
  }

  return 0;
}


static void UpdateSequenceFormMessage (ForM f, Int2 mssg)

{
  BaseFormPtr        bfp;
  StdEditorProcsPtr  sepp;

  bfp = (BaseFormPtr) GetObjectExtra (f);
  if (bfp == NULL) return;
  switch (mssg) {
    case VIB_MSG_CLOSE :
      Remove (f);
      break;
    case VIB_MSG_QUIT :
      QuitProc ();
      break;
    default :
      sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
      if (sepp != NULL && sepp->handleMessages != NULL) {
        sepp->handleMessages (f, mssg);
      }
      break;
  }
}

static void FreeUdpFields (UpsDataPtr udp)
{
  Uint2       entityID;

  udp->ovpict     = DeletePicture (udp->ovpict);
  udp->dtpict     = DeletePicture (udp->dtpict);
  udp->salp       = SeqAlignFree (udp->salp);
  entityID        = ObjMgrGetEntityIDForPointer (udp->newbsp);
  udp->newbsp     = ObjMgrFreeByEntityID (entityID);
  udp->seq1       = MemFree (udp->seq1);
  udp->seq2       = MemFree (udp->seq2);
  udp->aln1       = MemFree (udp->aln1);
  udp->aln1       = NULL;
  udp->aln2       = MemFree (udp->aln2);
  udp->aln2       = NULL;
  udp->indels     = ValNodeFree (udp->indels);
  udp->mismatches = ValNodeFree (udp->mismatches);
  udp->transl_except_list = ValNodeFree (udp->transl_except_list);
  udp->affected_variation_features = ValNodeFree (udp->affected_variation_features);
}

static void CleanupUpdateSequenceForm (GraphiC g, VoidPtr data)

{
  UpsDataPtr  udp;

  udp = (UpsDataPtr) data;
  if (udp != NULL)
    FreeUdpFields (udp);
  StdCleanupFormProc (g, data);
}

static CharPtr txt1 =
  "Sequence Relationship displays sequence lengths";

static CharPtr txt2 =
  "Alignment Details displays sequence positions";

static CharPtr txt3 =
  "Click above to scroll Alignment Text position";

/*------------------------------------------------------------------*/
/*                                                                  */
/* DetermineButtonState () -- Enable/disable buttons based on the   */
/*                            nature of the alignment.              */
/*                                                                  */
/*------------------------------------------------------------------*/

static void DetermineButtonState (UpsDataPtr   udp,
				  ButtoN PNTR  replaceButtonPtr,
				  ButtoN PNTR  extend5ButtonPtr,
				  ButtoN PNTR  extend3ButtonPtr,
				  ButtoN PNTR  patchButtonPtr)
{
  BioSourcePtr       biop1;
  BioSourcePtr       biop2;
  SeqMgrDescContext  dcontext;
  Uint2              entityID;
  SeqMgrFeatContext  fcontext;
  OrgRefPtr          orp1;
  OrgRefPtr          orp2;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;

  /* If no alignment then disable the patch button */

  if (udp->salp == NULL) {
    if (udp->do_update && patchButtonPtr != NULL && *patchButtonPtr != NULL) {
      Disable (*patchButtonPtr);
    }
    if (!udp->do_update)
    {
      Disable (*replaceButtonPtr);
    }
  }

  /* Extend 5' */

  else if (udp->new5 > udp->old5 && udp->new3 < udp->old3) {
    SetValue (udp->rmc, UPDATE_EXTEND5);
    Disable (*extend3ButtonPtr);
    udp->recomb2 = udp->aln_length;
    if (! udp->do_update) {
      Disable (*replaceButtonPtr);
    }
  }

  /* Extend 3' */

  else if (udp->new5 < udp->old5 && udp->new3 > udp->old3) {
    SetValue (udp->rmc, UPDATE_EXTEND3);
    Disable (*extend5ButtonPtr);
    udp->recomb1 = 0;
    if (! udp->do_update) {
      Disable (*replaceButtonPtr);
    }
  }

  /* Replace */

  else {
    SetValue (udp->rmc, UPDATE_REPLACE);
    Disable (*extend5ButtonPtr);
    Disable (*extend3ButtonPtr);
    udp->recomb1 = 0;
    udp->recomb2 = udp->aln_length;
  }
  
  switch (udp->rmcval)
  {
    case UPDATE_REPLACE:
      if (Enabled (*replaceButtonPtr))
      {
        SetValue (udp->rmc, UPDATE_REPLACE);
      }
      break;
    case UPDATE_EXTEND5:
      if (Enabled (*extend5ButtonPtr))
      {
        SetValue (udp->rmc, UPDATE_EXTEND5);
      }
      break;
    case UPDATE_EXTEND3:
      if (Enabled (*extend3ButtonPtr))
      {
        SetValue (udp->rmc, UPDATE_EXTEND3);
      }
      break;
    case UPDATE_PATCH:
      if (Enabled (*patchButtonPtr))
      {
        SetValue (udp->rmc, UPDATE_PATCH);
      }
      break;
  }

  /* If no features, must be sequence update only */

  entityID = ObjMgrGetEntityIDForPointer (udp->newbsp);
  if (! SeqMgrFeaturesAreIndexed (entityID))
    SeqMgrIndexFeatures (entityID, NULL);

  sfp = SeqMgrGetNextFeature (udp->newbsp, NULL, 0, 0, &fcontext);
  if (sfp == NULL) {
    SetValue (udp->sfb, UPDATE_SEQUENCE_ONLY);
    Disable (udp->sfb);
    Disable (udp->replace_all);
    Disable (udp->nobm);
  }

  /* If different organisms, must be feature propagation only */

  orp1 = NULL;
  orp2 = NULL;
  sdp = SeqMgrGetNextDescriptor (udp->oldbsp, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL) {
    biop1 = (BioSourcePtr) sdp->data.ptrvalue;
    if (biop1 != NULL) {
      orp1 = biop1->org;
    }
  }
  sdp = SeqMgrGetNextDescriptor (udp->newbsp, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL) {
    biop2 = (BioSourcePtr) sdp->data.ptrvalue;
    if (biop2 != NULL) {
      orp2 = biop2->org;
    }
  }
  if (orp1 != NULL && orp2 != NULL) {
    if (StringICmp (orp1->taxname, orp2->taxname) != 0) {
      if (sfp != NULL) {
        SetValue (udp->sfb, UPDATE_FEATURES_ONLY);
        Disable (udp->sfb);
        udp->diffOrgs = TRUE;
	if (FALSE == udp->isSet)
	  Message (MSG_OK, "Organisms are different, so features will"
		   " be propagated, but sequence will not be changed");
      } else {
        /* no features, cannot do anything */
        SetValue (udp->sfb, UPDATE_CHOICE_NOT_SET);
        Disable (udp->sfb);
        Disable (udp->replace_all);
        Disable (udp->nobm);
	if (FALSE == udp->isSet)
	  Message (MSG_OK, "Organisms are different, no features"
		   " to propagate, so nothing to do");
      }
      Disable (udp->rmc);
    }
  }

  /* If either sequence is not raw and not indexer version, only allow feature propagation */
  if (!indexerVersion &&
      (udp->oldbsp->repr != Seq_repr_raw || udp->newbsp->repr != Seq_repr_raw)) {
    SetValue (udp->sfb, UPDATE_FEATURES_ONLY);
    Disable (udp->sfb);
  }

  /* Disable accept button unless rmc and sfb are both preset */

  UpdateAccept (udp->rmc);

}

static void ChangeUpdateReplaceAll (ButtoN b)
{
  UpsDataPtr udp;
  
  if (b == NULL) return;
  udp = (UpsDataPtr) GetObjectExtra (b);
  if (udp == NULL) return;
 
  if (GetStatus (b))
  {
  	Disable (udp->nobm);
  }
  else
  {
  	Enable (udp->nobm);
  }
  /* update accept button */
  UpdateAccept (udp->rmc);
}

static void CancelUpdate (ButtoN b)
{
  UpsDataPtr udp;
  
  if (b == NULL) return;
  udp = (UpsDataPtr) GetObjectExtra (b);
  if (udp == NULL) return;
  CloseOutSequenceUpdateLog (udp);
  StdCancelButtonProc (b);
}

static GrouP CreateUpdateOperationsGroup (GrouP parent, UpsDataPtr udp)
{
  GrouP  g;
  GrouP  gp1, gp2, gp3;
  ButtoN b1 = NULL, b2 = NULL, b3 = NULL, b4 = NULL;
  
  if (udp == NULL) return NULL;
  
  g = HiddenGroup (parent, -1, 0, NULL);
  SetGroupSpacing (g, 5, 5);
  
  gp1 = NormalGroup (g, 4, 0, "Alignment Relationship", programFont, NULL);
  udp->rmc = HiddenGroup (gp1, 4, 0, UpdateButtons);
  SetObjectExtra (udp->rmc, (Pointer) udp, NULL);
  SetGroupSpacing (udp->rmc, 10, 5);
  if (udp->do_update) {
    b1 = RadioButton (udp->rmc, "Replace");
  } else {
    b1= RadioButton (udp->rmc, "Extend Both Ends");
  }
  b2 = RadioButton (udp->rmc, "Extend 5'");
  b3 = RadioButton (udp->rmc, "Extend 3'");
  if (udp->do_update) {
    b4 = RadioButton (udp->rmc, "Patch");
  } else {
    b4 = NULL;
  }

  if (udp->do_update) {
    gp2 = NormalGroup (g, 4, 0, "Update Operation", programFont, NULL);
    udp->sfb = HiddenGroup (gp2, 3, 0, UpdateAccept);
    SetObjectExtra (udp->sfb, (Pointer) udp, NULL);
    SetGroupSpacing (udp->sfb, 10, 5);
    RadioButton (udp->sfb, "Sequence");
    RadioButton (udp->sfb, "Features");
    RadioButton (udp->sfb, "Sequence + Features");
  
    udp->keepProteinIDs = CheckBox (g, "Keep Protein IDs", NULL);
  
    gp3 = NormalGroup (g, 1, 0, "Feature Policy", programFont, NULL);
    udp->replace_all = CheckBox (gp3, "Replace All Features", ChangeUpdateReplaceAll);
    SetObjectExtra (udp->replace_all, (Pointer) udp, NULL);
    SetValue (udp->replace_all, FALSE);
    
    udp->nobm = NormalGroup (gp3, 5, 0, "Duplicate Features Only", programFont, UpdateAccept);
    SetObjectExtra (udp->nobm, (Pointer) udp, NULL);
    SetGroupSpacing (udp->nobm, 10, 5);
    RadioButton (udp->nobm, "New");
    RadioButton (udp->nobm, "Old");
    RadioButton (udp->nobm, "Both");
    RadioButton (udp->nobm, "Merge");
    RadioButton (udp->nobm, "Replace");
  
    AlignObjects (ALIGN_CENTER, (HANDLE) gp1, (HANDLE) gp2, 
                  (HANDLE) udp->keepProteinIDs,
                  (HANDLE) gp3,  NULL);
  }
  /* Enable/disable buttons based on the nature of the alignment */
    
  DetermineButtonState (udp, &b1, &b2, &b3, &b4);
  return g;
  
}


static void ChangeProteinUpdateStatus (ButtoN b)
{
  UpsDataPtr udp;
  
  udp = (UpsDataPtr) GetObjectExtra (b);
  if (udp == NULL) return;
  SetProteinOptionsEnable (udp);
}

static void SkipUpdate (UpsDataPtr udp)
{
  Char     acc_str [256];
  SeqIdPtr sip, sip_next;
  
  if (udp == NULL) return;

  OpenSequenceUpdateLog (udp);
  if (udp->log_fp != NULL && udp->oldbsp != NULL && udp->oldbsp->id != NULL)
  {
    sip = SeqIdFindBest (udp->oldbsp->id, 0);
    if (sip != NULL)
    {
      sip_next = sip->next;
      sip->next = NULL;
      SeqIdWrite (sip, acc_str, PRINTID_REPORT, sizeof (acc_str));
      fprintf (udp->log_fp, "Skipped update for %s\n", acc_str);
	    udp->data_in_log = TRUE;
	    sip->next = sip_next;
    }
  }
  
  /* if we are updating a set from a SeqSub, we don't want to free the SeqSub yet */
  if (udp->seqsubsep != NULL)
  {
    udp->newbsp = NULL;
  }
  FreeUdpFields (udp);
  UpdateNextBioseqInFastaSet (udp); 
}

static GrouP CreateExtraUpdateOptionsGroup (GrouP g, UpsDataPtr udp)
{
  GrouP y, protein_options;
  
  if (udp == NULL || g == NULL) return NULL;

  y = HiddenGroup (g, -1, 0, NULL);
  
  udp->add_cit_subs = CheckBox (y, "Add Cit-subs for Updated Sequences", NULL);
  udp->update_quality_scores_btn = NULL;
  if (! ISA_aa (udp->oldbsp->mol) && udp->do_update)
  {
    udp->update_quality_scores_btn = CheckBox (y, "Replace Quality Scores", NULL);
    SetStatus (udp->update_quality_scores_btn, TRUE);

    udp->update_proteins = CheckBox (y, "Update Proteins for Updated Sequences", ChangeProteinUpdateStatus);
    SetObjectExtra (udp->update_proteins, (Pointer) udp, NULL);
    SetStatus (udp->update_proteins, FALSE);
    protein_options = HiddenGroup (y, 1, 0, NULL);
    udp->truncate_proteins_btn = CheckBox (protein_options,
                                     "Truncate retranslated proteins at stops",
                                           NULL);
    SetStatus (udp->truncate_proteins_btn,
               udp->truncate_proteins);
    udp->extend_proteins3_btn = CheckBox (protein_options,
                                  "Extend retranslated proteins without stops",
                                         NULL);
    udp->extend_proteins5_btn = CheckBox (protein_options,
                                  "Extend retranslated proteins without starts",
                                         NULL);
    udp->correct_cds_genes_btn = CheckBox (protein_options, "Correct CDS genes", NULL);
  
    SetStatus (udp->extend_proteins3_btn, udp->extend_proteins3);
    SetStatus (udp->extend_proteins5_btn, udp->extend_proteins5);
    SetStatus (udp->correct_cds_genes_btn, udp->correct_cds_genes);
  }
  AlignObjects (ALIGN_CENTER, (HANDLE) udp->add_cit_subs, 
                (HANDLE) udp->update_quality_scores_btn,
                (HANDLE) udp->update_proteins,
                (HANDLE) protein_options,
                NULL);
  return y;
}

static void SkipUpdateBtn (ButtoN b)
{
  UpsDataPtr udp;
  
  udp = (UpsDataPtr) GetObjectExtra (b);
  if (udp == NULL) return;
  SafeHide (udp->form);
  Remove (udp->form);

  SkipUpdate (udp);  
}

/*------------------------------------------------------------------*/
/*                                                                  */
/* UpdateSequenceForm () -- Compares two sequences and displays a   */
/*                          window giving the user options on how   */
/*                          to update one from the other.           */
/*                                                                  */
/*------------------------------------------------------------------*/

static ForM UpdateSequenceForm (UpsDataPtr udp)
{
  ButtoN             b;
  GrouP              c, g, y, k, x, z = NULL;
  Uint2              hgt;
  GrouP              ppt0, ppt1, ppt2, ppt3;
  RecT               r;
  BaR                sb;
  Int4               scaleX;
  SeqIdPtr           sip;
  Char               strid1 [MAX_ID_LEN], strid2 [MAX_ID_LEN], txt0 [256];
  CharPtr            title;
  WindoW             w;
  GrouP              misc_options;
  Int4               prompt_width = 400;
  GrouP              left_panel;
  GrouP              right_panel;

  /* Check parameters */

  if ((udp->oldbsp == NULL) || (udp->newbsp == NULL))
    return NULL;

  /* Create window */

  if (udp->do_update) {
    title = "Update Sequence";
  } else {
    title = "Extend Sequence";
  }
  w = FixedWindow (-50, -33, -10, -10, title, NULL);

  if (w == NULL)
    return NULL;
  
  if (FALSE == udp->isSet)
    SetObjectExtra (w, (Pointer) udp, CleanupUpdateSequenceForm);
  udp->form = (ForM) w;
  
  /* Get string IDs for the Bioseqs */

  sip = SeqIdFindWorst (udp->oldbsp->id);
  SeqIdWrite (sip, strid1, PRINTID_REPORT, sizeof (strid1) - 1);
  sip = SeqIdFindWorst (udp->newbsp->id);
  SeqIdWrite (sip, strid2, PRINTID_REPORT, sizeof (strid2) - 1);
  if (StringNICmp (strid2, "SequinUpdateSequence", 20) == 0 &&
      udp->newbsp->id->next != NULL) {
    sip = SeqIdFindWorst (udp->newbsp->id->next);
    SeqIdWrite (sip, strid2, PRINTID_REPORT, sizeof (strid2) - 1);
  }

  /* FIll in some of the data structure */
  /* for passing to the callbacks.      */

  udp->formmessage = UpdateSequenceFormMessage;

#ifdef WIN_MAC
  udp->activate = UpdateSequenceFormActivated;
  SetActivate (w, UpdateSequenceFormActivate);
#endif

  udp->diffOrgs = FALSE;

  CalculateOverhangs (udp);

  /* Display the sequences */

  sprintf (txt0,
	   "New sequence: %s - Length: %ld\nOld Sequence: %s - Length: %ld",
	   strid2, (long) udp->newbsp->length, strid1,
	   (long) udp->oldbsp->length);
  ppt0 = MultiLinePrompt (w, txt0, prompt_width, programFont);
  
  x = HiddenGroup (w, 2, 0, NULL);
  left_panel = HiddenGroup (x, -1, 0, NULL);
  y = left_panel;
  
  ppt1 = MultiLinePrompt (y, txt1, prompt_width, programFont);
  udp->overview = CreateViewer (y, prompt_width + Nlm_vScrollBarWidth, 100,
				FALSE, FALSE);
  
  ppt2 = MultiLinePrompt (y, txt2, prompt_width, programFont);
  udp->details = CreateViewer (y, prompt_width + Nlm_vScrollBarWidth, 80,
			       FALSE, FALSE);
  
  ppt3 = MultiLinePrompt (y, txt3, prompt_width, programFont);
    
#ifdef WIN_MAC
  hgt = 90;
#else
  hgt = 110;
#endif
  udp->letters = AutonomousPanel4 (y, prompt_width + Nlm_vScrollBarWidth, hgt,
				   LetDraw, NULL, LetScrl, 0, NULL, NULL);
  SetObjectExtra (udp->letters, (Pointer) udp, NULL);
  
  if (indexerVersion && shftKey) {
    ButtoN  b;
    z = HiddenGroup (y, 2, 0, NULL);
    SetGroupSpacing (z, 10, 3);
    b = PushButton (z, "Print Graphic", PrintGAln);
    SetObjectExtra (b, (Pointer) udp, NULL);
    b = PushButton (z, "Display Alignment", PrintTAln);
    SetObjectExtra (b, (Pointer) udp, NULL);
  }
  
  udp->ovpict = NULL;
  udp->dtpict = NULL;

  AlignObjects (ALIGN_CENTER, (HANDLE) udp->overview,
                (HANDLE) udp->details, (HANDLE) udp->letters,
                (HANDLE) ppt1, (HANDLE) ppt2,
                (HANDLE) ppt3, NULL);

  right_panel = HiddenGroup (x, -1, 0, NULL);
  y = right_panel;

  k = HiddenGroup (y, -1, 0, NULL);
  g = CreateUpdateOperationsGroup (k, udp);
  misc_options = CreateExtraUpdateOptionsGroup (k, udp);  
  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) misc_options, NULL);
  
  c = HiddenGroup (w, 5, 0, NULL);
  if (udp->isSet)
  {
    udp->accept = DefaultButton (c, "Accept", AcceptRMCOrExtendSet);
    SetObjectExtra (udp->accept, (Pointer) udp, NULL);
    udp->acceptAll = DefaultButton (c, "Accept All", AcceptRMCOrExtendAll);
    SetObjectExtra (udp->acceptAll, (Pointer) udp, NULL);
    b = PushButton (c, "Skip", SkipUpdateBtn);
    SetObjectExtra (b, (Pointer) udp, NULL);  
  }
  else
  {
    udp->accept = DefaultButton (c, "Accept", AcceptRMCOrExtend);
    SetObjectExtra (udp->accept, (Pointer) udp, NULL);
  }
  
  b = PushButton (c, "Cancel", CancelUpdate);
  SetObjectExtra (b, (Pointer) udp, NULL);  
  UpdateButtons (udp->rmc);

  AlignObjects (ALIGN_CENTER,     
                (HANDLE) ppt0,
                (HANDLE) x,
                (HANDLE) c, (HANDLE) z, NULL);
  RealizeWindow (w);
  
  udp->ovpict = MakeAlignPicture (udp, strid1, strid2, udp->salp);
  scaleX = CalculateBestScale (udp, udp->overview, udp->ovpict);
  AttachPicture (udp->overview, udp->ovpict, 0, 0, UPPER_LEFT,
		 scaleX, 1, FrameVwr);
  
  udp->dtpict = MakeAlignDetails (udp, strid1, strid2, udp->salp);
  scaleX = CalculateBestScale (udp, udp->details, udp->dtpict);
  udp->scaleX = scaleX;
  AttachPicture (udp->details, udp->dtpict, 0, 0, UPPER_LEFT,
		 scaleX, 1, FrameVwr);
  SetViewerData (udp->details, (Pointer) udp, NULL);
  SetViewerProcs (udp->details, DtlClck, NULL, NULL, NULL);
  
  udp->indels = ValNodeSort (udp->indels, SortVnpByInt);
  udp->mismatches = ValNodeSort (udp->mismatches, SortVnpByInt);
  
  SelectFont (SetSmallFont ());
  ObjectRect (udp->letters, &r);
  InsetRect (&r, 4, 4);
  udp->lineheight = LineHeight ();
  udp->charwidth = MaxCharWidth ();
  udp->maxchars = (r.right-r.left-2+udp->charwidth - 1) / udp->charwidth;
  SelectFont (systemFont);
  
  sb = GetSlateHScrollBar ((SlatE) udp->letters);
  SetBarMax (sb, udp->aln_length - (Int4) udp->maxchars);
  CorrectBarPage (sb, (Int4) udp->maxchars - 1, (Int4) udp->maxchars - 1);
  
  udp->recomb1 = -1;
  udp->recomb2 = -1;
  
  return (ForM) w;
}


/*=====================================================================*/
/*                                                                     */
/* PrepareToUpdateSequences ()                                         */
/*                                                                     */
/*=====================================================================*/

static Boolean PrepareToUpdateSequences (UpsDataPtr udp)
{
  ForM         f;

  if ( ! PrepareUpdatePtr (udp)) 
  {
    CloseOutSequenceUpdateLog (udp);
    return FALSE;	
  }

  if (TRUE == udp->useGUI)
  {
    if (udp->salp == NULL && udp->do_update &&
        (udp->no_aln_choice == UPDATE_REPLACE_THIS
         || udp->no_aln_choice == UPDATE_SKIP_THIS
         || udp->no_aln_choice == UPDATE_REPLACE_ALL
         || udp->no_aln_choice == UPDATE_REPLACE_THIS)) 
    {
      CalculateOverhangs (udp);
      DoAcceptRMCOrExtendSet (udp);
      UpdateNextBioseqInFastaSet (udp);
    }
    else
    {
      f = UpdateSequenceForm (udp);
      if (f == NULL) 
      {
        CloseOutSequenceUpdateLog (udp);	
        return FALSE;
      }
      Show (f);
      Select (f);
      SendHelpScrollMessage (helpForm, "Edit Menu", "Update Sequence");
    }
  }
  else {
    CalculateOverhangs (udp);
    DoAcceptRMCOrExtendSet (udp);
  }
  return TRUE;
}

/*=====================================================================*/
/*                                                                     */
/* FindMatchingBioseq () -- Callback function for exploring Bioseqs.   */
/*                          Finds the bioseq that matches a given      */
/*                          string ID.                                 */
/*                                                                     */
/*=====================================================================*/

static Boolean LIBCALLBACK FindMatchingBioseq (BioseqPtr bsp,
					SeqMgrBioseqContextPtr bContext)
{
  Char          currentId[MAX_ID_LEN];
  UpdateDataPtr pUpdateData;
  Int2          result;
  SeqIdPtr      sip, sip_next;

  pUpdateData = (UpdateDataPtr) bContext->userdata;

  /* Get the string IDs for the current Bioseq */
  
  for (sip = bsp->id; sip != NULL; sip = sip->next) 
  {
    sip_next = sip->next;
    sip->next = NULL;
    SeqIdWrite (sip, currentId, PRINTID_TEXTID_ACC_ONLY,
	            sizeof (currentId) - 1);

    /* Compare it to the string ID of the new Bioseq */

    result = StringICmp (pUpdateData->newId, currentId);

    /* if TEXTID_ACC_ONLY doesn't match, try PRINTID_REPORT */
    if (result != 0)
    {
      SeqIdWrite (sip, currentId, PRINTID_REPORT,
	                sizeof (currentId) - 1);

      /* Compare it to the string ID of the new Bioseq */

      result = StringICmp (pUpdateData->newId, currentId);
    }

	  sip->next = sip_next;
    /* If they match, save the Bioseq and quit searching */

    if (0 == result) {
      pUpdateData->matchingBsp = bsp;
      return FALSE;
    }
  }

  /* Else continue */

  return TRUE;
}

static Boolean SkipProteinInNucUpdate (SeqEntryPtr sep, UpsDataPtr udp)
{
  Char       newId[MAX_ID_LEN];
  SeqIdPtr   sip;
  BioseqPtr  bsp;
  Boolean    rval = FALSE;
  MsgAnswer  ans;
  
  if (sep == NULL || ! IS_Bioseq (sep) || sep->data.ptrvalue == NULL || udp == NULL)
  {
    return FALSE;
  }
   
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (ISA_na (bsp->mol))
  {
    return FALSE;
  }
  
  sip = SeqIdFindWorst (bsp->id);
  SeqIdWrite (sip, newId, PRINTID_REPORT, sizeof (newId) - 1);
  ans = Message (MSG_YN, "Found a protein (%s) in the update file, expecting "
         "only nucleotides.  Do you want to skip this sequence and continue?", newId);
  if (ans == ANS_YES)
  {
    rval = TRUE;
    fprintf (udp->log_fp, "Skipped protein Bioseq (%s) in nucleotide update\n", newId);
	  udp->data_in_log = TRUE;
  }
  return rval;
}

static void RemoveUpdateSet (UpsDataPtr udp)
{
#if 0
  ObjMgrPtr omp;
  Int2      entityID;
  
  if (udp == NULL || udp->seqsubsep == NULL)
  {
    return;
  }
  entityID = SeqMgrGetEntityIDForSeqEntry (udp->seqsubsep);  
  udp->seqsubsep = SeqEntryFree (udp->seqsubsep);
	omp = ObjMgrGet ();	
	ObjMgrReapOne (omp);
	ObjMgrFreeCache (0);
  ObjMgrSendMsg(OM_MSG_DEL, entityID, 0, 0);
#endif
}

/*=====================================================================*/
/*                                                                     */
/* UpdateNextBioseqInFastaSet () - Reads in one Bioseq from a FASTA set*/
/*                                 file and updates the corresponding  */
/*                                 Bioseq in memory.                   */
/*                                                                     */
/*=====================================================================*/

static Int2 UpdateNextBioseqInFastaSet (UpsDataPtr udp)
{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Pointer       dataptr;
  Uint2         datatype;
  Char          errMsg[256];
  BioseqPtr     nbsp;
  SeqEntryPtr   nwsep = NULL;
  SeqEntryPtr   sep = NULL;
  SeqIdPtr      sip;
  SeqSubmitPtr  ssp;
  UpdateData    updateData;
  SeqEntryPtr   nthBspSep;
  BioseqPtr     nthBsp;
  Boolean       skip_prot_in_nuc;

  sep = GetTopSeqEntryForEntityID (udp->input_entityID);
  bsp = GetBioseqGivenIDs (udp->input_entityID,
			   udp->input_itemID,
			   udp->input_itemtype);

  updateData.matchingBsp = NULL;
  /* keep reading file until we find a sequence that matches
   * one that we have.
   */
  while (updateData.matchingBsp == NULL)
  {
    skip_prot_in_nuc = FALSE;
    if (udp->seqsubsep == NULL)
    {
      /* Read in one sequence from the file */
      dataptr = ReadAsnFastaOrFlatFile (udp->fp, &datatype, NULL, FALSE, FALSE,
		                   	                  TRUE, FALSE);      

      if (NULL == dataptr) 
      {
        FileClose (udp->fp);
        CloseOutSequenceUpdateLog (udp);
        RemoveUpdateSet (udp); 
        return FASTA_READ_DONE;
      }

      /* Convert the file data to a SeqEntry */
  
      if (datatype == OBJ_SEQENTRY)
        nwsep = (SeqEntryPtr) dataptr;
      else if (datatype == OBJ_BIOSEQ || datatype == OBJ_BIOSEQSET)
        nwsep = SeqMgrGetSeqEntryForData (dataptr);
      else if (datatype == OBJ_SEQSUB) 
      {
        ssp = (SeqSubmitPtr) dataptr;
        if (ssp != NULL && ssp->datatype == 1)
        {
          nwsep = (SeqEntryPtr) ssp->data;
          udp->seqsubsep = nwsep;
          udp->seqsubpos = 1; 
        }
      }
  
      if (nwsep == NULL) 
      {
        FileClose (udp->fp);
        ErrPostEx (SEV_ERROR, 0, 0, "Unable to convert file data into SeqEntry.");
        CloseOutSequenceUpdateLog (udp);
        return FASTA_READ_ERROR;
      }  

      /* Use the new SeqEntry to get a Bioseq */
  
      if (ISA_na (bsp->mol))
      {
        nbsp = FindNucBioseq (nwsep);
      }
      else 
      {
        nwsep = FindNthBioseq (nwsep, 1);
        if (nwsep == NULL || nwsep->choice != 1)
        {
          CloseOutSequenceUpdateLog (udp);
          return FASTA_READ_ERROR;
        }
        nbsp = (BioseqPtr) nwsep->data.ptrvalue;
      }
      if (nbsp == NULL) 
      {
        if (ISA_na (bsp->mol))
        {
          skip_prot_in_nuc = SkipProteinInNucUpdate (nwsep, udp);
        }
        if (! skip_prot_in_nuc)
        {
          FileClose (udp->fp);
          ErrPostEx (SEV_ERROR, 0, 0, "Unable to convert file data into Bioseq.");
          CloseOutSequenceUpdateLog (udp);
          return FASTA_READ_ERROR;
        }
      }
    }
    else
    {
      /* get the next Bioseq from the record */
      udp->seqsubpos ++;
      nthBspSep = FindNthBioseq (udp->seqsubsep, udp->seqsubpos);
      nbsp = NULL;
      while (nthBspSep != NULL && nbsp == NULL)
      {
        if (!IS_Bioseq (nthBspSep))
        {
          udp->seqsubpos++;
          nthBspSep = FindNthBioseq (udp->seqsubsep, udp->seqsubpos);
        }
        else
        {
          nthBsp = nthBspSep->data.ptrvalue;
          if (ISA_na (bsp->mol) && ISA_na (nthBsp->mol))
          {
            nbsp = nthBsp;
          }
          else if (ISA_aa (bsp->mol) && ISA_aa (nthBsp->mol))
          {
            nbsp = nthBsp;
          }
          else
          {
            udp->seqsubpos++;
            nthBspSep = FindNthBioseq (udp->seqsubsep, udp->seqsubpos);
          }
        }
      }
      if (nthBspSep == NULL)
      {
        RemoveUpdateSet (udp); 
        return FASTA_READ_DONE;
      }
    }
  
    if (!skip_prot_in_nuc)
    {
      /* Get the string ID for the new Bioseq so that we */
      /* can find a matching ID among current Bioseqs.   */
  
      sip = SeqIdFindWorst (nbsp->id);
      SeqIdWrite (sip, updateData.newId, PRINTID_REPORT,
	                sizeof (updateData.newId) - 1);

      /* Find the matching bioseq in the current sequence set */

      updateData.matchingBsp = NULL;
      if (2 == sep->choice ) 
      {
        bssp = (BioseqSetPtr) sep->data.ptrvalue;
        SeqMgrExploreBioseqs (0, (Pointer) bssp, &updateData, FindMatchingBioseq,
			                     TRUE, TRUE, TRUE);
      }
      else 
      {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        SeqMgrExploreBioseqs (0, (Pointer) bsp, &updateData, FindMatchingBioseq,
			                        TRUE, TRUE, TRUE);
      }


      if (updateData.matchingBsp == NULL) 
      {
        OpenSequenceUpdateLog (udp);
        if (udp->log_fp != NULL)
        {
          fprintf (udp->log_fp, "No Bioseq found with ID matching that of the"
	                     " one in the file (%s)\n", updateData.newId);
	        udp->data_in_log = TRUE;
        }
        sprintf (errMsg, "No Bioseq found with ID matching that of the"
                       " one in the file (%s)", updateData.newId);
        ErrPostEx (SEV_ERROR, 0, 0, errMsg);
      }
    }
  }

  /* Do the updating of the sequences */
  
  udp->oldbsp = updateData.matchingBsp;
  udp->newbsp = nbsp;

  if (! PrepareToUpdateSequences (udp))
  {
    return FASTA_READ_DONE;
  }

  return FASTA_READ_OK;
}
  
/*=====================================================================*/
/*                                                                     */
/* UpdateFastaSet () - Updates a set of sequence from a FASTA file     */
/*                     containing a set of sequences.                  */
/*                                                                     */
/*=====================================================================*/

static void UpdateOrExtendFastaSet (IteM i, Boolean do_update)
{
  BaseFormPtr   bfp;
  FILE         *fp;
  Char          path [PATH_MAX];
  UpsDataPtr    udp;

  /* Check parameters */

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif

  if (bfp == NULL)
    return;

  /* Read in the update data from a file */

  if (! GetInputFileName (path, sizeof (path),"","TEXT"))
    return;
  fp = FileOpen (path, "r");
  if (fp == NULL)
    return;

  /* Create data ptr */

  udp = (UpsDataPtr) MemNew (sizeof (UpsData));
  if (udp == NULL)
    return;

  udp->input_entityID = bfp->input_entityID;
  udp->input_itemID   = bfp->input_itemID;
  udp->input_itemtype = bfp->input_itemtype;
  udp->fp             = fp;
  udp->useGUI         = TRUE;
  udp->isSet          = TRUE;
  udp->convertPubs    = CONVERTPUBS_NO; /* was CONVERTPUBS_NOT_SET */
  udp->do_update      = do_update;
  udp->suppress_continue_msg = FALSE;
  udp->suppress_instant_refresh = FALSE;
  udp->log_fp         = NULL;
  udp->data_in_log    = FALSE;
  udp->transl_except_list = NULL;
  udp->aln1           = NULL;
  udp->aln2           = NULL;

  /* Update one Bioseq from the file.  Note that this chains */
  /* to the processing of the Bioseq after that, so that     */
  /* actually all Bioseqs are processed by this call.        */

  UpdateNextBioseqInFastaSet (udp);

}

extern void UpdateFastaSet (IteM i)
{
  UpdateOrExtendFastaSet (i, TRUE);
}

extern void ExtendFastaSet (IteM i)
{
  UpdateOrExtendFastaSet (i, FALSE);
}

typedef struct extendsequences
{
  FEATURE_FORM_BLOCK

  BioseqPtr  newbsp;
  Boolean    add_cit_sub;
  Boolean    extend5;
  FILE *     log_fp;
  Char       log_path [PATH_MAX];
  ValNodePtr sequence_list;
  Boolean    data_in_log;
  LisT       sequence_list_ctrl;
  ButtoN     add_cit_sub_btn;
  GrouP      extend_end;
  
} ExtendSequencesData, PNTR ExtendSequencesPtr;

static void ExtendAllSequencesInSetCallback (BioseqPtr bsp, Pointer userdata)
{
  ExtendSequencesPtr    esp;
  SeqIdPtr              sip, id_next;
  Char                  acc_str [256];
  CharPtr               origSeqStr;
  CharPtr               newSeqStr;
  CharPtr               mergedSeqStr;
  Int4                  mergedLen;
  ByteStorePtr          mergedBS;
  Int4                  offset;
  SeqMgrFeatContext     context;
  SeqFeatPtr            sfp;
  CdRegionPtr           crp;
  CodeBreakPtr          cbp;
  RnaRefPtr             rrp;
  tRNAPtr               trp;
  
  if (bsp == NULL || userdata == NULL) return;
  esp = (ExtendSequencesPtr) userdata;
  
  if (bsp == esp->newbsp) return;
  
  /* Get original and new sequences */

  origSeqStr = GetSequenceByBsp (bsp);
  newSeqStr = GetSequenceByBsp (esp->newbsp);
  
  /* create string to hold extended sequence */
  mergedLen =  StringLen (newSeqStr) + StringLen (origSeqStr);
  mergedSeqStr = (CharPtr) MemNew (mergedLen + 1);
    
  if (esp->extend5)
  {
    /* prepend the new sequence */
    sprintf (mergedSeqStr, "%s%s", newSeqStr, origSeqStr);
  }
  else
  {
    /* append the new sequence */
    sprintf (mergedSeqStr, "%s%s", origSeqStr, newSeqStr);
  }

  /* Convert the new sequence into a ByteStore */

  mergedBS = BSNew (mergedLen);
  BSWrite (mergedBS, (VoidPtr) mergedSeqStr, mergedLen);

  /* Replace the original sequence with the */
  /* new concatenated sequence.             */

  bsp->seq_data      = BSFree (bsp->seq_data);
  bsp->seq_data      = mergedBS;
  bsp->seq_data_type = Seq_code_iupacna;
  bsp->length        = mergedLen;

  /* shift the features downstream for 5' extension */
  offset = StringLen (newSeqStr);
  sip = SeqIdFindBest (bsp->id, 0);
  if (sip == NULL) return;
  if (esp->extend5 && offset > 0)
  {
    sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &context);
    while (sfp != NULL) {
      OffsetLocation (sfp->location, offset, sip);
      switch (sfp->data.choice) {
        case SEQFEAT_CDREGION :
          crp = (CdRegionPtr) sfp->data.value.ptrvalue;
          if (crp != NULL) {
            for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
              OffsetLocation (cbp->loc, offset, sip);
            }
          }
          break;
        case SEQFEAT_RNA :
          rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
          if (rrp != NULL && rrp->ext.choice == 2) {
            trp = (tRNAPtr) rrp->ext.value.ptrvalue;
            if (trp != NULL && trp->anticodon != NULL) {
              OffsetLocation (trp->anticodon, offset, sip);
            }
          }
          break;
        default :
          break;
      }
      sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &context);
    }
  }
  
  if (esp->add_cit_sub)
  {
    AddCitSubToUpdatedSequence ( bsp, esp->input_entityID);
  }
  
  if (esp->log_fp != NULL)
  {
    id_next = sip->next;
    sip->next = NULL;
    SeqIdWrite (sip, acc_str, PRINTID_REPORT, sizeof (acc_str) - 1);
    sip->next = id_next;
    if (esp->extend5)
    {
      fprintf (esp->log_fp, "Extended %s at 5' end\n", acc_str);
    }
    else
    {
      fprintf (esp->log_fp, "Extended %s at 3' end\n", acc_str);
    }
    esp->data_in_log = TRUE;  
  }
}

static void DoExtendAllSequencesInSet (ButtoN b)
{
  ExtendSequencesPtr    esp;
  SeqEntryPtr           topsep;
  ValNodePtr            sip_list, vnp;
  SeqIdPtr              sip;
  BioseqPtr             bsp;

  esp = (ExtendSequencesPtr) GetObjectExtra (b);
  if (esp == NULL)
  {
    return;
  }
  
  esp->add_cit_sub = GetStatus (esp->add_cit_sub_btn);
  if (GetValue (esp->extend_end) == 1)
  {
    esp->extend5 = TRUE;
  }
  else
  {
    esp->extend5 = FALSE;
  }
  sip_list = GetSelectedSequenceList (esp->sequence_list_ctrl);

  /* create file for log */
  TmpNam (esp->log_path);
  esp->log_fp = FileOpen (esp->log_path, "wb");
  
  topsep = GetTopSeqEntryForEntityID (esp->input_entityID);
  if (topsep == NULL)
    return;
  
  for (vnp = sip_list; vnp != NULL; vnp = vnp->next)
  {
    sip = (SeqIdPtr) vnp->data.ptrvalue;
    bsp = BioseqFind (sip);
    if (bsp != NULL)
    {
      ExtendAllSequencesInSetCallback (bsp, esp);
    }
  }

  if (esp->log_fp != NULL) 
  {
    FileClose (esp->log_fp);
    esp->log_fp = NULL;
    if (esp->data_in_log) 
    {
      LaunchGeneralTextViewer (esp->log_path, "Extended Sequences");
      esp->data_in_log = FALSE;
    }
    FileRemove (esp->log_path);  	
  }
  ObjMgrSetDirtyFlag (esp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, esp->input_entityID, 0, 0);  
  Remove (esp->form);
  Update ();
}

static void SelectAllSequencesForExtend (ButtoN b)
{
  ExtendSequencesPtr    esp;

  esp = (ExtendSequencesPtr) GetObjectExtra (b);
  if (esp == NULL)
  {
    return;
  }
  SelectAllSequencesInListCtrl (esp->sequence_list_ctrl);  
}

static void UnSelectAllSequencesForExtend (ButtoN b)
{
  ExtendSequencesPtr    esp;

  esp = (ExtendSequencesPtr) GetObjectExtra (b);
  if (esp == NULL)
  {
    return;
  }
  UnSelectAllSequencesInListCtrl (esp->sequence_list_ctrl);  
}


extern void ExtendAllSequencesInSet (IteM i)
{
  BaseFormPtr        bfp;
  FILE              *fp;
  Char               path [PATH_MAX];
  Pointer            dataptr;
  Uint2              datatype;
  SeqEntryPtr        nwsep, topsep;
  SeqSubmitPtr       ssp;
  BioseqPtr          nbsp;
  BioseqPtr          bsp;
  ExtendSequencesPtr esp;
  WindoW             w;
  GrouP              h, g, c;
  ButtoN             b;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif

  if (bfp == NULL)
    return;
  
  topsep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (topsep == NULL)
    return;
  
  bsp = GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID,
			   bfp->input_itemtype);
  if (bsp == NULL)
  {
    Message (MSG_ERROR, "You must select a bioseq");
    return;
  }

  /* Read in the update data from a file */

  if (! GetInputFileName (path, sizeof (path),"","TEXT"))
    return;
  fp = FileOpen (path, "r");
  if (fp == NULL)
  {
    Message (MSG_ERROR, "Unable to open %s", path);
    return;
  }
    
  /* Create data ptr */
  dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, FALSE,
		                 	              TRUE, FALSE);

  FileClose (fp);
  if (NULL == dataptr) 
  {
    Message (MSG_ERROR, "Unable to read sequence data from %s", path);
    return;
  }

  /* Convert the file data to a SeqEntry */
  
  if (datatype == OBJ_SEQENTRY)
    nwsep = (SeqEntryPtr) dataptr;
  else if (datatype == OBJ_BIOSEQ || datatype == OBJ_BIOSEQSET)
    nwsep = SeqMgrGetSeqEntryForData (dataptr);
  else if (datatype == OBJ_SEQSUB) 
  {
    ssp = (SeqSubmitPtr) dataptr;
    if (ssp != NULL && ssp->datatype == 1)
    {
      nwsep = (SeqEntryPtr) ssp->data;
    }
  }
  
  if (nwsep == NULL) 
  {
    ErrPostEx (SEV_ERROR, 0, 0, "Unable to convert file data into SeqEntry.");
    return;
  }  

  /* Use the new SeqEntry to get a Bioseq */
  
  if (ISA_na (bsp->mol))
  {
    nbsp = FindNucBioseq (nwsep);
  }
  else 
  {
    nwsep = FindNthBioseq (nwsep, 1);
    if (nwsep == NULL || nwsep->choice != 1)
    {
      return;
    }
    nbsp = (BioseqPtr) nwsep->data.ptrvalue;
  }
  
  if (nbsp == NULL) 
  {
    ErrPostEx (SEV_ERROR, 0, 0, "Unable to convert file data into Bioseq.");
    return;
  }
  
  esp = (ExtendSequencesPtr) MemNew (sizeof (ExtendSequencesData));
  if (esp == NULL) return;
  esp->newbsp = nbsp;
  
  w = FixedWindow (-50, -33, -10, -10, "Extend Sequences", NULL);

  SetObjectExtra (w, esp, StdCleanupFormProc);
  esp->form = (ForM) w;
  esp->input_entityID = bfp->input_entityID;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  esp->sequence_list_ctrl = MakeSequenceListControl (h, topsep, NULL, NULL, TRUE, TRUE);
  g = HiddenGroup (h, 2, 0, NULL);
  b = PushButton (g, "Select All", SelectAllSequencesForExtend);
  SetObjectExtra (b, esp, NULL);
  b = PushButton (g, "Unselect All", UnSelectAllSequencesForExtend);
  SetObjectExtra (b, esp, NULL);
  
  esp->extend_end = HiddenGroup (h, 2, 0, NULL);
  RadioButton (esp->extend_end, "5' end");
  RadioButton (esp->extend_end, "3' end");
  SetValue (esp->extend_end, 1);
    
  esp->add_cit_sub_btn = CheckBox (h, "Add Cit Subs to extended sequences", NULL);


  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoExtendAllSequencesInSet);
  SetObjectExtra (b, esp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  AlignObjects (ALIGN_CENTER, (HANDLE) esp->sequence_list_ctrl,
                (HANDLE) esp->add_cit_sub_btn, (HANDLE) g,
                (HANDLE) esp->extend_end,
                (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
}



/*=====================================================================*/
/*                                                                     */
/* NewUpdateSequence () - Updates a sequence from a file.              */
/*                                                                     */
/*=====================================================================*/

static Boolean DeltaLitOnly (BioseqPtr bsp)

{
  ValNodePtr  vnp;

  if (bsp == NULL || bsp->repr != Seq_repr_delta) return FALSE;
  for (vnp = (ValNodePtr)(bsp->seq_ext); vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == 1) return FALSE;
  }
  return TRUE;
}

static void NewUpdateOrExtendSequence (IteM i, Boolean do_update)
{
  MsgAnswer     ans;
  BaseFormPtr   bfp;
  BioseqPtr     bsp, nbsp;
  Pointer       dataptr;
  Uint2         datatype;
  FILE          *fp;
  Char          path [PATH_MAX];
  SeqEntryPtr   sep, nwsep = NULL;
  SeqSubmitPtr  ssp;
  UpsDataPtr    udp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif

  /* Get the current Bioseq */

  if (bfp == NULL)
    return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL)
    return;
  bsp = GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID,
			   bfp->input_itemtype);
  if (bsp == NULL)
    return;

  /* Read in the update data from a file */

  if (! GetInputFileName (path, sizeof (path),"","TEXT"))
    return;
  fp = FileOpen (path, "r");
  if (fp == NULL)
  {
    Message (MSG_ERROR, "Unable to open %s", path);
    return;
  }
  if (bsp->repr == Seq_repr_delta)
  {
    nwsep = ImportOneGappedSequence (fp);
    if (nwsep != NULL && IS_Bioseq (nwsep))
    {
      nbsp = (BioseqPtr) nwsep->data.ptrvalue;
      if (nbsp->repr != Seq_repr_delta)
      {
        if (ANS_CANCEL == Message (MSG_OKC, "You are updating a delta sequence with a non-delta sequence. "
                          "If you choose replace, your delta sequence will no longer be a delta sequence. "
                          "Do you wish to continue?"))
        {
          SeqEntryFree (nwsep);
          return;
        }
      }
      SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) nwsep->data.ptrvalue, nwsep);
      SeqMgrAddToBioseqIndex (nwsep->data.ptrvalue);
      dataptr = nwsep;
      datatype = OBJ_SEQENTRY;
    }
  }
  else
  {
    dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, FALSE,
		  		    TRUE, FALSE);
  }
  FileClose (fp);
  if (dataptr == NULL)
    return;

  /* Get a pointer to the new SeqEntry */

  if (datatype == OBJ_SEQENTRY)
    nwsep = (SeqEntryPtr) dataptr;
  else if (datatype == OBJ_BIOSEQ || datatype == OBJ_BIOSEQSET)
    nwsep = SeqMgrGetSeqEntryForData (dataptr);
  else if (datatype == OBJ_SEQSUB) {
    ssp = (SeqSubmitPtr) dataptr;
    if (ssp != NULL && ssp->datatype == 1)
      nwsep = (SeqEntryPtr) ssp->data;
  }

  if (nwsep == NULL)
    return;

  /* Use the new SeqEntry to get a Bioseq */

  if (ISA_na (bsp->mol))
    nbsp = FindNucBioseq (nwsep);
  else {
    nwsep = FindNthBioseq (nwsep, 1);
    if (nwsep == NULL || nwsep->choice != 1) return;
    nbsp = (BioseqPtr) nwsep->data.ptrvalue;
  }

  if (nbsp == NULL)
    return;

  /* convert delta lit to raw so sequence can be updated */
  if (!indexerVersion)
  {
    /* If original sequence is a not a raw sequence then */
    /* ask user for advice on how to proceed.      */

    if (bsp->repr != Seq_repr_raw) {
      ans = Message (MSG_YN, "Only raw sequences can be updated."
		     " Do you wish to proceed for copying features?");
      if (ans == ANS_NO)
        return;
    }
  }

  /* Create data ptr */

  udp = (UpsDataPtr) MemNew (sizeof (UpsData));
  if (udp == NULL)
    return;

  udp->input_entityID = bfp->input_entityID;
  udp->input_itemID   = bfp->input_itemID;
  udp->input_itemtype = bfp->input_itemtype;
  udp->oldbsp         = bsp;
  udp->newbsp         = nbsp;
  udp->fp             = NULL;
  udp->isSet          = FALSE;
  udp->useGUI         = TRUE;
  udp->convertPubs    = CONVERTPUBS_NO; /* was CONVERTPUBS_NOT_SET */
  udp->do_update      = do_update;
  udp->suppress_continue_msg = FALSE;
  udp->suppress_instant_refresh = FALSE;
  udp->log_fp         = NULL;
  udp->data_in_log    = FALSE;
  udp->transl_except_list = NULL;
  udp->aln1           = NULL;
  udp->aln2           = NULL;

  /* Do the updating of the sequences */

  PrepareToUpdateSequences (udp);
}

extern void NewUpdateSequence (IteM i)

{
  NewUpdateOrExtendSequence (i, TRUE);
}

extern void NewExtendSequence (IteM i)

{
  NewUpdateOrExtendSequence (i, FALSE);
}

extern void UpdateSeqAfterDownload
(BaseFormPtr bfp,
 BioseqPtr oldbsp,
 BioseqPtr newbsp)
{
  MsgAnswer   ans;
  UpsDataPtr  udp;

  /* convert delta lit to raw so sequence can be updated */

  if (oldbsp->repr == Seq_repr_delta && DeltaLitOnly (oldbsp)) {
    if (indexerVersion) {
      SegOrDeltaBioseqToRaw (oldbsp);
      ObjMgrSetDirtyFlag (oldbsp->idx.entityID, TRUE);
    } else {
      ans = Message (MSG_YN, "Only raw sequences can be updated."
		     " Do you wish to convert this delta sequence to raw?");
      if (ans == ANS_YES) {
        SegOrDeltaBioseqToRaw (oldbsp);
        ObjMgrSetDirtyFlag (oldbsp->idx.entityID, TRUE);
      }
    }
  }

  if (newbsp->repr == Seq_repr_delta && DeltaLitOnly (newbsp)) {
    if (indexerVersion) {
      SegOrDeltaBioseqToRaw (newbsp);
      ObjMgrSetDirtyFlag (newbsp->idx.entityID, TRUE);
    } else {
      ans = Message (MSG_YN, "Only raw sequences can be updated."
		     " Do you wish to convert this delta sequence to raw?");
      if (ans == ANS_YES) {
        SegOrDeltaBioseqToRaw (newbsp);
        ObjMgrSetDirtyFlag (newbsp->idx.entityID, TRUE);
      }
    }
  }
  /* Create data ptr */

  udp = (UpsDataPtr) MemNew (sizeof (UpsData));
  if (udp == NULL)
    return;

  udp->input_entityID = bfp->input_entityID;
  udp->input_itemID   = bfp->input_itemID;
  udp->input_itemtype = bfp->input_itemtype;
  udp->oldbsp         = oldbsp;
  udp->newbsp         = newbsp;
  udp->fp             = NULL;
  udp->isSet          = FALSE;
  udp->useGUI         = TRUE;
  udp->convertPubs    = CONVERTPUBS_NO; /* was CONVERTPUBS_NOT_SET */
  udp->do_update      = TRUE;
  udp->suppress_continue_msg = FALSE;
  udp->suppress_instant_refresh = FALSE;
  udp->log_fp         = NULL;
  udp->data_in_log    = FALSE;
  udp->transl_except_list = NULL;
  udp->aln1           = NULL;
  udp->aln2           = NULL;

  /* Do the updating of the sequences */

  PrepareToUpdateSequences (udp);
}


extern void ExtendSeqAfterDownload 
(BaseFormPtr bfp,
 BioseqPtr oldbsp,
 BioseqPtr newbsp)

{
  MsgAnswer   ans;
  UpsDataPtr  udp;

  /* convert delta lit to raw so sequence can be updated */

  if (oldbsp->repr == Seq_repr_delta && DeltaLitOnly (oldbsp)) {
    if (indexerVersion) {
      SegOrDeltaBioseqToRaw (oldbsp);
      ObjMgrSetDirtyFlag (oldbsp->idx.entityID, TRUE);
    } else {
      ans = Message (MSG_YN, "Only raw sequences can be extended."
		     " Do you wish to convert this delta sequence to raw?");
      if (ans == ANS_YES) {
        SegOrDeltaBioseqToRaw (oldbsp);
        ObjMgrSetDirtyFlag (oldbsp->idx.entityID, TRUE);
      }
    }
  }

  if (newbsp->repr == Seq_repr_delta && DeltaLitOnly (newbsp)) {
    if (indexerVersion) {
      SegOrDeltaBioseqToRaw (newbsp);
      ObjMgrSetDirtyFlag (newbsp->idx.entityID, TRUE);
    } else {
      ans = Message (MSG_YN, "Only raw sequences can be extended."
		     " Do you wish to convert this delta sequence to raw?");
      if (ans == ANS_YES) {
        SegOrDeltaBioseqToRaw (newbsp);
        ObjMgrSetDirtyFlag (newbsp->idx.entityID, TRUE);
      }
    }
  }
  /* Create data ptr */

  udp = (UpsDataPtr) MemNew (sizeof (UpsData));
  if (udp == NULL)
    return;

  udp->input_entityID = bfp->input_entityID;
  udp->input_itemID   = bfp->input_itemID;
  udp->input_itemtype = bfp->input_itemtype;
  udp->oldbsp         = oldbsp;
  udp->newbsp         = newbsp;
  udp->fp             = NULL;
  udp->isSet          = FALSE;
  udp->useGUI         = TRUE;
  udp->convertPubs    = CONVERTPUBS_NO; /* was CONVERTPUBS_NOT_SET */
  udp->do_update      = FALSE;
  udp->suppress_continue_msg = FALSE;
  udp->suppress_instant_refresh = FALSE;
  udp->log_fp         = NULL;
  udp->data_in_log    = FALSE;
  udp->transl_except_list = NULL;
  udp->aln1           = NULL;
  udp->aln2           = NULL;

  /* Do the updating of the sequences */

  PrepareToUpdateSequences (udp);
}



/* NEW FEATURE PROPAGATION SECTION */


extern void FixCdsAfterPropagate (
  IteM i
)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitFeaturesInSep (sep, (Pointer) bfp, DoFixCDS);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}


extern void NewFeaturePropagate (
  IteM i
)

{
  BaseFormPtr        bfp;
  BioseqPtr          bsp;
  ForM               f;
  SeqMgrFeatContext  fcontext;
  Uint2              itemID = 0;
  SeqAlignPtr        salp;
  SeqFeatPtr         sfp;
  SelStructPtr       sel;
  SeqEntryPtr        sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  bsp = GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
  if (bsp == NULL) {
    Message (MSG_OK, "You must target a single sequence in order to propagate");
    return;
  }
  sfp = GetNextFeatureOnSegOrMaster (bsp, NULL, 0, 0, &fcontext);
  if (sfp == NULL)
  {
    Message (MSG_OK, "The sequence must have features in order to propagate");
    return;
  }
  
  salp = FindAlignmentsForBioseq (bsp);

  if (salp == NULL) {
    Message (MSG_OK, "The record must have an alignment in order to propagate");
    return;
  }

  sel = ObjMgrGetSelected ();
  if (sel != NULL && sel->entityID == bfp->input_entityID &&
      sel->next == NULL && sel->itemtype == OBJ_SEQFEAT) {
    sfp = SeqMgrGetDesiredFeature (bfp->input_entityID, NULL, sel->itemID, 0, NULL, &fcontext);
    if (sfp != NULL && fcontext.bsp == bsp) {
      itemID = sel->itemID;
    }
  }

  f = FeaturePropagateForm (bsp, salp, itemID);
  if (f == NULL) return;
  Show (f);
  Select (f);
  SendHelpScrollMessage (helpForm, "Edit Menu", "Feature Propagate");
}

/* taken from ripen.c */

static Boolean PropagateFromGenomicProductSet (SeqEntryPtr sep)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  SeqEntryPtr   seqentry;
  ValNodePtr    sourcedescr;

  if (sep != NULL) {
    if (sep->choice == 2 && sep->data.ptrvalue != NULL) {

      /* get descriptors packaged on genomic product set */

      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      sourcedescr = bssp->descr;
      if (sourcedescr == NULL) return FALSE;
      if (bssp->_class == BioseqseqSet_class_gen_prod_set) {
        seqentry = bssp->seq_set;
        while (seqentry != NULL) {

          /* copy descriptors onto contig and nuc-prot sets */

          if (seqentry->data.ptrvalue != NULL) {
            if (seqentry->choice == 1) {
              bsp = (BioseqPtr) seqentry->data.ptrvalue;
              ValNodeLink (&(bsp->descr),
                           AsnIoMemCopy ((Pointer) sourcedescr,
                                         (AsnReadFunc) SeqDescrAsnRead,
                                         (AsnWriteFunc) SeqDescrAsnWrite));
            } else if (seqentry->choice == 2) {
              bssp = (BioseqSetPtr) seqentry->data.ptrvalue;
              ValNodeLink (&(bssp->descr),
                           AsnIoMemCopy ((Pointer) sourcedescr,
                                         (AsnReadFunc) SeqDescrAsnRead,
                                         (AsnWriteFunc) SeqDescrAsnWrite));
            }
          }
          seqentry = seqentry->next;
        }

        /* and free descriptors from genomic product set */

        bssp = (BioseqSetPtr) sep->data.ptrvalue;
        bssp->descr = SeqDescrFree (bssp->descr);
        return TRUE;
      }
    }
  }
  return FALSE;
}


static void CopyUserObject (SeqFeatPtr sfp, Pointer userdata)

{
  BioseqPtr          bsp;
  SeqMgrDescContext  dcontext;
  SeqDescrPtr        sdp;
  UserObjectPtr      uop;

  if (sfp->idx.subtype != FEATDEF_mRNA) return;

  /* find product cdna of mrna feature */

  bsp = BioseqFindFromSeqLoc (sfp->product);
  if (bsp == NULL) return;

  /* get closest user object descriptor */

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &dcontext);

  /* make sure evidence user object is no higher than nuc-prot set */

  if (sdp == NULL || dcontext.level > 1) return;

  /* copy user object, place on mrna feature */

  uop = (UserObjectPtr) sdp->data.ptrvalue;
  if (uop == NULL) return;
  uop = AsnIoMemCopy (uop,
                       (AsnReadFunc) UserObjectAsnRead,
                       (AsnWriteFunc) UserObjectAsnWrite);
  if (uop == NULL) return;

  /* should not be a user object there, but use combine function just in case */

  sfp->ext = CombineUserObjects (sfp->ext, uop);
}

static void CopyGene (SeqFeatPtr sfp, Pointer userdata)

{
  BioseqPtr          bsp;
  SeqMgrFeatContext  gcontext;
  SeqFeatPtr         gene, copy, temp;

  /* input mrna features are multi-interval on contig */

  if (sfp->data.choice != SEQFEAT_RNA) return;

  /* overlapping gene should be single interval on contig */

  gene = SeqMgrGetOverlappingGene (sfp->location, &gcontext);
  if (gene == NULL) return;

  /* find cdna product of mrna */

  bsp = BioseqFindFromSeqLoc (sfp->product);
  if (bsp == NULL) return;

  /* make sure gene feature doesn't already exist on cDNA */

  temp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, 0, &gcontext);
  if (temp != NULL) return;

  /* copy gene feature fields to paste into new gene feature */

  temp = AsnIoMemCopy (gene,
                       (AsnReadFunc) SeqFeatAsnRead,
                       (AsnWriteFunc) SeqFeatAsnWrite);
  if (temp == NULL) return;

  /* make new gene feature on full-length of cdna */

  copy = CreateNewFeatureOnBioseq (bsp, SEQFEAT_GENE, NULL);
  if (copy == NULL) return;

  /* paste fields from temp copy of original gene */

  copy->data.value.ptrvalue = temp->data.value.ptrvalue;
  copy->partial = temp->partial;
  copy->excpt = temp->excpt;
  copy->comment = temp->comment;
  copy->qual = temp->qual;
  copy->title = temp->title;
  copy->ext = temp->ext;
  copy->cit = temp->cit;
  copy->exp_ev = temp->exp_ev;
  copy->xref = temp->xref;
  copy->dbxref = temp->dbxref;
  copy->pseudo = temp->pseudo;
  copy->except_text = temp->except_text;

  MemFree (temp); /* do not SeqFeatFree */
}

static void CopyCDS (SeqFeatPtr sfp, Pointer userdata)

{
  BioseqPtr          bsp, cdna;
  Boolean            madeMrnaProtLink = FALSE;
  SeqMgrFeatContext  mcontext;
  SeqFeatPtr         mrna, cds, temp;
  SeqLocPtr          slp, cbslp;
  UserObjectPtr      uop;
  CdRegionPtr        crp;
  CodeBreakPtr       cbp;

  /* input cds features are single interval on cdna in nuc-prot sets, non on contig */

  if (sfp->idx.subtype != FEATDEF_CDS) return;

  /* map cds location from cdna to contig via mrna feature intervals on contig */

  cdna = BioseqFindFromSeqLoc (sfp->location);
  if (cdna == NULL) return;
  mrna = SeqMgrGetRNAgivenProduct (cdna, &mcontext);
  if (mrna == NULL) return;
  slp = productLoc_to_locationLoc (mrna, sfp->location);
  if (slp == NULL) return;
  bsp = BioseqFindFromSeqLoc (slp);
  if (bsp == NULL) return;

  /* copy cds feature fields to paste into new cds feature */

  temp = AsnIoMemCopy (sfp,
                       (AsnReadFunc) SeqFeatAsnRead,
                       (AsnWriteFunc) SeqFeatAsnWrite);
  if (temp == NULL) return;

  cds = CreateNewFeatureOnBioseq (bsp, SEQFEAT_CDREGION, NULL);
  if (cds == NULL) return;

  /* replace cdna location with contig location */

  cds->location = SeqLocFree (cds->location);
  cds->location = slp;

  /* now two cds products point to same protein bioseq, okay for genomic product set */

  cds->product = AsnIoMemCopy (sfp->product,
                               (AsnReadFunc) SeqLocAsnRead,
                               (AsnWriteFunc) SeqLocAsnWrite);

  /* paste fields from temp copy of original cds */

  cds->data.value.ptrvalue = temp->data.value.ptrvalue;

  /* update code breaks */
  crp = (CdRegionPtr) cds->data.value.ptrvalue;
  if (crp != NULL) {
    if (mrna->excpt) {
      /* Exception, e.g. unclassified transcription discrepancy: gap */
      if (crp->code_break != NULL) {
	crp->code_break =  CodeBreakFree (crp->code_break); /* XXX Remove list of code breaks; any better fix?*/
      }
    } else {
      for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
	cbslp = productLoc_to_locationLoc (mrna, cbp->loc);
	assert (cbslp != NULL);
	cbp->loc = SeqLocFree (cbp->loc);
	cbp->loc = cbslp;
      }
    }
  }

  /*
  if (crp != NULL) {
    crp->frame = 0;
  }
  */
  cds->partial = temp->partial;
  cds->excpt = temp->excpt;
  cds->comment = temp->comment;
  cds->qual = temp->qual;
  cds->title = temp->title;
  cds->ext = temp->ext;
  cds->cit = temp->cit;
  cds->exp_ev = temp->exp_ev;
  cds->xref = temp->xref;
  cds->dbxref = temp->dbxref;
  cds->pseudo = temp->pseudo;
  cds->except_text = temp->except_text;

  MemFree (temp); /* do not SeqFeatFree */

  bsp = BioseqFindFromSeqLoc (sfp->product);
  if (bsp == NULL) return;

  if (mrna->excpt && StringICmp (mrna->except_text, "unclassified transcription discrepancy") == 0) {
    cds->excpt = TRUE;
    cds->except_text = StringSave ("unclassified translation discrepancy");
  }

  /* evidence user object and mrna-cds link user object combined onto mrna on contig */

  if (! madeMrnaProtLink) {
    uop = CreateMrnaProteinLinkUserObject (bsp);
    mrna->ext = CombineUserObjects (mrna->ext, uop);
  }
}

static void MakeProtRefXref (SeqFeatPtr sfp)

{
  BioseqPtr          bsp;
  SeqMgrFeatContext  pcontext;
  SeqFeatPtr         prot;
  ProtRefPtr         prp;
  SeqFeatXrefPtr     xref;

  for (xref = sfp->xref; xref != NULL; xref = xref->next) {
    if (xref->data.choice == SEQFEAT_PROT) return;
  }

  bsp = BioseqFindFromSeqLoc (sfp->product);
  if (bsp == NULL) return;

  prot = SeqMgrGetBestProteinFeature (bsp, &pcontext);
  if (prot != NULL && prot->data.choice == SEQFEAT_PROT) {
    prp = (ProtRefPtr) prot->data.value.ptrvalue;
    if (prp != NULL) {
      xref = SeqFeatXrefNew ();
      if (xref != NULL) {
        xref->data.choice = SEQFEAT_PROT;
        xref->data.value.ptrvalue = AsnIoMemCopy ((Pointer) prp,
                                                  (AsnReadFunc) ProtRefAsnRead,
                                                  (AsnWriteFunc) ProtRefAsnWrite);
        xref->next = sfp->xref;
        sfp->xref = xref;
      }
    }
  }

}

static void RemoveSfpTitle (SeqFeatPtr sfp, Pointer userdata)

{
  if (sfp->title == NULL) return;
  sfp->title = MemFree (sfp->title);
}

static void MakeRedundantGPS (
  IteM i, Boolean doprop, Boolean justprotxref
)

{
  BaseFormPtr        bfp;
  BioseqPtr          bsp;
  SeqMgrFeatContext  fcontext;
  SeqEntryPtr        sep;
  SeqFeatPtr         sfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  if (! justprotxref) {  	
    VisitFeaturesInSep (sep, NULL, CopyGene);
    VisitFeaturesInSep (sep, NULL, CopyUserObject);
  }

  /* find genomic sequence */
  bsp = FindNucBioseq (sep);
  if (bsp != NULL) {
    if (! justprotxref) {
      sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
      if (sfp == NULL) {
        /* only copy CDS features up if none already on genomic */
        VisitFeaturesInSep (sep, NULL, CopyCDS);
        /* reindex with new CDS features */
        SeqMgrIndexFeatures (bfp->input_entityID, NULL);
      }
    }
    sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
    while (sfp != NULL) {
      MakeProtRefXref (sfp);
      sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext);
    }
  }
  VisitFeaturesInSep (sep, NULL, RemoveSfpTitle);
  if (doprop) {
    PropagateFromGenomicProductSet (sep);
  }

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

extern void MakeRedundantGPSwithProp (
  IteM i
)

{
  MakeRedundantGPS (i, TRUE, FALSE);
}

extern void MakeRedundantGPSnoProp (
  IteM i
)

{
  MakeRedundantGPS (i, FALSE, FALSE);
}

extern void MakeRedundantGPSjustXref (IteM i)

{
  MakeRedundantGPS (i, FALSE, TRUE);
}

static void FuseFeatJoins (SeqFeatPtr sfp, Pointer userdata)

{
  BioseqPtr  bsp;
  Boolean    partial5;
  Boolean    partial3;
  SeqLocPtr  slp;

  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return;

  slp = SeqLocFindNext (sfp->location, NULL);
  if (slp == NULL) return;
  slp = SeqLocFindNext (sfp->location, slp);
  if (slp == NULL) return;

  slp = SeqLocMerge (bsp, sfp->location, NULL, FALSE, TRUE, FALSE);
  if (slp == NULL) return;
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  sfp->location = SeqLocFree (sfp->location);
  sfp->location = slp;
  SetSeqLocPartial (sfp->location, partial5, partial3);
}

extern void FuseSlpJoins (IteM i)

{
  BaseFormPtr        bfp;
  SeqEntryPtr        sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitFeaturesInSep (sep, NULL, FuseFeatJoins);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void DoAuthorityPrefix (BioSourcePtr biop, Pointer userdata)

{
  size_t      len;
  OrgModPtr   omp;
  OrgNamePtr  onp;
  OrgRefPtr   orp;
  CharPtr     str;

  if (biop == NULL) return;
  orp = biop->org;
  if (orp == NULL) return;
  if (StringHasNoText (orp->taxname)) return;
  len = StringLen (orp->taxname);
  onp = orp->orgname;
  if (onp == NULL) return;
  for (omp = onp->mod; omp != NULL; omp = omp->next) {
    if (omp->subtype != ORGMOD_authority) continue;
    if (StringNCmp (omp->subname, orp->taxname, len) == 0) continue;
    str = MemNew (StringLen (omp->subname) + len + 3);
    if (str == NULL) continue;
    StringCpy (str, orp->taxname);
    StringCat (str, " ");
    StringCat (str, omp->subname);
    omp->subname = MemFree (omp->subname);
    omp->subname = str;
  }
}

extern void PrefixAuthorityWithOrganism (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  VisitBioSourcesInSep (sep, NULL, DoAuthorityPrefix);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

typedef struct addtranslexceptdata {
  FEATURE_FORM_BLOCK

  TexT        cds_comment;
  ButtoN      strict_checking_btn;
  CharPtr     cds_comment_txt;
  Boolean     strict_checking;
} AddTranslExceptData, PNTR AddTranslExceptPtr; 

static void AddTranslExcept (SeqFeatPtr sfp, Pointer userdata)

{
  CdRegionPtr        crp;
  Boolean            partial5, partial3;
  Int4               dna_len;
  CharPtr            bases;
  Int4               except_len;
  Int4               total;
  TransTablePtr      tbl = NULL;
  Int2               state;
  CharPtr            codon_start;
  BioseqPtr          bsp;
  Int4               from, to;
  Boolean            changed;
  CodeBreakPtr       new_cbp, last_cbp;
  Boolean            table_is_local;
  CharPtr            cds_comment = NULL;
  CharPtr            new_comment;
  Int4               comment_len;
  AddTranslExceptPtr ap;
  Boolean            use_strict = FALSE;

  if (sfp == NULL 
      || sfp->idx.subtype != FEATDEF_CDS
      || (crp = (CdRegionPtr)sfp->data.value.ptrvalue)== NULL) {
    return;
  }
  
  ap = (AddTranslExceptPtr) userdata;
  if (ap != NULL)
  {
    cds_comment = ap->cds_comment_txt;
  	use_strict = ap->strict_checking;
  }

  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  if (partial3) return;

  dna_len = SeqLocLen (sfp->location);
  if (partial5 && crp->frame > 1) {
    except_len = (dna_len - crp->frame + 1) % 3;
  } else {
    except_len = dna_len % 3;
  }
  if (except_len == 0) return;

  /* don't add code break if one already exists */
  last_cbp = crp->code_break;
  if (last_cbp != NULL && last_cbp->aa.choice == 1 && last_cbp->aa.value.intvalue == 42) {
    return;
  }
  while (last_cbp != NULL && last_cbp->next != NULL) {
    if (last_cbp->aa.choice == 1 && last_cbp->aa.value.intvalue == 42) {
      return;
    }
    last_cbp = last_cbp->next;
  }

  bases = ReadCodingRegionBases (sfp->location, dna_len, crp->frame, &total);
  if (bases == NULL) return;

  /* don't add transl_except if cds has valid stop codon */
  state = 0;
  codon_start = bases + StringLen (bases) - 6;
  if (codon_start < bases) {
    MemFree (bases);
    return;
  }
  tbl = GetTranslationTable (crp, &table_is_local);
  state = 0;
  state = NextCodonState (tbl, state, (Uint1)*codon_start);
  state = NextCodonState (tbl, state, (Uint1)*(codon_start + 1));
  state = NextCodonState (tbl, state, (Uint1)*(codon_start + 2));
  if (IsOrfStop (tbl, state, TTBL_TOP_STRAND)) {
    MemFree (bases);
    if (table_is_local) {
      TransTableFree (tbl);
    }
    return;
  }
  if (table_is_local) {
    TransTableFree (tbl);
    tbl = NULL;
  }
  
  if (use_strict)
  {
  	if (except_len == 2)
  	{
  	  if (toupper (*(codon_start + 3)) != 'T' || toupper(*(codon_start + 4)) != 'A')
  	  {
  	  	MemFree (bases);
  	  	return;
  	  }
  	}
  	else
  	{
  	  if (toupper (*(codon_start + 3)) != 'T')
  	  {
  	  	MemFree (bases);
  	  	return;
  	  }
  	}
  
  }
  else
  {
    /* don't add transl_except if exception location does not start with 'T' or 'N' */
    if (*(codon_start + 3) != 'T' && *(codon_start + 3) != 't'
        && *(codon_start + 3) != 'N' && *(codon_start + 3) != 'n') {
      MemFree (bases);
      return;
    }
  }
  MemFree (bases);

  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return;
  from = GetOffsetInBioseq (sfp->location, bsp, SEQLOC_LEFT_END);
  to = GetOffsetInBioseq (sfp->location, bsp, SEQLOC_RIGHT_END);

  new_cbp = CodeBreakNew ();
  new_cbp->aa.choice = 1;
  new_cbp->aa.value.intvalue = 42;
  new_cbp->loc = SeqLocMerge (BioseqFindFromSeqLoc (sfp->location), sfp->location, NULL,
                                   FALSE, FALSE, FALSE);
  if (SeqLocStrand (new_cbp->loc) == Seq_strand_minus) {
    new_cbp->loc = SeqLocDelete (new_cbp->loc, SeqLocId (new_cbp->loc),  from + except_len, to,
                  FALSE, &changed);
  } else {
    new_cbp->loc = SeqLocDelete (new_cbp->loc, SeqLocId (new_cbp->loc), from, to - except_len,
                  FALSE, &changed);
  }

  /* add code break to end of list */
  if (last_cbp == NULL) {
    crp->code_break = new_cbp;
  } else {
    last_cbp->next = new_cbp;
  }
  
  /* add comment if there is one */
  if (cds_comment != NULL)
  {
  	if (StringHasNoText (sfp->comment))
  	{
  	  sfp->comment = MemFree (sfp->comment);
  	  sfp->comment = StringSave (cds_comment);
  	}
  	else
  	{
      comment_len = StringLen (sfp->comment) + StringLen (cds_comment) + 2;
      new_comment = (CharPtr) MemNew (sizeof (Char) * comment_len);
      if (new_comment != NULL)
      {
      	StringCpy (new_comment, sfp->comment);
      	StringCat (new_comment, ";");
      	StringCat (new_comment, cds_comment);
      	sfp->comment = MemFree (sfp->comment);
      	sfp->comment = new_comment;
      }
  	}
  }
}

extern void GlobalAddTranslExcept (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  VisitFeaturesInSep (sep, NULL, AddTranslExcept);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void DoAddTranslExceptWithComment (ButtoN b)
{
  AddTranslExceptPtr ap;
  SeqEntryPtr  sep;

  ap = (AddTranslExceptPtr) GetObjectExtra (b);
  if (ap == NULL) return;
  Hide (ap->form);
  sep = GetTopSeqEntryForEntityID (ap->input_entityID);
  if (sep == NULL) return;
  ap->cds_comment_txt = SaveStringFromText (ap->cds_comment);
  ap->strict_checking = GetStatus (ap->strict_checking_btn);
  if (StringHasNoText (ap->cds_comment_txt))
  {
  	ap->cds_comment_txt = MemFree (ap->cds_comment_txt);
  }
  VisitFeaturesInSep (sep, ap, AddTranslExcept);
  ObjMgrSetDirtyFlag (ap->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ap->input_entityID, 0, 0);
  ap->cds_comment_txt = MemFree (ap->cds_comment_txt);
  Remove (ap->form);
}

extern void AddTranslExceptWithComment (IteM i)
{
  BaseFormPtr        bfp;
  AddTranslExceptPtr ap;
  WindoW             w;
  GrouP              h, g, c;
  ButtoN             b;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  ap = (AddTranslExceptPtr) MemNew (sizeof (AddTranslExceptData));
  w = FixedWindow (-50, -33, -10, -10, "Add Translation Exception", NULL);

  SetObjectExtra (w, ap, StdCleanupFormProc);
  ap->form = (ForM) w;
  ap->input_entityID = bfp->input_entityID;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, 2, 0, NULL);
  StaticPrompt (g, "CDS comment", 0, 0, programFont, 'c');
  ap->cds_comment = DialogText (g, "", 20, NULL);
  ap->strict_checking_btn = CheckBox (h, "Overhang must be T or TA", NULL);

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoAddTranslExceptWithComment);
  SetObjectExtra (b, ap, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) ap->strict_checking_btn, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
}

typedef struct updatealignmentlengths
{
  Int4 aln_length; /* length of alignment */
  Int4 old5;       /* length of 5' end of original sequence that is not 
                    * included in the alignment.
                    */
  Int4 old3;       /* length of 3' end of original sequence that is not
                    * included in the alignment.
                    */
  Int4 olda;       /* length of update sequence that is included in the alignment. */
  Int4 new5;       /* length of 5' end of update sequence that is not 
                    * included in the alignment.
                    */
  Int4 new3;       /* length of 3' end of update sequence that is not
                    * included in the alignment.
                    */
  Int4 newa;       /* length of update sequence that is included in the alignment. */
  
  Int4 log10_aln_length; 
  
  Int4 recomb1;
  Int4 recomb2;
} UpdateAlignmentLengthsData, PNTR UpdateAlignmentLengthsPtr;

/* holds information about an update pair */
typedef struct updatepair 
{
  SeqAlignPtr salp;
  Boolean     revcomp;
  BioseqPtr   orig_bsp;
  BioseqPtr   update_bsp;
} UpdatePairData, PNTR UpdatePairPtr;

/* These structures hold information about how to perform the updates.
 * The dialogs for displaying and collecting these structures are farther
 * down in the code.
 */
typedef enum {
  eSequenceUpdateNoChange = 1,
  eSequenceUpdateReplace,
  eSequenceUpdatePatch,
  eSequenceUpdateExtend5,
  eSequenceUpdateExtend3
} ESequenceUpdateType;

typedef enum {
  eFeatureUpdateNoChange = 1,
  eFeatureUpdateAllExceptDups,
  eFeatureUpdateAllMergeDups,
  eFeatureUpdateAllReplaceDups,
  eFeatureUpdateAll
} EFeatureUpdateType;

typedef enum {
  eFeatureRemoveNone = 1,
  eFeatureRemoveAligned,
  eFeatureRemoveNotAligned,
  eFeatureRemoveAll
} EFeatureRemoveType;

typedef struct submitterupdateoptions
{
  ESequenceUpdateType   sequence_update_type;
  EFeatureUpdateType    feature_update_type;
  EFeatureRemoveType    feature_remove_type;
  Boolean               ignore_alignment;
} SubmitterUpdateOptionsData, PNTR SubmitterUpdateOptionsPtr;

typedef struct indexer_options
{
  Boolean keep_protein_ids;
  Boolean add_cit_subs;
  Boolean update_quality_scores;
  Boolean update_proteins;
  Boolean truncate_proteins;
  Boolean extend_proteins3;
  Boolean extend_proteins5;
  Boolean correct_cds_genes;
} IndexerOptionsData, PNTR IndexerOptionsPtr;

typedef struct updateoptions
{
  SubmitterUpdateOptionsPtr submitter_opts;
  IndexerOptionsPtr         indexer_opts;
} UpdateOptionsData, PNTR UpdateOptionsPtr;

static UpdateOptionsPtr UpdateOptionsFree (UpdateOptionsPtr uop)
{
  if (uop != NULL)
  {
  	uop->submitter_opts = MemFree (uop->submitter_opts);
  	uop->indexer_opts = MemFree (uop->indexer_opts);
  	uop = MemFree (uop);
  }
  return uop;
}

/* These functions remove features from portions of the updated sequence:
 * RemoveFeatsInAlignedRegion
 * RemoveFeatsNotInAlignedRegion
 */
static void 
RemoveFeatsInAlignedRegion 
(BioseqPtr                 orig_bsp,
 UpdateAlignmentLengthsPtr ualp)

{
  SeqMgrFeatContext  context;
  Int4               left, right;
  SeqFeatPtr         sfp;

  if (orig_bsp == NULL || ualp == NULL) return;

  left = ualp->old5;
  right = ualp->old5 + ualp->olda;

  sfp = SeqMgrGetNextFeature (orig_bsp, NULL, 0, 0, &context);

  while (sfp != NULL) {

    if (context.right >= left && context.left <= right) {
      sfp->idx.deleteme = TRUE;
      MarkProductForDeletion (sfp->product);
    }

    sfp = SeqMgrGetNextFeature (orig_bsp, sfp, 0, 0, &context);
  }
}

static void 
RemoveFeatsNotInAlignedRegion 
(BioseqPtr                 orig_bsp,
 UpdateAlignmentLengthsPtr ualp)

{
  SeqMgrFeatContext  context;
  Int4               left, right;
  SeqFeatPtr         sfp;

  if (orig_bsp == NULL || ualp == NULL) return;

  left = ualp->old5;
  right = ualp->old5 + ualp->olda;

  sfp = SeqMgrGetNextFeature (orig_bsp, NULL, 0, 0, &context);

  while (sfp != NULL) {

    if (context.right < left || context.left > right) {
      sfp->idx.deleteme = TRUE;
      MarkProductForDeletion (sfp->product);
    }

    sfp = SeqMgrGetNextFeature (orig_bsp, sfp, 0, 0, &context);
  }
}

/* This function adjusts an alignment based on updates to the sequence. */
static Boolean ShiftAlignmentForUpdate 
(SeqAlignPtr sap,
 UpdateAlignmentLengthsPtr ualp,
 ESequenceUpdateType choice)

{
  DenseSegPtr  dsp;
  Int2         j;

  if (ualp == NULL || sap == NULL) return FALSE;

  AMFreeAllIndexes (sap);

  if (sap->segtype == SAS_DENSEG) {
    dsp = (DenseSegPtr) sap->segs;

    switch (choice) {
      case eSequenceUpdateExtend5 :
        /* adjust alignment 5' */
        if (dsp != NULL && dsp->lens != NULL && dsp->numseg > 0) {
          dsp->lens [dsp->numseg - 1] += ualp->old3;
        }
        break;
      case eSequenceUpdateExtend3 :
        /* adjust alignment 3' */
        if (dsp != NULL && dsp->lens != NULL && dsp->starts != NULL && dsp->numseg > 0) {
          dsp->lens [0] += ualp->old5;
          dsp->starts [0] = 0;
          dsp->starts [1] = 0;
          for (j = 1; j < dsp->numseg; j++) {
            if (dsp->starts [1 + j * 2] != -1) {
              dsp->starts [1 + j * 2] += ualp->old5 - ualp->new5;
            }
          }
        }
        break;
      case eSequenceUpdatePatch :
        /* adjust alignment patch */
        if (dsp != NULL && dsp->lens != NULL && dsp->starts != NULL && dsp->numseg > 0) {
          dsp->lens [dsp->numseg - 1] += ualp->old3;
          dsp->lens [0] += ualp->old5;
          dsp->starts [0] = 0;
          dsp->starts [1] = 0;
          for (j = 1; j < dsp->numseg; j++) {
            if (dsp->starts [1 + j * 2] != -1) {
              dsp->starts [1 + j * 2] += ualp->old5 - ualp->new5;
            }
          }
        }
        break;
      default :
        break;
    }
  }

  AlnMgr2IndexSingleChildSeqAlign (sap);

  return TRUE;
}

/* This function shifts the position of features on the updated sequence based on
 * the number of nucleotides added upstream.
 */
static Boolean 
ShiftFeaturesForUpdate 
(BioseqPtr orig_bsp,
 BioseqPtr update_bsp, 
 Int4 offset)

{
  ByteStorePtr       bs;
  CodeBreakPtr       cbp;
  SeqMgrFeatContext  context;
  CdRegionPtr        crp;
  Int4               len;
  RnaRefPtr          rrp;
  Uint1              seq_data_type;
  SeqFeatPtr         sfp;
  SeqIdPtr           sip;
  tRNAPtr            trp;

  if (orig_bsp == NULL || update_bsp == NULL)
  {
    return FALSE;
  }

  sip = SeqIdFindBest (orig_bsp->id, 0);
  if (sip == NULL) return FALSE;

  if (offset > 0) {
    sfp = SeqMgrGetNextFeature (orig_bsp, NULL, 0, 0, &context);
    while (sfp != NULL) {
      OffsetLocation (sfp->location, offset, sip);
      switch (sfp->data.choice) {
        case SEQFEAT_CDREGION :
          crp = (CdRegionPtr) sfp->data.value.ptrvalue;
          if (crp != NULL) {
            for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
              OffsetLocation (cbp->loc, offset, sip);
            }
          }
          break;
        case SEQFEAT_RNA :
          rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
          if (rrp != NULL && rrp->ext.choice == 2) {
            trp = (tRNAPtr) rrp->ext.value.ptrvalue;
            if (trp != NULL && trp->anticodon != NULL) {
              OffsetLocation (trp->anticodon, offset, sip);
            }
          }
          break;
        default :
          break;
      }
      sfp = SeqMgrGetNextFeature (orig_bsp, sfp, 0, 0, &context);
    }
  }

  /* switch bioseqs to finish extension */

  bs = orig_bsp->seq_data;
  orig_bsp->seq_data = update_bsp->seq_data;
  update_bsp->seq_data = bs;
  len = orig_bsp->length;
  orig_bsp->length = update_bsp->length;
  update_bsp->length = len;
  seq_data_type = orig_bsp->seq_data_type;
  orig_bsp->seq_data_type = update_bsp->seq_data_type;
  update_bsp->seq_data_type = seq_data_type;

  return TRUE;
}

static Int4 ReadLongDataToString (SeqPortPtr spp, Int4 how_much, CharPtr str)
{
  Int4          ctr = 0, left_to_read, read_this;

  if (spp == NULL || str == NULL)
  {
    return 0;
  }
  left_to_read = how_much;
  while (left_to_read > 0)
  {
    read_this = SeqPortRead (spp, (Uint1Ptr)(str + ctr), MIN (left_to_read, INT2_MAX));
    left_to_read -= read_this;
    ctr += read_this;
  }
  return ctr;
}

static SeqLitPtr SeqLitFromBioseq (BioseqPtr bsp, Int4 start, Int4 stop)
{
  CharPtr    str;
  SeqLitPtr  slip;
  Uint1      seqcode;
  SeqPortPtr spp;
  Int4       ctr;
  
  if (bsp == NULL || bsp->length < 1 || start < 0 || stop >= bsp->length)
  {
    return NULL;
  }
  
  str = (CharPtr) MemNew (sizeof (Char) * (size_t) (stop - start + 5));
  if (str == NULL)
    return NULL;
  
  if (ISA_na (bsp->mol))
  {
    seqcode = Seq_code_iupacna;
  }
  else
  {
    seqcode = Seq_code_iupacaa;
  }  
    
  slip = SeqLitNew ();

	slip->seq_data_type = seqcode;
  slip->length = stop - start + 1;
	slip->seq_data = BSNew (slip->length);

  /* copy in the data */
  spp = SeqPortNew(bsp, start, stop, Seq_strand_plus, seqcode);
  ctr = ReadLongDataToString (spp, stop - start + 1, str);
  spp = SeqPortFree (spp);
  str[ctr] = '\0';
       
  BSWrite (slip->seq_data, (VoidPtr) str, stop - start + 1);
  return slip;  
}

static DeltaSeqPtr DeltaSeqListFree (DeltaSeqPtr list)
{
  DeltaSeqPtr dsp;
  SeqLitPtr   slip;
  
  dsp = list;
  while (dsp != NULL)
  {
    slip = (SeqLitPtr) (dsp->data.ptrvalue);
    SeqLitFree (slip);
    dsp = dsp->next;
  }
  ValNodeFree (list);
  return NULL;  
}

static Boolean 
UpdateSequenceExtendDelta5Prime 
(UpdatePairPtr             upp,
 UpdateAlignmentLengthsPtr ualp)

{
  Int4        currold_pos;
  SeqLitPtr   slip, slip_new;
  DeltaSeqPtr dspold, dspnew, dspold_prev;
  Int4        seqstart;
  DeltaSeqPtr new_list = NULL;
  
  if (upp == NULL 
      || upp->orig_bsp == NULL 
      || upp->update_bsp == NULL
      || (upp->orig_bsp->repr != Seq_repr_delta && upp->update_bsp->repr != Seq_repr_delta)
      || (upp->orig_bsp->repr == Seq_repr_delta && upp->orig_bsp->seq_ext_type != 4)
      || (upp->update_bsp->repr == Seq_repr_delta && upp->update_bsp->seq_ext_type != 4))
  {
    return FALSE;
  }

  /* copy in entire new sequence */
  if (upp->update_bsp->repr == Seq_repr_delta)
  {
    dspnew = (DeltaSeqPtr) upp->update_bsp->seq_ext;
    while (dspnew != NULL)
    {
      if (dspnew->data.ptrvalue == NULL || dspnew->choice != 2)
      {
        return FALSE;
      }
	    slip = (SeqLitPtr) (dspnew->data.ptrvalue);
		  slip_new = (SeqLitPtr) AsnIoMemCopy (slip, (AsnReadFunc) SeqLitAsnRead,
		                                       (AsnWriteFunc) SeqLitAsnWrite);
    		                                     
      ValNodeAddPointer (&new_list, 2, slip_new);
	    dspnew = dspnew->next;
    }
  } 
  else
  {
    slip = SeqLitFromBioseq (upp->update_bsp, 0, upp->update_bsp->length - 1);
    if (slip == NULL)
    {
      return FALSE;
    }
    ValNodeAddPointer (&new_list, 2, slip);
  }
  
  /* now copy in old sequence, all if no alignment, after alignment if present */
  if (upp->orig_bsp->repr == Seq_repr_delta)
  {
    if (upp->salp == NULL)
    {
      ValNodeLink (&new_list, upp->orig_bsp->seq_ext);
      upp->orig_bsp->seq_ext = new_list;
      upp->orig_bsp->length += upp->update_bsp->length;
    }
    else
    {
      /* skip over old 5 and aligned area */
      dspold = (DeltaSeqPtr) upp->orig_bsp->seq_ext;
      currold_pos = 0;
      dspold_prev = NULL;
      while (dspold != NULL && currold_pos < ualp->old5 + ualp->olda)
      {
        seqstart = currold_pos;
        if (dspold->data.ptrvalue == NULL || dspold->choice != 2)
        {
          return FALSE;
        }
        slip = (SeqLitPtr) (dspold->data.ptrvalue);
	      currold_pos += slip->length;
		    if (currold_pos > ualp->old5 + ualp->olda)
  		  {
          SplitDeltaSeq (dspold, ualp->old5 + ualp->olda - seqstart);
          slip = (SeqLitPtr) (dspold->data.ptrvalue);
          currold_pos = ualp->old5 + ualp->olda;		  
        }
  		  dspold_prev = dspold;
    	  dspold = dspold->next;
      }
      if (dspold_prev == NULL)
      {
        ValNodeLink (&new_list, upp->orig_bsp->seq_ext);
      }
      else
      {
        /* attach remainder of list to new list */
        ValNodeLink (&new_list, dspold);
        dspold_prev->next = NULL;
        
        /* free dsps in old sequence that are being replaced */
        upp->orig_bsp->seq_ext = DeltaSeqListFree (upp->orig_bsp->seq_ext);
      
        /* put new list in place */
        upp->orig_bsp->seq_ext = new_list;
      }
      upp->orig_bsp->length = ualp->new5 + ualp->newa + ualp->old3;
    }
  }
  else
  {
    if (upp->salp == NULL)
    {
      slip = SeqLitFromBioseq (upp->orig_bsp, 0, upp->orig_bsp->length - 1);
    }
    else
    {
      slip = SeqLitFromBioseq (upp->orig_bsp, ualp->old5 + ualp->olda, upp->orig_bsp->length - 1);
    }
    if (slip == NULL)
    {
      return FALSE;
    }
    ValNodeAddPointer (&new_list, 2, slip);
    upp->orig_bsp->seq_data = BSFree (upp->orig_bsp->seq_data);
    upp->orig_bsp->seq_ext = new_list;
    upp->orig_bsp->repr = Seq_repr_delta;
    
    upp->orig_bsp->length = ualp->new5 + ualp->newa + ualp->old3;
  }
    
  return TRUE;  
}

static Boolean 
UpdateSequenceExtend5Prime 
(UpdatePairPtr             upp,
 UpdateAlignmentLengthsPtr ualp)

{
  ByteStorePtr  bs;
  Int4          newlen;
  CharPtr       str;
  SeqPortPtr    spp;
  Uint1         seqcode;
  Int4          ctr;
  Boolean       rval;
  Int4          shift_len = 0;

  if (upp == NULL || upp->orig_bsp == NULL || upp->update_bsp == NULL || ualp == NULL)
  {
    return FALSE;
  }
  
  if (upp->orig_bsp->repr == Seq_repr_delta
      || upp->update_bsp->repr == Seq_repr_delta)
  {
    return UpdateSequenceExtendDelta5Prime (upp, ualp);
  }

  if (upp->salp == NULL)
  {
    /* add to 5' end */
    newlen = upp->update_bsp->length + upp->orig_bsp->length;
    shift_len = upp->update_bsp->length;
  }
  else
  {
    /* construct replacement sequence by recombining between old and overlap */
    newlen = ualp->new5 + ualp->newa + ualp->old3;
    shift_len = ualp->new5;
  }

  str = (CharPtr) MemNew (sizeof (Char) * (size_t) (newlen + 5));
  if (str == NULL)
    return FALSE;

  if (ISA_na (upp->orig_bsp->mol))
  {
    seqcode = Seq_code_iupacna;
  }
  else
  {
    seqcode = Seq_code_iupacaa;
  }

  if (upp->salp == NULL)
  {
    /* add new sequence */
    spp = SeqPortNew(upp->update_bsp, 0, upp->update_bsp->length - 1, Seq_strand_plus, seqcode);
    ctr = ReadLongDataToString (spp, upp->update_bsp->length, str);
    spp = SeqPortFree (spp);
    /* add old sequence */
    spp = SeqPortNew(upp->orig_bsp, 0,upp-> orig_bsp->length - 1, Seq_strand_plus, seqcode);
    ctr += ReadLongDataToString(spp, upp->orig_bsp->length, str + ctr);
    spp = SeqPortFree (spp);
  }
  else
  {
    /* take new 5' and aligned middle */
    spp = SeqPortNew(upp->update_bsp, 0, ualp->new5 + ualp->newa - 1, Seq_strand_plus, seqcode);
    ctr = ReadLongDataToString (spp, ualp->new5 + ualp->newa, str);
    spp = SeqPortFree (spp);
    /* take old 3' end */
    if (ualp->old3 > 0)
    {
      spp = SeqPortNew(upp->orig_bsp, ualp->old5 + ualp->olda, upp->orig_bsp->length - 1, Seq_strand_plus, seqcode);
      ctr += ReadLongDataToString (spp, ualp->old3, str + ctr);
      spp = SeqPortFree (spp);
    }
  }
  
  str[ctr] = '\0';
       
  bs = BSNew (newlen);
  BSWrite (bs, (VoidPtr) str, newlen);

  if (bs != NULL && BSLen (bs) < 1) {
    bs = BSFree (bs);
  }
  if (bs == NULL) return FALSE;

  /* overlap turned into replacement sequence */

  upp->update_bsp->seq_data = BSFree (upp->update_bsp->seq_data);
  upp->update_bsp->seq_data = bs;
  upp->update_bsp->seq_data_type = seqcode;
  upp->update_bsp->length = newlen;


  if (upp->salp == NULL)
  {
    /* then finish by replacing with new sequence */
    rval = ShiftFeaturesForUpdate (upp->orig_bsp, upp->update_bsp, shift_len);
  }
  else
  {
    /* adjust alignment and reindex */
    rval = ShiftAlignmentForUpdate (upp->salp, ualp, eSequenceUpdateExtend5);
    if (rval)
    {
      /* then finish by replacing with new sequence */
      ReplaceOneSequence (upp->salp, upp->orig_bsp, upp->update_bsp);
    }
  }

  return rval;
}

static Boolean 
UpdateSequenceExtendDelta3Prime 
(UpdatePairPtr             upp,
 UpdateAlignmentLengthsPtr ualp)

{
  Int4        currnew_pos = 0, currold_pos;
  SeqLitPtr   slip, slip_new;
  DeltaSeqPtr dspold, dspnew, dspold_prev;
  Int4        seqstart;
  DeltaSeqPtr new_list = NULL;
  
  if (upp == NULL 
      || upp->orig_bsp == NULL 
      || upp->update_bsp == NULL
      || (upp->orig_bsp->repr != Seq_repr_delta && upp->update_bsp->repr != Seq_repr_delta)
      || (upp->orig_bsp->repr == Seq_repr_delta && upp->orig_bsp->seq_ext_type != 4)
      || (upp->update_bsp->repr == Seq_repr_delta && upp->update_bsp->seq_ext_type != 4))
  {
    return FALSE;
  }

  /* retain 5' end of original sequence */
  if (upp->orig_bsp->repr == Seq_repr_delta)
  {
    if (upp->salp == NULL)
    {
      /* do nothing, keep the entire sequence */
    }
    else
    {
      /* discard the portion after the old 5' */
      dspold_prev = NULL;
      dspold = (DeltaSeqPtr) upp->orig_bsp->seq_ext;
      currold_pos = 0;
      while (dspold != NULL && currold_pos < ualp->old5)
      {
        seqstart = currold_pos;
        if (dspold->data.ptrvalue == NULL || dspold->choice != 2)
        {
          return FALSE;
        }
        slip = (SeqLitPtr) (dspold->data.ptrvalue);
	      currold_pos += slip->length;
		    if (currold_pos > ualp->old5)
  		  {
          SplitDeltaSeq (dspold, ualp->old5 - seqstart);
          slip = (SeqLitPtr) (dspold->data.ptrvalue);
          currold_pos = ualp->old5;		  
        }
  		  dspold_prev = dspold;
    	  dspold = dspold->next;
      }
      
      if (dspold_prev == NULL)
      {
        /* discard the entire sequence */
        upp->orig_bsp->seq_ext = DeltaSeqListFree (upp->orig_bsp->seq_ext); 
      }
      else
      {
        dspold_prev->next = DeltaSeqListFree (dspold_prev->next);
      }
    }
  }
  else
  {
    if (upp->salp == NULL)
    {
      slip = SeqLitFromBioseq (upp->orig_bsp, 0, upp->orig_bsp->length - 1);      
    }
    else
    {
      slip = SeqLitFromBioseq (upp->orig_bsp, 0, ualp->old5);
    }
    upp->orig_bsp->seq_data = BSFree (upp->orig_bsp->seq_data);
    new_list = upp->orig_bsp->seq_ext;
    ValNodeAddPointer (&new_list, 2, slip);
    upp->orig_bsp->seq_ext = new_list;
    upp->orig_bsp->repr = Seq_repr_delta;
  }
  
  
  /* now add in new sequence */
  if (upp->update_bsp->repr == Seq_repr_delta)
  {
    /* skip over new5 */
    dspnew = (DeltaSeqPtr) upp->update_bsp->seq_ext;
    currnew_pos = 0;
    if (upp->salp != NULL)
    {
      while (dspnew != NULL && currnew_pos < ualp->new5)
      {
        seqstart = currnew_pos;
        if (dspnew->data.ptrvalue == NULL || dspnew->choice != 2)
        {
          return FALSE;
        }
        slip = (SeqLitPtr) (dspnew->data.ptrvalue);
	      currnew_pos += slip->length;
		    if (currnew_pos > ualp->new5)
  		  {
          SplitDeltaSeq (dspnew, ualp->old5 - seqstart);
          slip = (SeqLitPtr) (dspnew->data.ptrvalue);
          currnew_pos = ualp->new5;		  
        }
    	  dspnew = dspnew->next;
      }
    }
    
    while (dspnew != NULL)
    {
      if (dspnew->data.ptrvalue == NULL || dspnew->choice != 2)
      {
        return FALSE;
      }
      slip = (SeqLitPtr) (dspnew->data.ptrvalue);
		  slip_new = (SeqLitPtr) AsnIoMemCopy (slip, (AsnReadFunc) SeqLitAsnRead,
		                                       (AsnWriteFunc) SeqLitAsnWrite);
      new_list = upp->orig_bsp->seq_ext;  		                                     
      ValNodeAddPointer (&new_list, 2, slip_new);
      upp->orig_bsp->seq_ext = new_list;
      dspnew = dspnew->next;
    }
  }
  else
  {
    if (upp->salp == NULL)
    {
      slip = SeqLitFromBioseq (upp->update_bsp, 0, upp->update_bsp->length - 1);      
    }
    else
    {
      slip = SeqLitFromBioseq (upp->orig_bsp, ualp->new5, upp->update_bsp->length - 1);
    }
    new_list = upp->orig_bsp->seq_ext;
    ValNodeAddPointer (&new_list, 2, slip);    
    upp->orig_bsp->seq_ext = new_list;
  }
  
  if (upp->salp == NULL)
  {
    upp->orig_bsp->length = upp->orig_bsp->length + upp->update_bsp->length;
  }
  else
  {
    upp->orig_bsp->length = ualp->old5 + ualp->newa + ualp->new3;
  }
    
  return TRUE;  
}

static Boolean 
UpdateSequenceExtend3Prime 
(UpdatePairPtr             upp,
 UpdateAlignmentLengthsPtr ualp)
{
  ByteStorePtr  bs;
  Int4          newlen;
  CharPtr       str;
  SeqPortPtr    spp;
  Uint1         seqcode;
  Int4          ctr;
  Boolean       rval;

  if (upp == NULL || upp->orig_bsp == NULL || upp->update_bsp == NULL || ualp == NULL)
  {
    return FALSE;
  }
  
  if (upp->orig_bsp->repr == Seq_repr_delta
      || upp->update_bsp->repr == Seq_repr_delta)
  {
    return UpdateSequenceExtendDelta3Prime (upp, ualp);
  }

  /* construct replacement sequence by recombining between old and overlap */
  if (upp->salp == NULL)
  {
    newlen = upp->orig_bsp->length + upp->update_bsp->length;
  }
  else
  {
    newlen = ualp->old5 + ualp->newa + ualp->new3;
  }
  
  str = (CharPtr) MemNew (sizeof (Char) * (size_t) (newlen + 5));
  if (str == NULL)
    return FALSE;

  if (ISA_na (upp->orig_bsp->mol))
  {
    seqcode = Seq_code_iupacna;
  }
  else
  {
    seqcode = Seq_code_iupacaa;
  }

  if (upp->salp == NULL)
  {
    /* add old sequence */
    spp = SeqPortNew(upp->orig_bsp, 0, upp->orig_bsp->length - 1, Seq_strand_plus, seqcode);
    ctr = ReadLongDataToString (spp, upp->orig_bsp->length, str);
    spp = SeqPortFree (spp);
    /* add new sequence */
    spp = SeqPortNew(upp->update_bsp, 0, upp->update_bsp->length - 1, Seq_strand_plus, seqcode);
    ctr += ReadLongDataToString (spp, upp->update_bsp->length, str + ctr);
    spp = SeqPortFree (spp);
  }
  else
  {
    /* take old 5'  */
    spp = SeqPortNew(upp->orig_bsp, 0, ualp->old5 - 1, Seq_strand_plus, seqcode);
    ctr = ReadLongDataToString(spp, ualp->old5, str);
    spp = SeqPortFree (spp);
    /* take aligned middle and new 3' end */
    spp = SeqPortNew(upp->update_bsp, ualp->new5, upp->update_bsp->length - 1, Seq_strand_plus, seqcode);
    ctr += ReadLongDataToString (spp, ualp->newa + ualp->new3, str + ctr);
    spp = SeqPortFree (spp);
  }
  
  str[ctr] = '\0';

  bs = BSNew (newlen);
  BSWrite (bs, (VoidPtr) str, newlen);

  if (bs != NULL && BSLen (bs) < 1) {
    bs = BSFree (bs);
  }
  if (bs == NULL) return FALSE;

  /* overlap turned into replacement sequence */

  upp->update_bsp->seq_data = BSFree (upp->update_bsp->seq_data);
  upp->update_bsp->seq_data = bs;
  upp->update_bsp->seq_data_type = Seq_code_iupacna;
  upp->update_bsp->length = newlen;

  if (upp->salp == NULL)
  {
    /* then finish by replacing with new sequence */
    rval = ShiftFeaturesForUpdate (upp->orig_bsp, upp->update_bsp, 0);
  }
  else
  {
    /* adjust alignment and reindex */
    rval = ShiftAlignmentForUpdate (upp->salp, ualp, eSequenceUpdateExtend3);
    if (rval)
    {
      /* then finish by replacing with new sequence */
      ReplaceOneSequence (upp->salp, upp->orig_bsp, upp->update_bsp);
    }
  }

  return rval;
}

static Boolean RawSequencePatchOk (BioseqPtr orig_bsp, BioseqPtr update_bsp)
{
  Boolean rval = TRUE;
  
  if (orig_bsp == NULL || update_bsp == NULL
      || orig_bsp->repr != Seq_repr_raw
      || update_bsp->repr != Seq_repr_raw)
  {
    rval = FALSE;
  }

  return rval;
}

/* this replaces just the middle */
static Boolean 
UpdateSequencePatchRaw
(UpdatePairPtr             upp,
 UpdateAlignmentLengthsPtr ualp)

{
  ByteStorePtr  bs;
  Int4          newlen;
  CharPtr       str;
  SeqPortPtr    spp;
  Uint1         seqcode;
  Int4          ctr;

  if (upp == NULL || upp->orig_bsp == NULL || upp->update_bsp == NULL 
      || upp->salp == NULL || ualp == NULL)
  {
    return FALSE;
  }

  newlen = ualp->old5 + ualp->newa + ualp->old3;
  str = (CharPtr) MemNew (sizeof (Char) * (size_t) (newlen + 5));
  if (str == NULL) return FALSE;

  if (ISA_na (upp->orig_bsp->mol))
  {
    seqcode = Seq_code_iupacna;
  }
  else
  {
    seqcode = Seq_code_iupacaa;
  }
  
  /* construct replacement sequence by double recombination */

  /* take old 5'  */
  spp = SeqPortNew(upp->orig_bsp, 0, ualp->old5 - 1, Seq_strand_plus, seqcode);
  ctr = SeqPortRead(spp, (Uint1Ptr)str, ualp->old5);
  spp = SeqPortFree (spp);
  /* take aligned middle */
  spp = SeqPortNew(upp->update_bsp, ualp->new5, ualp->new5 + ualp->newa - 1, Seq_strand_plus, seqcode);
  ctr += SeqPortRead(spp, (Uint1Ptr)(str + ctr), ualp->newa);
  spp = SeqPortFree (spp);
  /* take old 3' */
  spp = SeqPortNew(upp->orig_bsp, ualp->old5 + ualp->olda, upp->orig_bsp->length - 1, Seq_strand_plus, seqcode);
  ctr += SeqPortRead(spp, (Uint1Ptr)(str + ctr), ualp->old3);
  spp = SeqPortFree (spp);
  
  str[ctr] = '\0';
  
  bs = BSNew (newlen);
  BSWrite (bs, (VoidPtr) str, newlen);

  if (bs != NULL && BSLen (bs) < 1) {
    bs = BSFree (bs);
  }
  if (bs == NULL) return FALSE;

  /* overlap turned into replacement sequence */

  upp->update_bsp->seq_data = BSFree (upp->update_bsp->seq_data);
  upp->update_bsp->seq_data = bs;
  upp->update_bsp->seq_data_type = Seq_code_iupacna;
  upp->update_bsp->length = newlen;
  return TRUE;  
}

static Boolean DeltaSequencePatchOk (BioseqPtr orig_bsp, BioseqPtr update_bsp)
{
  Boolean rval = TRUE;
  
  if (orig_bsp == NULL || update_bsp == NULL
      || orig_bsp->repr != Seq_repr_delta || update_bsp->repr != Seq_repr_delta
      || orig_bsp->seq_ext_type != 4 || update_bsp->seq_ext_type != 4)
  {
    rval = FALSE;
  }

  return rval;
}

/* This function will patch a delta sequence with another delta sequence.
 * The pieces in the overlap from the old sequence will be replaced by pieces
 * in the overlap from the new sequence.
 */
static Boolean UpdateSequencePatchDelta 
(UpdatePairPtr             upp,
 UpdateAlignmentLengthsPtr ualp)

{
  Int4        currnew_pos = 0, currold_pos;
  SeqLitPtr   slip, slip_new;
  DeltaSeqPtr dspold, dspnew;
  Int4        seqstart;
  DeltaSeqPtr new_list = NULL;
  
  if (upp == NULL)
  {
    return FALSE;
  }
  if (! DeltaSequencePatchOk (upp->orig_bsp, upp->update_bsp)
      || upp->salp == NULL || ualp == NULL)
  {
    return FALSE;
  }

  /* keep old 5' end intact */
  currold_pos = 0;
  seqstart = 0;
  dspold = (DeltaSeqPtr) upp->orig_bsp->seq_ext;
  while (dspold != NULL && currold_pos < ualp->old5)
  {
    seqstart = currold_pos;
    if (dspold->data.ptrvalue == NULL || dspold->choice != 2)
    {
      return FALSE;
    }
    slip = (SeqLitPtr) (dspold->data.ptrvalue);
	  currold_pos += slip->length;
		if (currold_pos > ualp->old5)
		{
      SplitDeltaSeq (dspold, ualp->old5 - seqstart);
      slip = (SeqLitPtr) (dspold->data.ptrvalue);
      currold_pos = ualp->old5;		  
		}
		slip_new = (SeqLitPtr) AsnIoMemCopy (slip, (AsnReadFunc) SeqLitAsnRead,
		                                     (AsnWriteFunc) SeqLitAsnWrite);
		ValNodeAddPointer (&new_list, 2, slip_new);
	  dspold = dspold->next;
  }
  
  /* skip over new 5' end */
  currnew_pos = 0;
  seqstart = 0;
  dspnew = (DeltaSeqPtr) upp->update_bsp->seq_ext;
  while (dspnew != NULL && currnew_pos < ualp->new5)
  {
    seqstart = currold_pos;
    if (dspnew->data.ptrvalue == NULL || dspnew->choice != 2)
    {
      return FALSE;
    }
	  slip = (SeqLitPtr) (dspnew->data.ptrvalue);
	  currnew_pos += slip->length;
	  if (currnew_pos > ualp->new5)
	  {
      SplitDeltaSeq (dspnew, ualp->new5 - seqstart);
      currnew_pos = ualp->new5;
	  }
	  dspnew = dspnew->next;
  }
  
  /* copy in new overlap */
  while (dspnew != NULL && currnew_pos < ualp->new5 + ualp->newa)
  {
    seqstart = currold_pos;
    if (dspnew->data.ptrvalue == NULL || dspnew->choice != 2)
    {
      return FALSE;
    }
	  slip = (SeqLitPtr) (dspnew->data.ptrvalue);
	  currnew_pos += slip->length;
		if (currnew_pos > ualp->new5 + ualp->newa)
		{
      SplitDeltaSeq (dspnew, ualp->new5 + ualp->newa - seqstart);
      slip = (SeqLitPtr) (dspnew->data.ptrvalue);
      currnew_pos = ualp->new5 + ualp->newa;		  
		}
		slip_new = (SeqLitPtr) AsnIoMemCopy (slip, (AsnReadFunc) SeqLitAsnRead,
		                                     (AsnWriteFunc) SeqLitAsnWrite);
		ValNodeAddPointer (&new_list, 2, slip_new);
		dspnew = dspnew->next;
  }
  
  /* skip over old overlap */
  
  while (dspold != NULL && currold_pos < ualp->old5 + ualp->olda)
  {
    seqstart = currold_pos;
    if (dspold->data.ptrvalue == NULL || dspold->choice != 2)
    {
      return FALSE;
    }
    slip = (SeqLitPtr) (dspold->data.ptrvalue);
    currold_pos += slip->length;
    if (currold_pos > ualp->old5 + ualp->olda)
    {
      SplitDeltaSeq (dspold, ualp->new5 + ualp->newa - seqstart);
      currold_pos = ualp->old5 + ualp->olda;		        
    }
    dspold = dspold->next;
  }
  
  /* copy in old 3' */
  
  while (dspold != NULL)
  {
    if (dspold->data.ptrvalue == NULL || dspold->choice != 2)
    {
      return FALSE;
    }
    slip = (SeqLitPtr) (dspold->data.ptrvalue);
		slip_new = (SeqLitPtr) AsnIoMemCopy (slip, (AsnReadFunc) SeqLitAsnRead,
		                                     (AsnWriteFunc) SeqLitAsnWrite);
		ValNodeAddPointer (&new_list, 2, slip_new);
		dspold = dspold->next;
  }
  
  /* free newbsp's old SeqLit List */
  for (dspnew = (DeltaSeqPtr) upp->update_bsp->seq_ext;
       dspnew != NULL; 
       dspnew = dspnew->next)
  {
    slip = (SeqLitPtr) (dspnew->data.ptrvalue);
    SeqLitFree (slip);
  }
  upp->update_bsp->seq_ext = ValNodeFree (upp->update_bsp->seq_ext);
  upp->update_bsp->seq_ext = new_list;
  upp->update_bsp->length = ualp->old5 + ualp->newa + ualp->old3;
  return TRUE;  
}

static Boolean UpdateSequencePatch 
(UpdatePairPtr             upp,
 UpdateAlignmentLengthsPtr ualp)

{
  Boolean rval = FALSE;
  if (upp == NULL || upp->orig_bsp == NULL || upp->update_bsp == NULL
      || upp->salp == NULL)
  {
    return FALSE;
  }

  if (RawSequencePatchOk (upp->orig_bsp, upp->update_bsp))
  {
    rval = UpdateSequencePatchRaw (upp, ualp);
  }
  else if (DeltaSequencePatchOk (upp->orig_bsp, upp->update_bsp))
  {
    rval = UpdateSequencePatchDelta (upp, ualp);
  }
  
  if (!rval)
  {
    return rval;
  }

  /* adjust alignment and reindex */
  rval = ShiftAlignmentForUpdate (upp->salp, ualp, eSequenceUpdatePatch);
  if (rval)
  {
    /* then finish by replacing with new sequence */
    ReplaceOneSequence (upp->salp, upp->orig_bsp, upp->update_bsp);
  }

  return rval;
}

static Boolean DoSequencesHaveDifferentOrganisms (BioseqPtr bsp1, BioseqPtr bsp2)
{
  OrgRefPtr         orp1 = NULL, orp2 = NULL; 
  SeqDescrPtr       sdp;
  SeqMgrDescContext dcontext;
  BioSourcePtr      biop1, biop2;
  Boolean           rval = FALSE;

  sdp = SeqMgrGetNextDescriptor (bsp1, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL) {
    biop1 = (BioSourcePtr) sdp->data.ptrvalue;
    if (biop1 != NULL) {
      orp1 = biop1->org;
    }
  }
  sdp = SeqMgrGetNextDescriptor (bsp2, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL) {
    biop2 = (BioSourcePtr) sdp->data.ptrvalue;
    if (biop2 != NULL) 
    {
      orp2 = biop2->org;
    }
  }
  if (orp1 != NULL && orp2 != NULL
      && StringICmp (orp1->taxname, orp2->taxname) != 0) 
  {
    rval = TRUE;
  }
  return rval;
}

/*This function removes RefTrackDescriptors from the specified Bioseq */
static void RemoveRefTrackDescriptors (BioseqPtr bsp)
{
	UserObjectPtr uop;
  ObjectIdPtr		oip;
  SeqDescrPtr   desc, prev_desc = NULL, next_desc;
	
	for (desc=bsp->descr; desc != NULL; desc=next_desc) 
	{
	  next_desc = desc->next;
		if (desc->choice != Seq_descr_user || desc->data.ptrvalue == NULL) 
		{
		  prev_desc = desc;
		  continue;
		}
	  uop = desc->data.ptrvalue;
	  if ((oip = uop->type) == NULL || StringCmp(oip->str, "RefGeneTracking") != 0)
	  {
	    prev_desc = desc;
	  }
	  else
	  {
	    if (prev_desc == NULL)
	    {
	      bsp->descr = desc->next;
	    }
	    else
	    {
	      prev_desc->next = desc->next;
	    }
	    desc->next = NULL;
	    SeqDescrFree (desc);
	  }
	}
}

static void 
ImportFeatureProduct 
(SeqFeatPtr  dup, 
 Boolean     keepProteinIDs, 
 SeqEntryPtr top,
 SeqIdPtr    nuc_sip)
{
  BioseqPtr   bsp, newbsp;
  SeqEntryPtr prdsep, newsep;
  
  if (dup == NULL || dup->product == NULL || nuc_sip == NULL || top == NULL)
  {
    return;
  }
  
  SeqEntrySetScope (NULL);
  bsp = BioseqFindFromSeqLoc (dup->product);
  if (bsp != NULL) {
    prdsep = SeqMgrGetSeqEntryForData (bsp);
    if (prdsep != NULL) {
      newsep = AsnIoMemCopy ((Pointer) prdsep,
                             (AsnReadFunc) SeqEntryAsnRead,
                             (AsnWriteFunc) SeqEntryAsnWrite);
      if (newsep != NULL) {
        if (IS_Bioseq (newsep)) {
          newbsp = (BioseqPtr) newsep->data.ptrvalue;
          if (newbsp != NULL) {
            /* we do not want to import reftrack descriptors with the products */
            RemoveRefTrackDescriptors (newbsp);
            if (! keepProteinIDs) {
              newbsp->id = SeqIdSetFree (newbsp->id);
              newbsp->id = MakeNewProteinSeqId (NULL, nuc_sip);
              newbsp->hist = SeqHistFree (newbsp->hist);
              VisitFeaturesOnBsp (newbsp, (Pointer) newbsp->id, CorrectFeatureSeqIds);
              SetSeqFeatProduct (dup, newbsp);
            }
            SeqMgrReplaceInBioseqIndex (newbsp);
          }
        }
        AddSeqEntryToSeqEntry (top, newsep, TRUE);
      }
    }
  }
}

static Boolean 
ImportFeaturesWithOffset 
(UpdatePairPtr             upp,
 UpdateOptionsPtr          uop,
 UpdateAlignmentLengthsPtr ualp,
  Int4         offset,
  SeqAnnotPtr PNTR sapp)

{
  CodeBreakPtr       cbp;
  SeqMgrFeatContext  context;
  CdRegionPtr        crp;
  SeqFeatPtr         dup, sfp, last = NULL;
  Uint2              entityID;
  Boolean            keepProteinIDs = FALSE;
  SeqEntryPtr        top;
  RnaRefPtr          rrp;
  SeqAnnotPtr        sap = NULL, saptmp;
  SeqDescrPtr        sdp;
  SeqIdPtr           sip;
  tRNAPtr            trp;

  if (upp == NULL || uop == NULL || ualp == NULL) return FALSE;

  SeqEntrySetScope (NULL);

  sfp = SeqMgrGetNextFeature (upp->update_bsp, NULL, 0, 0, &context);
  if (sfp == NULL) return FALSE;

  if (uop->indexer_opts != NULL 
      && uop->indexer_opts->keep_protein_ids
      && ! DoSequencesHaveDifferentOrganisms (upp->orig_bsp, upp->update_bsp))
  {
    keepProteinIDs = TRUE;
  } 

  entityID = ObjMgrGetEntityIDForPointer (upp->orig_bsp);
  top = GetBestTopParentForData (entityID, upp->orig_bsp);

  sdp = ExtractBioSourceAndPubs (top);

  sip = SeqIdFindBest (upp->orig_bsp->id, 0);

  for (; sfp != NULL; sfp = SeqMgrGetNextFeature (upp->update_bsp, sfp, 0, 0, &context))
  {
    if (uop->submitter_opts->sequence_update_type == eSequenceUpdatePatch
        && (context.right < ualp->new5 || context.left > ualp->new5 + ualp->newa))
    {
      /* this was a patch operation, and feature is outside patch area */
      continue;
    }
    if (sfp->data.choice == SEQFEAT_USER)
    {
      /* this is where we will continue if this is a RefTrack object */
    }
    dup = AsnIoMemCopy ((Pointer) sfp,
                        (AsnReadFunc) SeqFeatAsnRead,
                        (AsnWriteFunc) SeqFeatAsnWrite);

    /* if this is the first feature we have imported, create an annotation to hold it.
     * we put the features on a separate annotation so that we will be able to resolve
     * duplicates later.
     */
    if (last == NULL) 
    {
      sap = SeqAnnotNew ();
      if (upp->orig_bsp->annot == NULL) 
      {
        upp->orig_bsp->annot = sap;
      } 
      else 
      {
        for (saptmp = upp->orig_bsp->annot; saptmp->next != NULL; saptmp = saptmp->next) continue;
        saptmp->next = sap;
      }
      sap->type = 1;
      sap->data = (Pointer) dup;
    } else {
      last->next = dup;
    }
    last = dup;

    OffsetLocation (dup->location, offset, sip);
    switch (dup->data.choice) {
      case SEQFEAT_CDREGION :
        crp = (CdRegionPtr) dup->data.value.ptrvalue;
        if (crp != NULL) {
          for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
            OffsetLocation (cbp->loc, offset, sip);
          }
        }
        break;
      case SEQFEAT_RNA :
        rrp = (RnaRefPtr) dup->data.value.ptrvalue;
        if (rrp != NULL && rrp->ext.choice == 2) {
          trp = (tRNAPtr) rrp->ext.value.ptrvalue;
          if (trp != NULL && trp->anticodon != NULL) {
            OffsetLocation (trp->anticodon, offset, sip);
          }
        }
        break;
      default :
        break;
    }
    ImportFeatureProduct (dup, keepProteinIDs, top, sip);
  }

  ReplaceBioSourceAndPubs (top, sdp);

  if (sapp != NULL) {
    *sapp = sap;
  }

  return TRUE;
}

static Int4 GetMaxPositionForUpdate 
(BioseqPtr                 orig_bsp,
 UpdateOptionsPtr          uop,
 UpdateAlignmentLengthsPtr ualp)
{
  Int4 new_len = 0;
  
  if (uop == NULL || uop->submitter_opts == NULL)
  {
    return 0;
  }
  
  switch (uop->submitter_opts->sequence_update_type)
  {
    case eSequenceUpdateNoChange:
      if (orig_bsp != NULL)
      {
        new_len = orig_bsp->length;
      }
      break;
    case eSequenceUpdatePatch:
      if (ualp != NULL)
      {
        new_len = ualp->old5 + ualp->newa + ualp->old3;
      }
      break;
    case eSequenceUpdateReplace:
      if (ualp != NULL)
      {
        new_len = ualp->new5 + ualp->newa + ualp->new3;
      }
      break;
    case eSequenceUpdateExtend5:
      if (ualp != NULL)
      {
        new_len = ualp->new5 + ualp->newa + ualp->old3;
      }
      break;
    case eSequenceUpdateExtend3:
      if (ualp != NULL)
      {
        new_len = ualp->old5 + ualp->newa + ualp->new3;
      }
      break;
  }
  return new_len;  
}

static Boolean ImportFeaturesViaAlignment 
(UpdatePairPtr             upp,
 UpdateOptionsPtr          uop,
 UpdateAlignmentLengthsPtr ualp,
 SeqAnnotPtr PNTR          sapp)

{
  CodeBreakPtr       cbp, prevcbp, nextcbp;
  SeqMgrFeatContext  context;
  CdRegionPtr        crp;
  SeqFeatPtr         dup, sfp, last = NULL;
  Uint2              entityID;
  Int4               from, to, max_position;
  Boolean            keepProteinIDs = FALSE;
  SeqLocPtr          newloc;
  SeqEntryPtr        top;
  RnaRefPtr          rrp;
  SeqAnnotPtr        sap = NULL, saptmp;
  SeqDescrPtr        sdp;
  SeqIdPtr           sip;
  Boolean            split;
  tRNAPtr            trp;
  Boolean            partial5, partial3;

  if (upp == NULL || uop == NULL || ualp == NULL) return FALSE;

  SeqEntrySetScope (NULL);

  sfp = SeqMgrGetNextFeature (upp->update_bsp, NULL, 0, 0, &context);
  if (sfp == NULL) return FALSE;

  if (uop->indexer_opts != NULL && uop->indexer_opts->keep_protein_ids)
  {
    keepProteinIDs = TRUE;
  }

  entityID = ObjMgrGetEntityIDForPointer (upp->orig_bsp);
  top = GetBestTopParentForData (entityID, upp->orig_bsp);

  sdp = ExtractBioSourceAndPubs (top);

  sip = SeqIdFindBest (upp->orig_bsp->id, 0);

  from = ualp->new5;
  to = ualp->new5 + ualp->newa;
  
  max_position = GetMaxPositionForUpdate (upp->orig_bsp, uop, ualp);

  while (sfp != NULL) {

    if (context.right >= from && context.left <= to) {
      split = FALSE;
      newloc = GetPropagatedLocation (sfp->location, upp->update_bsp, upp->orig_bsp, 
                                      max_position, upp->salp);
      if (newloc != NULL) {
        CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
        SetSeqLocPartial (newloc, partial5, partial3);
        dup = AsnIoMemCopy ((Pointer) sfp,
                            (AsnReadFunc) SeqFeatAsnRead,
                            (AsnWriteFunc) SeqFeatAsnWrite);
                            
        SeqLocFree (dup->location);
        dup->location = newloc;
        if (split) {
          dup->partial = TRUE;
        }
        dup->partial |= partial5;
        dup->partial |= partial3;

        if (last == NULL) {
          sap = SeqAnnotNew ();
          if (upp->orig_bsp->annot == NULL) {
            upp->orig_bsp->annot = sap;
          } else {
            for (saptmp = upp->orig_bsp->annot; saptmp->next != NULL; saptmp = saptmp->next) continue;
            saptmp->next = sap;
          }
          sap->type = 1;
          sap->data = (Pointer) dup;
        } else {
          last->next = dup;
        }
        last = dup;

        switch (dup->data.choice) {
          case SEQFEAT_CDREGION :
            crp = (CdRegionPtr) dup->data.value.ptrvalue;
            if (crp != NULL) {
              prevcbp = NULL;
              for (cbp = crp->code_break; cbp != NULL; cbp = nextcbp) {
                nextcbp = cbp->next;
                newloc = GetPropagatedLocation (cbp->loc, upp->update_bsp, upp->orig_bsp, 
                                                max_position, upp->salp);
                SeqLocFree (cbp->loc);
                cbp->loc = newloc;
                if (cbp->loc == NULL) {
                  if (prevcbp != NULL) {
                    prevcbp->next = nextcbp;
                  } else {
                    crp->code_break = nextcbp;
                  }
                  cbp->next = NULL;
                  CodeBreakFree (cbp);
                } else {
                  prevcbp = cbp;
                }
              }
            }
            break;
          case SEQFEAT_RNA :
            rrp = (RnaRefPtr) dup->data.value.ptrvalue;
            if (rrp != NULL && rrp->ext.choice == 2) {
              trp = (tRNAPtr) rrp->ext.value.ptrvalue;
              if (trp != NULL && trp->anticodon != NULL) {
                newloc = GetPropagatedLocation (trp->anticodon, upp->update_bsp, upp->orig_bsp, 
                                                max_position, upp->salp);
                SeqLocFree (trp->anticodon);
                trp->anticodon = newloc;
              }
            }
            break;
          default :
            break;
        }
        ImportFeatureProduct (dup, keepProteinIDs, top, sip);
      }
    }

    sfp = SeqMgrGetNextFeature (upp->update_bsp, sfp, 0, 0, &context);
  }

  ReplaceBioSourceAndPubs (top, sdp);

  if (sapp != NULL) {
    *sapp = sap;
  }

  return TRUE;
}

static Boolean AreSequenceResiduesIdentical (BioseqPtr bsp1, BioseqPtr bsp2)
{
  SeqPortPtr    spp1, spp2;
  Uint1         seqcode;
  Int4          buf_len = 255;
  Char          buf1[255], buf2[255];
  Int4          ctr1, ctr2, offset;
  Boolean       rval;
  
  if (bsp1 == NULL && bsp2 == NULL)
  {
    return TRUE;
  }
  else if (bsp1 == NULL || bsp2 == NULL)
  {
    return FALSE;
  }
  else if (bsp1->length != bsp2->length)
  {
    return FALSE;
  }
  else if (ISA_na (bsp1->mol) && ! ISA_na (bsp2->mol))
  {
    return FALSE;
  }
  else if (!ISA_na (bsp1->mol) && ISA_na (bsp2->mol))
  {
    return FALSE;
  }

  if (ISA_na (bsp1->mol))
  {
    seqcode = Seq_code_iupacna;
  }
  else
  {
    seqcode = Seq_code_iupacaa;
  }

  
  spp1 = SeqPortNew (bsp1, 0, bsp1->length - 1, Seq_strand_plus, seqcode);
  spp2 = SeqPortNew (bsp2, 0, bsp2->length - 1, Seq_strand_plus, seqcode);
  
  ctr1 = SeqPortRead (spp1, (Uint1Ptr)buf1, buf_len - 1);
  ctr2 = SeqPortRead (spp2, (Uint1Ptr)buf2, buf_len - 1);
  buf1 [ctr1] = 0;
  buf2 [ctr2] = 0;
  offset = ctr1;

  while (ctr1 == ctr2 && StringCmp (buf1, buf2) == 0 && offset < bsp1->length)
  {
    ctr1 = SeqPortRead (spp1, (Uint1Ptr)buf1, buf_len - 1);
    ctr2 = SeqPortRead (spp2, (Uint1Ptr)buf2, buf_len - 1);
    buf1 [ctr1] = 0;
    buf2 [ctr2] = 0;
    offset += ctr1;
  }
  
  if (ctr1 != ctr2 || StringCmp (buf1, buf2) != 0 || offset < bsp1->length)
  {
    rval = FALSE;
  }
  else
  {
    rval = TRUE;
  }
  
  spp1 = SeqPortFree (spp1);
  spp2 = SeqPortFree (spp2);
  
  return rval;
}

/* This function examines the contexts for two features.
 * If the featdeftypes are the same and the locations are identical,
 * the features are considered to be duplicates.
 */
static Boolean 
AreFeaturesDuplicates 
(SeqMgrFeatContextPtr context1, 
 SeqMgrFeatContextPtr context2)
{
  Boolean ivalssame = FALSE;
  Int4    i, j;
  
  if (context1 == NULL || context2 == NULL)
  {
    return FALSE;
  }
  
  if (context1->left == context2->left &&
      context1->right == context2->right &&
      context1->featdeftype == context2->featdeftype) 
  {
    if (context1->strand == context2->strand ||
        context2->strand == Seq_strand_unknown ||
        context1->strand == Seq_strand_unknown) 
    {
      ivalssame = TRUE;
      if (context1->numivals != context2->numivals ||
          context1->ivals == NULL ||
          context2->ivals == NULL) 
      {
        ivalssame = FALSE;
      } 
      else 
      {
        for (i = 0, j = 0; i < context2->numivals; i++, j += 2) 
        {
          if (context1->ivals [j] != context2->ivals [j]) 
          {
            ivalssame = FALSE;
          }
          if (context1->ivals [j + 1] != context2->ivals [j + 1]) 
          {
            ivalssame = FALSE;
          }
        }
      }
    }
  }
  return ivalssame;
}

static void ResolveDuplicateUpdateFeats 
(BioseqPtr        bsp,
 UpdateOptionsPtr uop,
 SeqAnnotPtr      newfeat_sap)

{
  SeqMgrFeatContext  context, lastcontext;
  SeqFeatPtr         lastsfp = NULL, sfp;

  if (bsp == NULL 
      || newfeat_sap == NULL 
      || uop == NULL
      || uop->submitter_opts == NULL
      || uop->submitter_opts->feature_update_type == eFeatureUpdateAll) 
  {
    return;
  }

  SeqMgrIndexFeatures (0, (Pointer) bsp);

  lastsfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &context);
  if (lastsfp == NULL) return;

  MemCopy ((Pointer) &lastcontext, (Pointer) &context, sizeof (SeqMgrFeatContext));

  sfp = SeqMgrGetNextFeature (bsp, lastsfp, 0, 0, &context);
  if (sfp == NULL) return;

  while (sfp != NULL) 
  {
    if (context.sap != lastcontext.sap 
        && (context.sap == newfeat_sap || lastcontext.sap == newfeat_sap)
        && AreFeaturesDuplicates (&context, &lastcontext))
    {
      if (uop->submitter_opts->feature_update_type == eFeatureUpdateAllReplaceDups) { /* keep new */
        if (context.sap == newfeat_sap) {
          lastsfp->idx.deleteme = TRUE;
          MarkProductForDeletion (lastsfp->product);
        } else if (lastcontext.sap == newfeat_sap) {
          sfp->idx.deleteme = TRUE;
          MarkProductForDeletion (sfp->product);
        }
      } else if (uop->submitter_opts->feature_update_type == eFeatureUpdateAllExceptDups) { /* keep old */
        if (context.sap == newfeat_sap) {
          sfp->idx.deleteme = TRUE;
          MarkProductForDeletion (sfp->product);
        } else if (lastcontext.sap == newfeat_sap) {
          lastsfp->idx.deleteme = TRUE;
          MarkProductForDeletion (lastsfp->product);
        }

      } else if (uop->submitter_opts->feature_update_type == eFeatureUpdateAllMergeDups) { /* merge */
        if (context.sap == newfeat_sap) {
          FuseFeatures (sfp, lastsfp);
          lastsfp->idx.deleteme = TRUE;
          MarkProductForDeletion (lastsfp->product);
        } else if (lastcontext.sap == newfeat_sap) {
          FuseFeatures (lastsfp, sfp);
          sfp->idx.deleteme = TRUE;
          MarkProductForDeletion (sfp->product);
        }
      }
    }

    lastsfp = sfp;
    MemCopy ((Pointer) &lastcontext, (Pointer) &context, sizeof (SeqMgrFeatContext));

    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &context);
  }
}



static void 
UpdateProteinForOneUpdatedCodingRegion 
(SeqFeatPtr        sfp,
 ValNodePtr        transl_except_list,
 IndexerOptionsPtr iop,
 Uint2             entityID,
 FILE              *log_fp,
 BoolPtr           data_in_log)
{
  SeqLocPtr     new_product;
  Int4          transl_except_len = 0;
  Boolean       fix_products = TRUE;
  BioseqPtr     protBsp = NULL;
  SeqAlignPtr   salp = NULL;
  BioseqPtr     newbsp = NULL;
  Boolean       rval;

  if (sfp == NULL 
      || sfp->idx.subtype != FEATDEF_CDS 
      || sfp->idx.deleteme)
  {
    return;
  }
    
  protBsp = BioseqFindFromSeqLoc (sfp->product);
  if (protBsp == NULL)
  {
    return;
  }

  transl_except_len = GetOriginalTranslExceptLen (sfp, transl_except_list);
  rval = PrepareUpdateAlignmentForProtein (sfp,
                                           protBsp,
                                           entityID,
                                           log_fp,
                                           iop == NULL ? FALSE : iop->truncate_proteins,
                                           iop == NULL ? FALSE : iop->extend_proteins5,
                                           iop == NULL ? FALSE : iop->extend_proteins3,
                                           iop == NULL ? FALSE : iop->correct_cds_genes,
                                           transl_except_len,
                                           data_in_log,
                                           &salp,
                                           &newbsp);
  if (!rval) return;
  
  if (protBsp->idx.deleteme)
  {
    fix_products = FALSE;
  }

  ReplaceOneSequence (salp, protBsp, newbsp);
  
  if (fix_products)
  {
    if (sfp->product->choice != SEQLOC_WHOLE) {
      new_product = SeqLocWholeNew (protBsp);
      if (new_product == NULL) return;
      SeqLocFree (sfp->product);
      sfp->product = new_product;
    }
    newbsp = BioseqFree (newbsp);
  }
}

static void 
UpdateProteinsForUpdatedCodingRegions 
(BioseqPtr         orig_bsp,
 ValNodePtr        transl_except_list,
 IndexerOptionsPtr iop,
 Uint2             entityID,
 FILE              *log_fp,
 BoolPtr           data_in_log)
{
  SeqFeatPtr        sfp;
  SeqMgrFeatContext context;
  
  sfp = SeqMgrGetNextFeature (orig_bsp, NULL, SEQFEAT_CDREGION, 0, &context);
  while (sfp != NULL)
  {
    UpdateProteinForOneUpdatedCodingRegion (sfp, transl_except_list, iop, 
                                            entityID, log_fp, data_in_log);    
    sfp = SeqMgrGetNextFeature (orig_bsp, sfp, SEQFEAT_CDREGION, 0, &context);
  }
}

static Boolean 
IsSequenceUpdateChoiceAllowed 
(ESequenceUpdateType       action,
 UpdatePairPtr             upp,
 UpdateAlignmentLengthsPtr ualp,
 Boolean                   is_indexer,
 Boolean                   ignore_alignment)
{

  if (upp == NULL || upp->orig_bsp == NULL || upp->update_bsp == NULL)
  {
    return FALSE;
  }
  else if (action == eSequenceUpdateNoChange)
  {
    if (ignore_alignment)
    {
      return FALSE;
    }
    else
    {
      return TRUE;
    }
  }
  else if (action != eSequenceUpdateNoChange 
           && !is_indexer
           && (upp->orig_bsp->repr != Seq_repr_raw || upp->update_bsp->repr != Seq_repr_raw))
  {
    /* If either sequence is not raw and not indexer version, do not allow sequence update */
    return FALSE;
  }  
  else if (action == eSequenceUpdatePatch
           && (upp->salp == NULL || ignore_alignment 
             || (! RawSequencePatchOk (upp->orig_bsp, upp->update_bsp)
                 && !DeltaSequencePatchOk (upp->orig_bsp, upp->update_bsp))))
  {
    /* If no alignment or patch not available then disable the patch button */
    return FALSE;
  }
  else if (ignore_alignment 
           && (action == eSequenceUpdateExtend3 
               || action == eSequenceUpdateExtend5))
  {
    return TRUE;
  }
  else if (ualp == NULL)
  {
    return FALSE;
  }
  else if (action == eSequenceUpdateExtend3 
           && (ualp->new3 < ualp->old3 || ualp->new3 == 0))
  {
    return FALSE;
  }
  else if (action == eSequenceUpdateExtend5 
           && (ualp->new5 < ualp->old5 || ualp->new5 == 0))
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}

static void 
ReportBadUpdateType 
(ESequenceUpdateType action,
 UpdatePairPtr       upp,
 FILE                *log_fp,
 BoolPtr             data_in_log)
{
  Char id_txt [MAX_ID_LEN];
  
  if (upp == NULL 
      || upp->orig_bsp == NULL 
      || log_fp == NULL 
      || data_in_log == NULL)
  {
    return;
  }
  
  SeqIdWrite (SeqIdFindBest (upp->orig_bsp->id, SEQID_GENBANK), id_txt, 
              PRINTID_REPORT, sizeof (id_txt) - 1);
  
  switch (action)
  {
    case eSequenceUpdateNoChange:
      /* do nothing */
      break;
    case eSequenceUpdateReplace:
      fprintf (log_fp, "Unable to replace sequence for %s\n", id_txt);
      *data_in_log = TRUE;
      break;
    case eSequenceUpdatePatch:
      fprintf (log_fp, "Unable to patch sequence for %s\n", id_txt);
      *data_in_log = TRUE;
      break;
    case eSequenceUpdateExtend5:
      fprintf (log_fp, "Unable to extend sequence on 5' end for %s\n", id_txt);
      *data_in_log = TRUE;
      break;
    case eSequenceUpdateExtend3:
      fprintf (log_fp, "Unable to extend sequence on 3' end for %s\n", id_txt);
      *data_in_log = TRUE;
      break;
    default:
      fprintf (log_fp, "Unknown update operation for %s\n", id_txt);
      *data_in_log = TRUE;
      break;
  }
}

static Boolean 
UpdateOrExtendOneSequence 
(UpdatePairPtr             upp,
 UpdateOptionsPtr          uop,
 UpdateAlignmentLengthsPtr ualp,
 Uint2                     entityID,
 FILE                      *log_fp,
 BoolPtr                   data_in_log)
{
  ValNodePtr   prot_feat_list = NULL;
  ValNodePtr   transl_except_list = NULL;
  Boolean      sequence_update_successful = FALSE;
  Boolean      feature_update_successful = FALSE;
  SeqEntryPtr  sep;
  Boolean      rval = FALSE;
  Int4         feature_update_offset;
  SeqAnnotPtr  newfeat_sap = NULL;
  Int4         orig_seq_len;
  Uint2        update_entityID;
  Char         id_txt [128];

  if (upp == NULL || upp->orig_bsp == NULL || upp->update_bsp == NULL
      || uop == NULL || uop->submitter_opts == NULL)
  {
    return FALSE;
  }
  
  update_entityID = upp->update_bsp->idx.entityID; 

  if (!IsSequenceUpdateChoiceAllowed (uop->submitter_opts->sequence_update_type,
                                      upp, ualp, indexerVersion,
                                      uop->submitter_opts->ignore_alignment))
  {
    ReportBadUpdateType (uop->submitter_opts->sequence_update_type,
                         upp, log_fp, data_in_log);
    return FALSE;
  }
  
  if (uop->submitter_opts->sequence_update_type == eSequenceUpdateReplace
      && upp->salp == NULL)
  {
    SeqIdWrite (SeqIdFindBest (upp->orig_bsp->id, SEQID_GENBANK), id_txt,
                PRINTID_REPORT, sizeof (id_txt) - 1);
    
    if (ANS_NO == Message (MSG_YN, 
                           "There is no alignment for %s.  Are you sure that you want to replace this sequence?",
                           id_txt))
    {
      return FALSE;
    }        
  }
  
  /* remove features */
  switch (uop->submitter_opts->feature_remove_type)
  {
    case eFeatureRemoveAll:
      RemoveOldFeats (upp->orig_bsp);
      DeleteMarkedObjects (entityID, 0, NULL);
      break;
    case eFeatureRemoveAligned:
      RemoveFeatsInAlignedRegion (upp->orig_bsp, ualp);
      DeleteMarkedObjects (entityID, 0, NULL);
      break;
    case eFeatureRemoveNotAligned:
      RemoveFeatsNotInAlignedRegion (upp->orig_bsp, ualp);
      DeleteMarkedObjects (entityID, 0, NULL);
      break;
    default:
      break;
  }
  
  if (ISA_na (upp->orig_bsp->mol) && uop->indexer_opts != NULL && uop->indexer_opts->update_proteins)
  {
    prot_feat_list = FindProductProtRefs (upp->orig_bsp);
    transl_except_list = FindTranslExceptCDSs (upp->orig_bsp);
  }

  /* get length of original sequence before updating */
  orig_seq_len = upp->orig_bsp->length;

  /* update sequence */  
  switch (uop->submitter_opts->sequence_update_type)
  {
    case eSequenceUpdateReplace:
      ReplaceOneSequence (upp->salp, upp->orig_bsp, upp->update_bsp);
      sequence_update_successful = TRUE;
      rval = sequence_update_successful;
      break;
    case eSequenceUpdatePatch:
      sequence_update_successful = UpdateSequencePatch (upp, ualp); 
      rval = sequence_update_successful;
      break;
    case eSequenceUpdateExtend5:
      sequence_update_successful = UpdateSequenceExtend5Prime (upp, ualp);    
      rval = sequence_update_successful;
      break;
    case eSequenceUpdateExtend3:
      sequence_update_successful = UpdateSequenceExtend3Prime (upp, ualp);    
      rval = sequence_update_successful;
      break;
    case eSequenceUpdateNoChange:
      rval = TRUE;
      break;    
  }

  if (log_fp != NULL && data_in_log != NULL && upp->revcomp)
  {
     SeqIdWrite (SeqIdFindBest (upp->orig_bsp->id, SEQID_GENBANK), id_txt,
                 PRINTID_REPORT, sizeof (id_txt) - 1);
     fprintf (log_fp, "Reverse complemented %s\n", id_txt);
     *data_in_log = TRUE;
  }
  
  /* update features */
  if (uop->submitter_opts->feature_update_type != eFeatureUpdateNoChange)
  {
    if (! SeqMgrFeaturesAreIndexed (update_entityID))
      SeqMgrIndexFeatures (update_entityID, NULL);
  
    if (sequence_update_successful)
    {
      feature_update_offset = 0;
      if (uop->submitter_opts->sequence_update_type == eSequenceUpdateExtend3
          || uop->submitter_opts->sequence_update_type == eSequenceUpdatePatch)
      {
        if (upp->salp == NULL)
        {
          feature_update_offset = orig_seq_len;
        }
        else
        {
          feature_update_offset = ualp->old5 - ualp->new5;
        }
      }
      feature_update_successful = ImportFeaturesWithOffset (upp, uop, ualp,
                                                            feature_update_offset,
                                                            &newfeat_sap);
      rval &= feature_update_successful;
    }
    else
    {
      feature_update_successful = ImportFeaturesViaAlignment (upp, uop, ualp, &newfeat_sap);
    }
  }

  if (uop->submitter_opts->sequence_update_type == eSequenceUpdateReplace
      && upp->revcomp 
      && sequence_update_successful)
  {
    ReverseComplementBioseqAndFeats (upp->orig_bsp, entityID);
  }


  if (uop->indexer_opts != NULL && uop->indexer_opts->add_cit_subs
      && (feature_update_successful
          || (sequence_update_successful 
              && ! AreSequenceResiduesIdentical(upp->orig_bsp, upp->update_bsp))))
     
  {
    AddCitSubToUpdatedSequence (upp->orig_bsp, entityID);
  }
  
  /* update proteins for coding regions on the updated sequence */
  if (sequence_update_successful
      && uop->indexer_opts != NULL
      && uop->indexer_opts->update_proteins)
  {
    SeqMgrClearFeatureIndexes (entityID, upp->orig_bsp);
    SeqMgrIndexFeatures (entityID, NULL);
    UpdateProteinsForUpdatedCodingRegions (upp->orig_bsp, transl_except_list,
                                           uop->indexer_opts, entityID,
                                           log_fp, data_in_log);
    FixProtRefPtrs (prot_feat_list);        
  }
  
  if (uop->indexer_opts != NULL && uop->indexer_opts->update_quality_scores
      && uop->submitter_opts->sequence_update_type == eSequenceUpdateReplace)
  {
    ReplaceQualityScores (upp->orig_bsp, upp->update_bsp, upp->salp, 
                          log_fp, data_in_log);
  }
 
  prot_feat_list = ValNodeFree (prot_feat_list);
  transl_except_list = ValNodeFree (transl_except_list);
  if (feature_update_successful
      && uop->submitter_opts->feature_update_type > eFeatureUpdateNoChange
      && uop->submitter_opts->feature_remove_type != eFeatureRemoveAll)
  {
    /* resolve features unless the policy was to remove all the old ones */
  
    sep = GetTopSeqEntryForEntityID (entityID);
    /* need to set scope to make sure we mark the right bioseq for deletion */
    SeqEntrySetScope (sep);
    ResolveDuplicateUpdateFeats (upp->orig_bsp, uop, newfeat_sap);
    SeqEntrySetScope (NULL);
    DeleteMarkedObjects (entityID, 0, NULL);
    SeqMgrClearFeatureIndexes (entityID, NULL);
  }

  /* Note - we do not remove the update sequence here - that will happen when
   * the UpdateSequence form is removed.
   */
  upp->update_bsp = NULL;

  return rval;
}

static void 
CalculateUpdateAlignmentLengths 
(SeqAlignPtr               salp,
 BioseqPtr                 orig_bsp,
 BioseqPtr                 update_bsp,
 UpdateAlignmentLengthsPtr ualp)
{
  Int4         stopold, startold, lenold, stopnew, startnew, lennew;
  Int4         aln_length;

  if (ualp == NULL || update_bsp == NULL || orig_bsp == NULL) return;
  
  MemSet (ualp, 0, sizeof (UpdateAlignmentLengthsData));
  
  if (salp == NULL)
  {
    return;
  }
    
  ualp->aln_length = AlnMgr2GetAlnLength (salp, FALSE);
  AlnMgr2GetNthSeqRangeInSA (salp, 1, &startold, &stopold);
  AlnMgr2GetNthSeqRangeInSA (salp, 2, &startnew, &stopnew);
  lenold = orig_bsp->length;
  lennew = update_bsp->length;

  ualp->old5 = startold;
  ualp->old3 = lenold - stopold - 1;
  ualp->olda = stopold - startold + 1;

  ualp->new5 = startnew;
  ualp->new3 = lennew - stopnew - 1;
  ualp->newa = stopnew - startnew + 1;

#if 0  
  ualp->startmax = MAX (startold, startnew);
  ualp->stopmax = MAX (aln_length + lenold - stopold, aln_length + lennew - stopnew);

  ualp->strandold = AlnMgr2GetNthStrand (sap, 1);
  ualp->strandnew = AlnMgr2GetNthStrand (sap, 2);
#endif  

  /* calculate logarithmic scale for length */
  ualp->log10_aln_length = 1;
  aln_length = ualp->aln_length;
  while (aln_length >= 10) {
    aln_length /= 10;
    (ualp->log10_aln_length)++;
  }
  
  /* calculate alignment recombination lengths */
  ualp->recomb1 = -1;
  ualp->recomb2 = -1;
  if (ualp->new5 > ualp->old5 && ualp->new3 < ualp->old3) {
    ualp->recomb2 = ualp->aln_length;
  }

  /* Extend 3' */

  else if (ualp->new5 < ualp->old5 && ualp->new3 > ualp->old3) {
    ualp->recomb1 = 0;
  }

  /* Replace */

  else {
    ualp->recomb1 = 0;
    ualp->recomb2 = ualp->aln_length;
  }
  

}

/* The following group of functions is used to generate the pictures for the
 * alignment viewer for the update sequence dialog.
 */

/* This function draws one rectangle for the new sequence and one rectangle for
 * the old sequence.  The aligned overlap is filled in with black, the unaligned
 * portions are drawn with outlines only.
 */
static SegmenT MakeAlignmentOverviewPicture (UpdateAlignmentLengthsPtr ualp)

{
  SegmenT   pict;
  Char      str [96];
  Int4      top, bottom;

  pict = CreatePicture ();
  if (ualp == NULL) return pict;
  
  top = 0;
  bottom = top - 10;

  DrawAlignBlock (pict, top, bottom, bottom, LOWER_CENTER, ualp->old5, ualp->olda, ualp->old3, ualp->aln_length);

  sprintf (str, "%ld", (long) ualp->aln_length);
  AddLabel (pict, ualp->aln_length / 2, 10, str, SMALL_TEXT, 5, MIDDLE_CENTER, 0);


  top = 30;
  bottom = top - 10;

  DrawAlignBlock (pict, top, bottom, top, UPPER_CENTER, ualp->new5, ualp->newa, ualp->new3, ualp->aln_length);

  return pict;
}

/* This function shows how the new sequence will be created when the sequence is
 * extended and the alignment is ignored.
 */
static SegmenT MakeExtensionPicture (Int4 len_old, Int4 len_new, Boolean extend5)

{
  SegmenT   pict;
  Char      str [96];
  Int4      top, bottom, label_loc;

  pict = CreatePicture ();
  
  top = 0;
  bottom = top - 10;
  
  label_loc = -(top + bottom) / 2;

  if (extend5)
  {
    sprintf (str, "New (%ld)", (long) len_new);
    AddLabel (pict, -len_new / 2, label_loc, str, SMALL_TEXT, 5, MIDDLE_LEFT, 0);
    AddRectangle (pict, -len_new, top, 0, bottom, NO_ARROW, FALSE, 0);
    sprintf (str, "Old (%ld)", (long) len_old);
    AddLabel (pict, len_old / 2, label_loc - 30, str, SMALL_TEXT, 5, MIDDLE_RIGHT, 0);
    AddRectangle (pict,  0, top, len_old, bottom, NO_ARROW, FALSE, 0);
  }
  else
  {
    sprintf (str, "Old (%ld)", (long) len_old);
    AddLabel (pict, -len_old / 2, label_loc - 30, str, SMALL_TEXT, 5, MIDDLE_LEFT, 0);
    AddRectangle (pict, -len_old, top, 0, bottom, NO_ARROW, FALSE, 0);
    sprintf (str, "New (%ld)", (long) len_new);
    AddLabel (pict, len_new / 2, label_loc, str, SMALL_TEXT, 5, MIDDLE_RIGHT, 0);
    AddRectangle (pict,  0, top, len_new, bottom, NO_ARROW, FALSE, 0);
  }

  return pict;  
}

static void DrawUpdateAlignmentBits (
  SegmenT         pict,
  Int4            top,
  Int4            bottom,
  Int4            row,
  Int4            pos1,
  Int4            pos2,
  SeqAlignPtr     sap,
  ValNodePtr PNTR indels /* pointer to ValNode list of integers which are the locations
                          * of insertions and deletions.
                          */
)

{
  AlnMsg2Ptr  amp;
  Int4       len, start, stop, from, to;
  Char       str [96];
  Boolean    wasgap;

  amp = AlnMsgNew2 ();
  if (amp == NULL) return;

  amp->from_aln = 0;
  amp->to_aln = -1;
  amp->row_num = row;

  start = 0;
  stop = 0;
  from = 0;
  to = 0;
  wasgap = FALSE;

  while (AlnMgr2GetNextAlnBit (sap, amp)) {
    len = amp->to_row - amp->from_row + 1;
    stop = start + len;
    if (amp->type == AM_GAP) {
      if (wasgap) {
        to = stop;
      } else {
        AddRectangle (pict, from, top, to, bottom, NO_ARROW, FALSE, 0);
        wasgap = TRUE;
        from = start;
        to = stop;
      }
    } else {
      if (wasgap) {

        /* record for accurate scrolling to text view */
        ValNodeAddInt (indels, 0, from);

        AddLine (pict, from, (top + bottom) / 2, to, (top + bottom) / 2, FALSE, 0);
        wasgap = FALSE;
        from = start;
        to = stop;
      } else {
        to = stop;
      }
    }
    start += len;
  }

  if (to > from) {
    if (wasgap) {

      /* record for accurate scrolling to text view */
      ValNodeAddInt (indels, 0, from);

      AddLine (pict, from, (top + bottom) / 2, to, (top + bottom) / 2, FALSE, 0);
    } else {
      AddRectangle (pict, from, top, to, bottom, NO_ARROW, FALSE, 0);
    }
  }

  AlnMsgFree2 (amp);

  sprintf (str, "%ld", (long) pos1);
  AddLabel (pict, 0, (top + bottom) / 2, str, SMALL_TEXT, 5, MIDDLE_LEFT, 0);

  sprintf (str, "%ld", (long) pos2);
  AddLabel (pict, to, (top + bottom) / 2, str, SMALL_TEXT, 5, MIDDLE_RIGHT, 0);
}

static void DrawUpdateAlignmentDiffs 
(SegmenT         pict,
 Int4            top,
 Int4            bottom,
 SeqAlignPtr     salp,
 BioseqPtr       orig_bsp,
 BioseqPtr       update_bsp,
 ValNodePtr PNTR mismatches)

{
  AlnMsg2Ptr  amp1, amp2;
  SegmenT    seg;
  Int4       len1, len2, i;
  Int4       seg_i, seg_n, seg_start, seg_stop;
  CharPtr      seq1, seq2;
  Uint2        entityID;
  SeqEntryPtr  sep;

  if (salp == NULL || orig_bsp == NULL || update_bsp == NULL 
      || pict == NULL || mismatches == NULL)
  {
    return;
  }

  entityID = ObjMgrGetEntityIDForPointer (orig_bsp);
  sep = GetTopSeqEntryForEntityID (entityID);
  SeqEntrySetScope (sep);
  seq1 = GetSequenceByBsp (orig_bsp);
  SeqEntrySetScope (NULL);

  entityID = ObjMgrGetEntityIDForPointer (update_bsp);
  sep = GetTopSeqEntryForEntityID (entityID);
  SeqEntrySetScope (sep);
  seq2 = GetSequenceByBsp (update_bsp);
  SeqEntrySetScope (NULL);

  if (seq1 == NULL || seq2 == NULL) 
  {
    seq1 = MemFree (seq1);
    seq2 = MemFree (seq2);
    return;
  }
  len1 = StringLen (seq1);
  len2 = StringLen (seq2);

  seg = CreateSegment (pict, 0, 0);
  AddAttribute (seg, COLOR_ATT, RED_COLOR, 0, 0, 0, 0);

  seg_n = AlnMgr2GetNumSegs(salp);
  for (seg_i = 1; seg_i<=seg_n; seg_i++) {
    AlnMgr2GetNthSegmentRange(salp, seg_i, &seg_start, &seg_stop);

    amp1 = AlnMsgNew2 ();
    amp2 = AlnMsgNew2 ();
    if (amp1 == NULL || amp2 == NULL) return;

    amp1->from_aln = seg_start;
    amp1->to_aln = seg_stop;
    amp1->row_num = 1;

    amp2->from_aln = seg_start;
    amp2->to_aln = seg_stop;
    amp2->row_num = 2;

    AlnMgr2GetNextAlnBit (salp, amp1);
    AlnMgr2GetNextAlnBit (salp, amp2);

    if (amp1->to_row - amp1->from_row == amp2->to_row - amp2->from_row &&
        amp1->type == AM_SEQ && amp2->type == AM_SEQ) {
      for (i=0; i<seg_stop-seg_start+1; i++) {
        if (seq1[amp1->from_row+i] != seq2[amp2->from_row+i]) {

          /* record for accurate scrolling to text view */
          ValNodeAddInt (mismatches, 0, i);

          AddLine (seg, seg_start+i, top, seg_start+i, bottom, FALSE, 0);
        }
      }
    }

    AlnMsgFree2 (amp1);
    AlnMsgFree2 (amp2);
  }
}

/* This function draws just the aligned region.
 */
static SegmenT MakeAlignmentDetailsPicture 
(SeqAlignPtr               salp,
 BioseqPtr                 orig_bsp,
 BioseqPtr                 update_bsp, 
 ValNodePtr PNTR           indels,
 ValNodePtr PNTR           mismatches,
 UpdateAlignmentLengthsPtr ualp,
 Boolean                   revcomp)

{
  Int4     aln_length;
  SegmenT  pict;
  Int4     top, bottom;

  pict = CreatePicture ();
  if (salp == NULL || ualp == NULL || indels == NULL) return pict;

  aln_length = ualp->aln_length;

  top = 0;
  bottom = top - 10;

  DrawUpdateAlignmentBits (pict, top, bottom, 1, ualp->old5 + 1, 
                           ualp->old5 + ualp->olda, salp, indels);

  top = 30;
  bottom = top - 10;

  if (revcomp) {
    DrawUpdateAlignmentBits (pict, top, bottom, 2, ualp->new3 + ualp->newa,
                             ualp->new3 + 1, salp, indels);
  } else {
    DrawUpdateAlignmentBits (pict, top, bottom, 2, ualp->new5 + 1, 
                             ualp->new5 + ualp->newa, salp, indels);
  }

  top = 15;
  bottom = top - 10;

  DrawUpdateAlignmentDiffs (pict, top, bottom, salp, orig_bsp, update_bsp, mismatches);

  return pict;
}



static Int4 CalculateBestUpdateAlignmentViewerScale (VieweR vwr, SegmenT pict)

{
  BoxInfo  box;
  Int2     i;
  Int4     max, worldwid, portwid;
  RecT     r;
  Int4     scaleX, oldscaleX;
  Int4     wid;

  ObjectRect (vwr, &r);
  InsetRect (&r, 4, 4);
  wid = (Int4) (r.right - r.left + 1);

  SegmentBox (pict, &box);
  oldscaleX = (box.right - box.left + wid - 1) / wid;
  RecalculateSegment (pict, oldscaleX, 1);
  SegmentBox (pict, &box);
  portwid = wid * oldscaleX;
  worldwid = box.right - box.left + 20 * oldscaleX + 1;
  max = MAX (worldwid, portwid);
  scaleX = (max + wid - 1) / wid;
  i = 0;
  while (i < 10 && (scaleX > oldscaleX || portwid < worldwid)) {
    oldscaleX = scaleX;
    RecalculateSegment (pict, oldscaleX, 1);
    SegmentBox (pict, &box);
    portwid = wid * oldscaleX;
    worldwid = box.right - box.left + 20 * oldscaleX + 1;
    max = MAX (worldwid, portwid);
    scaleX = (max + wid - 1) / wid;
    i++;
  }

  return scaleX;
}

static CharPtr MakeUpdateAlignmentSequenceString 
( SeqAlignPtr sap,
  Int4        aln_length,
  Int4        row,
  CharPtr     seq)

{
  CharPtr    aln;
  AlnMsg2Ptr  amp;
  Int4       len, lens, start, stop, from, to, i, j;

  if (sap == NULL || seq == NULL || aln_length < 1) return NULL;
  lens = StringLen (seq);

  aln = (CharPtr) MemNew (sizeof (Char) * (aln_length + 2));
  if (aln == NULL) return NULL;
  MemSet ((Pointer) aln, '-', aln_length);

  amp = AlnMsgNew2 ();
  if (amp == NULL) return aln;

  amp->from_aln = 0;
  amp->to_aln = -1;
  amp->row_num = row;

  start = 0;
  stop = 0;
  from = 0;
  to = 0;

  while (AlnMgr2GetNextAlnBit (sap, amp)) {
    len = amp->to_row - amp->from_row + 1;
    stop = start + len;

    if (amp->type == AM_SEQ) {
      for (i = start, j = amp->from_row; i < stop && j < lens; i++, j++) {
        aln [i] = seq [j];
      }
    }
    start += len;
  }

  AlnMsgFree2 (amp);

  return aln;
}

/* The alignment letters panel draws the individual letters of the sequence, showing
 * insertions, deletions, and mismatches.
 */
typedef struct alignmentletterspanel
{
  DIALOG_MESSAGE_BLOCK
  PaneL letters;

  Int4    lineheight;       /* height of lines in the letters panel */
  Int4    charwidth;        /* width of characters in the letters panel */
  Int4    maxchars;         /* maximum number of characters that can be 
                             * displayed in the letters panel. 
                             */
  
  UpdateAlignmentLengthsData uald; /* structure that holds the length of the alignment,
                                    * the lengths of the 3' and 5' overlaps, etc.
                                    */
  CharPtr            aln1;       /* holds the string representation of the alignment
                                  * for the original sequence 
                                  */
  CharPtr            aln2;       /* holds the string representation of the alignment
                                  * for the update sequence
                                  */              
                                  
  SeqAlignPtr        salp;       /* alignment used for update */
                                 /* This alignment should not be freed by this panel.
                                  */
  Boolean            revcomp;
} AlignmentLettersPanelData, PNTR AlignmentLettersPanelPtr;

static void 
PaintAlignmentString 
(CharPtr aln_str,
 Int4    char_offset,
 Int4    char_width,
 Int4    left_pos,
 Int4    top_pos,
 Int4    right_pos,
 Int4    maxchars,
 Int4    aln_length)
{
  Int2 i;
  Int4 j;
  Int4 pos_offset;
  
  if (aln_str != NULL)
  {
    MoveTo (left_pos, top_pos);
    pos_offset = left_pos;
    for (i = 0, j = char_offset, pos_offset = left_pos; 
         i < maxchars && j < aln_length && pos_offset + char_width <= right_pos;
         i++, j++, pos_offset += char_width) {
      PaintChar (aln_str [j]);
    }
  }
  
}


static void UpdateAlignmentLettersPanelOnDraw (PaneL pnl)

{
  Char        ch1, ch2;
  Int2        i, k, q, left, top, bottom, arrowwidth;
  size_t      len;
  Int4        offset, j, pos, realpos;
  RecT        r, x;
  BaR         sb;
  Char        str [32];
  AlignmentLettersPanelPtr  dlg;

  dlg = (AlignmentLettersPanelPtr) GetObjectExtra (pnl);
  if (dlg == NULL) return;

  ObjectRect (pnl, &r);
  InsetRect (&r, 4, 4);

  sb = GetSlateHScrollBar ((SlatE) pnl);
  offset = GetBarValue (sb);

  SelectFont (SetSmallFont ());

  /* draw top (new) letters */

  PaintAlignmentString (dlg->aln2, offset, dlg->charwidth,
                        r.left, r.top + 8 + 3 * dlg->lineheight,
                        r.right,
                        dlg->maxchars,
                        dlg->uald.aln_length);

  /* draw bottom (old) letters */

  PaintAlignmentString (dlg->aln1, offset, dlg->charwidth,
                        r.left, r.top + 8 + 5 * dlg->lineheight,
                        r.right,
                        dlg->maxchars,
                        dlg->uald.aln_length);

  /* draw recombination arrows */

  arrowwidth = MIN (6, dlg->charwidth);
  if (dlg->uald.recomb1 >= offset && dlg->uald.recomb1 <= offset + dlg->maxchars) {
    left = r.left + dlg->charwidth * (dlg->uald.recomb1 - offset);
    LoadRect (&x, left, r.top, left + arrowwidth, r.top + 6);
    CopyBits (&x, leftTriFillSym);
  }

  if (dlg->uald.recomb2 >= offset && dlg->uald.recomb2 <= offset + dlg->maxchars) {
    left = r.left + dlg->charwidth * (dlg->uald.recomb2 - offset - 1);
    LoadRect (&x, left, r.top, left + arrowwidth, r.top + 6);
    CopyBits (&x, rightTriFillSym);
  }

  if (dlg->aln1 == NULL || dlg->aln2 == NULL) 
  {
  	return;
  }
  /* draw red mismatch lines */

  Red ();
  top = r.top + 8 + 4 * dlg->lineheight - Ascent ();
  bottom = top + dlg->lineheight - 2;

  for (i = 0, j = offset; i < dlg->maxchars && j < dlg->uald.aln_length; i++, j++) {
    ch1 = dlg->aln1 [j];
    ch2 = dlg->aln2 [j];
    if (ch1 == ch2) {
    } else if (ch1 == '-' || ch2 == '-') {
    } else {
      left = r.left + i * dlg->charwidth + dlg->charwidth / 2 - 1;
      MoveTo (left, top);
      LineTo (left, bottom);
    }
  }
  Black ();

  /* draw top (new) tick marks and coordinates */

  bottom = r.top + 8 + 3 * dlg->lineheight - Ascent () - 2;
  top = bottom - 5;
  i = 0;
  j = offset;
  pos = AlnMgr2MapSeqAlignToBioseq (dlg->salp, j, 2);
  while (pos < 1 && i < dlg->maxchars && j < dlg->uald.aln_length) {
    i++;
    j++;
    pos = AlnMgr2MapSeqAlignToBioseq (dlg->salp, j, 2);
  }
  for (; i < dlg->maxchars + dlg->uald.log10_aln_length && j < dlg->uald.aln_length; i++, j++) {
    ch1 = dlg->aln2 [j];
    if (ch1 != '-') {
      if (dlg->revcomp) {
        realpos = (dlg->uald.new5 + dlg->uald.newa + dlg->uald.new3 - pos - 1);
      } else {
        realpos = pos;
      }
      if (((realpos + 1) % 10) == 0) {
        left = r.left + i * dlg->charwidth + dlg->charwidth / 2 - 1;
        if (i < dlg->maxchars && left <= r.right - dlg->charwidth) {
          MoveTo (left, top);
          LineTo (left, bottom);
        }
        sprintf (str, "%ld", (long) (realpos + 1));
        len = StringLen (str);
        if (len <= j + 1) {
          k = i - len + 1;
          q = 0;
          if (k < 0) {
            q -= k;
            k = 0;
          }
          if (q < len) {
            left = r.left + k * dlg->charwidth;
            MoveTo (left, r.top + 8 + dlg->lineheight);
            while (k < dlg->maxchars && q < len && left <= r.right - dlg->charwidth) {
              PaintChar (str [q]);
              k++;
              q++;
              left += dlg->charwidth;
            }
          }
        }
      } else if (((realpos + 1) % 5) == 0) {
        left = r.left + i * dlg->charwidth + dlg->charwidth / 2 - 1;
        if (i < dlg->maxchars && left <= r.right - dlg->charwidth) {
          MoveTo (left, top + 3);
          LineTo (left, bottom);
        }
      }
      pos++;
    }
  }

  /* draw bottom (old) tick marks and coordinates */

  top = r.top + 8 + 6 * dlg->lineheight - Ascent () + 2;
  bottom = top + 5;
  i = 0;
  j = offset;
  pos = AlnMgr2MapSeqAlignToBioseq (dlg->salp, j, 1);
  while (pos < 1 && i < dlg->maxchars && j < dlg->uald.aln_length) {
    i++;
    j++;
    pos = AlnMgr2MapSeqAlignToBioseq (dlg->salp, j, 1);
  }
  for (; i < dlg->maxchars + dlg->uald.log10_aln_length && j < dlg->uald.aln_length; i++, j++) {
    ch1 = dlg->aln1 [j];
    if (ch1 != '-') {
      if (((pos + 1) % 10) == 0) {
        left = r.left + i * dlg->charwidth + dlg->charwidth / 2 - 1;
        if (i < dlg->maxchars && left <= r.right - dlg->charwidth) {
          MoveTo (left, top);
          LineTo (left, bottom);
        }
        sprintf (str, "%ld", (long) (pos + 1));
        len = StringLen (str);
        if (len <= j + 1) {
          k = i - len + 1;
          q = 0;
          if (k < 0) {
            q -= k;
            k = 0;
          }
          if (q < len) {
            left = r.left + k * dlg->charwidth;
            MoveTo (left, r.top + 8 + 7 * dlg->lineheight);
            while (k < dlg->maxchars && q < len && left <= r.right - dlg->charwidth) {
              PaintChar (str [q]);
              k++;
              q++;
              left += dlg->charwidth;
            }
          }
        }
      } else if (((pos + 1) % 5) == 0) {
        left = r.left + i * dlg->charwidth + dlg->charwidth / 2 - 1;
        if (i < dlg->maxchars && left <= r.right - dlg->charwidth) {
          MoveTo (left, top);
          LineTo (left, bottom - 3);
        }
      }
      pos++;
    }
  }
  SelectFont (systemFont);
}

static void UpdateAlignmentLettersPanelOnScroll (
  BaR sb,
  SlatE slt,
  Int4 newval,
  Int4 oldval
)

{
  RecT                    r;
  AlignmentLettersPanelPtr  dlg;
  Int4                      scroll_distance;

  dlg = (AlignmentLettersPanelPtr) GetObjectExtra (slt);
  if (dlg == NULL) return;

  ObjectRect (dlg->letters, &r);
  InsetRect (&r, 4, 4);
  Select (dlg->letters);
  if (ABS (oldval - newval) < dlg->maxchars) {
    scroll_distance = (oldval - newval) * dlg->charwidth;
    ScrollRect (&r, scroll_distance, 0);
    if (scroll_distance > 0)
    {
      r.left += scroll_distance;
      r.right = r.left + dlg->charwidth;
    }
    else
    {
      r.right += scroll_distance;
      r.left = r.right - dlg->charwidth;
    }
    InvalRect (&r);
  } else {
    InsetRect (&r, -2, -2);
    InvalRect (&r);
  }
  Update ();
}

static void UpdatePairToAlignmentLettersPanel (DialoG d, Pointer data)
{
  AlignmentLettersPanelPtr  dlg;
  UpdatePairPtr             upp;
  CharPtr                   seq_orig, seq_update;
  BaR                       sb;
  RecT                      r;
  
  dlg = (AlignmentLettersPanelPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  dlg->aln1 = MemFree (dlg->aln1);
  dlg->aln2 = MemFree (dlg->aln2);
  dlg->salp = NULL;
  dlg->revcomp = FALSE;
  dlg->salp = NULL;
  MemSet (&(dlg->uald), 0, sizeof (UpdateAlignmentLengthsData));
  
  upp = (UpdatePairPtr) data;
  if (upp != NULL && upp->salp != NULL)
  {  
    dlg->salp = upp->salp;
    dlg->revcomp = upp->revcomp;
  
    CalculateUpdateAlignmentLengths (upp->salp, upp->orig_bsp, upp->update_bsp, &(dlg->uald));
  
    /* create alignment string for original */
    seq_orig = GetSequenceByBsp (upp->orig_bsp);
    dlg->aln1 = MakeUpdateAlignmentSequenceString (upp->salp, dlg->uald.aln_length,
                                                 1, seq_orig);
    seq_orig = MemFree (seq_orig);  

    /* create alignment string for update */
    seq_update = GetSequenceByBsp (upp->update_bsp);
    dlg->aln2 = MakeUpdateAlignmentSequenceString (upp->salp, dlg->uald.aln_length,
                                                 2, seq_update);
    seq_update = MemFree (seq_update);  
  
    sb = GetSlateHScrollBar ((SlatE) dlg->letters);
    SetBarMax (sb, dlg->uald.aln_length - (Int4) dlg->maxchars);
    CorrectBarPage (sb, (Int4) dlg->maxchars - 1, (Int4) dlg->maxchars - 1);
  }

  ObjectRect (dlg->letters, &r);
  InvalRect (&r);
  
}

static DialoG AlignmentLettersPanel (GrouP parent, Int4 prompt_width, Int4 hgt)
{
  AlignmentLettersPanelPtr  dlg;
  GrouP                     p;
  RecT                      r;

  dlg = (AlignmentLettersPanelPtr) MemNew (sizeof (AlignmentLettersPanelData));
  if (dlg == NULL)
  {
    return NULL;
  }
  
  p = HiddenGroup (parent, 1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  dlg->dialog = (DialoG) p;
  dlg->todialog = UpdatePairToAlignmentLettersPanel;
  
  dlg->aln1 = NULL;
  dlg->aln2 = NULL;
  dlg->salp = NULL;
  dlg->revcomp = FALSE;
  MemSet (&(dlg->uald), 0, sizeof (UpdateAlignmentLengthsData));
  
  dlg->letters = AutonomousPanel4 (p, prompt_width + Nlm_vScrollBarWidth, hgt,
                                   UpdateAlignmentLettersPanelOnDraw,
				                           NULL, UpdateAlignmentLettersPanelOnScroll, 0, NULL, NULL);
  SetObjectExtra (dlg->letters, (Pointer) dlg, NULL);  

  SelectFont (SetSmallFont ());
  ObjectRect (dlg->letters, &r);
  InsetRect (&r, 4, 4);
  dlg->lineheight = LineHeight ();
  dlg->charwidth = MaxCharWidth ();
  dlg->maxchars = (r.right-r.left-2+dlg->charwidth - 1) / dlg->charwidth;
  SelectFont (systemFont);
  return (DialoG) p;
}

static Int4 GetAlignmentLettersPanelMaxchars (DialoG d)
{
  AlignmentLettersPanelPtr  dlg;
  
  dlg = (AlignmentLettersPanelPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return 0;
  }
  else
  {
    return dlg->maxchars;
  }
}

static Int4 GetAlignmentLettersScrollPosition (DialoG d)
{
  AlignmentLettersPanelPtr  dlg;
  BaR                       sb;
  
  dlg = (AlignmentLettersPanelPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return 0;
  }

  sb = GetSlateHScrollBar ((SlatE) dlg->letters);

  return GetBarValue (sb);    
}

static void ScrollAlignmentLettersPanel (DialoG d, Int4 pos)
{
  AlignmentLettersPanelPtr  dlg;
  BaR                       sb;
  
  dlg = (AlignmentLettersPanelPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }

  sb = GetSlateHScrollBar ((SlatE) dlg->letters);

  SetBarValue (sb, pos);  
}

typedef struct updatetitlesdialog
{
  DIALOG_MESSAGE_BLOCK
  
  PrompT  new_sequence_id_ppt; /* indicate which sequence is being used for the update */
  PrompT  new_sequence_length_ppt;
  PrompT  old_sequence_id_ppt; /* indicates which sequence is being updated */
  PrompT  old_sequence_length_ppt;
                                    
  BioseqPtr          orig_bsp;   /* original bioseq */
  BioseqPtr          update_bsp; /* update bioseq */
  Uint2              entityID_orig;
  Uint2              entityID_update;
  Boolean            is_indexer;
  Boolean            is_update;
} UpdateTitlesDialogData, PNTR UpdateTitlesDialogPtr;

static Pointer UpdatePairFromUpdateTitlesDialog (DialoG d)
{
  UpdateTitlesDialogPtr dlg;
  UpdatePairPtr          upp;

  dlg = (UpdateTitlesDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  
  upp = (UpdatePairPtr) MemNew (sizeof (UpdatePairData));
  if (upp != NULL)
  {
    upp->orig_bsp = dlg->orig_bsp;
    upp->update_bsp = dlg->update_bsp;
    upp->revcomp = FALSE;
    upp->salp = NULL;
  }
  
  return upp;
}

static void UpdatePairToUpdateTitlesDialog (DialoG d, Pointer data)
{
  UpdateTitlesDialogPtr dlg;
  UpdatePairPtr         upp;
  Char                  ppt_txt [500];
  Char                  id_txt [MAX_ID_LEN];
  
  dlg = (UpdateTitlesDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  dlg->orig_bsp = NULL;
  dlg->update_bsp = NULL;
  
  upp = (UpdatePairPtr) data;
  
  if (upp == NULL)
  {
    SafeHide (dlg->new_sequence_id_ppt);
    SafeHide (dlg->new_sequence_length_ppt);
    SafeHide (dlg->old_sequence_id_ppt);    
    SafeHide (dlg->old_sequence_length_ppt);    
    return;
  }
  
  dlg->orig_bsp = upp->orig_bsp;
  dlg->update_bsp = upp->update_bsp;
  
  if (upp->orig_bsp == NULL)
  {
    SafeHide (dlg->new_sequence_id_ppt);
    SafeHide (dlg->new_sequence_length_ppt);
    SafeHide (dlg->old_sequence_id_ppt);    
    SafeHide (dlg->old_sequence_length_ppt);    
  }  
  /* show no update message */
  else if (upp->update_bsp == NULL)
  {
    SeqIdWrite (SeqIdFindBest (upp->orig_bsp->id, SEQID_GENBANK), id_txt, 
                PRINTID_REPORT, sizeof (id_txt) - 50);
    if (dlg->is_indexer)
    {
      sprintf (ppt_txt, "No update sequence for %s", id_txt);
      SetTitle (dlg->new_sequence_id_ppt, ppt_txt);
      sprintf (ppt_txt, "Length %ld", upp->orig_bsp->length);
      SetTitle (dlg->new_sequence_length_ppt, ppt_txt);
    }
    else
    {
      sprintf (ppt_txt, "No update sequence for %s (length %ld)", 
               id_txt, upp->orig_bsp->length);
      SetTitle (dlg->new_sequence_id_ppt, ppt_txt);
    }
    SafeHide (dlg->old_sequence_id_ppt);    
    SafeHide (dlg->old_sequence_length_ppt);    
  }  
  /* show no alignment message */
  else if (upp->salp == NULL && dlg->is_update)
  {
    SeqIdWrite (SeqIdFindBest (upp->orig_bsp->id, SEQID_GENBANK), id_txt, 
                PRINTID_REPORT, sizeof (id_txt) - 50);
    if (dlg->is_indexer)
    {
      sprintf (ppt_txt, "No alignment for %s", id_txt);
      SetTitle (dlg->new_sequence_id_ppt, ppt_txt);
      sprintf (ppt_txt, "Length: %ld", dlg->orig_bsp->length);
      SetTitle (dlg->new_sequence_length_ppt, ppt_txt);
    }
    else
    {
      sprintf (ppt_txt, "No alignment for %s Length: %ld", id_txt, dlg->orig_bsp->length);
      SetTitle (dlg->new_sequence_id_ppt, ppt_txt);
    }
    SafeHide (dlg->old_sequence_id_ppt);    
    SafeHide (dlg->old_sequence_length_ppt);    
  }
  else
  {
    /* update sequence title prompts */ 
    SeqIdWrite (SeqIdFindBest (dlg->update_bsp->id, SEQID_GENBANK), id_txt, 
                               PRINTID_REPORT, sizeof (id_txt) - 50);
    if (dlg->is_indexer)
    {
      sprintf (ppt_txt, "New sequence: %s", id_txt);
      SetTitle (dlg->new_sequence_id_ppt, ppt_txt);
      sprintf (ppt_txt, "Length: %ld", dlg->update_bsp->length);
      SetTitle (dlg->new_sequence_length_ppt, ppt_txt);
    }
    else
    {
      sprintf (ppt_txt, "New sequence: %s Length: %ld", id_txt, dlg->update_bsp->length);
      SetTitle (dlg->new_sequence_id_ppt, ppt_txt);
    }
  
    SeqIdWrite (SeqIdFindBest (dlg->orig_bsp->id, SEQID_GENBANK), id_txt, 
                               PRINTID_REPORT, sizeof (id_txt) - 50);
    if (dlg->is_indexer)
    {
      sprintf (ppt_txt, "Old sequence: %s", id_txt);
      SetTitle (dlg->old_sequence_id_ppt, ppt_txt);
      sprintf (ppt_txt, "Length: %ld", dlg->orig_bsp->length);
      SetTitle (dlg->old_sequence_length_ppt, ppt_txt);
    }
    else
    {
      sprintf (ppt_txt, "Old sequence: %s Length: %ld", id_txt, dlg->orig_bsp->length);
      SetTitle (dlg->old_sequence_id_ppt, ppt_txt);
    }
  
    SafeShow (dlg->new_sequence_id_ppt);
    SafeShow (dlg->new_sequence_length_ppt);
    SafeShow (dlg->old_sequence_id_ppt);  
    SafeShow (dlg->old_sequence_length_ppt);  
  }
}

static DialoG UpdateTitlesDialog (GrouP parent, Boolean is_indexer, Boolean is_update)
{
  UpdateTitlesDialogPtr dlg;
  GrouP                 p;
  Int4                  prompt_width = 450;
  
  dlg = (UpdateTitlesDialogPtr) MemNew (sizeof (UpdateTitlesDialogData));
  if (dlg == NULL)
  {
    return NULL;
  }
  
  p = HiddenGroup (parent, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 2, 2);
  
  dlg->dialog = (DialoG) p;
  
  dlg->todialog = UpdatePairToUpdateTitlesDialog;
  dlg->fromdialog = UpdatePairFromUpdateTitlesDialog;
  
  dlg->is_indexer = is_indexer;
  dlg->is_update = is_update;
  
  if (is_indexer)
  {
    prompt_width = 400;
  }
  else
  {
    prompt_width = 550;
  }
  
  dlg->new_sequence_id_ppt = StaticPrompt (p, "", prompt_width, 
                                           dialogTextHeight, programFont, 'l');
  if (is_indexer)
  {
    dlg->new_sequence_length_ppt = StaticPrompt (p, "", prompt_width, 
                                                 dialogTextHeight, programFont, 'l');
  }
  else
  {
    dlg->new_sequence_length_ppt = NULL;
  }
  
  dlg->old_sequence_id_ppt = StaticPrompt (p, "", prompt_width, dialogTextHeight, programFont, 'l');
  if (is_indexer)
  {
    dlg->old_sequence_length_ppt = StaticPrompt (p, "", prompt_width, 
                                                 dialogTextHeight, programFont, 'l');
  }
  
  return (DialoG) p;
}

typedef struct updatepreviewdialog
{
  DIALOG_MESSAGE_BLOCK

  GrouP   viewer_grp;       /* group that contains viewers and prompts
                             * this group is hidden when there is no alignment
                             */
  VieweR  overview;       /* shows sequence lengths/overlap */
  GrouP   details_grp;
  VieweR  details;        /* shows individual mismatches/deletions/insertions */
  GrouP   letters_grp;
  DialoG  letters_dlg;      /* shows sequence and gap characters and coordinates */

  SegmenT overview_picture; /* picture of sequence lengths/overlap */
  SegmenT details_picture;  /* picture of mismatches/deletions/insertions */
  Int4    details_scaleX;   /* scale at which the details picture is displayed */                                                                
  
  UpdateAlignmentLengthsData uald; /* structure that holds the length of the alignment,
                                    * the lengths of the 3' and 5' overlaps, etc.
                                    */
  ValNodePtr         indels;     /* ValNode structure that holds the positions
                                  * of insertions and deletions (in alignment
                                  * coordinates)
                                  */    
  ValNodePtr         mismatches; /* ValNode structure that holds the positions
                                  * of mismatches (in alignment coordinates)
                                  */    
                                  
  SeqAlignPtr        salp;       /* alignment used for update */
                                 /* The UpdatePreviewDialog is responsible for freeing
                                  * the alignment passed to it.
                                  */
  BioseqPtr          orig_bsp;   /* original bioseq */
  BioseqPtr          update_bsp; /* update bioseq */
  Uint2              entityID_orig;
  Uint2              entityID_update;
  Boolean            revcomp;              
  
  Boolean            ignore_alignment;
  Boolean            extend5;                       
} UpdatePreviewDialogData, PNTR UpdatePreviewDialogPtr;

/* if the user clicks on the details picture, scroll to the appropriate position
 * in the letters panel.
 */
static void UpdateAlignmentDetailsOnClick (VieweR vwr, SegmenT pict, PoinT pt)

{
  Int4                   goHere;
  Int4                   offset;
  Int4                   maxchars;
  Int4                   maxover2;
  PntInfo                pnt;
  UpdatePreviewDialogPtr dlg;
  ValNodePtr             vnp;

  dlg = (UpdatePreviewDialogPtr) GetViewerData (vwr);
  if (dlg == NULL) return;

  MapViewerToWorld (vwr, pt, &pnt);
  maxchars = GetAlignmentLettersPanelMaxchars (dlg->letters_dlg);
  maxover2 = maxchars / 2;
  if (pnt.x <= 0) {
    pnt.x = 0;
  } else if (pnt.x >= dlg->uald.aln_length) {
    pnt.x = dlg->uald.aln_length  - maxchars;
  } else if (pnt.x >= maxover2) {

    offset = GetAlignmentLettersScrollPosition (dlg->letters_dlg);

    /* look for clicks within 5 pixels of an indel start or a mismatch */

    goHere = -1;
    for (vnp = dlg->indels; vnp != NULL && goHere < 0; vnp = vnp->next) {
      if (ABS (pnt.x - vnp->data.intvalue) < dlg->details_scaleX * 5) {
        goHere = vnp->data.intvalue;
      }
    }
    for (vnp = dlg->mismatches; vnp != NULL && goHere < 0; vnp = vnp->next) {
      if (ABS (pnt.x - vnp->data.intvalue) < dlg->details_scaleX * 5) {
        goHere = vnp->data.intvalue;
      }
    }

    if (goHere >= 0) {
      pnt.x = goHere;
    } else {
      /* if already visible, no need to scroll */
      if (pnt.x - maxover2 > offset && pnt.x - maxover2 < offset + maxover2 - 5) return;
      if (pnt.x - maxover2 < offset && pnt.x - maxover2 > offset - maxover2 + 5) return;
    }

    /* go left 1/2 screen so desired point is in the middle */

    pnt.x -= maxover2;
  }

  ResetClip ();
  ScrollAlignmentLettersPanel (dlg->letters_dlg, pnt.x);
  Update ();
}

static Pointer UpdatePairFromUpdatePreviewDialog (DialoG d)
{
  UpdatePreviewDialogPtr dlg;
  UpdatePairPtr          upp;

  dlg = (UpdatePreviewDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  
  upp = (UpdatePairPtr) MemNew (sizeof (UpdatePairData));
  if (upp != NULL)
  {
    upp->orig_bsp = dlg->orig_bsp;
    upp->update_bsp = dlg->update_bsp;
    upp->revcomp = dlg->revcomp;
    upp->salp = dlg->salp;
  }
  
  return upp;
}

static void RedrawAlignmentPreview (DialoG d)
{
  UpdatePreviewDialogPtr dlg;
  UpdatePairData         upd;
  Int4                   scaleX;

  dlg = (UpdatePreviewDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  upd.orig_bsp = dlg->orig_bsp;
  upd.update_bsp = dlg->update_bsp;
  upd.salp = dlg->salp;
  upd.revcomp = dlg->revcomp;

  /* reverse the sequence for drawing purposes */
  if (upd.revcomp)
  {
    BioseqRevComp (upd.update_bsp);
    ReverseBioseqFeatureStrands (upd.update_bsp);
    SeqMgrReplaceInBioseqIndex (upd.update_bsp);
  }

  if (dlg->salp == NULL 
      || dlg->orig_bsp == NULL
      || dlg->update_bsp == NULL)
  {
    dlg->overview_picture = CreatePicture ();
    AttachPicture (dlg->overview, dlg->overview_picture, 0, 0, UPPER_LEFT,
		               1, 1, FrameVwr);
    dlg->details_picture = CreatePicture ();
    AttachPicture (dlg->details, dlg->details_picture, 0, 0, UPPER_LEFT,
		               1, 1, FrameVwr);
    MemSet (&(dlg->uald), 0, sizeof (UpdateAlignmentLengthsData));

    PointerToDialog (dlg->letters_dlg, NULL);
    
    Hide (dlg->viewer_grp);
    /* reverse the sequence for drawing purposes */
    if (upd.revcomp)
    {
      BioseqRevComp (upd.update_bsp);
      ReverseBioseqFeatureStrands (upd.update_bsp);
      SeqMgrReplaceInBioseqIndex (upd.update_bsp);
    }

    return;
  }
    
  CalculateUpdateAlignmentLengths (dlg->salp, dlg->orig_bsp, dlg->update_bsp, &(dlg->uald));

  if (dlg->ignore_alignment)
  {
    dlg->overview_picture = MakeExtensionPicture (dlg->orig_bsp->length, dlg->update_bsp->length, dlg->extend5);
  }
  else
  {
    dlg->overview_picture = MakeAlignmentOverviewPicture (&(dlg->uald));
  }

  scaleX = CalculateBestUpdateAlignmentViewerScale (dlg->overview, dlg->overview_picture);
  AttachPicture (dlg->overview, dlg->overview_picture, 0, 0, UPPER_LEFT,
		             scaleX, 1, FrameVwr);

  dlg->indels = ValNodeFree (dlg->indels);
  dlg->details_picture = MakeAlignmentDetailsPicture (dlg->salp, 
                                                      dlg->orig_bsp,
                                                      dlg->update_bsp,
                                                      &(dlg->indels),
                                                      &(dlg->mismatches),
                                                      &(dlg->uald), dlg->revcomp);
  dlg->details_scaleX = CalculateBestUpdateAlignmentViewerScale (dlg->details, dlg->details_picture);
  AttachPicture (dlg->details, dlg->details_picture, 0, 0, UPPER_LEFT,
		             dlg->details_scaleX, 1, FrameVwr);
  SetViewerData (dlg->details, (Pointer) dlg, NULL);
  SetViewerProcs (dlg->details, UpdateAlignmentDetailsOnClick, NULL, NULL, NULL);

  /* sort the insertions and deletions so they will be in order */
  dlg->indels = ValNodeSort (dlg->indels, SortVnpByInt);
  /* sort the matches so that they will be in order */
  dlg->mismatches = ValNodeSort (dlg->mismatches, SortVnpByInt);

  PointerToDialog (dlg->letters_dlg, &upd);
  Show (dlg->viewer_grp);
  if (dlg->ignore_alignment)
  {
    Hide (dlg->details_grp);
    Hide (dlg->letters_grp);
  }
  else
  {
    Show (dlg->details_grp);
    Show (dlg->letters_grp);
  }

  /* unreverse the sequence if it was reversed */
  if (upd.revcomp)
  {
    BioseqRevComp (upd.update_bsp);
    ReverseBioseqFeatureStrands (upd.update_bsp);
    SeqMgrReplaceInBioseqIndex (upd.update_bsp);
  }  
  
}

static void AlignmentToUpdatePreviewDialog (DialoG d, Pointer data)
{
  UpdatePreviewDialogPtr dlg;
  UpdatePairPtr          upp;
  
  dlg = (UpdatePreviewDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  dlg->overview_picture     = DeletePicture (dlg->overview_picture);
  dlg->details_picture     = DeletePicture (dlg->details_picture);
  dlg->salp = NULL;
  dlg->orig_bsp = NULL;
  dlg->update_bsp = NULL;
  dlg->revcomp = FALSE;
  dlg->salp = SeqAlignFree (dlg->salp);
  
  upp = (UpdatePairPtr) data;
  
  if (upp != NULL)
  {
    dlg->salp = upp->salp;
    dlg->revcomp = upp->revcomp;
    dlg->orig_bsp = upp->orig_bsp;
    dlg->update_bsp = upp->update_bsp;
    if (BioseqFind (dlg->update_bsp->id) != dlg->update_bsp)
    {
      SeqMgrReplaceInBioseqIndex (dlg->update_bsp);
    }
  }
  
  RedrawAlignmentPreview (d);
}

static void SetIgnoreAndExtendForPreviewPictures (DialoG d, Boolean ignore_alignment, Boolean extend5)
{
  UpdatePreviewDialogPtr dlg;
  
  dlg = (UpdatePreviewDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  dlg->ignore_alignment = ignore_alignment;
  dlg->extend5 = extend5;
  RedrawAlignmentPreview (d);
}

static void CleanupPreviewDialog (GraphiC g, Pointer data)
{
  UpdatePreviewDialogPtr dlg;

  dlg = (UpdatePreviewDialogPtr) data;
  if (dlg != NULL)
  {
    dlg->indels = ValNodeFree (dlg->indels);
    dlg->mismatches = ValNodeFree (dlg->mismatches);
    dlg->salp = SeqAlignFree (dlg->salp);
  }
  StdCleanupExtraProc (g, data);
}

static DialoG UpdatePreviewDialog (GrouP parent, Boolean is_indexer)
{
  UpdatePreviewDialogPtr dlg;
  GrouP                  p;
  Int4                   prompt_width = 400, hgt;
  GrouP                  ppt1, ppt2, ppt3;
  
  dlg = (UpdatePreviewDialogPtr) MemNew (sizeof (UpdatePreviewDialogData));
  if (dlg == NULL)
  {
    return NULL;
  }
  
  p = HiddenGroup (parent, -1, 0, NULL);
  SetObjectExtra (p, dlg, CleanupPreviewDialog);
  SetGroupSpacing (p, 2, 2);
  
  dlg->dialog = (DialoG) p;
  
  dlg->todialog = AlignmentToUpdatePreviewDialog;
  dlg->fromdialog = UpdatePairFromUpdatePreviewDialog;

  dlg->viewer_grp = HiddenGroup (p, -1, 0, NULL);
  SetGroupSpacing (dlg->viewer_grp, 2, 2);
  ppt1 = MultiLinePrompt (dlg->viewer_grp, txt1, prompt_width, programFont);
  dlg->overview = CreateViewer (dlg->viewer_grp, prompt_width + Nlm_vScrollBarWidth, 100,
				FALSE, FALSE);
  
  dlg->details_grp = HiddenGroup (dlg->viewer_grp, -1, 0, NULL);
  SetGroupSpacing (dlg->details_grp, 2, 2);
  ppt2 = MultiLinePrompt (dlg->details_grp, txt2, prompt_width, programFont);
  dlg->details = CreateViewer (dlg->details_grp, prompt_width + Nlm_vScrollBarWidth, 80,
			       FALSE, FALSE);
  AlignObjects (ALIGN_CENTER, (HANDLE) ppt2, (HANDLE) dlg->details, NULL);
  
  dlg->letters_grp = HiddenGroup (dlg->viewer_grp, -1, 0, NULL);
  SetGroupSpacing (dlg->letters_grp, 2, 2);
  ppt3 = MultiLinePrompt (dlg->letters_grp, txt3, prompt_width, programFont);
    
#ifdef WIN_MAC
  hgt = 90;
#else
  hgt = 110;
#endif

  dlg->letters_dlg = AlignmentLettersPanel (dlg->letters_grp, prompt_width, hgt);
  AlignObjects (ALIGN_CENTER, (HANDLE) ppt3, (HANDLE) dlg->letters_dlg, NULL);
   
  dlg->overview_picture = NULL;
  dlg->details_picture = NULL;
  dlg->indels = NULL;
  dlg->mismatches = NULL;
  MemSet (&(dlg->uald), 0, sizeof (UpdateAlignmentLengthsData));
  
  

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->overview,
                              (HANDLE) dlg->details_grp, 
                              (HANDLE) dlg->letters_grp,
                              (HANDLE) ppt1, 
                              NULL);
  
  return (DialoG) p;
}

typedef struct submitterupdateoptionsdialog
{
  DIALOG_MESSAGE_BLOCK
  GrouP                sequence_update_type;
  GrouP                feature_update_type;
  GrouP                feature_remove_type;
  
  ButtoN               sequence_update_btns [eSequenceUpdateExtend3];
  
  ButtoN               ignore_alignment;
  
  Boolean              do_update;
  
  Nlm_ChangeNotifyProc change_notify;
  Pointer              change_userdata;

} SubmitterUpdateOptionsDialogData, PNTR SubmitterUpdateOptionsDialogPtr;

static void SubmitterUpdateOptionsToDialog (DialoG d, Pointer data)
{
  SubmitterUpdateOptionsDialogPtr dlg;
  SubmitterUpdateOptionsPtr       suop;
  
  dlg = (SubmitterUpdateOptionsDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  suop = (SubmitterUpdateOptionsPtr) data;
  if (dlg->do_update)
  {
    if (suop == NULL)
    {
      SetValue (dlg->sequence_update_type, eSequenceUpdateNoChange);
      SetValue (dlg->feature_update_type, eFeatureUpdateNoChange);
      SetValue (dlg->feature_remove_type, eFeatureRemoveNone);
      SetStatus (dlg->ignore_alignment, FALSE);
    }
    else
    {
      SetValue (dlg->sequence_update_type, suop->sequence_update_type);
      SetValue (dlg->feature_update_type, suop->feature_update_type);
      SetValue (dlg->feature_remove_type, suop->feature_remove_type);
      SetStatus (dlg->ignore_alignment, suop->ignore_alignment);
    }
  }
  else
  {
  	if (suop == NULL)
  	{
  	  SetValue (dlg->sequence_update_type, 1);
      SetValue (dlg->feature_update_type, eFeatureUpdateNoChange);
      SetValue (dlg->feature_remove_type, eFeatureRemoveNone);  	  
  	}
  	else
  	{
  	  if (suop->sequence_update_type == eSequenceUpdateExtend3)
  	  {
        SetValue (dlg->sequence_update_type, 2);
  	  }
  	  else
  	  {
        SetValue (dlg->sequence_update_type, 1);
  	  }
  	  
      SetValue (dlg->feature_update_type, suop->feature_update_type);
      SetValue (dlg->feature_remove_type, suop->feature_remove_type);
  	}
  }
  
  if (dlg->change_notify != NULL)
  {
    (dlg->change_notify) (dlg->change_userdata);
  }
  
}

static Pointer SubmitterUpdateOptionsFromDialog (DialoG d)
{
  SubmitterUpdateOptionsDialogPtr dlg;
  SubmitterUpdateOptionsPtr       suop;
  
  dlg = (SubmitterUpdateOptionsDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  
  suop = (SubmitterUpdateOptionsPtr) MemNew (sizeof (SubmitterUpdateOptionsData));
  if (suop != NULL)
  {
    if (dlg->do_update)
    {
      suop->sequence_update_type = (ESequenceUpdateType) GetValue (dlg->sequence_update_type);
      if (suop->sequence_update_type > eSequenceUpdateNoChange
          && ! Enabled (dlg->sequence_update_btns [suop->sequence_update_type - 1]))
      {
        suop->sequence_update_type = eSequenceUpdateNoChange;
      }
    }
    else
    {
      if (GetValue (dlg->sequence_update_type) == 2)
      {
      	suop->sequence_update_type = eSequenceUpdateExtend3;
      }
      else
      {
      	suop->sequence_update_type = eSequenceUpdateExtend5;
      }
    }
    
    suop->feature_update_type = (EFeatureUpdateType) GetValue (dlg->feature_update_type);
    suop->feature_remove_type = (EFeatureRemoveType) GetValue (dlg->feature_remove_type);
    
    if (dlg->do_update)
    {
      if (Enabled (dlg->ignore_alignment)
          && (suop->sequence_update_type == eSequenceUpdateExtend5
              || suop->sequence_update_type == eSequenceUpdateExtend3))
      {
        suop->ignore_alignment = GetStatus (dlg->ignore_alignment);
      }
      else
      {
        suop->ignore_alignment = FALSE;
      }
    }
    else
    {
      suop->ignore_alignment = TRUE;
    }
  }
  return suop;
}

static void DisableSubmitterImportFeatureOptions (DialoG d)
{
  SubmitterUpdateOptionsDialogPtr dlg;
  
  dlg = (SubmitterUpdateOptionsDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  Disable (dlg->feature_update_type);
}

static void EnableSubmitterImportFeatureOptions (DialoG d)
{
  SubmitterUpdateOptionsDialogPtr dlg;
  
  dlg = (SubmitterUpdateOptionsDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  Enable (dlg->feature_update_type);
}

static void DisableSubmitterRemoveFeatureOptions (DialoG d)
{
  SubmitterUpdateOptionsDialogPtr dlg;
  
  dlg = (SubmitterUpdateOptionsDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  Disable (dlg->feature_remove_type);
}

static void EnableSubmitterRemoveFeatureOptions (DialoG d)
{
  SubmitterUpdateOptionsDialogPtr dlg;
  
  dlg = (SubmitterUpdateOptionsDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  Enable (dlg->feature_remove_type);
}

static void AdjustSequenceUpdateOptionsEnabled
(DialoG                    d, 
 UpdatePairPtr             upp,
 UpdateAlignmentLengthsPtr ualp,
 Boolean                   is_indexer)
{
  SubmitterUpdateOptionsDialogPtr dlg;
  ESequenceUpdateType             action;
  Boolean                         ignore_alignment;
  SubmitterUpdateOptionsPtr       suop;
  
  dlg = (SubmitterUpdateOptionsDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  if (!dlg->do_update)
  {
  	Enable (dlg->sequence_update_btns [0]);
  	Enable (dlg->sequence_update_btns [1]);
  	return;
  }
  
  suop = DialogToPointer (d);
  
  if (upp == NULL || upp->salp == NULL) 
  {
    Disable (dlg->ignore_alignment);
    ignore_alignment = TRUE;
  }
  else if (suop != NULL 
          && (suop->sequence_update_type == eSequenceUpdateNoChange
              || suop->sequence_update_type == eSequenceUpdatePatch))
  {
    Disable (dlg->ignore_alignment);
    ignore_alignment = FALSE;
  }
  else
  {
  
    Enable (dlg->ignore_alignment);
    ignore_alignment = GetStatus (dlg->ignore_alignment);
  }
  
  for (action = eSequenceUpdateNoChange;
       action <= eSequenceUpdateExtend3;
       action++)
  {
    if (IsSequenceUpdateChoiceAllowed (action, upp, ualp, is_indexer,
                                       ignore_alignment))
    {
      Enable (dlg->sequence_update_btns [action - 1]);
    }
    else
    {
      Disable (dlg->sequence_update_btns [action - 1]);
    }
  }
      
}

static void IgnoreAlignmentBtn (ButtoN b)
{
  SubmitterUpdateOptionsDialogPtr dlg;

  dlg = (SubmitterUpdateOptionsDialogPtr) GetObjectExtra (b);
  if (dlg != NULL && dlg->change_notify != NULL)
  {
    (dlg->change_notify) (dlg->change_userdata);
  }
}

static void UpdateTypeGroup (GrouP g)
{
  SubmitterUpdateOptionsDialogPtr dlg;

  dlg = (SubmitterUpdateOptionsDialogPtr) GetObjectExtra (g);
  if (dlg != NULL && dlg->change_notify != NULL)
  {
    (dlg->change_notify) (dlg->change_userdata);
  }
}

static DialoG 
SubmitterUpdateOptionsDialog 
(GrouP                parent,
 Boolean              is_indexer,
 Boolean              do_update,
 Nlm_ChangeNotifyProc change_notify,
 Pointer              change_userdata)
{
  SubmitterUpdateOptionsDialogPtr dlg;
  GrouP                  p, j, k, g;
  
  dlg = (SubmitterUpdateOptionsDialogPtr) MemNew (sizeof (SubmitterUpdateOptionsDialogData));
  if (dlg == NULL)
  {
    return NULL;
  }
  
  if (is_indexer)
  {
    p = HiddenGroup (parent, -1, 0, NULL);
  }
  else
  {
    p = HiddenGroup (parent, 1, 0, NULL);
  }
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = SubmitterUpdateOptionsToDialog;
  dlg->fromdialog = SubmitterUpdateOptionsFromDialog;
  
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;
  
  dlg->do_update = do_update;
  
  if (dlg->do_update)
  {
    j = NormalGroup (p, -1, 0, "Sequence Update", programFont, NULL);
    if (is_indexer)
    {
      dlg->sequence_update_type = HiddenGroup (j, 0, 1, UpdateTypeGroup);
    }
    else
    {
      SetGroupSpacing (j, 3, 10);
      dlg->sequence_update_type = HiddenGroup (j, 1, 0, UpdateTypeGroup);
      SetGroupSpacing (dlg->sequence_update_type, 3, 2);
    }
    SetObjectExtra (dlg->sequence_update_type, dlg, NULL);
    dlg->sequence_update_btns [0] = RadioButton (dlg->sequence_update_type, "No change");
    dlg->sequence_update_btns [1] = RadioButton (dlg->sequence_update_type, "Replace");
    dlg->sequence_update_btns [2] = RadioButton (dlg->sequence_update_type, "Patch");
    dlg->sequence_update_btns [3] = RadioButton (dlg->sequence_update_type, "Extend 5'");
    dlg->sequence_update_btns [4] = RadioButton (dlg->sequence_update_type, "Extend 3'");
    SetValue (dlg->sequence_update_type, eSequenceUpdateNoChange);
    dlg->ignore_alignment = CheckBox (j, "Ignore alignment", IgnoreAlignmentBtn);
    SetObjectExtra (dlg->ignore_alignment, dlg, NULL);
    if (is_indexer)
    {
      AlignObjects (ALIGN_CENTER, (HANDLE) dlg->sequence_update_type,
                                  (HANDLE) dlg->ignore_alignment,
                                  NULL);
    }
    else
    {
      AlignObjects (ALIGN_LEFT, (HANDLE) dlg->sequence_update_type,
                                (HANDLE) dlg->ignore_alignment,
                                NULL);
    }
  }
  else
  {
  	j = NormalGroup (p, -1, 0, "Extend Sequence", programFont, NULL);
    dlg->sequence_update_type = HiddenGroup (j, 0, 1, NULL);
    dlg->sequence_update_btns [0] = RadioButton (dlg->sequence_update_type, "Extend 5'");
    dlg->sequence_update_btns [1] = RadioButton (dlg->sequence_update_type, "Extend 3'");
    dlg->sequence_update_btns [2] = NULL;
    dlg->sequence_update_btns [3] = NULL;
    dlg->sequence_update_btns [4] = NULL;
    
    SetValue (dlg->sequence_update_type, 1);
    dlg->ignore_alignment = NULL;  	
  }

  if (is_indexer)
  {
    k = HiddenGroup (p, 2, 0, NULL);
  }
  else
  {
    k = p;
  }
  
  g = NormalGroup (k, 1, 0, "Existing Features", programFont, NULL);
  dlg->feature_remove_type = HiddenGroup (g, 0, 6, NULL);
  RadioButton (dlg->feature_remove_type, "Do not remove");
  RadioButton (dlg->feature_remove_type, "Remove in aligned area");
  RadioButton (dlg->feature_remove_type, "Remove outside aligned area");
  RadioButton (dlg->feature_remove_type, "Remove all");
  SetValue (dlg->feature_remove_type, eFeatureRemoveNone);
    
  g = NormalGroup (k, 1, 0, "Import Features", programFont, NULL);
  dlg->feature_update_type = HiddenGroup (g, 0, 6, NULL);
  RadioButton (dlg->feature_update_type, "Do not import features");
  RadioButton (dlg->feature_update_type, "Import all except duplicates");
  RadioButton (dlg->feature_update_type, "Import all, merge duplicates");
  RadioButton (dlg->feature_update_type, "Import all, replace duplicates");
  RadioButton (dlg->feature_update_type, "Import all, including duplicates");
  SetValue (dlg->feature_update_type, eFeatureUpdateNoChange);
  
  if (is_indexer)
  {
    AlignObjects (ALIGN_CENTER, (HANDLE) j,
                                (HANDLE) k,
                                NULL);
  }
  
  return (DialoG) p;
}

typedef struct indexeroptionsdialog
{
  DIALOG_MESSAGE_BLOCK
  ButtoN keep_protein_ids;
  ButtoN add_cit_subs;
  ButtoN update_quality_scores;
  ButtoN update_proteins;
  GrouP  protein_options;
  ButtoN truncate_proteins;
  ButtoN extend_proteins3;
  ButtoN extend_proteins5;
  ButtoN correct_cds_genes;
} IndexerOptionsDialogData, PNTR IndexerOptionsDialogPtr;

static void EnableIndexerOptions (ButtoN b)
{
  IndexerOptionsDialogPtr dlg;

  dlg = (IndexerOptionsDialogPtr) GetObjectExtra (b);
  if (dlg == NULL)
  {
    return;
  }
  
  if (GetStatus (dlg->update_proteins))
  {
    Enable (dlg->protein_options);
  }
  else
  {
    Disable (dlg->protein_options);
  }
}

static void IndexerOptionsToDialog (DialoG d, Pointer data)
{
  IndexerOptionsDialogPtr dlg;
  IndexerOptionsPtr       iop;
  
  dlg = (IndexerOptionsDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  iop = (IndexerOptionsPtr) data;
  
  if (iop == NULL)
  {
    SetStatus (dlg->add_cit_subs, FALSE);
    SetStatus (dlg->update_quality_scores, TRUE);
    SafeSetStatus (dlg->update_proteins, FALSE);
    SafeSetStatus (dlg->truncate_proteins, FALSE);
    SafeSetStatus (dlg->extend_proteins3, FALSE);
    SafeSetStatus (dlg->extend_proteins5, FALSE);
    SafeSetStatus (dlg->correct_cds_genes, FALSE);
    SafeSetStatus (dlg->keep_protein_ids, FALSE);
  }
  else
  {
    SetStatus (dlg->add_cit_subs, iop->add_cit_subs);
    SetStatus (dlg->update_quality_scores, iop->update_quality_scores);
    SafeSetStatus (dlg->update_proteins, iop->update_proteins);
    SafeSetStatus (dlg->keep_protein_ids, iop->keep_protein_ids);
    if (iop->update_proteins)
    {
      SafeSetStatus (dlg->truncate_proteins, iop->truncate_proteins);
      SafeSetStatus (dlg->extend_proteins3, iop->extend_proteins3);
      SafeSetStatus (dlg->extend_proteins5, iop->extend_proteins5);
      SafeSetStatus (dlg->correct_cds_genes, iop->correct_cds_genes);
    }
    else
    {
      SafeSetStatus (dlg->truncate_proteins, FALSE);
      SafeSetStatus (dlg->extend_proteins3, FALSE);
      SafeSetStatus (dlg->extend_proteins5, FALSE);
      SafeSetStatus (dlg->correct_cds_genes, FALSE);
    }
  }
  
  EnableIndexerOptions (dlg->update_proteins);
}

static Pointer IndexerOptionsFromDialog (DialoG d)
{
  IndexerOptionsDialogPtr dlg;
  IndexerOptionsPtr       iop;
  
  dlg = (IndexerOptionsDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  
  iop = (IndexerOptionsPtr) MemNew (sizeof (IndexerOptionsData));
  
  if (iop == NULL)
  {
    return NULL;
  }
  iop->add_cit_subs = GetStatus (dlg->add_cit_subs);
  iop->update_quality_scores = GetStatus (dlg->add_cit_subs);
  iop->keep_protein_ids = GetStatus (dlg->keep_protein_ids);
  if (dlg->protein_options == NULL || GetStatus (dlg->update_proteins) == FALSE)
  {
    iop->update_proteins = FALSE;
    iop->truncate_proteins = FALSE;
    iop->extend_proteins3 = FALSE;
    iop->extend_proteins5 = FALSE;
    iop->correct_cds_genes = FALSE;
  }
  else
  {
    iop->update_proteins = TRUE;
    iop->truncate_proteins = GetStatus (dlg->truncate_proteins);
    iop->extend_proteins3 = GetStatus (dlg->extend_proteins3);
    iop->extend_proteins5 = GetStatus (dlg->extend_proteins5);
    iop->correct_cds_genes = GetStatus (dlg->correct_cds_genes);
  }
  return iop;
}

static void DisableIndexerImportFeatureOptions (DialoG d)
{
  IndexerOptionsDialogPtr dlg;

  dlg = (IndexerOptionsDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  Disable (dlg->update_proteins);
  Disable (dlg->protein_options);
  SafeDisable (dlg->keep_protein_ids);
}

static void EnableIndexerImportFeatureOptions (DialoG d)
{
  IndexerOptionsDialogPtr dlg;

  dlg = (IndexerOptionsDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  SafeEnable (dlg->keep_protein_ids);
  Enable (dlg->update_proteins);
  EnableIndexerOptions (dlg->update_proteins);
}

static DialoG IndexerUpdateOptionsDialog (GrouP parent, Boolean is_nuc)
{
  GrouP p;
  IndexerOptionsDialogPtr dlg;
  
  dlg = (IndexerOptionsDialogPtr) MemNew (sizeof (IndexerOptionsDialogData));
  if (dlg == NULL)
  {
    return NULL;
  }
  p = HiddenGroup (parent, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 2, 2);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = IndexerOptionsToDialog;
  dlg->fromdialog = IndexerOptionsFromDialog;
  
  dlg->keep_protein_ids = CheckBox (p, "Keep protein IDs", NULL);
  dlg->add_cit_subs = CheckBox (p, "Add Cit-subs for Updated Sequences", NULL);
  dlg->update_quality_scores = CheckBox (p, "Replace Quality Scores", NULL);
  if (is_nuc)
  {
    dlg->update_proteins = CheckBox (p, "Update Proteins for Updated Sequences", EnableIndexerOptions);
    SetObjectExtra (dlg->update_proteins, dlg, NULL);
    dlg->protein_options = HiddenGroup (p, 1, 0, NULL);
    dlg->truncate_proteins = CheckBox (dlg->protein_options,
                                       "Truncate retranslated proteins at stops",
                                       NULL);
    dlg->extend_proteins3 = CheckBox (dlg->protein_options,
                                      "Extend retranslated proteins without stops",
                                      NULL);
    dlg->extend_proteins5 = CheckBox (dlg->protein_options,
                                      "Extend retranslated proteins without starts",
                                      NULL);
    dlg->correct_cds_genes = CheckBox (dlg->protein_options, "Correct CDS genes", NULL);
  }
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->keep_protein_ids,
                (HANDLE) dlg->add_cit_subs, 
                (HANDLE) dlg->update_quality_scores,
                (HANDLE) dlg->update_proteins,
                (HANDLE) dlg->protein_options,
                NULL);
  EnableIndexerOptions (dlg->update_proteins);
  return (DialoG) p;
}

typedef struct updateoptionsdialog
{
  DIALOG_MESSAGE_BLOCK
  DialoG submitter_opts;
  DialoG indexer_opts;
} UpdateOptionsDialogData, PNTR UpdateOptionsDialogPtr;

static void UpdateOptionsToDialog (DialoG d, Pointer data)
{
  UpdateOptionsDialogPtr dlg;
  UpdateOptionsPtr       uop;
  
  dlg = (UpdateOptionsDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  uop = (UpdateOptionsPtr) data;
  if (uop == NULL)
  {
    PointerToDialog (dlg->submitter_opts, NULL);
    PointerToDialog (dlg->indexer_opts, NULL);
  }
  else
  {
    PointerToDialog (dlg->submitter_opts, uop->submitter_opts);
    PointerToDialog (dlg->indexer_opts, uop->indexer_opts);
  }
}

static Pointer UpdateOptionsFromDialog (DialoG d)
{
  UpdateOptionsDialogPtr dlg;
  UpdateOptionsPtr       uop;
  
  dlg = (UpdateOptionsDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  uop = (UpdateOptionsPtr) MemNew (sizeof (UpdateOptionsData));
  if (uop != NULL)
  {
    uop->submitter_opts = DialogToPointer (dlg->submitter_opts);
    if (dlg->indexer_opts == NULL)
    {
      uop->indexer_opts = NULL;
    }
    else
    {
      uop->indexer_opts = DialogToPointer (dlg->indexer_opts);
    }
  }
  return uop;
}

static void DisableImportFeatureOptions (DialoG d)
{
  UpdateOptionsDialogPtr dlg;

  dlg = (UpdateOptionsDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  DisableSubmitterImportFeatureOptions (dlg->submitter_opts);
  DisableIndexerImportFeatureOptions (dlg->indexer_opts);
}

static void EnableImportFeatureOptions (DialoG d)
{
  UpdateOptionsDialogPtr dlg;

  dlg = (UpdateOptionsDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  EnableSubmitterImportFeatureOptions (dlg->submitter_opts);
  EnableIndexerImportFeatureOptions (dlg->indexer_opts);
}

static void DisableRemoveFeatureOptions (DialoG d)
{
  UpdateOptionsDialogPtr dlg;

  dlg = (UpdateOptionsDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  DisableSubmitterRemoveFeatureOptions (dlg->submitter_opts);
}

static void EnableRemoveFeatureOptions (DialoG d)
{
  UpdateOptionsDialogPtr dlg;

  dlg = (UpdateOptionsDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  EnableSubmitterRemoveFeatureOptions (dlg->submitter_opts);
}

static void 
AdjustUpdateOptionsEnabled 
(DialoG                    d, 
 UpdatePairPtr             upp,
 UpdateAlignmentLengthsPtr ualp,
 Boolean                   is_indexer)
{
  UpdateOptionsDialogPtr dlg;

  dlg = (UpdateOptionsDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }

  AdjustSequenceUpdateOptionsEnabled (dlg->submitter_opts, upp, ualp, is_indexer);    
}

static DialoG 
UpdateOptionsDialog 
(GrouP                parent,
 Boolean              is_nuc, 
 Boolean              indexer,
 Boolean              do_update,
 Nlm_ChangeNotifyProc change_notify,
 Pointer              change_userdata)
{
  UpdateOptionsDialogPtr dlg;
  GrouP                  p;
  
  dlg = (UpdateOptionsDialogPtr) MemNew (sizeof (UpdateOptionsDialogData));
  if (dlg == NULL)
  {
    return NULL;
  }
  
  p = HiddenGroup (parent, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);
  dlg->dialog = (DialoG) p;
  dlg->todialog = UpdateOptionsToDialog;
  dlg->fromdialog = UpdateOptionsFromDialog;
  
  dlg->submitter_opts = SubmitterUpdateOptionsDialog (p, indexer, do_update, change_notify, change_userdata);
  if (indexer)
  {
    dlg->indexer_opts = IndexerUpdateOptionsDialog (p, is_nuc);
  }
  else
  {
    dlg->indexer_opts = NULL;
  }
  
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->submitter_opts,
                              (HANDLE) dlg->indexer_opts,
                              NULL);
  
  return (DialoG) p;
}

static ValNodePtr GetNthValNode (ValNodePtr list, Int4 nth)
{
  if (nth < 0)
  {
    return NULL;
  }
  
  while (nth > 0 && list != NULL)
  {
    list = list->next;
    nth --;
  }
  return list;
}

static ValNodePtr ExtractNthValNode (ValNodePtr PNTR list, Int4 nth)
{
  ValNodePtr prev = NULL, this_vnp;
  if (nth < 0)
  {
    return NULL;
  }
  
  this_vnp = *list;
  while (nth > 0 && this_vnp != NULL)
  {
    prev = this_vnp;
    this_vnp = this_vnp->next;
    nth --;
  }
  
  if (this_vnp != NULL)
  {
    if (prev == NULL)
    {
      *list = (*list)->next;
    }
    else
    {
      prev->next = this_vnp->next;
    }
    this_vnp->next = NULL;
  }
  return this_vnp;
  
}

extern void 
ListBioseqsInSeqEntry 
(SeqEntryPtr     sep, 
 Boolean         is_na,
 Int4Ptr         seq_num, 
 ValNodePtr PNTR bioseq_list)
{
  BioseqPtr    bsp;
  BioseqSetPtr bssp;
  
  if (sep == NULL || bioseq_list == NULL || seq_num == NULL)
  {
    return;
  }
  if (IS_Bioseq (sep) && sep->data.ptrvalue != NULL)
  {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (ISA_na (bsp->mol))
    {
      if (is_na)
      {
        ValNodeAddPointer (bioseq_list, *seq_num, bsp);
        (*seq_num)++;
      }
    }
    else if (!is_na)
    {
      ValNodeAddPointer (bioseq_list, *seq_num, bsp);
      (*seq_num)++;
    } 
  }
  else if (IS_Bioseq_set (sep) && sep->data.ptrvalue != NULL)
  {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    ListBioseqsInSeqEntry (bssp->seq_set, is_na, seq_num, bioseq_list);
  }
  
  ListBioseqsInSeqEntry (sep->next, is_na, seq_num, bioseq_list);
}

/* This function compares the text from a local ID against the 
 * report string from non-local IDs in sip_list, useful when
 * comparing values from a file in which the user did not specify
 * the version or the gb| in the sequence ID.
 */
static Boolean RelaxedSeqIdIn (SeqIdPtr sip, SeqIdPtr sip_list)
{
  SeqIdPtr sip_next;
  Char     id_txt1 [128], id_txt2 [128];
  CharPtr  ptr;
  Int4     len;
  
  if (sip == NULL || sip_list == NULL || sip->choice != SEQID_LOCAL)
  {
    return FALSE;
  }
  
  SeqIdWrite (sip, id_txt1, PRINTID_REPORT, sizeof (id_txt1) - 1);
  
  while (sip_list != NULL)
  {
    if (sip_list->choice != SEQID_LOCAL)
    {
      sip_next = sip_list->next;
      sip_list->next = NULL;
      SeqIdWrite (sip_list, id_txt2, PRINTID_REPORT, sizeof (id_txt1) - 1);
      sip_list->next = sip_next;
      if (StringCmp (id_txt1, id_txt2) == 0)
      {
        return TRUE;
      }
      ptr = StringChr (id_txt2, '.');
      if (ptr != NULL)  /* ID in list has version */
      {
        len = StringLen (id_txt1);
        if (len == ptr - id_txt2 && StringNCmp (id_txt1, id_txt2, len) == 0)
        {
          return TRUE;
        }
      }
    }
    sip_list = sip_list->next;
  }
  return FALSE;
}

static BioseqPtr FindBioseqInList (ValNodePtr bioseq_list, SeqIdPtr sip, Int4Ptr position)
{
  ValNodePtr vnp;
  BioseqPtr  bsp = NULL;
  Int4       vnp_pos;
  
  if (position != NULL)
  {
    *position = -1;
  }
  if (bioseq_list == NULL)
  {
    return NULL;
  }
  
  for (vnp = bioseq_list, vnp_pos = 0;
       vnp != NULL && bsp == NULL;
       vnp = vnp->next, vnp_pos++)
  {
    bsp = (BioseqPtr) vnp->data.ptrvalue;
    if (SeqIdIn (sip, bsp->id) || RelaxedSeqIdIn (sip, bsp->id))
    {
      if (position != NULL)
      {
        *position = vnp_pos;
      }
    }
    else
    {
      bsp = NULL;
    }
  }
  return bsp;
}

typedef struct previewsequenceselectiondialog
{
  DIALOG_MESSAGE_BLOCK
  DoC                      doc;  
  Int4                     sequence_row;
  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
  ValNodePtr               sequence_list;
  ParData                  ParFmt;
  ColData                  ColFmt;
  
} PreviewSequenceSelectionDialogData, PNTR PreviewSequenceSelectionDialogPtr;

static void SelectionToPreviewSequenceSelectionDialog (DialoG d, Pointer userdata)
{
  PreviewSequenceSelectionDialogPtr dlg;
  SeqIdPtr                          sip;
  Int4                              seq_num, match_num;
  ValNodePtr                        vnp;

  dlg = (PreviewSequenceSelectionDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  sip = (SeqIdPtr) userdata;
  if (sip == NULL)
  {
    dlg->sequence_row = -1;  
  }
  else
  {
    match_num = -1;
    for (seq_num = 1, vnp = dlg->sequence_list;
         vnp != NULL && match_num < 0;
         seq_num++, vnp = vnp->next)
    {
      if (SeqIdComp (sip, vnp->data.ptrvalue))
      {
        match_num = seq_num;
      }
    }
    dlg->sequence_row = match_num;
  }
  InvalDocRows (dlg->doc, 0, 0, 0);
  
  UpdateDocument(dlg->doc, 0, 0);
  
  if (dlg->change_notify != NULL)
  {
    (dlg->change_notify) (dlg->change_userdata);
  }
}

static Pointer SelectionFromPreviewSequenceSelectionDialog (DialoG d)
{
  PreviewSequenceSelectionDialogPtr dlg;
  SeqIdPtr                          sip = NULL;
  ValNodePtr                        vnp;

  dlg = (PreviewSequenceSelectionDialogPtr) GetObjectExtra (d);
  if (dlg != NULL && dlg->sequence_row > 0)
  {
    vnp = GetNthValNode (dlg->sequence_list, dlg->sequence_row - 1);
    if (vnp != NULL)
    {
      sip = vnp->data.ptrvalue;
    }
  }
  return sip;
}

static void ResetSequenceList (DialoG d, ValNodePtr sequence_list)
{
  PreviewSequenceSelectionDialogPtr dlg;
  Char                              id_txt [MAX_ID_LEN];
  SeqIdPtr                          sip;
  
  dlg = (PreviewSequenceSelectionDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  dlg->sequence_list = ValNodeFree (dlg->sequence_list);
  dlg->sequence_list = sequence_list;
  Reset (dlg->doc);
  
  while (sequence_list != NULL)
  {
    sip = (SeqIdPtr) sequence_list->data.ptrvalue;
    if (sip != NULL)
    {
      /* add to sequence_selector doc */
      SeqIdWrite (SeqIdFindBest (sip, SEQID_GENBANK), id_txt, PRINTID_REPORT, sizeof (id_txt) - 1);
  	  AppendText (dlg->doc, id_txt, &(dlg->ParFmt), &(dlg->ColFmt), programFont);  	  
    }
    sequence_list = sequence_list->next;
  }
  InvalDocRows (dlg->doc, 0, 0, 0);
  dlg->sequence_row = -1;  
}

static void SelectPreviewSequence (DoC d, PoinT pt)
{
  Int2      item, row, prevrow;
  PreviewSequenceSelectionDialogPtr dlg;
  
  dlg = (PreviewSequenceSelectionDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  
  MapDocPoint (d, pt, &item, &row, NULL, NULL);
  if (item > 0 && row > 0) {
    prevrow = dlg->sequence_row;
    dlg->sequence_row = item;
    if (item != prevrow)
    {
      if (prevrow != -1)
      {
        InvalDocRows (d, prevrow, 1, 1);
      }
      InvalDocRows (d, item, 1, 1);
      if (dlg->change_notify != NULL)
      {
        (dlg->change_notify) (dlg->change_userdata);
      }
    }
  }
}

static Boolean PreviewSequenceHighlight (DoC doc, Int2 item, Int2 row, Int2 col)
{
  PreviewSequenceSelectionDialogPtr dlg;
  
  dlg = (PreviewSequenceSelectionDialogPtr) GetObjectExtra (doc);
  if (dlg == NULL) return FALSE;
  
  if (item == dlg->sequence_row) return TRUE;
  return FALSE;
}

static void CleanupPreviewSequenceSelectionDialog (GraphiC g, Pointer data)
{
  PreviewSequenceSelectionDialogPtr dlg;

  dlg = (PreviewSequenceSelectionDialogPtr) data;
  if (dlg != NULL)
  {
    dlg->sequence_list = ValNodeFree (dlg->sequence_list);
  }
  StdCleanupExtraProc (g, data);
}

static void SetPreviewSequenceSelectionByPosition (DialoG d, Int4 seq_pos)
{
  PreviewSequenceSelectionDialogPtr dlg;
  Int4                              prevrow, curr_scroll, max_scroll;
  BaR                               sb;

  dlg = (PreviewSequenceSelectionDialogPtr) GetObjectExtra (d);
  if (dlg == NULL || seq_pos < 0)
  {
    return;
  }
  
  prevrow = dlg->sequence_row;
  dlg->sequence_row = seq_pos + 1;
  if (dlg->sequence_row != prevrow)
  {
    if (prevrow != -1)
    {
      InvalDocRows (dlg->doc, prevrow, 1, 1);
    }
    InvalDocRows (dlg->doc, seq_pos + 1, 1, 1);
    sb = GetSlateVScrollBar( (SlatE) dlg->doc);
    if (sb != NULL)
    {
      curr_scroll = GetBarValue (sb);
      max_scroll = GetBarMax (sb);
      if (max_scroll > seq_pos)
      {
        SetBarValue (sb, seq_pos);
      }
      else
      {
        SetBarValue (sb, max_scroll);
      }
    }
    if (dlg->change_notify != NULL)
    {
      (dlg->change_notify) (dlg->change_userdata);
    }
  }
}

static DialoG 
PreviewSequenceSelectionDialog 
(GrouP                    parent, 
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  PreviewSequenceSelectionDialogPtr dlg;
  GrouP                             p;
  RecT                              r;
  
  dlg = (PreviewSequenceSelectionDialogPtr) MemNew (sizeof (PreviewSequenceSelectionDialogData));
  if (dlg == NULL)
  {
    return NULL;
  }
  
  p = HiddenGroup (parent, 1, 0, NULL);
  SetObjectExtra (p, dlg, CleanupPreviewSequenceSelectionDialog);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = SelectionToPreviewSequenceSelectionDialog;
  dlg->fromdialog = SelectionFromPreviewSequenceSelectionDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;
  dlg->sequence_list = NULL;
  
  dlg->doc = DocumentPanel (p, stdCharWidth * 10, stdLineHeight * 5);
  SetObjectExtra (dlg->doc, dlg, NULL);
  SetDocProcs (dlg->doc, SelectPreviewSequence, NULL, NULL, NULL);
  SetDocShade (dlg->doc, NULL, NULL, PreviewSequenceHighlight, NULL);
  
  /* initialize document paragraph format */
  dlg->ParFmt.openSpace = FALSE;
  dlg->ParFmt.keepWithNext = FALSE;
  dlg->ParFmt.keepTogether = FALSE;
  dlg->ParFmt.newPage = FALSE;
  dlg->ParFmt.tabStops = FALSE;
  dlg->ParFmt.minLines = 0;
  dlg->ParFmt.minHeight = 0;
  
  /* initialize document column format */
  ObjectRect (dlg->doc, &r);
  InsetRect (&r, 4, 4);
  dlg->ColFmt.pixWidth = r.right - r.left;
  dlg->ColFmt.pixInset = 0;
  dlg->ColFmt.charWidth = 80;
  dlg->ColFmt.charInset = 0;
  dlg->ColFmt.font = NULL;
  dlg->ColFmt.just = 'l';
  dlg->ColFmt.wrap = TRUE;
  dlg->ColFmt.bar = FALSE;
  dlg->ColFmt.underline = FALSE;
  dlg->ColFmt.left = FALSE;
  dlg->ColFmt.last = TRUE;

  
  return (DialoG) p;
}

typedef struct unmatchedsequencedialog
{
  DIALOG_MESSAGE_BLOCK
  DoC      doc;
  ParData  ParFmt;
  ColData  ColFmt;
  Int4     dlg_height;
} UnmatchedSequenceDialogData, PNTR UnmatchedSequenceDialogPtr;

static void ListToUnmatchedSequenceDialog (DialoG d, Pointer userdata)
{
  UnmatchedSequenceDialogPtr dlg;
  ValNodePtr                 bioseq_list;
  BioseqPtr                  bsp;
  Char                       id_txt [MAX_ID_LEN];
  
  dlg = (UnmatchedSequenceDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  Reset (dlg->doc);
  
  bioseq_list = (ValNodePtr) userdata;
  if (bioseq_list == NULL)
  {
    return;
  }

  while (bioseq_list != NULL)
  {
    bsp = (BioseqPtr) bioseq_list->data.ptrvalue;
    if (bsp != NULL)
    {
      /* add to sequence_selector doc */
      SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_txt, PRINTID_REPORT, sizeof (id_txt) - 1);
  	  AppendText (dlg->doc, id_txt, &(dlg->ParFmt), &(dlg->ColFmt), programFont);  	  
    }
    bioseq_list = bioseq_list->next;
  }
  InvalDocRows (dlg->doc, 0, 0, 0);
  
}

static DialoG 
UnmatchedSequenceDialog 
(GrouP           parent, 
 Nlm_BtnActnProc map_btn_proc,
 Pointer         change_userdata)
{
  UnmatchedSequenceDialogPtr dlg;
  GrouP                      p;
  RecT                       r;
  PrompT                     ppt = NULL;
  ButtoN                     b = NULL;
  
  dlg = (UnmatchedSequenceDialogPtr) MemNew (sizeof (UnmatchedSequenceDialogData));
  if (dlg == NULL)
  {
    return NULL;
  }
  
  p = HiddenGroup (parent, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = ListToUnmatchedSequenceDialog;
  
  if (map_btn_proc == NULL)
  {
  	ppt = StaticPrompt (p, "No Updates", 0, 0, programFont, 'l');
  }
  else
  {
    ppt = StaticPrompt (p, "Unmatched Sequences", 0, 0, programFont, 'l');
  }
  dlg->doc = DocumentPanel (p, stdCharWidth * 10, stdLineHeight * 5);

  /* initialize document paragraph format */
  dlg->ParFmt.openSpace = FALSE;
  dlg->ParFmt.keepWithNext = FALSE;
  dlg->ParFmt.keepTogether = FALSE;
  dlg->ParFmt.newPage = FALSE;
  dlg->ParFmt.tabStops = FALSE;
  dlg->ParFmt.minLines = 0;
  dlg->ParFmt.minHeight = 0;
  
  /* initialize document column format */
  ObjectRect (dlg->doc, &r);
  dlg->dlg_height = r.bottom - r.top;
  InsetRect (&r, 4, 4);
  dlg->ColFmt.pixWidth = r.right - r.left;
  dlg->ColFmt.pixInset = 0;
  dlg->ColFmt.charWidth = 80;
  dlg->ColFmt.charInset = 0;
  dlg->ColFmt.font = NULL;
  dlg->ColFmt.just = 'l';
  dlg->ColFmt.wrap = TRUE;
  dlg->ColFmt.bar = FALSE;
  dlg->ColFmt.underline = FALSE;
  dlg->ColFmt.left = FALSE;
  dlg->ColFmt.last = TRUE;


  if (map_btn_proc != NULL)
  {
    /* add button for loading map */
    b = PushButton (p, "Load Map", map_btn_proc);
    SetObjectExtra (b, change_userdata, NULL);
    ObjectRect (b, &r);
    dlg->dlg_height += 10 + r.bottom - r.top;
  }
  
  AlignObjects (ALIGN_CENTER, (HANDLE) ppt, (HANDLE) dlg->doc, (HANDLE) b, NULL);
  
  return (DialoG) p;
}

static void 
ChangePreviewSequenceSelectionDialogHeights 
(DialoG orig_sel,
 DialoG unmatched,
 DialoG no_updates,
 ButtoN remove_identical_btn,
 Boolean has_identical,
 WindoW match_win,
 Boolean show_unmatched,
 Boolean show_no_updates)
{
  PreviewSequenceSelectionDialogPtr orig_sel_dlg;
  UnmatchedSequenceDialogPtr        unmatched_dlg, no_updates_dlg;
  RecT                              pictures_r, unmatched_r, list_r;
  RecT                              no_updates_r, rem_ident_r;
  Int4                              no_updates_height;
  BaR                               sb;
  Int4                              updates_list_bottom;
  Int4                              spacing = 0, rem_ident_height = 0;

  orig_sel_dlg = (PreviewSequenceSelectionDialogPtr) GetObjectExtra (orig_sel);
  unmatched_dlg = (UnmatchedSequenceDialogPtr) GetObjectExtra (unmatched);
  no_updates_dlg = (UnmatchedSequenceDialogPtr) GetObjectExtra (no_updates);
  if (orig_sel_dlg == NULL || unmatched_dlg == NULL 
      || no_updates_dlg == NULL)
  {
    return;
  }

  ObjectRect (orig_sel_dlg->doc, &list_r);
  ObjectRect (unmatched_dlg->doc, &unmatched_r);
  ObjectRect (no_updates_dlg->doc, &no_updates_r);
  if (match_win == NULL)
  {
    pictures_r.top = list_r.top;
    pictures_r.left = list_r.right;
    pictures_r.right = list_r.right;
    pictures_r.bottom = no_updates_r.bottom;
  }
  else
  {
    ObjectRect (match_win, &pictures_r);
  }
  
  if (remove_identical_btn != NULL && has_identical)
  {
    ObjectRect (remove_identical_btn, &rem_ident_r);  
    rem_ident_height = rem_ident_r.bottom - rem_ident_r.top;
    spacing = 5;
  }
  
  /* updates list is on top, 
   * followed by remove_identical_btn,
   * followed by unmatched list, 
   * followed by no updates list */
  updates_list_bottom = pictures_r.bottom;

  if (show_no_updates)
  {
    /* align no updates list with bottom of preview pictures */
    no_updates_height = no_updates_r.bottom - no_updates_r.top;
    no_updates_r.bottom = updates_list_bottom;
    no_updates_r.top = no_updates_r.bottom - no_updates_height;
    sb = GetSlateVScrollBar( (SlatE) no_updates_dlg->doc);
    if (sb != NULL)
    {
      /* because the doc has a vertical scroll bar and the set position subtracts the
       * width of the scroll bar before positioning the list, must add the width of
       * the scroll bar to the rightt.
       */
      no_updates_r.right += Nlm_vScrollBarWidth;
    }
    SetPosition (no_updates, &no_updates_r);
    updates_list_bottom -= (no_updates_dlg->dlg_height + 10);
    Show (no_updates);
  }
  else
  {
  	Hide (no_updates);
  }

  if (show_unmatched)
  {
    /* align unmatched list with bottom of updates list or preview pictures */
    unmatched_r.bottom = updates_list_bottom ;
    unmatched_r.top = unmatched_r.bottom - unmatched_dlg->dlg_height;
  
    sb = GetSlateVScrollBar( (SlatE) unmatched_dlg->doc);
    if (sb != NULL)
    {
      /* because the doc has a vertical scroll bar and the set position subtracts the
       * width of the scroll bar before positioning the list, must add the width of
       * the scroll bar to the rightt.
       */
      unmatched_r.right += Nlm_vScrollBarWidth;
    }
    SetPosition (unmatched, &unmatched_r);
    updates_list_bottom -= (unmatched_dlg->dlg_height + 10);
    Show (unmatched);
  }
  else
  {
  	Hide (unmatched);
  }
  
  /* set position of sequence list - align top with pictures,
   * align bottom with top of unmatched (if shown) or
   * bottom of pictures
   */
  list_r.top = pictures_r.top;
  list_r.bottom = updates_list_bottom - spacing - rem_ident_height;
    
  sb = GetSlateVScrollBar( (SlatE) orig_sel_dlg->doc);
  if (sb != NULL)
  {
    /* because the doc has a vertical scroll bar and the set position subtracts the
     * width of the scroll bar before positioning the list, must add the width of
     * the scroll bar to the rightt.
     */
    list_r.right += Nlm_vScrollBarWidth;
  }
  SetPosition (orig_sel_dlg->doc, &list_r); 

  if (remove_identical_btn != NULL)
  {
    if (has_identical)
    {
      /* set position of remove identical btn */
      rem_ident_r.top = updates_list_bottom - rem_ident_height;
      rem_ident_r.bottom = rem_ident_r.top + rem_ident_height;
      SetPosition (remove_identical_btn, &rem_ident_r);
    }
    else
    {
      /* hide button */
      SafeHide (remove_identical_btn);
    }
  }
}

typedef struct multisequenceupdatepreviewdialog
{
  DIALOG_MESSAGE_BLOCK
  
  DialoG               update_list_dlg;
  DialoG               preview_dlg;
  DialoG               unmatched_list_dlg;
  DialoG               no_updates_list_dlg;
  DialoG               extend_preview_dlg;
  ButtoN               remove_identical_btn;
  WindoW               match_win;
  ValNodePtr           orig_bioseq_list;
  ValNodePtr           update_bioseq_list;
  Boolean              is_na;
  Nlm_ChangeNotifyProc change_notify;
  Pointer              change_userdata;
  
} MultiSequenceUpdatePreviewDialogData, PNTR MultiSequenceUpdatePreviewDialogPtr;

typedef struct multisequenceupdate
{
  ValNodePtr  orig_bioseq_list;
  ValNodePtr  update_bioseq_list;
  ValNodePtr  unmatched_updates_list;
  ValNodePtr  no_updates_list;
} MultiSequenceUpdateData, PNTR MultiSequenceUpdatePtr;

static void SelectSequenceForUpdatePreview (Pointer userdata)
{
  MultiSequenceUpdatePreviewDialogPtr dlg;
  UpdatePairData                      upd;
  ValNodePtr                          update_vnp;
  SeqIdPtr                            sip;
  Int4                                orig_pos = -1;
  
  dlg = (MultiSequenceUpdatePreviewDialogPtr) userdata;
  if (dlg == NULL)
  {
    return;
  }

  sip = (SeqIdPtr) DialogToPointer (dlg->dialog);
  
  upd.orig_bsp = FindBioseqInList (dlg->orig_bioseq_list, sip, &orig_pos);
  
  update_vnp = GetNthValNode (dlg->update_bioseq_list, orig_pos);
  if (update_vnp != NULL)
  {
    upd.update_bsp = update_vnp->data.ptrvalue;
  }
  else
  {
    upd.update_bsp = NULL;
  }                                         

  SeqMgrReplaceInBioseqIndex (upd.orig_bsp);

  upd.revcomp = FALSE;
  upd.salp = Sequin_GlobalAlign2Seq (upd.orig_bsp, upd.update_bsp, &(upd.revcomp));
  if (upd.revcomp)
  {
    BioseqRevComp (upd.update_bsp);
    ReverseBioseqFeatureStrands (upd.update_bsp);
    SeqMgrReplaceInBioseqIndex (upd.update_bsp);
  }
    
  PointerToDialog (dlg->preview_dlg, &upd);  
  PointerToDialog (dlg->extend_preview_dlg, &upd);

  SeqMgrReplaceInBioseqIndex (upd.orig_bsp);
  
  if (dlg->change_notify != NULL)
  {
    (dlg->change_notify) (dlg->change_userdata);
  }
}

static ValNodePtr SeqIdListFromBioseqList (ValNodePtr bioseq_list)
{
  ValNodePtr sip_list = NULL;
  Int4       seq_num = 0;
  BioseqPtr  bsp;
  
  while (bioseq_list != NULL)
  {
    if (bioseq_list->data.ptrvalue == NULL)
    {
      ValNodeAddPointer (&sip_list, seq_num, NULL);
    }
    else
    {
      bsp = (BioseqPtr) bioseq_list->data.ptrvalue;
      ValNodeAddPointer (&sip_list, seq_num, bsp->id);
    }
    seq_num++;
    bioseq_list = bioseq_list->next;
  }
  return sip_list;
}

static Boolean HasIdenticalUpdates 
(ValNodePtr orig_list,
 ValNodePtr update_list)
{
  BioseqPtr  orig_bsp, update_bsp;
  
  while (orig_list != NULL && update_list != NULL)
  {
    orig_bsp = orig_list->data.ptrvalue;
    update_bsp = update_list->data.ptrvalue;
    
    if (orig_bsp != NULL && update_bsp != NULL 
        && AreSequenceResiduesIdentical (orig_bsp, update_bsp))
    {
      return TRUE;
    }
    
    orig_list = orig_list->next;
    update_list = update_list->next;
  }
  
  return FALSE;
}

static void DataToMultiSequenceUpdatePreview (DialoG d, Pointer data)
{
  MultiSequenceUpdatePreviewDialogPtr dlg;
  MultiSequenceUpdatePtr              msup;
  BioseqPtr                           first_bsp;
  ValNodePtr                          unmatched_updates_list = NULL; 
  ValNodePtr                          no_updates_list = NULL; 
  Boolean                             has_identical = FALSE;
  
  dlg = (MultiSequenceUpdatePreviewDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
    
  msup = (MultiSequenceUpdatePtr) data;
  if (msup == NULL)
  {
    dlg->orig_bioseq_list = NULL;
    dlg->update_bioseq_list = NULL;
  }
  else
  {
    dlg->orig_bioseq_list = msup->orig_bioseq_list;
    dlg->update_bioseq_list = msup->update_bioseq_list;
    unmatched_updates_list = msup->unmatched_updates_list;
    no_updates_list = msup->no_updates_list;
  }

  ResetSequenceList (dlg->update_list_dlg, SeqIdListFromBioseqList(dlg->orig_bioseq_list));
  PointerToDialog (dlg->unmatched_list_dlg, unmatched_updates_list);
  PointerToDialog (dlg->no_updates_list_dlg, no_updates_list);
  
  if (dlg->remove_identical_btn != NULL)
  {
    has_identical = HasIdenticalUpdates (msup->orig_bioseq_list, msup->update_bioseq_list);
  }
 
  ChangePreviewSequenceSelectionDialogHeights (dlg->update_list_dlg, 
                                               dlg->unmatched_list_dlg,
                                               dlg->no_updates_list_dlg,
                                               dlg->remove_identical_btn,
                                               has_identical,
                                               dlg->match_win,
                                               (msup->unmatched_updates_list != NULL),
                                               (msup->no_updates_list != NULL));

  /* select first sequence */  
  if (dlg->orig_bioseq_list != NULL && dlg->orig_bioseq_list->data.ptrvalue != NULL)
  {
    first_bsp = (BioseqPtr) dlg->orig_bioseq_list->data.ptrvalue;

    PointerToDialog (dlg->update_list_dlg, first_bsp->id);
    SelectSequenceForUpdatePreview (dlg); 
  }
  else
  {
    if (dlg->change_notify != NULL)
    {
      (dlg->change_notify) (dlg->change_userdata);
    }    
  }
}

/* returns SeqID currently selected in preview list */
static Pointer SelectionFromMultiSequenceUpdatePreview (DialoG d)
{
  MultiSequenceUpdatePreviewDialogPtr dlg;
  ValNodePtr                          orig_vnp, update_vnp;
  BioseqPtr                           orig_bsp;
  SeqIdPtr                            return_sip = NULL;
  
  dlg = (MultiSequenceUpdatePreviewDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  
  if (dlg->update_list_dlg == NULL)
  {
    for (orig_vnp = dlg->orig_bioseq_list, update_vnp = dlg->update_bioseq_list;
         orig_vnp != NULL && update_vnp != NULL && return_sip == NULL;
         orig_vnp = orig_vnp->next, update_vnp = update_vnp->next)
    {
      if (orig_vnp->data.ptrvalue == NULL || update_vnp->data.ptrvalue == NULL)
      {
        continue;
      }
      orig_bsp = (BioseqPtr) orig_vnp->data.ptrvalue;
      return_sip = SeqIdDup (orig_bsp->id);
    }
  }
  else
  {
    return_sip = DialogToPointer (dlg->update_list_dlg);
  }
  return return_sip;
}

static void SelectUpdatePreviewSequenceByPosition (DialoG d, Int4 seq_pos)
{
  MultiSequenceUpdatePreviewDialogPtr dlg; 

  dlg = (MultiSequenceUpdatePreviewDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  SetPreviewSequenceSelectionByPosition (dlg->update_list_dlg, seq_pos);  
}

static UpdatePairPtr GetCurrentUpdatePair (DialoG d)
{
  MultiSequenceUpdatePreviewDialogPtr dlg; 

  dlg = (MultiSequenceUpdatePreviewDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  
  return (UpdatePairPtr) DialogToPointer (dlg->preview_dlg);
}

static void AddExtendDialogToMultiSequenceUpdatePreview (DialoG msup, DialoG ext)
{
  MultiSequenceUpdatePreviewDialogPtr dlg; 

  dlg = (MultiSequenceUpdatePreviewDialogPtr) GetObjectExtra (msup);
  if (dlg != NULL)
  {
    dlg->extend_preview_dlg = ext;
  }
}

static void RemoveIdenticalUpdates (ButtoN b);

static void SetMultiSequenceUpdatePreviewMatchWindow (DialoG d, WindoW w)
{
  MultiSequenceUpdatePreviewDialogPtr dlg; 

  dlg = (MultiSequenceUpdatePreviewDialogPtr) GetObjectExtra (d);
  if (dlg != NULL)
  {
    dlg->match_win = w;
  }
}

static void SetIgnoreAndExtendForMultiSequencePreview (DialoG d, Boolean ignore_alignment, Boolean extend5)
{
  MultiSequenceUpdatePreviewDialogPtr dlg; 

  dlg = (MultiSequenceUpdatePreviewDialogPtr) GetObjectExtra (d);
  if (dlg != NULL)
  {
    SetIgnoreAndExtendForPreviewPictures (dlg->preview_dlg, ignore_alignment, extend5);
  }
}

static DialoG 
MultiSequenceUpdatePreview 
(GrouP                parent, 
 Boolean              is_na,
 Boolean              do_update,
 Nlm_ChangeNotifyProc change_notify,
 Pointer              change_userdata,
 Nlm_BtnActnProc      map_btn_proc,
 Boolean              is_indexer)
{
  MultiSequenceUpdatePreviewDialogPtr dlg; 
  GrouP                               p, k;
  
  dlg = (MultiSequenceUpdatePreviewDialogPtr) MemNew (sizeof (MultiSequenceUpdatePreviewDialogData));
  if (dlg == NULL)
  {
    return NULL;
  }

  if (do_update)
  {
    p = HiddenGroup (parent, 0, 1, NULL);
  }
  else
  {
    p = HiddenGroup (parent, -1, 0, NULL);
  }
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = DataToMultiSequenceUpdatePreview;
  dlg->fromdialog = SelectionFromMultiSequenceUpdatePreview;
  
  dlg->is_na = is_na;
  
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  if (is_indexer)
  {
    k = HiddenGroup (p, 0, 4, NULL);
    /* note - the ValNodeSelectionDialog will free the seq_id list when done */  
    dlg->update_list_dlg = PreviewSequenceSelectionDialog (k, SelectSequenceForUpdatePreview,
                                                           dlg);
    if (do_update && is_indexer)
    {
      dlg->remove_identical_btn = PushButton (k, "Remove Identical Updates", RemoveIdenticalUpdates);
      SetObjectExtra (dlg->remove_identical_btn, change_userdata, NULL);      
    }

    dlg->unmatched_list_dlg = UnmatchedSequenceDialog (k, map_btn_proc, change_userdata);
    
    dlg->no_updates_list_dlg = UnmatchedSequenceDialog (k, NULL, NULL);
  }
  
  dlg->extend_preview_dlg = NULL;
  if (do_update)
  {
    dlg->preview_dlg = UpdatePreviewDialog (p, is_indexer);
    dlg->match_win = (WindoW) dlg->preview_dlg;
  }
  else
  {
    dlg->preview_dlg = NULL;
    dlg->match_win = NULL;
  }
  
  return (DialoG) p;
}

static void AddUniqueUpdateSequenceIDs (SeqEntryPtr sep)
{
  BioseqPtr bsp;
  BioseqSetPtr bssp;
  
  if (sep == NULL)
  {
    return;
  }
  else if (IS_Bioseq (sep))
  {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp != NULL && bsp->id == NULL)
    {
      bsp->id = MakeUniqueSeqID ("UpdateSequence");
    }
  }
  else if (IS_Bioseq_set (sep))
  {
    /* we could add IDs to segmented sets, but maybe we should just remove them? */
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL)
    {
      AddUniqueUpdateSequenceIDs (bssp->seq_set);
    }
  }
  
  AddUniqueUpdateSequenceIDs (sep->next);
}

static Boolean SeqIdListsOverlap (SeqIdPtr sip1, SeqIdPtr sip2)
{ 
  SeqIdPtr sip_tmp, sip_next;
  Char     tmp_id_str [MAX_ID_LEN + 5];
  Boolean  rval = FALSE;
  
  while (sip1 != NULL && !rval)
  {
    if (SeqIdIn (sip1, sip2))
    {
      rval = TRUE;
    }
    else if (sip1->choice == SEQID_LOCAL)
    {
      /* check to see if user just forgot to put "gb|" at the front of the IDs */
      sip_next = sip1->next;
      sip1->next = NULL;
      StringCpy (tmp_id_str, "gb|");
      SeqIdWrite (sip1, tmp_id_str + 3, PRINTID_REPORT, sizeof (tmp_id_str) - 4);     
      sip_tmp = MakeSeqID (tmp_id_str);
      if (SeqIdIn (sip_tmp, sip2))
      {
        rval = TRUE;
      }
      sip_tmp = SeqIdFree (sip_tmp);
      sip1->next = sip_next;
    }
    sip1 = sip1->next;
  }
  return rval;
}

static ValNodePtr ShuffleUpdateBioseqList (ValNodePtr PNTR update_bioseq_list, ValNodePtr orig_bioseq_list)
{
  ValNodePtr unmatched_list = NULL;
  ValNodePtr orig_vnp, update_vnp;
  ValNodePtr new_update_list = NULL;
  Int4       bsp_pos = 0, update_pos, pos;
  BioseqPtr  orig_bsp, update_bsp;
  
  if (update_bioseq_list == NULL || *update_bioseq_list == NULL)
  {
    return NULL;
  }
  else if (orig_bioseq_list == NULL)
  {
    unmatched_list = *update_bioseq_list;
    *update_bioseq_list = NULL;
  }

  for (orig_vnp = orig_bioseq_list; orig_vnp != NULL; orig_vnp = orig_vnp->next)
  {
    if (orig_vnp->data.ptrvalue == NULL)
    {
      ValNodeAddPointer (&new_update_list, bsp_pos, NULL);
      bsp_pos ++;
      continue;
    }
    orig_bsp = (BioseqPtr) orig_vnp->data.ptrvalue;
    update_pos = -1;
    for (update_vnp = *update_bioseq_list, pos = 0;
         update_vnp != NULL && update_pos < 0; 
         update_vnp = update_vnp->next, pos++)
    {
      if (update_vnp->data.ptrvalue != NULL)
      {
        update_bsp = (BioseqPtr) update_vnp->data.ptrvalue;
        if (SeqIdListsOverlap (update_bsp->id, orig_bsp->id))
        {
          update_pos = pos;
        }
      }
    }
    if (update_pos >= 0)
    {
      update_vnp = ExtractNthValNode (update_bioseq_list, update_pos);
      update_vnp->choice = bsp_pos;
      ValNodeLink (&new_update_list, update_vnp);
    }
    else
    {
      ValNodeAddPointer (&new_update_list, bsp_pos, NULL);
    }
    bsp_pos++;
  }
  
  unmatched_list = *update_bioseq_list;
  *update_bioseq_list = new_update_list;
  
  /* renumber unmatched_list */
  for (update_vnp = unmatched_list, update_pos = 0;
       update_vnp != NULL;
       update_vnp = update_vnp->next, update_pos++)
  {
    update_vnp->choice = update_pos;
  }
  
  return unmatched_list;
}

static ValNodePtr ExtractSequencesWithoutUpdates (ValNodePtr PNTR orig_bioseq_list, ValNodePtr PNTR update_bioseq_list)
{
  ValNodePtr orig_prev = NULL, update_prev = NULL;
  ValNodePtr orig_next = NULL, update_next = NULL;
  ValNodePtr orig_vnp, update_vnp;
  ValNodePtr no_update_list = NULL;
  Int4       seq_num;
  
  if (orig_bioseq_list == NULL || update_bioseq_list == NULL
      || *orig_bioseq_list == NULL || *update_bioseq_list == NULL)
  {
    return NULL;
  }
  
  orig_vnp = *orig_bioseq_list;
  update_vnp = *update_bioseq_list;
  
  while (orig_vnp != NULL && update_vnp != NULL)
  {
    orig_next = orig_vnp->next;
    update_next = update_vnp->next;
    if (orig_vnp->data.ptrvalue == NULL || update_vnp->data.ptrvalue == NULL)
    {
      if (orig_prev == NULL || update_prev == NULL)
      {
        *orig_bioseq_list = orig_vnp->next;
        *update_bioseq_list = update_vnp->next;
      }
      else
      {
        orig_prev->next = orig_vnp->next;
        update_prev->next = update_vnp->next;
      }
      orig_vnp->next = NULL;
      update_vnp->next = NULL;
      if (orig_vnp->data.ptrvalue == NULL)
      {
        ValNodeFree (orig_vnp);
      }
      else
      {
      	ValNodeLink (&no_update_list, orig_vnp);
      }

      ValNodeFree (update_vnp);
    }
    else
    {
      orig_prev = orig_vnp;
      update_prev = update_vnp;
    }
    
    orig_vnp = orig_next;
    update_vnp = update_next;
  }
  
  for (orig_vnp = *orig_bioseq_list, update_vnp = *update_bioseq_list, seq_num = 0;
       orig_vnp != NULL && update_vnp != NULL;
       orig_vnp = orig_vnp->next, update_vnp = update_vnp->next, seq_num++)
  {
    orig_vnp->choice = seq_num;
    update_vnp->choice = seq_num;
  }
  return no_update_list;
}

static void RemoveSequencesWithoutUpdates (ValNodePtr PNTR orig_bioseq_list, ValNodePtr PNTR update_bioseq_list)
{
  ValNodePtr orig_prev = NULL, update_prev = NULL;
  ValNodePtr orig_next = NULL, update_next = NULL;
  ValNodePtr orig_vnp, update_vnp;
  Int4       seq_num;
  
  if (orig_bioseq_list == NULL || update_bioseq_list == NULL
      || *orig_bioseq_list == NULL || *update_bioseq_list == NULL)
  {
    return;
  }
  
  orig_vnp = *orig_bioseq_list;
  update_vnp = *update_bioseq_list;
  
  while (orig_vnp != NULL && update_vnp != NULL)
  {
    orig_next = orig_vnp->next;
    update_next = update_vnp->next;
    if (orig_vnp->data.ptrvalue == NULL || update_vnp->data.ptrvalue == NULL)
    {
      if (orig_prev == NULL || update_prev == NULL)
      {
        *orig_bioseq_list = orig_vnp->next;
        *update_bioseq_list = update_vnp->next;
      }
      else
      {
        orig_prev->next = orig_vnp->next;
        update_prev->next = update_vnp->next;
      }
      orig_vnp->next = NULL;
      update_vnp->next = NULL;
      ValNodeFree (orig_vnp);
      ValNodeFree (update_vnp);
    }
    else
    {
      orig_prev = orig_vnp;
      update_prev = update_vnp;
    }
    
    orig_vnp = orig_next;
    update_vnp = update_next;
  }
  
  for (orig_vnp = *orig_bioseq_list, update_vnp = *update_bioseq_list, seq_num = 0;
       orig_vnp != NULL && update_vnp != NULL;
       orig_vnp = orig_vnp->next, update_vnp = update_vnp->next, seq_num++)
  {
    orig_vnp->choice = seq_num;
    update_vnp->choice = seq_num;
  }
}

/* This function should find all update Bioseqs that have colliding sequence IDs and
 * replace the colliding IDs with new sequence IDs.
 */
static void ReplaceCollidingUpdateIDs (ValNodePtr update_bioseq_list, ValNodePtr orig_bioseq_list)
{
  ValNodePtr vnp, orig_vnp;
  SeqIdPtr   replace_sip, sip;
  BioseqPtr  bsp;
  Char       id_txt [128];
  Int4       orig_pos;
  
  if (update_bioseq_list == NULL || orig_bioseq_list == NULL)
  {
    return;
  }
  
  for (vnp = update_bioseq_list; vnp != NULL; vnp = vnp->next)
  {
    bsp = (BioseqPtr) vnp->data.ptrvalue;
    if (bsp != NULL)
    {
      if (FindBioseqInList (orig_bioseq_list, bsp->id, &orig_pos))
      {
        orig_vnp = GetNthValNode (orig_bioseq_list, orig_pos);
        if (orig_vnp != NULL && orig_vnp->data.ptrvalue != NULL)
        {
          replace_sip = SeqIdFindBest (bsp->id, SEQID_GENBANK);
          SeqIdWrite (replace_sip, id_txt, PRINTID_REPORT, sizeof (id_txt) - 1);
          StringCat (id_txt, "_update");
          sip = MakeUniqueSeqID (id_txt);
          SeqMgrDeleteFromBioseqIndex (bsp);
          BioseqReplaceID (bsp, sip);
          sip = SeqIdFree (sip);   
        
          SeqMgrReplaceInBioseqIndex (orig_vnp->data.ptrvalue);
          SeqMgrReplaceInBioseqIndex (bsp);
        }
      }
    }
  }
  
}

static SeqEntryPtr ReadASNUpdateSequences (FILE *fp, BoolPtr chars_stripped)
{
  Pointer      dataptr;
  Uint2        datatype;
  SeqEntryPtr  sep;
  SeqSubmitPtr ssp;
  
  /* Read in one sequence from the file */
  dataptr = ReadAsnFastaOrFlatFileEx (fp, &datatype, NULL, FALSE, FALSE,
		                   	                  TRUE, FALSE, chars_stripped);      

  if (NULL == dataptr) 
  {
    return NULL;
  }

  /* Convert the file data to a SeqEntry */
  
  if (datatype == OBJ_SEQENTRY)
    sep = (SeqEntryPtr) dataptr;
  else if (datatype == OBJ_BIOSEQ || datatype == OBJ_BIOSEQSET)
    sep = SeqMgrGetSeqEntryForData (dataptr);
  else if (datatype == OBJ_SEQSUB) 
  {
    ssp = (SeqSubmitPtr) dataptr;
    if (ssp != NULL && ssp->datatype == 1)
    {
      sep = (SeqEntryPtr) ssp->data;
    }
  }
  return sep;  
}

static SeqEntryPtr ReadUpdateSequences (Boolean is_na)
{
  FILE          *fp;
  Char          path [PATH_MAX];
  SeqEntryPtr   sep_list;
  ValNodePtr    err_msg_list = NULL;
  BioseqSetPtr  top_bssp;
  BioseqPtr     bsp;
  Uint2         entityID;
  Boolean       chars_stripped = FALSE;
  
  if (! GetInputFileName (path, sizeof (path),"","TEXT"))
    return NULL;
  fp = FileOpen (path, "r");
  if (fp == NULL)
  {
    Message (MSG_ERROR, "Unable to open %s", path);
    return NULL;
  }
  sep_list = ImportSequencesFromFile (fp, NULL, is_na, TRUE, NULL, &err_msg_list, &chars_stripped);
  ValNodeFreeData (err_msg_list);
  AddUniqueUpdateSequenceIDs (sep_list);
  FileClose (fp);
  
  if (sep_list == NULL)
  {
    fp = FileOpen (path, "r");
    sep_list = ReadASNUpdateSequences (fp, &chars_stripped);
    FileClose (fp);
    
    if (sep_list == NULL)
    {
      return NULL;
    }
    else if (chars_stripped)
    {
      if (ANS_CANCEL == Message (MSG_OKC, "Illegal characters will be stripped from your sequence data.  Do you want to continue?"))
      {
        sep_list = SeqEntryFree (sep_list);
        return NULL;
      }
    }
    
    if (sep_list->choice == 1)
    {
      bsp = (BioseqPtr) sep_list->data.ptrvalue;
      entityID = ObjMgrGetEntityIDForPointer (bsp);
    }
    else
    {
      top_bssp = (BioseqSetPtr) sep_list->data.ptrvalue;
      entityID = ObjMgrGetEntityIDForPointer (top_bssp);
    }
    AddUniqueUpdateSequenceIDs (sep_list);
  }
  else if (sep_list != NULL)
  {
    if (chars_stripped)
    {
      if (ANS_CANCEL == Message (MSG_OKC, "Illegal characters will be stripped from your sequence data.  Do you want to continue?"))
      {
        sep_list = SeqEntryFree (sep_list);
        return NULL;
      }
    }  
    top_bssp = BioseqSetNew ();
    top_bssp->_class = BioseqseqSet_class_genbank;
    top_bssp->seq_set = sep_list;
    sep_list = SeqEntryNew ();
    sep_list->choice = 2;
    sep_list->data.ptrvalue = top_bssp;
    entityID = ObjMgrGetEntityIDForPointer (top_bssp);
  }
  
  AssignIDsInEntityEx (entityID, 0, NULL, NULL);
  SeqMgrIndexFeatures (entityID, NULL);

  return sep_list;
}

typedef struct updatemultisequenceform
{
  FORM_MESSAGE_BLOCK
  DialoG     update_preview;
  DialoG     options_dialog;
  DialoG     titles_dlg;
  ButtoN     update_this;
  ButtoN     skip_this;
  ValNodePtr orig_bioseq_list;
  ValNodePtr update_bioseq_list;
  ValNodePtr unmatched_updates_list;
  ValNodePtr no_updates_list;
  Boolean    is_na;
  LogInfoPtr lip;
  Int4       num_successful;
  Int4       num_failed;
  Int4       num_skipped;
  Char       undo_file [PATH_MAX];
} UpdateMultiSequenceFormData, PNTR UpdateMultiSequenceFormPtr;

static void DoTestUpdateOneSequence (ButtoN b)
{
  UpdateMultiSequenceFormPtr usfp;
  SeqIdPtr              sip;
  Char                  id_txt [MAX_ID_LEN];
  Int4                  orig_pos;
  ValNodePtr            orig_vnp, update_vnp;
  MultiSequenceUpdateData  msud;
  Boolean                  update_successful = FALSE;
  UpdateAlignmentLengthsData uald;
  UpdateOptionsPtr           uop;
  UpdatePairData             upd;
  Uint2                      update_entityID = 0;
  
  usfp = (UpdateMultiSequenceFormPtr) GetObjectExtra (b);
  if (usfp == NULL)
  {
    return;
  }
  
  sip = DialogToPointer (usfp->update_preview);
  if (sip == NULL)
  {
    Message (MSG_ERROR, "No sequence selected!");
    return;
  }
  
  uop = DialogToPointer (usfp->options_dialog);   
  if (uop == NULL || uop->submitter_opts == NULL
      || (uop->submitter_opts->sequence_update_type == eSequenceUpdateNoChange
          && uop->submitter_opts->feature_update_type == eFeatureUpdateNoChange
          && uop->submitter_opts->feature_remove_type == eFeatureRemoveNone))
  {
    Message (MSG_ERROR, "Invalid options selected!");
    uop = UpdateOptionsFree (uop);
    return;
  }

  
  WatchCursor ();
  Update ();
  
  SeqIdWrite (sip, id_txt, PRINTID_REPORT, sizeof (id_txt) - 1);
  
  upd.orig_bsp = FindBioseqInList (usfp->orig_bioseq_list, sip, &orig_pos);
  if (upd.orig_bsp != NULL)
  {
    upd.update_bsp = NULL;
    update_vnp = GetNthValNode (usfp->update_bioseq_list, orig_pos);
    if (update_vnp != NULL)
    {
      upd.update_bsp = update_vnp->data.ptrvalue;
	    update_entityID = upd.update_bsp->idx.entityID;
    }
    
    /* if we are going to ignore the alignment, don't calculate it */
    if (uop->submitter_opts->ignore_alignment)
    {
      upd.revcomp = FALSE;
      upd.salp = NULL;
    }
    else
    {
      upd.revcomp = FALSE;
      upd.salp = Sequin_GlobalAlign2Seq (upd.orig_bsp, upd.update_bsp, &(upd.revcomp)); 
    }
    
    CalculateUpdateAlignmentLengths (upd.salp, upd.orig_bsp, upd.update_bsp, &uald); 
    update_successful = UpdateOrExtendOneSequence (&upd, uop, &uald, 
                                                   usfp->input_entityID,
                                                   usfp->lip == NULL ? NULL : usfp->lip->fp,
                                                   usfp->lip == NULL ? NULL : &(usfp->lip->data_in_log));
	  upd.salp = SeqAlignFree (upd.salp);
  }
  uop = UpdateOptionsFree (uop);
 
  if (update_successful)
  {
    /* remove sequence and its pair from list */
    orig_vnp = ExtractNthValNode (&(usfp->orig_bioseq_list), orig_pos);
    orig_vnp = ValNodeFree (orig_vnp);
    update_vnp = ExtractNthValNode (&(usfp->update_bioseq_list), orig_pos);
    update_vnp = ValNodeFree (update_vnp);

    /* renumber valnode lists */
    for (orig_vnp = usfp->orig_bioseq_list; orig_vnp != NULL; orig_vnp = orig_vnp->next)
    {
      if (orig_vnp->choice > orig_pos)
      {
        orig_vnp->choice --;
      }
    }
    for (update_vnp = usfp->update_bioseq_list; update_vnp != NULL; update_vnp = update_vnp->next)
    {
      if (update_vnp->choice > orig_pos)
      {
        update_vnp->choice --;
      }
    }
    
    /* only update the preview if we have any sequences left to update
     * otherwise we get a flash of a strange result
     */
    if (usfp->orig_bioseq_list != NULL)
    {
      msud.orig_bioseq_list = usfp->orig_bioseq_list;
      msud.update_bioseq_list = usfp->update_bioseq_list;  
      msud.unmatched_updates_list = usfp->unmatched_updates_list;
      msud.no_updates_list = usfp->no_updates_list;  
      PointerToDialog (usfp->update_preview, &msud);    

      /* maintain list position */
      if (orig_pos >= ValNodeLen (msud.orig_bioseq_list))
      {
        orig_pos--;
      }
      SelectUpdatePreviewSequenceByPosition (usfp->update_preview, orig_pos);
    }
    usfp->num_successful ++;                                           
  }
  
  SeqMgrClearFeatureIndexes (usfp->input_entityID, NULL);
  ObjMgrSetDirtyFlag (usfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, usfp->input_entityID, 0, 0);  
 
  /* close window after updating the last sequence */
  if (usfp->orig_bioseq_list == NULL)
  {
    Remove (usfp->form);
  }         
  ArrowCursor ();
  Update ();                                          
}

static void DoNotTestUpdateOneSequence (ButtoN b)
{
  UpdateMultiSequenceFormPtr usfp;
  SeqIdPtr              sip;
  Char                  id_txt [MAX_ID_LEN];
  BioseqPtr                bsp, update_bsp;
  Int4                     orig_pos;
  ValNodePtr               orig_vnp, update_vnp;
  MultiSequenceUpdateData  msud;
  Uint2                    update_entityID;
  
  usfp = (UpdateMultiSequenceFormPtr) GetObjectExtra (b);
  if (usfp == NULL)
  {
    return;
  }
    
  sip = DialogToPointer (usfp->update_preview);
  if (sip == NULL)
  {
    Message (MSG_ERROR, "No sequence selected!");
    return;
  }
  
  SeqIdWrite (SeqIdFindBest (sip, 0), id_txt, PRINTID_REPORT, sizeof (id_txt) - 1);
  if (usfp->lip != NULL && usfp->lip->fp != NULL)
  {
    fprintf (usfp->lip->fp, "Skipped %s\n", id_txt);
    usfp->lip->data_in_log = TRUE;
  }
  
  usfp->num_skipped++;
  bsp = FindBioseqInList (usfp->orig_bioseq_list, sip, &orig_pos);
  if (bsp != NULL)
  {
    /* remove sequence and its pair from list */
    orig_vnp = ExtractNthValNode (&(usfp->orig_bioseq_list), orig_pos);
    orig_vnp = ValNodeFree (orig_vnp);
    update_vnp = ExtractNthValNode (&(usfp->update_bioseq_list), orig_pos);

    /* remove update sequence from Desktop */    
    if (update_vnp != NULL && update_vnp->data.ptrvalue != NULL)
    {
      update_bsp = (BioseqPtr) update_vnp->data.ptrvalue;
      update_bsp->idx.deleteme = TRUE;
      update_entityID = update_bsp->idx.entityID;
      DeleteMarkedObjects (update_entityID, 0, NULL);
    }
    
    update_vnp = ValNodeFree (update_vnp);
    /* renumber valnode lists */
    for (orig_vnp = usfp->orig_bioseq_list; orig_vnp != NULL; orig_vnp = orig_vnp->next)
    {
      if (orig_vnp->choice > orig_pos)
      {
        orig_vnp->choice --;
      }
    }
    for (update_vnp = usfp->update_bioseq_list; update_vnp != NULL; update_vnp = update_vnp->next)
    {
      if (update_vnp->choice > orig_pos)
      {
        update_vnp->choice --;
      }
    }
    msud.orig_bioseq_list = usfp->orig_bioseq_list;
    msud.update_bioseq_list = usfp->update_bioseq_list;   
    msud.unmatched_updates_list = usfp->unmatched_updates_list; 
    msud.no_updates_list = usfp->no_updates_list;                                           
    PointerToDialog (usfp->update_preview, &msud);
    /* maintain list position */
    if (orig_pos >= ValNodeLen (msud.orig_bioseq_list))
    {
      orig_pos--;
    }
    SelectUpdatePreviewSequenceByPosition (usfp->update_preview, orig_pos);
    
  }
  /* close window after skipping the last sequence */
  if (usfp->orig_bioseq_list == NULL)
  {
    Remove (usfp->form);
  }                                                   
}

static void UpdateAllSequences (ButtoN b)
{
  UpdateMultiSequenceFormPtr usfp;
  ValNodePtr                 orig_vnp, update_vnp;
  Char                       id_txt [MAX_ID_LEN];
  UpdatePairData             upd;
  Boolean                    update_successful;
  UpdateOptionsPtr           uop;
  UpdateAlignmentLengthsData uald;

  usfp = (UpdateMultiSequenceFormPtr) GetObjectExtra (b);
  if (usfp == NULL)
  {
    return;
  }

  WatchCursor ();
  Update ();
    
  for (orig_vnp = usfp->orig_bioseq_list, update_vnp = usfp->update_bioseq_list;
       orig_vnp != NULL && update_vnp != NULL;
       orig_vnp = orig_vnp->next, update_vnp = update_vnp->next)
  {
    if (orig_vnp->data.ptrvalue == NULL
        || update_vnp->data.ptrvalue == NULL)
    {
      continue;
    }
    upd.orig_bsp = (BioseqPtr)(orig_vnp->data.ptrvalue);
    upd.update_bsp = (BioseqPtr) (update_vnp->data.ptrvalue);
    /* Get Update Options */
    uop = DialogToPointer (usfp->options_dialog);   
    /* if we are going to ignore the alignment, don't calculate it */
    if (uop->submitter_opts->ignore_alignment)
    {
      upd.revcomp = FALSE;
      upd.salp = NULL;
    }
    else
    {
      upd.revcomp = FALSE;
      upd.salp = Sequin_GlobalAlign2Seq (upd.orig_bsp, upd.update_bsp, &(upd.revcomp)); 
    }
    CalculateUpdateAlignmentLengths (upd.salp, upd.orig_bsp, upd.update_bsp, &uald); 
    update_successful = UpdateOrExtendOneSequence (&upd, uop, &uald, 
                                                   usfp->input_entityID,
                                                   usfp->lip == NULL ? NULL : usfp->lip->fp,
                                                   usfp->lip == NULL ? NULL : &(usfp->lip->data_in_log));
    upd.salp = SeqAlignFree (upd.salp);
    uop = UpdateOptionsFree (uop);
    
    SeqIdWrite (upd.orig_bsp->id, id_txt, PRINTID_REPORT, sizeof (id_txt) - 1);
    if (usfp->lip != NULL && usfp->lip->fp != NULL && ! update_successful)
    {
      fprintf (usfp->lip->fp, "Failed to update %s\n", id_txt);
      usfp->lip->data_in_log = TRUE;
    }  
    if (update_successful)
    {
      usfp->num_successful++;
      /* remove update sequence from list */
      update_vnp->data.ptrvalue = NULL;
    }
    else
    {
      usfp->num_failed++;
    }
    
  }
  
  
  SeqMgrClearFeatureIndexes (usfp->input_entityID, NULL);
  ObjMgrSetDirtyFlag (usfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, usfp->input_entityID, 0, 0);  
    
  /* close window */
  if (usfp->num_failed == 0)
  {
    /* launch log viewer */
    CloseLog (usfp->lip);
    Remove (usfp->form);
  }
  ArrowCursor ();
  Update ();
}

static void UpdateSequenceSelectionChange (Pointer userdata)
{
  UpdateMultiSequenceFormPtr usfp;
  SeqIdPtr                   sip;
  BioseqPtr                  orig_bsp;
  Int4                       orig_pos;
  ValNodePtr                 update_vnp;
  UpdatePairPtr              upp;
  UpdateAlignmentLengthsData uald;
  Boolean                    is_indexer = indexerVersion;
  UpdateOptionsPtr           uop;

  usfp = (UpdateMultiSequenceFormPtr) userdata;
  if (usfp == NULL)
  {
    return;
  }
  
  sip = (SeqIdPtr) DialogToPointer (usfp->update_preview);
  if (sip == NULL)
  {
    Disable (usfp->update_this);
    Disable (usfp->skip_this);
    upp = NULL;
  }
  else
  {
    orig_bsp = FindBioseqInList (usfp->orig_bioseq_list, sip, &orig_pos);
    if (orig_bsp == NULL)
    {
      Disable (usfp->update_this);
      Disable (usfp->skip_this);
    }
    else
    {
      update_vnp = GetNthValNode (usfp->update_bioseq_list, orig_pos);
      if (update_vnp == NULL || update_vnp->data.ptrvalue == NULL)
      {
        Disable (usfp->update_this);
        Enable (usfp->skip_this);
      }
      else
      {
        Enable (usfp->update_this);
        Enable (usfp->skip_this);
      }
    }
    upp = GetCurrentUpdatePair (usfp->update_preview);
  }
  
  if (upp != NULL)
  CalculateUpdateAlignmentLengths (upp->salp,
                                   upp->orig_bsp, 
                                   upp->update_bsp, 
                                   &uald); 
  AdjustUpdateOptionsEnabled (usfp->options_dialog, upp, &uald, is_indexer);
  
  PointerToDialog (usfp->titles_dlg, upp);

  uop = (UpdateOptionsPtr) DialogToPointer (usfp->options_dialog);
  if (uop != NULL && uop->submitter_opts != NULL)
  {
    if (uop->submitter_opts->ignore_alignment)
    {
      if (uop->submitter_opts->sequence_update_type == eSequenceUpdateExtend5)
      {
        SetIgnoreAndExtendForMultiSequencePreview (usfp->update_preview, TRUE, TRUE);
      }
      else if (uop->submitter_opts->sequence_update_type == eSequenceUpdateExtend3)
      {
        SetIgnoreAndExtendForMultiSequencePreview (usfp->update_preview, TRUE, FALSE);
      }
      else
      {
        SetIgnoreAndExtendForMultiSequencePreview (usfp->update_preview, FALSE, FALSE);        
      }
    }
    else
    {
      SetIgnoreAndExtendForMultiSequencePreview (usfp->update_preview, FALSE, FALSE);
    }
  }  
  else
  {
    SetIgnoreAndExtendForMultiSequencePreview (usfp->update_preview, FALSE, FALSE);
  }
  uop = UpdateOptionsFree (uop);
}

static Boolean BioseqHasFeatures (BioseqPtr bsp)
{
  SeqFeatPtr        sfp;
  SeqMgrFeatContext context;
  Uint2             entityID;
  
  if (bsp == NULL)
  {
    return FALSE;
  }

  entityID = bsp->idx.entityID;
  if (entityID == 0)
  {
    entityID = ObjMgrGetEntityIDForPointer (bsp);
  }
  if (! SeqMgrFeaturesAreIndexed (entityID))
  {
    SeqMgrIndexFeatures (entityID, NULL);
  }
  
  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &context);
  if (sfp == NULL)
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}

static Boolean BioseqListHasFeatures (ValNodePtr bioseq_list)
{
  Boolean has_features = FALSE;
  
  while (bioseq_list != NULL && !has_features)
  {
    has_features = BioseqHasFeatures (bioseq_list->data.ptrvalue);
    bioseq_list = bioseq_list->next;
  }
  return has_features;
}



static void CleanupUpdateMultiSequence (GraphiC g, Pointer data)
{
  UpdateMultiSequenceFormPtr usfp;
  ValNodePtr                 vnp;
  BioseqPtr                  update_bsp;
  Uint2                      update_entityID;
  ObjMgrPtr  omp;
  
  usfp = (UpdateMultiSequenceFormPtr) data;
  if (usfp != NULL)
  {
    /* report successes, skips, and failures */
    if (usfp->num_successful > 0)
    {
      if (usfp->num_skipped > 0)
      {
        if (usfp->num_failed > 0)
        {
          Message (MSG_OK, "%d succeeded, %d failed, %d skipped",
                           usfp->num_successful,
                           usfp->num_failed,
                           usfp->num_skipped);
        }
        else
        {
          Message (MSG_OK, "%d succeeded, %d skipped",
                           usfp->num_successful,
                           usfp->num_skipped);
        }
      }
      else if (usfp->num_failed > 0)
      {
        Message (MSG_OK, "%d succeeded, %d failed",
                         usfp->num_successful,
                         usfp->num_failed);
      }
      else
      {
        Message (MSG_OK, "%d succeeded",
                         usfp->num_successful);
      }
    }
    else
    {
      if (usfp->num_skipped > 0)
      {
        if (usfp->num_failed > 0)
        {
          Message (MSG_OK, "%d failed, %d skipped",
                           usfp->num_failed,
                           usfp->num_skipped);
        }
        else
        {
          Message (MSG_OK, "%d skipped",
                           usfp->num_skipped);
        }
      }
      else if (usfp->num_failed > 0)
      {
        Message (MSG_OK, "%d failed",
                         usfp->num_failed);
      }
    }
      
    /* we don't need to free the data for the original bioseq list */
    usfp->orig_bioseq_list = ValNodeFree (usfp->orig_bioseq_list);
    
    /* remove unused update sequences */
    for (vnp = usfp->update_bioseq_list; vnp != NULL; vnp = vnp->next)
    {
      if (vnp->data.ptrvalue != NULL)
      {
        update_bsp = (BioseqPtr) vnp->data.ptrvalue;
        update_bsp->idx.deleteme = TRUE;
        update_entityID = update_bsp->idx.entityID;
        DeleteMarkedObjects (update_entityID, 0, NULL);
      }
    }
    usfp->update_bioseq_list = ValNodeFree (usfp->update_bioseq_list);
    
    /* remove unmatched update sequences */
    for (vnp = usfp->unmatched_updates_list; vnp != NULL; vnp = vnp->next)
    {
      if (vnp->data.ptrvalue != NULL)
      {
        update_bsp = (BioseqPtr) vnp->data.ptrvalue;
        update_bsp->idx.deleteme = TRUE;
        update_entityID = update_bsp->idx.entityID;
        DeleteMarkedObjects (update_entityID, 0, NULL);
      }
    }
    usfp->unmatched_updates_list = ValNodeFree (usfp->unmatched_updates_list);
    
    CloseLog (usfp->lip);
    usfp->lip = FreeLog (usfp->lip);
    
    FileRemove (usfp->undo_file);
    
    /* clear indexes */
    omp = ObjMgrGet ();
    ObjMgrReapOne (omp);
    ObjMgrFreeCache (0);
    FreeSeqIdGiCache ();
  }
  
  StdCleanupExtraProc (g, data);
}

static void LoadUpdateSequenceMapFile (UpdateMultiSequenceFormPtr usfp)
{
  Char                       path [PATH_MAX];
  ReadBufferData             rbd;
  CharPtr                    line, ptr;
  SeqIdPtr                   sip_1, sip_2, sip_orig, sip_update;
  BioseqPtr                  bsp_orig, bsp_update;
  Int4                       orig_position, update_position;
  ValNodePtr                 update_vnp;
  MultiSequenceUpdateData    msud;
  
  if (usfp == NULL)
  {
    return;
  }
  
  if (! GetInputFileName (path, sizeof (path),"","TEXT"))
    return;
  
  rbd.fp = FileOpen (path, "r");
  if (rbd.fp == NULL)
  {
    Message (MSG_ERROR, "Unable to open %s", path);
    return;
  }

  rbd.current_data = NULL;
  line = AbstractReadFunction (&rbd);
  while (line != NULL)
  {
    ptr = line;
    /* skip over initial white space if any*/
    ptr += StringSpn (ptr, " \t");
    /* skip over first sequence ID */ 
    ptr += StringCSpn (ptr, "\t ");
    if (*ptr != 0 && ptr != line)
    {
      /* truncate string after first ID */
      *ptr = '\0';
      ptr++;
      /* skip over any white space before second ID */
      ptr += StringSpn (ptr, "\t ");
      /* truncate trailing white space */
      while (*(ptr + StringLen (ptr) - 1) == ' ' || *(ptr + StringLen (ptr) - 1) == '\t')
      {
        *(ptr + StringLen (ptr) - 1) = 0;
      }
      
      sip_1 = MakeSeqID (line);
      sip_2 = MakeSeqID (ptr);
      
      /* try to determine which is original and which is update */
      bsp_orig = FindBioseqInList (usfp->no_updates_list, sip_1, &orig_position);
      bsp_update = FindBioseqInList (usfp->unmatched_updates_list, sip_2, &update_position);
      if (bsp_orig == NULL && bsp_update == NULL)
      {
        bsp_orig = FindBioseqInList (usfp->no_updates_list, sip_2, &orig_position);
        bsp_update = FindBioseqInList (usfp->unmatched_updates_list, sip_1, &update_position);
        sip_orig = sip_2;
        sip_update = sip_1;
      }
      else
      {
        sip_orig = sip_1;
        sip_update = sip_2;
      }

      if (bsp_orig != NULL && bsp_update != NULL)
      {
        /* remove Bioseq from no update list and add to list of originals */
        update_vnp = ExtractNthValNode (&(usfp->no_updates_list), orig_position);
        ValNodeLink (&(usfp->orig_bioseq_list), update_vnp);
        /* remove Bioseq from unmatched list and add to list of updates */
        update_vnp = ExtractNthValNode (&(usfp->unmatched_updates_list), update_position);
        ValNodeLink (&(usfp->update_bioseq_list), update_vnp);
      }
      sip_orig = SeqIdFree (sip_orig);
      sip_update = SeqIdFree (sip_update);      
    }
    line = MemFree (line);
  	line = AbstractReadFunction (&rbd);
  }
  FileClose (rbd.fp);
  
  msud.orig_bioseq_list = usfp->orig_bioseq_list;
  msud.update_bioseq_list = usfp->update_bioseq_list;     
  msud.unmatched_updates_list = usfp->unmatched_updates_list;
  msud.no_updates_list = usfp->no_updates_list;                                          
  PointerToDialog (usfp->update_preview, &msud);                                                 
}

static void 
ReportIdenticalUpdateSequences 
(ValNodePtr orig_list,
 ValNodePtr update_list)
{
  Char       id_str [128];
  BioseqPtr  orig_bsp, update_bsp;
  LogInfoPtr lip;
  Int4       num_identical = 0, num_total = 0;
  
  lip = OpenLog ("Sequences Identical to Updates");
  if (lip == NULL)
  {
    return;
  }
  WatchCursor ();
  Update ();
  while (orig_list != NULL && update_list != NULL)
  {
    orig_bsp = orig_list->data.ptrvalue;
    update_bsp = update_list->data.ptrvalue;
    
    if (orig_bsp != NULL)
    {
      num_total++;
    }
    
    if (orig_bsp != NULL && update_bsp != NULL 
        && AreSequenceResiduesIdentical (orig_bsp, update_bsp))
    {
      SeqIdWrite (SeqIdFindBest (orig_bsp->id, SEQID_GENBANK), id_str, PRINTID_REPORT, sizeof (id_str) - 1);
      fprintf (lip->fp, "%s\n", id_str);
      lip->data_in_log = TRUE;
      num_identical++;
    }
    
    orig_list = orig_list->next;
    update_list = update_list->next;
  }
  
  if (num_identical == 1)
  {
    Message (MSG_OK, "One sequence out of %d is identical to its update sequence.", num_total);
  } 
  else if (num_identical > 1)
  {
    Message (MSG_OK, "%d sequences out of %d are identical to their update sequences.", num_identical, num_total);
  }
  CloseLog (lip);
  lip = FreeLog (lip);  
  ArrowCursor ();
  Update ();
}

static void RemoveIdenticalUpdates (ButtoN b)
{
  UpdateMultiSequenceFormPtr usfp;
  SeqIdPtr                   sip;
  Int4                       orig_pos = 0, rem_pos = 0;
  ValNodePtr                 orig_list, update_list;
  ValNodePtr                 prev_orig_list = NULL, prev_update_list = NULL;
  ValNodePtr                 next_orig_list = NULL, next_update_list = NULL;
  Uint2                      update_entityID = 0;
  MultiSequenceUpdateData    msud;
  BioseqPtr                  orig_bsp, update_bsp;
    
  usfp = (UpdateMultiSequenceFormPtr) GetObjectExtra (b);
  if (usfp == NULL)
  {
    return;
  }
  
  sip = (SeqIdPtr) DialogToPointer (usfp->update_preview);
  
  FindBioseqInList (usfp->orig_bioseq_list, sip, &orig_pos);

  
  orig_list = usfp->orig_bioseq_list;
  update_list = usfp->update_bioseq_list;

  while (orig_list != NULL && update_list != NULL)
  {
    next_orig_list = orig_list->next;
    next_update_list = update_list->next;
    orig_bsp = orig_list->data.ptrvalue;
    update_bsp = update_list->data.ptrvalue;
    
    if (orig_bsp != NULL && update_bsp != NULL 
        && AreSequenceResiduesIdentical (orig_bsp, update_bsp))
    {
      /* remove update Bioseq from desktop */
      if (update_entityID == 0)
      {
        update_entityID = update_bsp->idx.entityID;
      }
      update_bsp->idx.deleteme = TRUE;
      
      /* remove pair from list */
      if (prev_orig_list == NULL)
      {
        usfp->orig_bioseq_list = orig_list->next;
      }
      else
      {
        prev_orig_list->next = orig_list->next;
      }
      if (prev_update_list == NULL)
      {
        usfp->update_bioseq_list = update_list->next;
      }
      else
      {
        prev_update_list->next = update_list->next;
      }
      
      orig_list->next = NULL;
      orig_list = ValNodeFree (orig_list);
      update_list->next = NULL;
      update_list = ValNodeFree (update_list);
      if (rem_pos < orig_pos)
      {
        orig_pos --;
      }
    }
    else
    {
      prev_orig_list = orig_list;
      prev_update_list = update_list;
      rem_pos++;
    }
    
    orig_list = next_orig_list;
    update_list = next_update_list;
  }
  
  DeleteMarkedObjects (update_entityID, 0, NULL);
  
  ObjMgrSetDirtyFlag (update_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, update_entityID, 0, 0);
  
  msud.orig_bioseq_list = usfp->orig_bioseq_list;
  msud.update_bioseq_list = usfp->update_bioseq_list;   
  msud.unmatched_updates_list = usfp->unmatched_updates_list; 
  msud.no_updates_list = usfp->no_updates_list;                                           
  PointerToDialog (usfp->update_preview, &msud);
  
  /* maintain list position */
  if (orig_pos >= ValNodeLen (msud.orig_bioseq_list))
  {
    orig_pos--;
  }
  SelectUpdatePreviewSequenceByPosition (usfp->update_preview, orig_pos);
    
  /* close window after skipping the last sequence */
  if (usfp->orig_bioseq_list == NULL)
  {
    Message (MSG_OK, "All sequences were identical!");
    Remove (usfp->form);
  }                                                   
  
  
  Update ();    
}

static void LoadUpdateSequenceMapFileBtn (ButtoN b)
{
  UpdateMultiSequenceFormPtr usfp;
  
  usfp = (UpdateMultiSequenceFormPtr) GetObjectExtra (b);
  if (usfp == NULL) {
    return;
  }
  LoadUpdateSequenceMapFile (usfp); 
  ReportIdenticalUpdateSequences (usfp->orig_bioseq_list, usfp->update_bioseq_list);
  if (BioseqListHasFeatures (usfp->update_bioseq_list)) {
    EnableImportFeatureOptions(usfp->options_dialog);
  } else {
    DisableImportFeatureOptions (usfp->options_dialog);        
  }  
  if (BioseqListHasFeatures (usfp->orig_bioseq_list)) {
    EnableRemoveFeatureOptions (usfp->options_dialog);
  } else {
    DisableRemoveFeatureOptions (usfp->options_dialog);
  }
}

extern SeqEntryPtr RestoreFromFile (CharPtr path);

static void UndoUpdates (ButtoN b)
{
  UpdateMultiSequenceFormPtr usfp;
  SeqEntryPtr                oldsep, currsep;
  
  usfp = (UpdateMultiSequenceFormPtr) GetObjectExtra (b);
  if (usfp == NULL)
  {
    return;
  }
  WatchCursor ();
  Update ();

  SeqEntrySetScope (NULL);
  oldsep = RestoreFromFile (usfp->undo_file);
  currsep = GetTopSeqEntryForEntityID (usfp->input_entityID);
  ReplaceSeqEntryWithSeqEntry (currsep, oldsep, TRUE);
  SeqEntrySetScope (NULL);
  usfp->input_entityID = ObjMgrGetEntityIDForChoice (currsep);
  ObjMgrSetDirtyFlag (usfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, usfp->input_entityID, 0, 0);

  usfp->num_successful = 0;
  usfp->num_failed = 0;
  usfp->num_skipped = 0;
  usfp->lip = FreeLog (usfp->lip);

  Remove (usfp->form);
  ArrowCursor ();
  Update ();  
}

static void TestUpdateSequenceSet (IteM i, Boolean is_indexer, Boolean do_update)
{
  BaseFormPtr                bfp;
  WindoW                     w;
  CharPtr                    title;
  GrouP                      h, k, ext_grp;
  BioseqPtr                  orig_bsp;
  SeqEntryPtr                sep;
  SeqEntryPtr                update_list;
  GrouP                      c;
  ButtoN                     b;
  Boolean                    is_na;
  Int4                       num_orig = 0, num_update = 0;
  UpdateMultiSequenceFormPtr usfp;
  MultiSequenceUpdateData    msud;
  UpdateOptionsPtr           uop;
  DialoG                     ext_dlg;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  /* for test purposes, need to load sequence and update sequence */
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL)
    return;
  
  orig_bsp = GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID,
			   bfp->input_itemtype);
  if (orig_bsp == NULL)
  {
    is_na = TRUE;
  }
  else
  {
    is_na = ISA_na (orig_bsp->mol);
  }

  /* Read in the update data from a file */
  /* for now, just handling FASTA */
  update_list = ReadUpdateSequences (is_na);
  if (update_list == NULL)
  {
    return;
  }

  usfp = (UpdateMultiSequenceFormPtr) MemNew (sizeof (UpdateMultiSequenceFormData));
  if (usfp == NULL)
  {
    return;
  }
  usfp->input_entityID = bfp->input_entityID;
  usfp->is_na = is_na;
  usfp->orig_bioseq_list = NULL;
  usfp->update_bioseq_list = NULL;
  
  usfp->num_successful = 0;
  usfp->num_failed = 0;
  usfp->num_skipped = 0;
  
  usfp->lip = OpenLog ("Update Sequence Log");
  
  /* save the current SeqEntry to a file, so we can restore it if we need to */
  TmpNam (usfp->undo_file);
  write_out_put(sep, usfp->undo_file, FALSE);
    
  ListBioseqsInSeqEntry (sep, usfp->is_na, &num_orig, &usfp->orig_bioseq_list);
  ListBioseqsInSeqEntry (update_list, usfp->is_na, &num_update, &usfp->update_bioseq_list);

  usfp->unmatched_updates_list = ShuffleUpdateBioseqList (&usfp->update_bioseq_list, 
                                                          usfp->orig_bioseq_list);

  usfp->no_updates_list = ExtractSequencesWithoutUpdates (&usfp->orig_bioseq_list, &usfp->update_bioseq_list);  
  if (!is_indexer)
  {
    if (usfp->unmatched_updates_list != NULL)
    {
      if (ANS_YES == Message (MSG_YN, "Some sequences in your update file are not present in the record to be updated.  Would you like to load a tab-delimited file to map the update files to the record files?"))
      {
        LoadUpdateSequenceMapFile (usfp);
      }
    }
    RemoveSequencesWithoutUpdates (&usfp->orig_bioseq_list, &usfp->update_bioseq_list);
    usfp->no_updates_list = ValNodeFree (usfp->no_updates_list);
  }
  
  ReplaceCollidingUpdateIDs (usfp->update_bioseq_list, usfp->orig_bioseq_list);
  
  /* Create window */

  if (do_update) {
    title = "Update Sequence";
  } else {
    title = "Extend Sequence";
  }
  w = MovableModalWindow (-50, -33, -10, -10, title, NULL);

  if (w == NULL)
    return;
  
  SetObjectExtra (w, usfp, CleanupUpdateMultiSequence);
  usfp->form = (ForM) w;
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  usfp->titles_dlg = UpdateTitlesDialog (h, is_indexer, do_update);
  
  k = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (k, 10, 10);
  
  usfp->update_preview = MultiSequenceUpdatePreview (k, usfp->is_na, do_update,
                                                     UpdateSequenceSelectionChange,
                                                     usfp,
                                                     LoadUpdateSequenceMapFileBtn,
                                                     is_indexer);

  if (! do_update)
  {
    ext_grp = HiddenGroup (k, -1, 0, NULL);
    
    ext_dlg = UpdateTitlesDialog (ext_grp, is_indexer, FALSE);
    AddExtendDialogToMultiSequenceUpdatePreview (usfp->update_preview, ext_dlg);
                                                 ;
    
    usfp->options_dialog = UpdateOptionsDialog (ext_grp, usfp->is_na, is_indexer, do_update,
                                                UpdateSequenceSelectionChange, usfp);
    AlignObjects (ALIGN_CENTER, (HANDLE) ext_dlg, (HANDLE) usfp->options_dialog, NULL);
    SetMultiSequenceUpdatePreviewMatchWindow (usfp->update_preview, (WindoW) ext_grp);
  }
  else
  {
    usfp->options_dialog = UpdateOptionsDialog (k, usfp->is_na, is_indexer, do_update,
                                                UpdateSequenceSelectionChange, usfp);
  }


  c = HiddenGroup (h, 6, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  usfp->update_this = PushButton (c, "Update This Sequence", DoTestUpdateOneSequence);
  SetObjectExtra (usfp->update_this, usfp, NULL);
  usfp->skip_this = PushButton (c, "Skip This Sequence", DoNotTestUpdateOneSequence);
  SetObjectExtra (usfp->skip_this, usfp, NULL);
  b = PushButton (c, "Update All Sequences", UpdateAllSequences);
  SetObjectExtra (b, usfp, NULL);
  
  b = PushButton (c, "Stop Updating", StdCancelButtonProc);
  
  b = PushButton (c, "Cancel", UndoUpdates);
  SetObjectExtra (b, usfp, NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) usfp->titles_dlg, (HANDLE) k, (HANDLE) c, NULL);
 
  msud.orig_bioseq_list = usfp->orig_bioseq_list;
  msud.update_bioseq_list = usfp->update_bioseq_list;     
  msud.unmatched_updates_list = usfp->unmatched_updates_list; 
  msud.no_updates_list = usfp->no_updates_list;                                         
  PointerToDialog (usfp->update_preview, &msud);                                               

  if (!BioseqListHasFeatures (usfp->update_bioseq_list))
  {
    DisableImportFeatureOptions (usfp->options_dialog);        
  }

  if (BioseqListHasFeatures (usfp->orig_bioseq_list)) {
    EnableRemoveFeatureOptions (usfp->options_dialog);
  } else {
    DisableRemoveFeatureOptions (usfp->options_dialog);
  }

  uop = DialogToPointer (usfp->options_dialog);
  uop->submitter_opts->sequence_update_type = eSequenceUpdateReplace;
  PointerToDialog (usfp->options_dialog, uop);
  uop = UpdateOptionsFree (uop);

  Show (w);
  Update ();
  ReportIdenticalUpdateSequences (usfp->orig_bioseq_list, usfp->update_bioseq_list);
}

extern void TestUpdateSequenceSetSubmitter (IteM i)
{
  TestUpdateSequenceSet (i, FALSE, TRUE);
}

extern void TestExtendSequenceSetSubmitter (IteM i)
{
  TestUpdateSequenceSet (i, FALSE, FALSE);
}

extern void TestUpdateSequenceSetIndexer (IteM i)
{
  TestUpdateSequenceSet (i, TRUE, TRUE);
}

extern void TestExtendSequenceSetIndexer (IteM i)
{
  TestUpdateSequenceSet (i, TRUE, FALSE);
}

typedef struct singlesequenceupdateform
{
  FORM_MESSAGE_BLOCK
  DialoG     preview_dlg;
  DialoG     options_dlg;
  DialoG     titles_dlg;
  Boolean    is_na;
  LogInfoPtr lip;
  UpdatePairData update_pair;
  SeqEntryPtr    update_list;
  Uint2          update_entity_ID;
} SingleSequenceUpdateFormData, PNTR SingleSequenceUpdateFormPtr;

static void MarkUpdateListForDeletion (SeqEntryPtr sep)
{
  BioseqPtr    bsp;
  BioseqSetPtr bssp;
  
  if (sep == NULL)
  {
    return;
  }
  MarkUpdateListForDeletion (sep->next);
  if (IS_Bioseq (sep))
  {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp != NULL)
    {
      bsp->idx.deleteme = TRUE;
    }
  }
  else if (IS_Bioseq_set (sep))
  {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL)
    {
      bssp->idx.deleteme = TRUE;
      MarkUpdateListForDeletion (bssp->seq_set);
    }
  }
}

static void CleanupSingleSequenceUpdateForm (GraphiC g, Pointer data)
{
  SingleSequenceUpdateFormPtr ssufp;
  
  ssufp = (SingleSequenceUpdateFormPtr) data;
  if (ssufp != NULL)
  {
    CloseLog (ssufp->lip);
    ssufp->lip = FreeLog (ssufp->lip);
    
    if (ssufp->update_pair.update_bsp != NULL)
    {
      ssufp->update_pair.update_bsp->idx.deleteme = TRUE;
    }
    MarkUpdateListForDeletion (ssufp->update_list);
    DeleteMarkedObjects (ssufp->update_entity_ID, 0, NULL);
  }
  StdCleanupExtraProc (g, data);
}

static void DoUpdateSingleSequence (ButtoN b)
{
  SingleSequenceUpdateFormPtr ssufp;
  Char                        id_txt [MAX_ID_LEN];
  UpdateAlignmentLengthsData  uald;
  UpdateOptionsPtr            uop;
  Boolean                     update_successful = FALSE;

  ssufp = (SingleSequenceUpdateFormPtr) GetObjectExtra (b);

  WatchCursor ();
  Update ();
  
  if (ssufp == NULL || ssufp->update_pair.orig_bsp == NULL)
  {
    return;
  }

  uop = DialogToPointer (ssufp->options_dlg);   
  if (uop == NULL || uop->submitter_opts == NULL
      || (uop->submitter_opts->sequence_update_type == eSequenceUpdateNoChange
          && uop->submitter_opts->feature_update_type == eFeatureUpdateNoChange
          && uop->submitter_opts->feature_remove_type == eFeatureRemoveNone))
  {
    Message (MSG_ERROR, "Invalid options selected!");
    uop = UpdateOptionsFree (uop);
    return;
  }
            
  /* if we are going to ignore the alignment, remove it */
  if (uop->submitter_opts->ignore_alignment)
  {
    ssufp->update_pair.salp = NULL;
    ssufp->update_pair.revcomp = FALSE;
  }
  else
  {
    /* reverse the sequence prior to the update if necessary */
    if (ssufp->update_pair.revcomp)
    {
      BioseqRevComp (ssufp->update_pair.update_bsp);
      ReverseBioseqFeatureStrands (ssufp->update_pair.update_bsp);
      SeqMgrReplaceInBioseqIndex (ssufp->update_pair.update_bsp);
    }
  }
 
  CalculateUpdateAlignmentLengths (ssufp->update_pair.salp,
                                   ssufp->update_pair.orig_bsp, 
                                   ssufp->update_pair.update_bsp, 
                                   &uald); 

  update_successful = UpdateOrExtendOneSequence (&(ssufp->update_pair), 
                                                 uop, 
                                                 &uald,
                                                 ssufp->input_entityID, 
                                                 ssufp->lip->fp, &(ssufp->lip->data_in_log));
                                                 
  if (!update_successful && ssufp->lip != NULL && ssufp->lip->fp != NULL)
  {
    SeqIdWrite (ssufp->update_pair.orig_bsp->id, id_txt, PRINTID_REPORT, sizeof (id_txt) - 1);
    fprintf (ssufp->lip->fp, "Failed to update %s\n", id_txt);
    ssufp->lip->data_in_log = TRUE;
  }
  uop = UpdateOptionsFree (uop);
  
    
  SeqMgrClearFeatureIndexes (ssufp->input_entityID, NULL);
  ObjMgrSetDirtyFlag (ssufp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ssufp->input_entityID, 0, 0);  

  Remove (ssufp->form);
  ArrowCursor ();
  Update ();
}

static void UpdateSingleSequenceOptions (Pointer userdata)
{
  SingleSequenceUpdateFormPtr ssufp;
  UpdateAlignmentLengthsData  uald;
  Boolean                     is_indexer = indexerVersion;
  UpdateOptionsPtr            uop;

  ssufp = (SingleSequenceUpdateFormPtr) userdata;
  if (ssufp == NULL)
  {
    return;
  }
  
  CalculateUpdateAlignmentLengths (ssufp->update_pair.salp,
                                   ssufp->update_pair.orig_bsp, 
                                   ssufp->update_pair.update_bsp, 
                                   &uald); 
  
  AdjustUpdateOptionsEnabled (ssufp->options_dlg, &(ssufp->update_pair), &uald, is_indexer);
  uop = (UpdateOptionsPtr) DialogToPointer (ssufp->options_dlg);
  if (uop != NULL && uop->submitter_opts != NULL)
  {
    if (uop->submitter_opts->ignore_alignment)
    {
      if (uop->submitter_opts->sequence_update_type == eSequenceUpdateExtend5)
      {
        SetIgnoreAndExtendForPreviewPictures (ssufp->preview_dlg, TRUE, TRUE);
      }
      else if (uop->submitter_opts->sequence_update_type == eSequenceUpdateExtend3)
      {
        SetIgnoreAndExtendForPreviewPictures (ssufp->preview_dlg, TRUE, FALSE);
      }
      else
      {
        SetIgnoreAndExtendForPreviewPictures (ssufp->preview_dlg, FALSE, FALSE);        
      }
    }
    else
    {
      SetIgnoreAndExtendForPreviewPictures (ssufp->preview_dlg, FALSE, FALSE);
    }
  }  
  else
  {
    SetIgnoreAndExtendForPreviewPictures (ssufp->preview_dlg, FALSE, FALSE);
  }
  uop = UpdateOptionsFree (uop);  
}

static void 
UpdateSingleSequence 
(BioseqPtr   orig_bsp,
 SeqEntryPtr update_list,
 Boolean     is_indexer,
 Boolean     do_update,
 Uint2       entityID)
{
  WindoW             w;
  CharPtr            title;
  GrouP              h, k;
  GrouP              c;
  ButtoN             b;
  Boolean            is_na;
  Int4               num_update = 0;
  SingleSequenceUpdateFormPtr ssufp;
  ValNodePtr                  orig_bioseq_list = NULL, update_bioseq_list = NULL;
  UpdateAlignmentLengthsData  uald;
  UpdateOptionsPtr            uop;

  if (orig_bsp == NULL || update_list == NULL)
  {
    return;
  }

  
  is_na = ISA_na (orig_bsp->mol);

  ssufp = (SingleSequenceUpdateFormPtr) MemNew (sizeof (SingleSequenceUpdateFormData));
  if (ssufp == NULL)
  {
    return;
  }
  ssufp->is_na = is_na;
  ssufp->lip = OpenLog ("Update Sequence Log");
  
  ssufp->update_pair.orig_bsp = orig_bsp;
  
  ssufp->update_list = update_list;
  
  ValNodeAddPointer (&orig_bioseq_list, 0, orig_bsp);
  ListBioseqsInSeqEntry (update_list, ssufp->is_na, &num_update, &update_bioseq_list);

  ReplaceCollidingUpdateIDs (update_bioseq_list, orig_bioseq_list);
  
  ssufp->update_pair.update_bsp = update_bioseq_list->data.ptrvalue;
  update_bioseq_list = ValNodeFree (update_bioseq_list);

  if (ssufp->update_pair.update_bsp != NULL)
  {
    ssufp->update_entity_ID = ssufp->update_pair.update_bsp->idx.entityID;    
  }
  
  if (AreSequenceResiduesIdentical (ssufp->update_pair.orig_bsp,
                                    ssufp->update_pair.update_bsp))
  {
    Message (MSG_OK, "Sequence is identical to update sequence!");
  }
  
  /* Create window */

  if (do_update) {
    title = "Update Sequence";
  } else {
    title = "Extend Sequence";
  }
  w = FixedWindow (-50, -33, -10, -10, title, NULL);

  if (w == NULL)
    return;
  
  SetObjectExtra (w, ssufp, CleanupSingleSequenceUpdateForm);
  ssufp->form = (ForM) w;
  ssufp->input_entityID = entityID;
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  ssufp->titles_dlg = UpdateTitlesDialog (h, is_indexer, do_update);
  
  k = HiddenGroup (h, 2, 0, NULL);
  
  SetGroupSpacing (k, 10, 10);
  if (do_update)
  {
    ssufp->preview_dlg = UpdatePreviewDialog (k, is_indexer);
  }
  
  ssufp->update_pair.revcomp = FALSE;
  ssufp->update_pair.salp = Sequin_GlobalAlign2Seq (ssufp->update_pair.orig_bsp,
                                                    ssufp->update_pair.update_bsp, 
                                                    &(ssufp->update_pair.revcomp));
  if (ssufp->update_pair.revcomp)
  {
    BioseqRevComp (ssufp->update_pair.update_bsp);
    ReverseBioseqFeatureStrands (ssufp->update_pair.update_bsp);
    SeqMgrReplaceInBioseqIndex (ssufp->update_pair.update_bsp);
  }
  
  ssufp->options_dlg = UpdateOptionsDialog (k, ssufp->is_na, is_indexer, do_update,
                                            UpdateSingleSequenceOptions, ssufp);
  PointerToDialog (ssufp->preview_dlg, &(ssufp->update_pair));  
  PointerToDialog (ssufp->titles_dlg, &(ssufp->update_pair));
  
  c = HiddenGroup (h, 4, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  b = PushButton (c, "Update Sequence", DoUpdateSingleSequence);
  SetObjectExtra (b, ssufp, NULL);
  b = PushButton (c, "Cancel", StdCancelButtonProc);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) ssufp->titles_dlg, (HANDLE) k, (HANDLE) c, NULL);
  
  if (BioseqHasFeatures (ssufp->update_pair.update_bsp)) {
    EnableImportFeatureOptions(ssufp->options_dlg);
    uop = DialogToPointer (ssufp->options_dlg);
    /* for public version only, if there are features, default to
     * import all except duplicates
     */
    if (!is_indexer) {
      uop->submitter_opts->feature_update_type = eFeatureUpdateAllExceptDups;
      PointerToDialog (ssufp->options_dlg, uop);
      uop = UpdateOptionsFree (uop);
    }
  } else {
    DisableImportFeatureOptions (ssufp->options_dlg);            
  }

  if (BioseqHasFeatures (ssufp->update_pair.orig_bsp)) {
    EnableRemoveFeatureOptions (ssufp->options_dlg);
  } else {
    DisableRemoveFeatureOptions (ssufp->options_dlg);
  }
  
  CalculateUpdateAlignmentLengths (ssufp->update_pair.salp,
                                   ssufp->update_pair.orig_bsp, 
                                   ssufp->update_pair.update_bsp, 
                                   &uald); 
  AdjustUpdateOptionsEnabled (ssufp->options_dlg, &(ssufp->update_pair), &uald, is_indexer);

  uop = DialogToPointer (ssufp->options_dlg);
  if (is_indexer)
  {
    uop->submitter_opts->sequence_update_type = eSequenceUpdateReplace;
  }
  else
  {
    if (ssufp->update_pair.salp == NULL) 
    {
      /* no alignment - default is to replace */
      uop->submitter_opts->sequence_update_type = eSequenceUpdateReplace;
    }
    else if (uald.new5 > uald.old5 && uald.new3 < uald.old3) 
    {
      /* new 5' end is bigger, new 3' end is smaller,
       * default is extend 5' */
      uop->submitter_opts->sequence_update_type = eSequenceUpdateExtend5;
    }
    else if (uald.new5 < uald.old5 && uald.new3 > uald.old3)
    {
      /* new 3' end is bigger, new 5' end is smaller,
       * default is extend 3' */
      uop->submitter_opts->sequence_update_type = eSequenceUpdateExtend3;
    }
    else
    {
      uop->submitter_opts->sequence_update_type = eSequenceUpdateReplace;
    }
  }
  PointerToDialog (ssufp->options_dlg, uop);
  uop = UpdateOptionsFree (uop);
  
  Show (w);
  Update ();  
}

static void TestUpdateSequence (IteM i, Boolean is_indexer, Boolean do_update)
{
  BaseFormPtr        bfp;
  BioseqPtr          orig_bsp;
  SeqEntryPtr        sep;
  SeqEntryPtr        update_list;
  Boolean            is_na;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  /* for test purposes, need to load sequence and update sequence */
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL)
    return;
  orig_bsp = GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID,
			   bfp->input_itemtype);
  if (orig_bsp == NULL)
    return;
  
  is_na = ISA_na (orig_bsp->mol);

  /* Read in the update data from a file */
  /* for now, just handling FASTA */
  update_list = ReadUpdateSequences (is_na);
  if (update_list == NULL)
  {
    return;
  }
  
  UpdateSingleSequence (orig_bsp, update_list, is_indexer, do_update, bfp->input_entityID);
}

extern void TestUpdateSequenceIndexer (IteM i)
{
  TestUpdateSequence (i, TRUE, TRUE);
}

extern void TestExtendSequenceIndexer (IteM i)
{
  TestUpdateSequence (i, TRUE, FALSE);
}

extern void TestUpdateSequenceSubmitter (IteM i)
{
  TestUpdateSequence (i, FALSE, TRUE);
}

extern void TestExtendSequenceSubmitter (IteM i)
{
  TestUpdateSequence (i, FALSE, FALSE);
}

typedef struct seqentrydownload
{
  FORM_MESSAGE_BLOCK
  GrouP           accntype;
  TexT            accession;
  ButtoN          accept;
} SeqEntryDownloadData, PNTR SeqEntryDownloadPtr;

static void TypeAccessionProc (TexT t)

{
  Boolean       alldigits;
  Char          ch;
  SeqEntryDownloadPtr sedp;
  CharPtr       ptr;
  Char          str [32];

  sedp = (SeqEntryDownloadPtr) GetObjectExtra (t);
  if (sedp == NULL) return;
  GetTitle (t, str, sizeof (str));
  if (StringHasNoText (str)) {
    SafeDisable (sedp->accept);
  } else {
    SafeEnable (sedp->accept);
    TrimSpacesAroundString (str);
    alldigits = TRUE;
    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      if (! IS_DIGIT (ch)) {
        alldigits = FALSE;
      }
      ptr++;
      ch = *ptr;
    }
    if (alldigits) {
      SafeSetValue (sedp->accntype, 2);
    } else {
      SafeSetValue (sedp->accntype, 1);
    }
  }
}

static SeqEntryPtr DownloadUpdateSequence (void)
{
  GrouP                 c;
  SeqEntryDownloadData  sedd;
  GrouP                 g;
  WindoW                w;
  ModalAcceptCancelData acd;
  ButtoN                b;
  SeqEntryPtr           fetched_sep = NULL;
  Char                  str [32];
  Int4                  uid;  
  SeqIdPtr              sip;
  Uint2                 entityID;
  BioseqPtr             bsp;
  BioseqSetPtr          bssp;
  
  w = MovableModalWindow (-50, -33, -10, -10, "Download From Entrez", NULL);
  SetObjectExtra (w, &sedd, NULL);
  sedd.form = (ForM) w;
  SetGroupSpacing (w, 10, 10);

  g = HiddenGroup (w, -3, 0, NULL);
  StaticPrompt (g, "Type", 0, stdLineHeight, programFont, 'l');
  sedd.accntype = HiddenGroup (g, 4, 0, NULL);
  RadioButton (sedd.accntype, "Accession");
  RadioButton (sedd.accntype, "GI");
  SetValue (sedd.accntype, 1);
  sedd.accession = DialogText (g, "", 6, TypeAccessionProc);
  SetObjectExtra (sedd.accession, &sedd, NULL);

  c = HiddenGroup (w, 4, 0, NULL);
  SetGroupSpacing (c, 10, 2);
  sedd.accept = DefaultButton (c, "Retrieve", ModalAcceptButton);
  SetObjectExtra (sedd.accept, &acd, NULL);
  Disable (sedd.accept);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
  RealizeWindow (w);

  Select (sedd.accession);
  Show (w);
  Select (w);
  Update ();

  acd.accepted = FALSE;
  acd.cancelled = FALSE;
  while (!acd.accepted && ! acd.cancelled)
  {
    ProcessExternalEvent ();
    Update ();
  }
  ProcessAnEvent ();

  if (acd.accepted)
  {
    Hide (w);
    WatchCursor ();
    Update ();
    uid = 0;
    
    GetTitle (sedd.accession, str, sizeof (str));
    if (!StringHasNoText (str)) 
    {
      uid = 0;
      if (GetValue (sedd.accntype) == 1) 
      {
        sip = SeqIdFromAccessionDotVersion (str);
        if (sip != NULL)
        {
          uid = GetGIForSeqId (sip);
          sip = SeqIdFree (sip);
        }
      }
      else 
      {
        if (! StrToLong (str, &uid)) 
        {
          uid = 0;
        }
      }
      ArrowCursor ();
      Update ();
      if (uid > 0) 
      {
        fetched_sep = PubSeqSynchronousQuery (uid, 0, -1);
      }
      if (fetched_sep == NULL) 
      {
        Message (MSG_OK, "Unable to find this record in the database.");
      }
      else
      {
        if (IS_Bioseq (fetched_sep))
        {
          bsp = (BioseqPtr) fetched_sep->data.ptrvalue;
          entityID = ObjMgrGetEntityIDForPointer (bsp);
        }
        else
        {
          bssp = (BioseqSetPtr) fetched_sep->data.ptrvalue;
          entityID = ObjMgrGetEntityIDForPointer (bssp);
        }
        AssignIDsInEntityEx (entityID, 0, NULL, NULL);
      }
    }
  }
  Remove (w);
  return fetched_sep;
}


static void UpdateSequenceViaDownload (IteM i, Boolean is_indexer)
{
  BaseFormPtr bfp;
  SeqEntryPtr sep;
  BioseqPtr   orig_bsp;
  Boolean     is_na;
  SeqEntryPtr update_sep;
  
#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  /* for test purposes, need to load sequence and update sequence */
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL)
    return;
  orig_bsp = GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID,
			   bfp->input_itemtype);
  if (orig_bsp == NULL)
    return;
  
  is_na = ISA_na (orig_bsp->mol);

  update_sep = DownloadUpdateSequence ();
  if (update_sep == NULL)
  {
    return;
  }
 
  UpdateSingleSequence (orig_bsp, update_sep, is_indexer, TRUE, bfp->input_entityID);
}

extern void UpdateSequenceViaDownloadIndexer (IteM i)
{
  UpdateSequenceViaDownload (i, TRUE);
}

extern void UpdateSequenceViaDownloadSubmitter (IteM i)
{
  UpdateSequenceViaDownload (i, FALSE);
}

