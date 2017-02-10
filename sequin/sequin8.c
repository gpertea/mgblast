/*   sequin8.c
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
* File Name:  sequin8.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   2/3/98
*
* $Revision: 6.356 $
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

#include "sequin.h"
#include <objsub.h>
#include <valid.h>
#include <cdrgn.h>
#include <suggslp.h>
#include <toasn3.h>
#include <subutil.h>
#include <explore.h>
#include <medarch.h>
#include <medutil.h>
#include <tofasta.h>
#include <asn2gnbk.h>
#include <alignmgr2.h>
#include <spidey.h>
#include <blast.h>
#include <salpanel.h>
#include <edutil.h>
#include <asn2gnbp.h> /* included for discrepancy report */


#define DEFLINE_MAX_LEN          380
#define TEXT_MAX_LEN             64
#define DEFLINE_MAX_GENENAME_LEN 64
#define ALL_FEATURES             255

typedef struct evidenceformdata {
  FEATURE_FORM_BLOCK

  LisT           objlist;
  TexT           findthis;
  Uint2          itemtype;
  Uint2          subtype;
  PopuP          evdence;
  Uint2          exp_ev;
  ValNodePtr     head;
  Boolean        stringfound;
  Char           findStr [128];
  ButtoN         case_insensitive;
  ButtoN         when_string_not_present;
} EvidenceFormData, PNTR 

EvidenceFormPtr;

typedef struct codebreakformdata {
  FEATURE_FORM_BLOCK
  PopuP  aminoAcidPopup;
  Char   currentCodonStr [4];
  TexT   codonText;
  ButtoN acceptButton;
} CodeBreakFormData, PNTR CodeBreakFormPtr;

static Boolean IsRealImpFeat (Uint2 subtype)

{
  if (subtype >= FEATDEF_allele && subtype <= FEATDEF_site_ref) return TRUE;
  if (subtype == FEATDEF_oriT) return TRUE;
  return FALSE;
}



static void BreakIntoAGroup (BioseqSetPtr parent, Uint1 _class, SeqEntryPtr list)

{
  BioseqSetPtr  bssp;
  Int2          count;
  SeqEntryPtr   sep;
  SeqEntryPtr   tmp;

  while (list != NULL) {
    bssp = BioseqSetNew ();
    if (bssp == NULL) return;
    bssp->_class = _class;
    sep = SeqEntryNew ();
    if (sep == NULL) return;
    sep->choice = 2;
    sep->data.ptrvalue = (Pointer) bssp;
    if (parent->seq_set != NULL) {
      tmp = parent->seq_set;
      while (tmp->next != NULL) {
        tmp = tmp->next;
      }
      tmp->next = sep;
    } else {
      parent->seq_set = sep;
    }
    bssp->seq_set = list;
    for (tmp = list, count = 0; tmp != NULL && count < 99; tmp = tmp->next, count++) continue;
    if (tmp != NULL) {
      list = tmp->next;
      tmp->next = NULL;
    } else {
      list = NULL;
    }
  }
}

extern void MakeGroupsOf200 (IteM i)

{
  AsnIoPtr       aip;
  BaseFormPtr    bfp;
  BioseqSetPtr   bssp;
  Uint1          _class;
  Int2           count;
  Char           file [FILENAME_MAX];
  SeqEntryPtr    list;
  SeqEntryPtr    next;
  ObjMgrDataPtr  omdptop;
  ObjMgrData     omdata;
  Char           output [PATH_MAX];
  Uint2          parenttype;
  Pointer        parentptr;
  Char           path [PATH_MAX];
  CharPtr        ptr;
  SeqEntryPtr    sep;
  SeqSubmitPtr   ssp;
  Char           str [FILENAME_MAX];
  SeqEntryPtr    tmp;
#ifdef WIN_MAC
  FILE           *f;
#endif

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  if (! IS_Bioseq_set (sep)) return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp == NULL) return;
  _class = bssp->_class;
  if (_class != 7 && _class != 13 && _class != 14 &&
      _class != 15 && _class != 16 && _class != 18) return;

  SaveSeqEntryObjMgrData (sep, &omdptop, &omdata);
  GetSeqEntryParent (sep, &parentptr, &parenttype);

  list = bssp->seq_set;
  bssp->seq_set = NULL;
  bssp->_class = 7;
  BreakIntoAGroup (bssp, _class, list);

  SeqMgrLinkSeqEntry (sep, parenttype, parentptr);
  RestoreSeqEntryObjMgrData (sep, omdptop, &omdata);
  PropagateFromGenBankBioseqSet (sep, TRUE);

  if (parenttype == OBJ_SEQSUB) {
    if (GetOutputFileName (path, sizeof (path), "")) {
      ssp = (SeqSubmitPtr) parentptr;
      if (ssp != NULL && ssp->datatype == 1) {
        sep = (SeqEntryPtr) ssp->data;
        ptr = StringRChr (path, DIRDELIMCHR);
        if (ptr != NULL) {
          ptr++;
          StringNCpy_0 (file, ptr, sizeof (file));
          *ptr = '\0';
          tmp = bssp->seq_set;
          count = 0;
          while (tmp != NULL) {
            next = tmp->next;
            tmp->next = NULL;
            ssp->data = (Pointer) tmp;
            StringCpy (output, path);
            count++;
            if (count < 10) {
              sprintf (str, "%s0%1d", file, (int) count);
            } else {
              sprintf (str, "%s%2d", file, (int) count);
            }
            FileBuildPath (output, NULL, str);
#ifdef WIN_MAC
            f = FileOpen (output, "r");
            if (f != NULL) {
              FileClose (f);
            } else {
              FileCreate (output, "TEXT", "ttxt");
            }
#endif
            aip = AsnIoOpen (output, "w");
            if (aip != NULL) {
              SeqSubmitAsnWrite (ssp, aip, NULL);
            }
            AsnIoClose (aip);
            tmp->next = next;
            tmp = next;
          }
          ssp->data = (Pointer) sep;
        }
      }
    }
  }

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();
}

extern void ParseInNucUpdates (IteM i)

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
  Message (MSG_OK, "Not yet implemented");
}

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

typedef struct recompdata {
  Int4        count;
  MonitorPtr  mon;
  BioseqPtr   batchbsp;
  Boolean     include_stop;
  Boolean     no_stop_at_end_of_complete_cds;
  Boolean     fix_genes;
  Uint2       entityID;
} RecompData, PNTR RecompDataPtr;

static Int2 GeneticCodeFromCrp (CdRegionPtr crp)

{
  Int2            code;
  GeneticCodePtr  gcp;
  Char            name [256];
  ValNodePtr      tmp;

  code = 0;
  name [0] = '\0';
  gcp = crp->genetic_code;
  if (gcp != NULL) {
    tmp = (ValNodePtr) gcp->data.ptrvalue;
    for (tmp = (ValNodePtr) gcp->data.ptrvalue; tmp != NULL; tmp = tmp->next) {
      switch (tmp->choice) {
        case 1 :
          if (name [0] == '\0') {
            StringNCpy_0 (name, (CharPtr) tmp->data.ptrvalue, sizeof (name));
          }
          break;
        case 2 :
          code = tmp->data.intvalue;
          break;
        default :
          break;
      }
    }
    if (code == 0) {
      gcp = GeneticCodeFind (code, name);
      if (gcp != NULL) {
        for (tmp = (ValNodePtr) gcp->data.ptrvalue; tmp != NULL; tmp = tmp->next) {
          switch (tmp->choice) {
            case 2 :
              code = tmp->data.intvalue;
              break;
            default :
              break;
          }
        }
      }
    }
  }
  return code;
}

extern void ExtendSeqLocToPosition (SeqLocPtr slp, Boolean end5, Int4 pos)
{
  Uint1          strand;
  SeqLocPtr      slp_to_change, slp_index;
  Int4           extent_to_change;
  Int4           start, stop;
  SeqIdPtr       sip;
  BioseqPtr      bsp;
  
  if (slp == NULL || pos < 0) return;
  
  bsp = BioseqFindFromSeqLoc (slp);
  if (bsp == NULL) return;

  slp_to_change = NULL;
  strand = SeqLocStrand (slp);
  switch (slp->choice)
  {
    case SEQLOC_INT:
      slp_to_change = slp;
      break;
    case SEQLOC_MIX:
  	case SEQLOC_PACKED_INT:
      sip = SeqLocId (slp);
      if (sip == NULL) return; /* can only process if all on one bioseq */
      slp_to_change = NULL;
      if ((strand == Seq_strand_minus && end5)
        || (strand != Seq_strand_minus && !end5))
      {
        extent_to_change = 0;
        for (slp_index = (SeqLocPtr)slp->data.ptrvalue; slp_index != NULL; slp_index = slp_index->next)
        {
          stop = GetOffsetInBioseq (slp_index, bsp, SEQLOC_STOP);
          if (stop > extent_to_change)
          {
            slp_to_change = slp_index;
            extent_to_change = stop;
          }
        }
      }
      else
      {
        extent_to_change = bsp->length;
        for (slp_index = (SeqLocPtr)slp->data.ptrvalue; slp_index != NULL; slp_index = slp_index->next)
        {
          start = GetOffsetInBioseq (slp_index, bsp, SEQLOC_START);
          if (start < extent_to_change)
          {
            slp_to_change = slp_index;
            extent_to_change = start;
          }
        }
      }
      break;
  }

  if (slp_to_change != NULL)
  {
    if ((strand == Seq_strand_minus && end5)
      || (strand != Seq_strand_minus && !end5))
    {
      start = GetOffsetInBioseq (slp_to_change, bsp, SEQLOC_START);
      stop = pos;
    }
    else
    {
      start = pos;
      stop = GetOffsetInBioseq (slp_to_change, bsp, SEQLOC_STOP);
    }
    if (start < 0 
        || stop > bsp->length - 1
        || start > stop)
    {
      return;
    }
    expand_seq_loc (start, stop, strand, slp_to_change);
  }
}

extern Int4 ExtendSeqLocToEnd (SeqLocPtr slp, BioseqPtr bsp, Boolean end5)
{
  Uint1          strand;
  SeqLocPtr      slp_to_change, slp_index;
  Int4           extent_to_change;
  Int4           start, stop;
  SeqIdPtr       sip;
  Int4           start_diff = 0;
  
  if (slp == NULL || bsp == NULL) return;

  slp_to_change = NULL;
  strand = SeqLocStrand (slp);
  switch (slp->choice)
  {
    case SEQLOC_INT:
      slp_to_change = slp;
      break;
    case SEQLOC_MIX:
  	case SEQLOC_PACKED_INT:
      sip = SeqLocId (slp);
      if (sip == NULL) return; /* can only process if all on one bioseq */
      slp_to_change = NULL;
      if ((strand == Seq_strand_minus && end5)
        || (strand != Seq_strand_minus && !end5))
      {
        extent_to_change = 0;
        for (slp_index = (SeqLocPtr)slp->data.ptrvalue; slp_index != NULL; slp_index = slp_index->next)
        {
          stop = GetOffsetInBioseq (slp_index, bsp, SEQLOC_STOP);
          if (stop > extent_to_change)
          {
            slp_to_change = slp_index;
            extent_to_change = stop;
          }
        }
      }
      else
      {
        extent_to_change = bsp->length;
        for (slp_index = (SeqLocPtr)slp->data.ptrvalue; slp_index != NULL; slp_index = slp_index->next)
        {
          start = GetOffsetInBioseq (slp_index, bsp, SEQLOC_START);
          if (start < extent_to_change)
          {
            slp_to_change = slp_index;
            extent_to_change = start;
          }
        }
      }
      break;
  }

  if (slp_to_change != NULL)
  {
    if ((strand == Seq_strand_minus && end5)
      || (strand != Seq_strand_minus && !end5))
    {
      start = GetOffsetInBioseq (slp_to_change, bsp, SEQLOC_START);
      stop = bsp->length - 1;
    }
    else
    {
      start = 0;
      stop = GetOffsetInBioseq (slp_to_change, bsp, SEQLOC_STOP);
    }
    if (end5) {
        if (strand == Seq_strand_minus) {
            start_diff = bsp->length - 1 - GetOffsetInBioseq(slp_to_change, bsp, SEQLOC_START);
        } else {
            start_diff = GetOffsetInBioseq(slp_to_change, bsp, SEQLOC_START);
        }
    }
    
    expand_seq_loc (start, stop, strand, slp_to_change);
  }
  return start_diff;
}

static void ExtendOnePartialFeature (SeqFeatPtr sfp, Pointer userdata)
{
  BioseqPtr   bsp;
  Boolean     partial3, partial5;
  Int4        start_diff;
  CdRegionPtr crp;

  if (sfp == NULL) return;
  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return;
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  if (partial5)
  {
    start_diff = ExtendSeqLocToEnd (sfp->location, bsp, TRUE);
    if (start_diff > 0 && sfp->data.choice == SEQFEAT_CDREGION) {
      crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      if (crp != NULL) {
          if (crp->frame == 0) {
              crp->frame = 1;
          }
          crp->frame = (crp->frame + start_diff - 1) % 3 + 1;
      }
    }
  }
  if (partial3)
  {
    ExtendSeqLocToEnd (sfp->location, bsp, FALSE);
  }
}

extern void ExtendPartialFeatures (IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep, old_scope;
  SelStructPtr      sel;
  SeqFeatPtr        sfp;
  SeqMgrFeatContext fcontext;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  sel = ObjMgrGetSelected ();
  WatchCursor ();
  Update ();
  old_scope = SeqEntrySetScope(sep);
  if (sel == NULL)
  {
    VisitFeaturesInSep (sep, NULL, ExtendOnePartialFeature);
  }
  else
  {
    while (sel != NULL)
    {
      if (sel->entityID == bfp->input_entityID
        && sel->itemtype == OBJ_SEQFEAT)
      {
        sfp = SeqMgrGetDesiredFeature (bfp->input_entityID, NULL, sel->itemID, 0, NULL, &fcontext);
        if (sfp != NULL)
        {
          ExtendOnePartialFeature (sfp, NULL);
        }
      }
      sel = sel->next;
    }
  }
  SeqEntrySetScope(old_scope);
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static Boolean HasValidStartCodon (SeqFeatPtr cds)
{
  ByteStorePtr  bs;
  CharPtr       prot;

  if (cds == NULL) return FALSE;

  bs = ProteinFromCdRegionEx (cds, TRUE, FALSE);

  if (bs == NULL) return FALSE;
  prot = BSMerge (bs, NULL);
  bs = BSFree (bs);
  if (prot == NULL) return FALSE;
  if (prot [0] != 'M') return FALSE;
  return TRUE;
}

static void FixReadingFrame (
  CdRegionPtr crp,
  SeqLocPtr slp,
  BioseqPtr bsp,
  Int4 start
)
{
  Int4 offset;

  if (crp == NULL || slp == NULL) return;

  if (SeqLocStrand (slp) == Seq_strand_minus)
  {
    offset = bsp->length - start;
  }
  else
  {
    offset = start;
  }
  start = start % 3;
  start = start + 1;
  crp->frame = start;
}

extern void RecomputeSuggestedIntervalsForCDS 
(Uint2          entityID,
 BioseqPtr PNTR batchbsp,
 Int4Ptr        count,
 MonitorPtr     mon,
 SeqFeatPtr     sfp)
{
  Int2           code;
  CdRegionPtr    crp;
  BioseqPtr      nucbsp;
  BioseqPtr      protbsp;
  SeqIdPtr       sip;
  SeqLocPtr      slp;
  Char           str [256];
  Char           tmp [256];
  Boolean        partial3, partial5;
  Int4           start;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return;
  crp = (CdRegionPtr) sfp->data.value.ptrvalue;
  if (crp == NULL) return;

  code = GeneticCodeFromCrp (crp);

  nucbsp = GetBioseqGivenSeqLoc (sfp->location, entityID);
  if (nucbsp != NULL && batchbsp != NULL && *batchbsp != NULL 
      && nucbsp != *batchbsp) {
    ClearBatchSuggestNucleotide ();
    *batchbsp = nucbsp;
    SetBatchSuggestNucleotide (*batchbsp, code);
/*    Message (MSG_POSTERR, "Recompute Suggest is reverting to slower processing"); */
  }
  sip = SeqLocId (sfp->product);
  if (sip != NULL) {
    protbsp = BioseqFind (sip);
    if (nucbsp != NULL && protbsp != NULL &&
        ISA_na (nucbsp->mol) && ISA_aa (protbsp->mol) &&
        nucbsp->length > 0 && protbsp->length > 0) {
      str [0] = '\0';
      tmp [0] = '\0';
      sip = SeqIdFindWorst (protbsp->id);
      SeqIdWrite (sip, tmp, PRINTID_REPORT, sizeof (tmp));
      if (count != NULL)
      {
        (*count) ++;
        if (mon != NULL)
        {
          sprintf (str, "Processing sequence %d [%s]", *count, tmp);
          MonitorStrValue (mon, str);
          Update ();
        }
      }
      slp = PredictCodingRegion (nucbsp, protbsp, code);
      if (slp == NULL) return;

      /* correct for partial conditions */
      CheckSeqLocForPartial (sfp->location, &partial5, &partial3);

      sfp->location = SeqLocFree (sfp->location);
      sfp->location = slp;

      /* if no valid start codon, cds is 5' partial */
      if (! HasValidStartCodon (sfp))
      {
        partial5 = TRUE;
      }

      if (partial5 || partial3)
      {
        if (partial5)
        {
          start = GetOffsetInBioseq (sfp->location, nucbsp, SEQLOC_START);
          ExtendSeqLocToEnd (sfp->location, nucbsp, TRUE);
          FixReadingFrame (crp, sfp->location, nucbsp, start);
        }
        if (partial3)
        {
          ExtendSeqLocToEnd (sfp->location, nucbsp, FALSE);
        }
        SetSeqLocPartial (sfp->location, partial5, partial3);
      }
      sfp->partial = LocationHasNullsBetween (sfp->location);
      sfp->partial |= partial5 || partial3;
    }
  }
}

static void RecomputeIntervalsForOneCDS (SeqFeatPtr sfp, RecompDataPtr rdp)
{
  SeqFeatPtr gene_to_update = NULL;
  SeqLocPtr      orig_loc = NULL;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION || rdp == NULL)
  {
    return;
  }

  if (rdp->fix_genes)
  {
    gene_to_update = SeqMgrGetOverlappingGene (sfp->location, NULL);
    orig_loc = (SeqLocPtr) AsnIoMemCopy (sfp->location, 
                                         (AsnReadFunc) SeqLocAsnRead,
                                         (AsnWriteFunc) SeqLocAsnWrite);
  }
  RecomputeSuggestedIntervalsForCDS (rdp->entityID, &(rdp->batchbsp),
                                     &(rdp->count), rdp->mon, sfp);

  if (gene_to_update != NULL)
  {
    UpdateGeneLocation (gene_to_update, orig_loc,
                        sfp->location, rdp->entityID);
  }
  orig_loc = SeqLocFree (orig_loc);
  
}

static Boolean RecomputeSuggCallback (GatherContextPtr gcp)
{
  RecompDataPtr  rdp;
  SeqFeatPtr     sfp;

  if (gcp == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  rdp = (RecompDataPtr) gcp->userdata;
  if (rdp == NULL) return TRUE;
  sfp = (SeqFeatPtr) gcp->thisitem;
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return TRUE;

  RecomputeIntervalsForOneCDS (sfp, rdp);
  return TRUE;
}

extern void RecomputeSuggestEx (Uint2 entityID, Boolean fix_genes, Boolean recompute_all)

{
  Int2              code;
  GatherScope       gs;
  SeqEntryPtr       nucsep;
  RecompData        rd;
  SeqEntryPtr       sep;
  SelStructPtr      sel;
  SeqFeatPtr        sfp;
  SeqMgrFeatContext fcontext;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return;
  sel = ObjMgrGetSelected ();
  WatchCursor ();
  Update ();
  rd.count = 0;
  rd.mon = MonitorStrNewEx ("Correcting Coding Regions", 20, FALSE);
  rd.batchbsp = NULL;
  rd.no_stop_at_end_of_complete_cds = FALSE;
  rd.fix_genes = fix_genes;
  rd.entityID = entityID;
  nucsep = FindNucSeqEntry (sep);
  if (nucsep != NULL && IS_Bioseq (nucsep)) {
    rd.batchbsp = (BioseqPtr) nucsep->data.ptrvalue;
  }
  if (rd.batchbsp != NULL) {
    code = SeqEntryToGeneticCode (sep, NULL, NULL, 0);
    SetBatchSuggestNucleotide (rd.batchbsp, code);
  }
  if (sel == NULL || recompute_all)
  {
    MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
    gs.seglevels = 1;
    gs.get_feats_location = FALSE;
    MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
    gs.ignore[OBJ_BIOSEQ] = FALSE;
    gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
    gs.ignore[OBJ_SEQFEAT] = FALSE;
    gs.ignore[OBJ_SEQANNOT] = FALSE;
    GatherEntity (entityID, (Pointer) (&rd), RecomputeSuggCallback, &gs);
  }
  else
  {
    while (sel != NULL)
    {
      if (sel->entityID == entityID
        && sel->itemtype == OBJ_SEQFEAT)
      {
        sfp = SeqMgrGetDesiredFeature (entityID, NULL, sel->itemID, 0, NULL, &fcontext);
        if (sfp != NULL && sfp->idx.subtype == FEATDEF_CDS)
        {
          RecomputeIntervalsForOneCDS (sfp, &rd);
        }
      }
      sel = sel->next;
    }
  }
  MonitorFree (rd.mon);
  if (rd.batchbsp != NULL) {
    ClearBatchSuggestNucleotide ();
  }
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
}

extern void RecomputeSuggest (IteM i)
{
  BaseFormPtr       bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL)
  {
    return;
  }
  RecomputeSuggestEx (bfp->input_entityID, FALSE, FALSE);
}

extern void RecomputeSuggestFixGenes (IteM i)
{
  BaseFormPtr       bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL)
  {
    return;
  }
  RecomputeSuggestEx (bfp->input_entityID, TRUE, FALSE);
}


extern Boolean RetranslateOneCDS 
( SeqFeatPtr sfp,
  Uint2 entityID,
  Boolean include_stop,
  Boolean no_stop_at_end_of_complete_cds)

{
  SeqFeatPtr    bestprot;
  ByteStorePtr  bs;
  BioseqPtr     bsp;
  Char          ch;
  SeqFeatPtr    gene;
  GeneRefPtr    grp;
  SeqEntryPtr   master;
  MolInfoPtr    mip;
  SeqEntryPtr   old;
  Boolean       partial5;
  Boolean       partial3;
  CharPtr       prot;
  CharPtr       ptr;
  SeqEntryPtr   sep;
  SeqIdPtr      sip;
  ValNodePtr    vnp;
  ProtRefPtr    prp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return TRUE;

  /* bail on pseudo CDS */

  if (sfp->pseudo) return TRUE;
  grp = SeqMgrGetGeneXref (sfp);
  if (grp != NULL) {
    if (grp->pseudo) return TRUE;
  } else {
    gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
    if (gene != NULL) {
      if (gene->pseudo) return TRUE;
      grp = (GeneRefPtr) gene->data.value.ptrvalue;
      if (grp != NULL && grp->pseudo) return TRUE;
    }
  }

  if (sfp->location == NULL) return TRUE;
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);

  if (sfp->product == NULL) {
    master = NULL;
    old = NULL;
    bsp = GetBioseqGivenSeqLoc (sfp->location, entityID);
    if (bsp != NULL) {
      master = GetBestTopParentForData (entityID, bsp);
    }
    bsp = BioseqNew ();
    if (bsp != NULL) {
      bsp->mol = Seq_mol_aa;
      bsp->repr = Seq_repr_raw;
      bsp->seq_data_type = Seq_code_ncbieaa;
      bsp->length = 0;
      bsp->seq_data = BSNew (0);
      if (master != NULL) {
        old = SeqEntrySetScope (master);
      }
      bsp->id = MakeNewProteinSeqId (sfp->location, NULL);
      SeqMgrAddToBioseqIndex (bsp);
      if (master != NULL) {
        SeqEntrySetScope (old);
      }
      sep = SeqEntryNew ();
      if (sep != NULL) {
        sep->choice = 1;
        sep->data.ptrvalue = (Pointer) bsp;
        SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);
      }
      SetSeqFeatProduct (sfp, bsp);
      if (master != NULL && sep != NULL) {
        AddSeqEntryToSeqEntry (master, sep, TRUE);
      }
    }
  }

  sip = SeqLocId (sfp->product);
  if (sip != NULL) {
    bsp = BioseqFind (sip);
    if (bsp != NULL && ISA_aa (bsp->mol) && bsp->repr == Seq_repr_raw) {
      bestprot = FindBestProtein (entityID, sfp->product);
      bs = ProteinFromCdRegionExWithTrailingCodonHandling (sfp,
                                              include_stop,
                                              FALSE,
                                              no_stop_at_end_of_complete_cds );
      if (bs == NULL) return TRUE;
      prot = BSMerge (bs, NULL);
      bs = BSFree (bs);
      if (prot == NULL) return TRUE;
      ptr = prot;
      ch = *ptr;
      while (ch != '\0') {
        *ptr = TO_UPPER (ch);
        ptr++;
        ch = *ptr;
      }
      bs = BSNew (1000);
      if (bs != NULL) {
        ptr = prot;
        /*
        if (prot [0] == '-') {
          ptr++;
        }
        */
        BSWrite (bs, (VoidPtr) ptr, (Int4) StringLen (ptr));
      }
      bsp->repr = Seq_repr_raw;
      bsp->mol = Seq_mol_aa;
      bsp->seq_data_type = Seq_code_ncbieaa;
      bsp->seq_data = BSFree (bsp->seq_data);
      bsp->seq_data = bs;
      bsp->length = BSLen (bs);
      sep = SeqMgrGetSeqEntryForData (bsp);
      if (sep == NULL) return TRUE;
      if (bestprot == NULL)
      {
        bestprot = CreateNewFeature (sep, NULL, SEQFEAT_PROT, NULL);
        prp = ProtRefNew ();
        bestprot->data.value.ptrvalue = prp;
      }
      if (bestprot != NULL) {
        bestprot->location = SeqLocFree (bestprot->location);
        bestprot->location = CreateWholeInterval (sep);
        SetSeqLocPartial (bestprot->location, partial5, partial3);
        bestprot->partial = (partial5 || partial3);
      }
      vnp = SeqEntryGetSeqDescr (sep, Seq_descr_molinfo, NULL);
      if (vnp == NULL) {
        vnp = CreateNewDescriptor (sep, Seq_descr_molinfo);
        if (vnp != NULL) {
          mip = MolInfoNew ();
          vnp->data.ptrvalue = (Pointer) mip;
          if (mip != NULL) {
            mip->biomol = 8;
            mip->tech = 13;
          }
        }
      }
      if (vnp != NULL) {
        mip = (MolInfoPtr) vnp->data.ptrvalue;
        if (mip != NULL) {
          if (partial5 && partial3) {
            mip->completeness = 5;
          } else if (partial5) {
            mip->completeness = 3;
          } else if (partial3) {
            mip->completeness = 4;
          /*
          } else if (partial) {
            mip->completeness = 2;
          */
          } else {
            mip->completeness = 0;
          }
        }
      }
    }
  }
  return TRUE;
}

static Boolean RetranslateCDSCallback (GatherContextPtr gcp)

{
  RecompDataPtr  rdp;
  SeqFeatPtr     sfp;

  if (gcp == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  rdp = (RecompDataPtr) gcp->userdata;
  if (rdp == NULL) return TRUE;
  sfp = (SeqFeatPtr) gcp->thisitem;
  return RetranslateOneCDS (sfp, gcp->entityID, rdp->include_stop,
                            rdp->no_stop_at_end_of_complete_cds);
}

extern void RetranslateCdRegionsEx (
  Uint2   entityID,
  Boolean include_stop,
  Boolean no_stop_at_end_of_complete_cds )

{
  GatherScope  gs;
  RecompData   rd;
  SeqEntryPtr  sep;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return;
  WatchCursor ();
  Update ();
  rd.count = 0;
  rd.mon = MonitorStrNewEx ("Correcting Coding Regions", 20, FALSE);
  rd.include_stop = include_stop;
  rd.no_stop_at_end_of_complete_cds = no_stop_at_end_of_complete_cds;
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = FALSE;
  MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  GatherEntity (entityID, (Pointer) (&rd), RetranslateCDSCallback, &gs);
  MonitorFree (rd.mon);
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
}

static void RetranslateCdRegions (
  IteM i,
  Boolean include_stop,
  Boolean no_stop_at_end_of_complete_cds )

{
  BaseFormPtr  bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  RetranslateCdRegionsEx (bfp->input_entityID, include_stop, no_stop_at_end_of_complete_cds);
}

extern void RetranslateCdRegionsNoStop (IteM i)

{
  RetranslateCdRegions (i, FALSE, FALSE);
}

extern void RetranslateCdRegionsDoStop (IteM i)

{
  RetranslateCdRegions (i, TRUE, FALSE);
}

extern void RetranslateCdRegionsNoStopExceptEndCompleteCDS (IteM i)
{
  RetranslateCdRegions (i, TRUE, TRUE);
}


static void DoReprocessPeptides (SeqFeatPtr sfp, Pointer userdata)

{
  SeqFeatPtr    bestprot;
  ByteStorePtr  bs;
  BioseqPtr     bsp;
  Char          ch;
  MolInfoPtr    mip;
  Boolean       partial5;
  Boolean       partial3;
  CharPtr       prot;
  ProtRefPtr    prp;
  CharPtr       ptr;
  SeqEntryPtr   sep;
  SeqIdPtr      sip;
  ValNodePtr    vnp;

  if (sfp->data.choice != SEQFEAT_PROT) return;
  if (sfp->product == NULL) return;
  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  if (prp == NULL) return;
  if (prp->processed < 1 || prp->processed > 4) return;
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  sip = SeqLocId (sfp->product);
  if (sip == NULL) return;
  bsp = BioseqFind (sip);
  if (bsp != NULL && ISA_aa (bsp->mol) && bsp->repr == Seq_repr_raw) {
    bestprot = FindBestProtein (sfp->idx.entityID, sfp->product);
    prot = GetSequenceByFeature (sfp);
    if (prot == NULL) return;
    ptr = prot;
    ch = *ptr;
    while (ch != '\0') {
      *ptr = TO_UPPER (ch);
      ptr++;
      ch = *ptr;
    }
    bs = BSNew (1000);
    if (bs != NULL) {
      ptr = prot;
      BSWrite (bs, (VoidPtr) ptr, (Int4) StringLen (ptr));
    }
    bsp->repr = Seq_repr_raw;
    bsp->mol = Seq_mol_aa;
    bsp->seq_data_type = Seq_code_ncbieaa;
    bsp->seq_data = BSFree (bsp->seq_data);
    bsp->seq_data = bs;
    bsp->length = BSLen (bs);
    sep = SeqMgrGetSeqEntryForData (bsp);
    if (sep == NULL) return;
    if (bestprot != NULL) {
      bestprot->location = SeqLocFree (bestprot->location);
      bestprot->location = CreateWholeInterval (sep);
      SetSeqLocPartial (bestprot->location, partial5, partial3);
      bestprot->partial = (partial5 || partial3);
    }
    vnp = SeqEntryGetSeqDescr (sep, Seq_descr_molinfo, NULL);
    if (vnp == NULL) {
      vnp = CreateNewDescriptor (sep, Seq_descr_molinfo);
      if (vnp != NULL) {
        mip = MolInfoNew ();
        vnp->data.ptrvalue = (Pointer) mip;
        if (mip != NULL) {
          mip->biomol = 8;
          mip->tech = 13;
        }
      }
    }
    if (vnp != NULL) {
      mip = (MolInfoPtr) vnp->data.ptrvalue;
      if (mip != NULL) {
        if (partial5 && partial3) {
          mip->completeness = 5;
        } else if (partial5) {
          mip->completeness = 3;
        } else if (partial3) {
          mip->completeness = 4;
        /*
        } else if (partial) {
          mip->completeness = 2;
        */
        } else {
          mip->completeness = 0;
        }
      }
    }
  }
}

extern void ReprocessPeptideProducts (IteM i);
extern void ReprocessPeptideProducts (IteM i)

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
  VisitFeaturesInSep (sep, NULL, DoReprocessPeptides);
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void DoReprocessMrnas (SeqFeatPtr sfp, Pointer userdata)

{
  ByteStorePtr  bs;
  BioseqPtr     bsp;
  Char          ch;
  MolInfoPtr    mip;
  Boolean       partial5;
  Boolean       partial3;
  CharPtr       prot;
  CharPtr       ptr;
  RnaRefPtr     rrp;
  SeqEntryPtr   sep;
  SeqIdPtr      sip;
  ValNodePtr    vnp;

  if (sfp->data.choice != SEQFEAT_RNA) return;
  if (sfp->product == NULL) return;
  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp == NULL) return;
  if (rrp->type != 2) return;
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  sip = SeqLocId (sfp->product);
  if (sip == NULL) return;
  bsp = BioseqFind (sip);
  if (bsp != NULL && ISA_na (bsp->mol) && bsp->repr == Seq_repr_raw) {
    prot = GetSequenceByFeature (sfp);
    if (prot == NULL) return;
    ptr = prot;
    ch = *ptr;
    while (ch != '\0') {
      *ptr = TO_UPPER (ch);
      ptr++;
      ch = *ptr;
    }
    bs = BSNew (1000);
    if (bs != NULL) {
      ptr = prot;
      BSWrite (bs, (VoidPtr) ptr, (Int4) StringLen (ptr));
    }
    bsp->repr = Seq_repr_raw;
    bsp->mol = Seq_mol_na;
    bsp->seq_data_type = Seq_code_iupacna;
    bsp->seq_data = BSFree (bsp->seq_data);
    bsp->seq_data = bs;
    bsp->length = BSLen (bs);
    sep = SeqMgrGetSeqEntryForData (bsp);
    if (sep == NULL) return;
    vnp = SeqEntryGetSeqDescr (sep, Seq_descr_molinfo, NULL);
    if (vnp == NULL) {
      vnp = CreateNewDescriptor (sep, Seq_descr_molinfo);
      if (vnp != NULL) {
        mip = MolInfoNew ();
        vnp->data.ptrvalue = (Pointer) mip;
        if (mip != NULL) {
          mip->biomol = 8;
          mip->tech = 13;
        }
      }
    }
    if (vnp != NULL) {
      mip = (MolInfoPtr) vnp->data.ptrvalue;
      if (mip != NULL) {
        if (partial5 && partial3) {
          mip->completeness = 5;
        } else if (partial5) {
          mip->completeness = 3;
        } else if (partial3) {
          mip->completeness = 4;
        /*
        } else if (partial) {
          mip->completeness = 2;
        */
        } else {
          mip->completeness = 0;
        }
      }
    }
  }
}

extern void ReprocessmRNAProducts (IteM i);
extern void ReprocessmRNAProducts (IteM i)

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
  VisitFeaturesInSep (sep, NULL, DoReprocessMrnas);
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

/*
 * ApplyCodeBreakToCDS
 */
 
static Boolean ApplyCodeBreakToCDS (SeqFeatPtr sfp, CharPtr codonStr, Int2 aaNum)
{
  Uint1            aaChar;
  SeqLocPtr        aaSlp;
  Int4             aaPosition;
  SeqPntPtr        aaSpp;
  CharPtr          basePtr;
  CharPtr          bases;
  CodeBreakPtr     cbp;
  CdRegionPtr      crp;
  Int4             dnaLen;
  SeqLocPtr        dnaSlp;
  CodeBreakPtr     lastCbp;
  SeqCodeTablePtr  sctp;
  Int4             total;
  Boolean          added_code_breaks = FALSE;
  
  if (sfp == NULL || codonStr == NULL)
  {
    return FALSE;
  }

  /* Get the nucleotide sequence */

  dnaLen = SeqLocLen (sfp->location);
  if (dnaLen < 1)
    return FALSE;

  crp = (CdRegionPtr) sfp->data.value.ptrvalue;

  bases = ReadCodingRegionBases (sfp->location, dnaLen, crp->frame, &total);

  /* Search for the selected codon in the */
  /* nucleotide sequence.  If found, add  */
  /* it as a codebreak.                   */

  basePtr = bases;
  aaPosition = 0;
  while (basePtr[0] != '\0') {
    if (StringNCmp (basePtr, codonStr, 3) == 0) {

      /* Create a new seq point object with the aa location */

      aaSpp = SeqPntNew ();
      aaSpp->point  = aaPosition;
      aaSpp->strand = Seq_strand_plus;
      aaSpp->id = SeqLocId (sfp->product);
      aaSpp->fuzz   = NULL;

      /* Make a SeqLoc using the seq point */

      aaSlp = (SeqLocPtr) ValNodeNew (NULL);
      aaSlp->choice = SEQLOC_PNT;
      aaSlp->data.ptrvalue = (Pointer) aaSpp;

      /* Convert the seqloc to a DNA location */

      dnaSlp = aaLoc_to_dnaLoc (sfp, aaSlp);

      /* Create the code break using the DNA location */

      cbp = CodeBreakNew ();
      cbp->loc = dnaSlp;
      sctp = SeqCodeTableFind (Seq_code_ncbieaa);
      aaChar = (Uint1) GetSymbolForResidue (sctp, aaNum);
      cbp->aa.value.intvalue = aaChar;
      cbp->aa.choice = 1; /* ncbieaa */

      /* Insert the code break into the CDS's */
      /* existing list of code breaks.        */

      lastCbp = crp->code_break;
      if (lastCbp == NULL)
      {
	      crp->code_break = cbp;        
      }
      else 
      {
        while (lastCbp->next != NULL)
        {
          lastCbp = lastCbp->next;          
        }
	      lastCbp->next = cbp;
	      cbp->next = NULL;
      }

      added_code_breaks = TRUE;
    }
    basePtr += 3;
    aaPosition++;
  }

  return added_code_breaks;
}

/*---------------------------------------------------------------------*/
/*                                                                     */
/* ApplyCodeBreak_FeatureCallback () -- Called for each CDS feature in */
/*                                      a Bioseq.  Checks for any      */
/*                                      nucleotide triplets that match */
/*                                      the one in the given code      */
/*                                      break and sets a code break    */
/*                                      for each one that is found.    */
/*                                                                     */
/*---------------------------------------------------------------------*/

static Boolean LIBCALLBACK ApplyCodeBreak_FeatureCallback (SeqFeatPtr sfp,
					SeqMgrFeatContextPtr fcontext)
{
  Int2             aaNum;
  CodeBreakFormPtr cbfp;
  Char             codonStr [4];
  Int2             i;

  cbfp = (CodeBreakFormPtr) fcontext->userdata;


  /* Get the selected Amino Acid and codon triplet */

  GetTitle (cbfp->codonText, codonStr, sizeof (codonStr));
  for (i = 0; i < 3; i++)
    codonStr [i] = TO_UPPER (codonStr [i]);

  aaNum = GetValue (cbfp->aminoAcidPopup);
  aaNum += 63;

  if (aaNum >= 74)
  {
  	aaNum++;
  }
  if (aaNum >= 79)
  {
  	aaNum++;
  }

  if (ApplyCodeBreakToCDS (sfp, codonStr, aaNum))
  {
    /* Retranslate the CDS */

    RetranslateOneCDS (sfp, fcontext->entityID, TRUE, FALSE);
    
  }
  
  /* Return TRUE to continue on to the next CDS feature */

  return TRUE;
}

/*---------------------------------------------------------------------*/
/*                                                                     */
/* ApplyCodeBreak_BioseqCallback () -- Called by SeqMgrExploreBioseqs  */
/*                                     for each Bioseq.  Searches the  */
/*                                     Bioseq for CDS features and adds*/
/*                                     the given code break to any     */
/*                                     found.                          */
/*                                                                     */
/*---------------------------------------------------------------------*/

static Boolean LIBCALLBACK ApplyCodeBreak_BioseqCallback (BioseqPtr bsp,
					 SeqMgrBioseqContextPtr bcontext)
{
  Boolean featureFilterArray [SEQFEAT_MAX];

  /* Set up to explore only CDS features */

  MemSet ((Pointer) (featureFilterArray),
	  (int) FALSE,
	  SEQFEAT_MAX);

  featureFilterArray[SEQFEAT_CDREGION] = TRUE;

  /* Explore the Bioseq's CDS features, marking the */
  /* ones with internal stop codons as pseudo.      */

  SeqMgrExploreFeatures (bsp, bcontext->userdata,
			 ApplyCodeBreak_FeatureCallback, NULL,
			 featureFilterArray, NULL);

  /* Return TRUE to continue on to the next Bioseq */

  return TRUE;
}

/*---------------------------------------------------------------------*/
/*                                                                     */
/* DoAddCodeBreak_Callback () -- Called when the 'Apply' button is     */
/*                               pressed in the "Add Code Break"       */
/*                               window.  Adds the entered code break  */
/*                               to all CDS features.                  */
/*                                                                     */
/*---------------------------------------------------------------------*/

static void DoAddCodeBreak_Callback (ButtoN b)
{
  CodeBreakFormPtr  cbfp;

  cbfp = (CodeBreakFormPtr) GetObjectExtra (b);

  /* Change to the "working" cursor */

  Hide (cbfp->form);
  WatchCursor ();
  Update ();

  /* Visit all the Bioseqs, where we will */
  /* then explore their CDS features.     */

  SeqMgrExploreBioseqs (cbfp->input_entityID, NULL, (Pointer) cbfp,
			ApplyCodeBreak_BioseqCallback, TRUE, FALSE, TRUE);

  /* Restore the cursor and force an update */

  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (cbfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, cbfp->input_entityID, 0, 0);
  Remove (cbfp->form);
}

/*---------------------------------------------------------------------*/
/*                                                                     */
/* PopulateAAPopup () -- Creates a popup list of amino acids.          */
/*                                                                     */
/*     NOTE : This function is identical to (and identically named as) */
/*            a function in cdrgn.c                                    */
/*                                                                     */
/*---------------------------------------------------------------------*/

static void PopulateAAPopup (PopuP AAitem)

{
  Char             ch;
  Uint1            first;
  Uint1            i;
  Char             item [77];
  Uint1            last;
  SeqCodeTablePtr  sctp;
  CharPtr          str;

  sctp = SeqCodeTableFind (Seq_code_ncbieaa);
  first = FirstResidueInCode (sctp);
  last = LastResidueInCode (sctp);
  PopupItem (AAitem, " ");
  for (i = 65; i <= last; i++) {
    if (i == 74 || i == 79) {
      continue;
    }
    ch = GetSymbolForResidue (sctp, i);
    str = (CharPtr) GetNameForResidue (sctp, i);
    sprintf (item, "%c    %s", ch, str);
    PopupItem (AAitem, item);
  } 
  SetValue (AAitem, 1); 
}

/*---------------------------------------------------------------------*/
/*                                                                     */
/* IsLegalCodon () - Determines if a three base string is a legal      */
/*                   codon.                                            */
/*                                                                     */
/*---------------------------------------------------------------------*/

static Boolean IsLegalCodon (CharPtr codonStr)
{
  Int2 i;
  Char baseChar;

  /* Only allow three characters */

  if (StringLen (codonStr) > 3)
    return FALSE;

  /* Allow only the character A,C,G,T,U and */
  /* convert the U to a T.                  */

  i = 0;
  while (i < 3) {

    if (codonStr [i] == '\0')
      break;

    baseChar = codonStr [i];
    
    if (StringChr ("acgtuACGTU", baseChar) == NULL)
      return FALSE;
    if ('U' == baseChar)
      codonStr [i] = 'T';
    else if ('u' == baseChar)
      codonStr [i] = 't';
    
    i++;
  }

  /* If we made it this far, it's a valid codon */

  return TRUE;
}

/*---------------------------------------------------------------------*/
/*                                                                     */
/* CodonText_Callback () -- Called whenever a keystoke is entered in   */
/*                          Codon text field.  Validates to see if the */
/*                          keystroke should be allowed.               */
/*                                                                     */
/*---------------------------------------------------------------------*/

static void CodonText_Callback (TexT codonText)

{
  CodeBreakFormPtr  cbfp;
  Int2              aaNum;
  Char              newCodonStr [5];

  /* Get the currect code break data */

  cbfp = (CodeBreakFormPtr) GetObjectExtra (codonText);
  if (cbfp == NULL)
    return;

  /* If the new codon string is not legal */
  /* then reset to the previous text.     */

  GetTitle (codonText, newCodonStr, sizeof (newCodonStr));

  if (!IsLegalCodon (newCodonStr))
    StringCpy (newCodonStr, cbfp->currentCodonStr);
  else
    StringCpy (cbfp->currentCodonStr, newCodonStr);

  SafeSetTitle (cbfp->codonText, newCodonStr);

  /* Only enable the accept button if */
  /* we have a full codon.            */

  if (StringLen (newCodonStr) != 3) {
    SafeDisable (cbfp->acceptButton);
    return;
  }

  /* See if an amino acid has been selected yet */

  aaNum = GetValue (cbfp->aminoAcidPopup);
  if (aaNum <= 1) {
    SafeDisable (cbfp->acceptButton);
    return;
  }

  /* If we made it this far then we have both a codon and */
  /* an amino acid, so enable the accept button.          */
 
  SafeEnable (cbfp->acceptButton);
}

/*---------------------------------------------------------------------*/
/*                                                                     */
/* SelectAminoAcid_Callback () -- Called whenever a new amino acid is  */
/*                                selected in the Amino Acid Popup.    */
/*                                Toggles 'Accept' button base on      */
/*                                current state.                       */
/*                                                                     */
/*---------------------------------------------------------------------*/

static void SelectAminoAcid_Callback (PopuP p)
{
  CodeBreakFormPtr  cbfp;
  Char              codonStr [4];
  Int2              aaNum;

  /* Get the currect code break data */

  cbfp = (CodeBreakFormPtr) GetObjectExtra (p);
  if (cbfp == NULL)
    return;

  /* Only enable the accept button if */
  /* we have a full codon.            */

  GetTitle (cbfp->codonText, codonStr, sizeof (codonStr));
  if (StringLen (codonStr) != 3) {
    SafeDisable (cbfp->acceptButton);
    return;
  }

  /* Get the newly selected amino acid */

  aaNum = GetValue (cbfp->aminoAcidPopup);

  /* If an amino acid is selected then */
  /* enable the accept button.         */

  if (aaNum > 1)
    SafeEnable (cbfp->acceptButton);
  else
    SafeDisable (cbfp->acceptButton);
}

/*---------------------------------------------------------------------*/
/*                                                                     */
/* AddGlobalCodeBreak () -- Gets a nucleotide triplet and an amino     */
/*                          acid from the user and adds them as        */
/*                          codebreaks for all CDS features.           */
/*                                                                     */
/*---------------------------------------------------------------------*/

extern void AddGlobalCodeBreak (IteM i);
extern void AddGlobalCodeBreak (IteM i)

{
  BaseFormPtr      bfp;
  WindoW           breakWin;
  GrouP            mainGroup;
  GrouP            buttGroup;
  CodeBreakFormPtr cbfp;

  /* Get the current state of things */

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL)
    return;

  cbfp = (CodeBreakFormPtr) MemNew (sizeof (EvidenceFormData));

  /* Create a window to get the codon and */
  /* the Amino acid from the user.        */

  breakWin = FixedWindow (-50, -33, -10, -10,
			  "Add Code Break", StdCloseWindowProc);
  SetObjectExtra (breakWin, cbfp, StdCleanupFormProc);
  cbfp->form = (ForM) breakWin;
  cbfp->formmessage = NULL;
  cbfp->input_entityID = bfp->input_entityID;

  mainGroup = HiddenGroup (breakWin, -2, 0, NULL);

  /* Create a text entry box for the nucl. codon */

  StaticPrompt (mainGroup, "Triplet Codon", 0, popupMenuHeight,
		programFont, 'l');
  cbfp->codonText = DialogText (mainGroup, "", 3, CodonText_Callback);
  SetObjectExtra (cbfp->codonText, cbfp, NULL);
  cbfp->currentCodonStr [0] = '\0';

  /* Add a Popup list of Amino Acids */

  StaticPrompt (mainGroup, "Amino Acid", 0, popupMenuHeight,
		programFont, 'l');
  cbfp->aminoAcidPopup = PopupList (mainGroup, TRUE,
				    SelectAminoAcid_Callback);
  PopulateAAPopup (cbfp->aminoAcidPopup);
  SetObjectExtra (cbfp->aminoAcidPopup, cbfp, NULL);

  /* Add Accept and Cancel buttons */

  buttGroup = HiddenGroup (breakWin, 2, 0, NULL);
  cbfp->acceptButton = DefaultButton (buttGroup, "Accept",
				   DoAddCodeBreak_Callback);
  SetObjectExtra (cbfp->acceptButton, cbfp, NULL);
  SafeDisable (cbfp->acceptButton);
  PushButton (buttGroup, "Cancel", StdCancelButtonProc);

  /* Line things up and display the window */

  AlignObjects (ALIGN_CENTER, (HANDLE) mainGroup, (HANDLE) buttGroup, NULL);
  RealizeWindow (breakWin);
  Show (breakWin);
  Update ();

}

static void ParseCodonQualToCodeBreakCallback (SeqFeatPtr sfp, Pointer userdata)
{
  SeqCodeTablePtr sctp;
  GBQualPtr       gqual, prev_qual = NULL, next_qual;
  CharPtr         cp;
  Char            codon_text[4];
  Char            symbol_text [4];
  Int2            aaNum;
  Int4            i;
  Boolean         converted_qual;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION || userdata == NULL)
  {
    return;
  }
  
  sctp = (SeqCodeTablePtr) userdata;
  
  for (gqual = sfp->qual; gqual != NULL; gqual = next_qual)
  {
    next_qual = gqual->next;
    converted_qual = FALSE;
    if (StringCmp (gqual->qual, "codon") == 0)
    {
      cp = StringSearch (gqual->val, "seq:\"");
      if (cp != NULL)
      {
        cp += 5;
        StringNCpy (codon_text, cp, 3);
        codon_text [3] = 0;
        for (i = 0; i < 3; i++)
          codon_text [i] = TO_UPPER (codon_text [i]);

        cp = StrChr (cp, ':');
        if (cp != NULL)
        {
          cp++;
          StringNCpy (symbol_text, cp, 3);
          symbol_text [3] = 0;
          aaNum = FindResidueByName (symbol_text, sctp);         
          if (ApplyCodeBreakToCDS (sfp, codon_text, aaNum))
          {
            /* Retranslate the CDS */

            RetranslateOneCDS (sfp, sfp->idx.entityID, TRUE, FALSE);
            
            /* remove the codon qual */
            if (prev_qual == NULL)
            {
              sfp->qual = gqual->next;
            }
            else
            {
              prev_qual->next = gqual->next;
            }
            gqual->next = NULL;
            GBQualFree (gqual);
            converted_qual = TRUE;
          }
        }
      }
    }
    if (!converted_qual)
    {
      prev_qual = gqual;
    }
  }
}

extern void ParseCodonQualToCodeBreak (IteM i)
{
  BaseFormPtr      bfp;
  SeqEntryPtr      sep;
  SeqCodeTablePtr  sctp;
  
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
  
  sctp = SeqCodeTableFind (Seq_code_ncbieaa);
  
  VisitFeaturesInSep (sep, sctp, ParseCodonQualToCodeBreakCallback);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update ();
}

static Boolean CorrectGenCodeCallback (GatherContextPtr gcp)

{
  CdRegionPtr     crp;
  GeneticCodePtr  gc;
  Int2            genCode;
  SeqEntryPtr     sep;
  SeqFeatPtr      sfp;
  ValNodePtr      vnp;

  if (gcp == NULL) return TRUE;
  sep = (SeqEntryPtr) gcp->userdata;
  if (sep == NULL ) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  sfp = (SeqFeatPtr) gcp->thisitem;
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return TRUE;
  crp = (CdRegionPtr) sfp->data.value.ptrvalue;
  if (crp == NULL) return TRUE;
  genCode = SeqEntryToGeneticCode (sep, NULL, NULL, 0);
  gc = GeneticCodeNew ();
  if (gc == NULL) return TRUE;
  crp->genetic_code = GeneticCodeFree (crp->genetic_code);
  vnp = ValNodeNew (NULL);
  gc->data.ptrvalue = vnp;
  if (vnp != NULL) {
    vnp->choice = 2;
    vnp->data.intvalue = (Int4) genCode;
  }
  crp->genetic_code = gc;
  return TRUE;
}

extern void CorrectGenCodes (SeqEntryPtr sep, Uint2 entityID)

{
  BioseqSetPtr  bssp;
  GatherScope   gs;

  if (sep == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 ||
                         (IsPopPhyEtcSet (bssp->_class)))) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        CorrectGenCodes (sep, entityID);
      }
      return;
    }
  }
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = FALSE;
  MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.scope = sep;
  GatherEntity (entityID, (Pointer) sep, CorrectGenCodeCallback, &gs);
}

extern void CorrectCDSGenCodes (IteM i)

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
  CorrectGenCodes (sep, bfp->input_entityID);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static BioseqPtr SqnGetBioseqGivenSeqLoc (SeqLocPtr slp, Uint2 entityID)

{
  BioseqPtr    bsp;
  SeqEntryPtr  sep;
  SeqIdPtr     sip;
  SeqLocPtr    tmp;

  if (slp == NULL) return NULL;
  bsp = NULL;
  sip = SeqLocId (slp);
  if (sip != NULL) {
    bsp = BioseqFind (sip);
  } else {
    tmp = SeqLocFindNext (slp, NULL);
    if (tmp != NULL) {
      sip = SeqLocId (tmp);
      if (sip != NULL) {
        bsp = BioseqFind (sip);
        if (bsp != NULL) {
          sep = SeqMgrGetSeqEntryForData (bsp);
          entityID = ObjMgrGetEntityIDForChoice (sep);
          bsp = GetBioseqGivenSeqLoc (slp, entityID);
        }
      }
    }
  }
  return bsp;
}

static BioseqPtr GetBioseqReferencedByAnnot (SeqAnnotPtr sap, Uint2 entityID)

{
  SeqAlignPtr   align;
  BioseqPtr     bsp;
  DenseDiagPtr  ddp;
  DenseSegPtr   dsp;
  SeqFeatPtr    feat;
  SeqGraphPtr   graph;
  SeqIdPtr      sip;
  SeqLocPtr     slp;
  StdSegPtr     ssp;
  SeqLocPtr     tloc;

  if (sap == NULL) return NULL;
  switch (sap->type) {
    case 1 :
      feat = (SeqFeatPtr) sap->data;
      while (feat != NULL) {
        slp = feat->location;
        if (slp != NULL) {
          bsp = SqnGetBioseqGivenSeqLoc (slp, entityID);
          if (bsp != NULL) return bsp;
        }
        feat = feat->next;
      }
      break;
    case 2 :
      align = (SeqAlignPtr) sap->data;
      while (align != NULL) {
        if (align->segtype == 1) {
          ddp = (DenseDiagPtr) align->segs;
          if (ddp != NULL) {
            for (sip = ddp->id; sip != NULL; sip = sip->next) {
              bsp = BioseqFind (sip);
              if (bsp != NULL) return bsp;
            }
          }
        } else if (align->segtype == 2) {
          dsp = (DenseSegPtr) align->segs;
          if (dsp != NULL) {
            for (sip = dsp->ids; sip != NULL; sip = sip->next) {
              bsp = BioseqFind (sip);
              if (bsp != NULL) return bsp;
            }
          }
        } else if (align->segtype == 3) {
          ssp = (StdSegPtr) align->segs;
          if (ssp != NULL && ssp->loc != NULL) {
            for (tloc = ssp->loc; tloc != NULL; tloc = tloc->next) {
              bsp = BioseqFind (SeqLocId (tloc));
              if (bsp != NULL) return bsp;
            }
          }
        }
        align = align->next;
      }
      break;
    case 3 :
      graph = (SeqGraphPtr) sap->data;
      while (graph != NULL) {
        slp = graph->loc;
        if (slp != NULL) {
          bsp = SqnGetBioseqGivenSeqLoc (slp, entityID);
          if (bsp != NULL) return bsp;
        }
        graph = graph->next;
      }
      break;
    default :
      break;
  }
  return NULL;
}

static Int4 GetScore (ScorePtr score)

{
  ObjectIdPtr  id;

  while (score != NULL) {
    id = score->id;
    if (id != NULL) {
      if (StringICmp (id->str, "score") == 0) {
        if (score->choice == 1) {
          return (score->value.intvalue);
        }
      }
    }
    score = score->next;
  }
  return 0;
}

static Int4 FindScore (SeqAlignPtr align)

{
  if (align == NULL) return 0;
  if (align->score != NULL) {
    return GetScore (align->score);
  }
  return 0;
}

static int LIBCALLBACK SortByScoreCallback (VoidPtr ptr1, VoidPtr ptr2)

{
  SeqAlignPtr   sap1;
  SeqAlignPtr   sap2;
  Int4          score1;
  Int4          score2;

  if (ptr1 != NULL && ptr2 != NULL) {
    sap1 = *((SeqAlignPtr PNTR) ptr1);
    sap2 = *((SeqAlignPtr PNTR) ptr2);
    if (sap1 != NULL && sap2 != NULL) {
      score1 = FindScore (sap1);
      score2 = FindScore (sap2);
      if (score1 < score2) {
        return 1;
      } else if (score1 > score2) {
        return -1;
      } else {
        return 0;
      }
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}

static SeqAlignPtr SortBySeqAlignScore (SeqAlignPtr list)

{
  SeqAlignPtr  align;
  Int4         count, i;
  SeqAlignPtr  PNTR head;

  if (list == NULL) return 0;
  count = 0;
  for (align = list; align != NULL; align = align->next) {
    count++;
  }
  head = MemNew (sizeof (SeqAlignPtr) * (size_t) (count + 1));
  if (head == NULL) return 0;
  for (align = list, i = 0; align != NULL && i < count; i++) {
    head [i] = align;
    align = align->next;
  }
  HeapSort (head, (size_t) count, sizeof (SeqAlignPtr), SortByScoreCallback);
  for (i = 0; i < count; i++) {
    align = head [i];
    align->next = head [i + 1];
  }
  list = head [0];
  MemFree (head);
  return list;
}

static void TakeTop10Alignments (SeqAnnotPtr sap)

{
  SeqAlignPtr  align;
  MsgAnswer    ans;
  Int2         count;
  SeqAlignPtr  next;

  if (sap == NULL || sap->type != 2 || sap->data == NULL) return;
  count = 0;
  for (align = (SeqAlignPtr) sap->data; align != NULL; align = align->next) {
    count++;
  }
  if (count <= 10) return;
  ans = Message (MSG_YN, "Do you want to take only the top 10 (out of %d) alignments?", (int) count);
  if (ans == ANS_NO) return;
  sap->data = SortBySeqAlignScore ((SeqAlignPtr) sap->data);
  for (align = (SeqAlignPtr) sap->data, count = 0; align != NULL && count < 10; align = align->next) {
    count++;
  }
  next = align->next;
  align->next = NULL;
  align = next;
  while (align != NULL) {
    next = align->next;
    align->next = NULL;
    SeqAlignFree (align);
    align = next;
  }
}

static void DoOnePub (PubdescPtr pdp)

{
  ValNodePtr    citartptr = NULL;
  Int4          muid = 0;
  Int4          pmid = 0;
  ValNodePtr    tmp = NULL;
  ValNodePtr    vnp;

  if (pdp != NULL) {
    for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
      if (vnp->choice == PUB_Muid) {
        muid = vnp->data.intvalue;
      } else if (vnp->choice == PUB_PMid) {
        pmid = vnp->data.intvalue;
      } else if (vnp->choice == PUB_Article) {
        citartptr = vnp;
      }
    }
    if (pmid != 0) {
      tmp = MedArchGetPubPmId (pmid);
      muid = MedArchPm2Mu (pmid);
    } else if (muid != 0) {
      tmp = MedArchGetPub (muid);
      pmid = MedArchMu2Pm (muid);
    } else if (citartptr != NULL) {
      muid = MedArchCitMatch (citartptr);
      if (muid != 0) {
        tmp = MedArchGetPub (muid);
        pmid = MedArchMu2Pm (muid);
      }
    }
    if (tmp != NULL) {
      MedlineToISO (tmp);
      if (pmid != 0) {
        ValNodeAddInt (&tmp, PUB_PMid, pmid);
      }
      if (muid != 0) {
        ValNodeAddInt (&tmp, PUB_Muid, muid);
      }
      pdp->pub = PubEquivFree (pdp->pub);
      pdp->pub = tmp;
    }
  }
}

static void DoLookupPub (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  PubdescPtr    pdp;
  SeqAnnotPtr   sap;
  ValNodePtr    sdp;
  SeqFeatPtr    sfp;

  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
    sdp = bsp->descr;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
    sdp = bssp->descr;
  } else return;
  while (sap != NULL) {
    if (sap->type == 1) {
      for (sfp = (SeqFeatPtr) sap->data; sfp != NULL; sfp = sfp->next) {
        if (sfp->data.choice == SEQFEAT_PUB) {
          pdp = (PubdescPtr) sfp->data.value.ptrvalue;
          DoOnePub (pdp);
        }
      }
    }
    sap = sap->next;
  }
  while (sdp != NULL) {
    if (sdp->choice == Seq_descr_pub) {
      pdp = (PubdescPtr) sdp->data.ptrvalue;
      DoOnePub (pdp);
    }
    sdp = sdp->next;
  }
}

extern void LookupAllPubs (IteM i);
extern void LookupAllPubs (IteM i)

{
  BaseFormPtr  bfp;
  MonitorPtr   mon = NULL;
  SeqEntryPtr  sep;
  ErrSev       sev;


  if (! useMedarch) return;
#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  sev = ErrSetMessageLevel (SEV_FATAL);
  WatchCursor ();
  mon = MonitorStrNewEx ("Processing Publications", 40, FALSE);
  MonitorStrValue (mon, "Connecting to MedArch");
  Update ();
  if (! MedArchInit ()) {
    MonitorFree (mon);
    ArrowCursor ();
    Update ();
    Message (MSG_POST, "Unable to connect to MedArch");
    return;
  }
  SeqEntryExplore (sep, NULL, DoLookupPub);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  MonitorStrValue (mon, "Closing MedArch");
  Update ();
  MedArchFini ();
  MonitorFree (mon);
  ArrowCursor ();
  Update ();
  ErrSetMessageLevel (sev);
  ErrClear ();
  ErrShow ();
}

static void LookupPublications (SeqAnnotPtr sap)

{
  MonitorPtr  mon = NULL;
  PubdescPtr  pdp;
  SeqFeatPtr  sfp;
  ValNodePtr  tmp;
  Int4        uid;
  Boolean     usingMedarch = FALSE;
  ValNodePtr  vnp;

  if (! useMedarch) return;
  if (sap == NULL || sap->type != 1) return;
  for (sfp = (SeqFeatPtr) sap->data; sfp != NULL; sfp = sfp->next) {
    if (sfp->data.choice == SEQFEAT_PUB) {
      pdp = (PubdescPtr) sfp->data.value.ptrvalue;
      if (pdp != NULL) {
        vnp = pdp->pub;
        if (vnp != NULL && vnp->next == NULL) {
          if (vnp->choice == PUB_Muid || vnp->choice == PUB_PMid) {
            if (! usingMedarch) {
              WatchCursor ();
              mon = MonitorStrNewEx ("Processing Publications", 40, FALSE);
              MonitorStrValue (mon, "Connecting to MedArch");
              Update ();
              if (MedArchInit ()) {
                usingMedarch = TRUE;
              } else {
                MonitorFree (mon);
                ArrowCursor ();
                Update ();
                Message (MSG_POST, "Unable to connect to MedArch");
                return;
              }
            }
          }
          tmp = NULL;
          if (vnp->choice == PUB_Muid) {
            uid = vnp->data.intvalue;
            tmp = MedArchGetPub (uid);
          } else if (vnp->choice == PUB_PMid) {
            uid = vnp->data.intvalue;
            tmp = MedArchGetPubPmId (uid);
          }
          if (tmp != NULL) {
            MedlineToISO (tmp);
            tmp->next = vnp;
            pdp->pub = tmp;
          }
        }
      }
    }
  }
  if (usingMedarch) {
    MonitorStrValue (mon, "Closing MedArch");
    Update ();
    MedArchFini ();
    MonitorFree (mon);
    ArrowCursor ();
    Update ();
  }
}

static void PromotePubs (SeqFeatPtr first, BioseqPtr bsp, Uint2 entityID)

{
  MsgAnswer    ans;
  Boolean      asked = FALSE;
  PubdescPtr   pdp;
  SeqDescrPtr  sdp;
  SeqEntryPtr  sep;
  SeqFeatPtr   sfp;
  ValNode      vn;

  MemSet ((Pointer) &vn, 0, sizeof (ValNode));
  vn.choice = SEQLOC_WHOLE;
  vn.data.ptrvalue = (Pointer) SeqIdFindBest (bsp->id, 0);
  vn.next = NULL;

  for (sfp = first; sfp != NULL; sfp = sfp->next) {
    if (sfp->data.choice == SEQFEAT_PUB) {
      if (SeqLocCompare (sfp->location, &vn) == SLC_A_EQ_B) {
        if (! asked) {
          ans = Message (MSG_YN, "Do you wish to convert full-length publication features to descriptors?");
          if (ans == ANS_NO) return;
          asked = TRUE;
        }
      }
    }
  }

  sep = GetBestTopParentForData (entityID, bsp);
  for (sfp = first; sfp != NULL; sfp = sfp->next) {
    if (sfp->data.choice == SEQFEAT_PUB) {
      if (SeqLocCompare (sfp->location, &vn) == SLC_A_EQ_B) {
        sfp->idx.deleteme = TRUE;
        sfp->data.choice = SEQFEAT_COMMENT;
        pdp = (PubdescPtr) sfp->data.value.ptrvalue;
        sfp->data.value.ptrvalue = NULL;
        sdp = CreateNewDescriptor (sep, Seq_descr_pub);
        if (sdp != NULL) {
          sdp->data.ptrvalue = (Pointer) pdp;
        }
      }
    }
  }

  DeleteMarkedObjects (entityID, 0, NULL);
}

extern Uint2 SmartAttachSeqAnnotToSeqEntry (Uint2 entityID, SeqAnnotPtr sap)

{
  BioseqPtr      bsp;
  Int2           genCode;
  SeqEntryPtr    oldscope;
  OMProcControl  ompc;
  SeqEntryPtr    sep;
  SeqFeatPtr     sfp = NULL;

  if (sap == NULL) return entityID;
  bsp = GetBioseqReferencedByAnnot (sap, entityID);
  if (bsp == NULL) {
    oldscope = SeqEntrySetScope (NULL);
    if (oldscope != NULL) {
      bsp = GetBioseqReferencedByAnnot (sap, entityID);
      SeqEntrySetScope (oldscope);
    }
  }
  if (bsp != NULL) {
    sep = SeqMgrGetSeqEntryForData (bsp);
    entityID = ObjMgrGetEntityIDForChoice (sep);
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      sep = GetBestTopParentForData (entityID, bsp);
      genCode = SeqEntryToGeneticCode (sep, NULL, NULL, 0);
      SetEmptyGeneticCodes (sap, genCode);
      LookupPublications (sap);
    } else if (sap->type == 2) {
      TakeTop10Alignments (sap);
    }
    MemSet ((Pointer) &ompc, 0, sizeof (OMProcControl));
    ompc.input_entityID = entityID;
    ompc.input_itemID = GetItemIDGivenPointer (entityID, OBJ_BIOSEQ, (Pointer) bsp);
    ompc.input_itemtype = OBJ_BIOSEQ;
    ompc.output_itemtype = OBJ_SEQANNOT;
    ompc.output_data = (Pointer) sap;
    if (! AttachDataForProc (&ompc, FALSE)) {
      Message (MSG_ERROR, "SmartAttachSeqAnnotToSeqEntry failed");
    } else if (sfp != NULL) {
      PromoteXrefs (sfp, bsp, entityID);
      PromotePubs (sfp, bsp, entityID);
    }
  } else {
    Message (MSG_ERROR, "Feature table identifiers do not match record");
  }
  return entityID;
}

typedef struct removeformdata {
  FEATURE_FORM_BLOCK

  Boolean        is_feature;
  LisT           objlist;
  TexT           findthis;
  TexT           fromTxt;
  TexT           toTxt;
  Uint2          itemtype;
  Uint2          subtype;
  ValNodePtr     head;
  Boolean        stringfound;
  Char           findStr [128];
  Boolean        take_action_when_string_present;
  GrouP          string_constraint_type;
  ButtoN         case_insensitive;
  Int4           from;
  Int4           to;
  ValNodePtr     bsplist;
  ValNodePtr     bssplist;
} RemoveFormData, PNTR RemoveFormPtr;

static Boolean ObjectInRange (SeqFeatPtr sfp, Int4 from, Int4 to)

{
  SeqMgrFeatContext  context;

  if (sfp == NULL || from < 0 || to < 0) return TRUE;
  if (SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, 0, 0, sfp, &context) == sfp) {
    if (context.left > to) return FALSE;
    if (context.right < from) return FALSE;
  }
  return TRUE;
}


static void RemoveFeatureCallback (SeqFeatPtr sfp, Pointer userdata)
{
  RemoveFormPtr rfp;
  SeqIdPtr      sip;
  BioseqPtr     productbsp, productcdna;
  BioseqSetPtr  productnps;
  
  if (sfp == NULL || userdata == NULL) return;

  rfp = (RemoveFormPtr) userdata;
  if (rfp == NULL) return;
  if (sfp->idx.subtype == rfp->subtype ||
      (rfp->subtype == FEATDEF_IMP && IsRealImpFeat (sfp->idx.subtype)) ||
      rfp->subtype == ALL_FEATURES) 
  {
    if ((rfp->from == -1 && rfp->to == -1) || ObjectInRange (sfp, rfp->from, rfp->to)) 
    {
      if (sfp->data.choice == SEQFEAT_CDREGION) 
      {
        if (sfp->product != NULL) 
        {
          sip = SeqLocId (sfp->product);
          if (sip != NULL) 
          {
            productbsp = BioseqFind (sip);
            if (productbsp != NULL) 
            {
              ValNodeAddPointer (&(rfp->bsplist), 0, (Pointer) productbsp);
            }
          }
        }
      } 
      else if (sfp->data.choice == SEQFEAT_RNA) 
      {
        if (sfp->product != NULL) 
        {
          sip = SeqLocId (sfp->product);
          if (sip != NULL) 
          {
            productcdna = BioseqFind (sip);
            if (productcdna != NULL && productcdna->idx.parenttype == OBJ_BIOSEQSET) 
            {
              productnps = (BioseqSetPtr) productcdna->idx.parentptr;
              if (productnps != NULL && productnps->_class == BioseqseqSet_class_nuc_prot) 
              {
                ValNodeAddPointer (&(rfp->bssplist), 0, (Pointer) productnps);
              }
            }
          }
        }
      }
      sfp->idx.deleteme = TRUE;  
    }
  }
}

static void RemoveFeatures (SeqEntryPtr sep, RemoveFormPtr rfp)
{
  FeaturesWithTextData fd;
  Char           str [32];
  Int4           swap;
  long int       val;
  
  GetTitle (rfp->findthis, rfp->findStr, sizeof (rfp->findStr) - 1);
  fd.search_text = rfp->findStr;
  fd.no_text = StringHasNoText (rfp->findStr);
  fd.seqFeatChoice = 0;
  fd.featDefChoice = 0;
  fd.case_insensitive = GetStatus (rfp->case_insensitive);
  fd.whole_word = FALSE;
  fd.act_when_string_not_present = ! rfp->take_action_when_string_present;
  fd.userdata = rfp;
  fd.callback = RemoveFeatureCallback;
  GetTitle (rfp->fromTxt, str, sizeof (str) - 1);
  if ((! StringHasNoText (str)) && sscanf (str, "%ld", &val) == 1 && val >= 0) {
    rfp->from = (Int4) val;
  } else {
    rfp->from = -1;
  }
  GetTitle (rfp->toTxt, str, sizeof (str) - 1);
  if ((! StringHasNoText (str)) && sscanf (str, "%ld", &val) == 1 && val >= 0) {
    rfp->to = (Int4) val;
  } else {
    rfp->to = -1;
  }
  if (rfp->from > rfp->to) {
    swap = rfp->from;
    rfp->from = rfp->to;
    rfp->to = swap;
  }
  OperateOnSeqEntryFeaturesWithText (sep, &fd);
  DeleteMarkedObjects (rfp->input_entityID, OBJ_SEQENTRY, (Pointer) sep);
}

static void RemoveDescriptorCallback (SeqDescrPtr sdp, Pointer userdata)
{
  ObjValNodePtr ovp;
  RemoveFormPtr rfp;
  
  if (sdp == NULL || userdata == NULL || sdp->extended == 0) return;
  rfp = (RemoveFormPtr) userdata;
  if (rfp == NULL) return;
  
  ovp = (ObjValNodePtr) sdp;
    
  if (sdp->choice == rfp->subtype) 
  {
    ovp->idx.deleteme = TRUE;	
  }
}

static void RemoveDescriptors (SeqEntryPtr sep, RemoveFormPtr rfp)
{
  DescriptorsWithTextData dd;
  Char                    str [32];
  Int4                    swap;
  long int                val;
  
  GetTitle (rfp->findthis, rfp->findStr, sizeof (rfp->findStr) - 1);
  dd.search_text = rfp->findStr;
  dd.no_text = StringHasNoText (rfp->findStr);
  dd.case_insensitive = GetStatus (rfp->case_insensitive);
  dd.whole_word = FALSE;
  dd.act_when_string_not_present = ! rfp->take_action_when_string_present;
  dd.userdata = rfp;
  dd.callback = RemoveDescriptorCallback;
  GetTitle (rfp->fromTxt, str, sizeof (str) - 1);
  if ((! StringHasNoText (str)) && sscanf (str, "%ld", &val) == 1 && val >= 0) {
    rfp->from = (Int4) val;
  } else {
    rfp->from = -1;
  }
  GetTitle (rfp->toTxt, str, sizeof (str) - 1);
  if ((! StringHasNoText (str)) && sscanf (str, "%ld", &val) == 1 && val >= 0) {
    rfp->to = (Int4) val;
  } else {
    rfp->to = -1;
  }
  if (rfp->from > rfp->to) {
    swap = rfp->from;
    rfp->from = rfp->to;
    rfp->to = swap;
  }
  OperateOnSeqEntryDescriptorsWithText (sep, &dd);
  DeleteMarkedObjects (rfp->input_entityID, OBJ_SEQENTRY, (Pointer) sep);
}

static void DoRemoveAsnObject (ButtoN b)

{
  MsgAnswer      ans;
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  Uint2          itemID;
  OMProcControl  ompc;
  RemoveFormPtr  rfp;
  SeqEntryPtr    sep;
  ValNodePtr     tmp;
  Int2           val;
  ValNodePtr     vnp;
  Boolean        removed_some_features;

  rfp = GetObjectExtra (b);
  if (rfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (rfp->input_entityID);
  if (sep == NULL) return;
  Hide (rfp->form);
  WatchCursor ();
  Update ();
  if (rfp->is_feature) {
    rfp->itemtype = OBJ_SEQFEAT;
  } else {
    rfp->itemtype = OBJ_SEQDESC;
  }

  if (rfp->itemtype == 0) return;

  removed_some_features = FALSE;
  
  if (GetValue (rfp->string_constraint_type) == 1)
  {
  	rfp->take_action_when_string_present = TRUE;
  }
  else
  {
    rfp->take_action_when_string_present = FALSE;
  }

  val = 1;
  for (vnp = rfp->head; vnp != NULL; vnp = vnp->next)
  {
    if (GetItemStatus (rfp->objlist, val))
    {
      rfp->subtype = vnp->choice;
      if (rfp->subtype != 0) {
        if (rfp->is_feature) {
          RemoveFeatures (sep, rfp);
          removed_some_features = TRUE;
        } else {
          RemoveDescriptors (sep, rfp);
        }
      }
    }
    val ++;
  }

  if (removed_some_features) {
    if (rfp->bsplist != NULL) {
      ans = Message (MSG_YN, "Remove protein products?");
      if (ans == ANS_YES) {
        for (tmp = rfp->bsplist; tmp != NULL; tmp = tmp->next) {
          bsp = (BioseqPtr) tmp->data.ptrvalue;
          itemID = GetItemIDGivenPointer (rfp->input_entityID, OBJ_BIOSEQ, (Pointer) bsp);
          if (itemID > 0) {
            MemSet ((Pointer) (&ompc), 0, sizeof (OMProcControl));
            ompc.do_not_reload_from_cache = TRUE;
            ompc.input_entityID = rfp->input_entityID;
            ompc.input_itemID = itemID;
            ompc.input_itemtype = OBJ_BIOSEQ;
            if (! DetachDataForProc (&ompc, FALSE)) {
              Message (MSG_POSTERR, "DetachDataForProc failed");
            }
            SeqMgrDeleteFromBioseqIndex (bsp);
          }
        }
        ans = Message (MSG_YN, "Renormalize Nuc-Prot sets?");
        if (ans == ANS_YES)
        {
          RemoveOrphanProteins (rfp->input_entityID, sep);
          RenormalizeNucProtSets (sep, TRUE);   	
        }
      }
    }
    if (rfp->bssplist != NULL) {
      ans = Message (MSG_YN, "Remove cDNA nuc-prot products?");
      if (ans == ANS_YES) {
        for (tmp = rfp->bssplist; tmp != NULL; tmp = tmp->next) {
          bssp = (BioseqSetPtr) tmp->data.ptrvalue;
          itemID = GetItemIDGivenPointer (rfp->input_entityID, OBJ_BIOSEQSET, (Pointer) bssp);
          if (itemID > 0) {
            MemSet ((Pointer) (&ompc), 0, sizeof (OMProcControl));
            ompc.do_not_reload_from_cache = TRUE;
            ompc.input_entityID = rfp->input_entityID;
            ompc.input_itemID = itemID;
            ompc.input_itemtype = OBJ_BIOSEQSET;
            if (! DetachDataForProc (&ompc, FALSE)) {
              Message (MSG_POSTERR, "DetachDataForProc failed");
            }
          }
        }
      }
    }
  }
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (rfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, rfp->input_entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  Remove (rfp->form);
}

static void RemoveDefLinesCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  ValNodePtr     nextsdp;
  Pointer PNTR   prevsdp;
  ValNodePtr     sdp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sdp = bsp->descr;
    prevsdp = (Pointer PNTR) &(bsp->descr);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sdp = bssp->descr;
    prevsdp = (Pointer PNTR) &(bssp->descr);
  } else return;

  while (sdp != NULL) {
    nextsdp = sdp->next;
    if (sdp->choice == Seq_descr_title)
    {
      *(prevsdp) = sdp->next;
      sdp->next = NULL;
      SeqDescFree (sdp);
    } else {
      prevsdp = (Pointer PNTR) &(sdp->next);
    }
    sdp = nextsdp;
  }
}

extern void RemoveDefLinesToolBtn (ButtoN b)
{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  
  WatchCursor ();
  Update ();

  SeqEntryExplore (sep, NULL, RemoveDefLinesCallback);

  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  CommonApplyToAllProc (bfp, ADD_TITLE);
}

int LIBCALLBACK SortByVnpChoice (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr   vnp1;
  ValNodePtr   vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      if (vnp1->choice > vnp2->choice) {
        return 1;
      } else if (vnp1->choice < vnp2->choice) {
        return -1;
      } else {
        return 0;
      }
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}

static void RemoveMessageProc (ForM f, Int2 mssg)

{
  RemoveFormPtr  rfp;

  rfp = (RemoveFormPtr) GetObjectExtra (f);
  if (rfp != NULL) {
    if (rfp->appmessage != NULL) {
      rfp->appmessage (f, mssg);
    }
  }
}

static void CleanupRemovePage (GraphiC g, VoidPtr data)

{
  RemoveFormPtr  rfp;

  rfp = (RemoveFormPtr) data;
  if (rfp != NULL) {
    ValNodeFreeData (rfp->head);
    ValNodeFree (rfp->bsplist);
    ValNodeFree (rfp->bssplist);
  }
  StdCleanupFormProc (g, data);
}

static CharPtr descNames [] = {
  " ", " ", " ", " ", "Name",
  "Title", " ", "Comment", "Numbering",
  "MapLoc", "PIR", "GenBank", "Publication",
  "Region", "User", "SWISS-PROT", "dbXREF",
  "EMBL", "Create Date", "Update Date", "PRF",
  "PDB", "Heterogen", "BioSource", "MolInfo", NULL
};

/*
#ifdef INTERNAL_NCBI_SEQUIN
#define LISTHEIGHT 16
#else
#define LISTHEIGHT 8
#endif
*/

CharPtr MostUsedDescriptorList[] = { "Title" };

static Boolean isMostUsedDescriptor (CharPtr descname)
{
  Int2 i;

  if (descname == NULL) return FALSE;

  for (i=0; i < sizeof (MostUsedDescriptorList) / sizeof (CharPtr); i++)
  {
    if (StringCmp (descname, MostUsedDescriptorList[i]) == 0)
      return TRUE;
  }
  return FALSE;
}

static int LIBCALLBACK SortMostUsedDescriptorsFirst (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr   vnp1;
  ValNodePtr   vnp2;
  CharPtr      str1;
  CharPtr      str2;
  Boolean      str1_is_most_used;
  Boolean      str2_is_most_used;

  /* Check parameters */

  if ((NULL == ptr1) || (NULL == ptr2))
    return 0;

  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if ((NULL == vnp1) || (NULL == vnp2))
    return 0;

  str1 = (CharPtr) vnp1->data.ptrvalue;
  str2 = (CharPtr) vnp2->data.ptrvalue;
  if ((NULL == str1) || (NULL == str2))
    return 0;

  str1_is_most_used = isMostUsedDescriptor (str1);
  str2_is_most_used = isMostUsedDescriptor (str2);

  if ((str1_is_most_used && str2_is_most_used)
    || (!str1_is_most_used && !str2_is_most_used))
  {
    return SortVnpByString (ptr1, ptr2);
  }
  else if (str1_is_most_used)
  {
    return -1;
  }
  else
  {
    return 1;
  }
}

extern ValNodePtr BuildDescriptorValNodeList (void)
{
  Int4 j;
  ValNodePtr vnp;
  ValNodePtr head = NULL;
  
  for (j = 1; descNames [j] != NULL; j++) {
    if (StringHasNoText (descNames [j])) continue;
    vnp = ValNodeNew (head);
    if (head == NULL) {
      head = vnp;
    }
    if (vnp != NULL) {
      vnp->choice = j;
      vnp->data.ptrvalue = StringSave (descNames [j]);
    }
  }
  head = SortValNode (head, SortMostUsedDescriptorsFirst);
  return head;
}

static void RemoveAsnObject (IteM i, Boolean feature)

{
  BaseFormPtr        bfp;
  ButtoN             b;
  GrouP              c;
  GrouP              g;
  GrouP              h;
  ValNodePtr         head;
  GrouP              k;
  Int2               listHeight;
  GrouP              m;
  RemoveFormPtr      rfp;
  SeqEntryPtr        sep;
  StdEditorProcsPtr  sepp;
  CharPtr            title;
  ValNodePtr         vnp;
  WindoW             w;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  rfp = (RemoveFormPtr) MemNew (sizeof (RemoveFormData));
  if (rfp == NULL) return;
  if (feature) {
    title = "Feature Removal";
  } else {
    title = "Descriptor Removal";
  }
  w = FixedWindow (-50, -33, -10, -10, title, StdCloseWindowProc);
  SetObjectExtra (w, rfp, CleanupRemovePage);
  rfp->form = (ForM) w;
  rfp->formmessage = RemoveMessageProc;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    rfp->appmessage = sepp->handleMessages;
  }

  rfp->input_entityID = bfp->input_entityID;
  rfp->input_itemID = bfp->input_itemID;
  rfp->input_itemtype = bfp->input_itemtype;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, 0, 2, NULL);
  rfp->is_feature = feature;
  if (feature) {
    StaticPrompt (g, "Feature", 0, 0, programFont, 'c');
  } else {
    StaticPrompt (g, "Descriptor", 0, 0, programFont, 'c');
  }
  if (indexerVersion) {
    listHeight = 16;
  } else {
    listHeight = 8;
  }
  rfp->objlist = MultiList (g, 16, listHeight, NULL);
  head = NULL;
  if (feature) {
    head = BuildFeatureValNodeList (TRUE, "All", ALL_FEATURES, TRUE, FALSE);
  } else {
    head = BuildDescriptorValNodeList();
  }
  if (head != NULL) {

    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      ListItem (rfp->objlist, (CharPtr) vnp->data.ptrvalue);
    }
  }
  rfp->head = head;
  rfp->bsplist = NULL;
  rfp->bssplist = NULL;

  k = NormalGroup (h, 0, 3, "Optional string constraint", NULL, NULL);
  rfp->string_constraint_type = HiddenGroup (k, 0, 2, NULL);
  RadioButton (rfp->string_constraint_type, "Remove when text is present");
  RadioButton (rfp->string_constraint_type, "Remove when text is not present");
  SetValue (rfp->string_constraint_type, 1);
  rfp->findthis = DialogText (k, "", 14, NULL);
  rfp->case_insensitive = CheckBox (k, "Case Insensitive", NULL);

  m = NULL;
  if (feature) {
    m = HiddenGroup (h, 4, 0, NULL);
    StaticPrompt (m, "From", 0, dialogTextHeight, programFont, 'l');
    rfp->fromTxt = DialogText (m, "", 6, NULL);
    StaticPrompt (m, "To", 0, dialogTextHeight, programFont, 'l');
    rfp->toTxt = DialogText (m, "", 6, NULL);
  }

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoRemoveAsnObject);
  SetObjectExtra (b, rfp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) k, (HANDLE) c, (HANDLE) m, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
}

extern void RemoveDescriptor (IteM i)

{
  RemoveAsnObject (i, FALSE);
}

#define SLCT_FEAT    1
#define SLCT_DESC    2
#define SLCT_BIOSEQ  3
#define SLCT_PUB     4

typedef struct selectformdata {
  FEATURE_FORM_BLOCK

  Int2           type;
  LisT           objlist;
  TexT           findthis;
  Uint2          itemtype;
  Uint2          subtype;
  ObjMgrPtr      omp;
  ObjMgrTypePtr  omtp;
  ValNodePtr     head;
  Boolean        stringfound;
  Char           findStr [128];
  ButtoN         when_string_not_present;
  ButtoN         case_insensitive;
} SelectFormData, PNTR SelectFormPtr;

static void FeatureSelectCallback (SeqFeatPtr sfp, Pointer userdata)
{
  Uint1Ptr subtype;
  if (sfp == NULL) return;
  
  if (userdata != NULL)
  {
  	subtype = (Uint1Ptr) userdata;
  	if (*subtype != sfp->idx.subtype) return;
  }
  ObjMgrAlsoSelect (sfp->idx.entityID, sfp->idx.itemID, OBJ_SEQFEAT, 0, NULL);
}

static void DescriptorSelectCallback (SeqDescrPtr sdp, Pointer userdata)
{
  ObjValNodePtr ovp;
  Uint1Ptr subtype;
  
  if (sdp == NULL || sdp->extended == 0) return;
  
  ovp = (ObjValNodePtr) sdp;
  if (userdata != NULL)
  {
  	subtype = (Uint1Ptr) userdata;
  	if (*subtype != ovp->idx.subtype) return;
  }
    
  ObjMgrAlsoSelect (ovp->idx.entityID, ovp->idx.itemID, OBJ_SEQDESC, 0, NULL);
}

static void BioseqSelectCallback (BioseqPtr bsp, Pointer userdata)
{
  Uint1Ptr subtype;
  
  if (bsp == NULL) return;
  
  if (userdata != NULL)
  {
  	subtype = (Uint1Ptr) userdata;
  	if (*subtype != bsp->idx.subtype) return;
  }
    
  ObjMgrAlsoSelect (bsp->idx.entityID, bsp->idx.itemID, bsp->idx.itemtype, 0, NULL);
}


static void DoSelectAsnObject (ButtoN b)

{
  SelectFormPtr           selfp;
  SeqEntryPtr             sep;
  Int2                    val;
  ValNodePtr              vnp;
  FeaturesWithTextData    fd;
  DescriptorsWithTextData dd;
  Uint1                   bioseq_choice = Seq_repr_raw;
  Uint1                   pub_choice = FEATDEF_PUB;

  selfp = GetObjectExtra (b);
  if (selfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (selfp->input_entityID);
  if (sep == NULL) return;
  Hide (selfp->form);
  
  vnp = NULL;
  if (selfp->type == SLCT_FEAT || selfp->type == SLCT_DESC)
  {
    val = GetValue (selfp->objlist);
    if (val > 0) {
      vnp = selfp->head;
      while (vnp != NULL && val > 1) {
        val--;
        vnp = vnp->next;
      }
    }
  }

  switch (selfp->type) {
    case SLCT_FEAT :
      GetTitle (selfp->findthis, selfp->findStr, sizeof (selfp->findStr) - 1);
      fd.search_text = selfp->findStr;
      fd.no_text = StringHasNoText (selfp->findStr);
      fd.seqFeatChoice = 0;
      fd.featDefChoice = 0;
      fd.case_insensitive = GetStatus (selfp->case_insensitive);
      fd.whole_word = FALSE;
      fd.act_when_string_not_present = GetStatus (selfp->when_string_not_present);
      fd.callback = FeatureSelectCallback;
      if (vnp == NULL)
      {
      	fd.userdata = NULL;
      }
      else
      {
        fd.userdata = (Pointer) &(vnp->choice);
      }
      OperateOnSeqEntryFeaturesWithText (sep, &fd); 	
      break;
    case SLCT_DESC :
      GetTitle (selfp->findthis, selfp->findStr, sizeof (selfp->findStr) - 1);
      dd.search_text = selfp->findStr;
      dd.no_text = StringHasNoText (selfp->findStr);
      dd.case_insensitive = GetStatus (selfp->case_insensitive);
      dd.whole_word = FALSE;
      dd.act_when_string_not_present = GetStatus (selfp->when_string_not_present);
      dd.callback = DescriptorSelectCallback;
      if (vnp == NULL)
      {
      	dd.userdata = NULL;
      }
      else
      {
        dd.userdata = (Pointer) &(vnp->choice);
      }
      OperateOnSeqEntryDescriptorsWithText (sep, &dd); 	
      break;
    case SLCT_BIOSEQ :
 	  VisitBioseqsInSep (sep, (Pointer) &bioseq_choice, BioseqSelectCallback);	
      break;
  	case SLCT_PUB:
      GetTitle (selfp->findthis, selfp->findStr, sizeof (selfp->findStr) - 1);
      fd.search_text = selfp->findStr;
      fd.no_text = StringHasNoText (selfp->findStr);
      fd.seqFeatChoice = 0;
      fd.featDefChoice = 0;
      fd.case_insensitive = GetStatus (selfp->case_insensitive);
      fd.whole_word = FALSE;
      fd.act_when_string_not_present = GetStatus (selfp->when_string_not_present);
      fd.callback = FeatureSelectCallback;
      fd.userdata = (Pointer) &pub_choice;
      OperateOnSeqEntryFeaturesWithText (sep, &fd); 
      dd.search_text = fd.search_text;
      dd.no_text = fd.no_text;
      dd.case_insensitive = fd.case_insensitive;
      dd.whole_word = fd.whole_word;
      dd.act_when_string_not_present = fd.act_when_string_not_present;
      dd.callback = DescriptorSelectCallback;
      dd.userdata = fd.userdata;	
      OperateOnSeqEntryDescriptorsWithText (sep, &dd); 	
      break;  	
    default :
      Remove (selfp->form);
      Update ();
      return;
  }
  WatchCursor ();
  Update ();

  ArrowCursor ();
  Update ();
  /* ObjMgrSendMsg (OM_MSG_UPDATE, selfp->input_entityID, 0, 0); */
  Remove (selfp->form);
}

static void SelectMessageProc (ForM f, Int2 mssg)

{
  SelectFormPtr  selfp;

  selfp = (SelectFormPtr) GetObjectExtra (f);
  if (selfp != NULL) {
    if (selfp->appmessage != NULL) {
      selfp->appmessage (f, mssg);
    }
  }
}

static void CleanupSelectPage (GraphiC g, VoidPtr data)

{
  SelectFormPtr  selfp;

  selfp = (SelectFormPtr) data;
  if (selfp != NULL) {
    ValNodeFreeData (selfp->head);
  }
  StdCleanupFormProc (g, data);
}

static void SelectAsnObject (IteM i, Int2 type)

{
  BaseFormPtr        bfp;
  ButtoN             b;
  GrouP              c;
  GrouP              g;
  GrouP              h;
  GrouP              k, m;
  ValNodePtr         head;
  Int2               j;
  Int2               listHeight;
  SelectFormPtr      selfp;
  SeqEntryPtr        sep;
  StdEditorProcsPtr  sepp;
  CharPtr            title;
  ValNodePtr         vnp;
  WindoW             w;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  selfp = (SelectFormPtr) MemNew (sizeof (SelectFormData));
  if (selfp == NULL) return;
  switch (type) {
    case SLCT_FEAT :
      title = "Feature Selection";
      break;
    case SLCT_DESC :
      title = "Descriptor Selection";
      break;
    case SLCT_BIOSEQ :
      title = "Sequence Selection";
      break;
  	case SLCT_PUB:
  	  title = "Publication Selection";
  	  break;
    default :
      title = "? Selection";
      break;
  }
  w = FixedWindow (-50, -33, -10, -10, title, StdCloseWindowProc);
  SetObjectExtra (w, selfp, CleanupSelectPage);
  selfp->form = (ForM) w;
  selfp->formmessage = SelectMessageProc;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    selfp->appmessage = sepp->handleMessages;
  }

  selfp->input_entityID = bfp->input_entityID;
  selfp->input_itemID = bfp->input_itemID;
  selfp->input_itemtype = bfp->input_itemtype;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, 0, 2, NULL);
  selfp->type = type;
  switch (type) {
    case SLCT_FEAT :
      StaticPrompt (g, "Feature", 0, 0, programFont, 'c');
      break;
    case SLCT_DESC :
      StaticPrompt (g, "Descriptor", 0, 0, programFont, 'c');
      break;
    case SLCT_BIOSEQ :
      StaticPrompt (g, "Sequence", 0, 0, programFont, 'c');
      break;
  	case SLCT_PUB:
  	  StaticPrompt (g, "Publication", 0, 0, programFont, 'c');
  	  break;
    default :
      break;
  }
  if (indexerVersion) {
    listHeight = 16;
  } else {
    listHeight = 8;
  }
  
  if (type != SLCT_PUB)
  {
    selfp->objlist = SingleList (g, 16, listHeight, NULL);
  }
  head = NULL;
  if (type == SLCT_FEAT) {
    head = BuildFeatureValNodeList (TRUE, NULL, ALL_FEATURES, TRUE, FALSE);
  } else if (type == SLCT_DESC) {
    for (j = 1; descNames [j] != NULL; j++) {
      if (StringHasNoText (descNames [j])) continue;
      vnp = ValNodeNew (head);
      if (head == NULL) {
        head = vnp;
      }
      if (vnp != NULL) {
        vnp->choice = j;
        vnp->data.ptrvalue = StringSave (descNames [j]);
      }
    }
  }
  if (head != NULL) {
    if (type != SLCT_FEAT) {
      head = SortValNode (head, SortByVnpChoice);
    }
    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      ListItem (selfp->objlist, (CharPtr) vnp->data.ptrvalue);
    }
  }
  selfp->head = head;

  if (selfp->type == SLCT_FEAT || selfp->type == SLCT_DESC || selfp->type == SLCT_PUB)
  {
    k = HiddenGroup (h, 0, 3, NULL);
    StaticPrompt (k, "Optional string constraint", 0, dialogTextHeight, programFont, 'c');
    selfp->findthis = DialogText (k, "", 14, NULL);
    m = HiddenGroup (k, 2, 0, NULL);
    selfp->case_insensitive = CheckBox (m, "Case Insensitive", NULL);
    selfp->when_string_not_present = CheckBox (m, "When String Not Present", NULL);
  }

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoSelectAsnObject);
  SetObjectExtra (b, selfp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  if (selfp->type == SLCT_FEAT || selfp->type == SLCT_DESC)
  {
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) k, (HANDLE) c, NULL);
  }
  else
  {
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
  }
  RealizeWindow (w);
  if (type == SLCT_BIOSEQ) {
    DoSelectAsnObject (b);
    Update ();
    return;
  }
  Show (w);
  Update ();
}

extern void SelectDescriptor (IteM i)

{
  SelectAsnObject (i, SLCT_DESC);
}

extern void SelectBioseq (IteM i)

{
  SelectAsnObject (i, SLCT_BIOSEQ);
}

extern void SelectPubs (IteM i)

{
  SelectAsnObject (i, SLCT_PUB);
}

typedef struct fuseformdata {
  FEATURE_FORM_BLOCK

  LisT           objlist;
  Uint2          subtype;
  ValNodePtr     head;
} FuseFormData, PNTR FuseFormPtr;

static SeqLocPtr FuseTwoLocations (Uint2 entityID, SeqLocPtr slp1, SeqLocPtr slp2)
{
  Boolean           slp1_partial5, slp1_partial3;
  Boolean           slp2_partial5, slp2_partial3;
  Int4              slp1_start_pos, slp1_stop_pos;
  Int4              slp2_start_pos, slp2_stop_pos;
  SeqLocPtr         slp_result = NULL, slp2_copy, slp1_copy, slp_list;
  BioseqPtr         bsp1, bsp2;

  if (slp1 == NULL || slp2 == NULL)
  {
    return NULL;
  }
  
  bsp1 = GetBioseqGivenSeqLoc (slp1, entityID);
  bsp2 = GetBioseqGivenSeqLoc (slp2, entityID);
  
  /* preserve partialness of ends */
  CheckSeqLocForPartial (slp1, &slp1_partial5, &slp1_partial3);
  slp1_start_pos = SeqLocStart (slp1);
  slp1_stop_pos = SeqLocStop (slp1);
  CheckSeqLocForPartial (slp2, &slp2_partial5, &slp2_partial3);
  slp2_start_pos = SeqLocStart (slp2);
  slp2_stop_pos = SeqLocStop (slp2);
  if (slp1_start_pos > slp2_start_pos)
  {
    slp1_partial5 = slp2_partial5;
  }
  if (slp1_stop_pos < slp2_stop_pos)
  {
    slp1_partial3 = slp2_partial3;
  }

  if (bsp1 == bsp2)
  {
    slp_result = SeqLocMerge (bsp1, slp2, slp1, FALSE, TRUE, FALSE);
    SetSeqLocPartial (slp_result, slp1_partial5, slp1_partial3);
  }
  else 
  {
    /* we are dealing with a segmented set, and both locations are segments */
    slp1_copy = (SeqLocPtr) AsnIoMemCopy (slp1, (AsnReadFunc) SeqLocAsnRead,
                                          (AsnWriteFunc) SeqLocAsnWrite);
    slp2_copy = (SeqLocPtr) AsnIoMemCopy (slp2, (AsnReadFunc) SeqLocAsnRead,
                                          (AsnWriteFunc) SeqLocAsnWrite);
    if (slp1_copy != NULL && slp2_copy != NULL)
    {
      /* if the second location is already a mix, don't create a nested
       * mixed location.
       */
      if (slp2_copy->choice == SEQLOC_MIX)
      {
        slp_list = (SeqLocPtr) slp2_copy->data.ptrvalue;
        slp2_copy->data.ptrvalue = NULL;
        slp2_copy = SeqLocFree (slp2_copy);
        slp2_copy = slp_list;
      }
    
      if (slp1_copy->choice == SEQLOC_MIX)
      {
        slp_list = slp1_copy->data.ptrvalue;
        while (slp_list != NULL && slp_list->next != NULL)
        {
          slp_list = slp_list->next;
        }
        
        
        if (slp_list == NULL)
        {
          slp1_copy->data.ptrvalue = slp2_copy;
        }
        else
        {
          slp_list->next = slp2_copy;
        }
        slp_result = slp1_copy;
      }
      else
      {
        slp_result = ValNodeNew (NULL);
        if (slp_result != NULL)
        {
          slp_result->choice = SEQLOC_MIX;
          slp_result->data.ptrvalue = slp1_copy;
          slp1_copy->next = slp2_copy;
        }
      }
    }
  }
  
  return slp_result;
}

static void CombineProductFeatures (BioseqPtr pbsp1, BioseqPtr pbsp2, Int4 old_len)
{
  SeqAnnotPtr sap1, sap2, last_sap1 = NULL;
  SeqFeatPtr  sfp2, sfp_new, last_sfp1 = NULL, main_prot = NULL;
  Boolean     partial5_orig = TRUE, partial3_orig = TRUE;
  Boolean     partial5_new, partial3_new;
  SeqEntryPtr psep;
  ProtRefPtr  prp;
  
  if (pbsp1 == NULL || pbsp2 == NULL)
  {
    return;
  }
  
  sap2 = pbsp2->annot;
  while (sap2 != NULL && sap2->type != 1)
  {
    sap2 = sap2->next;
  }
  if (sap2 == NULL || sap2->data == NULL)
  {
    /* second sequence has no features */
    return;
  }
  
  sap1 = pbsp1->annot;
  while (sap1 != NULL && sap1->type != 1)
  {
    last_sap1 = sap1;
    sap1 = sap1->next;
  }
  if (sap1 == NULL)
  {
    sap1 = SeqAnnotNew();
    if (sap1 == NULL)
    {
      return;
    }
    sap1->type = 1;
    sap1->data = NULL;
    if (last_sap1 == NULL)
    {
      pbsp1->annot = sap1;
    }
    else
    {
      last_sap1->next = sap1;
    }
  }
  
  last_sfp1 = sap1->data;
  if (last_sfp1 != NULL && last_sfp1->idx.subtype == FEATDEF_PROT)
  {
    main_prot = last_sfp1;
    CheckSeqLocForPartial (main_prot->location, &partial5_orig, &partial3_orig);
  }
  while (last_sfp1 != NULL && last_sfp1->next != NULL)
  {
    if (main_prot == NULL && last_sfp1->idx.subtype == FEATDEF_PROT)
    {
      main_prot = last_sfp1;
      CheckSeqLocForPartial (main_prot->location, &partial5_orig, &partial3_orig);
    }
    last_sfp1 = last_sfp1->next;
  }
  if (last_sfp1 != NULL && main_prot == NULL && last_sfp1->idx.subtype == FEATDEF_PROT)
  {
    main_prot = last_sfp1;
    CheckSeqLocForPartial (main_prot->location, &partial5_orig, &partial3_orig);
  }
  
  partial3_new = partial3_orig;
  
  while (sap2 != NULL)
  {
    if (sap2->type == 1)
    {
      for (sfp2 = sap2->data; sfp2 != NULL; sfp2 = sfp2->next)
      {  
        if (sfp2->idx.subtype == FEATDEF_PROT)
        {
          CheckSeqLocForPartial (sfp2->location, &partial5_new, &partial3_new);
          /* do not create additional full-length protein features */
          continue;
        }
        sfp_new = (SeqFeatPtr) AsnIoMemCopy (sfp2, (AsnReadFunc) SeqFeatAsnRead,
                                               (AsnWriteFunc) SeqFeatAsnWrite);
        if (sfp_new != NULL)
        {
          OffsetLocation (sfp_new->location, old_len, pbsp1->id);

          if (last_sfp1 == NULL)
          {
            sap1->data = sfp_new;
          }
          else
          {
            last_sfp1->next = sfp_new;
          }
          last_sfp1 = sfp_new;
        }
      }
    }
    sap2 = sap2->next;
  }

  /* make sure there is one full-length protein feature */
  if (main_prot == NULL)
  {
    psep = SeqMgrGetSeqEntryForData (pbsp1);
    main_prot = CreateNewFeature (psep, NULL, SEQFEAT_PROT, NULL);
    if (main_prot != NULL) 
    {
      prp = ProtRefNew ();
      main_prot->data.value.ptrvalue = (Pointer) prp;
    }
  }
  if (main_prot != NULL)
  {
    if (main_prot->location == NULL || main_prot->location->choice != SEQLOC_INT)
    {
      main_prot->location = SeqLocFree (main_prot->location);
      main_prot->location = SeqLocIntNew (0, pbsp1->length - 1, 
                                          Seq_strand_plus,
                                          SeqIdDup (pbsp1->id));
    }
    SetSeqLocPartial (main_prot->location, partial5_orig, partial3_new);
    main_prot->partial = (Boolean) (partial5_orig || partial3_new);
  }
}

static void FuseTwoProducts (SeqFeatPtr sfp1, SeqFeatPtr sfp2, Uint2 entityID)
{
  BioseqPtr    pbsp1, pbsp2;
  CharPtr      pstr1, pstr2;
  ByteStorePtr byte_store, bs2 = NULL;
  Int4         old_length, added_length;
  SeqIdPtr     sip1;
  
  if (sfp1 == NULL || sfp2 == NULL 
      || sfp1->idx.subtype != FEATDEF_CDS
      || sfp2->idx.subtype != FEATDEF_CDS)
  {
    return;
  }
  pbsp1 = BioseqFindFromSeqLoc (sfp1->product);
  pbsp2 = BioseqFindFromSeqLoc (sfp2->product);
  
  if (pbsp1 == NULL)
  {
    sip1 = SeqLocId (sfp1->product);
    pbsp1 = BioseqFind (sip1);
    if (pbsp1 == NULL)
    {
      RetranslateOneCDS (sfp1, entityID, TRUE, FALSE);
      pbsp1 = BioseqFindFromSeqLoc (sfp1->product);
      if (pbsp1 == NULL)
      {
        sip1 = SeqLocId (sfp1->product);
        pbsp1 = BioseqFind (sip1);
      }
    }
  }
  if (pbsp1 == NULL)
  {
    return;
  }
  pstr1 = BSMerge ((ByteStorePtr)(pbsp1->seq_data), NULL);
  old_length = pbsp1->length;
  
  if (pbsp2 == NULL)
  {
    bs2 = ProteinFromCdRegionEx (sfp2, TRUE, FALSE);
    pstr2 = BSMerge (bs2, NULL);
    added_length = BSLen (bs2);
  }
  else
  {
    pstr2 = BSMerge ((ByteStorePtr)(pbsp2->seq_data), NULL);
    added_length = pbsp2->length;
  }
  
  byte_store = BSNew (old_length + added_length);
  if (byte_store != NULL)
  {
    BSWrite (byte_store, pstr1, StringLen (pstr1));
    BSWrite (byte_store, pstr2, StringLen (pstr2));
    pbsp1->seq_data = BSFree (pbsp1->seq_data);
    pbsp1->seq_data = byte_store;
    pbsp1->length += added_length;
      
    /* now copy features from the second protein to the first */
    CombineProductFeatures (pbsp1, pbsp2, old_length);     
      
    /* remove unused protein */
    if (pbsp2 != NULL)
    {
      pbsp2->idx.deleteme = TRUE;
    }
  }
  bs2 = BSFree (bs2);
}

static void FuseFeatureCallback (BioseqPtr bsp, Pointer userdata)
{
  FuseFormPtr       ffp;
  SeqFeatPtr        first = NULL, sfp = NULL;
  SeqMgrFeatContext context;
  SeqLocPtr         slp;
  
  if (bsp == NULL || userdata == NULL)
  {
    return;
  }
  
  ffp = (FuseFormPtr) userdata;
  
  sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &context);
  while (sfp != NULL)
  {
    if (sfp->idx.subtype == ffp->subtype ||
           (ffp->subtype == FEATDEF_IMP &&
            IsRealImpFeat (sfp->idx.subtype))) 
    {
      if (first == NULL)
      {
        first = sfp;
      }
      else
      {
        slp = FuseTwoLocations (ffp->input_entityID, first->location, sfp->location);
        first->location = SeqLocFree (first->location);
        first->location = slp;
        first->partial = CheckSeqLocForPartial (slp, NULL, NULL);
        sfp->idx.deleteme = TRUE;
        if (sfp->idx.subtype == FEATDEF_CDS)
        {
          FuseTwoProducts (first, sfp, ffp->input_entityID);
        }
      }
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &context);
  }  
}


static void DoFuseFeature (ButtoN b)

{
  FuseFormPtr  ffp;
  SeqEntryPtr  sep;
  Int2         val;
  ValNodePtr   vnp;

  ffp = (FuseFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  sep = GetTopSeqEntryForEntityID (ffp->input_entityID);
  if (sep == NULL) return;
  Hide (ffp->form);
  WatchCursor ();
  Update ();

  vnp = NULL;
  val = GetValue (ffp->objlist);
  if (val > 0) {
    vnp = ffp->head;
    while (vnp != NULL && val > 1) {
      val--;
      vnp = vnp->next;
    }
  }
  if (vnp != NULL) {
    ffp->subtype = vnp->choice;
    VisitBioseqsInSep (sep, ffp, FuseFeatureCallback);
    DeleteMarkedObjects (ffp->input_entityID, 0, NULL);
  }

  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (ffp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ffp->input_entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  Remove (ffp->form);
}

static void FuseMessageProc (ForM f, Int2 mssg)

{
  FuseFormPtr  ffp;

  ffp = (FuseFormPtr) GetObjectExtra (f);
  if (ffp != NULL) {
    if (ffp->appmessage != NULL) {
      ffp->appmessage (f, mssg);
    }
  }
}

static void CleanupFusePage (GraphiC g, VoidPtr data)

{
  FuseFormPtr  ffp;

  ffp = (FuseFormPtr) data;
  if (ffp != NULL) {
    ValNodeFreeData (ffp->head);
  }
  StdCleanupFormProc (g, data);
}

extern void FuseFeature (IteM i)

{
  BaseFormPtr        bfp;
  ButtoN             b;
  GrouP              c;
  FeatDefPtr         curr;
  FuseFormPtr        ffp;
  GrouP              g;
  GrouP              h;
  ValNodePtr         head;
  Uint1              key;
  CharPtr            label = NULL;
  Int2               listHeight;
  SeqEntryPtr        sep;
  StdEditorProcsPtr  sepp;
  Uint2              subtype;
  ValNodePtr         vnp;
  WindoW             w;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  ffp = (FuseFormPtr) MemNew (sizeof (FuseFormData));
  if (ffp == NULL) return;
  w = FixedWindow (-50, -33, -10, -10, "Fuse Feature", StdCloseWindowProc);
  SetObjectExtra (w, ffp, CleanupFusePage);
  ffp->form = (ForM) w;
  ffp->formmessage = FuseMessageProc;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    ffp->appmessage = sepp->handleMessages;
  }

  ffp->input_entityID = bfp->input_entityID;
  ffp->input_itemID = bfp->input_itemID;
  ffp->input_itemtype = bfp->input_itemtype;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, 0, 2, NULL);
  StaticPrompt (g, "Feature", 0, 0, programFont, 'c');
  if (indexerVersion) {
    listHeight = 16;
  } else {
    listHeight = 8;
  }
  ffp->objlist = SingleList (g, 16, listHeight, NULL);
  head = NULL;
  curr = FeatDefFindNext (NULL, &key, &label, FEATDEF_ANY, TRUE);
  while (curr != NULL) {
    if (key != FEATDEF_BAD) {
      subtype = curr->featdef_key;
      if (subtype != FEATDEF_misc_RNA &&
          subtype != FEATDEF_precursor_RNA &&
          subtype != FEATDEF_mat_peptide &&
          subtype != FEATDEF_sig_peptide &&
          subtype != FEATDEF_transit_peptide &&
          subtype != FEATDEF_Imp_CDS) {
        vnp = ValNodeNew (head);
        if (head == NULL) {
          head = vnp;
        }
        if (vnp != NULL) {
          vnp->choice = subtype;
          vnp->data.ptrvalue = StringSave (curr->typelabel);
        }
      }
    }
    curr = FeatDefFindNext (curr, &key, &label, FEATDEF_ANY, TRUE);
  }
  if (head != NULL) {
    head = SortValNode (head, SortByVnpChoice);
    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      ListItem (ffp->objlist, (CharPtr) vnp->data.ptrvalue);
    }
  }
  ffp->head = head;

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoFuseFeature);
  SetObjectExtra (b, ffp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
}





typedef struct genomeprojiduserdialog {
  DIALOG_MESSAGE_BLOCK
  DialoG        ids;
} GenomeprojidUserDialog, PNTR GenomeprojidUserDialogPtr;

typedef struct genomeprojiduserform {
  FEATURE_FORM_BLOCK
  SeqEntryPtr   sep;
} GenomeprojidUserForm, PNTR GenomeprojidUserFormPtr;

static void UserObjectPtrToGenomeprojidDialog (DialoG d, Pointer data)

{
  Char                       buf [64];
  UserFieldPtr               curr;
  ValNodePtr                 head = NULL;
  ObjectIdPtr                oip;
  Int4                       parentID;
  Int4                       projectID;
  GenomeprojidUserDialogPtr  rdp;
  UserObjectPtr              uop;
  Int4                       val;

  rdp = (GenomeprojidUserDialogPtr) GetObjectExtra (d);
  if (rdp == NULL) return;
  uop = (UserObjectPtr) data;
  if (uop == NULL || uop->type == NULL || StringICmp (uop->type->str, "GenomeProjectsDB") != 0) {
    PointerToDialog (rdp->ids, NULL);
    return;
  }
  projectID = 0;
  parentID = 0;
  for (curr = uop->data; curr != NULL; curr = curr->next) {
    oip = curr->label;
    if (oip == NULL) continue;
    if (StringICmp (oip->str, "ProjectID") == 0) {
      if (curr->choice == 2) {
        val = (Int4) curr->data.intvalue;
        if (projectID > 0) {
          sprintf (buf, "%ld\t%ld", (long) projectID, (long) parentID);
          ValNodeCopyStr (&head, 0, buf);
          parentID = 0;
        }
        projectID = val;
      }
    } else if (StringICmp (oip->str, "ParentID") == 0) {
      if (curr->choice == 2) {
        val = (Int4) curr->data.intvalue;
        parentID = val;
      }
    }
  }
  if (projectID > 0) {
    sprintf (buf, "%ld\t%ld", (long) projectID, (long) parentID);
    ValNodeCopyStr (&head, 0, buf);
  }

  PointerToDialog (rdp->ids, (Pointer) head);
  ValNodeFreeData (head);
}

static Pointer GenomeprojidDialogToUserObjectPtr (DialoG d)

{
  Char                       buf [64];
  ValNodePtr                 head;
  Int4                       parentID;
  Int4                       projectID;
  CharPtr                    ptr1;
  CharPtr                    ptr2;
  GenomeprojidUserDialogPtr  rdp;
  CharPtr                    str;
  UserObjectPtr              uop;
  long int                   val;
  ValNodePtr                 vnp;

  rdp = (GenomeprojidUserDialogPtr) GetObjectExtra (d);
  if (rdp == NULL) return NULL;

  uop = CreateGenomeProjectsDBUserObject ();
  if (uop == NULL) return NULL;

  head = (ValNodePtr) DialogToPointer (rdp->ids);
  if (head == NULL) return NULL;

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    projectID = 0;
    parentID = 0;
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    StringNCpy_0 (buf, str, sizeof (buf));
    ptr1 = StringChr (buf, '\t');
    if (ptr1 != NULL) {
      *ptr1 = '\0';
      ptr1++;
      ptr2 = StringChr (ptr1, '\n');
      if (ptr2 == NULL) {
        ptr2 = StringChr (ptr1, '\t');
      }
      if (ptr2 != NULL) {
        *ptr2 = '\0';
      }
      if (sscanf (buf, "%ld", &val) == 1 && val > 0) {
        projectID = (Int4) val;
        if (sscanf (ptr1, "%ld", &val) == 1 && val > 0) {
          parentID = (Int4) val;
        }
        AddIDsToGenomeProjectsDBUserObject (uop, projectID, parentID);
      }
    }
  }

  ValNodeFreeData (head);

  return uop;
}

static void ValNodePtrToGenomeprojidDialog (DialoG d, Pointer data)

{
  ValNodePtr   head;
  Int2         j;
  ValNodePtr   list;
  CharPtr      str;
  TagListPtr   tlp;
  ValNodePtr   vnp;

  tlp = (TagListPtr) GetObjectExtra (d);
  list = (ValNodePtr) data;
  if (tlp != NULL) {
    head = NULL;
    while (list != NULL) {
      vnp = ValNodeNew (head);
      if (head == NULL) {
        head = vnp;
      }
      if (vnp != NULL) {
        str = MemNew (StringLen ((CharPtr) list->data.ptrvalue) + 3);
        if (str != NULL) {
          StringCpy (str, (CharPtr) list->data.ptrvalue);
          StringCat (str, "\n");
        }
        vnp->data.ptrvalue = str;
      }
      list = list->next;
    }
    SendMessageToDialog (tlp->dialog, VIB_MSG_RESET);
    tlp->vnp = head;
    SendMessageToDialog (tlp->dialog, VIB_MSG_REDRAW);
    for (j = 0, vnp = tlp->vnp; vnp != NULL; j++, vnp = vnp->next) {
    }
    tlp->max = MAX ((Int2) 0, (Int2) (j - tlp->rows + 1));
    CorrectBarMax (tlp->bar, tlp->max);
    CorrectBarPage (tlp->bar, tlp->rows - 1, tlp->rows - 1);
  }
}

static Pointer GenomeprojidDialogToValNodePtr (DialoG d)

{
  Char         ch;
  ValNodePtr   head;
  Int2         j;
  Int2         len;
  ValNodePtr   list;
  Boolean      okay;
  CharPtr      str;
  TagListPtr   tlp;
  ValNodePtr   vnp;

  head = NULL;
  tlp = (TagListPtr) GetObjectExtra (d);
  if (tlp != NULL && tlp->vnp != NULL) {
    list = NULL;
    for (vnp = tlp->vnp; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      okay = FALSE;
      len = StringLen (str);
      for (j = 0; j < len; j++) {
        ch = str [j];
        if (ch != ' ' && ch != '\t' && ch != '\n') {
          okay = TRUE;
        }
      }
      if (okay) {
        list = ValNodeNew (list);
        if (head == NULL) {
          head = list;
        }
        if (list != NULL) {
          list->choice = 0;
          list->data.ptrvalue = StringSave (str);
        }
      }
    }
  }
  return (Pointer) head;
}

Uint2 genproj_types [] = {
  TAGLIST_TEXT, TAGLIST_TEXT
};

Uint2 genproj_widths [] = {
  10, 10, 0
};

static DialoG CreateGenomeProjectsDBDialog (GrouP g)

{
  GrouP                      p;
  GenomeprojidUserDialogPtr  rdp;
  GrouP                      x;
  GrouP                      y;

  p = HiddenGroup (g, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  rdp = (GenomeprojidUserDialogPtr) MemNew (sizeof (GenomeprojidUserDialog));
  if (rdp == NULL) return NULL;

  SetObjectExtra (p, rdp, NULL);
  rdp->dialog = (DialoG) p;
  rdp->todialog = UserObjectPtrToGenomeprojidDialog;
  rdp->fromdialog = GenomeprojidDialogToUserObjectPtr;

  x = HiddenGroup (p, 0, 2, NULL);
  y = HiddenGroup (x, 3, 0, NULL);
  StaticPrompt (y, "Project ID", 10 * stdCharWidth, 0, programFont, 'c');
  StaticPrompt (y, "Parent ID", 10 * stdCharWidth, 0, programFont, 'c');

  rdp->ids = CreateTagListDialog (x, 3, 2, -1,
                                  genproj_types, genproj_widths, NULL,
                                  ValNodePtrToGenomeprojidDialog,
                                  GenomeprojidDialogToValNodePtr);

  return (DialoG) p;
}

static void GenomeProjectsDBUserFormMessage (ForM f, Int2 mssg)

{
  GenomeprojidUserFormPtr  rfp;

  rfp = (GenomeprojidUserFormPtr) GetObjectExtra (f);
  if (rfp != NULL) {
    switch (mssg) {
      case VIB_MSG_CLOSE :
        Remove (f);
        break;
      case VIB_MSG_CUT :
        StdCutTextProc (NULL);
        break;
      case VIB_MSG_COPY :
        StdCopyTextProc (NULL);
        break;
      case VIB_MSG_PASTE :
        StdPasteTextProc (NULL);
        break;
      case VIB_MSG_DELETE :
        StdDeleteTextProc (NULL);
        break;
      default :
        if (rfp->appmessage != NULL) {
          rfp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static ForM CreateGenomeProjectsDBDescForm (Int2 left, Int2 top, Int2 width,
                                           Int2 height, CharPtr title, ValNodePtr sdp,
                                           SeqEntryPtr sep, FormActnFunc actproc)

{
  ButtoN                   b;
  GrouP                    c;
  GrouP                    g;
  GenomeprojidUserFormPtr  rfp;
  StdEditorProcsPtr        sepp;
  WindoW                   w;

  w = NULL;
  rfp = (GenomeprojidUserFormPtr) MemNew (sizeof (GenomeprojidUserForm));
  if (rfp != NULL) {
    w = FixedWindow (left, top, width, height, title, StdCloseWindowProc);
    SetObjectExtra (w, rfp, StdDescFormCleanupProc);
    rfp->form = (ForM) w;
    rfp->actproc = actproc;
    rfp->formmessage = GenomeProjectsDBUserFormMessage;

    rfp->sep = sep;

#ifndef WIN_MAC
    CreateStdEditorFormMenus (w);
#endif
    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      SetActivate (w, sepp->activateForm);
      rfp->appmessage = sepp->handleMessages;
    }

    g = HiddenGroup (w, -1, 0, NULL);
    rfp->data = CreateGenomeProjectsDBDialog (g);

    c = HiddenGroup (w, 2, 0, NULL);
    b = DefaultButton (c, "Accept", StdAcceptFormButtonProc);
    SetObjectExtra (b, rfp, NULL);
    PushButton (c, "Cancel", StdCancelButtonProc);
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
    RealizeWindow (w);
  }
  return (ForM) w;
}

extern Int2 LIBCALLBACK GenomeProjectsDBUserGenFunc (Pointer data);
extern Int2 LIBCALLBACK GenomeProjectsDBUserGenFunc (Pointer data)

{
  ObjectIdPtr              oip;
  OMProcControlPtr         ompcp;
  OMUserDataPtr            omudp;
  ObjMgrProcPtr            proc;
  GenomeprojidUserFormPtr  rfp;
  ValNodePtr               sdp;
  SeqEntryPtr              sep;
  UserObjectPtr            uop;
  WindoW                   w;

  ompcp = (OMProcControlPtr) data;
  w = NULL;
  sdp = NULL;
  sep = NULL;
  uop = NULL;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  proc = ompcp->proc;
  switch (ompcp->input_itemtype) {
    case OBJ_SEQDESC :
      sdp = (ValNodePtr) ompcp->input_data;
      if (sdp != NULL && sdp->choice != Seq_descr_user) {
        return OM_MSG_RET_ERROR;
      }
      uop = (UserObjectPtr) sdp->data.ptrvalue;
      break;
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      break;
    default :
      return OM_MSG_RET_ERROR;
  }
  omudp = ItemAlreadyHasEditor (ompcp->input_entityID, ompcp->input_itemID,
                                ompcp->input_itemtype, ompcp->proc->procid);
  if (omudp != NULL) {
    if (StringCmp (proc->procname, "Edit GenomeProjectsDB User Desc") == 0) {
      rfp = (GenomeprojidUserFormPtr) omudp->userdata.ptrvalue;
      if (rfp != NULL) {
        Select (rfp->form);
      }
      return OM_MSG_RET_DONE;
    } else {
      return OM_MSG_RET_OK; /* not this type, check next registered user object editor */
    }
  }
  if (uop != NULL) {
    oip = uop->type;
    if (oip == NULL || oip->str == NULL) return OM_MSG_RET_OK;
    if (StringCmp (oip->str, "GenomeProjectsDB") != 0) return OM_MSG_RET_OK;
  }
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  w = (WindoW) CreateGenomeProjectsDBDescForm (-50, -33, -10, -10,
                                               "Genome Projects DB", sdp, sep,
                                               StdDescFormActnProc);
  rfp = (GenomeprojidUserFormPtr) GetObjectExtra (w);
  if (rfp != NULL) {
    rfp->input_entityID = ompcp->input_entityID;
    rfp->input_itemID = ompcp->input_itemID;
    rfp->input_itemtype = ompcp->input_itemtype;
    rfp->this_itemtype = OBJ_SEQDESC;
    rfp->this_subtype = Seq_descr_user;
    rfp->procid = ompcp->proc->procid;
    rfp->proctype = ompcp->proc->proctype;
    rfp->userkey = OMGetNextUserKey ();
    omudp = ObjMgrAddUserData (ompcp->input_entityID, ompcp->proc->procid,
	                           OMPROC_EDIT, rfp->userkey);
    if (omudp != NULL) {
      omudp->userdata.ptrvalue = (Pointer) rfp;
      omudp->messagefunc = StdVibrantEditorMsgFunc;
    }
    SendMessageToForm (rfp->form, VIB_MSG_INIT);
    if (sdp != NULL) {
      PointerToDialog (rfp->data, (Pointer) sdp->data.ptrvalue);
      SetClosestParentIfDuplicating ((BaseFormPtr) rfp);
    }
  }
  Show (w);
  Select (w);
  return OM_MSG_RET_DONE;
}



extern Int2 LIBCALLBACK RefGeneUserGenFunc (Pointer data);

#define REFGENE_ASSEMBLY   1
#define REFGENE_RELATED    2
#define REFGENE_SPLICEVAR  3
#define REFGENE_RELDREK    4
#define REFGENE_REJECT     5
#define REFGENE_UNKNOWN    6

typedef struct refgeneuserdialog {
  DIALOG_MESSAGE_BLOCK
  GrouP         status;
  TexT          curator;
  TexT          source;
  Int2          indexer;
  DialoG        fields;
  ButtoN        pipebtn;
} RefgeneUserDialog, PNTR RefgeneUserDialogPtr;

typedef struct refgeneuserform {
  FEATURE_FORM_BLOCK
  SeqEntryPtr   sep;
} RefgeneUserForm, PNTR RefgeneUserFormPtr;

static ENUM_ALIST(changeflags_alist)
  {" ",           0},
  {"Sequence",    1},
  {"Annotation",  2},
  {"Both",        3},
END_ENUM_ALIST

static ENUM_ALIST(refgene_alist)
  {" ",          0},
  {"Assembly",    REFGENE_ASSEMBLY},
  {"Related",     REFGENE_RELATED},
  {"SpliceVar",   REFGENE_SPLICEVAR},
  {"RelatedDrek", REFGENE_RELDREK},
  {"Reject",      REFGENE_REJECT},
  {"Unknown",     REFGENE_UNKNOWN},
END_ENUM_ALIST

static Uint2 refgene_types [] = {
  TAGLIST_TEXT, TAGLIST_TEXT, TAGLIST_TEXT, TAGLIST_POPUP, TAGLIST_POPUP
};

static Uint2 refgene_widths [] = {
  9, 7, 15, 0, 0
};

static EnumFieldAssocPtr refgene_popups [] = {
  NULL, NULL, NULL, changeflags_alist, refgene_alist
};

static CharPtr refgene_labels [] = {
  "", "Assembly", "Related", "SpliceVar", "RelatedDrek", "Reject", "Unknown", NULL
};

static CharPtr refgene_fields [] = {
  "Accession", "GI", "Comment", "Change", "Type", NULL
};

static void AccessionUserFieldPtrToVisStringDialog (DialoG d, Pointer data)

{
  CharPtr       accession;
  Boolean       annotChange;
  CharPtr       comment;
  UserFieldPtr  curr;
  UserFieldPtr  entry;
  Int2          field;
  Int2          flags;
  Int4          from;
  Int4          gi;
  ValNodePtr    head;
  Int2          i;
  Int2          j;
  CharPtr       name;
  ObjectIdPtr   oip;
  Boolean       seqChange;
  CharPtr       str;
  TagListPtr    tlp;
  Int4          to;
  UserFieldPtr  ufp;
  ValNodePtr    vnp;

  tlp = (TagListPtr) GetObjectExtra (d);
  if (tlp == NULL) return;
  str = MemNew (sizeof (Char) * 1024);
  head = NULL;
  curr = (UserFieldPtr) data;
  while (curr != NULL) {
    oip = curr->label;
    if (oip != NULL) {
      field = 0;
      for (i = REFGENE_ASSEMBLY; i <= REFGENE_UNKNOWN; i++) {
        if (StringICmp (oip->str, refgene_labels [i]) == 0 && curr->choice == 11) {
          field = i;
        }
      }
      if (field > 0) {
        entry = (UserFieldPtr) curr->data.ptrvalue;
        while (entry != NULL && entry->choice == 11) {
          accession = NULL;
          comment = NULL;
          name = NULL;
          gi = 0;
          from = 0;
          to = 0;
          annotChange = FALSE;
          seqChange = FALSE;
          ufp = (UserFieldPtr) entry->data.ptrvalue;
          while (ufp != NULL) {
            oip = ufp->label;
            if (oip != NULL && oip->str != NULL) {
              if (StringICmp (oip->str, "accession") == 0 && ufp->choice == 1) {
                accession = (CharPtr) ufp->data.ptrvalue;
              } else if (StringICmp (oip->str, "gi") == 0 && ufp->choice == 2) {
                gi = ufp->data.intvalue;
              } else if (StringICmp (oip->str, "from") == 0 && ufp->choice == 2) {
                from = ufp->data.intvalue;
              } else if (StringICmp (oip->str, "to") == 0 && ufp->choice == 2) {
                to = ufp->data.intvalue;
              } else if (StringICmp (oip->str, "sequenceChange") == 0 && ufp->choice == 4) {
                seqChange = ufp->data.boolvalue;
              } else if (StringICmp (oip->str, "annotationChange") == 0 && ufp->choice == 4) {
                annotChange = ufp->data.boolvalue;
              } else if (StringICmp (oip->str, "comment") == 0 && ufp->choice == 1) {
                comment = (CharPtr) ufp->data.ptrvalue;
              } else if (StringICmp (oip->str, "name") == 0 && ufp->choice == 1) {
                name = (CharPtr) ufp->data.ptrvalue;
              }
            }
            ufp = ufp->next;
          }
          if (accession != NULL) {
            if (comment == NULL) {
              comment = "";
            }
            flags = 0;
            if (seqChange) {
              flags++;
            }
            if (annotChange) {
              flags += 2;
            }
            if (gi != 0) {
              sprintf (str, "%s\t%ld\t%s\t%d\t%d\n", accession, (long) gi, comment, (int) flags, (int) field);
            } else {
              sprintf (str, "%s\t\t%s\t%d\t%d\n", accession, comment, (int) flags, (int) field);
            }
            vnp = ValNodeNew (head);
            if (head == NULL) {
              head = vnp;
            }
            if (vnp != NULL) {
              vnp->data.ptrvalue = StringSave (str);
            }
          } else if (name != NULL) {
            sprintf (str, "\t\t%s\t0\t%d\n", name, (int) field);
            vnp = ValNodeNew (head);
            if (head == NULL) {
              head = vnp;
            }
            if (vnp != NULL) {
              vnp->data.ptrvalue = StringSave (str);
            }
          }
          entry = entry->next;
        }
      }
    }
    curr = curr->next;
  }
  MemFree (str);
  SendMessageToDialog (tlp->dialog, VIB_MSG_RESET);
  tlp->vnp = head;
  SendMessageToDialog (tlp->dialog, VIB_MSG_REDRAW);
  for (j = 0, vnp = tlp->vnp; vnp != NULL; j++, vnp = vnp->next) {
  }
  tlp->max = MAX ((Int2) 0, (Int2) (j - tlp->rows + 1));
  CorrectBarMax (tlp->bar, tlp->max);
  CorrectBarPage (tlp->bar, tlp->rows - 1, tlp->rows - 1);
}

static Pointer VisStringDialogToUserFieldPtr (DialoG d)

{
  return NULL;
}

static void UserObjectPtrToRefGeneDialog (DialoG d, Pointer data)

{
  UserFieldPtr          curr;
  ObjectIdPtr           oip;
  RefgeneUserDialogPtr  rdp;
  Int2                  status = 0;
  CharPtr               str;
  UserObjectPtr         uop;

  rdp = (RefgeneUserDialogPtr) GetObjectExtra (d);
  if (rdp == NULL) return;
  uop = (UserObjectPtr) data;
  if (uop == NULL || uop->type == NULL || StringICmp (uop->type->str, "RefGeneTracking") != 0) {
    SetValue (rdp->status, 0);
    PointerToDialog (rdp->fields, NULL);
    return;
  }
  PointerToDialog (rdp->fields, uop->data);
  for (curr = uop->data; curr != NULL; curr = curr->next) {
    oip = curr->label;
    if (oip != NULL && StringICmp (oip->str, "Status") == 0) {
      break;
    }
  }
  if (curr != NULL && curr->choice == 1) {
    str = (CharPtr) curr->data.ptrvalue;
    if (StringICmp (str, "Inferred") == 0) {
      status = 1;
    } else if (StringICmp (str, "Predicted") == 0) {
      status = 2;
    } else if (StringICmp (str, "Provisional") == 0) {
      status = 3;
    } else if (StringICmp (str, "Validated") == 0) {
      status = 4;
    } else if (StringICmp (str, "Reviewed") == 0) {
      status = 5;
    } else if (StringICmp (str, "Model") == 0) {
      status = 6;
    } else if (StringICmp (str, "WGS") == 0) {
      status = 7;
    } else if (StringICmp (str, "Pipeline") == 0) {
      status = 8;
      SafeEnable (rdp->pipebtn);
    }
  }
  for (curr = uop->data; curr != NULL; curr = curr->next) {
    oip = curr->label;
    if (oip != NULL && StringICmp (oip->str, "Collaborator") == 0) {
      break;
    }
  }
  if (curr != NULL && curr->choice == 1) {
    str = (CharPtr) curr->data.ptrvalue;
    SetTitle (rdp->curator, str);
  }
  for (curr = uop->data; curr != NULL; curr = curr->next) {
    oip = curr->label;
    if (oip != NULL && StringICmp (oip->str, "GenomicSource") == 0) {
      break;
    }
  }
  if (curr != NULL && curr->choice == 1) {
    str = (CharPtr) curr->data.ptrvalue;
    SetTitle (rdp->source, str);
  }
  SetValue (rdp->status, status);
  for (curr = uop->data; curr != NULL; curr = curr->next) {
    oip = curr->label;
    if (oip != NULL && StringICmp (oip->str, "Indexer") == 0) {
      break;
    }
  }
  if (curr != NULL && curr->choice == 2) {
    rdp->indexer = (Int2) curr->data.intvalue;
  }
}

static void AddIndexerToRefGeneTrackUserObject (UserObjectPtr uop, Int2 indexer)

{
  UserFieldPtr  curr;
  ObjectIdPtr   oip;

  if (uop == NULL || indexer < 1) return;
  oip = uop->type;
  if (oip == NULL || StringICmp (oip->str, "RefGeneTracking") != 0) return;

  for (curr = uop->data; curr != NULL; curr = curr->next) {
    oip = curr->label;
    if (oip != NULL && StringICmp (oip->str, "Indexer") == 0) {
      break;
    }
  }

  if (curr == NULL) {
    curr = UserFieldNew ();
    oip = ObjectIdNew ();
    oip->str = StringSave ("Indexer");
    curr->label = oip;
    curr->choice = 2; /* integer */

    /* link indexer at beginning of list */

    curr->next = uop->data;
    uop->data = curr;
  }

  if (curr == NULL || curr->choice != 2) return;

  /* replace any existing indexer indication */

  curr->data.intvalue = (Int4) indexer;
}

static Pointer RefGeneDialogToUserObjectPtr (DialoG d)

{
  Boolean               annotChange;
  Char                  ch;
  Char                  curator [256];
  Int2                  i;
  Int2                  j;
  size_t                len;
  Int4                  num [5];
  Boolean               okay;
  RefgeneUserDialogPtr  rdp;
  Boolean               seqChange;
  Char                  source [64];
  Int2                  status;
  CharPtr               str;
  TagListPtr            tlp;
  CharPtr               txt [5];
  UserObjectPtr         uop;
  long int              val;
  ValNodePtr            vnp;

  rdp = (RefgeneUserDialogPtr) GetObjectExtra (d);
  if (rdp == NULL) return NULL;

  uop = CreateRefGeneTrackUserObject ();
  if (uop == NULL) return NULL;

  status = GetValue (rdp->status);
  if (status == 1) {
    AddStatusToRefGeneTrackUserObject (uop, "Inferred");
  } else if (status == 2) {
    AddStatusToRefGeneTrackUserObject (uop, "Predicted");
  } else if (status == 3) {
    AddStatusToRefGeneTrackUserObject (uop, "Provisional");
  } else if (status == 4) {
    AddStatusToRefGeneTrackUserObject (uop, "Validated");
  } else if (status == 5) {
    AddStatusToRefGeneTrackUserObject (uop, "Reviewed");
  } else if (status == 6) {
    AddStatusToRefGeneTrackUserObject (uop, "Model");
  } else if (status == 7) {
    AddStatusToRefGeneTrackUserObject (uop, "WGS");
  } else if (status == 8) {
    AddStatusToRefGeneTrackUserObject (uop, "Pipeline");
  }

  GetTitle (rdp->source, source, sizeof (source));
  if (! StringHasNoText (source)) {
    AddSourceToRefGeneTrackUserObject (uop, source);
  }

  GetTitle (rdp->curator, curator, sizeof (curator));
  if (! StringHasNoText (curator)) {
    AddCuratorToRefGeneTrackUserObject (uop, curator);
  }

  if (rdp->indexer > 0) {
    AddIndexerToRefGeneTrackUserObject (uop, rdp->indexer);
  }

  tlp = (TagListPtr) GetObjectExtra (rdp->fields);
  if (tlp != NULL && tlp->vnp != NULL) {
    for (vnp = tlp->vnp; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      okay = FALSE;
      len = StringLen (str);
      for (j = 0; j < len; j++) {
        ch = str [j];
        if (ch != ' ' && ch != '\t' && ch != '\n') {
          okay = TRUE;
        }
      }
      if (okay) {
        for (j = 0; j < 5; j++) {
          txt [j] = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, j);
          num [j] = 0;
        }
        for (j = 1; j < 5; j++) {
          if (j != 2) {
            num [j] = 0;
            if (txt [j] != NULL && sscanf (txt [j], "%ld", &val) == 1) {
              num [j] = val;
            }
          }
        }
        annotChange = FALSE;
        seqChange = FALSE;
        if (num [3] >= 2) {
          annotChange = TRUE;
          (num [3]) -= 2;
        }
        if (num [3] > 0) {
          seqChange = TRUE;
        }
        i = num [4];
        if (i >= REFGENE_ASSEMBLY && i <= REFGENE_UNKNOWN) {
          if (! StringHasNoText (txt [0])) {
            AddAccessionToRefGeneTrackUserObject (uop, refgene_labels [i],
                                                  txt [0], num [1],
                                                  seqChange, annotChange,
                                                  txt [2]);
          } else if (! StringHasNoText (txt [2])) {
            /* comment by itself goes into name */
            AddAccessionToRefGeneTrackUserObject (uop, refgene_labels [i],
                                                  NULL, num [1],
                                                  seqChange, annotChange,
                                                  txt [2]);
          }
        }
        for (j = 0; j < 5; j++) {
          txt [j] = MemFree (txt [j]);
        }
      }
    }
  }

  return uop;
}

static DialoG CreateRefGeneDialog (GrouP g)

{
  Int2                  i;
  PrompT                lastppt;
  GrouP                 p;
  PrompT                ppt;
  GrouP                 q;
  RefgeneUserDialogPtr  rdp;
  TagListPtr            tlp;
  GrouP                 x;
  GrouP                 y;

  p = HiddenGroup (g, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  rdp = (RefgeneUserDialogPtr) MemNew (sizeof (RefgeneUserDialog));
  if (rdp == NULL) return NULL;

  SetObjectExtra (p, rdp, NULL);
  rdp->dialog = (DialoG) p;
  rdp->todialog = UserObjectPtrToRefGeneDialog;
  rdp->fromdialog = RefGeneDialogToUserObjectPtr;

  x = HiddenGroup (p, 4, 0, NULL);
  /* StaticPrompt (x, "Status", 0, stdLineHeight, programFont, 'l'); */
  rdp->status = HiddenGroup (x, 8, 0, NULL);
  SetObjectExtra (rdp->status, rdp, NULL);
  RadioButton (rdp->status, "Inferred");
  RadioButton (rdp->status, "Predicted");
  RadioButton (rdp->status, "Provisional");
  RadioButton (rdp->status, "Validated");
  RadioButton (rdp->status, "Reviewed");
  RadioButton (rdp->status, "Model");
  RadioButton (rdp->status, "WGS");
  rdp->pipebtn = RadioButton (rdp->status, "Pipeline");
  Disable (rdp->pipebtn);

  y = HiddenGroup (p, 4, 0, NULL);
  StaticPrompt (y, "Curator", 0, dialogTextHeight, programFont, 'l');
  rdp->curator = DialogText (y, "", 14, NULL);
  StaticPrompt (y, "Genomic Source", 0, dialogTextHeight, programFont, 'l');
  rdp->source = DialogText (y, "", 7, NULL);

  rdp->indexer = 0;

  q = HiddenGroup (p, -6, 0, NULL);
  lastppt = NULL;
  ppt = NULL;
  for (i = 0; i < 5; i++) {
    lastppt = ppt;
    ppt = StaticPrompt (q, refgene_fields [i], refgene_widths [i] * stdCharWidth, 0, systemFont, 'c');
  }
  rdp->fields = CreateTagListDialog (p, 6, 5, STD_TAG_SPACING,
                                     refgene_types, refgene_widths, refgene_popups,
                                     AccessionUserFieldPtrToVisStringDialog,
                                     VisStringDialogToUserFieldPtr);

  tlp = (TagListPtr) GetObjectExtra (rdp->fields);
  if (tlp != NULL) {
    AlignObjects (ALIGN_JUSTIFY, (HANDLE) tlp->control [3], (HANDLE) lastppt, NULL);
    AlignObjects (ALIGN_JUSTIFY, (HANDLE) tlp->control [4], (HANDLE) ppt, NULL);
  }

  AlignObjects (ALIGN_CENTER, (HANDLE) x, (HANDLE) y, (HANDLE) q, (HANDLE) rdp->fields, NULL);
  return (DialoG) p;
}

static void RefgeneUserFormMessage (ForM f, Int2 mssg)

{
  RefgeneUserFormPtr  rfp;

  rfp = (RefgeneUserFormPtr) GetObjectExtra (f);
  if (rfp != NULL) {
    switch (mssg) {
      case VIB_MSG_CLOSE :
        Remove (f);
        break;
      case VIB_MSG_CUT :
        StdCutTextProc (NULL);
        break;
      case VIB_MSG_COPY :
        StdCopyTextProc (NULL);
        break;
      case VIB_MSG_PASTE :
        StdPasteTextProc (NULL);
        break;
      case VIB_MSG_DELETE :
        StdDeleteTextProc (NULL);
        break;
      default :
        if (rfp->appmessage != NULL) {
          rfp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static ForM CreateRefGeneDescForm (Int2 left, Int2 top, Int2 width,
                                   Int2 height, CharPtr title, ValNodePtr sdp,
                                   SeqEntryPtr sep, FormActnFunc actproc)

{
  ButtoN              b;
  GrouP               c;
  GrouP               g;
  RefgeneUserFormPtr  rfp;
  StdEditorProcsPtr   sepp;
  WindoW              w;

  w = NULL;
  rfp = (RefgeneUserFormPtr) MemNew (sizeof (RefgeneUserForm));
  if (rfp != NULL) {
    w = FixedWindow (left, top, width, height, title, StdCloseWindowProc);
    SetObjectExtra (w, rfp, StdDescFormCleanupProc);
    rfp->form = (ForM) w;
    rfp->actproc = actproc;
    rfp->formmessage = RefgeneUserFormMessage;

    rfp->sep = sep;

#ifndef WIN_MAC
    CreateStdEditorFormMenus (w);
#endif
    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      SetActivate (w, sepp->activateForm);
      rfp->appmessage = sepp->handleMessages;
    }

    g = HiddenGroup (w, -1, 0, NULL);
    rfp->data = CreateRefGeneDialog (g);

    c = HiddenGroup (w, 2, 0, NULL);
    b = DefaultButton (c, "Accept", StdAcceptFormButtonProc);
    SetObjectExtra (b, rfp, NULL);
    PushButton (c, "Cancel", StdCancelButtonProc);
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
    RealizeWindow (w);
  }
  return (ForM) w;
}

extern Int2 LIBCALLBACK RefGeneUserGenFunc (Pointer data)

{
  ObjectIdPtr         oip;
  OMProcControlPtr    ompcp;
  OMUserDataPtr       omudp;
  ObjMgrProcPtr       proc;
  RefgeneUserFormPtr  rfp;
  ValNodePtr          sdp;
  SeqEntryPtr         sep;
  UserObjectPtr       uop;
  WindoW              w;

  ompcp = (OMProcControlPtr) data;
  w = NULL;
  sdp = NULL;
  sep = NULL;
  uop = NULL;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  proc = ompcp->proc;
  switch (ompcp->input_itemtype) {
    case OBJ_SEQDESC :
      sdp = (ValNodePtr) ompcp->input_data;
      if (sdp != NULL && sdp->choice != Seq_descr_user) {
        return OM_MSG_RET_ERROR;
      }
      uop = (UserObjectPtr) sdp->data.ptrvalue;
      break;
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      break;
    default :
      return OM_MSG_RET_ERROR;
  }
  omudp = ItemAlreadyHasEditor (ompcp->input_entityID, ompcp->input_itemID,
                                ompcp->input_itemtype, ompcp->proc->procid);
  if (omudp != NULL) {
    if (StringCmp (proc->procname, "Edit RefGene UserTrack Desc") == 0) {
      rfp = (RefgeneUserFormPtr) omudp->userdata.ptrvalue;
      if (rfp != NULL) {
        Select (rfp->form);
      }
      return OM_MSG_RET_DONE;
    } else {
      return OM_MSG_RET_OK; /* not this type, check next registered user object editor */
    }
  }
  if (uop != NULL) {
    oip = uop->type;
    if (oip == NULL || oip->str == NULL) return OM_MSG_RET_OK;
    if (StringCmp (oip->str, "RefGeneTracking") != 0) return OM_MSG_RET_OK;
  }
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  w = (WindoW) CreateRefGeneDescForm (-50, -33, -10, -10,
                                      "Reference Gene Tracking", sdp, sep,
                                      StdDescFormActnProc);
  rfp = (RefgeneUserFormPtr) GetObjectExtra (w);
  if (rfp != NULL) {
    rfp->input_entityID = ompcp->input_entityID;
    rfp->input_itemID = ompcp->input_itemID;
    rfp->input_itemtype = ompcp->input_itemtype;
    rfp->this_itemtype = OBJ_SEQDESC;
    rfp->this_subtype = Seq_descr_user;
    rfp->procid = ompcp->proc->procid;
    rfp->proctype = ompcp->proc->proctype;
    rfp->userkey = OMGetNextUserKey ();
    omudp = ObjMgrAddUserData (ompcp->input_entityID, ompcp->proc->procid,
	                           OMPROC_EDIT, rfp->userkey);
    if (omudp != NULL) {
      omudp->userdata.ptrvalue = (Pointer) rfp;
      omudp->messagefunc = StdVibrantEditorMsgFunc;
    }
    SendMessageToForm (rfp->form, VIB_MSG_INIT);
    if (sdp != NULL) {
      PointerToDialog (rfp->data, (Pointer) sdp->data.ptrvalue);
      SetClosestParentIfDuplicating ((BaseFormPtr) rfp);
    }
  }
  Show (w);
  Select (w);
  return OM_MSG_RET_DONE;
}

/*
static void TestGeneRefStuff (void)

{
  UserObjectPtr uop;
  ValNodePtr    sdp;

  uop = CreateRefGeneTrackUserObject ();
  AddAccessionToRefGeneTrackUserObject (uop, "Assembly", "U12345", 57, 29, 1995);
  AddAccessionToRefGeneTrackUserObject (uop, "Assembly", "L97531", 142, 66, 963);
  AddAccessionToRefGeneTrackUserObject (uop, "Assembly", "M66778", 823, 7677, 343);
  AddAccessionToRefGeneTrackUserObject (uop, "Related", "P34345", 445, 0, 0);
  AddAccessionToRefGeneTrackUserObject (uop, "Reject", "S19635", 1765, 0, 0);
  AddAccessionToRefGeneTrackUserObject (uop, "Related", "Q14884", 664, 35, 97);
  sdp = ValNodeNew (NULL);
  sdp->choice = Seq_descr_user;
  sdp->data.ptrvalue = (Pointer) uop;
  if (! ObjMgrRegister (OBJ_SEQDESC, (Pointer) sdp)) {
     ErrPostEx (SEV_ERROR, 0, 0, "ObjMgrRegister failed.");
  }
}
*/

#define CKA_GAPLEN  50 /* max allowed unaligned gap size */

typedef struct cka_acc {
   CharPtr      accession;
   SeqIdPtr     sip_whole;
   SeqAlignPtr  sap;
   Int4         start_acc;
   Int4         stop_acc;
   Int4         start_seq;
   Int4         stop_seq;
   Uint1        strand;
   Int4         num;
   struct cka_acc PNTR next;
} CKA_Acc, PNTR CKA_AccPtr;

static Int4     CKA_blast_wordsize;
static FloatHi  CKA_blast_expect_value;
static Boolean  CKA_blast_allow_repeats;
static Int4     CKA_blast_detailed_wordsize;
static FloatHi  CKA_blast_detailed_expect_value;
static Boolean  CKA_blast_detailed_allow_repeats;

static SeqAlignPtr CKA_MakeAlign(BioseqPtr bsp, CKA_AccPtr acc_head);

static Boolean SPI_GetAccessionFromSeqId(SeqIdPtr sip, Int4Ptr gi, CharPtr PNTR id)
{
   Boolean numeric_id_type = FALSE;
   Int2 id_len;
   GiimPtr gip;
   ObjectIdPtr oip;
   TextSeqIdPtr textsip;
   DbtagPtr dbtag;
   PatentSeqIdPtr psip;
   PDBSeqIdPtr pdbsip;

   *id = NULL;
   *gi = 0;

   switch (sip->choice) {
   case SEQID_GI: case SEQID_GIBBSQ: case SEQID_GIBBMT:
      *gi = sip->data.intvalue;
      numeric_id_type = TRUE;
      break;
   case SEQID_GIIM:
      gip = (GiimPtr) sip->data.ptrvalue;
      *gi = gip->id;
      numeric_id_type = TRUE;
      break;
   case SEQID_LOCAL:
      oip = (ObjectIdPtr) sip->data.ptrvalue;

      if (oip->str) {
         id_len = StringLen(oip->str);
         *id = (CharPtr) MemNew(id_len+1);
         sprintf(*id, "%s", oip->str);
      } else {
         *id = (CharPtr) MemNew(6);
         sprintf(*id, "%d", oip->id);
      }
      break;
   case SEQID_GENBANK: case SEQID_EMBL: case SEQID_PIR: case SEQID_TPG: case SEQID_TPE: case SEQID_TPD:
   case SEQID_SWISSPROT: case SEQID_DDBJ: case SEQID_PRF:
   case SEQID_OTHER: case SEQID_GPIPE:
      textsip = (TextSeqIdPtr)sip->data.ptrvalue;
      id_len = StringLen(textsip->accession);
      *id = (CharPtr) MemNew(id_len+1);
      if (textsip->version > 0)
         sprintf(*id, "%s.%d", textsip->accession, textsip->version);
      else
         sprintf(*id, "%s", textsip->accession);
      break;
   case SEQID_GENERAL:
      dbtag = (DbtagPtr) sip->data.ptrvalue;
      if (dbtag->tag->str == NULL) {
         numeric_id_type = TRUE;
         *gi = dbtag->tag->id;
      } else {
         id_len = StringLen(dbtag->tag->str);
         *id = (CharPtr) MemNew(id_len+1);
         sprintf(*id, "%s", dbtag->tag->str);
      }
      break;
   case SEQID_PATENT:
      psip = (PatentSeqIdPtr) sip->data.ptrvalue;
      *gi = (Int4) psip->seqid;
      numeric_id_type = TRUE;
      break;
   case SEQID_PDB:
      pdbsip = (PDBSeqIdPtr) sip->data.ptrvalue;
      id_len = StringLen(pdbsip->mol);
      *id = (CharPtr) MemNew(id_len+4);
      sprintf(*id, "%s%d", pdbsip->mol, pdbsip->chain);
      break;
   default: break;
   }

   return numeric_id_type;
}

static void CKA_FindAllTpaDescr(SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
   CKA_AccPtr         acc;
   CKA_AccPtr         PNTR acc_head;
   CKA_AccPtr         acc_prev;
   BioseqPtr          bsp;
   SeqMgrDescContext  context;
   UserFieldPtr       curr;
   ObjectIdPtr        oip;
   SeqDescrPtr        sdp;
   UserFieldPtr       ufp;
   UserObjectPtr      uop;

   acc_head = (CKA_AccPtr PNTR)data;
   acc_prev = *acc_head;
   while (acc_prev != NULL && acc_prev->next != NULL)
   {
      acc_prev = acc_prev->next;
   }
   sdp = NULL;
   if (IS_Bioseq(sep))
   {
      bsp = (BioseqPtr)sep->data.ptrvalue;
      if (ISA_na(bsp->mol))
      {
         while ((sdp = SeqMgrGetNextDescriptor(bsp, sdp, Seq_descr_user, &context)) != NULL)
         {
            uop = (UserObjectPtr)sdp->data.ptrvalue;
            if (!StringICmp(uop->type->str, "TpaAssembly"))
            {
               for (curr = uop->data; curr != NULL; curr = curr->next)
               {
                  if (curr->choice != 11) continue;
                  
                  acc = (CKA_AccPtr)MemNew(sizeof(CKA_Acc));
                  acc->sip_whole = SeqIdSetDup(bsp->id);
                  /* will use these to mark the span for blast2seq */
                  acc->start_acc = acc->stop_acc = -1;
                  if (acc_prev == NULL)
                    *acc_head = acc_prev = acc;
                  else {
                    acc_prev->next = acc;
                    acc_prev = acc;
                  }
                  
                  for (ufp = curr->data.ptrvalue; ufp != NULL; ufp = ufp->next) {
                    oip = ufp->label;
                    if (oip == NULL) continue;
                    if (StringICmp (oip->str, "accession") == 0 && ufp->choice == 1) {
                      acc->accession = StringSave((CharPtr)ufp->data.ptrvalue);
                    } else if (StringICmp (oip->str, "from") == 0 && ufp->choice == 2) {
                      acc->start_acc = (Int4) ufp->data.intvalue;
                    } else if (StringICmp (oip->str, "to") == 0 && ufp->choice == 2) {
                      acc->stop_acc = (Int4) ufp->data.intvalue;
                    }
                  }
               }
            }
         }
      }
   }
}

static int LIBCALLBACK CKA_SortAccs(VoidPtr ptr1, VoidPtr ptr2)
{
   CKA_AccPtr  acc1;
   CKA_AccPtr  acc2;

   acc1 = *((CKA_AccPtr PNTR)ptr1);
   acc2 = *((CKA_AccPtr PNTR)ptr2);
   if (acc1->start_seq < acc2->start_seq)
      return -1;
   else if (acc1->start_seq > acc2->start_seq)
      return 1;
   else if (acc1->stop_seq < acc2->stop_seq)
      return -1;
   else if (acc1->stop_seq > acc2->stop_seq)
      return 1;
   else
      return 0; /* no alignment */
}

static SeqIdPtr SqnSeqIdFindBestAccession (SeqIdPtr sip)
{
	Uint1 order[NUM_SEQID];

	if (sip == NULL)
		return NULL;
	SeqIdBestRank(order, NUM_SEQID);
        order[SEQID_GI]=order[SEQID_LOCAL]+2;
        order[SEQID_PATENT]=order[SEQID_LOCAL]+1;
	return SeqIdSelect (sip, order, NUM_SEQID);
}

static Boolean CKA_ValidateSeqAlign(SeqAlignPtr sap, CKA_AccPtr acc_head, Int4 bioseqlen)
{
   CKA_AccPtr        acc;
   CKA_AccPtr        PNTR accarray;
   AMAlignIndex2Ptr  amaip;
   Int4              first, first_align;
   Boolean           found;
   Int4              gi;
   Int4              i;
   Int4              j;
   Int4              k;
   Int4              last;
   Int4              longest;
   Int4              max;
   Int4              n;
   FILE              PNTR ofp;
   Int4              prev;
   Boolean           retval;
   Char              string[300];
   CharPtr           textid;
   Char              textid2[42];

   if (sap == NULL || sap->saip == NULL || sap->saip->indextype != INDEX_PARENT)
      return FALSE;
   amaip = (AMAlignIndex2Ptr)(sap->saip);
   for (i=0; i<amaip->numsaps; i++)
   {
      acc = acc_head;
      found = FALSE;
      while (acc != NULL && !found)
      {
         if (amaip->saps[i] == acc->sap)
            found = TRUE;
         if (!found)
            acc = acc->next;
      }
      if (!found) /* big error */
         return FALSE;
      acc->num = i+1;
      AlnMgr2GetNthSeqRangeInSA(amaip->saps[i], 1, &acc->start_seq, &acc->stop_seq);
      AlnMgr2GetNthSeqRangeInSA(amaip->saps[i], 2, &acc->start_acc, &acc->stop_acc);
      acc->strand = AlnMgr2GetNthStrand(amaip->saps[i], 2);
      acc->start_seq++;
      acc->stop_seq++;
      acc->start_acc++;
      acc->stop_acc++;
   }
   acc = acc_head;
   i = 0;
   while (acc != NULL)
   {
      if (acc->start_seq == 0 && acc->stop_seq == 0)
      {
         AlnMgr2GetNthSeqRangeInSA(acc->sap, 1, &acc->start_seq, &acc->stop_seq);
         AlnMgr2GetNthSeqRangeInSA(acc->sap, 2, &acc->start_acc, &acc->stop_acc);
         acc->strand = AlnMgr2GetNthStrand(acc->sap, 2);
         acc->start_seq++;
         acc->stop_seq++;
         acc->start_acc++;
         acc->stop_acc++;
      }
      if (acc->num == 0)
         acc->num = amaip->numsaps; /* sort these guys all to the end */
      i++;
      acc = acc->next;
   }
   accarray = (CKA_AccPtr PNTR)MemNew(i*sizeof(CKA_AccPtr));
   i = 0;
   acc = acc_head;
   while (acc != NULL)
   {
      accarray[i] = acc;
      i++;
      acc = acc->next;
   }
   HeapSort(accarray, i, sizeof(CKA_AccPtr), CKA_SortAccs);
   n=0;
   while (accarray[n]->sap == NULL && n < i)
   {
      n++;
   }
   SPI_GetAccessionFromSeqId(SqnSeqIdFindBestAccession(accarray[0]->sip_whole), &gi, &textid);
   if (textid == NULL)
   {
      sprintf(textid2, "%d", gi);
      textid = textid2;
   }
   first = last = -1;
   prev = -1;
   retval = TRUE;
   for (j=0; j<i /*&& first <=0*/ ; j++)
   {
      acc = accarray[j];
      if (acc->sap != NULL)
      {
        if (first == -1) {
          first = acc->start_seq;
          first_align = j;
        }
        last = MAX(last, acc->stop_seq);
      } else {
        continue;
      }
      if (prev != -1)
      {
         if (acc->start_seq > prev + CKA_GAPLEN)
         {
            sprintf(string, "Primary accessions do not completely cover the bioseq %s:\n %s aligns to %d-%d but next aln is %s to %d-%d", textid, accarray[j-1]->accession, accarray[j-1]->start_seq, accarray[j-1]->stop_seq, acc->accession, acc->start_seq, acc->stop_seq);
            ofp = FileOpen("cka_validate.txt", "a");
            fprintf(ofp, "%s", string);
            FileClose(ofp);
            Message(MSG_ERROR, "The following message was written to cka_validate.txt: %s", string);
            retval = FALSE;
         } else if (acc->start_seq > prev)
         {
            sprintf(string, "Primary accessions do not completely cover the bioseq %s:\n %s aligns to %d-%d but the next aln is %s to %d-%d;\n the gap is less than %d and is acceptable.", textid, accarray[j-1]->accession, accarray[j-1]->start_seq, accarray[j-1]->stop_seq, acc->accession, acc->start_seq, acc->stop_seq, CKA_GAPLEN);
            ofp = FileOpen("cka_validate.txt", "a");
            fprintf(ofp, "%s", string);
            FileClose(ofp);
            Message(MSG_ERROR, "The following message was written to cka_validate.txt: %s", string);
            retval = FALSE;
         }
      }
      prev = acc->stop_seq+1;
   }
   if (first != 1 || last != bioseqlen)
   {
      if (first > CKA_GAPLEN)
      {
         sprintf(string, "Primary accessions do not completely cover the bioseq %s:\n %s (the first aln) starts at position %d", textid, accarray[first_align]->accession, accarray[first_align]->start_seq);
         ofp = FileOpen("cka_validate.txt", "a");
         fprintf(ofp, "%s\n", string);
         FileClose(ofp);
         Message(MSG_ERROR, "The following message was written to cka_validate.txt: %s", string);
         retval = FALSE;
      } else if (first != 1)
      {
         Message(MSG_OK, "Primary accessions do not completely cover the bioseq %s:\n %s (the first alignment) starts at position %d, but the gap is less than %d and is acceptable.", textid, accarray[first_align]->accession, accarray[first_align]->start_seq, CKA_GAPLEN);
      }
      max = 0;
      for (k=0; k<i; k++)
      {
         if (accarray[k]->stop_seq > max)
         {
            max = accarray[k]->stop_seq;
            longest = k;
         }
      }
      if (accarray[longest]->stop_seq < bioseqlen-CKA_GAPLEN)
      {
         sprintf (string, "Primary accessions do not completely cover the bioseq %s:\n %s (the last aln) goes to %d, bioseq length is %d", textid, accarray[longest]->accession, accarray[longest]->stop_seq, bioseqlen);
         ofp = FileOpen("cka_validate.txt", "a");
         fprintf(ofp, "%s", string);
         FileClose(ofp);
         Message(MSG_ERROR, "The following message was written to cka_validate.txt: %s", string);
         retval = FALSE;
      } else if (accarray[longest]->stop_seq < bioseqlen)
      {
         Message(MSG_OK, "Primary accessions do not completely cover the bioseq %s:\n %s (the last alignment) goes to %d, bioseq length is %d, but the gap is less than %d and is acceptable", textid, accarray[longest]->accession, accarray[longest]->stop_seq, bioseqlen, CKA_GAPLEN);
      }
      MemFree(accarray);
      accarray = NULL;
   }
   MemFree(accarray);
   return retval;
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

static void CKA_ShowAln(SeqAlignPtr sap, CKA_AccPtr acc_head)
{
   CKA_AccPtr   acc;
   BioseqPtr    bsp;
   DenseSegPtr  dsp;
   Boolean      found;
   GrouP        g;
   Int4         gi;
   Int4         i;
   Int4         len;
   Int4         numsaps;
   SegmenT      picture;
   SeqAlignPtr  salp;
   SeqIdPtr     sip;
   Int4         start;
   Int4         start_r;
   Int4         stop;
   Int4         stop_r;
   Char         tmp[42];
   CharPtr      textid;
   Char         textid2[42];
   VieweR       v;
   WindoW       w;

   w = FixedWindow(-1, -1, -1, -1, "TPA display", StdCloseWindowProc);
   g = HiddenGroup(w, 1, 0, NULL);
   v = CreateViewer(g, 750, 300, TRUE, TRUE);
   picture = CreatePicture();
   salp = (SeqAlignPtr)(sap->segs);
   numsaps = 0;
   while (salp != NULL)
   {
      numsaps++;
      salp = salp->next;
   }
   salp = (SeqAlignPtr)(sap->segs);
   numsaps++;
   dsp = (DenseSegPtr)(salp->segs);
   sip = dsp->ids;
   SPI_GetAccessionFromSeqId(SqnSeqIdFindBestAccession(sip), &gi, &textid);
   if (textid == NULL)
   {
      sprintf(textid2, "%d", gi);
      textid = textid2;
   }
   bsp = BioseqLockById(sip);
   len = bsp->length;
   AddRectangle(picture, 0, numsaps*10, (bsp->length*680)/len, numsaps*10-7, 0, TRUE, 0);
   sprintf(tmp, "1");
   AddLabel(picture, 0-10, numsaps*10-3, tmp, 0, 0, MIDDLE_LEFT, 0);
   sprintf(tmp, "%d  %s", bsp->length, textid);
   AddLabel(picture, ((bsp->length+10)*680)/len, numsaps*10-3, tmp, 0, 0, MIDDLE_RIGHT, 0);
   BioseqUnlock(bsp);
   i = numsaps-1;
   while (salp != NULL)
   {
      acc = acc_head;
      found = FALSE;
      while (acc != NULL && !found)
      {
         if (acc->sap == salp)
            found = TRUE;
         else
            acc = acc->next;
      }
      AlnMgr2GetNthSeqRangeInSA(salp, 1, &start, &stop);
      start_r = (start*680)/len;
      stop_r = (stop*680)/len;
      AddRectangle(picture, start_r, i*10, stop_r, i*10-7, 0, TRUE, 0);
      dsp = (DenseSegPtr)(salp->segs);
      sprintf(tmp, "%d", start+1);
      AddLabel(picture, start_r-10, i*10-3, tmp, 0, 0, MIDDLE_LEFT, 0);
      sprintf(tmp, "%d  %s", stop+1, acc->accession);
      AddLabel(picture, stop_r+10, i*10-3, tmp, 0, 0, MIDDLE_RIGHT, 0);
      salp = salp->next;
      i--;
   }
   AttachPicture(v, picture, 0, 0, UPPER_LEFT, 1, 1, FrameVwr);
   Show(w);
}

static void CKA_RunChecker(SeqEntryPtr sep)
{
   CKA_AccPtr   acc;
   CKA_AccPtr   acc_head;
   CKA_AccPtr   acc_head_next;
   CKA_AccPtr   acc_head_prev;
   CKA_AccPtr   acc_head_real;
   CKA_AccPtr   acc_head_tmp;
   BioseqPtr    bsp;
   Boolean      found;
   Int4         gi;
   SeqIdPtr     lastid;
   SeqAlignPtr  sap;
   SeqHistPtr   shp;
   CharPtr      textid;
   Char         textid2[42];

   if (sep == NULL)
   {
      Message(MSG_ERROR, "Null SeqEntry passed to CKA_RunChecker");
      return;
   }
   acc_head = NULL;
   SeqEntryExplore(sep, &acc_head, CKA_FindAllTpaDescr);
   lastid = NULL;
   if (acc_head == NULL)
   {
      Message(MSG_ERROR, "No Tpa features found in SeqEntry.");
      return;
   }
   acc_head_real = acc_head;
   while (acc_head != NULL)
   {
      lastid = acc_head->sip_whole;
      acc_head_prev = acc_head;
      acc_head_tmp = acc_head->next;
      found = FALSE;
      while (!found && acc_head_tmp != NULL)
      {
         if (SeqIdComp(lastid, acc_head_tmp->sip_whole) != SIC_YES)
            found = TRUE;
         else
         {
            acc_head_prev = acc_head_tmp;
            acc_head_tmp = acc_head_tmp->next;
         }
      }
      acc_head_next = acc_head_prev->next;
      acc_head_prev->next = NULL;
      bsp = BioseqLockById(acc_head->sip_whole);
      if (ISA_na(bsp->mol))
      {
         sap = CKA_MakeAlign(bsp, acc_head);
         acc = acc_head;
         while (acc != NULL && acc->sap == NULL)
         {
            acc = acc->next;
         }
         SPI_GetAccessionFromSeqId(SqnSeqIdFindBestAccession(acc_head->sip_whole), &gi, &textid);
         if (textid == NULL)
         {
            sprintf(textid2, "%d", gi);
            textid = textid2;
         }

         /* report each accession that does not align to the bioseq */
         acc = acc_head;
         while (acc != NULL) {
           if (acc->sap == NULL) {
             Message (MSG_OK, "Accession %s does not align to the bioseq %s.", acc->accession, textid);
           }
           acc = acc->next;
         }

	 if (CKA_ValidateSeqAlign(sap, acc_head, bsp->length))
	    Message(MSG_OK, "Alignments were successfully created and are being added to %s.", textid);
	 else if (sap != NULL)
	    Message(MSG_OK, "Alignments were created but are not valid. They are being added to %s for review.", textid);
         else if (sap == NULL) {
            Message(MSG_ERROR, "No alignments could be created for %s.", textid);
            acc_head = acc_head_next;
            continue;
         }
         AlnMgr2IndexLite(sap);
         AlnMgr2SortAlnSetByNthRowPos(sap, 1);
	 /* make seq-hist and add it to record */
	 if (bsp->hist != NULL)
	 {
	    shp = bsp->hist;
	    if (shp->assembly != NULL)
	       SeqAlignSetFree(shp->assembly);
	    shp->assembly = (SeqAlignPtr)(sap->segs);
	 } else
	 {
	    shp = SeqHistNew();
	    shp->assembly = (SeqAlignPtr)(sap->segs);
	    bsp->hist = shp;
	 }
         sap->segs = NULL;
         SeqAlignFree(sap);
	     } else
	 Message(MSG_ERROR, "%s is annotated on a non-nucleotide bioseq.", acc_head->accession);
    /* while (acc_head != NULL)
     {
        acc_head_tmp = acc_head->next;
        acc_head->next = NULL;
        MemFree(acc_head->accession);
        SeqIdFree(acc_head->sip_whole);
        MemFree(acc_head);
        acc_head = acc_head_tmp;
     }*/
     acc_head = acc_head_next;
  }
  /*CKA_ShowAln(sap, acc_head_real);*/
  while (acc_head_real != NULL)
  {
     acc_head_tmp = acc_head_real->next;
     MemFree(acc_head_real->accession);
     SeqIdFree(acc_head_real->sip_whole);
     MemFree(acc_head_real);
     acc_head_real = acc_head_tmp;
  }
}

static void CKA_FindNuc(SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
   BioseqPtr      bsp;
   BioseqPtr      PNTR bspptr;

   bspptr = (BioseqPtr PNTR)data;
   if (IS_Bioseq(sep))
   {
      bsp = (BioseqPtr)sep->data.ptrvalue;
      if (ISA_na(bsp->mol))
      {
         *bspptr = bsp;
      }
   }
}

static int LIBCALLBACK CKA_CompareAlns(VoidPtr ptr1, VoidPtr ptr2)
{
   Int4         len1;
   Int4         len2;
   SeqAlignPtr  sap1;
   SeqAlignPtr  sap2;

   sap1 = *((SeqAlignPtr PNTR) ptr1);
   sap2 = *((SeqAlignPtr PNTR) ptr2);
   if (sap1 == NULL || sap2 == NULL)
      return 0;
   len1 = AlnMgr2GetAlnLength(sap1, FALSE);
   len2 = AlnMgr2GetAlnLength(sap2, FALSE);
   if (len1 < len2)
      return 1;
   if (len1 > len2)
      return -1;
   return 0;
}

static void amconssetfree(AMConsSetPtr acp)
{
   AMConsSetPtr  acp_next;

   while (acp != NULL)
   {
      acp_next = acp->next;
      MemFree(acp->starts);
      MemFree(acp->stops);
      MemFree(acp->strands);
      MemFree(acp);
      acp = acp_next;
   }
}

static int LIBCALLBACK CKA_SortForConsistent(VoidPtr ptr1, VoidPtr ptr2)
{
   AMConsSetPtr  acp1;
   AMConsSetPtr  acp2;
   FloatHi       bitscore;
   FloatHi       evalue;
   Int4          number;
   SAIndex2Ptr   saip1;
   SAIndex2Ptr   saip2;

   acp1 = *((AMConsSetPtr PNTR)ptr1);
   acp2 = *((AMConsSetPtr PNTR)ptr2);
   saip1 = (SAIndex2Ptr)(acp1->sap->saip);
   saip2 = (SAIndex2Ptr)(acp2->sap->saip);
   if (saip1->score == 0)
      GetScoreAndEvalue(acp1->sap, &saip1->score, &bitscore, &evalue, &number);
   if (saip2->score == 0)
      GetScoreAndEvalue(acp2->sap, &saip2->score, &bitscore, &evalue, &number);
   if (saip1->score > saip2->score)
      return -1;
   else if (saip1->score < saip2->score)
      return 1;
   else
      return 0;
}

static void CKA_RemoveInconsistentAlnsFromSet(SeqAlignPtr sap_head, Int4 fuzz)
{
   AMConsSetPtr  acp;
   AMConsSetPtr  acp_head;
   AMConsSetPtr  acp_prev;
   AMConsSetPtr  PNTR acparray;
   DenseSegPtr   dsp;
   Int4          i;
   Int4          j;
   Int4          k;
   Int4          lfuzz;
   SeqAlignPtr   newsap;
   Int4          numrows;
   Int4          numsaps;
   Int4          orientation;
   Int4          row;
   SAIndex2Ptr   saip;
   SeqAlignPtr   salp_head;
   SeqAlignPtr   salp_prev;
   SeqAlignPtr   sap;
   SeqAlignPtr   sapnext;
   Int4          score;
   SeqIdPtr      sip;
   SeqIdPtr      sip_head;
   Uint1         strand;

   lfuzz = fuzz;
   if (fuzz < 0)
      fuzz = 1;
   sap = (SeqAlignPtr)(sap_head->segs);
   if (sap->next == NULL)
      return;
   dsp = (DenseSegPtr)(sap->segs);
   sip_head = dsp->ids;
   numrows = AlnMgr2GetNumRows(sap);
   acp_head = NULL;
   strand = AlnMgr2GetNthStrand(sap, 1);
   numsaps = 0;
   while (sap != NULL)
   {
      if (AlnMgr2GetNumRows(sap) != numrows)
      {
         amconssetfree(acp_head);
         return;
      }
      numsaps++;
      acp = (AMConsSetPtr)MemNew(sizeof(AMConsSet));
      acp->starts = (Int4Ptr)MemNew(numrows*sizeof(Int4));
      acp->stops = (Int4Ptr)MemNew(numrows*sizeof(Int4));
      acp->strands = (Uint1Ptr)MemNew(numrows*sizeof(Uint1));
      acp->which = (Int4Ptr)MemNew(numrows*sizeof(Int4));
      acp->sap = sap;
      if (acp_head != NULL)
      {
         acp_prev->next = acp;
         acp_prev = acp;
      } else
         acp_head = acp_prev = acp;
      sip = sip_head;
      row = AlnMgr2GetFirstNForSip(sap, sip);
      if (row <= 0)
      {
         amconssetfree(acp_head);
         return;
      }
      if (acp->strands[row] != strand)
      {
         sapnext = acp->sap->next;
         acp->sap->next = NULL;
         score = ((SAIndex2Ptr)(acp->sap->saip))->score;
         SeqAlignListReverseStrand(acp->sap);
         AMAlignIndexFreeEitherIndex(acp->sap);
         AlnMgr2IndexSingleChildSeqAlign(acp->sap);
         saip = (SAIndex2Ptr)(acp->sap->saip);
         saip->score = score;
         acp->strands[row] = strand;
         acp->sap->next = sapnext;
      }
      for (i=0; i<numrows; i++)
      {
         acp->which[i] = row;
         AlnMgr2GetNthSeqRangeInSA(sap, i+1, &acp->starts[i], &acp->stops[i]);
         acp->strands[i] = AlnMgr2GetNthStrand(sap, i+1);
      }
      sap = sap->next;
   }
   acparray = (AMConsSetPtr PNTR)MemNew(numsaps*sizeof(AMConsSetPtr));
   acp = acp_head;
   i = 0;
   while (acp != NULL)
   {
      acparray[i] = acp;
      acp = acp->next;
      i++;
   }
   HeapSort(acparray, numsaps, sizeof(AMConsSetPtr), CKA_SortForConsistent);
   /* orientation -1 means that ith is before jth in ALL rows, 1 means ith is after jth in ALL rows */
   for (i=0; i<numsaps; i++)
   {
      if (acparray[i]->used != -1)
      {
         for (j=i+1; j<numsaps; j++)
         {
            orientation = 0;
            for (k=0; acparray[j]->used != -1 && k<numrows; k++)
            {
               if (acparray[i]->starts[k] - fuzz < acparray[j]->starts[k])
               {
                  if (acparray[i]->stops[k] - fuzz < acparray[j]->starts[k])
                  {
                     if (orientation == 0)
                     {
                        if (acparray[i]->strands[k] == Seq_strand_minus)
                           orientation = 1;
                        else
                           orientation = -1;
                     }
                  } else
                  {
                     if (lfuzz >= 0) /* just mark it for deletion */
                        acparray[j]->used = -1;
                     else /* truncate it */
                     {
                        if (acparray[j]->stops[k] >
                            acparray[i]->stops[k] + CKA_blast_wordsize)
                        {
                           newsap = AlnMgr2GetSubAlign(acparray[j]->sap, acparray[i]->stops[k]+1,
 acparray[j]->stops[k], k+1, TRUE);
                           AlnMgr2IndexSingleChildSeqAlign(newsap);
                           SeqAlignFree(acparray[j]->sap);
                           acparray[j]->sap = newsap;
                           acparray[j]->starts[k] = acparray[i]->stops[k]+1;
                        } else
                           acparray[j]->used = -1;
                     }
                  }
               } else if (acparray[i]->starts[k] - fuzz > acparray[j]->starts[k])
               {
                 if (acparray[i]->starts[k] + fuzz > acparray[j]->stops[k])
                  {
                     if (orientation == 0)
                     {
                        if (acparray[i]->strands[k] == Seq_strand_minus)
                           orientation = -1;
                        else
                           orientation = 1;
                     }
                  } else
                  {
                     if (lfuzz >= 0) /* mark for deletion */
                        acparray[j]->used = -1;
                     else /* truncate */
                     {
                        if (acparray[j]->starts[k] <
                            acparray[i]->starts[k] - CKA_blast_wordsize)
                        {
                           newsap = AlnMgr2GetSubAlign(acparray[j]->sap, acparray[j]->starts[k], acparray[i]->starts[k]-1, k+1, TRUE);
                           AlnMgr2IndexSingleChildSeqAlign(newsap);
                           SeqAlignFree(acparray[j]->sap);
                           acparray[j]->sap = newsap;
                           acparray[j]->stops[k] = acparray[i]->starts[k]-1;
                        } else
                           acparray[j]->used = -1;
                     }
                  }
               } else
                  acparray[j]->used = -1;
            }
         }
      }
   }
   /* now free all the unused ones, stick the rest back together, reindex, and return */
   salp_head = salp_prev = NULL;
   for (i=0; i<numsaps; i++)
   {
      if (acparray[i]->used == -1)
      {
         SeqAlignFree(acparray[i]->sap);
         acparray[i]->sap = NULL;
      } else
      {
         if (salp_head != NULL)
         {
            salp_prev->next = acparray[i]->sap;
            salp_prev = acparray[i]->sap;
            salp_prev->next = NULL;
         } else
         {
            salp_head = salp_prev = acparray[i]->sap;
            salp_prev->next = NULL;
         }
      }
   }
   amconssetfree(acp_head);
   MemFree(acparray);
   sap_head->segs = (Pointer)(salp_head);
   AMAlignIndex2Free2(sap_head->saip);
   AlnMgr2IndexLite(sap_head);
}

static BioseqPtr ReadFromTraceDb (CharPtr number)

{
  BioseqPtr    bsp = NULL;
  CONN         conn;
  time_t       currtime, starttime;
  FILE         *fp;
  Int2         max = 0;
  size_t       n_written;
  Char         path [PATH_MAX];
  Char         query [64];
  SeqEntryPtr  sep = NULL;
  EIO_Status   status;
  STimeout     timeout;
  long int     val;

  if (StringHasNoText (number)) return NULL;
  if (sscanf (number, "%ld", &val) != 1) return NULL;
  sprintf (query, "cmd=raw&query=retrieve+fasta+%ld", (long) val);
  conn = QUERY_OpenUrlQuery ("www.ncbi.nlm.nih.gov", 80, "/Traces/trace.cgi",
                             query, "Sequin", 30, eMIME_T_NcbiData,
                             eMIME_Fasta, eENCOD_None, 0);
  if (conn == NULL) return NULL;
  status = CONN_Write (conn, (const void *) query, StringLen (query),
                       &n_written, eIO_WritePersist);
  if (status != eIO_Success) return NULL;
  QUERY_SendQuery (conn);

#ifdef OS_MAC 
  timeout.sec = 0;
  timeout.usec = 0;
#else
  timeout.sec = 100;
  timeout.usec = 0;
#endif

  starttime = GetSecs ();
  while ((status = CONN_Wait (conn, eIO_Read, &timeout)) != eIO_Success && max < 300) {
    currtime = GetSecs ();
    max = currtime - starttime;
  }

  if (status == eIO_Success) {
    TmpNam (path);
    fp = FileOpen (path, "w");
    QUERY_CopyResultsToFile (conn, fp);
    FileClose (fp);
    /*
    LaunchGeneralTextViewer (path, "QueueFastaQueryToURL results");
    */
    fp = FileOpen (path, "r");
    sep = FastaToSeqEntry (fp, TRUE);
    FileClose (fp);
    FileRemove (path);
    if (sep != NULL) {
      bsp = FindNucBioseq (sep);
    }
  }
  CONN_Close (conn);

  return bsp;
}


static SeqAlignPtr CKA_MakeAlign(BioseqPtr bsp, CKA_AccPtr acc_head)
{
   CKA_AccPtr           acc;
   CKA_AccPtr           acc_new;
   CKA_AccPtr           acc_new_head_head;
   CKA_AccPtr           acc_new_head;
   CKA_AccPtr           acc_new_prev;
   SeqAlignPtr          allsap;
   SeqAlignPtr          allsap_prev;
   AMAlignIndex2Ptr     amaip;
   BioseqPtr            bsp_tmp;
   Int4                 i;
   BLAST_OptionsBlkPtr  options;
   SeqAlignPtr          sap_new;
   SeqAlignPtr          sap_tmp;
   SeqAlignPtr          sap_tmp_next;
   SeqIdPtr             sip;
   Uint1                strand;

   if (bsp == NULL || acc_head == NULL)
      return NULL;
   acc = acc_head;
   allsap = NULL;
   allsap_prev = NULL;
   acc_new_head_head = acc_new_prev = NULL;
   while (acc != NULL)
   {
      options = BLASTOptionNew("blastn", TRUE);
      /*if (bsp->length > 10000)*/
      {
/*    faster with no filtering */
         options->is_megablast_search = TRUE;
         options->gap_open = options->gap_extend = 0;
         options->wordsize = CKA_blast_wordsize;
      }/* else
         options->filter_string = StringSave("m L;R");*/
      options->expect_value = CKA_blast_expect_value;
      bsp_tmp = NULL;
      if (StringNICmp (acc->accession, "ti", 2) == 0) {
        bsp_tmp = ReadFromTraceDb (acc->accession + 2);
      } else {
        sip = SeqIdFromAccessionDotVersion(acc->accession);
        bsp_tmp = BioseqLockById(sip);
      }
      if (bsp_tmp == NULL) break;
      if (bsp_tmp->id->next) {
        /* find the best accession */
        SeqIdPtr sip = SeqIdDup(SqnSeqIdFindBestAccession(bsp_tmp->id));
        bsp_tmp->id = SeqIdSetFree(bsp_tmp->id);
        bsp_tmp->id = sip;
      }
      if (!ISA_na(bsp_tmp->mol))
      {
         BioseqUnlock(bsp_tmp);
         Message(MSG_ERROR, "%s is not a nucleotide bioseq.", acc->accession);
         break;
      }
      WatchCursor();
      if (acc->start_acc >=0 && acc->stop_acc >=0 &&
          acc->start_acc < bsp_tmp->length &&
          acc->start_acc < bsp_tmp->length) {
        SeqLocPtr slp1, slp2;
        if (acc->start_acc <= acc->stop_acc) {
          slp1 = SeqLocIntNew
            (acc->start_acc, acc->stop_acc, Seq_strand_plus, bsp_tmp->id);
        } else {
          slp1 = SeqLocIntNew
            (acc->stop_acc, acc->start_acc, Seq_strand_minus, bsp_tmp->id);
        }
        slp2 = SeqLocIntNew(0, bsp->length-1, Seq_strand_plus, bsp->id);
        acc->sap = BlastTwoSequencesByLoc(slp1, slp2, "blastn", options);
        SeqLocFree(slp1);
        SeqLocFree(slp2);
      } else {
        acc->sap = BlastTwoSequences(bsp_tmp, bsp, "blastn", options);
      }
      BLASTOptionDelete(options);
      ArrowCursor();
      acc->start_acc = acc->stop_acc = 0; /* reset, for later usage */
      if (acc->sap != NULL)
         SPI_flip_sa_list(acc->sap);
      acc_new_head = NULL;
      if (acc->sap != NULL && acc->sap->next != NULL)
      {
         AlnMgr2IndexLite(acc->sap);
         if (!CKA_blast_allow_repeats) {
           CKA_RemoveInconsistentAlnsFromSet(acc->sap, -1);
         }
         sap_tmp = acc->sap;
         acc->sap = (SeqAlignPtr)(acc->sap->segs);
         sap_tmp->segs = NULL;
         SeqAlignFree(sap_tmp);
         sap_tmp = acc->sap->next;
         acc->sap->next = NULL;
         while (sap_tmp != NULL)
         {
            AlnMgr2IndexSingleChildSeqAlign(sap_tmp);
            sap_tmp_next = sap_tmp->next;
            sap_tmp->next = NULL;
            acc_new = (CKA_AccPtr)MemNew(sizeof(CKA_Acc));
            acc_new->accession = StringSave(acc->accession);
            acc_new->sip_whole = SeqIdDup(acc->sip_whole);
            acc_new->sap = sap_tmp;
            sap_tmp = sap_tmp_next;
            if (!acc_new_head) {
              acc_new_head = acc_new;
            }
            if (acc_new_prev) {
              acc_new_prev->next = acc_new;
            }
            acc_new_prev = acc_new;
         }
      } else if (acc->sap != NULL)
         AlnMgr2IndexSingleChildSeqAlign(acc->sap);
      if (acc->sap != NULL)
      {
         strand = AlnMgr2GetNthStrand(acc->sap, 1);
         if (strand == Seq_strand_minus)
         {
            SeqAlignListReverseStrand(acc->sap);
            SAIndex2Free2(acc->sap->saip);
            acc->sap->saip = NULL;
         }
      }
      if (allsap != NULL && acc->sap != NULL)
      {
         allsap_prev->next = acc->sap;
         allsap_prev = allsap_prev->next;;
      } else if (acc->sap != NULL)
         allsap_prev = allsap = (acc->sap);
      acc_new = acc_new_head;
      while (acc_new != NULL)
      {
         strand = AlnMgr2GetNthStrand(acc_new->sap, 1);
         if (strand == Seq_strand_minus)
         {
            SeqAlignListReverseStrand(acc_new->sap);
            SAIndex2Free2(acc_new->sap->saip);
            acc_new->sap->saip = NULL;
         }
         if (allsap != NULL)
         {
            allsap_prev->next = acc_new->sap;
            allsap_prev = allsap_prev->next;;
         } else
            allsap_prev = allsap = acc_new->sap;
         acc_new = acc_new->next;
      }
      if (allsap_prev != NULL)
      {
         while (allsap_prev->next != NULL)
         {
            allsap_prev = allsap_prev->next;
         }
      }  
      BioseqUnlock(bsp_tmp);
      acc = acc->next;
      if (!acc_new_head_head) {
        acc_new_head_head = acc_new_head;
      }
   }
   acc = acc_head;
   while (acc->next != NULL)
   {
      acc = acc->next;
   }
   acc->next = acc_new_head_head;
   if (allsap == NULL)
      return NULL;
   sap_new = SeqAlignNew();
   sap_new->segtype = SAS_DISC;
   sap_new->segs = (Pointer)(allsap);
   allsap = sap_new;
   AlnMgr2IndexLite(allsap);
   AlnMgr2SortAlnSetByNthRowPos(allsap, 1);
   amaip = (AMAlignIndex2Ptr)(allsap->saip);
   for (i=0; i<amaip->numsaps-1; i++)
   {
      amaip->saps[i]->next = amaip->saps[i+1];
   }
   amaip->saps[amaip->numsaps-1]->next = NULL;
   allsap->segs = (Pointer)(amaip->saps[0]);
   return allsap;
}

static void DoCreateSeqHistTPA (IteM i)
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

  CKA_RunChecker(sep);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

extern void CreateSeqHistTPA (IteM i);
extern void CreateSeqHistTPA (IteM i)
{
  CKA_blast_wordsize = 28;
  CKA_blast_expect_value = 0.000001;
  CKA_blast_allow_repeats = FALSE;
  DoCreateSeqHistTPA(i);
}

static SeqAlignPtr DeltaSeq2SeqAlign(BioseqPtr bsp)
{
  DeltaSeqPtr  deltasp;
  SeqLocPtr    slp;
  SeqAlignPtr  sap, sap_head = NULL;
  DenseSegPtr  dsp;
  SeqIntPtr    intp;
  SeqLitPtr    litp;
  int          curr_start = 0;

  if (bsp == NULL || bsp->repr != Seq_repr_delta || bsp->seq_ext_type != 4) {
    return NULL;
  }
  if (!(deltasp = (DeltaSeqPtr) bsp->seq_ext)) {
    return NULL;
  }

  while (deltasp) {
    if (deltasp->choice == 1) { 
      slp = (SeqLocPtr) deltasp->data.ptrvalue;
      if (sap_head) {
        sap = sap->next = SeqAlignNew();
      } else {
        sap_head = sap = SeqAlignNew();
      }
      dsp = DenseSegNew();
      
      sap->type = SAT_PARTIAL;
      sap->segtype = SAS_DENSEG;
      sap->dim = 2;
      sap->segs = (Pointer)(dsp);
      dsp->dim = 2;
      dsp->numseg = 1;
      dsp->lens = (Int4Ptr)MemNew((dsp->numseg)*sizeof(Int4));
      dsp->starts = (Int4Ptr)MemNew((dsp->numseg)*(dsp->dim)*sizeof(Int4));
      dsp->strands = (Uint1Ptr)MemNew((dsp->numseg)*(dsp->dim)*sizeof(Int4));
      
      dsp->ids = SeqIdDup(bsp->id);
      if (dsp->ids->next) {
        /* Dense-seg ids do not support lists, only 1 id per sequence */
        SeqIdFree(dsp->ids->next);
        dsp->ids->next = NULL;
      }
      switch (slp->choice) {
      case SEQLOC_INT:
        intp = (SeqIntPtr) slp->data.ptrvalue;
        dsp->starts[0] = curr_start;
        dsp->starts[1] = intp->from;
        curr_start += dsp->lens[0] = intp->to - intp->from + 1;
        dsp->strands[0] = Seq_strand_plus;
        dsp->strands[1] = intp->strand;
        dsp->ids->next = SeqIdDup(intp->id);
        break;
      default:
        /* exception */
        break;
      }
    } else if (deltasp->choice == 2) { 
      litp = (SeqLitPtr) deltasp->data.ptrvalue;
      if (litp != NULL) {
        curr_start += litp->length;
      }
    }
    deltasp = deltasp->next;
  }
  return sap_head;
}

static void DoDeltaHist (BioseqPtr bsp, Pointer userdata)

{
  SeqAlignPtr  salp;
  SeqHistPtr   shp;

  if (bsp == NULL || bsp->repr != Seq_repr_delta || bsp->seq_ext_type != 4) return;
  shp = bsp->hist;
  if (shp != NULL && shp->assembly != NULL) return;
  salp = DeltaSeq2SeqAlign (bsp);
  if (salp == NULL) return;
  if (shp == NULL) {
    shp = SeqHistNew ();
    bsp->hist = shp;
  }
  if (shp == NULL) return;
  shp->assembly = salp;
}

extern void CreateSeqHistDelta (IteM i);
extern void CreateSeqHistDelta (IteM i)
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

  VisitBioseqsInSep (sep, NULL, DoDeltaHist);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static TexT blast_wordsize_text = NULL;
static TexT blast_expect_value_text = NULL;
static ButtoN blast_allow_repeats_button = NULL;
static IteM blast_i;

static void DoAcceptBlastOptions (ButtoN b)

{
   Char    buf [64];
   long    val1;
   FloatHi val2;

   GetTitle (blast_wordsize_text, buf, sizeof (buf));
   if (sscanf (buf, "%ld", &val1) == 1) {
     CKA_blast_wordsize = CKA_blast_detailed_wordsize = (Int4) val1;
   }
   GetTitle (blast_expect_value_text, buf, sizeof (buf));
   if (sscanf (buf, "%lf", &val2) == 1) {
     CKA_blast_expect_value = CKA_blast_detailed_expect_value = (FloatHi) val2;
   }
   CKA_blast_allow_repeats = CKA_blast_detailed_allow_repeats =
     (Boolean) GetStatus (blast_allow_repeats_button);
   Remove (ParentWindow (b));
   DoCreateSeqHistTPA(blast_i);
}

extern void CreateSeqHistTPADetailed (IteM i);
extern void CreateSeqHistTPADetailed (IteM i)
{
   GrouP   c;
   GrouP   g;
   GrouP   h;
   WindoW  w;
   Char    buf[64];

   blast_i = i;

   w = FixedWindow (-50, -33, -10, -10, "Blast Options", NULL);
   h = HiddenGroup (w, -1, 0, NULL);
   SetGroupSpacing (h, 3, 2);
   g = HiddenGroup (h, 2, 0, NULL);

   StaticPrompt (g, "Word Size", 0, dialogTextHeight, programFont, 'l');
   if (CKA_blast_detailed_wordsize <= 0) {
     CKA_blast_detailed_wordsize = 14;
   }
   CKA_blast_wordsize = CKA_blast_detailed_wordsize;
   sprintf(buf, "%d", CKA_blast_wordsize);
   blast_wordsize_text = DialogText (g, buf, 10, NULL);

   StaticPrompt (g, "Expect Value", 0, dialogTextHeight, programFont, 'l');
   if (CKA_blast_detailed_expect_value <= 0.0) {
     CKA_blast_detailed_expect_value = 0.001;
   }
   CKA_blast_expect_value = CKA_blast_detailed_expect_value;
   sprintf(buf, "%f", CKA_blast_expect_value);
   blast_expect_value_text = DialogText (g, buf, 10, NULL);

   blast_allow_repeats_button = CheckBox (g, "Allow Repeats", NULL);
   if (CKA_blast_detailed_allow_repeats) {
     SetStatus(blast_allow_repeats_button, TRUE);
   }

   c = HiddenGroup (w, 2, 0, NULL);
   SetGroupSpacing (c, 5, 5);
   DefaultButton (c, "Accept", DoAcceptBlastOptions);
   PushButton (c, "Cancel", StdCancelButtonProc);
   AlignObjects (ALIGN_CENTER, (HANDLE) h, (HANDLE) c, NULL);
   RealizeWindow (w);
   Show (w);
   Select (w);
   Select (blast_wordsize_text);
}

extern Int2 LIBCALLBACK AssemblyUserGenFunc (Pointer data);

typedef struct assemblyuserdialog {
  DIALOG_MESSAGE_BLOCK
  DialoG        accns;
} AssemblyUserDialog, PNTR AssemblyUserDialogPtr;

typedef struct assemblyuserform {
  FEATURE_FORM_BLOCK
  SeqEntryPtr   sep;
} AssemblyUserForm, PNTR AssemblyUserFormPtr;

static void UserObjectPtrToAssemblyDialog (DialoG d, Pointer data)

{
  AssemblyUserDialogPtr  adp;
  Char                   buf [16];
  UserFieldPtr           curr;
  Int4                   from;
  ValNodePtr             head = NULL;
  ObjectIdPtr            oip;
  CharPtr                str;
  CharPtr                tmp;
  Int4                   to;
  UserFieldPtr           ufp;
  UserObjectPtr          uop;

  adp = (AssemblyUserDialogPtr) GetObjectExtra (d);
  if (adp == NULL) return;

  uop = (UserObjectPtr) data;
  if (uop == NULL || uop->type == NULL || StringICmp (uop->type->str, "TpaAssembly") != 0) {
    PointerToDialog (adp->accns, NULL);
    return;
  }

  for (curr = uop->data; curr != NULL; curr = curr->next) {
    if (curr->choice != 11) continue;
    str = NULL;
    from = 0;
    to = 0;
    for (ufp = curr->data.ptrvalue; ufp != NULL; ufp = ufp->next) {
      oip = ufp->label;
      if (oip == NULL) continue;
      if (StringICmp (oip->str, "accession") == 0 && ufp->choice == 1) {
        str = (CharPtr) ufp->data.ptrvalue;
      } else if (StringICmp (oip->str, "from") == 0 && ufp->choice == 2) {
        from = (Int4) ufp->data.intvalue;
      } else if (StringICmp (oip->str, "to") == 0 && ufp->choice == 2) {
        to = (Int4) ufp->data.intvalue;
      }
    }
    if (StringHasNoText (str)) continue;
    tmp = MemNew (StringLen (str) + 32);
    StringCpy (tmp, str);
    StringCat (tmp, "\t");
    if (from > 0 || to > 0) {
      sprintf (buf, "%ld", (long) (from + 1));
      StringCat (tmp, buf);
      StringCat (tmp, "\t");
      sprintf (buf, "%ld", (long) (to + 1));
      StringCat (tmp, buf);
      StringCat (tmp, "\t");
    } else {
      StringCat (tmp, "\t\t");
    }
    ValNodeAddStr (&head, 0, (Pointer) tmp);
  }

  PointerToDialog (adp->accns, (Pointer) head);
  ValNodeFreeData (head);
}

static void DoAddAccessionToTpa (UserObjectPtr uop, CharPtr last)

{
  Int4      from = 0;
  CharPtr   ptr1, ptr2, ptr3;
  Int4      to = 0;
  long int  val;

  ptr1 = StringChr (last, '\t');
  if (ptr1 != NULL) {
    *ptr1 = '\0';
    ptr1++;
    ptr2 = StringChr (ptr1, '\t');
    if (ptr2 != NULL) {
      *ptr2 = '\0';
      ptr2++;
      ptr3 = StringChr (ptr2, '\t');
      if (ptr3 != NULL) {
        *ptr3 = '\0';
      }
    }
    if (sscanf (ptr1, "%ld", &val) == 1 && val > 0) {
      from = val - 1;
      if (sscanf (ptr2, "%ld", &val) == 1 && val > 0) {
        to = val - 1;
      } else {
        from = 0;
        to = 0;
      }
    }
  }
  AddAccessionToTpaAssemblyUserObject (uop, last, from, to);
}

static Pointer AssemblyDialogToUserObjectPtr (DialoG d)

{
  AssemblyUserDialogPtr  adp;
  Char                   ch;
  ValNodePtr             head;
  CharPtr                last;
  UserObjectPtr          uop;
  CharPtr                ptr;
  CharPtr                str;
  CharPtr                tmp;
  ValNodePtr             vnp;

  adp = (AssemblyUserDialogPtr) GetObjectExtra (d);
  if (adp == NULL) return NULL;

  uop = CreateTpaAssemblyUserObject ();
  if (uop == NULL) return NULL;

  head = (ValNodePtr) DialogToPointer (adp->accns);
  if (head == NULL) return NULL;

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    tmp = StringSave (str);
    last = tmp;
    ptr = last;
    ch = *ptr;
    while (ch != '\0') {
      if (ch == ',' || ch == ' ') {
        *ptr = '\0';
        TrimSpacesAroundString (last);
        DoAddAccessionToTpa (uop, last);
        ptr++;
        last = ptr;
        ch = *ptr;
      } else {
        ptr++;
        ch = *ptr;
      }
    }
    if (! StringHasNoText (last)) {
      TrimSpacesAroundString (last);
      DoAddAccessionToTpa (uop, last);
    }
    MemFree (tmp);
  }

  ValNodeFreeData (head);

  return uop;
}

static void ValNodePtrToAssemblyDialog (DialoG d, Pointer data)

{
  ValNodePtr   head;
  Int2         j;
  ValNodePtr   list;
  CharPtr      str;
  TagListPtr   tlp;
  ValNodePtr   vnp;

  tlp = (TagListPtr) GetObjectExtra (d);
  list = (ValNodePtr) data;
  if (tlp != NULL) {
    head = NULL;
    while (list != NULL) {
      vnp = ValNodeNew (head);
      if (head == NULL) {
        head = vnp;
      }
      if (vnp != NULL) {
        str = MemNew (StringLen ((CharPtr) list->data.ptrvalue) + 3);
        if (str != NULL) {
          StringCpy (str, (CharPtr) list->data.ptrvalue);
          StringCat (str, "\n");
        }
        vnp->data.ptrvalue = str;
      }
      list = list->next;
    }
    SendMessageToDialog (tlp->dialog, VIB_MSG_RESET);
    tlp->vnp = head;
    SendMessageToDialog (tlp->dialog, VIB_MSG_REDRAW);
    for (j = 0, vnp = tlp->vnp; vnp != NULL; j++, vnp = vnp->next) {
    }
    tlp->max = MAX ((Int2) 0, (Int2) (j - tlp->rows + 1));
    CorrectBarMax (tlp->bar, tlp->max);
    CorrectBarPage (tlp->bar, tlp->rows - 1, tlp->rows - 1);
  }
}

static Pointer AssemblyDialogToValNodePtr (DialoG d)

{
  Char         ch;
  ValNodePtr   head;
  Int2         j;
  Int2         len;
  ValNodePtr   list;
  Boolean      okay;
  CharPtr      str;
  TagListPtr   tlp;
  ValNodePtr   vnp;

  head = NULL;
  tlp = (TagListPtr) GetObjectExtra (d);
  if (tlp != NULL && tlp->vnp != NULL) {
    list = NULL;
    for (vnp = tlp->vnp; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      okay = FALSE;
      len = StringLen (str);
      for (j = 0; j < len; j++) {
        ch = str [j];
        if (ch != ' ' && ch != '\t' && ch != '\n') {
          okay = TRUE;
        }
      }
      if (okay) {
        list = ValNodeNew (list);
        if (head == NULL) {
          head = list;
        }
        if (list != NULL) {
          list->choice = 0;
          list->data.ptrvalue = StringSave ((CharPtr) vnp->data.ptrvalue);
        }
      }
    }
  }
  return (Pointer) head;
}

Uint2 assmbly_types [] = {
  TAGLIST_TEXT, TAGLIST_TEXT, TAGLIST_TEXT
};

Uint2 assmbly_widths [] = {
  16, 8, 8, 0
};

static DialoG CreateAssemblyDialog (GrouP g)

{
  AssemblyUserDialogPtr  adp;
  GrouP                  p;
  GrouP                  x;
  GrouP                  y;

  p = HiddenGroup (g, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  adp = (AssemblyUserDialogPtr) MemNew (sizeof (AssemblyUserDialog));
  if (adp == NULL) return NULL;

  SetObjectExtra (p, adp, NULL);
  adp->dialog = (DialoG) p;
  adp->todialog = UserObjectPtrToAssemblyDialog;
  adp->fromdialog = AssemblyDialogToUserObjectPtr;

  x = HiddenGroup (p, 0, 2, NULL);
  y = HiddenGroup (x, 3, 0, NULL);
  StaticPrompt (y, "Accessions", 16 * stdCharWidth, 0, programFont, 'c');
  StaticPrompt (y, "From", 8 * stdCharWidth, 0, programFont, 'c');
  StaticPrompt (y, "To", 8 * stdCharWidth, 0, programFont, 'c');
  adp->accns = CreateTagListDialog (x, 3, 3, -1,
                                    assmbly_types, assmbly_widths, NULL,
                                    ValNodePtrToAssemblyDialog,
                                    AssemblyDialogToValNodePtr);

  return (DialoG) p;
}

static void AssemblyUserFormMessage (ForM f, Int2 mssg)

{
  AssemblyUserFormPtr  afp;

  afp = (AssemblyUserFormPtr) GetObjectExtra (f);
  if (afp != NULL) {
    switch (mssg) {
      case VIB_MSG_CLOSE :
        Remove (f);
        break;
      case VIB_MSG_CUT :
        StdCutTextProc (NULL);
        break;
      case VIB_MSG_COPY :
        StdCopyTextProc (NULL);
        break;
      case VIB_MSG_PASTE :
        StdPasteTextProc (NULL);
        break;
      case VIB_MSG_DELETE :
        StdDeleteTextProc (NULL);
        break;
      default :
        if (afp->appmessage != NULL) {
          afp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static ForM CreateAssemblyDescForm (Int2 left, Int2 top, Int2 width,
                                   Int2 height, CharPtr title, ValNodePtr sdp,
                                   SeqEntryPtr sep, FormActnFunc actproc)

{
  AssemblyUserFormPtr  afp;
  ButtoN               b;
  GrouP                c;
  GrouP                g;
  StdEditorProcsPtr    sepp;
  WindoW               w;

  w = NULL;
  afp = (AssemblyUserFormPtr) MemNew (sizeof (AssemblyUserForm));
  if (afp != NULL) {
    w = FixedWindow (left, top, width, height, title, StdCloseWindowProc);
    SetObjectExtra (w, afp, StdDescFormCleanupProc);
    afp->form = (ForM) w;
    afp->actproc = actproc;
    afp->formmessage = AssemblyUserFormMessage;

    afp->sep = sep;

#ifndef WIN_MAC
    CreateStdEditorFormMenus (w);
#endif
    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      SetActivate (w, sepp->activateForm);
      afp->appmessage = sepp->handleMessages;
    }

    g = HiddenGroup (w, -1, 0, NULL);
    afp->data = CreateAssemblyDialog (g);

    c = HiddenGroup (w, 2, 0, NULL);
    b = DefaultButton (c, "Accept", StdAcceptFormButtonProc);
    SetObjectExtra (b, afp, NULL);
    PushButton (c, "Cancel", StdCancelButtonProc);
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
    RealizeWindow (w);
  }
  return (ForM) w;
}

extern Int2 LIBCALLBACK AssemblyUserGenFunc (Pointer data)

{
  AssemblyUserFormPtr  afp;
  ObjectIdPtr          oip;
  OMProcControlPtr     ompcp;
  OMUserDataPtr        omudp;
  ObjMgrProcPtr        proc;
  ValNodePtr           sdp;
  SeqEntryPtr          sep;
  UserObjectPtr        uop;
  WindoW               w;

  ompcp = (OMProcControlPtr) data;
  w = NULL;
  sdp = NULL;
  sep = NULL;
  uop = NULL;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  proc = ompcp->proc;
  switch (ompcp->input_itemtype) {
    case OBJ_SEQDESC :
      sdp = (ValNodePtr) ompcp->input_data;
      if (sdp != NULL && sdp->choice != Seq_descr_user) {
        return OM_MSG_RET_ERROR;
      }
      uop = (UserObjectPtr) sdp->data.ptrvalue;
      break;
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      break;
    default :
      return OM_MSG_RET_ERROR;
  }
  omudp = ItemAlreadyHasEditor (ompcp->input_entityID, ompcp->input_itemID,
                                ompcp->input_itemtype, ompcp->proc->procid);
  if (omudp != NULL) {
    if (StringCmp (proc->procname, "Edit Assembly User Desc") == 0) {
      afp = (AssemblyUserFormPtr) omudp->userdata.ptrvalue;
      if (afp != NULL) {
        Select (afp->form);
      }
      return OM_MSG_RET_DONE;
    } else {
      return OM_MSG_RET_OK; /* not this type, check next registered user object editor */
    }
  }
  if (uop != NULL) {
    oip = uop->type;
    if (oip == NULL || oip->str == NULL) return OM_MSG_RET_OK;
    if (StringCmp (oip->str, "TpaAssembly") != 0) return OM_MSG_RET_OK;
  }
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  w = (WindoW) CreateAssemblyDescForm (-50, -33, -10, -10,
                                       "Assembly Tracking", sdp, sep,
                                       StdDescFormActnProc);
  afp = (AssemblyUserFormPtr) GetObjectExtra (w);
  if (afp != NULL) {
    afp->input_entityID = ompcp->input_entityID;
    afp->input_itemID = ompcp->input_itemID;
    afp->input_itemtype = ompcp->input_itemtype;
    afp->this_itemtype = OBJ_SEQDESC;
    afp->this_subtype = Seq_descr_user;
    afp->procid = ompcp->proc->procid;
    afp->proctype = ompcp->proc->proctype;
    afp->userkey = OMGetNextUserKey ();
    omudp = ObjMgrAddUserData (ompcp->input_entityID, ompcp->proc->procid,
	                           OMPROC_EDIT, afp->userkey);
    if (omudp != NULL) {
      omudp->userdata.ptrvalue = (Pointer) afp;
      omudp->messagefunc = StdVibrantEditorMsgFunc;
    }
    SendMessageToForm (afp->form, VIB_MSG_INIT);
    if (sdp != NULL) {
      PointerToDialog (afp->data, (Pointer) sdp->data.ptrvalue);
      SetClosestParentIfDuplicating ((BaseFormPtr) afp);
    }
  }
  Show (w);
  Select (w);
  return OM_MSG_RET_DONE;
}


typedef struct historyformdata {
  FEATURE_FORM_BLOCK

  BioseqPtr      bsp;
  DialoG         replace_date;
  DialoG         replace_ids;
  ButtoN         secondary_on_part;
  DialoG         replaced_by_date;
  DialoG         replaced_by_ids;
  ButtoN         deleted;
  DialoG         deleted_date;
} HistoryFormData, PNTR HistoryFormPtr;

static SeqIdPtr VisStrDialogToSeqIdSet (DialoG d)

{
  long        gi;
  SeqIdPtr    head = NULL;
  ValNodePtr  list;
  SeqIdPtr    sip;
  CharPtr     str;
  ValNodePtr  vnp;

  if (d == NULL) return NULL;
  list = DialogToPointer (d);
  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (str != NULL) {
      if (sscanf (str, "%ld", &gi)) {
        /*
        ValNodeAddInt (&head, SEQID_GI, (Int4) gi);
        */
      } else {
        sip = SeqIdFromAccessionDotVersion (str);
        if (sip != NULL) {
          ValNodeLink (&head, sip);
        }
      }
    }
  }
  ValNodeFreeData (list);
  return head;
}

static int LIBCALLBACK SortByName (VoidPtr ptr1, VoidPtr ptr2)

{
  CharPtr     str1;
  CharPtr     str2;
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      str1 = (CharPtr) vnp1->data.ptrvalue;
      str2 = (CharPtr) vnp2->data.ptrvalue;
      if (str1 != NULL && str2 != NULL) {
        return StringICmp (str1, str2);
      } else {
        return 0;
      }
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}

static ValNodePtr GetStringsForSeqIDs (SeqIdPtr sip)

{
  Char          buf [40];
  ValNodePtr    head = NULL;
  TextSeqIdPtr  tsip;

  if (sip == NULL) return NULL;
  while (sip != NULL) {
    buf [0] = '\0';
    switch (sip->choice) {
      case SEQID_GENBANK :
      case SEQID_EMBL :
      case SEQID_DDBJ :
      case SEQID_OTHER :
        tsip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tsip != NULL && (! StringHasNoText (tsip->accession))) {
          StringNCpy_0 (buf, tsip->accession, sizeof (buf));
        }
        break;
      case SEQID_GI :
        /*
        gi = sip->data.intvalue;
        if (gi > 0) {
          sprintf (buf, "%ld", (long) gi);
        }
        */
        break;
      default :
        break;
    }
    if (! StringHasNoText (buf)) {
      ValNodeCopyStr (&head, 0, buf);
    }
    sip = sip->next;
  }
  return head;
}

static void AddGenBankBlockToBioseq (BioseqPtr bsp, ValNodePtr head1, ValNodePtr head2)

{
  GBBlockPtr       gbp = NULL;
  CharPtr          last = NULL;
  ValNodePtr       next;
  ValNodePtr PNTR  prev;
  ValNodePtr       sdp;
  SeqEntryPtr      sep;
  CharPtr          str1, str2;
  ValNodePtr       vnp, vnp1, vnp2;

  sdp = BioseqGetSeqDescr (bsp, Seq_descr_genbank, NULL);
  if (sdp != NULL) {
    gbp = (GBBlockPtr) sdp->data.ptrvalue;
    if (gbp != NULL) {
      for (vnp1 = head1; vnp1 != NULL; vnp1 = vnp1->next) {
        str1 = (CharPtr) vnp1->data.ptrvalue;
        if (str1 != NULL) {
          for (vnp2 = gbp->extra_accessions; vnp2 != NULL; vnp2 = vnp2->next) {
            str2 = (CharPtr) vnp2->data.ptrvalue;
            if (StringICmp (str1, str2) == 0) {
              vnp2->data.ptrvalue = MemFree (vnp2->data.ptrvalue);
            }
          }
        }
      }
    }
  }
  if (sdp == NULL) {
    sep = SeqMgrGetSeqEntryForData (bsp);
    sdp = CreateNewDescriptor (sep, Seq_descr_genbank);
    if (sdp != NULL) {
      sdp->data.ptrvalue = GBBlockNew ();
    }
  }
  if (sdp != NULL) {
    gbp = (GBBlockPtr) sdp->data.ptrvalue;
    if (gbp != NULL) {
      while (head2 != NULL) {
        ValNodeCopyStr (&(gbp->extra_accessions), 0, (CharPtr) head2->data.ptrvalue);
        head2 = head2->next;
      }
      /*
      ValNodeLink (&(gbp->extra_accessions), head2);
      head2 = NULL;
      */
      gbp->extra_accessions = SortValNode (gbp->extra_accessions, SortByName);
      prev = &(gbp->extra_accessions);
      vnp = gbp->extra_accessions;
      last = NULL;
      while (vnp != NULL) {
        next = vnp->next;
        str2 = (CharPtr) vnp->data.ptrvalue;
        if (str2 == NULL || StringHasNoText (str2) || StringICmp (last, str2) == 0) {
          *prev = next;
          vnp->next = NULL;
          MemFree (vnp);
          vnp = next;
        } else {
          last = str2;
          prev = &(vnp->next);
          vnp = next;
        }
      }
    }
  }

}

static void DoChangeHistory (ButtoN b)

{
  MsgAnswer       ans;
  BioseqPtr       bsp;
  ValNodePtr      head1 = NULL, head2 = NULL;
  HistoryFormPtr  hfp;
  SeqHistPtr      hist;
  BioseqPtr       pbsp;
  SeqEntryPtr     sep;
  SeqLocPtr       slp;
  CharPtr         str1, str2;
  ValNodePtr      vnp1, vnp2;

  hfp = (HistoryFormPtr) GetObjectExtra (b);
  if (hfp == NULL) return;
  ans = Message (MSG_OKC, "Are you sure you want to edit the history?");
  if (ans == ANS_CANCEL) {
    return;
  }
  Hide (hfp->form);
  Update ();
  bsp = hfp->bsp;
  hist = bsp->hist;
  if (hist == NULL) {
    hist = SeqHistNew ();
    bsp->hist = hist;
  }
  if (hist != NULL) {

    hist->replace_date = DateFree (hist->replace_date);
    hist->replace_date = DialogToPointer (hfp->replace_date);
    head1 = GetStringsForSeqIDs (hist->replace_ids);
    hist->replace_ids = SeqIdSetFree (hist->replace_ids);
    hist->replace_ids = VisStrDialogToSeqIdSet (hfp->replace_ids);
    head2 = GetStringsForSeqIDs (hist->replace_ids);

    hist->replaced_by_date = DateFree (hist->replaced_by_date);
    hist->replaced_by_date = DialogToPointer (hfp->replaced_by_date);
    hist->replaced_by_ids = SeqIdSetFree (hist->replaced_by_ids);
    hist->replaced_by_ids = VisStrDialogToSeqIdSet (hfp->replaced_by_ids);

    hist->deleted = GetStatus (hfp->deleted);
    hist->deleted_date = DateFree (hist->deleted_date);
    hist->deleted_date = DialogToPointer (hfp->deleted_date);
  }

  if (hist->assembly == NULL &&
      hist->replace_date == NULL && hist->replace_ids == NULL &&
      hist->replaced_by_date == NULL && hist->replaced_by_ids == NULL &&
      (! hist->deleted) && hist->deleted_date == NULL) {
    bsp->hist = SeqHistFree (bsp->hist);
  }

  head1 = SortValNode (head1, SortByName);
  head2 = SortValNode (head2, SortByName);
  for (vnp1 = head1; vnp1 != NULL; vnp1 = vnp1->next) {
    str1 = (CharPtr) vnp1->data.ptrvalue;
    for (vnp2 = head2; vnp2 != NULL; vnp2 = vnp2->next) {
      str2 = (CharPtr) vnp2->data.ptrvalue;
      if (StringICmp (str1, str2) == 0) {
        vnp1->data.ptrvalue = MemFree (vnp1->data.ptrvalue);
      }
    }
  }

  AddGenBankBlockToBioseq (bsp, head1, head2);

  if (GetStatus (hfp->secondary_on_part)) {
    if (bsp->repr == Seq_repr_seg) {
      for (slp = (SeqLocPtr) bsp->seq_ext; slp != NULL; slp = slp->next) {
        pbsp = BioseqFind (SeqLocId (slp));
        if (pbsp != NULL) {
          AddGenBankBlockToBioseq (pbsp, head1, head2);
        }
      }
    }
  }

  ValNodeFreeData (head1);
  ValNodeFreeData (head2);

  sep = GetTopSeqEntryForEntityID (hfp->input_entityID);
  EntryCheckGBBlock (sep);

  Update ();
  ObjMgrSetDirtyFlag (hfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, hfp->input_entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  Remove (hfp->form);
}

static void SeqIdSetToVisStrDialog (DialoG d, SeqIdPtr sip)

{
  ValNodePtr  head = NULL;

  if (d == NULL || sip == NULL) return;
  head = GetStringsForSeqIDs (sip);
  PointerToDialog (d, head);
  ValNodeFreeData (head);
}

static void HistoryMessageProc (ForM f, Int2 mssg)

{
  HistoryFormPtr  hfp;

  hfp = (HistoryFormPtr) GetObjectExtra (f);
  if (hfp != NULL) {
    if (hfp->appmessage != NULL) {
      hfp->appmessage (f, mssg);
    }
  }
}

static void CleanupHistoryPage (GraphiC g, VoidPtr data)

{
  HistoryFormPtr  hfp;

  hfp = (HistoryFormPtr) data;
  if (hfp != NULL) {
  }
  StdCleanupFormProc (g, data);
}

extern void EditSequenceHistory (IteM i)

{
  ButtoN             b;
  BaseFormPtr        bfp;
  BioseqPtr          bsp;
  GrouP              c;
  GrouP              g;
  GrouP              h;
  HistoryFormPtr     hfp;
  SeqHistPtr         hist;
  GrouP              j;
  GrouP              k;
  PrompT             ppt1, ppt2, ppt3;
  SeqEntryPtr        sep;
  StdEditorProcsPtr  sepp;
  WindoW             w;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  bsp = GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
  if (bsp == NULL) return;
  hist = bsp->hist;

  hfp = (HistoryFormPtr) MemNew (sizeof (HistoryFormData));
  if (hfp == NULL) return;
  w = FixedWindow (-50, -33, -10, -10, "Sequence History", StdCloseWindowProc);
  SetObjectExtra (w, hfp, CleanupHistoryPage);
  hfp->form = (ForM) w;
  hfp->formmessage = HistoryMessageProc;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    hfp->appmessage = sepp->handleMessages;
  }

  hfp->input_entityID = bfp->input_entityID;
  hfp->input_itemID = bfp->input_itemID;
  hfp->input_itemtype = bfp->input_itemtype;

  hfp->bsp = bsp;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  ppt1 = StaticPrompt (h, "Replaces", 0, 0, systemFont, 'c');

  g = HiddenGroup (h, -1, 0, NULL);
  hfp->replace_ids = CreateVisibleStringDialog (g, 4, -1, 10);
  hfp->secondary_on_part = CheckBox (g, "Secondary on Parts", NULL);
  hfp->replace_date = CreateDateDialog (g, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) hfp->replace_ids, (HANDLE) hfp->secondary_on_part, (HANDLE) hfp->replace_date, NULL);

  ppt2 = StaticPrompt (h, "Replaced By", 0, 0, systemFont, 'c');

  j = HiddenGroup (h, -1, 0, NULL);
  hfp->replaced_by_ids = CreateVisibleStringDialog (j, 4, -1, 10);
  hfp->replaced_by_date = CreateDateDialog (j, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) hfp->replaced_by_ids, (HANDLE) hfp->replaced_by_date, NULL);

  ppt3 = StaticPrompt (h, "Status", 0, 0, systemFont, 'c');

  k = HiddenGroup (h, -1, 0, NULL);
  hfp->deleted = CheckBox (k, "Deleted", NULL);
  hfp->deleted_date = CreateDateDialog (k, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) hfp->deleted, (HANDLE) hfp->deleted_date, NULL);

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoChangeHistory);
  SetObjectExtra (b, hfp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) j, (HANDLE) k, (HANDLE) c,
                (HANDLE) ppt1, (HANDLE) ppt2, (HANDLE) ppt3, NULL);
  RealizeWindow (w);

  if (bsp->repr != Seq_repr_seg) {
    Disable (hfp->secondary_on_part);
  } else {
    SetStatus (hfp->secondary_on_part, TRUE);
  }

  if (hist != NULL) {
    PointerToDialog (hfp->replace_date, hist->replace_date);
    SeqIdSetToVisStrDialog (hfp->replace_ids, hist->replace_ids);
    PointerToDialog (hfp->replaced_by_date, hist->replaced_by_date);
    SeqIdSetToVisStrDialog (hfp->replaced_by_ids, hist->replaced_by_ids);
    PointerToDialog (hfp->deleted_date, hist->deleted_date);
    SetStatus (hfp->deleted, hist->deleted);
  }

  Show (w);
  Update ();
}

/* automatic defline generator */

typedef struct deffeats {
  SeqFeatPtr  sfp;
  SeqFeatPtr  gene;
  SeqFeatPtr  prot;
  CharPtr     genename;
  CharPtr     allelename;
  CharPtr     protname;
  Boolean     alreadyTrimmed;
  Uint2       entityID;
  Uint2       itemID;
  Uint2       subtype;
  Boolean     isDNA;
  Boolean     isAlleleGroup;
  Boolean     lastInString;
  Boolean     lastInGroup;
  Boolean     lastInType;
  Boolean     lastInPenultimate;
  Boolean     pseudo;
  Boolean     ignore;
  Boolean     suppressprefix;
  Int2        altSplices;
  Int2        numUnknown;
} DefFeatsData, PNTR DefFeatsPtr;

static Boolean GetMolBioFeatsGatherFunc (GatherContextPtr gcp, Boolean getGene, Boolean getSnoRNA)

{
  DefFeatsPtr  dfp;
  RnaRefPtr    rrp;
  SeqFeatPtr   sfp;
  CharPtr      str;
  Uint1        type;
  ValNodePtr   PNTR vnpp;

  if (gcp == NULL || gcp->thisitem == NULL || gcp->userdata == NULL)
    return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  vnpp = (ValNodePtr PNTR) gcp->userdata;
  sfp = (SeqFeatPtr) gcp->thisitem;
  switch (sfp->data.choice) {
    case SEQFEAT_GENE :
      if (getGene) {
        dfp = MemNew (sizeof (DefFeatsData));
        if (dfp == NULL) return TRUE;
        dfp->entityID = gcp->entityID;
        dfp->itemID = gcp->itemID;
        dfp->sfp = sfp;
        dfp->subtype = FEATDEF_GENE;
        ValNodeAddPointer (vnpp, 0, (Pointer) dfp);
      }
      break;
    case SEQFEAT_CDREGION :
      dfp = MemNew (sizeof (DefFeatsData));
      if (dfp == NULL) return TRUE;
      dfp->entityID = gcp->entityID;
      dfp->itemID = gcp->itemID;
      dfp->sfp = sfp;
      dfp->subtype = FEATDEF_CDS;
      dfp->altSplices = 1;
      ValNodeAddPointer (vnpp, 0, (Pointer) dfp);
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp == NULL) return TRUE;
      switch (rrp->type) {
        case 3 :
          dfp = MemNew (sizeof (DefFeatsData));
          if (dfp == NULL) return TRUE;
          dfp->entityID = gcp->entityID;
          dfp->itemID = gcp->itemID;
          dfp->sfp = sfp;
          dfp->subtype = FEATDEF_tRNA;
          ValNodeAddPointer (vnpp, 0, (Pointer) dfp);
          break;
        case 4 :
          dfp = MemNew (sizeof (DefFeatsData));
          if (dfp == NULL) return TRUE;
          dfp->entityID = gcp->entityID;
          dfp->itemID = gcp->itemID;
          dfp->sfp = sfp;
          dfp->subtype = FEATDEF_rRNA;
          ValNodeAddPointer (vnpp, 0, (Pointer) dfp);
          break;
         case 5 :
          if (getSnoRNA) {
            dfp = MemNew (sizeof (DefFeatsData));
            if (dfp == NULL) return TRUE;
            dfp->entityID = gcp->entityID;
            dfp->itemID = gcp->itemID;
            dfp->sfp = sfp;
            dfp->subtype = FEATDEF_snRNA;
            ValNodeAddPointer (vnpp, 0, (Pointer) dfp);
          }
          break;
       case 7 :
          if (getSnoRNA) {
            dfp = MemNew (sizeof (DefFeatsData));
            if (dfp == NULL) return TRUE;
            dfp->entityID = gcp->entityID;
            dfp->itemID = gcp->itemID;
            dfp->sfp = sfp;
            dfp->subtype = FEATDEF_snoRNA;
            ValNodeAddPointer (vnpp, 0, (Pointer) dfp);
          }
          break;
        case 255 :
          if (rrp->ext.choice == 1) {
            str = (CharPtr) rrp->ext.value.ptrvalue;
            if (StringICmp (str, "internal transcribed spacer 1") == 0 ||
                StringICmp (str, "internal transcribed spacer 2") == 0 ||
                StringICmp (str, "internal transcribed spacer 3") == 0 ||
                StringICmp (str, "internal transcribed spacer ITS1") == 0 ||
                StringICmp (str, "internal transcribed spacer ITS2") == 0 ||
                StringICmp (str, "internal transcribed spacer ITS3") == 0 ||
                StringICmp (str, "ITS1") == 0 ||
                StringICmp (str, "ITS2") == 0 ||
                StringICmp (str, "ITS3") == 0) {
              dfp = MemNew (sizeof (DefFeatsData));
              if (dfp == NULL) return TRUE;
              dfp->entityID = gcp->entityID;
              dfp->itemID = gcp->itemID;
              dfp->sfp = sfp;
              dfp->subtype = FEATDEF_otherRNA;
              ValNodeAddPointer (vnpp, 0, (Pointer) dfp);
            }
          }
          break;
        default :
          break;
      }
      break;
    case SEQFEAT_IMP :
      type = FindFeatDefType (sfp);
      if (type == FEATDEF_LTR || type == FEATDEF_exon) {
        dfp = MemNew (sizeof (DefFeatsData));
        if (dfp == NULL) return TRUE;
        dfp->entityID = gcp->entityID;
        dfp->itemID = gcp->itemID;
        dfp->sfp = sfp;
        dfp->subtype = type;
        ValNodeAddPointer (vnpp, 0, (Pointer) dfp);
      }
      break;
    default :
      break;
  }
  return TRUE;
}

static Boolean GetCDStRNArRNAGatherFunc (GatherContextPtr gcp)

{
  return GetMolBioFeatsGatherFunc (gcp, FALSE, FALSE);
}


extern EnumFieldAssoc  orgmod_subtype_alist [];
extern EnumFieldAssoc  subsource_subtype_alist [];

extern void StringToLower (CharPtr str)

{
  Char  ch;

  if (str == NULL) return;
  ch = *str;
  while (ch != '\0') {
    *str = TO_LOWER (ch);
    str++;
    ch = *str;
  }
}

static CharPtr molinfo_tech_list [] = {
  "?", "standard", "EST", "STS", "survey", "genetic map", "physical map",
  "derived", "concept-trans", "seq-pept", "both", "seq-pept-overlap",
  "seq-pept-homol", "concept-trans-a", "htgs 1", "htgs 2", "htgs 3",
  "fli cDNA", "htgs 0", "htc", "wgs", "barcode", "composite-wgs-htgs", NULL
};

static void MakeNucleotideTitlesInSequinStyle (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  EnumFieldAssocPtr  ap;
  BioseqContextPtr   bcp;
  BioSourcePtr       biop;
  BioseqPtr          bsp;
  MolInfoPtr         mip;
  OrgModPtr          mod;
  OrgNamePtr         onp;
  OrgRefPtr          orp;
  SeqDescrPtr        sdp;
  SubSourcePtr       ssp;
  CharPtr            str;
  Uint1              tech = 0;
  Char               text [256];
  ValNodePtr         ttl;
  ValNodePtr         vnp;

  if (sep == NULL) return;
  if (! IS_Bioseq (sep)) return;
  bsp = sep->data.ptrvalue;
  if (bsp == NULL) return;
  if (! ISA_na (bsp->mol)) return;
  bcp = BioseqContextNew (bsp);
  sdp = BioseqContextGetSeqDescr (bcp, Seq_descr_source, NULL, NULL);
  BioseqContextFree (bcp);
  if (sdp == NULL) return;
  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop == NULL) return;
  if (bsp->descr != NULL) {
    vnp = ValNodeExtract (&(bsp->descr), Seq_descr_title);
    vnp = ValNodeFreeData (vnp);
  }
  bcp = BioseqContextNew (bsp);
  sdp = BioseqContextGetSeqDescr (bcp, Seq_descr_molinfo, NULL, NULL);
  if (sdp != NULL) {
    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip != NULL) {
      switch (mip->tech) {
        case MI_TECH_est :
        case MI_TECH_sts :
        case MI_TECH_survey :
        case MI_TECH_htgs_1 :
        case MI_TECH_htgs_2 :
        case MI_TECH_htgs_3 :
        case MI_TECH_fli_cdna :
        case MI_TECH_htgs_0 :
        case MI_TECH_htc :
        case MI_TECH_wgs :
          tech = mip->tech;
          break;
        default :
          break;
      }
    }
  }
  BioseqContextFree (bcp);
  str = MemNew (2000);

  orp = biop->org;
  if (orp != NULL) {
    StringCpy (text, "[organism=");
    StringCat (text, orp->taxname);
    StringCat (text, "] ");
    StringCat (str, text);
  }

  ssp = biop->subtype;
  while (ssp != NULL) {
    for (ap = subsource_subtype_alist; ap->name != NULL; ap++) {
      if (ssp->subtype == ap->value) {
        StringCpy (text, "[");
        if (ap->value != 255) {
          StringCat (text, ap->name);
        } else {
          StringCat (text, "subsource");
        }
        StringToLower (text);
        StringCat (text, "=");
        StringCat (text, ssp->name);
        StringCat (text, "] ");
        StringCat (str, text);
      }
    }
    ssp = ssp->next;
  }
  orp = biop->org;
  if (orp != NULL) {
    onp = orp->orgname;
    if (onp != NULL) {
      mod = onp->mod;
      while (mod != NULL) {
        for (ap = orgmod_subtype_alist; ap->name != NULL; ap++) {
          if (mod->subtype == ap->value) {
            StringCpy (text, "[");
            if (ap->value != 255) {
              StringCat (text, ap->name);
            } else {
              StringCat (text, "note");
            }
            StringToLower (text);
            StringCat (text, "=");
            StringCat (text, mod->subname);
            StringCat (text, "] ");
            StringCat (str, text);
          }
        }
        mod = mod->next;
      }
    }
  }

  if (tech > 0) {
    StringCpy (text, "[tech=");
    StringCat (text, molinfo_tech_list [tech]);
    StringCat (text, "] ");
    StringCat (str, text);
  }

  TrimSpacesAroundString (str);
  if (! StringHasNoText (str)) {
    ttl = CreateNewDescriptor (sep, Seq_descr_title);
    if (ttl != NULL) {
      ttl->data.ptrvalue = StringSave (str);
    }
  }
  MemFree (str);
}

extern Int2 LIBCALLBACK MakeSequinProteinTitles (Pointer data);
extern Int2 LIBCALLBACK MakeSequinNucleotideTitles (Pointer data);
extern Int2 LIBCALLBACK MakeSequinFeatureTable (Pointer data);
extern Int2 LIBCALLBACK MakeContigBuildTable (Pointer data);

extern Int2 LIBCALLBACK MakeSequinNucleotideTitles (Pointer data)

{
  OMProcControlPtr  ompcp;
  SeqEntryPtr       sep;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      return OM_MSG_RET_ERROR;
    default :
      return OM_MSG_RET_ERROR;
  }
  if (ompcp->input_data == NULL) return OM_MSG_RET_ERROR;
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  SeqEntryExplore (sep, NULL, MakeNucleotideTitlesInSequinStyle);
  ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, 0, 0);
  return OM_MSG_RET_DONE;
}

static Boolean StringHasEqualSignOrBrackets (CharPtr str)

{
  Char  ch;

  if (StringHasNoText (str)) return FALSE;
  ch = *str;
  while (ch != '\0') {
    if (ch == '=' || ch == '[' || ch == ']') return TRUE;
    str++;
    ch = *str;
  }
  return FALSE;
}

static void MakeProteinTitlesInSequinStyle (Uint2 entityID, SeqEntryPtr sep)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  DefFeatsPtr   dfp;
  Char          quot [4];
  GeneRefPtr    grp;
  GatherScope   gs;
  Boolean       has_equal_or_brackets;
  ValNodePtr    head;
  SeqEntryPtr   nsep;
  ProtRefPtr    prp;
  SeqEntryPtr   psep;
  Char          str [256];
  Char          text [256];
  ValNodePtr    ttl;
  ValNodePtr    vnp;

  if (sep == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    if (bssp->_class == 7 ||
        (IsPopPhyEtcSet (bssp->_class))) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        MakeProteinTitlesInSequinStyle (entityID, sep);
      }
      return;
    }
  }
  nsep = FindNucSeqEntry (sep);
  if (nsep == NULL) return;

  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = TRUE;
  MemSet ((Pointer) (gs.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.scope = sep;
  gs.target = NULL;
  head = NULL;
  GatherEntity (entityID, (Pointer) (&head), GetCDStRNArRNAGatherFunc, &gs);
  /* head = SortValNode (head, SortCDStRNArRNAByLocation); */
  if (head == NULL) return;

  quot [0] = '"';
  quot [1] = '\0';

  vnp = head;
  while (vnp != NULL) {
    dfp = (DefFeatsPtr) vnp->data.ptrvalue;
    if (dfp != NULL && dfp->sfp != NULL && dfp->subtype == FEATDEF_CDS) {
      FindGeneAndProtForCDS (entityID, dfp->sfp, &(dfp->gene), &(dfp->prot));
      bsp = GetBioseqGivenSeqLoc (dfp->sfp->product, entityID);
      if (bsp != NULL) {
        str [0] = '\0';
        if (dfp->gene != NULL) {
          grp = (GeneRefPtr) dfp->gene->data.value.ptrvalue;
          if (grp != NULL) {
            StringNCpy_0 (text, (CharPtr) grp->locus, sizeof (text));
            if (! StringHasNoText (text)) {
              StringCat (str, "[gene=");
              StringCat (str, text);
              StringCat (str, "]");
            }
            if (grp->syn != NULL) {
              StringNCpy_0 (text, (CharPtr) grp->syn->data.ptrvalue, sizeof (text));
              if (! StringHasNoText (text)) {
                if (str [0] != '\0') {
                  StringCat (str, " ");
                }
                StringCat (str, "[gene_syn=");
                StringCat (str, text);
                StringCat (str, "]");
              }
            }
          }
        }
        if (dfp->prot != NULL) {
          prp = (ProtRefPtr) dfp->prot->data.value.ptrvalue;
          if (prp != NULL) {
            if (prp->name != NULL) {
              StringNCpy_0 (text, (CharPtr) prp->name->data.ptrvalue, sizeof (text));
              if (! StringHasNoText (text)) {
                if (str [0] != '\0') {
                  StringCat (str, " ");
                }
                StringCat (str, "[protein=");
                has_equal_or_brackets = StringHasEqualSignOrBrackets (text);
                if (has_equal_or_brackets) {
                  StringCat (str, quot);
                }
                StringCat (str, text);
                if (has_equal_or_brackets) {
                  StringCat (str, quot);
                }
                StringCat (str, "]");
              }
            }
            StringNCpy_0 (text, (CharPtr) prp->desc, sizeof (text));
            if (! StringHasNoText (text)) {
              StringCat (str, "[prot_desc=");
              has_equal_or_brackets = StringHasEqualSignOrBrackets (text);
              if (has_equal_or_brackets) {
                StringCat (str, quot);
              }
              StringCat (str, text);
              if (has_equal_or_brackets) {
                StringCat (str, quot);
              }
              StringCat (str, "]");
            }
          }
        }
        if (! StringHasNoText (str)) {
          psep = SeqMgrGetSeqEntryForData (bsp);
          if (psep != NULL) {
            ttl = CreateNewDescriptor (psep, Seq_descr_title);
            if (ttl != NULL) {
              ttl->data.ptrvalue = StringSave (str);
            }
          }
        }
      }
    }
    vnp = vnp->next;
  }
  ValNodeFreeData (head);
}

extern Int2 LIBCALLBACK MakeSequinProteinTitles (Pointer data)

{
  OMProcControlPtr  ompcp;
  SeqEntryPtr       sep;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      return OM_MSG_RET_ERROR;
    default :
      return OM_MSG_RET_ERROR;
  }
  if (ompcp->input_data == NULL) return OM_MSG_RET_ERROR;
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  MakeProteinTitlesInSequinStyle (ompcp->input_entityID, sep);
  ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, 0, 0);
  return OM_MSG_RET_DONE;
}

static Boolean LIBCALLBACK SequinFTableBioseq (BioseqPtr bsp, SeqMgrBioseqContextPtr context)

{
  FILE      *fp;

  if (bsp == NULL) return TRUE;
  fp = (FILE *) context->userdata;
  BioseqToGnbk (bsp, NULL, FTABLE_FMT, DUMP_MODE, NORMAL_STYLE, 0, 0, 0, NULL, fp);
  return TRUE;
}

extern Int2 LIBCALLBACK MakeSequinFeatureTable (Pointer data)

{
  FILE              *fp;
  OMProcControlPtr  ompcp;
  Char              path [PATH_MAX];
  SeqEntryPtr       sep;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      return OM_MSG_RET_ERROR;
    default :
      return OM_MSG_RET_ERROR;
  }
  if (ompcp->input_data == NULL) return OM_MSG_RET_ERROR;
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  TmpNam (path);
  fp = FileOpen (path, "w");
  if (fp == NULL) return OM_MSG_RET_ERROR;
  SeqMgrExploreBioseqs (ompcp->input_entityID, NULL, (Pointer) fp, SequinFTableBioseq, TRUE, FALSE, FALSE);
  FileClose (fp);
  LaunchGeneralTextViewer (path, "Gene - CDS Feature Table");
  FileRemove (path);
  return OM_MSG_RET_DONE;
}

/* the following two functions are modified from PrintGenome in wprint.c */

static void SqLitPrintGenome(SeqLitPtr slp, FILE *fp)
{
	if (slp->seq_data != NULL)         /* not a gap */
	{
		if (slp->length == 0)  /* unknown length */
		{
			fprintf(fp, "gap\t0\t0\t0\n");
		} else {
/* don't know what to do here */
		}
	} else {                  /* gap length was set */
			fprintf(fp, "gap\t0\t0\t%ld\n", (long) slp->length);
	}
}

static void SqLocPrintGenome(SeqLocPtr slp_head, FILE *fp)
{
	SeqLocPtr	slp;
	static Char		buf[11];
	SeqIdPtr	sid, newid;
	Int4 		start, stop;
		
	for (slp = slp_head; slp; slp = slp->next) {
		sid = SeqLocId(slp);
		if (slp->choice == SEQLOC_INT || slp->choice == SEQLOC_WHOLE) {
			start = SeqLocStart(slp);
			stop = SeqLocStop(slp);
		} else if (slp->choice == SEQLOC_NULL){
			fprintf(fp, "gap\t0\t0\t0\n");
			continue;
		} else {
			continue;
		}
		if (sid->choice == SEQID_GI) {
			newid = GetSeqIdForGI(sid->data.intvalue);
		} else if (sid->choice == SEQID_GENERAL) {
			newid = sid;
		} else {
			newid = sid;
		}
		SeqIdWrite(newid, buf, PRINTID_TEXTID_ACCESSION, 10);
		fprintf(fp, "%s\t%ld\t%ld\n", buf, (long) start+1, (long) stop+1);
	}
}

static Boolean DeltaLitOnly (BioseqPtr bsp)

{
  ValNodePtr  vnp;

  if (bsp == NULL || bsp->repr != Seq_repr_delta) return FALSE;
  for (vnp = (ValNodePtr)(bsp->seq_ext); vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == 1) return FALSE;
  }
  return TRUE;
}

extern Int2 LIBCALLBACK MakeContigBuildTable (Pointer data)

{
  BioseqPtr         bsp;
  DeltaSeqPtr       dsp;
  FILE              *fp;
  SeqLitPtr 	    litp;
  OMProcControlPtr  ompcp;
  Char              path [PATH_MAX];
  SeqEntryPtr       sep;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  bsp = NULL;
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQ :
      bsp = (BioseqPtr) ompcp->input_data;
      break;
    default :
      return OM_MSG_RET_ERROR;
  }
  if (bsp == NULL) return OM_MSG_RET_ERROR;
  sep = GetBestTopParentForData (ompcp->input_entityID, bsp);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  if (IsAGenomeRecord (sep) ||
      IsSegmentedBioseqWithoutParts (sep)) {
  } else if (IsADeltaBioseq (sep) && (! DeltaLitOnly (bsp))) {
  } else return OM_MSG_RET_ERROR;
  TmpNam (path);
  fp = FileOpen (path, "w");
  if (fp == NULL) return OM_MSG_RET_ERROR;
/* the following code is modified from PrintGenome in wprint.c */
	if (bsp->seq_ext_type == 1) {
		SqLocPrintGenome((SeqLocPtr) bsp->seq_ext, fp);
	} else if (bsp->seq_ext_type == 4) {
		for (dsp = (DeltaSeqPtr) bsp->seq_ext; dsp; dsp=dsp->next) {
			if (dsp->choice == 1) {  /* SeqLoc */
				SqLocPrintGenome((SeqLocPtr)(dsp->data.ptrvalue), fp);
			} else {
				litp = (SeqLitPtr)(dsp->data.ptrvalue);
				if (litp == NULL) continue;
				SqLitPrintGenome(litp, fp);
			}
		}
	}
  FileClose (fp);
  LaunchGeneralTextViewer (path, "Contig Build Table");
  FileRemove (path);
  return OM_MSG_RET_DONE;
}

static void InsertGeneLocusTagPrefixCallback (SeqFeatPtr sfp, Pointer userdata)
{
  CharPtr    prefix;
  GeneRefPtr grp;
  CharPtr    new_locus_tag = NULL;
  Int4       new_locus_tag_len = 0;
  
  if (sfp == NULL || userdata == NULL || sfp->data.choice != SEQFEAT_GENE)
  {
    return;
  }
  
  prefix = (CharPtr) userdata;
  
  if (StringHasNoText (prefix))
  {
    return;
  }
  
  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (grp == NULL)
  {
    grp = GeneRefNew();
    sfp->data.value.ptrvalue = grp;
  }
  if (grp == NULL)
  {
    return;
  }
  
  if (StringHasNoText (grp->locus_tag))
  {
    grp->locus_tag = MemFree (grp->locus_tag);
    grp->locus_tag = StringSave (prefix);
  }
  else
  {
    new_locus_tag_len = StringLen (prefix) + StringLen (grp->locus_tag) + 1;
    new_locus_tag = (CharPtr) MemNew (sizeof (Char) * new_locus_tag_len);
    if (new_locus_tag != NULL)
    {
      StringCpy (new_locus_tag, prefix);
      StringCat (new_locus_tag, grp->locus_tag);
      grp->locus_tag = MemFree (grp->locus_tag);
      grp->locus_tag = new_locus_tag;
    }
  }
} 

typedef struct locustagprefix
{
  FEATURE_FORM_BLOCK

  TexT prefix_txt;  
} LocusTagPrefixData, PNTR LocusTagPrefixPtr;

static void InsertGeneLocusTagPrefixButton (ButtoN b)
{
  LocusTagPrefixPtr ltpp;
  CharPtr           prefix;
  SeqEntryPtr       sep;
  
  ltpp = (LocusTagPrefixPtr) GetObjectExtra (b);
  if (ltpp == NULL)
  {
    return;
  }
  
  prefix = SaveStringFromText (ltpp->prefix_txt);
  if (!StringHasNoText (prefix))
  {
    sep = GetTopSeqEntryForEntityID (ltpp->input_entityID);
    if (sep != NULL)
    {
      VisitFeaturesInSep (sep, prefix, InsertGeneLocusTagPrefixCallback); 
    }
    
  }
  prefix = MemFree (prefix);
  ObjMgrSetDirtyFlag (ltpp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ltpp->input_entityID, 0, 0);
  Remove (ltpp->form);
  Update ();
}

extern void InsertGeneLocusTagPrefix (IteM i)
{
  BaseFormPtr        bfp;
  WindoW             w;
  LocusTagPrefixPtr  ltpp;
  GrouP              h, g, c;
  ButtoN             b;
  
#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  
  ltpp = (LocusTagPrefixPtr) MemNew (sizeof (LocusTagPrefixData));
  if (ltpp == NULL) return;
  ltpp->input_entityID = bfp->input_entityID;

  w = FixedWindow (-50, -33, -10, -10, "Feature Evidence", StdCloseWindowProc);
  SetObjectExtra (w, ltpp, StdCleanupFormProc);
  ltpp->form = (ForM) w;
  h = HiddenGroup (w, -1, 0, NULL);
  g = HiddenGroup (h, 2, 0, NULL);
  StaticPrompt (g, "Prefix for Gene Locus Tag", 0, 0, programFont, 'c');
  ltpp->prefix_txt = DialogText (g, "", 14, NULL);
  c = HiddenGroup (h, 2, 0, NULL);
  b = PushButton (c, "Accept", InsertGeneLocusTagPrefixButton);
  SetObjectExtra (b, ltpp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();    
}


static CharPtr MakeFlybaseTagString (CharPtr locus_tag)
{
  CharPtr new_str;
  CharPtr cp;
  
  if (locus_tag == NULL) return NULL;
  cp = locus_tag;
  while (*cp != 0 && !isdigit(*cp))
  {
    cp++;
  }
  if (StringLen (cp) > 7)
  {
    cp += StringLen (cp) - 7;
  }
  
  new_str = (CharPtr) MemNew (12 * sizeof (Char));
  if (new_str != NULL)
  {
    sprintf (new_str, "FBti0000000");
    StringCpy (new_str + 4 + 7 - StringLen (cp), cp);
  }
  return new_str;
}

static void 
ReplaceRepeatRegionLocusTagWithDbxrefCallback
(SeqFeatPtr sfp,
 Pointer userdata)
{
  GBQualPtr       gbqual, prev_qual = NULL, next_qual;
  ValNodePtr      vnp;
  Boolean         has_dbxref = FALSE;
  DbtagPtr        tag;
  CharPtr         new_string = NULL;
  CharPtr         new_comment;
	SeqFeatXrefPtr 	xrp, prev_xrp = NULL, next_xrp;
	GeneRefPtr      grp;
  
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_repeat_region)
  {
    return;
  }
  
	for (xrp = sfp->xref; xrp; xrp = next_xrp) 
	{
	  next_xrp = xrp->next;
		if (xrp->data.choice == SEQFEAT_GENE) 
		{
			grp = (GeneRefPtr) xrp->data.value.ptrvalue;
			if (grp != NULL && !StringHasNoText (grp->locus_tag))
			{
			  new_string = StringSave (grp->locus_tag);
			  if (prev_xrp == NULL)
			  {
			    sfp->xref = xrp->next;
			  }
			  else
			  {
			    prev_xrp->next = xrp->next;
			  }
			  xrp->next = NULL;
			  SeqFeatXrefFree (xrp);
			}
			else
			{
			  prev_xrp = xrp;
			}
		}
		else
		{
		  prev_xrp = xrp;
		}
	}
  
  for (gbqual = sfp->qual; 
       gbqual != NULL && new_string == NULL;
       gbqual = next_qual)
  {
    next_qual = gbqual->next;
    if (StringCmp (gbqual->qual, "locus_tag") == 0)
    {
      new_string = StringSave (gbqual->val);
      if (prev_qual == NULL)
      {
        sfp->qual = gbqual->next;
      }
      else
      {
        prev_qual->next = gbqual->next;
      }
      gbqual->next = NULL;
      GBQualFree (gbqual);
    }
    else
    {
      prev_qual = gbqual;
    }
  }
  if (new_string == NULL)
  {
    return;
  }
  
  for (vnp = sfp->dbxref; vnp != NULL && !has_dbxref; vnp = vnp->next)
  {
    tag = (DbtagPtr) vnp->data.ptrvalue;
    if (tag != NULL && StringCmp (tag->db, "FLYBASE") == 0)
    {
      has_dbxref = TRUE;
    }
  }
  if (!has_dbxref)
  {
    tag = DbtagNew ();
    if (tag != NULL)
    {
      tag->db = StringSave ("FLYBASE");
      tag->tag = ObjectIdNew ();
      tag->tag->str = MakeFlybaseTagString (new_string);
      ValNodeAddPointer (&(sfp->dbxref), 0, tag);
    }
  }
  
  if (StringHasNoText (sfp->comment))
  {
    sfp->comment = MemFree (sfp->comment);
    sfp->comment = StringSave (new_string);
  }
  else
  {
    new_comment = (CharPtr) MemNew (sizeof (Char) * (StringLen (sfp->comment) + StringLen (new_string) + 2));
    if (new_comment != NULL)
    {
      StringCpy (new_comment, new_string);
      StringCat (new_comment, ";");
      StringCat (new_comment, sfp->comment);
      sfp->comment = MemFree (sfp->comment);
      sfp->comment = new_comment;
    }
  }
  
  MemFree (new_string);
}

extern void ReplaceRepeatRegionLocusTagWithDbxref (IteM i)
{
  BaseFormPtr        bfp;
  SeqEntryPtr       sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitFeaturesInSep (sep, NULL, ReplaceRepeatRegionLocusTagWithDbxrefCallback);
    
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update (); 
}

static void SetPrimerBindPairStrandsCallback (BioseqPtr bsp, Pointer userdata)
{
  SeqFeatPtr        primer_1 = NULL, primer_2 = NULL, sfp;
  SeqMgrFeatContext context;
  Uint1             first_strand, second_strand;
  
  if (bsp == NULL)
  {
    return;
  }
  
  /* must have exactly two primer_bind features */
  primer_1 = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_primer_bind, &context);
  if (primer_1 == NULL)
  {
    return;
  }
  primer_2 = SeqMgrGetNextFeature (bsp, primer_1, 0, FEATDEF_primer_bind, &context);
  if (primer_2 == NULL)
  {
    return;
  }
  
  /* if there are three, must abandon */
  sfp = SeqMgrGetNextFeature (bsp, primer_1, 0, FEATDEF_primer_bind, &context);
  if (sfp != NULL)
  {
    return;
  }
  
  first_strand = SeqLocStrand (primer_1->location);
  second_strand = SeqLocStrand (primer_2->location);
  
  if (first_strand == Seq_strand_minus)
  {
    if (second_strand == Seq_strand_minus)
    {
      SetSeqLocStrand (primer_2->location, Seq_strand_plus);
    }
  }
  else
  {
    if (second_strand != Seq_strand_minus)
    {
      SetSeqLocStrand (primer_2->location, Seq_strand_minus);
    }
  }
}

extern void SetPrimerBindPairStrands (IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitBioseqsInSep (sep, NULL, SetPrimerBindPairStrandsCallback);
    
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update (); 
}

/* This function returns the length change, so that the start
 * position for subsequent features can be adjusted.
 */
static Int4 
FindAndConvertGapFeat 
(BioseqPtr  bsp,
 Int4       start,
 Boolean    make_known)
{
  ValNodePtr vnp;
  SeqLocPtr  slp;
  SeqLitPtr  litp;
  Int4       currpos = 0, len_diff;
  
  if (bsp == NULL || bsp->repr != Seq_repr_delta
      || start < 0)
  {
    return 0;
  }
  for (vnp = (ValNodePtr)(bsp->seq_ext); 
       vnp != NULL && currpos < start; 
       vnp = vnp->next) 
  {
    if (vnp->choice == 1) 
    {
      slp = (SeqLocPtr) vnp->data.ptrvalue;
      if (slp == NULL) continue;
      currpos += SeqLocLen (slp);
    }
    else if (vnp->choice == 2) 
    {
      litp = (SeqLitPtr) vnp->data.ptrvalue;
      if (litp == NULL) continue;
      currpos += litp->length;
    }
  }
  
  if (currpos < start || vnp == NULL || vnp->choice != 2)
  {
    return 0;
  }
  
  litp = (SeqLitPtr) vnp->data.ptrvalue;
  if (litp == NULL)
  {
    return 0;
  }
  
  if (make_known)
  {
    litp->fuzz = IntFuzzFree (litp->fuzz);
    len_diff = 0;
  }
  else if (litp->fuzz == NULL)
  {
    litp->fuzz = IntFuzzNew ();
    litp->fuzz->choice = 4;
    len_diff = litp->length - 100;
    litp->length = 100;
    bsp->length -= len_diff;
  }
  return len_diff;
}

typedef struct gapadjust 
{
  Int4 adjust_start;
  Int4 adjust_len;
} GapAdjustData, PNTR GapAdjustPtr;

static Int4 GetAdjustedGapStart (Int4 feat_left, ValNodePtr prev_adjust)
{
  GapAdjustPtr p;
  Int4         new_feat_left;
  
  new_feat_left = feat_left;
  while (prev_adjust != NULL)
  {
    p = (GapAdjustPtr) prev_adjust->data.ptrvalue;
    if (p != NULL && feat_left > p->adjust_start)
    {
      new_feat_left -= p->adjust_len;
    }
    prev_adjust = prev_adjust->next;
  }
  return new_feat_left;
}

static void ConvertSelectedGapFeatures (IteM i, Boolean to_known)
{
  BaseFormPtr       bfp;
  SelStructPtr      sel;
  SeqMgrFeatContext fcontext;
  BioseqPtr         bsp;
  SeqFeatPtr        sfp;
  Int4              adjusted_start, len_diff;
  ValNodePtr        adjustment_list = NULL;
  GapAdjustPtr      p;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sel = ObjMgrGetSelected ();
  if (sel == NULL)
  {
    Message (MSG_ERROR, "Must select gap features to convert!");
    return;
  }
  WatchCursor ();
  Update ();
  while (sel != NULL)
  {
    if (sel->entityID == bfp->input_entityID
        && sel->itemtype == OBJ_SEQFEAT)
    {
      sfp = SeqMgrGetDesiredFeature (bfp->input_entityID, NULL, sel->itemID, 0, NULL, &fcontext);
      if (sfp != NULL && sfp->idx.subtype == FEATDEF_gap)
      {
        bsp = BioseqFindFromSeqLoc (sfp->location);
        if (bsp != NULL && bsp->repr == Seq_repr_delta)
        {
          adjusted_start = GetAdjustedGapStart (fcontext.left, adjustment_list);
          len_diff = FindAndConvertGapFeat (bsp, adjusted_start, to_known);
          if (len_diff != 0)
          {
            p = (GapAdjustPtr) MemNew (sizeof (GapAdjustData));
            if (p != NULL)
            {
              p->adjust_start = fcontext.left;
              p->adjust_len = len_diff;
              ValNodeAddPointer (&adjustment_list, 0, p);
            }
          }
        }
      }
    }
    sel = sel->next;
  }
  
  adjustment_list = ValNodeFreeData (adjustment_list);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update (); 
}


typedef struct convertgaptounknown
{
  FEATURE_FORM_BLOCK
  
  TexT start_unknown_txt;
  Int4 start_unknown;
} ConvertGapToUnknownData, PNTR ConvertGapToUnknownPtr;


static void FixDeltaFeatures (BioseqPtr bsp, Int4 offset, Int4 len_diff)
{
  SeqFeatPtr             sfp;
  SeqMgrFeatContext      fcontext;
  
  if (len_diff == 0)
  {
    return;
  }
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext))
  {
    AdjustFeatureForGapChange (sfp, bsp, offset, len_diff);
  }  
}

static void ConvertGapFeaturesToUnknownCallback (BioseqPtr bsp, Pointer userdata)
{
  ConvertGapToUnknownPtr cgtup;
  DeltaSeqPtr            dsp;
  SeqLitPtr              slip;
  SeqLocPtr              slp;
  Int4                   len_diff;
  Int4                   offset = 0, add_len;
  
  if (bsp == NULL || bsp->repr != Seq_repr_delta 
      || bsp->seq_ext_type != 4 || bsp->seq_ext == NULL
      || userdata == NULL)
  {
    return;
  }
  
  cgtup = (ConvertGapToUnknownPtr) userdata;
  
  dsp = (DeltaSeqPtr) bsp->seq_ext;
  while (dsp != NULL)
  {
    add_len = 0;
    if (dsp->choice == 1 && dsp->data.ptrvalue != NULL)
    {
      slp = (SeqLocPtr) dsp->data.ptrvalue;
      add_len = SeqLocLen (slp);
    }
    else if (dsp->choice == 2 && dsp->data.ptrvalue != NULL)
    {
      slip = (SeqLitPtr) (dsp->data.ptrvalue);
      add_len = slip->length;
      if (slip->seq_data == NULL
          && slip->length >= cgtup->start_unknown)
      {
        len_diff = slip->length - 100;
        slip->length = 100;
        if (slip->fuzz != NULL)
        {
          slip->fuzz = IntFuzzFree (slip->fuzz);
        }
        slip->fuzz = IntFuzzNew();
        slip->fuzz->choice = 4;
        
        if (len_diff > 0) {
            FixDeltaFeatures (bsp, offset, len_diff);
        }
        add_len -= len_diff;      
      }
    }
    offset += add_len;
    dsp = dsp->next;
  }
  bsp->length = offset;
}

static void ConvertGapFeaturesToUnknownButton (ButtoN b)
{
  ConvertGapToUnknownPtr cgtup;
  SeqEntryPtr            sep;
  CharPtr                str;
  
  cgtup = (ConvertGapToUnknownPtr) GetObjectExtra (b);
  if (cgtup == NULL)
  {
    return;
  }
  
  str = SaveStringFromText (cgtup->start_unknown_txt);
  if (StringHasNoText (str))
  {
    str = MemFree (str);
    return;
  }
  
  cgtup->start_unknown = atoi (str);
  str = MemFree (str);
  if (cgtup->start_unknown <= 0)
  {
    return;
  }
  
  WatchCursor();
  Update();
  
  sep = GetTopSeqEntryForEntityID (cgtup->input_entityID);
  
  VisitBioseqsInSep (sep, cgtup, ConvertGapFeaturesToUnknownCallback);
  
  ObjMgrSetDirtyFlag (cgtup->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, cgtup->input_entityID, 0, 0);
  Remove (cgtup->form);
  ArrowCursor ();
  Update (); 
}

extern void ConvertGapFeaturesToUnknown (IteM i)
{
  BaseFormPtr            bfp;
  ConvertGapToUnknownPtr cgtup;
  WindoW                 w;
  GrouP                  h, g, c;
  ButtoN                 b;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  cgtup = (ConvertGapToUnknownPtr) MemNew (sizeof (ConvertGapToUnknownData));
  if (cgtup == NULL) return;
  w = FixedWindow (-50, -33, -10, -10, "Convert Known Length Gaps to Unknown", StdCloseWindowProc);
  SetObjectExtra (w, cgtup, StdCleanupFormProc);
  cgtup->form = (ForM) w;
  cgtup->input_entityID = bfp->input_entityID;
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  g = HiddenGroup (h, 2, 0, NULL);
  StaticPrompt (g, "Convert gaps longer or equal to", 0, popupMenuHeight, programFont, 'r');
  cgtup->start_unknown_txt = DialogText (g, "100", 10, NULL);
  
  c = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  b = PushButton (c, "Accept", ConvertGapFeaturesToUnknownButton);
  SetObjectExtra (b, cgtup, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
  
  Show (w);
}


typedef struct changegaplen
{
  FEATURE_FORM_BLOCK
  
  TexT   length_txt;
  Int4   length;
} ChangeGapLenData, PNTR ChangeGapLenPtr;


static void ChangeOneGapLength (SeqFeatPtr sfp, Int4 new_length)
{
  BioseqPtr   bsp;
  DeltaSeqPtr dsp;
  SeqLitPtr              slip;
  SeqLocPtr              slp;
  Int4                   len_diff;
  Int4                   offset = 0, gap_start;
  
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_gap || new_length < 0
      || new_length == SeqLocLen (sfp->location))
  {
    return;
  }
  
  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL || bsp->repr != Seq_repr_delta 
      || bsp->seq_ext_type != 4 || bsp->seq_ext == NULL)
  {
    return;
  }
  
  gap_start = SeqLocStart (sfp->location);
  
  dsp = (DeltaSeqPtr) bsp->seq_ext;
  while (dsp != NULL && offset < gap_start)
  {
    if (dsp->choice == 1 && dsp->data.ptrvalue != NULL)
    {
      slp = (SeqLocPtr) dsp->data.ptrvalue;
      offset += SeqLocLen (slp);
    }
    else if (dsp->choice == 2 && dsp->data.ptrvalue != NULL)
    {
      slip = (SeqLitPtr) (dsp->data.ptrvalue);
      offset += slip->length;
    }
    dsp = dsp->next;
  }
  
  if (offset == gap_start && dsp != NULL && dsp->choice == 2 && dsp->data.ptrvalue != NULL)
  {
    slip = (SeqLitPtr) (dsp->data.ptrvalue);
    if (slip->seq_data == NULL && slip->fuzz == NULL)
    {
      len_diff = slip->length - new_length;
      slip->length = new_length;
      FixDeltaFeatures (bsp, offset, len_diff);
    }
  }
  
  
  bsp->length -= len_diff;
}


static void ChangeGapLength (ButtoN b)
{
  ChangeGapLenPtr   cglp;
  SelStructPtr      sel;
  SeqFeatPtr        sfp;
  SeqMgrFeatContext fcontext;
  CharPtr           str;

  cglp = (ChangeGapLenPtr) GetObjectExtra (b);
  if (cglp == NULL)
  {
    return;
  }
  
  cglp->length = 0;
  str = SaveStringFromText (cglp->length_txt);
  if (!StringHasNoText (str))
  {
    cglp->length = atoi (str);
  }
  str = MemFree (str);
  if (cglp->length < 1)
  {
    Message (MSG_ERROR, "Must select a gap size greater than zero!");
    return;
  }
  
  sel = ObjMgrGetSelected ();
  if (sel == NULL)
  {
    Message (MSG_ERROR, "No gaps selected!");
    return;
  }
  WatchCursor ();
  Update ();
  while (sel != NULL)
  {
    if (sel->entityID == cglp->input_entityID
        && sel->itemtype == OBJ_SEQFEAT)
    {
      sfp = SeqMgrGetDesiredFeature (cglp->input_entityID, NULL, sel->itemID, 0, NULL, &fcontext);
      if (sfp != NULL && sfp->idx.subtype == FEATDEF_gap)
      {
        ChangeOneGapLength (sfp, cglp->length);
      }
    }
    sel = sel->next;
  }
  
  ObjMgrSetDirtyFlag (cglp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, cglp->input_entityID, 0, 0);
  ArrowCursor ();
  Update (); 
  
  if (!GetStatus (cglp->leave_dlg_up))
  {
    Remove (cglp->form);
  }
}

extern void ChangeKnownGapLength (IteM i)
{
  BaseFormPtr     bfp;
  ChangeGapLenPtr cglp;
  WindoW          w;
  GrouP           h, g, c;
  ButtoN          b;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  cglp = (ChangeGapLenPtr) MemNew (sizeof (ChangeGapLenData));
  if (cglp == NULL) return;
  w = FixedWindow (-50, -33, -10, -10, "Change Length of Selected Known Length Gaps", StdCloseWindowProc);
  SetObjectExtra (w, cglp, StdCleanupFormProc);
  cglp->form = (ForM) w;
  cglp->input_entityID = bfp->input_entityID;
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  g = HiddenGroup (h, 2, 0, NULL);
  StaticPrompt (g, "Change length to", 0, popupMenuHeight, programFont, 'r');
  cglp->length_txt = DialogText (g, "100", 10, NULL);
  
  c = HiddenGroup (h, 3, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  b = PushButton (c, "Accept", ChangeGapLength);
  SetObjectExtra (b, cglp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  cglp->leave_dlg_up = CheckBox (c, "Leave Dialog Up", NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
  
  Show (w);
}

extern void ConvertSelectedGapFeaturesToKnown (IteM i)
{
  ConvertSelectedGapFeatures (i, TRUE);  
}

extern void ConvertSelectedGapFeaturesToUnknown (IteM i)
{
  ConvertSelectedGapFeatures (i, FALSE);
}


static Int4 CountNsAtEndOfSeqLit (SeqLitPtr slip, Uint2 which_end)
{
  Int2         residue;
  Int4         change_len = 0;
  
  if (slip == NULL || slip->seq_data == NULL)
  {
    return 0;
  }
  
  if (slip->seq_data_type != Seq_code_iupacna)
  {
    slip->seq_data = BSConvertSeq(slip->seq_data, Seq_code_iupacna, 
                          slip->seq_data_type, 
                          slip->length);
    slip->seq_data_type = Seq_code_iupacna;
  }
  
  if (which_end == SEQLOC_LEFT_END)
  {
    BSSeek(slip->seq_data, 0, SEEK_SET);
  }
  else
  {
    BSSeek(slip->seq_data, slip->length - 1, SEEK_SET);    
  }
  
  residue = BSGetByte(slip->seq_data);
  while (change_len < slip->length && residue == 'N')
  {
    change_len++;
    if (which_end != SEQLOC_LEFT_END)
    {
      BSSeek (slip->seq_data, slip->length - change_len - 1, SEEK_SET);
    }
    residue = BSGetByte (slip->seq_data);
  }

  return change_len;
}

static void RemoveSeqLitEnd (SeqLitPtr slip, Int4 change_len, Uint2 which_end)
{
  ByteStorePtr bs_new;

  if (slip == NULL || change_len < 1)
  {
    return;
  }
  if (slip->seq_data_type != Seq_code_iupacna)
  {
    slip->seq_data = BSConvertSeq(slip->seq_data, Seq_code_iupacna, 
                                  slip->seq_data_type, 
                                  slip->length);
    slip->seq_data_type = Seq_code_iupacna;
  }
  bs_new = BSNew (slip->length - change_len);
  if (which_end == SEQLOC_LEFT_END)
  {
    BSSeek (slip->seq_data, change_len, SEEK_SET);
  }
  else
  {
    BSSeek(slip->seq_data, 0, SEEK_SET);
  }
  BSInsertFromBS (bs_new, slip->seq_data, slip->length - change_len);
  slip->seq_data = BSFree (slip->seq_data);
  slip->seq_data = bs_new;
  slip->length -= change_len;
}

static void ExpandGapsToIncludeFlankingNs (BioseqPtr bsp, Pointer userdata)
{
  DeltaSeqPtr  dsp;
  Int4         offset = 0, last_length = 0, change_len, add_offset;
  SeqLocPtr    slp;
  SeqLitPtr    slip, prev_slip = NULL, next_slip;
  
  if (bsp == NULL || bsp->repr != Seq_repr_delta 
      || bsp->seq_ext_type != 4 || bsp->seq_ext == NULL)
  {
    return;
  }
  
  for (dsp = bsp->seq_ext; dsp != NULL; dsp = dsp->next)
  {
    add_offset = 0;
    if (dsp->choice == 1)
    {
      slp = (SeqLocPtr) dsp->data.ptrvalue;
      add_offset = SeqLocLen (slp);
      last_length = 0;
      prev_slip = NULL;
    }
    else if (dsp->choice == 2)
    {
      slip = (SeqLitPtr) dsp->data.ptrvalue;
      add_offset = slip->length;
      if (slip->seq_data == NULL)
      {
        if (slip->fuzz == NULL)
        {
          /* check for Ns before gap of known length */
          if (prev_slip != NULL && last_length > 0)
          {
            change_len = CountNsAtEndOfSeqLit (prev_slip, SEQLOC_RIGHT_END);

            if (change_len > 0)
            {
              RemoveSeqLitEnd (prev_slip, change_len, SEQLOC_RIGHT_END);
              slip->length += change_len;
            }
          }
          /* check for Ns after gap of known length */
          if (dsp->next != NULL && dsp->next->choice == 2 && dsp->next->data.ptrvalue != NULL)
          {
            next_slip = (SeqLitPtr) dsp->next->data.ptrvalue;
            change_len = CountNsAtEndOfSeqLit (next_slip, SEQLOC_LEFT_END);
            
            if (change_len > 0)
            {
              RemoveSeqLitEnd (next_slip, change_len, SEQLOC_LEFT_END);
              slip->length += change_len;
              add_offset += change_len;
            }
          }
        }
        prev_slip = NULL;
        last_length = 0;
      }
      else
      {
        prev_slip = slip;
        last_length = slip->length;
      }
    }
    else
    {
      prev_slip = NULL;
      last_length = 0;
    }
    offset += add_offset;
  }
  BioseqPack (bsp);
}

extern void AddFlankingNsToKnownLengthGaps (IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  
  WatchCursor ();
  Update ();
  VisitBioseqsInSep (sep, NULL, ExpandGapsToIncludeFlankingNs);
  
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update (); 
}


static void CombineAdjacentGapsCallback (BioseqPtr bsp, Pointer userdata)
{
  ValNodePtr vnp, prev_vnp = NULL, next_vnp;
  SeqLitPtr  litp, prev_litp = NULL;
  SeqLocPtr  slp;
  Int4       offset = 0;
  Int4       len_diff;
  SeqFeatPtr sfp;
  SeqMgrFeatContext fcontext;
  
  if (bsp == NULL || bsp->repr != Seq_repr_delta)
  {
    return;
  }
  
  for (vnp = (ValNodePtr)(bsp->seq_ext); 
       vnp != NULL; 
       vnp = next_vnp) 
  {
    next_vnp = vnp->next;
    if (vnp->choice == 1) 
    {
      prev_vnp = NULL;
      prev_litp = NULL;
      slp = (SeqLocPtr) vnp->data.ptrvalue;
      offset += SeqLocLen (slp);
    }
    else if (vnp->choice == 2) 
    {
      litp = (SeqLitPtr) vnp->data.ptrvalue;
      if (litp == NULL || litp->seq_data != NULL)
      {
        prev_vnp = NULL;
        prev_litp = NULL;
        offset += litp->length;
      }
      else if (prev_vnp != NULL && prev_litp != NULL)
      {
        /* combine two adjacent gaps */
        len_diff = 0;
        offset = offset - prev_litp->length;
        if (prev_litp->fuzz != NULL || litp->fuzz != NULL)
        {
          /* at least one is unknown length, so result will be unknown */
          len_diff = prev_litp->length + litp->length - 100;
          
          prev_litp->length = 100;
          if (prev_litp->fuzz == NULL)
          {
            prev_litp->fuzz = IntFuzzNew ();
            prev_litp->fuzz->choice = 4;
          }          
        }
        else
        {
          prev_litp->length += litp->length;
        }
        prev_vnp->next = vnp->next;
        vnp->next = NULL;
        litp = SeqLitFree (litp);
        vnp = ValNodeFree (vnp);
        if (len_diff != 0)
        {
          FixDeltaFeatures (bsp, offset, len_diff);
        }
        offset += prev_litp->length;
      }
      else
      {
        prev_vnp = vnp;
        prev_litp = litp;
        offset += litp->length;
      }
    }
  }
  bsp->length = offset;
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, FEATDEF_CDS, &fcontext);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, FEATDEF_CDS, &fcontext))
  {
    AdjustCDSLocationsForGapsCallback (sfp, NULL);
  }
}

extern void CombineAdjacentGaps (IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  
  WatchCursor ();
  Update ();
  VisitBioseqsInSep (sep, NULL, CombineAdjacentGapsCallback);
  
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update (); 
}

static void MarkPseudoGenesCallback (SeqFeatPtr sfp, Pointer userdata)
{
  SeqFeatPtr gene;
  
  if (sfp == NULL || ! sfp->pseudo 
      || sfp->data.choice == SEQFEAT_GENE
      || SeqMgrGetGeneXref (sfp) != NULL)
  {
    return;
  }
  
  gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
  if (gene != NULL) 
  {
    gene->pseudo = TRUE;
  }
}

extern void MarkGenesWithPseudoFeaturesPseudo (IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  
  WatchCursor ();
  Update ();
  VisitFeaturesInSep (sep, NULL, MarkPseudoGenesCallback);
  
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update (); 
  
}

static void RemoveOneNomenclature (UserObjectPtr PNTR puop)
{
  UserObjectPtr  uop, obj, prev_obj, next_obj;
  UserFieldPtr   prev_ufp, next_ufp, ufp;
  ObjectIdPtr    oip;

  if (puop == NULL || *puop == NULL) return;
  uop = *puop;
  
  for (ufp = uop->data, prev_ufp = NULL; 
       ufp != NULL; 
       ufp = next_ufp) {
    next_ufp = ufp->next;
    if (ufp->choice == 6) {
      obj = (UserObjectPtr) ufp->data.ptrvalue;
      RemoveOneNomenclature (&obj);
      ufp->data.ptrvalue = obj;
    } else if (ufp->choice == 12) {
      for (obj = (UserObjectPtr) ufp->data.ptrvalue, prev_obj = NULL;
           obj != NULL;
           obj = next_obj) {
        next_obj = obj->next;
        RemoveOneNomenclature (&obj);
        if (obj == NULL)
        {
          if (prev_obj == NULL)
          {
            ufp->data.ptrvalue = next_obj;
          }
          else
          {
            prev_obj->next = next_obj;
          }
          obj = UserObjectFree (obj);
        }
        else
        {
          prev_obj = obj;
        }
      }
    }
    if ((ufp->choice == 6 || ufp->choice == 12) && ufp->data.ptrvalue == NULL)
    {
      if (prev_ufp == NULL)
      {
        uop->data = ufp->next;
      }
      else
      {
        prev_ufp->next = ufp->next;
      }
      ufp = UserFieldFree (ufp);
    }
    else
    {
      prev_ufp = ufp;
    }
  }
  
  oip = uop->type;
  if (oip != NULL && StringCmp (oip->str, "OfficialNomenclature") == 0)
  {
    uop = UserObjectFree (uop);
    *puop = uop;
  }
}

static void RemoveNomenclatureCallback (SeqFeatPtr sfp, Pointer userdata)
{
  UserObjectPtr uop;
  
  if (sfp != NULL)
  {
    uop = sfp->ext;
    RemoveOneNomenclature (&uop);
    sfp->ext = uop;
  }
}

extern void RemoveNomenclature (IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  
  WatchCursor ();
  Update ();
  
  VisitFeaturesInSep (sep, NULL, RemoveNomenclatureCallback);
  
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update (); 
}


static void RemoveUnindexedFeaturesInSeqEntry (SeqEntryPtr sep, Uint2 entityID)
{
  BioseqPtr         bsp;
  BioseqSetPtr      bssp;
  SeqAnnotPtr       sap = NULL;
  SeqFeatPtr        sfp;
  SeqMgrFeatContext context;

  if (sep == NULL || sep->data.ptrvalue == NULL)
  {
	return;
  }

  if (IS_Bioseq(sep))
  {
	bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
  }
  else if (IS_Bioseq_set (sep))
  {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;

	sep = bssp->seq_set;
	while (sep != NULL)
	{
      RemoveUnindexedFeaturesInSeqEntry (sep, entityID);
	  sep = sep->next;
	}
    sap = bssp->annot;
  }

  while (sap != NULL)
  {
    if (sap->type == 1)
	{
      sfp = (SeqFeatPtr) sap->data;
	  while (sfp != NULL)
	  {
	    if (SeqMgrGetDesiredFeature (entityID, NULL, 0, 0, sfp, &context) == NULL)
		{
		  sfp->idx.deleteme = TRUE;
		}
		sfp = sfp->next;
	  }
	}
	sap = sap->next;
  }
}


extern void RemoveUnindexedFeatures (IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  
  WatchCursor ();
  Update ();
  
  RemoveUnindexedFeaturesInSeqEntry (sep, bfp->input_entityID);

  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update (); 
}


static void CopyLocusToLocusTagCallback (SeqFeatPtr sfp, Pointer userdata)
{
  GeneRefPtr grp;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_GENE || sfp->data.value.ptrvalue == NULL)
  {
    return;
  }
  
  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (!StringHasNoText (grp->locus) && StringHasNoText (grp->locus_tag))
  {
    grp->locus_tag = MemFree (grp->locus_tag);
    grp->locus_tag = StringSave (grp->locus);
  }
}

extern void CopyLocusToLocusTag (IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  
  WatchCursor ();
  Update ();
  
  VisitFeaturesInSep (sep, NULL, CopyLocusToLocusTagCallback);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update (); 
}


/* Discrepancy Report */

typedef void (*DiscrepancyCallback) (ValNodePtr item_list, Pointer userdata);
typedef void (*DiscrepancyCallbackDataFree) (Pointer userdata);

/* To add a new type of test, do ALL Of the following:
 * 1. add an item to the DiscrepancyType enum (this will fill the disc_type value)
 * 2. add a collection function and declare it with the others
 * 3. add an item to discrepancy_info_list that corresponds with the position of the
 *    new enum value.  If you are combining multiple types in one collection function,
 *    be sure to list them together.
 */

typedef enum {
  DISC_GENE_MISSING = 0,
  DISC_SUPERFLUOUS_GENE,
  DISC_GENE_MISSING_LOCUS_TAG,
  DISC_GENE_DUPLICATE_LOCUS_TAG,
  DISC_GENE_LOCUS_TAG_BAD_FORMAT,
  DISC_GENE_LOCUS_TAG_INCONSISTENT_PREFIX,
  DISC_NON_GENE_LOCUS_TAG,
  DISC_MISSING_PROTEIN_ID,
  DISC_INCONSISTENT_PROTEIN_ID_PREFIX,
  DISC_GENE_CDS_mRNA_LOCATION_CONFLICT,
  DISC_GENE_PRODUCT_CONFLICT,
  DISC_GENE_DUPLICATE_LOCUS,
  DISC_EC_NUMBER_NOTE,
  DISC_PSEUDO_MISMATCH,
  DISC_JOINED_FEATURES,
  DISC_OVERLAPPING_GENES,
  DISC_OVERLAPPING_CDS,
  DISC_SHORT_CONTIG,
  DISC_INCONSISTENT_BIOSRC,
  DISC_SUSPECT_PRODUCT_NAME,
  DISC_INCONSISTENT_BIOSRC_DEFLINE,
  MAX_DISC_TYPE
} DiscrepancyType;

typedef void (*PerformDiscrepancyTest) PROTO ((SeqEntryPtr, ValNodePtr PNTR));

typedef struct discrepancyinfo 
{
  CharPtr                conf_name;
  CharPtr                setting_name;
  PerformDiscrepancyTest test_func;
} DiscrepancyInfoData, PNTR DiscrepancyInfoPtr;

/* declarations for discrepancy tests */
static void AddMissingAndSuperfluousGeneDiscrepancies (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list);
static void AddDiscrepanciesForMissingOrNonUniqueGeneLocusTags (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list);
static void AddDiscrepanciesForNonGeneLocusTags (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list);
static void FindMissingProteinIDs (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list);
static void FindCDSmRNAGeneLocationDiscrepancies (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list);
static void FindCDSGeneProductConflicts (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list);
static void FindDuplicateGeneLocus (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list);
static void AddECNumberNoteDiscrepancies (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list);
static void FindPseudoDiscrepancies (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list);
static void AddJoinedFeatureDiscrepancies (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list);
static void AddOverlappingGeneDiscrepancies (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list);
static void AddOverlappingCodingRegionDiscrepancies (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list);
static void FindShortContigs (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list);
static void FindNonmatchingContigSources (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list);
static void FindSuspectProductNames (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list);
static void FindInconsistentSourceAndDefline (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list);

const DiscrepancyInfoData discrepancy_info_list[] = 
{
  { "Missing Genes", "MISSING_GENES", AddMissingAndSuperfluousGeneDiscrepancies },
  { "Extra Genes", "EXTRA_GENES", AddDiscrepanciesForMissingOrNonUniqueGeneLocusTags },
  { "Missing Locus Tags", "MISSING_LOCUS_TAGS", AddDiscrepanciesForMissingOrNonUniqueGeneLocusTags },
  { "Duplicate Locus Tags", "DUPLICATE_LOCUS_TAGS", AddDiscrepanciesForMissingOrNonUniqueGeneLocusTags },
  { "Bad Locus Tag Format", "BAD_LOCUS_TAG_FORMAT", AddDiscrepanciesForMissingOrNonUniqueGeneLocusTags },
  { "Inconsistent Locus Tag Prefix", "INCONSISTENT_LOCUS_TAG_PREFIX", AddDiscrepanciesForMissingOrNonUniqueGeneLocusTags },
  { "Nongene Locus Tag", "NON_GENE_LOCUS_TAG", AddDiscrepanciesForNonGeneLocusTags },
  { "Missing Protein ID", "MISSING_PROTEIN_ID", FindMissingProteinIDs },
  { "Inconsistent Protein ID", "INCONSISTENT_PROTEIN_ID", FindMissingProteinIDs },
  { "Feature Location Conflict", "FEATURE_LOCATION_CONFLICT", FindCDSmRNAGeneLocationDiscrepancies },
  { "Gene Product Conflict", "GENE_PRODUCT_CONFLICT", FindCDSGeneProductConflicts },
  { "Duplicate Gene Locus", "DUPLICATE_GENE_LOCUS", FindDuplicateGeneLocus },
  { "EC Number Note", "EC_NUMBER_NOTE", AddECNumberNoteDiscrepancies },
  { "Pseudo Mismatch", "PSEUDO_MISMATCH", FindPseudoDiscrepancies },
  { "Joined Features", "JOINED_FEATURES", AddJoinedFeatureDiscrepancies },
  { "Overlapping Genes", "OVERLAPPING_GENES", AddOverlappingGeneDiscrepancies },
  { "Overlapping CDS", "OVERLAPPING_CDS", AddOverlappingCodingRegionDiscrepancies },
  { "Short Contig", "SHORT_CONTIG", FindShortContigs },
  { "Inconsistent BioSource", "INCONSISTENT_BIOSOURCE", FindNonmatchingContigSources },
  { "Suspect Product Name", "SUSPECT_PRODUCT_NAMES", FindSuspectProductNames },
  { "Inconsistent Source And Definition Line", "INCONSISTENT_SOURCE_DEFLINE", FindInconsistentSourceAndDefline }
    
};


typedef struct discrepancyitem 
{
  DiscrepancyType             disc_type;
  CharPtr                     description;
  ValNodePtr                  item_list;
  DiscrepancyCallback         callback_func; 
  DiscrepancyCallbackDataFree datafree_func; 
  Pointer                     callback_data;
  Boolean                     chosen;
  ValNodePtr                  subcategories;
  Boolean                     expanded;
  Int4                        level;
} DiscrepancyItemData, PNTR DiscrepancyItemPtr;

typedef struct discrepancyconfig
{
  Boolean conf_list[MAX_DISC_TYPE];
  Boolean use_feature_table_format;
} DiscrepancyConfigData, PNTR DiscrepancyConfigPtr;

typedef struct discrepancyreportform 
{
  FORM_MESSAGE_BLOCK
  SeqEntryPtr     sep;
  ValNodePtr      discrepancy_list;
  DoC             doc;
  DoC             item_list;
  ButtoN          recheck_btn;
  Int2            selected;
  Int2            item_selected;
  Int2            clicked;
  Boolean         dblClick;
  BaseFormPtr     bfp;
  
  Int4            num_levels;
  Nlm_ColPtr PNTR col_fmt_array_array;
  DiscrepancyConfigPtr dcp;
} DiscrepancyReportFormData, PNTR DiscrepancyReportFormPtr;


static DiscrepancyConfigPtr DiscrepancyConfigFree (DiscrepancyConfigPtr dcp)
{
  return MemFree (dcp);  
}

static DiscrepancyConfigPtr DiscrepancyConfigNew (void)
{
  DiscrepancyConfigPtr dcp;
  Int4                 i;
  
  dcp = (DiscrepancyConfigPtr) MemNew (sizeof (DiscrepancyConfigData));
  for (i = 0; i < MAX_DISC_TYPE; i++)
  {
    dcp->conf_list[i] = TRUE;
  }
  dcp->use_feature_table_format = FALSE;
  return dcp;
}

static DiscrepancyConfigPtr ReadDiscrepancyConfig (void)
{
  DiscrepancyConfigPtr dcp;
  Int4                 i;
  Char                 str[20];
  
  dcp = DiscrepancyConfigNew();
  if (dcp != NULL)
  {
    for (i = 0; i < MAX_DISC_TYPE; i++)
    {
      if (GetSequinAppParam ("DISCREPANCY_REPORT", discrepancy_info_list[i].setting_name, NULL, str, sizeof (str))) {
        if (StringICmp (str, "FALSE") == 0) {
          dcp->conf_list[i] = FALSE;
        }
      }
    }
    if (GetSequinAppParam ("DISCREPANCY_REPORT", "USE_FEATURE_TABLE_FORMAT", NULL, str, sizeof (str))) {
      if (StringICmp (str, "TRUE") == 0) {
        dcp->use_feature_table_format = TRUE;
      }
    }
  }
  return dcp;
}

static void SaveDiscrepancyConfig (DiscrepancyConfigPtr dcp)
{
  Int4 i;
  
  if (dcp == NULL)
  {
    return;
  }
  
  for (i = 0; i < MAX_DISC_TYPE; i++)
  {
    if (dcp->conf_list[i])
    {
      SetAppParam ("SEQUINCUSTOM", "DISCREPANCY_REPORT", discrepancy_info_list[i].setting_name, "TRUE");
    }
    else
    {
      SetAppParam ("SEQUINCUSTOM", "DISCREPANCY_REPORT", discrepancy_info_list[i].setting_name, "FALSE");
    }
  }
  if (dcp->use_feature_table_format)
  {
    SetAppParam ("SEQUINCUSTOM", "DISCREPANCY_REPORT", "USE_FEATURE_TABLE_FORMAT", "TRUE");
  }
  else
  {
    SetAppParam ("SEQUINCUSTOM", "DISCREPANCY_REPORT", "USE_FEATURE_TABLE_FORMAT", "FALSE");
  }
}


static void UnselectDiscrepancyList(ButtoN b)
{
  ButtoN *test_options;
  Int4    i;
  
  test_options = (ButtoN *) GetObjectExtra (b);
  if (test_options != NULL)
  {
    for (i = 0; i < MAX_DISC_TYPE; i++)
    {
      SetStatus (test_options[i], FALSE);
    }
  }
}

static void SelectDiscrepancyList(ButtoN b)
{
  ButtoN *test_options;
  Int4    i;
  
  test_options = (ButtoN *) GetObjectExtra (b);
  if (test_options != NULL)
  {
    for (i = 0; i < MAX_DISC_TYPE; i++)
    {
      SetStatus (test_options[i], TRUE);
    }
  }
}

/* This function returns TRUE if there was a change to the discrepancy config,
 * FALSE otherwise.
 */
static Boolean EditDiscrepancyConfig (DiscrepancyConfigPtr dcp)
{
  WindoW                w;
  GrouP                 h, g, k, c;
  ButtoN                b, use_feature_table_format_btn;
  ModalAcceptCancelData acd;
  Int4                  i;
  ButtoN                test_options[MAX_DISC_TYPE];
  Boolean               rval = FALSE;
  
  if (dcp == NULL)
  {
    return rval;
  }
  
  acd.accepted = FALSE;
  acd.cancelled = FALSE;
  
  w = ModalWindow(-20, -13, -10, -10, NULL);
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  g = NormalGroup (h, 0, 10, "Discrepancy Tests to Run", programFont, NULL);
  SetGroupSpacing (g, 10, 10);
  for (i = 0; i < MAX_DISC_TYPE; i++)
  {
    test_options[i] = CheckBox (g, discrepancy_info_list[i].conf_name, NULL);
    SetStatus (test_options[i], dcp->conf_list[i]);
  }
  
  use_feature_table_format_btn = CheckBox (h, "Use feature table format for features in report", NULL);
  SetStatus (use_feature_table_format_btn, dcp->use_feature_table_format);
  
  k = HiddenGroup (h, 2, 0, NULL);
  b = PushButton (k, "Select All", SelectDiscrepancyList);
  SetObjectExtra (b, test_options, NULL);
  b = PushButton (k, "Unselect All", UnselectDiscrepancyList);
  SetObjectExtra (b, test_options, NULL);  
  
  c = HiddenGroup (h, 3, 0, NULL);
  b = PushButton (c, "Accept", ModalAcceptButton);
  SetObjectExtra (b, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) use_feature_table_format_btn, (HANDLE) k, (HANDLE) c, NULL);
  
  Show(w); 
  Select (w);
  while (!acd.accepted && ! acd.cancelled)
  {
    ProcessExternalEvent ();
    Update ();
  }
  ProcessAnEvent ();
  if (acd.accepted)
  {
    for (i = 0; i < MAX_DISC_TYPE; i++)
    {
      dcp->conf_list [i] = GetStatus (test_options[i]);
    }
    dcp->use_feature_table_format = GetStatus (use_feature_table_format_btn);
    rval = TRUE;
    SaveDiscrepancyConfig (dcp);
  }

  Remove (w);
  return rval;
}


static WindoW discrepancyReportWindow = NULL;


static ValNodePtr FreeDiscrepancyList (ValNodePtr list)
{
  DiscrepancyItemPtr dip;
  
  if (list == NULL)
  {
    return list;
  }
  FreeDiscrepancyList (list->next);
  list->next = NULL;
  dip = (DiscrepancyItemPtr) list->data.ptrvalue;
  if (dip != NULL)
  {
    dip->description = MemFree (dip->description);
    if (dip->datafree_func != NULL)
    {
      (dip->datafree_func) (dip->callback_data);
    }
    dip->item_list = ValNodeFree (dip->item_list);
    
    dip->subcategories = FreeDiscrepancyList (dip->subcategories);
    dip = MemFree (dip);
    list->data.ptrvalue = NULL;
  }
  list = ValNodeFree (list);
  return list;
}


static DiscrepancyItemPtr 
NewDiscrepancyItem 
(DiscrepancyType disc_type,
 CharPtr         description_fmt,
 ValNodePtr      item_list)
{
  DiscrepancyItemPtr dip;
  
  dip = (DiscrepancyItemPtr) MemNew (sizeof (DiscrepancyItemData));
  if (dip != NULL)
  {
    dip->disc_type = disc_type;
    dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (description_fmt) + 15));
    sprintf (dip->description, description_fmt, ValNodeLen (item_list));
    dip->callback_func = NULL;
    dip->datafree_func = NULL;
    dip->callback_data = NULL;
    dip->item_list = item_list;
    dip->subcategories = NULL;
    dip->expanded = FALSE;
    dip->level = 0;
  }
  return dip;  
}


static void 
RemoveUnwantedDiscrepancyItems 
(ValNodePtr PNTR      discrepancy_list,
 DiscrepancyConfigPtr dcp)
{
  ValNodePtr         vnp, prev = NULL, vnp_next;
  DiscrepancyItemPtr dip;
  
  if (dcp == NULL || discrepancy_list == NULL || *discrepancy_list == NULL)
  {
    return;
  }
  
  for (vnp = *discrepancy_list; vnp != NULL; vnp = vnp_next)
  {
    vnp_next = vnp->next;
    dip = (DiscrepancyItemPtr) vnp->data.ptrvalue;
    if (dip == NULL || ! dcp->conf_list[dip->disc_type])
    {
      if (prev == NULL)
      {
        *discrepancy_list = vnp_next;
      }
      else
      {
        prev->next = vnp_next;
      }
      vnp->next = NULL;
      vnp = FreeDiscrepancyList (vnp);
    }
    else
    {
      prev = vnp;
    }
  }
  
}


static void ValNodeLinkCopy (ValNodePtr PNTR list1, ValNodePtr list2)
{
  if (list1 == NULL) return;
  while (list2 != NULL)
  {
    ValNodeAddPointer (list1, list2->choice, list2->data.ptrvalue);
    list2 = list2->next;
  }
}


static Boolean ValNodeStringListMatch (ValNodePtr vnp1, ValNodePtr vnp2)
{
  if (vnp1 == NULL && vnp2 == NULL)
  {
    return TRUE;
  }
  else if (vnp1 == NULL || vnp2 == NULL)
  {
    return FALSE;
  }
  else if (StringCmp (vnp1->data.ptrvalue, vnp2->data.ptrvalue) != 0)
  {
    return FALSE;
  }
  else
  {
    return ValNodeStringListMatch (vnp1->next, vnp2->next);
  }
}


static Boolean GeneRefMatch (GeneRefPtr grp1, GeneRefPtr grp2)
{
  if (grp1 == NULL && grp2 == NULL)
  {
    return TRUE;
  }
  else if (grp1 == NULL || grp2 == NULL)
  {
    return FALSE;
  }
  else if (StringCmp (grp1->locus, grp2->locus) != 0
           || StringCmp (grp1->allele, grp2->allele) != 0
           || StringCmp (grp1->desc, grp2->desc) != 0
           || StringCmp (grp1->maploc, grp2->maploc) != 0
           || StringCmp (grp1->locus_tag, grp2->locus_tag) != 0
           || (grp1->pseudo && !grp2->pseudo)
           || (!grp1->pseudo && grp2->pseudo)
           || !ValNodeStringListMatch (grp1->db, grp2->db)
           || !ValNodeStringListMatch (grp1->syn, grp2->syn))
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}


static void ExtractGeneFromListByGeneRef (ValNodePtr PNTR list, GeneRefPtr grp)
{
  ValNodePtr prev = NULL, this_vnp, next_vnp;
  
  if (list == NULL || grp == NULL)
  {
    return;
  }
  
  this_vnp = *list;
  while (this_vnp != NULL)
  {
    next_vnp = this_vnp->next;
    if (GeneRefMatch (this_vnp->data.ptrvalue, grp))
    {
      if (prev == NULL)
      {
        *list = next_vnp;
      }
      else
      {
        prev->next = next_vnp;
      }
      this_vnp->next = NULL;
      ValNodeFree (this_vnp);
    }
    else
    {
      prev = this_vnp;
    }
    this_vnp = next_vnp;    
  }
}


static void ExtractGeneFromListByGene (ValNodePtr PNTR list, SeqFeatPtr gene)
{
  ValNodePtr prev = NULL, this_vnp, next_vnp;
  
  if (list == NULL || gene == NULL)
  {
    return;
  }
  
  this_vnp = *list;
  while (this_vnp != NULL)
  {
    next_vnp = this_vnp->next;
    if (this_vnp->data.ptrvalue == gene)
    {
      if (prev == NULL)
      {
        *list = next_vnp;
      }
      else
      {
        prev->next = next_vnp;
      }
      this_vnp->next = NULL;
      ValNodeFree (this_vnp);
    }
    else
    {
      prev = this_vnp;
    }
    this_vnp = next_vnp;    
  }
}


static void 
CheckGenesForFeatureType 
(ValNodePtr PNTR features_without_genes,
 ValNodePtr PNTR superfluous_genes,
 BioseqPtr  bsp,
 Uint2      feature_type,
 Uint2      feature_subtype,
 Boolean    makes_gene_not_superfluous)
{
  SeqFeatPtr         sfp, gene_sfp;
  GeneRefPtr         grp;
  SeqMgrFeatContext  context;
  
  if (features_without_genes == NULL
      || superfluous_genes == NULL
      || bsp == NULL)
  {
    return;
  }
  
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, feature_type, feature_subtype, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, feature_type, feature_subtype, &context))
  {
    /* check for gene xref */
    grp = SeqMgrGetGeneXref (sfp);
    if (grp != NULL)
    {
      if (SeqMgrGeneIsSuppressed (grp))
      {
        ValNodeAddPointer (features_without_genes, OBJ_SEQFEAT, sfp);
      }
      else
      {
        ExtractGeneFromListByGeneRef (superfluous_genes, grp);
      }
    }
    else
    {
      gene_sfp = SeqMgrGetOverlappingGene (sfp->location, NULL);
      if (gene_sfp == NULL)
      {
        ValNodeAddPointer (features_without_genes, OBJ_SEQFEAT, sfp);
      }
      else if (makes_gene_not_superfluous)
      {
        ExtractGeneFromListByGene (superfluous_genes, gene_sfp);
      }
    }  
  }  
}

typedef struct misssupergenes
{
  ValNodePtr missing_list;
  ValNodePtr super_list;
} MissSuperGenesData, PNTR MissSuperGenesPtr;


static void FindMissingGenes (BioseqPtr bsp, Pointer userdata)
{
  SeqFeatPtr         sfp;
  SeqMgrFeatContext  context;
  ValNodePtr         features_without_genes = NULL;
  ValNodePtr         superfluous_genes = NULL;
  MissSuperGenesPtr  msgp;
  
  if (bsp == NULL || userdata == NULL)
  {
    return;
  }
  
  msgp = (MissSuperGenesPtr) userdata;
  
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, FEATDEF_GENE, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_GENE, FEATDEF_GENE, &context))
  {
    ValNodeAddPointer (&superfluous_genes, OBJ_SEQFEAT, sfp);
  }
  
  CheckGenesForFeatureType (&features_without_genes, &superfluous_genes, bsp, 
                            SEQFEAT_CDREGION, 0, TRUE);
  CheckGenesForFeatureType (&features_without_genes, &superfluous_genes, bsp, 
                            SEQFEAT_RNA, 0, TRUE);
  CheckGenesForFeatureType (&features_without_genes, &superfluous_genes, bsp, 
                            SEQFEAT_IMP, FEATDEF_RBS, FALSE);

  ValNodeLink (&(msgp->missing_list), features_without_genes);
  ValNodeLink (&(msgp->super_list), superfluous_genes);
}


static void 
GetPseudoAndNonPseudoGeneList 
(ValNodePtr      super_list,
 ValNodePtr PNTR pseudo_list, 
 ValNodePtr PNTR non_pseudo_list)
{
  ValNodePtr vnp;
  SeqFeatPtr gene;
  GeneRefPtr grp;
  
  if (pseudo_list == NULL || non_pseudo_list == NULL)
  {
    return;
  }
  *pseudo_list = NULL;
  *non_pseudo_list = NULL;
  
  for (vnp = super_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == OBJ_SEQFEAT)
    {
      gene = (SeqFeatPtr) vnp->data.ptrvalue;
      if (gene != NULL && gene->data.choice == SEQFEAT_GENE)
      {
        grp = (GeneRefPtr) gene->data.value.ptrvalue;
        if (gene->pseudo || (grp != NULL && grp->pseudo))
        {
          ValNodeAddPointer (pseudo_list, OBJ_SEQFEAT, gene);
        }
        else
        {
          ValNodeAddPointer (non_pseudo_list, OBJ_SEQFEAT, gene);
        }
      }
    }
  }
}


static void AddMissingAndSuperfluousGeneDiscrepancies (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list)
{
  DiscrepancyItemPtr dip, pseudo_dip, non_pseudo_dip;
  CharPtr            missing_genes_fmt = "%d features have no genes.";
  CharPtr            extra_genes_fmt = "%d gene features are not associated with a CDS or RNA feature.";
  CharPtr            pseudo_extra_genes_fmt = "%d pseudo gene features are not associated with a CDS or RNA feature.";
  CharPtr            non_pseudo_extra_genes_fmt = "%d non-pseudo gene features are not associated with a CDS or RNA feature.";
  MissSuperGenesData msgd;
  ValNodePtr         non_pseudo_list = NULL, pseudo_list = NULL;

  if (discrepancy_list == NULL)
  {
    return;
  }
  
  msgd.missing_list = NULL;
  msgd.super_list = NULL;
  
  VisitBioseqsInSep (sep, &msgd, FindMissingGenes);
  
  if (msgd.missing_list != NULL)
  {
    dip = NewDiscrepancyItem (DISC_GENE_MISSING, missing_genes_fmt, msgd.missing_list);
    if (dip != NULL)
    {
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
  
  if (msgd.super_list != NULL)
  {
    dip = NewDiscrepancyItem (DISC_SUPERFLUOUS_GENE, extra_genes_fmt, msgd.super_list);
    if (dip != NULL)
    {
      ValNodeAddPointer (discrepancy_list, 0, dip);      
      GetPseudoAndNonPseudoGeneList (msgd.super_list, &pseudo_list, &non_pseudo_list);
      non_pseudo_dip = NewDiscrepancyItem (DISC_SUPERFLUOUS_GENE, non_pseudo_extra_genes_fmt, non_pseudo_list);
      non_pseudo_dip->level = 1;
      pseudo_dip = NewDiscrepancyItem (DISC_SUPERFLUOUS_GENE, pseudo_extra_genes_fmt, pseudo_list);
      pseudo_dip->level = 1;
      ValNodeAddPointer (&(dip->subcategories), 0, non_pseudo_dip);
      ValNodeAddPointer (&(dip->subcategories), 0, pseudo_dip);
    }
  }  
}


typedef struct prefixcheck 
{
  CharPtr prefix;
  ValNodePtr feature_list;
} PrefixCheckData, PNTR PrefixCheckPtr;


static ValNodePtr FreePrefixCheckList (ValNodePtr prefix_list)
{
  PrefixCheckPtr pcp;
  
  if (prefix_list == NULL)
  {
    return NULL;
  }
  
  prefix_list->next = FreePrefixCheckList (prefix_list->next);
  
  pcp = (PrefixCheckPtr) prefix_list->data.ptrvalue;
  if (pcp != NULL)
  {
    pcp->prefix = MemFree (pcp->prefix);
    pcp->feature_list = ValNodeFree (pcp->feature_list);
    pcp = MemFree (pcp);
  }
  prefix_list = ValNodeFree (prefix_list);
  return NULL;
}


static DiscrepancyItemPtr InconsistentPrefix (PrefixCheckPtr pcp, CharPtr bad_fmt, DiscrepancyType disc_type)
{
  DiscrepancyItemPtr dip = NULL;

  if (pcp == NULL || StringHasNoText (pcp->prefix) || pcp->feature_list == NULL)
  {
    return NULL;
  }
  dip = (DiscrepancyItemPtr) MemNew (sizeof (DiscrepancyItemData));
  if (dip != NULL)
  {
    dip->disc_type = DISC_GENE_LOCUS_TAG_INCONSISTENT_PREFIX;
    dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + StringLen (pcp->prefix)+ 15));
    sprintf (dip->description, bad_fmt, ValNodeLen (pcp->feature_list), pcp->prefix);
    dip->callback_func = NULL;
    dip->datafree_func = NULL;
    dip->callback_data = NULL;
    dip->item_list = pcp->feature_list;
    pcp->feature_list = NULL;
  }      
  return dip;
}


static DiscrepancyItemPtr InconsistentLocusTagPrefix (PrefixCheckPtr pcp)
{
  CharPtr            bad_fmt = "%d features have locus tag prefix %s."; 

  return InconsistentPrefix (pcp, bad_fmt, DISC_GENE_LOCUS_TAG_INCONSISTENT_PREFIX);
}


static DiscrepancyItemPtr InconsistentProteinIDPrefix (PrefixCheckPtr pcp)
{
  CharPtr            bad_fmt = "%d sequences have protein ID prefix %s."; 

  return InconsistentPrefix (pcp, bad_fmt, DISC_INCONSISTENT_PROTEIN_ID_PREFIX);
}


typedef struct missinconstprotids 
{
  ValNodePtr missing_list;
  ValNodePtr inconsistent_list;
} MissInconstProtIDsData, PNTR MissInconstProtIDsPtr;


static void FindMissingAndInconsistentProteinIDsCallback (BioseqPtr bsp, Pointer userdata)
{
  SeqIdPtr              sip;
  MissInconstProtIDsPtr mipip;
  DbtagPtr              dbt = NULL;
  PrefixCheckPtr        pcp;
  ValNodePtr            vnp;
  
  if (bsp == NULL || ! ISA_aa (bsp->mol) || userdata == NULL)
  {
    return;
  }
  
  mipip = (MissInconstProtIDsPtr) userdata;
  
  for (sip = bsp->id; sip != NULL && dbt == NULL; sip = sip->next)
  {
    if (sip->choice == SEQID_GENERAL)
    {
      dbt = (DbtagPtr) sip->data.ptrvalue;
      if (dbt != NULL && StringICmp (dbt->db, "TMSMART") == 0)
      {
        dbt = NULL;
      }
    }
  }
  if (dbt == NULL)
  {
    ValNodeAddPointer (&(mipip->missing_list), OBJ_BIOSEQ, bsp);
  }
  else
  {  
    /* look for inconsistent prefixes */
    pcp = NULL;
    for (vnp = mipip->inconsistent_list; vnp != NULL && pcp == NULL; vnp = vnp->next)
    {
      pcp = (PrefixCheckPtr) vnp->data.ptrvalue;
      if (pcp != NULL && StringCmp (pcp->prefix, dbt->db) == 0)
      {
        ValNodeAddPointer (&pcp->feature_list, OBJ_BIOSEQ, bsp);
      }
      else
      {
        pcp = NULL;
      }
    }
    if (pcp == NULL)
    {
      pcp = (PrefixCheckPtr) MemNew (sizeof (PrefixCheckData));
      if (pcp != NULL)
      {
        pcp->prefix = StringSave (dbt->db);
        pcp->feature_list = ValNodeNew (NULL);
        pcp->feature_list->choice = OBJ_BIOSEQ;
        pcp->feature_list->data.ptrvalue = bsp;
        pcp->feature_list->next = NULL;
        ValNodeAddPointer (&(mipip->inconsistent_list), 0, pcp);
      }
    }
  }  
}


static void FindMissingProteinIDs (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list)
{
  CharPtr                bad_fmt = "%d proteins have invalid IDs.";
  DiscrepancyItemPtr     dip;
  MissInconstProtIDsData mipid;
  ValNodePtr             vnp;
  
  if (discrepancy_list == NULL) return;
  
  mipid.missing_list = NULL;
  mipid.inconsistent_list = NULL;
  
  VisitBioseqsInSep (sep, &mipid, FindMissingAndInconsistentProteinIDsCallback);
  
  if (mipid.missing_list != NULL)
  {
    dip = (DiscrepancyItemPtr) MemNew (sizeof (DiscrepancyItemData));
    if (dip != NULL)
    {
      dip->disc_type = DISC_MISSING_PROTEIN_ID;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (mipid.missing_list));
      dip->callback_func = NULL;
      dip->callback_data = NULL;
      dip->datafree_func = NULL;
      dip->item_list = mipid.missing_list;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
  

  if (mipid.inconsistent_list != NULL)
  {
    if (mipid.inconsistent_list->next != NULL)
    {
      for (vnp = mipid.inconsistent_list; vnp != NULL; vnp = vnp->next)
      {
        dip = InconsistentProteinIDPrefix (vnp->data.ptrvalue);
        if (dip != NULL)
        {
          ValNodeAddPointer (discrepancy_list, 0, dip);
        }
      }
    }
    mipid.inconsistent_list = FreePrefixCheckList (mipid.inconsistent_list);
  }
    
}

typedef struct locustagcheck
{
  ValNodePtr locus_tags_list;
  ValNodePtr missing_list;
  ValNodePtr duplicate_list;
  ValNodePtr bad_format_list;
  ValNodePtr inconsistent_prefix;
} LocusTagCheckData, PNTR LocusTagCheckPtr;

static void GeneLocusTagDiscrepancyCallback (ValNodePtr item_list, Pointer userdata)
{
  Message (MSG_OK, "I could launch the editor for the individual gene...");
}

static Boolean IsLocusTagFormatBad (CharPtr locus_tag)
{
  CharPtr cp;
  Boolean after_underscore = FALSE;
  
  if (StringHasNoText (locus_tag))
  {
    return FALSE;
  }
  
  cp = locus_tag;
  if (!isalpha (*cp))
  {
    return TRUE;
  }
  cp++;
  while (*cp != 0)
  {
    if (*cp == '_')
    {
      if (after_underscore)
      {
        return TRUE;
      }
      else
      {
        after_underscore = TRUE;
        if (*(cp + 1) == 0)
        {
          return TRUE;
        }
      }
    }
    else if (!isalpha (*cp) && !isdigit (*cp))
    {
      return TRUE;
    }
    cp++;
  }
  if (after_underscore)
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}

static CharPtr GetLocusTagPrefix (CharPtr locus_tag)
{
  Int4    prefix_len;
  CharPtr prefix;
  
  if (StringHasNoText (locus_tag))
  {
    return NULL;
  }
  
  prefix_len = StringCSpn (locus_tag, "_");
  if (prefix_len == 0 || prefix_len == StringLen (locus_tag))
  {
    return NULL;
  }
  else
  {
    prefix = (CharPtr) MemNew ((prefix_len + 1) * sizeof (Char));
    StringNCpy (prefix, locus_tag, prefix_len);
    prefix [prefix_len] = 0;
    return prefix;
  }
}

static void CheckGeneLocusTag (SeqFeatPtr sfp, Pointer userdata)
{
  GeneRefPtr         grp;
  LocusTagCheckPtr   ltcp;
  ValNodePtr         vnp;
  Boolean            found_duplicate;
  SeqFeatPtr         dup_sfp = NULL;
  CharPtr            prefix;
  PrefixCheckPtr     pcp;
  
  if (sfp == NULL || userdata == NULL || sfp->data.choice != SEQFEAT_GENE || sfp->data.value.ptrvalue == NULL)
  {
    return;
  }
  
  ltcp = (LocusTagCheckPtr) userdata;
  
  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (StringHasNoText (grp->locus_tag))
  {
    ValNodeAddPointer (&(ltcp->missing_list), OBJ_SEQFEAT, sfp);
  }
  else
  {
    /* look for badly formatted locus tags */
    if (IsLocusTagFormatBad (grp->locus_tag))
    {
      ValNodeAddPointer (&(ltcp->bad_format_list), OBJ_SEQFEAT, sfp);
    }
    
    /* look for inconsistent locus tag prefixes */
    prefix = GetLocusTagPrefix (grp->locus_tag);
    if (prefix != NULL)
    {
      pcp = NULL;
      for (vnp = ltcp->inconsistent_prefix; vnp != NULL && pcp == NULL; vnp = vnp->next)
      {
        pcp = (PrefixCheckPtr) vnp->data.ptrvalue;
        if (pcp != NULL && StringCmp (pcp->prefix, prefix) == 0)
        {
          ValNodeAddPointer (&pcp->feature_list, OBJ_SEQFEAT, sfp);
          prefix = MemFree (prefix);
        }
        else
        {
          pcp = NULL;
        }
      }
      if (pcp == NULL)
      {
        pcp = (PrefixCheckPtr) MemNew (sizeof (PrefixCheckData));
        if (pcp != NULL)
        {
          pcp->prefix = prefix;
          pcp->feature_list = ValNodeNew (NULL);
          pcp->feature_list->choice = OBJ_SEQFEAT;
          pcp->feature_list->data.ptrvalue = sfp;
          pcp->feature_list->next = NULL;
          ValNodeAddPointer (&(ltcp->inconsistent_prefix), 0, pcp);
        }
      }
    }
    
  
    dup_sfp = NULL;
    for (vnp = ltcp->locus_tags_list, found_duplicate = FALSE;
         vnp != NULL && !found_duplicate;
         vnp = vnp->next)
    { 
      dup_sfp = (SeqFeatPtr) vnp->data.ptrvalue;
      if (dup_sfp != NULL && dup_sfp->data.choice == SEQFEAT_GENE)
      {
        grp = (GeneRefPtr) dup_sfp->data.value.ptrvalue;
        if (grp != NULL && StringCmp (vnp->data.ptrvalue, grp->locus_tag) == 0)
        {
          found_duplicate = TRUE;
          if (vnp->choice != 0)
          {
            dup_sfp = NULL;
          }
          else
          {
            vnp->choice = 1;
          }
        }
        else
        {
          dup_sfp = NULL;
        }
      }
      else
      {
        dup_sfp = NULL;
      }
    }
    if (!found_duplicate)
    {
      ValNodeAddPointer (&(ltcp->locus_tags_list), 0, sfp);
    }
    else
    {
      if (dup_sfp != NULL)
      {
        ValNodeAddPointer (&(ltcp->duplicate_list), OBJ_SEQFEAT, dup_sfp);
      }
      ValNodeAddPointer (&(ltcp->duplicate_list), OBJ_SEQFEAT, sfp);
    }
  }
}

static void AddDiscrepanciesForMissingOrNonUniqueGeneLocusTags (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list)
{
  LocusTagCheckData  ltcd;
  DiscrepancyItemPtr dip = NULL;
  CharPtr            missing_fmt = "%d genes have no locus tags.";
  CharPtr            duplicate_fmt = "%d genes have duplicate locus tags.";
  CharPtr            bad_fmt = "%d locus tags are incorrectly formatted.";
  ValNodePtr         vnp;
  
  if (discrepancy_list == NULL) return;
  ltcd.locus_tags_list = NULL;
  ltcd.missing_list = NULL;
  ltcd.duplicate_list = NULL;
  ltcd.bad_format_list = NULL;
  ltcd.inconsistent_prefix = NULL;
  
  VisitFeaturesInSep (sep, &ltcd, CheckGeneLocusTag);
  ltcd.locus_tags_list = ValNodeFree (ltcd.locus_tags_list);
  
  if (ltcd.missing_list != NULL)
  {
    dip = (DiscrepancyItemPtr) MemNew (sizeof (DiscrepancyItemData));
    if (dip != NULL)
    {
      dip->disc_type = DISC_GENE_MISSING_LOCUS_TAG;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (missing_fmt) + 15));
      sprintf (dip->description, missing_fmt, ValNodeLen (ltcd.missing_list));
      dip->callback_func = GeneLocusTagDiscrepancyCallback;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = ltcd.missing_list;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
  
  if (ltcd.duplicate_list != NULL)
  {
    dip = (DiscrepancyItemPtr) MemNew (sizeof (DiscrepancyItemData));
    if (dip != NULL)
    {
      dip->disc_type = DISC_GENE_DUPLICATE_LOCUS_TAG;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (duplicate_fmt) + 15));
      sprintf (dip->description, duplicate_fmt, ValNodeLen (ltcd.duplicate_list));
      dip->callback_func = GeneLocusTagDiscrepancyCallback;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = ltcd.duplicate_list;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  } 
  
  if (ltcd.bad_format_list != NULL)
  {
    dip = (DiscrepancyItemPtr) MemNew (sizeof (DiscrepancyItemData));
    if (dip != NULL)
    {
      dip->disc_type = DISC_GENE_LOCUS_TAG_BAD_FORMAT;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (ltcd.bad_format_list));
      dip->callback_func = NULL;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = ltcd.bad_format_list;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
  
  if (ltcd.inconsistent_prefix != NULL)
  {
    if (ltcd.inconsistent_prefix->next != NULL)
    {
      for (vnp = ltcd.inconsistent_prefix; vnp != NULL; vnp = vnp->next)
      {
        dip = InconsistentLocusTagPrefix (vnp->data.ptrvalue);
        if (dip != NULL)
        {
          ValNodeAddPointer (discrepancy_list, 0, dip);
        }
      }
    }
    ltcd.inconsistent_prefix = FreePrefixCheckList (ltcd.inconsistent_prefix);
  }
}

static void NonGeneLocusTagDiscrepancyCallback (ValNodePtr item_list, Pointer userdata)
{
  Message (MSG_OK, "I could launch the editor for the individual gene, or I could remove all non-gene locus tags, or I could convert all non-gene locus tags to notes...");
}

static void AddDiscrepancyForNonGeneLocusTag (SeqFeatPtr sfp, Pointer userdata)
{
  ValNodePtr PNTR    locus_tag_list;
  GBQualPtr          qual;
  
  if (sfp == NULL || userdata == NULL || sfp->data.choice == SEQFEAT_GENE)
  {
    return;
  }
  
  locus_tag_list = (ValNodePtr PNTR) userdata;
  
  for (qual = sfp->qual; qual != NULL; qual = qual->next)
  {
    if (StringICmp(qual->qual, "locus_tag") == 0) 
    {
      ValNodeAddPointer (locus_tag_list, OBJ_SEQFEAT, sfp);
      return;
    }
  }
}

static void AddDiscrepanciesForNonGeneLocusTags (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list)
{
  ValNodePtr locus_tag_list = NULL;
  CharPtr    bad_fmt = "%d non-gene features have locus tags.";
  DiscrepancyItemPtr dip;
  
  VisitFeaturesInSep (sep, &locus_tag_list, AddDiscrepancyForNonGeneLocusTag);
  if (locus_tag_list != NULL)
  {
    dip = (DiscrepancyItemPtr) MemNew (sizeof (DiscrepancyItemData));
    if (dip != NULL)
    {
      dip->disc_type = DISC_NON_GENE_LOCUS_TAG;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (locus_tag_list));
      dip->callback_func = NonGeneLocusTagDiscrepancyCallback;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = locus_tag_list;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }  
}


static Boolean 
IsGeneLocationOk 
(SeqMgrFeatContextPtr feat_context, 
 SeqMgrFeatContextPtr gene_context,
 BioseqPtr            bsp)
{
  SeqFeatPtr        rbs_sfp;
  SeqMgrFeatContext rbs_context;
  
  if (feat_context == NULL || gene_context == NULL)
  {
    return FALSE;
  }  
  else if (feat_context->strand != gene_context->strand)
  {
    return FALSE;
  }
  else if (gene_context->left == feat_context->left && gene_context->right == feat_context->right)
  {
    return TRUE;
  }
  else if ((gene_context->strand == Seq_strand_minus && gene_context->left == feat_context->left)
           || (gene_context->strand != Seq_strand_minus && gene_context->right == feat_context->right))
  {
    /* find RBS to extend gene on 5' end */
    for (rbs_sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_RBS, &rbs_context);
         rbs_sfp != NULL;
         rbs_sfp = SeqMgrGetNextFeature (bsp, rbs_sfp, 0, FEATDEF_RBS, &rbs_context))
    {
      if (rbs_context.strand != gene_context->strand)
      {
        continue;
      }
      if (rbs_context.strand == Seq_strand_minus)
      {
        if (rbs_context.right == gene_context->right 
            && rbs_context.left >= feat_context->right)
        {
          return TRUE;
        }
      }
      else
      {
        if (rbs_context.left == gene_context->left
            && rbs_context.right <= feat_context->left)
        {
          return  TRUE;
        }
      }
    }
  }
  return FALSE;
}

static void
CheckFeatureTypeForLocationDiscrepancies 
(BioseqPtr       bsp, 
 Uint2           feature_type,
 ValNodePtr PNTR feature_list)
{
  SeqMgrFeatContext context, gene_context;
  GeneRefPtr        grp;
  SeqFeatPtr        sfp, gene_sfp;
  ValNodePtr        found_genes = NULL;
  Boolean           found_match;
  
  if (bsp == NULL || ISA_aa (bsp->mol) || feature_list == NULL)
  {
    return;
  }

  for (sfp = SeqMgrGetNextFeature (bsp, NULL, feature_type, 0, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, feature_type, 0, &context))
  {
    grp = SeqMgrGetGeneXref (sfp);
    if (grp == NULL)
    {
      gene_sfp = SeqMgrGetOverlappingGene (sfp->location, &gene_context);
      if (gene_sfp != NULL && !IsGeneLocationOk (&context, &gene_context, bsp))
      {
        ValNodeAddPointer (feature_list, OBJ_SEQFEAT, sfp);
        ValNodeAddPointer (feature_list, OBJ_SEQFEAT, gene_sfp);
      }
    }
    else if (!SeqMgrGeneIsSuppressed (grp))
    {
      found_match = FALSE;
      for (gene_sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, FEATDEF_GENE, &gene_context);
           gene_sfp != NULL && ! found_match;
           gene_sfp = SeqMgrGetNextFeature (bsp, gene_sfp, SEQFEAT_GENE, FEATDEF_GENE, &gene_context))
      {
        if (GeneRefMatch (gene_sfp->data.value.ptrvalue, grp) && gene_context.strand == context.strand)
        {
          if (IsGeneLocationOk (&context, &gene_context, bsp))
          {
            found_match = TRUE;
          }
          else
          {
            ValNodeAddPointer (&found_genes, OBJ_SEQFEAT, gene_sfp);
          }
        }
      }
      if (found_match || found_genes == NULL)
      {
        found_genes = ValNodeFree (found_genes);
        continue;
      }
      ValNodeAddPointer (feature_list, OBJ_SEQFEAT, sfp);
      ValNodeLink (feature_list, found_genes);
    }
  }
  
}

static void CDSmRNAGeneLocationDiscrepanciesCallback (BioseqPtr bsp, Pointer userdata)
{
  ValNodePtr PNTR   feature_list;

  if (bsp == NULL || ! ISA_na (bsp->mol) || userdata == NULL)
  {
    return;
  }
  
  feature_list = (ValNodePtr PNTR) userdata;
  
  CheckFeatureTypeForLocationDiscrepancies (bsp, SEQFEAT_CDREGION, feature_list);
  CheckFeatureTypeForLocationDiscrepancies (bsp, SEQFEAT_RNA, feature_list);
}


static void FindCDSmRNAGeneLocationDiscrepancies (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list)
{
  ValNodePtr feature_list = NULL;
  CharPtr            bad_fmt = "%d coding regions, mRNAs, and genes have inconsistent locations.";
  DiscrepancyItemPtr dip;
  
  VisitBioseqsInSep (sep, &feature_list, CDSmRNAGeneLocationDiscrepanciesCallback);
  
  if (feature_list != NULL)
  {
    dip = (DiscrepancyItemPtr) MemNew (sizeof (DiscrepancyItemData));
    if (dip != NULL)
    {
      dip->disc_type = DISC_GENE_CDS_mRNA_LOCATION_CONFLICT;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (feature_list));
      dip->callback_func = NULL;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = feature_list;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }  
}


typedef struct cdsgeneproduct 
{
  SeqFeatPtr cds;
  CharPtr    gene_locus;
  CharPtr    product_name;
} CDSGeneProductData, PNTR CDSGeneProductPtr;


static CharPtr GetGeneLabel (SeqFeatPtr sfp)
{
  GeneRefPtr grp;
  SeqFeatPtr gene_sfp;
  
  grp = SeqMgrGetGeneXref (sfp);
  if (grp == NULL)
  {
    gene_sfp = SeqMgrGetOverlappingGene (sfp->location, NULL);
    if (gene_sfp != NULL)
    {
      grp = gene_sfp->data.value.ptrvalue;
    }
  }
  if (grp != NULL)
  {
    if (!StringHasNoText (grp->locus))
    {
      return grp->locus;
    }
  }
  return NULL;
}

static void FindCDSGeneProductConflictsCallback (SeqFeatPtr sfp, Pointer userdata)
{
  CDSGeneProductPtr cgpp, cgpp_compare;
  SeqMgrFeatContext context;
  ValNodePtr PNTR   cds_list;
  ValNodePtr        prev = NULL, vnp;
  Boolean           found_match = FALSE;
  Uint1             new_choice = 0;
  CharPtr           gene_label;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION || userdata == NULL)
  {
    return;
  }
  
  sfp = SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, 0, 0, sfp, &context);
  if (sfp == NULL)
  {
    return;
  }
  
  gene_label = GetGeneLabel (sfp);
  if (StringHasNoText (gene_label)) return;
  
  cgpp = (CDSGeneProductPtr) MemNew (sizeof (CDSGeneProductData));
  if (cgpp != NULL)
  {
    cgpp->cds = sfp;
    cgpp->gene_locus = gene_label;
    cgpp->product_name = StringSave (context.label);
    
    cds_list = (ValNodePtr PNTR) userdata;
    if (*cds_list == NULL)
    {
      ValNodeAddPointer (cds_list, 0, cgpp);
    }
    else
    {
      vnp = *cds_list;
      while (vnp != NULL && !found_match)
      {
        cgpp_compare = (CDSGeneProductPtr) vnp->data.ptrvalue;
        if (cgpp_compare != NULL 
            && StringCmp (cgpp_compare->gene_locus, cgpp->gene_locus) == 0
            && StringCmp (cgpp_compare->product_name, cgpp->product_name) != 0)
        {
          found_match = TRUE;
          vnp->choice = 1;
          new_choice = 1;
        }
        prev = vnp;
        vnp = vnp->next;
      }
      
      if (found_match)
      {
        vnp = prev;
        /* insert at end of matches */
        while (found_match && vnp != NULL)
        {
          cgpp_compare = (CDSGeneProductPtr) vnp->data.ptrvalue;
          if (cgpp_compare != NULL 
            && StringCmp (cgpp_compare->gene_locus, cgpp->gene_locus) != 0)
          {
            found_match = FALSE;
          }
          else
          {
            prev = vnp;
          }
          vnp = vnp->next;
        }
      }

      /* add to list */
      vnp = ValNodeNew (NULL);
      vnp->choice = new_choice;
      vnp->data.ptrvalue = cgpp;
      vnp->next = prev->next;
      prev->next = vnp;
      
    }
  }
}

static void FindCDSGeneProductConflicts (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list)
{
  ValNodePtr         cds_list = NULL, non_conflict = NULL, vnp;
  CDSGeneProductPtr  cgpp;
  CharPtr            bad_fmt = "%d coding regions have the same gene name as another coding region but a different product.";
  DiscrepancyItemPtr dip;
  
  VisitFeaturesInSep (sep, &cds_list, FindCDSGeneProductConflictsCallback);

  /* remove CDSs without conflicts */
  non_conflict = ValNodeExtractList (&cds_list, 0);
  non_conflict = ValNodeFree (non_conflict);
  
  /* for each item, replace structure used for search with just the feature */
  for (vnp = cds_list; vnp != NULL; vnp = vnp->next)
  {
    cgpp = (CDSGeneProductPtr) vnp->data.ptrvalue;
    if (cgpp != NULL)
    {
      vnp->data.ptrvalue = cgpp->cds;
      vnp->choice = OBJ_SEQFEAT;
      cgpp->product_name = MemFree (cgpp->product_name);
      /* note - we are not freeing gene_locus because we didn't make a copy */
      cgpp = MemFree (cgpp);
    }
  }
    
  if (cds_list != NULL)
  {
    dip = (DiscrepancyItemPtr) MemNew (sizeof (DiscrepancyItemData));
    if (dip != NULL)
    {
      dip->disc_type = DISC_GENE_PRODUCT_CONFLICT;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (cds_list));
      dip->callback_func = NULL;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = cds_list;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }  
}


static void DuplicateGeneLocusCallback (SeqFeatPtr sfp, Pointer userdata)
{
  ValNodePtr PNTR gene_list;
  SeqFeatPtr      sfp_compare;
  GeneRefPtr      grp, grp_compare;
  ValNodePtr      prev = NULL, vnp;
  Boolean         found_match = FALSE;
  Uint1           new_choice = 0;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_GENE || sfp->data.value.ptrvalue == NULL
      || userdata == NULL)
  {
    return;
  }
  
  gene_list = (ValNodePtr PNTR) userdata;
  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (StringHasNoText (grp->locus))
  {
    return;
  }
  
  if (*gene_list == NULL)
  {
    ValNodeAddPointer (gene_list, 0, sfp);
  }
  else
  {
    vnp = *gene_list;
    while (vnp != NULL && !found_match)
    {
      sfp_compare = (SeqFeatPtr) vnp->data.ptrvalue;
      grp_compare = (GeneRefPtr) sfp_compare->data.value.ptrvalue;
      if (StringCmp (grp_compare->locus, grp->locus) == 0)
      {
        found_match = TRUE;
        vnp->choice = OBJ_SEQFEAT;
        new_choice = OBJ_SEQFEAT;
      }
      prev = vnp;
      vnp = vnp->next;
    }
    
    if (found_match)
    {
      vnp = prev;
      /* insert at end of matches */
      while (found_match && vnp != NULL)
      {
        sfp_compare = (SeqFeatPtr) vnp->data.ptrvalue;
        grp_compare = (GeneRefPtr) sfp_compare->data.value.ptrvalue;
        if (StringCmp (grp_compare->locus, grp->locus) != 0)
        {
          found_match = FALSE;
        }
        else
        {
          prev = vnp;
        }
        vnp = vnp->next;
      }
    }

    /* add to list */
    vnp = ValNodeNew (NULL);
    vnp->choice = new_choice;
    vnp->data.ptrvalue = sfp;
    vnp->next = prev->next;
    prev->next = vnp;
      
  }
  
}


static void FindDuplicateGeneLocus (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list)
{
  ValNodePtr         gene_list = NULL, non_conflict = NULL;
  CharPtr            bad_fmt = "%d genes have the same locus as another gene.";
  DiscrepancyItemPtr dip;
  
  VisitFeaturesInSep (sep, &gene_list, DuplicateGeneLocusCallback);

  /* remove Genes without conflicts */
  non_conflict = ValNodeExtractList (&gene_list, 0);
  non_conflict = ValNodeFree (non_conflict);
  
  if (gene_list != NULL)
  {
    dip = (DiscrepancyItemPtr) MemNew (sizeof (DiscrepancyItemData));
    if (dip != NULL)
    {
      dip->disc_type = DISC_GENE_DUPLICATE_LOCUS;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (gene_list));
      dip->callback_func = NULL;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = gene_list;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }  
  
}


static void ECNumberNoteDiscrepancyCallback (ValNodePtr item_list, Pointer userdata)
{
  Message (MSG_OK, "I could launch the editor for the individual feature, or I could convert the note to an EC number qual");
}


static Boolean ECNumberInText (CharPtr text) 
{
  CharPtr cp_num, cp_dash, cp_start;
  Int4    nums_found = 0;
  
  if (StringHasNoText (text)) {
      return FALSE;
  }
  
  cp_num = text + StrCSpn (text, "1234567890");
  cp_dash = StringStr (text, "-");
  if (cp_dash == NULL || cp_dash > cp_num) {
    cp_start = cp_num;
  } else {
    cp_start = cp_dash;
  }
  
  while (*cp_start != 0) {
    nums_found++;
    while (isdigit (*cp_start) || *cp_start == '-') {
      cp_start++;
    }
    if (*cp_start != '.') {
      if (nums_found == 4) {
        return TRUE;
      } else {
        return ECNumberInText (cp_start);
      }
    } else {
      cp_start ++;
    }
  }    
  if (nums_found == 4) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static void FindECNumberNotes (SeqFeatPtr sfp, Pointer userdata)
{
  ValNodePtr PNTR    ec_number_features;
  BioseqPtr          prot_bsp;
  SeqMgrFeatContext  fcontext;
  SeqFeatPtr         prot_sfp;
  ProtRefPtr         prp;
  ValNodePtr         vnp;
  
  if (sfp == NULL || userdata == NULL || StringHasNoText (sfp->comment))
  {
    return;
  }
  
  ec_number_features = (ValNodePtr PNTR) userdata;
  
  if (ECNumberInText (sfp->comment))
  {
    ValNodeAddPointer (ec_number_features, OBJ_SEQFEAT, sfp);
  }
  else if (sfp->data.choice == SEQFEAT_CDREGION && sfp->product != NULL) 
  {
    prot_bsp = BioseqFindFromSeqLoc(sfp->product);
    prot_sfp = SeqMgrGetNextFeature(prot_bsp, NULL, SEQFEAT_PROT, FEATDEF_PROT, &fcontext);
    if (prot_sfp != NULL && prot_sfp->data.value.ptrvalue != NULL) {
      prp = (ProtRefPtr) prot_sfp->data.value.ptrvalue;
      for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
        if (ECNumberInText (vnp->data.ptrvalue)) {
          ValNodeAddPointer (ec_number_features, OBJ_SEQFEAT, sfp);
          return;
        }
      }
      if (ECNumberInText (prp->desc)) {
        ValNodeAddPointer (ec_number_features, OBJ_SEQFEAT, sfp);
        return;
      }
    }
  }  
}

static void AddECNumberNoteDiscrepancies (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list)
{
  ValNodePtr ec_number_features = NULL;
  DiscrepancyItemPtr dip;
  CharPtr            bad_fmt = "%d features have EC numbers in notes or products.";
  
  if (discrepancy_list == NULL) return;
  
  VisitFeaturesInSep (sep, &ec_number_features, FindECNumberNotes);
  
  if (ec_number_features != NULL)
  {
    dip = (DiscrepancyItemPtr) MemNew (sizeof (DiscrepancyItemData));
    if (dip != NULL)
    {
      dip->disc_type = DISC_EC_NUMBER_NOTE;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (ec_number_features));
      dip->callback_func = ECNumberNoteDiscrepancyCallback;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = ec_number_features;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
  
}


static void FindPseudoDiscrepanciesCallback (SeqFeatPtr sfp, Pointer userdata)
{
  ValNodePtr PNTR pseudo_features;
  GeneRefPtr      grp;
  SeqFeatPtr      gene_sfp = NULL;
  Boolean         found_conflict = FALSE;
  
  if (sfp == NULL || (sfp->data.choice != SEQFEAT_CDREGION && sfp->data.choice != SEQFEAT_RNA)
      || userdata == NULL)
  {
    return;
  }
  
  grp = SeqMgrGetGeneXref (sfp);
  if (grp != NULL)
  {
    return;
  }
  
  gene_sfp = SeqMgrGetOverlappingGene (sfp->location, NULL);
  if (gene_sfp == NULL)
  {
    return;
  }
  
  if ((sfp->pseudo && ! gene_sfp->pseudo)
      || (!sfp->pseudo && gene_sfp->pseudo))
  {
    pseudo_features = (ValNodePtr PNTR) userdata;
    ValNodeAddPointer (pseudo_features, OBJ_SEQFEAT, sfp);
    if (gene_sfp != NULL)
    {
      ValNodeAddPointer (pseudo_features, OBJ_SEQFEAT, gene_sfp);
    }
  }
}


static void FindPseudoDiscrepancies (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list)
{
  ValNodePtr pseudo_features = NULL;
  DiscrepancyItemPtr dip;
  CharPtr            bad_fmt = "%d CDSs, RNAs, and genes have mismatching pseudos.";
  
  if (discrepancy_list == NULL) return;
  
  VisitFeaturesInSep (sep, &pseudo_features, FindPseudoDiscrepanciesCallback);
  
  if (pseudo_features != NULL)
  {
    dip = (DiscrepancyItemPtr) MemNew (sizeof (DiscrepancyItemData));
    if (dip != NULL)
    {
      dip->disc_type = DISC_PSEUDO_MISMATCH;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (pseudo_features));
      dip->callback_func = NULL;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = pseudo_features;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
  
}


static void FindJoinedLocations (SeqFeatPtr sfp, Pointer userdata)
{
  ValNodePtr PNTR joined_features;
  
  if (sfp == NULL || userdata == NULL || sfp->location == NULL)
  {
    return;
  }
  
  joined_features = (ValNodePtr PNTR) userdata;
  if (sfp->location->choice == SEQLOC_MIX)
  {
    ValNodeAddPointer (joined_features, OBJ_SEQFEAT, sfp);
  }
}

static void AddJoinedFeatureDiscrepancies (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list)
{
  ValNodePtr joined_features = NULL;
  
  DiscrepancyItemPtr dip;
  CharPtr            bad_fmt = "%d features have joined locations.";
  
  if (discrepancy_list == NULL) return;
  
  VisitFeaturesInSep (sep, &joined_features, FindJoinedLocations);
  
  if (joined_features != NULL)
  {
    dip = (DiscrepancyItemPtr) MemNew (sizeof (DiscrepancyItemData));
    if (dip != NULL)
    {
      dip->disc_type = DISC_JOINED_FEATURES;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (joined_features));
      dip->callback_func = NULL;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = joined_features;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
}


static void FindOverlappingGenes (BioseqPtr bsp, Pointer userdata)
{
  SeqFeatPtr         sfp, sfp_compare;
  SeqMgrFeatContext  context;
  ValNodePtr PNTR    overlapping_genes = NULL, non_overlap;
  ValNodePtr         gene_list = NULL, vnp, vnp_next;
  
  if (bsp == NULL || userdata == NULL)
  {
    return;
  }
  
  overlapping_genes = (ValNodePtr PNTR) userdata;
  
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, FEATDEF_GENE, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_GENE, FEATDEF_GENE, &context))
  {
    ValNodeAddPointer (&gene_list, 0, sfp);
  }
  
  for (vnp = gene_list; vnp != NULL && vnp->next != NULL; vnp = vnp->next)
  {
    sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    for (vnp_next = vnp->next; vnp_next != NULL; vnp_next = vnp_next->next)
    {
      sfp_compare = (SeqFeatPtr) vnp_next->data.ptrvalue;
      
      if (SeqLocStrand (sfp->location) != SeqLocStrand (sfp_compare->location))
      {
        continue;
      }
      
      if (SeqLocCompare (sfp->location, sfp_compare->location) != SLC_NO_MATCH)
      {
        vnp->choice = OBJ_SEQFEAT;
        vnp_next->choice = OBJ_SEQFEAT;
      }
    }
  }
  
  non_overlap = ValNodeExtractList (&gene_list, 0);
  non_overlap = ValNodeFree (non_overlap);
  ValNodeLink (overlapping_genes, gene_list);
  
}

static void AddOverlappingGeneDiscrepancies (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list)
{
  DiscrepancyItemPtr dip;
  CharPtr            bad_fmt = "%d genes overlap another gene on the same strand.";
  ValNodePtr         overlapping_genes = NULL;

  if (discrepancy_list == NULL)
  {
    return;
  }
  
  VisitBioseqsInSep (sep, &overlapping_genes, FindOverlappingGenes);
  
  if (overlapping_genes != NULL)
  {
    dip = (DiscrepancyItemPtr) MemNew (sizeof (DiscrepancyItemData));
    if (dip != NULL)
    {
      dip->disc_type = DISC_OVERLAPPING_GENES;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (overlapping_genes));
      dip->callback_func = NULL;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = overlapping_genes;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
}


typedef struct cdsoverlap 
{
  CharPtr    product_name;
  SeqFeatPtr sfp;  
  Int4       left;
  Int4       right;
} CDSOverlapData, PNTR CDSOverlapPtr;


static CDSOverlapPtr CDSOverlapNew (SeqFeatPtr sfp, CharPtr product_name, Int4 left, Int4 right)
{
  CDSOverlapPtr cop;
  
  cop = (CDSOverlapPtr) MemNew (sizeof (CDSOverlapData));
  if (cop != NULL)
  {
    cop->product_name = StringSave (product_name);
    cop->sfp = sfp;
    cop->left = left;
    cop->right = right;
  }
  return cop;
}


static ValNodePtr FreeCDSOverlapList (ValNodePtr vnp)
{
  CDSOverlapPtr cop;
  
  if (vnp != NULL)  
  {
    vnp->next = FreeCDSOverlapList (vnp->next);
    cop = (CDSOverlapPtr) vnp->data.ptrvalue;
    if (cop != NULL)
    {
      cop->product_name = MemFree (cop->product_name);
      cop = MemFree (cop);
      vnp->data.ptrvalue = NULL;
    }
    vnp = ValNodeFree (vnp);
  }
  return vnp;
}


static ValNodePtr FeatureListFromOverlapList (ValNodePtr vnp)
{
  ValNodePtr     feat_list = NULL;
  CDSOverlapPtr cop;
  
  while (vnp != NULL)
  {
    if (vnp->choice != 0 && vnp->data.ptrvalue != NULL)
    {
      cop = (CDSOverlapPtr) vnp->data.ptrvalue;
      ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, cop->sfp);
    }
    vnp = vnp->next;
  }
  return feat_list;
}


static CharPtr similar_product_words[] = 
{ "transposase",
  "integrase"
};

const int num_similar_product_words = sizeof (similar_product_words) / sizeof (CharPtr);

static CharPtr ignore_similar_product_words[] = 
{ "hypothetical protein",
  "phage"
};

const int num_ignore_similar_product_words = sizeof (ignore_similar_product_words) / sizeof (CharPtr);


static Boolean OverlappingProductNameSimilar (CharPtr str1, CharPtr str2)
{
  Int4 i;
  Boolean str1_has_similarity_word = FALSE, str2_has_similarity_word = FALSE;
  
  if (StringHasNoText (str1) && StringHasNoText (str2))
  {
    return TRUE;
  }
  else if (StringHasNoText (str1) || StringHasNoText (str2))
  {
    return FALSE;
  }
  
  /* if both product names contain one of the special case similarity words,
   * the product names are similar. */
  for (i = 0; i < num_similar_product_words; i++)
  {
    if (StringISearch (str1, similar_product_words [i]) != NULL)
    {
      str1_has_similarity_word = TRUE;
    }
    if (StringISearch (str2, similar_product_words [i]) != NULL)
    {
      str2_has_similarity_word = TRUE;
    }
  }
  if (str1_has_similarity_word && str2_has_similarity_word)
  {
    return TRUE;
  }
  
  /* otherwise, if one of the product names contains one of special ignore similarity
   * words, the product names are not similar.
   */
  for (i = 0; i < num_ignore_similar_product_words; i++)
  {
    if (StringISearch (str1, ignore_similar_product_words[i]) != NULL
        || StringISearch (str2, ignore_similar_product_words[i]) != NULL)
    {
      return FALSE;
    }
  }
  
  if (StringICmp (str1, str2) == 0)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}


static void FindOverlappingCDSs (BioseqPtr bsp, Pointer userdata)
{
  SeqFeatPtr         sfp;
  SeqMgrFeatContext  context;
  ValNodePtr PNTR    overlapping_cds = NULL, cds_list;
  ValNodePtr         overlap_list = NULL, vnp, vnp_next;
  CDSOverlapPtr      cop, cop_compare;
  
  if (bsp == NULL || userdata == NULL)
  {
    return;
  }
  
  overlapping_cds = (ValNodePtr PNTR) userdata;
  
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, FEATDEF_CDS, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, FEATDEF_CDS, &context))
  {
    ValNodeAddPointer (&overlap_list, 0, CDSOverlapNew (sfp, context.label, context.left, context.right));
  }
  
  for (vnp = overlap_list; vnp != NULL && vnp->next != NULL; vnp = vnp->next)
  {
    cop = (CDSOverlapPtr) vnp->data.ptrvalue;
    if (cop == NULL)
    {
      continue;
    }
    for (vnp_next = vnp->next; vnp_next != NULL; vnp_next = vnp_next->next)
    {
      cop_compare = (CDSOverlapPtr) vnp_next->data.ptrvalue;
      if (cop_compare == NULL)
      {
        continue;
      }
      else if (cop_compare->left > cop->right)
      {
        break;
      }
      if (!OverlappingProductNameSimilar (cop->product_name, cop_compare->product_name))
      {
        continue;
      }
      if (SeqLocStrand (cop->sfp->location) != SeqLocStrand (cop_compare->sfp->location))
      {
        continue;
      }
      
      if (SeqLocCompare (cop->sfp->location, cop_compare->sfp->location) != SLC_NO_MATCH)
      {
        vnp->choice = OBJ_SEQFEAT;
        vnp_next->choice = OBJ_SEQFEAT;
      }
    }
  }
  
  cds_list = FeatureListFromOverlapList(overlap_list);
  if (cds_list != NULL)
  {
    ValNodeLink (overlapping_cds, cds_list);
  }
  overlap_list = FreeCDSOverlapList (overlap_list);
}


static void AddOverlappingCodingRegionDiscrepancies (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list)
{
  DiscrepancyItemPtr dip;
  CharPtr            bad_fmt = "%d coding regions overlap another coding region with a similar or identical name.";
  ValNodePtr         overlapping_cds = NULL;

  if (discrepancy_list == NULL)
  {
    return;
  }
  
  VisitBioseqsInSep (sep, &overlapping_cds, FindOverlappingCDSs);
  
  if (overlapping_cds != NULL)
  {
    dip = NewDiscrepancyItem (DISC_OVERLAPPING_CDS, bad_fmt, overlapping_cds);

    if (dip != NULL)
    {
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
}


static void FindShortContigsCallback (BioseqPtr bsp, Pointer userdata)
{
  ValNodePtr PNTR bioseq_list;
  
  if (bsp == NULL || !ISA_na (bsp->mol) || userdata == NULL || bsp->length >= 200)
  {
    return;
  }
  
  bioseq_list = (ValNodePtr PNTR) userdata;
  
  ValNodeAddPointer (bioseq_list, OBJ_BIOSEQ, bsp);
}

static void FindShortContigs (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list)
{
  DiscrepancyItemPtr dip;
  CharPtr            bad_fmt = "%d contigs are shorter than 200 nt.";
  ValNodePtr         bioseq_list = NULL;
  
  if (discrepancy_list == NULL) return;
  
  VisitBioseqsInSep (sep, &bioseq_list, FindShortContigsCallback);
  
  if (bioseq_list != NULL)
  {
    dip = (DiscrepancyItemPtr) MemNew (sizeof (DiscrepancyItemData));
    if (dip != NULL)
    {
      dip->disc_type = DISC_SHORT_CONTIG;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (bioseq_list));
      dip->callback_func = NULL;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = bioseq_list;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
}

typedef struct biosrccheck 
{
  BioSourcePtr biop;
  ValNodePtr   sdp_list;
} BioSrcCheckData, PNTR BioSrcCheckPtr;

static ValNodePtr FreeBioSrcCheckList (ValNodePtr biosrc_list)
{
  BioSrcCheckPtr  bscp;
  
  if (biosrc_list == NULL)
  {
    return NULL;
  }
  
  biosrc_list->next = FreeBioSrcCheckList (biosrc_list->next);
  
  bscp = (BioSrcCheckPtr) biosrc_list->data.ptrvalue;
  if (bscp != NULL)
  {
    bscp->sdp_list = ValNodeFree (bscp->sdp_list);
    bscp = MemFree (bscp);
  }
  biosrc_list = ValNodeFree (biosrc_list);
  return NULL;
}


static void FindInconsistentSourcesCallback (BioseqPtr bsp, Pointer userdata)
{
  ValNodePtr PNTR biosrc_list, vnp;
  SeqDescrPtr     sdp;
  BioSrcCheckPtr  bscp;
  Boolean         found = FALSE;
  SeqMgrDescContext context;
  
  if (bsp == NULL || !ISA_na (bsp->mol) || userdata == NULL)
  {
    return;
  }
  
  biosrc_list = (ValNodePtr PNTR) userdata;
  
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
  if (sdp != NULL)
  {
    for (vnp = *biosrc_list; vnp != NULL && !found; vnp = vnp->next)
    {
      bscp = (BioSrcCheckPtr) vnp->data.ptrvalue;
      if (bscp != NULL && BioSourceMatch (sdp->data.ptrvalue, bscp->biop))
      {
        ValNodeAddPointer (&(bscp->sdp_list), OBJ_SEQDESC, sdp);
        found = TRUE;
      }
    }
    if (!found)
    {
      bscp = (BioSrcCheckPtr) MemNew (sizeof (BioSrcCheckData));
      if (bscp != NULL)
      {
        bscp->biop = sdp->data.ptrvalue;
        ValNodeAddPointer (&(bscp->sdp_list), OBJ_SEQDESC, sdp);
        ValNodeAddPointer (biosrc_list, 0, bscp);
      }
    }
  }
}


static DiscrepancyItemPtr InconsistentBiosrc (BioSrcCheckPtr bscp)
{
  DiscrepancyItemPtr dip = NULL;
  CharPtr            bad_fmt = "%d contigs have identical sources that do not match another contig source.";

  if (bscp == NULL || bscp->sdp_list == NULL)
  {
    return NULL;
  }
  
  dip = (DiscrepancyItemPtr) MemNew (sizeof (DiscrepancyItemData));
  if (dip != NULL)
  {
    dip->disc_type = DISC_INCONSISTENT_BIOSRC;
    dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
    sprintf (dip->description, bad_fmt, ValNodeLen (bscp->sdp_list));
    dip->callback_func = NULL;
    dip->datafree_func = NULL;
    dip->callback_data = NULL;
    dip->item_list = bscp->sdp_list;
    bscp->sdp_list = NULL;
  }      
  return dip;
}


static void FindNonmatchingContigSources (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list)
{
  DiscrepancyItemPtr dip;
  ValNodePtr         biosrc_list = NULL, vnp;
  
  if (discrepancy_list == NULL) return;
  
  VisitBioseqsInSep (sep, &biosrc_list, FindInconsistentSourcesCallback);
  
  if (biosrc_list != NULL && biosrc_list->next != NULL)
  {
    for (vnp = biosrc_list; vnp != NULL; vnp = vnp->next)
    {
      dip = InconsistentBiosrc (vnp->data.ptrvalue);
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
  biosrc_list = FreeBioSrcCheckList (biosrc_list);
}


static CharPtr suspect_product_names[] = 
{
"Similar to",
"Related to",
"interrupt",
"Homolog",
"Homologue",
"Fragment",
"Frameshift",
"Intein",
"COG",
"N-terminal",
"N terminal",
"C-terminal",
"C terminal",
"N-terminus",
"C-terminus",
"Chloroplast",
"Mitochondrial",
"Brackets or parenthesis [] ()",
"ending with period or comma",
"ortholog",
"orthologue",
"paralog",
"paralogue",
"bifunctional",
"pseudogene",
"frame shift",
"protien",
"partial",
"B.subtilis",
"E.coli",
"Escherichia",
"coli",
"Bacillus",
"Subtilis",
"Staphlococcus",
"aureus",
"Salmonella",
"Streptococcus",
"Staphlococcal",
"streptococcal",
"Helicobacter",
"pylori",
"Campylobacter",
"Jejuni",
"Yersinia",
"Pestis",
"Rhodobacter",
"sphaeroides",
};

const int num_suspect_product_names = sizeof (suspect_product_names) / sizeof (CharPtr);

const int brackets_name = 17;
const int end_with_punct_name = 18;


static void FindSuspectProductNamesCallback (SeqFeatPtr sfp, Pointer userdata)
{
  ValNodePtr PNTR feature_list;
  Int4            k, len;
  ProtRefPtr      prp;
  ValNodePtr      vnp;
  CharPtr         str;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_PROT || sfp->data.value.ptrvalue == NULL
      || userdata == NULL)
  {
    return;
  }
  
  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  feature_list = (ValNodePtr PNTR) userdata;
  
  for (k = 0; k < num_suspect_product_names; k++)
  {
    if (k == brackets_name)
    {
      for (vnp = prp->name; vnp != NULL; vnp = vnp->next)
      {
        if (StringChr (vnp->data.ptrvalue, '[') != NULL
            || StringChr (vnp->data.ptrvalue, ']') != NULL
            || StringChr (vnp->data.ptrvalue, '(') != NULL
            || StringChr (vnp->data.ptrvalue, ')') != NULL)
        {
          ValNodeAddPointer (&(feature_list[k]), OBJ_SEQFEAT, sfp);
          break;
        }
      }
    }
    else if (k == end_with_punct_name)
    {
      for (vnp = prp->name; vnp != NULL; vnp = vnp->next)
      {
        len = StringLen (vnp->data.ptrvalue);
        str = (CharPtr) vnp->data.ptrvalue;
        if (str[len - 1] == '.' || str[len - 1] == ',')
        {
          ValNodeAddPointer (&(feature_list[k]), OBJ_SEQFEAT, sfp);
          break;
        }
      }
    }
    else
    {
      if (StringLen (suspect_product_names [k]) < 4)
      {
        for (vnp = prp->name; vnp != NULL; vnp = vnp->next)
        {
          if (StringSearch(vnp->data.ptrvalue, suspect_product_names[k]) != NULL)
          {
            ValNodeAddPointer (&(feature_list[k]), OBJ_SEQFEAT, sfp);
            break;
          }
        }
      }
      else
      {
        for (vnp = prp->name; vnp != NULL; vnp = vnp->next)
        {
          if (StringISearch(vnp->data.ptrvalue, suspect_product_names[k]) != NULL)
          {
            ValNodeAddPointer (&(feature_list[k]), OBJ_SEQFEAT, sfp);
            break;
          }
        }
      }
    }
  }
  
}


static DiscrepancyItemPtr SuspectProductName (CharPtr product_name, ValNodePtr feature_list)
{
  DiscrepancyItemPtr dip = NULL;
  CharPtr            bad_fmt = "%d product names contain %s";

  if (feature_list == NULL || StringHasNoText (product_name))
  {
    return NULL;
  }
  
  dip = (DiscrepancyItemPtr) MemNew (sizeof (DiscrepancyItemData));
  if (dip != NULL)
  {
    dip->disc_type = DISC_SUSPECT_PRODUCT_NAME;
    dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + StringLen (product_name) + 15));
    sprintf (dip->description, bad_fmt, ValNodeLen (feature_list), product_name);
    dip->callback_func = NULL;
    dip->datafree_func = NULL;
    dip->callback_data = NULL;
    dip->item_list = feature_list;
  }      
  return dip;
}


static void FindSuspectProductNames (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list)
{
  ValNodePtr PNTR   feature_list = NULL;
  ValNodePtr         master_list = NULL;
  Int4               k;
  DiscrepancyItemPtr dip;
  ValNodePtr         subcategories = NULL;
  
  if (discrepancy_list == NULL) return;

  feature_list = (ValNodePtr PNTR) MemNew (sizeof (ValNodePtr) * num_suspect_product_names);
  if (feature_list == NULL) return;
  
  /* initialize array for suspicious product names */
  for (k = 0; k < num_suspect_product_names; k++)
  {
    feature_list[k] = NULL;
  }
  
  VisitFeaturesInSep (sep, feature_list, FindSuspectProductNamesCallback);
  
  for (k = 0; k < num_suspect_product_names; k++)
  {
    if (feature_list[k] != NULL)
    {
      dip = SuspectProductName (suspect_product_names[k], feature_list[k]);
      if (dip != NULL)
      {
        ValNodeAddPointer (&subcategories, 0, dip);
      }
      ValNodeLinkCopy (&master_list, feature_list[k]);
    }
  }
  
  if (master_list != NULL)
  {
    dip = SuspectProductName ("suspect phrase or characters", master_list);
    if (dip != NULL)
    {
      dip->subcategories = subcategories;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }

  MemFree (feature_list);
}


static DiscrepancyItemPtr InconsistentSourceDefline (SeqDescrPtr biop_sdp, SeqDescrPtr title_sdp)
{
  DiscrepancyItemPtr dip = NULL;
  CharPtr            bad_fmt = "Organism description not found in definition line: %s.";
  BioSourcePtr       biop;
  CharPtr            desc = NULL;

  if (biop_sdp == NULL || title_sdp == NULL)
  {
    return NULL;
  }
  
  biop = (BioSourcePtr) biop_sdp->data.ptrvalue;
  if (biop != NULL && biop->org != NULL && !StringHasNoText (biop->org->taxname))
  {
    desc = biop->org->taxname;
  }
  else
  {
    desc = title_sdp->data.ptrvalue;
  }
  if (StringHasNoText (desc)) {
    return NULL;
  }
  
  dip = (DiscrepancyItemPtr) MemNew (sizeof (DiscrepancyItemData));
  if (dip != NULL)
  {
    dip->disc_type = DISC_INCONSISTENT_BIOSRC_DEFLINE;
    dip->description = (CharPtr)MemNew (StringLen (bad_fmt) + StringLen (desc));
    sprintf (dip->description, bad_fmt, desc);
    dip->callback_func = NULL;
    dip->datafree_func = NULL;
    dip->callback_data = NULL;
    dip->item_list = NULL;
    ValNodeAddPointer (&(dip->item_list), OBJ_SEQDESC, biop_sdp);
    ValNodeAddPointer (&(dip->item_list), OBJ_SEQDESC, title_sdp);
  }      
  return dip;
}


static void FindInconsistentSourceAndDeflineCallback (BioseqPtr bsp, Pointer userdata)
{
  DiscrepancyItemPtr dip;
  ValNodePtr PNTR discrepancy_list;
  SeqDescrPtr        biop_sdp, title_sdp;
  SeqMgrDescContext  context;
  BioSourcePtr       biop;
  
  discrepancy_list = (ValNodePtr PNTR) userdata;
  if (bsp == NULL || discrepancy_list == NULL) return;
  
  biop_sdp = SeqMgrGetNextDescriptor(bsp, NULL, Seq_descr_source, &context);
  if (biop_sdp == NULL || biop_sdp->data.ptrvalue == NULL)
  {
    return;
  }
  biop = (BioSourcePtr) biop_sdp->data.ptrvalue;
  if (biop->org == NULL)
  {
    return;
  }
  if (StringHasNoText (biop->org->taxname)) 
  {
    return;
  }
  
  title_sdp = SeqMgrGetNextDescriptor(bsp, NULL, Seq_descr_title, &context);
  if (title_sdp == NULL) return;
  
  if (StringStr (title_sdp->data.ptrvalue, biop->org->taxname) == NULL)
  {
    dip = InconsistentSourceDefline (biop_sdp, title_sdp);
    if (dip != NULL)
    {
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
}


static void FindInconsistentSourceAndDefline (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list)
{  
  ValNodePtr disc_pairs = NULL;
  CharPtr    bad_fmt = "%d sources do not match definition lines.";
  DiscrepancyItemPtr dip;
  
  if (discrepancy_list == NULL) return;

  VisitBioseqsInSep (sep, &disc_pairs, FindInconsistentSourceAndDeflineCallback);

  if (disc_pairs == NULL) 
  {
    return;
  }
  else if (disc_pairs->next == NULL)
  {
    ValNodeLink (discrepancy_list, disc_pairs);
  }
  else
  {
    dip = NewDiscrepancyItem (DISC_INCONSISTENT_BIOSRC_DEFLINE, bad_fmt, disc_pairs);
    dip->item_list = NULL;
    dip->subcategories = disc_pairs;
    
    ValNodeAddPointer (discrepancy_list, 0, dip);
  }
}


static void SetDiscrepancyLevels (ValNodePtr discrepancy_list, Int4 level)
{
  DiscrepancyItemPtr dip;
  
  while (discrepancy_list != NULL)
  {
    dip = (DiscrepancyItemPtr) discrepancy_list->data.ptrvalue;
    if (dip != NULL)
    {
      dip->level = level;
      SetDiscrepancyLevels (dip->subcategories, level + 1);
    }
    discrepancy_list = discrepancy_list->next;
  }
}


/* Note that this function contains a hack - it assumes that all of the
 * test types that use the same collection function are listed together.
 */
static ValNodePtr CollectDiscrepancies (SeqEntryPtr sep, DiscrepancyConfigPtr dcp)
{
  ValNodePtr             discrepancy_list = NULL;
  Int4                   i;
  PerformDiscrepancyTest last_test_func = NULL;

  for (i = 0; i < MAX_DISC_TYPE; i++)
  {
    if ((dcp == NULL || dcp->conf_list[i])
        && discrepancy_info_list[i].test_func != NULL
        && discrepancy_info_list[i].test_func != last_test_func)
    {
      discrepancy_info_list[i].test_func (sep, &discrepancy_list);
      last_test_func = discrepancy_info_list[i].test_func;
    }
  }
  
  /* because some tests are run together, need to remove unwanted results */
  RemoveUnwantedDiscrepancyItems (&discrepancy_list, dcp);

  /* normalize the discrepancy levels so that they will be correctly displayed */
  SetDiscrepancyLevels (discrepancy_list, 0);
  return discrepancy_list;  
}

static Nlm_ParData discParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static Nlm_ColData discColFmt[2] = {{16, 0, 0, 0, NULL, 'l', 0,0,0,0, FALSE},
                                    {1000, 0, 0, 0, NULL, 'l', 1,0,0,0, TRUE}};


static Int4 CountChosenDiscrepancies (ValNodePtr discrepancy_list, Boolean count_all)
{
  Int4               num_chosen = 0;
  DiscrepancyItemPtr dip;
  
  while (discrepancy_list != NULL)
  {
    dip = (DiscrepancyItemPtr) discrepancy_list->data.ptrvalue;
    if (dip != NULL)
    {
      if (dip->chosen || count_all)
      {
        if (dip->expanded)
        {
          num_chosen += CountChosenDiscrepancies (dip->subcategories, TRUE);
        }
        else
        {
          num_chosen ++;
        }
      }
      else if (dip->expanded)
      {
        num_chosen += CountChosenDiscrepancies (dip->subcategories, FALSE);
      }
    }
    discrepancy_list = discrepancy_list->next;
  }
  return num_chosen;
}

static Int4 CountLevels (ValNodePtr discrepancy_list)
{
  Int4       num_levels = 1, num, num_sublevels = 0;
  ValNodePtr vnp;
  DiscrepancyItemPtr dip;
  
  if (discrepancy_list == NULL) 
  {
    return 0;
  }
  
  for (vnp = discrepancy_list; vnp != NULL; vnp = vnp->next)
  {
    dip = (DiscrepancyItemPtr) vnp->data.ptrvalue;
    if (dip == NULL || dip->subcategories == NULL || !dip->expanded)
    {
      continue;
    }
    num = CountLevels (dip->subcategories);
    if (num > num_sublevels) num_sublevels = num;
  }
  
  /* one level for the top plus levels for the subcategories */

  return 1 + num_sublevels;
}


static Nlm_ColPtr PNTR GetColumnFormatArrays (Int4 num_levels, DoC doc)
{
  Int4               n, k;
  Nlm_ColPtr PNTR    col_fmt_array_array = NULL;
  RecT               r;
  Int4               doc_width;
    
  if (num_levels == 0)
  {
    return NULL;
  }
  
  ObjectRect (doc, &r);
  InsetRect (&r, 4, 4);
  doc_width = r.right - r.left;
  
  col_fmt_array_array = (Nlm_ColPtr PNTR) MemNew (sizeof (Nlm_ColPtr) * num_levels);
  for (n = 0; n < num_levels; n++)
  {
    col_fmt_array_array[n] = (Nlm_ColPtr) MemNew (sizeof (Nlm_ColData) * (n + 3));
    for (k = 0; k < n + 2; k++)
    {
      col_fmt_array_array[n][k].pixWidth = 16;
      col_fmt_array_array[n][k].pixInset = 0;
      col_fmt_array_array[n][k].charWidth = 0;
      col_fmt_array_array[n][k].charInset = 0;
      col_fmt_array_array[n][k].font = NULL;
      col_fmt_array_array[n][k].just = 'l';
      col_fmt_array_array[n][k].wrap = 0;
      col_fmt_array_array[n][k].bar = 0;
      col_fmt_array_array[n][k].underline = 0;
      col_fmt_array_array[n][k].left = 0;
      col_fmt_array_array[n][k].last = 0;
    }
    col_fmt_array_array[n][k].pixWidth = doc_width - ((n + 2) * 16);
    col_fmt_array_array[n][k].pixInset = 0;
    col_fmt_array_array[n][k].charWidth = 0;
    col_fmt_array_array[n][k].charInset = 0;
    col_fmt_array_array[n][k].font = NULL;
    col_fmt_array_array[n][k].just = 'l';
    col_fmt_array_array[n][k].wrap = 1;
    col_fmt_array_array[n][k].bar = 0;
    col_fmt_array_array[n][k].underline = 0;
    col_fmt_array_array[n][k].left = 0;
    col_fmt_array_array[n][k].last = 1;
  }
  return col_fmt_array_array;
}


static Nlm_ColPtr PNTR FreeColumnFormatArrays (Nlm_ColPtr PNTR col_fmt_array_array, Int4 num_levels)
{
  Int4 n;
  
  if (col_fmt_array_array == NULL || num_levels < 1)
  {
    return NULL;
  }
  for (n = 0; n < num_levels; n++)
  {
    col_fmt_array_array [n] = MemFree (col_fmt_array_array [n]);
  }
  col_fmt_array_array = MemFree (col_fmt_array_array);
  return col_fmt_array_array;
}

static void AddDiscrepancyItem (DiscrepancyReportFormPtr drfp, DiscrepancyItemPtr dip, Int4 level)
{
  CharPtr            item_text;
  ValNodePtr         vnp;
  Int4               n;

  if (dip == NULL)
  {
    return;
  }
  item_text = (CharPtr) MemNew (sizeof (Char) * (StringLen (dip->description) + 6 + level));
  for (n = 0; n < level; n++)
  {
    StringCat (item_text, "\t");
  }
  StringCat (item_text, " \t \t");
  StringCat (item_text, dip->description);
  StringCat (item_text, "\n");
  AppendText (drfp->doc, item_text, &discParFmt, drfp->col_fmt_array_array [level], programFont);
  if (dip->expanded)
  {
    for (vnp = dip->subcategories; vnp != NULL; vnp = vnp->next)
    {
      AddDiscrepancyItem (drfp, vnp->data.ptrvalue, level + 1);
    }
  }
}


static void PopulateDiscrepancyList (DiscrepancyReportFormPtr drfp, ValNodePtr discrepancy_list)
{
  Int2               numItems;
  Int4               num_levels;
  
  if (drfp == NULL || drfp->doc == NULL) 
  {
    return;
  }
  
  Reset (drfp->doc);
  
  num_levels = CountLevels (drfp->discrepancy_list);
  if (num_levels != drfp->num_levels)
  {
    drfp->col_fmt_array_array = FreeColumnFormatArrays (drfp->col_fmt_array_array, drfp->num_levels);
    drfp->num_levels = num_levels;
    drfp->col_fmt_array_array = GetColumnFormatArrays (drfp->num_levels, drfp->doc);
  }
  
  while (discrepancy_list != NULL)
  {
    AddDiscrepancyItem (drfp, discrepancy_list->data.ptrvalue, 0);
    discrepancy_list = discrepancy_list->next;
  }
  GetDocParams (drfp->doc, &numItems, NULL);
  UpdateDocument (drfp->doc, 0, numItems);

}

static Nlm_ParData discItemParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static Nlm_ColData discItemColFmt [3] = {{0, 5, 10, 0, NULL, 'l', 1,0,0,0, FALSE},
                                         {0, 0, 10, 0, NULL, 'l', 1,0,0,0, FALSE},
                                         {0, 0, 10, 0, NULL, 'l', 1,0,0,0, TRUE}};


static CharPtr GetDiscrepancyItemText (ValNodePtr vnp)
{
  CharPtr           row_text = NULL;
  SeqFeatPtr        sfp;
  BioseqPtr         bsp;
  SeqMgrFeatContext context;
  CharPtr           location;
  CharPtr           label;
  SeqDescrPtr       sdp;
  
  if (vnp == NULL)
  {
    return NULL;
  }
  if (vnp->choice == OBJ_SEQFEAT)
  {
    sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    if (sfp != NULL)
    {
      sfp = SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, 0, 0, sfp, &context);
      if (sfp != NULL)
      {
        location = SeqLocPrintUseBestID (sfp->location);
        label = (CharPtr) FeatDefTypeLabel(sfp);
        row_text = (CharPtr) MemNew (sizeof (Char) * 
                                     (StringLen (label) 
                                      + StringLen (context.label) 
                                      + StringLen (location) 
                                      + 4));
        sprintf (row_text, "%s\t%s\t%s\n", label, context.label, location);
        location = MemFree (location);
      }
    }
  }
  else if (vnp->choice == OBJ_BIOSEQ)
  {
    bsp = (BioseqPtr) vnp->data.ptrvalue;
    if (bsp != NULL)
    {
      row_text = (CharPtr) MemNew (sizeof (Char) * 41);
      SeqIdWrite (bsp->id, row_text, PRINTID_REPORT, 39);
      StringCat (row_text, "\n");
    }
  }
  else if (vnp->choice == OBJ_SEQDESC)
  {
    sdp = (SeqDescrPtr) vnp->data.ptrvalue;
    if (sdp != NULL)
    {
      row_text = (CharPtr) MemNew (sizeof (Char) * 61);
      SeqDescLabel (sdp, row_text, 59, TRUE);
      StringCat (row_text, "\n");
    }
  }
  
    
  return row_text;
}

static void PopulateDiscrepancyItemList (DoC doc, DiscrepancyItemPtr dip)
{
  ValNodePtr        vnp;
  Int2              numItems;
  CharPtr           row_text;
  RecT              r;
  
  if (doc == NULL)
  {
    return;
  }
  Reset (doc);
  
  if (dip == NULL)
  {
    return;
  }
  
  if (dip->item_list == NULL)
  {
    AppendText (doc, "No items listed", NULL, NULL, programFont);
  }
  
  ObjectRect (doc, &r);
  InsetRect (&r, 4, 4);
  
  discItemColFmt[0].pixWidth = 5 * stdCharWidth;
  discItemColFmt[1].pixWidth = (r.right - r.left - discItemColFmt[0].pixWidth) / 2;
  discItemColFmt[2].pixWidth = (r.right - r.left - discItemColFmt[0].pixWidth) / 2;
  
  vnp = dip->item_list;
  
  while (vnp != NULL)
  {
    row_text = GetDiscrepancyItemText (vnp);
    if (row_text != NULL)
    {
      if (vnp->choice == OBJ_SEQFEAT)
      {
        AppendText (doc, row_text, &discItemParFmt, discItemColFmt, programFont);
      }
      else
      {
        AppendText (doc, row_text, &discItemParFmt, NULL, programFont);
      }
      row_text = MemFree (row_text);
    }
    vnp = vnp->next;
  }
  GetDocParams (doc, &numItems, NULL);
  UpdateDocument (doc, 0, numItems);  
}

static DiscrepancyItemPtr GetSubDiscrepancy (ValNodePtr discrepancy_list, Int2Ptr pitem)
{
  DiscrepancyItemPtr dip = NULL;

  if (discrepancy_list == NULL || pitem == NULL)
  {
    return NULL;
  }
  while (*pitem > 0 && discrepancy_list != NULL)
  {
    (*pitem)--;
    dip = (DiscrepancyItemPtr) discrepancy_list->data.ptrvalue;
    if (*pitem > 0)
    {
      if (dip != NULL && dip->expanded)
      {
        dip = GetSubDiscrepancy (dip->subcategories, pitem);
      }
    }
    discrepancy_list = discrepancy_list->next;
  }
  if (*pitem > 0)
  {
    dip = NULL;
  }
  return dip;
}

static DiscrepancyItemPtr GetSelectedDiscrepancy (ValNodePtr discrepancy_list, Int2 item)
{
  DiscrepancyItemPtr dip = NULL;
  
  dip = GetSubDiscrepancy (discrepancy_list, &item);

  return dip;
}


static void ClickDiscrepancy (DoC d, PoinT pt)

{
  Int2                      item, numItems;
  Int2                      row;
  Int2                      col;
  DiscrepancyReportFormPtr  drfp;
  DiscrepancyItemPtr        dip;

  drfp = GetObjectExtra (d);
  if (drfp != NULL) {
    MapDocPoint (d, pt, &item, &row, &col, NULL);
    if (item > 0 && row > 0 && drfp->clicked == item) {
      drfp->dblClick = dblClick;
    } else {
      drfp->dblClick = FALSE;
    }
    drfp->clicked = 0;
    if (item > 0 && row > 0) {
      drfp->clicked = item;
    }
    if (item > 0 && row > 0 && !dblClick)
    {
      dip = GetSelectedDiscrepancy (drfp->discrepancy_list, item);
      if (dip != NULL)
      {
        if (col == dip->level + 1)
        {
          dip->chosen = !dip->chosen;
          GetDocParams (d, &numItems, NULL);
          UpdateDocument (d, 0, numItems);
        }
        else if (col == dip->level + 2)
        {
          dip->expanded = !dip->expanded;
          PopulateDiscrepancyList (drfp, drfp->discrepancy_list);
        }
      }
    }
  }
}


static void ActOnDiscrepancy (ValNodePtr discrepancy_list, Int2 item)
{
  DiscrepancyItemPtr dip;
  
  dip = GetSelectedDiscrepancy (discrepancy_list, item);
  if (dip != NULL && dip->callback_func != NULL)
  {
    (dip->callback_func) (dip->item_list, dip->callback_data);
  }
}

static void InvalBorder (DoC d, Int2 item)

{
  Int2  bottom;
  RecT  r;
  Int2  top;

  ObjectRect (d, &r);
  InsetRect (&r, 4, 4);
  if (ItemIsVisible (d, item, &top, &bottom, NULL)) {
    r.top = top;
    r.bottom = bottom;
    r.right = r.left + 4;
    InsetRect (&r, -1, -1);
    InvalRect (&r);
  }
}

static void ReleaseDiscrepancy (DoC d, PoinT pt)

{
  Int2           item;
  Int2           old;
  Int2           row;
  DiscrepancyReportFormPtr  drfp;

  drfp = GetObjectExtra (d);
  if (drfp != NULL) {
    ResetClip ();
    MapDocPoint (d, pt, &item, &row, NULL, NULL);
    if (item > 0 && row > 0) {
      if (item == drfp->clicked) {
        old = drfp->selected;
        drfp->selected = item;
        if (old != item) {
          if (old == 0) {
            UpdateDocument (d, item, item);
          } else {
            UpdateDocument (d, old, old);
            UpdateDocument (d, item, item);
          }
          Update ();
        }
      }
    } else if (drfp->clicked == 0) {
      if (drfp->selected != 0) {
        old = drfp->selected;
        drfp->selected = 0;
        InvalBorder (d, old);
      }
      Update ();
    }
    if (drfp->selected > 0 && drfp->dblClick)
    {
      ActOnDiscrepancy (drfp->discrepancy_list, drfp->selected);
    }
    else if (drfp->selected > 0)
    {
      drfp->item_selected = 0;
      PopulateDiscrepancyItemList (drfp->item_list, 
                                   GetSelectedDiscrepancy (drfp->discrepancy_list,
                                                           drfp->selected));
      
    }
  }
}


static void DrawDiscrepancy (DoC d, RectPtr r, Int2 item, Int2 firstLine)

{
  DiscrepancyReportFormPtr drfp;
  RecT                     rct;
  DiscrepancyItemPtr       dip;
  Int4                     level_offset;

  drfp = (DiscrepancyReportFormPtr) GetObjectExtra (d);
  if (drfp != NULL && r != NULL && item > 0 && firstLine == 0) {
    rct = *r;
  
    dip = GetSelectedDiscrepancy (drfp->discrepancy_list, item);
    if (dip != NULL)
    {
      level_offset = dip->level * 16;
      rct.left += level_offset;
      rct.right += level_offset;
    }
  
    /* draw selection */
    if (item == drfp->selected) {
      rct.right = rct.left + 4;
      PaintRect (&rct);
    }

    /* draw chosen checkboxes */
    rct.left += 5;
    rct.right = rct.left + 10;
    rct.bottom = rct.top + (rct.right - rct.left);
    FrameRect (&rct);
    
    if (dip != NULL && dip->chosen) {
      MoveTo (rct.left, rct.top);
      LineTo (rct.right - 1, rct.bottom - 1);
      MoveTo (rct.left, rct.bottom - 1);
      LineTo (rct.right - 1, rct.top);
    }
    
    /* draw open/closed checkboxes */
    if (dip!= NULL && dip->subcategories != NULL)
    {
      rct.left += 10;
      rct.right = rct.left + 10;
      rct.bottom = rct.top + (rct.right - rct.left);
      FrameRect (&rct);
      MoveTo (rct.left, (rct.top + rct.bottom) / 2);
      LineTo (rct.right - 1, (rct.top + rct.bottom) / 2);
      if (!dip->expanded)
      {
        MoveTo ((rct.left + rct.right) / 2, rct.top);
        LineTo ((rct.left + rct.right) / 2, rct.bottom - 1);
      }
    }
  }
}



static void DrawDiscrepancyItem (DoC d, RectPtr r, Int2 item, Int2 firstLine)

{
  DiscrepancyReportFormPtr drfp;
  RecT                     rct;

  drfp = (DiscrepancyReportFormPtr) GetObjectExtra (d);
  if (drfp != NULL && r != NULL && item > 0 && firstLine == 0) {
    rct = *r;
  
    /* draw selection */
    if (item == drfp->item_selected) {
      rct = *r;
      rct.right = rct.left + 4;
      PaintRect (&rct);
    }
  }
}


static void EditDiscrepancyItem (ValNodePtr vnp)
{
  SeqFeatPtr sfp;
  BioseqPtr  bsp;
  SeqDescrPtr sdp;
  ObjValNodePtr ovp;
  
  if (vnp == NULL)
  {
    return;
  }
  if (vnp->choice == OBJ_SEQFEAT)
  {
    sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    if (sfp != NULL)
    {
      GatherProcLaunch (OMPROC_EDIT, FALSE, sfp->idx.entityID, sfp->idx.itemID,
                        OBJ_SEQFEAT, 0, 0, OBJ_SEQFEAT, 0);
    }
  }
  else if (vnp->choice == OBJ_BIOSEQ)
  {
    bsp = (BioseqPtr) vnp->data.ptrvalue;
    if (bsp != NULL)
    {
      GatherProcLaunch (OMPROC_EDIT, FALSE, bsp->idx.entityID, bsp->idx.itemID,
                         OBJ_BIOSEQ, 0, 0, OBJ_BIOSEQ, 0);
    }
  }
  else if (vnp->choice == OBJ_SEQDESC)
  {
    sdp = (SeqDescrPtr) (vnp->data.ptrvalue);
    if (sdp != NULL && sdp->extended != 0)
    {
      ovp = (ObjValNodePtr) sdp;
      GatherProcLaunch (OMPROC_EDIT, FALSE, ovp->idx.entityID, ovp->idx.itemID,
                         OBJ_SEQDESC, 0, 0, OBJ_SEQDESC, 0);
    }
  }

}


static void SetBioseqViewTargetByBioseq (BaseFormPtr bfp, BioseqPtr bsp)
{
  Char       id_text [41];
  
  if (bsp != NULL && bfp != NULL)
  {
    SeqIdWrite (SeqIdFindBest (bsp->id, 0), id_text, PRINTID_REPORT, sizeof (id_text));
    SetBioseqViewTarget (bfp, id_text);
  }
}


static BioseqPtr GetFirstBioseqInSeqEntry (SeqEntryPtr sep)
{
  BioseqPtr    bsp = NULL;
  BioseqSetPtr bssp;
  
  if (sep == NULL || sep->data.ptrvalue == NULL)
  {
    return NULL;
  }
  else if (IS_Bioseq (sep))
  {
    bsp = sep->data.ptrvalue;
  }
  else if (IS_Bioseq_set (sep))
  {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    for (sep = bssp->seq_set; sep != NULL && bsp == NULL; sep = sep->next)
    {
      bsp = GetFirstBioseqInSeqEntry (sep);
    }
  }
  return bsp;
}


static BioseqPtr GetBioseqForDescriptor (ObjValNodePtr ovp)
{
  BioseqPtr    bsp = NULL;
  BioseqSetPtr bssp;
  SeqEntryPtr  sep;
  
  if (ovp == NULL || ovp->idx.parentptr == NULL)
  {
    return NULL;
  }
  else if (ovp->idx.parenttype == OBJ_BIOSEQ) {
    bsp = (BioseqPtr) ovp->idx.parentptr;
  } else if (ovp->idx.parenttype == OBJ_BIOSEQSET) {
    bssp = (BioseqSetPtr) ovp->idx.parentptr;
    for (sep = bssp->seq_set; sep != NULL && bsp == NULL; sep = sep->next)
    {
      bsp = GetFirstBioseqInSeqEntry (sep);
    }  
  }
  return bsp;
}

static void ScrollToDiscrepancyItem (ValNodePtr vnp, BaseFormPtr bfp)
{
  SeqFeatPtr    sfp;
  BioseqPtr     bsp;
  SeqDescrPtr   sdp;
  ObjValNodePtr ovp;
  
  if (vnp == NULL || bfp == NULL)
  {
    return;
  }
  if (vnp->choice == OBJ_SEQFEAT)
  {
    sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    if (sfp != NULL)
    {
      /* need to scroll to item */
      bsp = BioseqFindFromSeqLoc (sfp->location);
      SetBioseqViewTargetByBioseq (bfp, bsp);
      ObjMgrSelect (sfp->idx.entityID, sfp->idx.itemID, OBJ_SEQFEAT, 0, NULL);
    }
  }
  else if (vnp->choice == OBJ_BIOSEQ)
  {
    bsp = (BioseqPtr) vnp->data.ptrvalue;
    SetBioseqViewTargetByBioseq (bfp, bsp);
  }
  else if (vnp->choice == OBJ_SEQDESC)
  {
    sdp = (SeqDescrPtr) (vnp->data.ptrvalue);
    if (sdp != NULL && sdp->extended != 0)
    {
      ovp = (ObjValNodePtr) sdp;      
      bsp = GetBioseqForDescriptor (ovp);
      SetBioseqViewTargetByBioseq (bfp, bsp);
      ObjMgrSelect (ovp->idx.entityID, ovp->idx.itemID, OBJ_SEQDESC, 0, NULL);
    }
  }
}

static void ClickDiscrepancyItem (DoC d, PoinT pt)

{
  Int2                      item, last_selected, numItems;
  Int2                      row;
  DiscrepancyReportFormPtr  drfp;
  DiscrepancyItemPtr        dip;
  ValNodePtr                vnp;

  drfp = GetObjectExtra (d);
  if (drfp != NULL) {
    MapDocPoint (d, pt, &item, &row, NULL, NULL);
    if (item > 0 && row > 0) {  
      dip = GetSelectedDiscrepancy (drfp->discrepancy_list, drfp->selected);
      if (dip != NULL && dip->item_list != NULL)
      {
        vnp = dip->item_list;
        
        last_selected = drfp->item_selected;
        drfp->item_selected = item;
        
        if (item != last_selected)
        {
          GetDocParams (d, &numItems, NULL);
          UpdateDocument (d, 0, numItems);
        }
    
        /* find item in list */
        while (item > 1 && vnp != NULL)
        {
          vnp = vnp->next;
          item--;
        }
        
        if (dblClick)
        {
          EditDiscrepancyItem (vnp);
        }
        else
        {
          ScrollToDiscrepancyItem (vnp, drfp->bfp);
        } 
         
      }
    }
  }
}


static void CleanupDiscrepancyReportForm (GraphiC g, VoidPtr data)

{
  DiscrepancyReportFormPtr drfp;

  drfp = (DiscrepancyReportFormPtr) data;
  if (drfp != NULL) {
    drfp->discrepancy_list = FreeDiscrepancyList (drfp->discrepancy_list);
    drfp->dcp = DiscrepancyConfigFree (drfp->dcp);
    drfp->col_fmt_array_array = FreeColumnFormatArrays (drfp->col_fmt_array_array, drfp->num_levels);
    ObjMgrFreeUserData (drfp->input_entityID, drfp->procid, drfp->proctype, drfp->userkey);
    discrepancyReportWindow = NULL;
  }
  StdCleanupFormProc (g, data);
}

static void RecheckDiscrepancyProc (ButtoN b)
{
  DiscrepancyReportFormPtr drfp;

  drfp = (DiscrepancyReportFormPtr) GetObjectExtra (b);
  if (drfp != NULL)
  {
    drfp->discrepancy_list = FreeDiscrepancyList (drfp->discrepancy_list);
    drfp->discrepancy_list = CollectDiscrepancies (drfp->sep, drfp->dcp);
    PopulateDiscrepancyList (drfp, drfp->discrepancy_list);
    if (drfp->discrepancy_list != NULL)
    {
      drfp->selected = 1;
      drfp->item_selected = 0;
      PopulateDiscrepancyItemList (drfp->item_list, 
                                   GetSelectedDiscrepancy (drfp->discrepancy_list,
                                                           drfp->selected));
    }
    else
    {
      Reset (drfp->item_list);
    }
  }
}


static ValNodePtr ValNodePointerDup (ValNodePtr vnp)
{
  ValNodePtr vnp_new = NULL;
  
  if (vnp != NULL)
  {
    vnp_new = ValNodeNew (NULL);
    vnp_new->choice = vnp->choice;
    vnp_new->data.ptrvalue = vnp->data.ptrvalue;
    vnp_new->next = ValNodePointerDup (vnp->next);
  }
  return vnp_new;
}

static ValNodePtr ReplaceDiscrepancyItemWithFeatureTableStrings (ValNodePtr feat_list)
{
  BioseqPtr       bsp;
  CstType         custom_flags = 0;
  Asn2gbJobPtr    ajp;
  BaseBlockPtr    bbp;
  Int4            index;
  SeqFeatPtr      sfp;
  ValNodePtr      vnp, list_copy = NULL, list_vnp;
  
  if (feat_list == NULL) return NULL;
  
  list_copy = ValNodePointerDup (feat_list);
  for (vnp = list_copy; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == OBJ_SEQFEAT)
    {
      sfp = (SeqFeatPtr) vnp->data.ptrvalue;
      bsp = BioseqFindFromSeqLoc (sfp->location);
      ajp = asn2gnbk_setup (bsp, NULL, NULL, FTABLE_FMT, DUMP_MODE, NORMAL_STYLE,
                             0, 0, custom_flags, NULL);
      if (ajp != NULL)
      {
        for (index = 0; index < ajp->numParagraphs; index++) 
        {
          bbp = ajp->paragraphArray [index];
          for (list_vnp = vnp; list_vnp != NULL; list_vnp = list_vnp->next)
          {
            if (list_vnp->choice == OBJ_SEQFEAT)
            {
              sfp = (SeqFeatPtr) list_vnp->data.ptrvalue;
              if (sfp != NULL 
                  && bbp->entityID == sfp->idx.entityID
                  && bbp->itemtype == sfp->idx.itemtype
                  && bbp->itemID == sfp->idx.itemID)
              {
                /* replace list feature with description, change choice */
                list_vnp->choice = 0;
                list_vnp->data.ptrvalue = asn2gnbk_format (ajp, (Int4) index);
              }
            }
          }
        }
        asn2gnbk_cleanup (ajp);
      }
    }
  }
  return list_copy;
}


static void WriteDiscrepancy (FILE *fp, DiscrepancyItemPtr dip, Boolean use_feature_table_fmt)
{
  ValNodePtr vnp, list_copy = NULL;
  CharPtr    row_text;
  
  if (fp == NULL || dip == NULL)
  {
    return;
  }
  
  fprintf (fp, "%s\n", dip->description);
  vnp = dip->item_list;
  
  if (use_feature_table_fmt)
  {
    list_copy = ReplaceDiscrepancyItemWithFeatureTableStrings (vnp);
    vnp = list_copy;
  }

  while (vnp != NULL)
  {
    if (vnp->choice == 0 && use_feature_table_fmt)
    {
      row_text = vnp->data.ptrvalue;
    }
    else
    {
      row_text = GetDiscrepancyItemText (vnp);
    }
    if (row_text != NULL)
    {
      fprintf (fp, row_text);
      row_text = MemFree (row_text);
    }
    vnp = vnp->next;
  }
  
  fprintf (fp, "\n");
}

static void 
WriteDiscrepancyReport 
(FILE       *fp,
 ValNodePtr discrepancy_list, 
 Boolean    show_all,
 Boolean    use_feature_table_fmt)
{
  DiscrepancyItemPtr       dip;
  ValNodePtr               vnp;
  Int4                     num_chosen;

  if (fp == NULL || discrepancy_list == NULL)
  {
    return;
  }
  for (vnp = discrepancy_list; vnp != NULL; vnp = vnp->next)
  {
    dip = (DiscrepancyItemPtr) vnp->data.ptrvalue;
    if (dip != NULL)
    {
      if (dip->expanded)
      {
        num_chosen = CountChosenDiscrepancies (dip->subcategories, show_all | dip->chosen);
        if (num_chosen > 0)
        {
          if (dip->chosen || show_all)
          {
            fprintf (fp, "%s\n", dip->description);
          }
          WriteDiscrepancyReport (fp, dip->subcategories, show_all | dip->chosen, use_feature_table_fmt);          
        }
      }
      else if (dip->chosen || show_all)
      {
        WriteDiscrepancy (fp, dip, use_feature_table_fmt);
      }
    }
  }
}

static Boolean DiscrepancyReportExportProc (ForM f, CharPtr filename)

{
  FILE           *fp;
  Char           path [PATH_MAX];
  DiscrepancyReportFormPtr drfp;
  Int4                     num_disc = 0;
  Boolean                  show_all = FALSE;

  drfp = (DiscrepancyReportFormPtr) GetObjectExtra (f);
  if (drfp == NULL) 
  {
    return FALSE;
  }
  
  num_disc = CountChosenDiscrepancies (drfp->discrepancy_list, FALSE);

  if (num_disc == 0) 
  {
    if (ANS_CANCEL == Message (MSG_OKC, "No discrepancies selected!  Export all?"))
    {
      return FALSE;
    }
    else
    {
      show_all = TRUE;
    }
  }
  
  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  if (path [0] != '\0' || GetOutputFileName (path, sizeof (path), NULL)) {
#ifdef WIN_MAC
    fp = FileOpen (path, "r");
    if (fp != NULL) {
      FileClose (fp);
    } else {
      FileCreate (path, "TEXT", "ttxt");
    }
#endif
    fp = FileOpen (path, "w");
    if (fp != NULL) {
      WriteDiscrepancyReport (fp, drfp->discrepancy_list, show_all, 
                              (Boolean)(drfp->dcp != NULL && drfp->dcp->use_feature_table_format));
      FileClose (fp);
      return TRUE;
    }
  }
  return FALSE;
}


static void GenerateDiscrepancyReport (ButtoN b)
{
  DiscrepancyReportFormPtr drfp;
  Char                     path [PATH_MAX];

  drfp = (DiscrepancyReportFormPtr) GetObjectExtra (b);
  if (drfp == NULL)
  {
    return;
  }

  TmpNam (path);  
  if (DiscrepancyReportExportProc (drfp->form, path))
  {
    LaunchGeneralTextViewer (path, "Discrepancy Report");
  }
  FileRemove (path);  
}


static void CopyDiscrepancyReportToClipboard (DiscrepancyReportFormPtr drfp)

{
  FILE         *fp;
  Char         path [PATH_MAX];

  if (drfp == NULL || drfp->discrepancy_list == NULL) return;
  TmpNam (path);
  fp = FileOpen (path, "w");
  if (fp != NULL) {
    WriteDiscrepancyReport (fp, drfp->discrepancy_list, TRUE, 
                            (Boolean)(drfp->dcp != NULL && drfp->dcp->use_feature_table_format));
    FileClose (fp);
    FileToClipboard (path);
  }
  FileRemove (path);
}


static void DiscrepancyReportFormMessage (ForM f, Int2 mssg)

{
  DiscrepancyReportFormPtr drfp;

  drfp = (DiscrepancyReportFormPtr) GetObjectExtra (f);
  if (drfp != NULL) {
    switch (mssg) {
      case VIB_MSG_EXPORT :
        DiscrepancyReportExportProc (f, NULL);
        break;
      case VIB_MSG_PRINT :
        break;
      case VIB_MSG_CLOSE :
        Remove (f);
        break;
      case VIB_MSG_CUT :
        CopyDiscrepancyReportToClipboard (drfp);
        break;
      case VIB_MSG_COPY :
        CopyDiscrepancyReportToClipboard (drfp);
        break;
      case VIB_MSG_PASTE :
        break;
      case VIB_MSG_DELETE :
        drfp->discrepancy_list = ValNodeFree (drfp->discrepancy_list);
        Reset (drfp->doc);
        Reset (drfp->item_list);
        break;
      default :
        if (drfp->appmessage != NULL) {
          drfp->appmessage (f, mssg);
        }
        break;
    }
  }
}


static void ReactivateDiscrepancyReport (BaseFormPtr bfp)
{
  DiscrepancyReportFormPtr drfp;

  if (discrepancyReportWindow == NULL) 
  {
    CreateDiscrepancyReportWindow (bfp);
  }
  
  drfp = (DiscrepancyReportFormPtr) GetObjectExtra (discrepancyReportWindow);
  if (drfp == NULL)
  {
    Remove (discrepancyReportWindow);
    discrepancyReportWindow = NULL;
    CreateDiscrepancyReportWindow (bfp);
  }
  
  drfp->bfp = bfp;
  drfp->input_entityID = bfp->input_entityID;
  drfp->sep = GetTopSeqEntryForEntityID (drfp->input_entityID);
  
  /* populate discrepancy lists */
  RecheckDiscrepancyProc (drfp->recheck_btn);
  Show (discrepancyReportWindow);  
}


static void EditDiscrepancyConfigBtn (ButtoN b)
{
  DiscrepancyReportFormPtr drfp;
  
  drfp = (DiscrepancyReportFormPtr) GetObjectExtra (discrepancyReportWindow);
  if (drfp == NULL) return;
  
  if (EditDiscrepancyConfig (drfp->dcp))
  {
    RecheckDiscrepancyProc (b);
  }
}


#ifndef WIN_MAC
extern void CreateStdValidatorFormMenus (WindoW w);
#endif


static Int2 LIBCALLBACK DiscrepancyReportMsgFunc (OMMsgStructPtr ommsp)
{
  WindoW                   currentport,
                           temport;
  OMUserDataPtr            omudp;
  DiscrepancyReportFormPtr drfp = NULL;
  
  omudp = (OMUserDataPtr)(ommsp->omuserdata);
  if (omudp == NULL) return OM_MSG_RET_ERROR;
  drfp = (DiscrepancyReportFormPtr) omudp->userdata.ptrvalue;
  if (drfp == NULL) return OM_MSG_RET_ERROR;

  currentport = ParentWindow (drfp->form);
  temport = SavePort (currentport);
  UseWindow (currentport);
  Select (drfp->form);
  switch (ommsp->message) 
  {
      case OM_MSG_UPDATE:
          break;
      case OM_MSG_DESELECT:
          break;

      case OM_MSG_SELECT: 
          break;
      case OM_MSG_DEL:
          Remove (drfp->form);
          break;
      case OM_MSG_HIDE:
          break;
      case OM_MSG_SHOW:
          break;
      case OM_MSG_FLUSH:
          Remove (drfp->form);	
          break;
      default:
          break;
  }
  RestorePort (temport);
  UseWindow (temport);
  return OM_MSG_RET_OK;
}


extern void CreateDiscrepancyReportWindow (BaseFormPtr bfp)
{
  DiscrepancyReportFormPtr drfp;
  GrouP                    h, g;
  ButtoN                   b;
  GrouP                    c;
  RecT                     r;
  WindoW                   w;
  OMUserDataPtr            omudp;

  if (bfp == NULL) return;
  
  if (discrepancyReportWindow != NULL)
  {
    ReactivateDiscrepancyReport (bfp);
    return; 
  }
  
  drfp = (DiscrepancyReportFormPtr) MemNew (sizeof (DiscrepancyReportFormData));
  if (drfp == NULL)
  {
    return;
  }
  
  drfp->bfp = bfp;
  drfp->input_entityID = bfp->input_entityID;
  drfp->sep = GetTopSeqEntryForEntityID (drfp->input_entityID);
  w = FixedWindow (-50, -33, -10, -10, "Discrepancy Report", StdCloseWindowProc);
  SetObjectExtra (w, drfp, CleanupDiscrepancyReportForm);
  drfp->form = (ForM) w;
  drfp->formmessage = DiscrepancyReportFormMessage;
  drfp->exportform = DiscrepancyReportExportProc;
  
  /* read in config file */
  drfp->dcp = ReadDiscrepancyConfig();
  
  /* register to receive update messages */
  drfp->userkey = OMGetNextUserKey ();
  drfp->procid = 0;
  drfp->proctype = OMPROC_EDIT;
  omudp = ObjMgrAddUserData (drfp->input_entityID, drfp->procid, drfp->proctype, drfp->userkey);
  if (omudp != NULL) {
    omudp->userdata.ptrvalue = (Pointer) drfp;
    omudp->messagefunc = DiscrepancyReportMsgFunc;
  }


#ifndef WIN_MAC
  CreateStdValidatorFormMenus (w);
#endif

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  g = HiddenGroup (h, 2, 0, NULL);
  StaticPrompt (g, "Discrepancies", 0, popupMenuHeight, programFont, 'c');
  StaticPrompt (g, "Affected Items", 0, popupMenuHeight, programFont, 'c');
  
  drfp->doc = DocumentPanel (g, stdCharWidth * 30, stdLineHeight * 20);
  SetObjectExtra (drfp->doc, drfp, NULL);
  SetDocAutoAdjust (drfp->doc, FALSE);
  SetDocProcs (drfp->doc, ClickDiscrepancy, NULL, ReleaseDiscrepancy, NULL);
  SetDocShade (drfp->doc, DrawDiscrepancy, NULL, NULL, NULL);
  
  drfp->item_list = DocumentPanel (g, stdCharWidth * 30 + 5, stdLineHeight * 20);
  SetObjectExtra (drfp->item_list, drfp, NULL);
  SetDocAutoAdjust (drfp->item_list, FALSE);
  SetDocProcs (drfp->item_list, ClickDiscrepancyItem, NULL, NULL, NULL);
  SetDocShade (drfp->item_list, DrawDiscrepancyItem, NULL, NULL, NULL);

  c = HiddenGroup (h, 4, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  b = PushButton (c, "Generate Report", GenerateDiscrepancyReport);
  SetObjectExtra (b, drfp, NULL);
  drfp->recheck_btn = PushButton (c, "Recheck", RecheckDiscrepancyProc);
  SetObjectExtra (drfp->recheck_btn, drfp, NULL);
  
  b = PushButton (c, "Configure", EditDiscrepancyConfigBtn);
  SetObjectExtra (b, drfp, NULL);
  
  PushButton (c, "Dismiss", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);

  RealizeWindow (w);
  
  /* adjust column width for discrepancy list */
  ObjectRect (drfp->doc, &r);
  InsetRect (&r, 4, 4);
  discColFmt[1].pixWidth = r.right - r.left - discColFmt[0].pixWidth;

  /* populate discrepancy lists */
  RecheckDiscrepancyProc (drfp->recheck_btn);
  Show (w);
  discrepancyReportWindow = w;
}

