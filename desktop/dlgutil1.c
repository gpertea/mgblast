/*   dlgutil1.c
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
* File Name:  dlgutil1.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   1/22/95
*
* $Revision: 6.97 $
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

#include <dlogutil.h>
#include <document.h>
#include <gather.h>
#include <subutil.h>
#include <objfdef.h>
#include <gbfeat.h>
#include <gbftdef.h>
#include <edutil.h>
#include <explore.h>
#include <sqnutils.h>
#include <alignmgr2.h>
#include <toasn3.h>

#define NUMBER_OF_SUFFIXES    8

static CharPtr name_suffix_labels [] = {
  " ", "Jr.", "Sr.", "II", "III", "IV", "V", "VI", NULL
};

static ENUM_ALIST(name_suffix_alist)
  {" ",    0},
  {"Jr.",  1},
  {"Sr.",  2},
  {"II",   3},
  {"III",  4},
  {"IV",   5},
  {"V",    6},
  {"VI",   7},
END_ENUM_ALIST

Uint2 author_types [] = {
  TAGLIST_TEXT, TAGLIST_TEXT, TAGLIST_TEXT, TAGLIST_POPUP
};

Uint2 std_author_widths [] = {
  8, 4, 9, 0
};

static EnumFieldAssocPtr author_popups [] = {
  NULL, NULL, NULL, name_suffix_alist
};

Uint2 str_author_widths [] = {
  24
};

typedef struct authordialog {
  DIALOG_MESSAGE_BLOCK
  DialoG             stdAuthor;
  DialoG             strAuthor;
  GrouP              stdGrp;
  GrouP              strGrp;
  Uint1              type;
} AuthorDialog, PNTR AuthorDialogPtr;

StdPrintOptionsPtr  spop = NULL;

extern Boolean SetupPrintOptions (void)

{
  if (spop == NULL) {
    spop = StdPrintOptionsNew (NULL);
    if (spop != NULL) {
      spop->newline = "\r";
      spop->indent = "";
    } else {
      Message (MSG_FATAL, "StdPrintOptionsNew failed");
    }
  }
  return (Boolean) (spop != NULL);
}

extern void FreePrintOptions (void)

{
  spop = StdPrintOptionsFree (spop);
}

ENUM_ALIST(months_alist)
  {" ",     0},
  {"Jan",   1},
  {"Feb",   2},
  {"Mar",   3},
  {"Apr",   4},
  {"May",   5},
  {"Jun",   6},
  {"Jul",   7},
  {"Aug",   8},
  {"Sep",   9},
  {"Oct",  10},
  {"Nov",  11},
  {"Dec",  12},
END_ENUM_ALIST

extern void SetDescriptorPropagate (BioseqSetPtr bssp)
{
  BioseqPtr         bsp;
  SeqEntryPtr       seqentry;
  ValNodePtr        sourcedescr;

  if (bssp != NULL) {
    sourcedescr = bssp->descr;
    if (sourcedescr != NULL) {
      bssp->descr = NULL;
      seqentry = bssp->seq_set;
      while (seqentry != NULL) {
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
      SeqDescrFree (sourcedescr);
    }
  }
}

extern Int2 LIBCALLBACK DescriptorPropagate (Pointer data)

{
  BioseqSetPtr      bssp = NULL;
  OMProcControlPtr  ompcp;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->input_entityID == 0) {
    Message (MSG_ERROR, "Please select a BioseqSet");
    return OM_MSG_RET_ERROR;
  }
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQSET :
      bssp = (BioseqSetPtr) ompcp->input_data;
      break;
    case 0 :
      Message (MSG_ERROR, "Please select a BioseqSet");
      return OM_MSG_RET_ERROR;
    default :
      Message (MSG_ERROR, "Please select a BioseqSet");
      return OM_MSG_RET_ERROR;
  }

  SetDescriptorPropagate (bssp);

  ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, 0, 0);
  return OM_MSG_RET_DONE;
}

extern Boolean DescFormReplaceWithoutUpdateProc (ForM f)

{
  MsgAnswer          ans;
  DescriptorFormPtr  dfp;
  Int4Ptr            intptr;
  OMProcControl      ompc;
  Boolean            rsult;
  ValNodePtr         sdp;
  SeqEntryPtr        sep;
  BioseqSetPtr       bssp;

  rsult = FALSE;
  dfp = (DescriptorFormPtr) GetObjectExtra (f);
  if (dfp != NULL) {
    MemSet ((Pointer) &ompc, 0, sizeof (OMProcControl));
    ompc.input_entityID = dfp->input_entityID;
    ompc.input_itemID = dfp->input_itemID;
    ompc.input_itemtype = dfp->input_itemtype;
    ompc.output_itemtype = dfp->input_itemtype;
    sdp = SeqDescrNew (NULL);
    if (sdp != NULL) {
      sdp->choice = dfp->this_subtype;
      switch (sdp->choice) {
        case Seq_descr_mol_type :
        case Seq_descr_method :
          intptr = (Int4Ptr) DialogToPointer (dfp->data);
          if (intptr != NULL) {
            sdp->data.intvalue = *intptr;
          }
          break;
        default :
          sdp->data.ptrvalue = DialogToPointer (dfp->data);
          break;
      }
      ompc.output_data = (Pointer) sdp;
      if (ompc.input_entityID == 0) {
        if (! ObjMgrRegister (OBJ_SEQDESC, (Pointer) sdp)) {
          Message (MSG_ERROR, "ObjMgrRegister failed");
        }
      } else if (ompc.input_itemtype != OBJ_SEQDESC) {
        ompc.output_itemtype = OBJ_SEQDESC;
        if (! AttachDataForProc (&ompc, FALSE)) {
          Message (MSG_ERROR, "AttachDataForProc failed");
        }
        rsult = TRUE;
      } else {
        if (! ReplaceDataForProc (&ompc, FALSE)) {
          Message (MSG_ERROR, "ReplaceDataForProc failed");
        }
        rsult = TRUE;
      }
    }

    /* If the descriptor was added to a GenBank set then*/
    /* optionally propagate it to the set's Bioseqs.    */

    if (ompc.input_itemtype == OBJ_BIOSEQSET) {
      sep = (SeqEntryPtr) ompc.input_choice;
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      if (bssp->_class == BioseqseqSet_class_genbank) {
	ans = Message (MSG_YN, "Do you wish to propagate the descriptor to "
		       "the set's Bioseqs?");
	if (ANS_YES == ans)
	  DescriptorPropagate (&ompc);
      }
    }
    
  }
  return rsult;
}

extern void StdDescFormActnProc (ForM f)

{
  DescriptorFormPtr  dfp;

  if (DescFormReplaceWithoutUpdateProc (f)) {
    dfp = (DescriptorFormPtr) GetObjectExtra (f);
    if (dfp != NULL) {
      GetRidOfEmptyFeatsDescStrings (dfp->input_entityID, NULL);
      if (GetAppProperty ("InternalNcbiSequin") != NULL) {
        ExtendGeneFeatIfOnMRNA (dfp->input_entityID, NULL);
      }
      ObjMgrSetDirtyFlag (dfp->input_entityID, TRUE);
      ObjMgrSendMsg (OM_MSG_UPDATE, dfp->input_entityID,
                     dfp->input_itemID, dfp->input_itemtype);
    }
  }
}

extern void StdDescFormCleanupProc (GraphiC g, VoidPtr data)

{
  DescriptorFormPtr  dfp;
  Uint2              userkey;

  dfp = (DescriptorFormPtr) data;
  if (dfp != NULL) {
    if (dfp->input_entityID > 0 && dfp->userkey > 0) {
      userkey = dfp->userkey;
      dfp->userkey = 0;
      ObjMgrFreeUserData (dfp->input_entityID, dfp->procid, dfp->proctype, userkey);
    }
  }
  StdCleanupExtraProc (g, data);
}

extern OMUserDataPtr ItemAlreadyHasEditor (Uint2 entityID, Uint2 itemID, Uint2 itemtype, Uint2 procid)

{
  BaseFormPtr    bfp;
  Int2           j;
  Int2           num;
  ObjMgrPtr      omp;
  ObjMgrDataPtr  PNTR omdpp;
  OMUserDataPtr  omudp;
  ObjMgrDataPtr  tmp;

  if (entityID == 0 || itemID == 0 || itemtype == 0 || procid == 0) return NULL;
  omp = ObjMgrGet ();
  if (omp == NULL) return NULL;
  num = omp->currobj;
  for (j = 0, omdpp = omp->datalist; j < num && omdpp != NULL; j++, omdpp++) {
    tmp = *omdpp;
    if (tmp->parentptr == NULL && tmp->EntityID == entityID) {

      for (omudp = tmp->userdata; omudp != NULL; omudp = omudp->next) {
        if (omudp->proctype == OMPROC_EDIT && omudp->procid == procid) {
          bfp = (BaseFormPtr) omudp->userdata.ptrvalue;
          if (bfp != NULL) {
            if (bfp->input_itemID == itemID && bfp->input_itemtype == itemtype) {
              return omudp;
            }
          }
        }
      }
    }
  }
  return NULL;
}

typedef struct genegatherlist {
  FeatureFormPtr  ffp;
  ObjMgrPtr       omp;
  SeqLocPtr       slp;
  GeneRefPtr      genexref;
  Boolean         xrefmatch;
  Int2            idx;
  Int2            val;
  Int4            min;
  Uint2           geneEntityID;
  Uint2           geneItemID;
  Uint2           geneItemtype;
  Boolean         geneFound;
  SeqLocPtr       old_feature_location;
  ValNodePtr      lastgene;
} GeneGatherList, PNTR GeneGatherPtr;

static Boolean GeneFindFunc (GatherContextPtr gcp)

{
  GeneGatherPtr  ggp;
  ObjMgrTypePtr  omtp;
  SeqFeatPtr     sfp;
  Char           thislabel [41];

  if (gcp == NULL) return TRUE;

  ggp = (GeneGatherPtr) gcp->userdata;
  if (ggp == NULL ) return TRUE;

  thislabel [0] = '\0';

  if (gcp->thistype == OBJ_SEQFEAT) {
    sfp = (SeqFeatPtr) gcp->thisitem;
    if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE) {
      omtp = ObjMgrTypeFind (ggp->omp, gcp->thistype, NULL, NULL);
      if (omtp == NULL) {
        return TRUE;
      }
      if (omtp->labelfunc != NULL) {
        (*(omtp->labelfunc)) (gcp->thisitem, thislabel, 40, OM_LABEL_CONTENT);
      }
      if (thislabel [0] != '\0') {
        ggp->idx++;
        if (ggp->idx == ggp->val) {
          ggp->geneEntityID = gcp->entityID;
          ggp->geneItemID = gcp->itemID;
          ggp->geneItemtype = gcp->thistype;
          ggp->geneFound = TRUE;
          return FALSE;
        }
      }
    }
  }
  return TRUE;
}

extern void Nlm_LaunchGeneFeatEd (ButtoN b);
extern void Nlm_LaunchGeneFeatEd (ButtoN b)

{
  FeatureFormPtr  ffp;
  GeneGatherList  ggl;
  GatherScope     gs;
  Int2            handled;
  Int2            val;

  ffp = (FeatureFormPtr) GetObjectExtra (b);
  if (ffp != NULL && ffp->gene != NULL && GetValue (ffp->useGeneXref) == 1) {
    val = GetValue (ffp->gene);
    if (val > 2) {
      ggl.ffp = ffp;
      ggl.omp = ObjMgrGet ();
      ggl.idx = 2;
      ggl.val = val;
      ggl.min = INT4_MAX;
      ggl.geneFound = FALSE;
      ggl.geneEntityID = 0;
      ggl.geneItemID = 0;
      ggl.geneItemtype = 0;
      MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
      gs.seglevels = 1;
      gs.get_feats_location = TRUE;
	  MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
	  gs.ignore[OBJ_BIOSEQ] = FALSE;
	  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
	  gs.ignore[OBJ_SEQFEAT] = FALSE;
	  gs.ignore[OBJ_SEQANNOT] = FALSE;
      gs.scope = GetBestTopParentForItemID (ffp->input_entityID,
                                            ffp->input_itemID,
                                            ffp->input_itemtype);
      GatherEntity (ffp->input_entityID, (Pointer) &ggl, GeneFindFunc, &gs);
      if (ggl.geneFound) {
        WatchCursor ();
        Update ();
        handled = GatherProcLaunch (OMPROC_EDIT, FALSE, ggl.geneEntityID, ggl.geneItemID,
                                    ggl.geneItemtype, 0, 0, ggl.geneItemtype, 0);
        ArrowCursor ();
        Update ();
        if (handled != OM_MSG_RET_DONE || handled == OM_MSG_RET_NOPROC) {
          Message (MSG_ERROR, "Unable to launch editor on gene feature.");
        }
      }
    }
  }
}

extern void UpdateGeneLocation 
(SeqFeatPtr gene,
 SeqLocPtr  old_feat_loc,
 SeqLocPtr  new_feat_loc,
 Uint2      entityID)
{
  Uint1          strandfeat, strandgene, strandold;
  BioseqPtr      bsp;
  SeqLocPtr      tmpslp, slp;
  Boolean        hasNulls;
  Boolean        noLeft;
  Boolean        noRight;
  Boolean        noLeftFeat;
  Boolean        noLeftGene;
  Boolean        noRightFeat;
  Boolean        noRightGene;

  if (gene == NULL || new_feat_loc == NULL)
  {
    return;
  }
  
  strandfeat = SeqLocStrand (new_feat_loc);
  strandgene = SeqLocStrand (gene->location); 
  if (old_feat_loc == NULL)
  {
    strandold = strandfeat;
  }
  else
  {
    strandold = SeqLocStrand (old_feat_loc);
  }
          
  /* only correct gene location if gene is on same strand as old feature
   * location and contained in new feature location (on either strand).
   */
  if (SeqLocAinB (new_feat_loc, gene->location) <= 0
      && ((strandold == Seq_strand_minus && strandgene == Seq_strand_minus)
          || (strandold != Seq_strand_minus && strandgene != Seq_strand_minus))) 
  {
    bsp = GetBioseqGivenSeqLoc (gene->location, entityID);
    if (bsp != NULL) {
      hasNulls = LocationHasNullsBetween (gene->location);
      if ((strandfeat == Seq_strand_minus && strandgene != Seq_strand_minus)
          || (strandfeat != Seq_strand_minus && strandgene == Seq_strand_minus))
      {
        tmpslp = SeqLocCopy (gene->location);
        SeqLocRevCmp (tmpslp);
        slp = SeqLocMerge (bsp, tmpslp, new_feat_loc, TRUE, FALSE, hasNulls);
        tmpslp = SeqLocFree (tmpslp);
      }
      else
      {
        slp = SeqLocMerge (bsp, gene->location, new_feat_loc, TRUE, FALSE, hasNulls);
      }

      if (slp != NULL) {
        CheckSeqLocForPartial (gene->location, &noLeftGene, &noRightGene);
        gene->location = SeqLocFree (gene->location);
        gene->location = slp;
        CheckSeqLocForPartial (new_feat_loc, &noLeftFeat, &noRightFeat);
        if (bsp->repr == Seq_repr_seg) {
          slp = SegLocToPartsEx (bsp, gene->location, TRUE);
          gene->location = SeqLocFree (gene->location);
          gene->location = slp;
          hasNulls = LocationHasNullsBetween (gene->location);
          gene->partial = (gene->partial || hasNulls);
        }
        FreeAllFuzz (gene->location);
        noLeft = (noLeftFeat || noLeftGene);
        noRight = (noRightFeat || noRightGene);
        SetSeqLocPartial (gene->location, noLeft, noRight);
        gene->partial = (gene->partial || noLeft || noRight);
      }
    }
  }
}

static Boolean GeneUpdateFunc (GatherContextPtr gcp)

{
  GeneGatherPtr  ggp;
  ObjMgrTypePtr  omtp;
  SeqFeatPtr     sfp;
  Char           thislabel [41];

  if (gcp == NULL) return TRUE;

  ggp = (GeneGatherPtr) gcp->userdata;
  if (ggp == NULL ) return TRUE;

  thislabel [0] = '\0';

  if (gcp->thistype == OBJ_SEQFEAT) {
    sfp = (SeqFeatPtr) gcp->thisitem;
    if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE) {
      omtp = ObjMgrTypeFind (ggp->omp, gcp->thistype, NULL, NULL);
      if (omtp == NULL) {
        return TRUE;
      }
      if (omtp->labelfunc != NULL) {
        (*(omtp->labelfunc)) (gcp->thisitem, thislabel, 40, OM_LABEL_CONTENT);
      }
      if (thislabel [0] != '\0') {
        ggp->idx++;
        if (ggp->idx == ggp->val) {
          UpdateGeneLocation (sfp, ggp->old_feature_location, ggp->slp, gcp->entityID);
          return FALSE;
        }
      }
    }
  }
  return TRUE;
}

static Boolean GeneGatherFunc (GatherContextPtr gcp)

{
  FeatureFormPtr  ffp;
  GeneGatherPtr   ggp;
  GeneRefPtr      grp;
  ObjMgrTypePtr   omtp;
  SeqFeatPtr      sfp;
  Char            thislabel [41];
  ValNodePtr      vnp;

  if (gcp == NULL) return TRUE;

  ggp = (GeneGatherPtr) gcp->userdata;
  if (ggp == NULL || ggp->ffp == NULL) return TRUE;

  thislabel [0] = '\0';

  if (gcp->thistype == OBJ_SEQFEAT) {
    sfp = (SeqFeatPtr) gcp->thisitem;
    if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE) {
      omtp = ObjMgrTypeFind (ggp->omp, gcp->thistype, NULL, NULL);
      if (omtp == NULL) {
        return TRUE;
      }
      if (omtp->labelfunc != NULL) {
        (*(omtp->labelfunc)) (gcp->thisitem, thislabel, 40, OM_LABEL_CONTENT);
      }
      if (thislabel [0] != '\0') {
        ffp = (FeatureFormPtr) ggp->ffp;
        grp = (GeneRefPtr) sfp->data.value.ptrvalue;
        if (grp != NULL 
            && (grp->locus != NULL || grp->locus_tag != NULL || grp->desc != NULL))
        {
          vnp = ValNodeNew (ggp->lastgene);
          if (ffp->geneNames == NULL) {
            ffp->geneNames = vnp;
          }
          ggp->lastgene = vnp;
          if (vnp != NULL) {
            vnp->data.ptrvalue = StringSave (thislabel);
            if (grp->locus != NULL) {
              vnp->choice = 1;
            } else if (grp->desc != NULL) {
              vnp->choice = 2;
            } else if (grp->locus_tag != NULL) {
              vnp->choice = 3;
            }
          }
        }
      }
    }
  }

  return TRUE;
}

extern void PopulateGenePopup (FeatureFormPtr ffp)

{
  Int2            count;
  GeneGatherList  ggl;
  GatherScope     gs;
  CharPtr         str;
  Boolean         usePopupForGene;
  ValNodePtr      vnp;

  if (ffp != NULL && ffp->genePopup != NULL && ffp->geneList != NULL) {
    ggl.ffp = ffp;
    ggl.omp = ObjMgrGet ();
    ggl.slp = NULL;
    ggl.genexref = NULL;
    ggl.xrefmatch = FALSE;
    ggl.idx = 0;
    ggl.val = 0;
    ggl.min = 0;
    ggl.lastgene = NULL;
    MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
    gs.seglevels = 1;
    gs.get_feats_location = TRUE;
	MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
	gs.ignore[OBJ_BIOSEQ] = FALSE;
	gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
	gs.ignore[OBJ_SEQFEAT] = FALSE;
	gs.ignore[OBJ_SEQANNOT] = FALSE;
    gs.scope = GetBestTopParentForItemID (ffp->input_entityID,
                                          ffp->input_itemID,
                                          ffp->input_itemtype);
    GatherEntity (ffp->input_entityID, (Pointer) &ggl, GeneGatherFunc, &gs);
    count = 0;
    for (vnp = ffp->geneNames; vnp != NULL; vnp = vnp->next) {
      count++;
    }
    if (count < 32) {
      usePopupForGene = TRUE;
      ffp->gene = ffp->genePopup;
    } else {
      usePopupForGene = FALSE;
      ffp->gene = ffp->geneList;
    }
    for (vnp = ffp->geneNames; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) {
        str = "??";
      }
      if (usePopupForGene) {
        PopupItem (ffp->gene, str);
      } else {
        ListItem (ffp->gene, str);
      }
    }
    Show (ffp->gene);
  }
}

static Boolean GeneMatchFunc (GatherContextPtr gcp)

{
  Int4            diff;
  FeatureFormPtr  ffp;
  GeneRefPtr      genexref;
  GeneGatherPtr   ggp;
  GeneRefPtr      grp;
  ObjMgrTypePtr   omtp;
  SeqFeatPtr      sfp;
  Char            thislabel [41];

  if (gcp == NULL) return TRUE;

  ggp = (GeneGatherPtr) gcp->userdata;
  if (ggp == NULL || ggp->ffp == NULL) return TRUE;

  thislabel [0] = '\0';

  if (gcp->thistype == OBJ_SEQFEAT) {
    sfp = (SeqFeatPtr) gcp->thisitem;
    if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE && sfp->data.value.ptrvalue != NULL) {
      omtp = ObjMgrTypeFind (ggp->omp, gcp->thistype, NULL, NULL);
      if (omtp == NULL) {
        return TRUE;
      }
      if (omtp->labelfunc != NULL) {
        (*(omtp->labelfunc)) (gcp->thisitem, thislabel, 40, OM_LABEL_CONTENT);
      }
      if (thislabel [0] != '\0') {
        ffp = (FeatureFormPtr) ggp->ffp;
        ggp->idx++;
        genexref = ggp->genexref;
        if (genexref != NULL) {
          grp = (GeneRefPtr) sfp->data.value.ptrvalue;
          if (! StringHasNoText (genexref->locus)) {
            if (StringICmp (genexref->locus, grp->locus) == 0) {
              ggp->val = ggp->idx;
              ggp->xrefmatch = TRUE;
              if (ffp != NULL) {
                SetValue (ffp->useGeneXref, 2);
              }
            }
          } else if (! StringHasNoText (genexref->locus_tag)) {
            if (StringICmp (genexref->locus_tag, grp->locus_tag) == 0) {
              ggp->val = ggp->idx;
              ggp->xrefmatch = TRUE;
              if (ffp != NULL) {
                SetValue (ffp->useGeneXref, 2);
              }
            }
          } else if (! StringHasNoText (genexref->desc)) {
            if (StringICmp (genexref->desc, grp->desc) == 0) {
              ggp->val = ggp->idx;
              ggp->xrefmatch = TRUE;
              if (ffp != NULL) {
                SetValue (ffp->useGeneXref, 2);
              }
            }
          }
        }
        diff = SeqLocAinB (ggp->slp, sfp->location);
        if (diff >= 0) {
          if (diff < ggp->min) {
            ggp->min = diff;
            if (! ggp->xrefmatch) {
              ggp->val = ggp->idx;
            }
          }
        }
      }
    }
  }

  return TRUE;
}

static void SaveGoTermsInSfp (UserObjectPtr uop, Pointer userdata)

{
  FeatureFormPtr  ffp;
  ObjectIdPtr     oip;

  if (uop == NULL || userdata == NULL) return;
  oip = uop->type;
  if (oip == NULL) return;
  if (StringCmp (oip->str, "GeneOntology") == 0) {
    ffp = (FeatureFormPtr) userdata;
    ffp->goTermUserObj = AsnIoMemCopy (uop, (AsnReadFunc) UserObjectAsnRead,
                                       (AsnWriteFunc) UserObjectAsnWrite);
  }
}

static void FeatIDtoText (TexT t, ChoicePtr cp)

{
  Char         buf [32];
  ObjectIdPtr  oip;

  if (t == NULL) return;
  if (cp == NULL) {
    SetTitle (t, "");
    return;
  }

  if (cp->choice == 3) {
    oip = (ObjectIdPtr) cp->value.ptrvalue;
    if (oip != NULL) {
      if (StringDoesHaveText (oip->str)) {
        SetTitle (t, oip->str);
        return;
      } else {
        sprintf (buf, "%ld", (long) oip->id);
        SetTitle (t, buf);
        return;
      }
    }
  }

  SetTitle (t, "");
}

static void TextToFeatID (TexT t, ChoicePtr cp)

{
  Boolean      all_digits = TRUE;
  Char         buf [128];
  Char         ch;
  ObjectIdPtr  oip;
  CharPtr      str;
  long int     val;

  if (t == NULL || cp == NULL) return;

  GetTitle (t, buf, sizeof (buf) - 1);
  if (StringHasNoText (buf)) {
    SeqFeatIdFree (cp);
    cp->choice = 0;
    return;
  }

  oip = ObjectIdNew ();
  if (oip == NULL) return;

  str = buf;
  ch = *str;
  while (ch != '\0') {
    if (! IS_DIGIT (ch)) {
      all_digits = FALSE;
    }
    str++;
    ch = *str;
  }

  if (all_digits && sscanf (buf, "%ld", &val) == 1) {
    oip->id = (Int4) val;
  } else {
    oip->str = StringSave (buf);
  }
  SeqFeatIdFree (cp);
  cp->choice = 3;
  cp->value.ptrvalue = (Pointer) oip;
}

static void FeatXreftoText (TexT t, SeqFeatPtr sfp)

{
  Char            buf [32];
  ObjectIdPtr     oip;
  SeqFeatXrefPtr  xref;

  if (t == NULL) return;
  if (sfp == NULL) {
    SetTitle (t, "");
    return;
  }

  for (xref = sfp->xref; xref != NULL; xref = xref->next) {
    if (xref->id.choice != 3) continue;
    oip = (ObjectIdPtr) xref->id.value.ptrvalue;
    if (oip != NULL) {
      if (StringDoesHaveText (oip->str)) {
        SetTitle (t, oip->str);
        return;
      } else {
        sprintf (buf, "%ld", (long) oip->id);
        SetTitle (t, buf);
        return;
      }
    }
  }

  SetTitle (t, "");
}

static void TextToFeatXref (TexT t, SeqFeatPtr sfp)

{
  Boolean         all_digits = TRUE;
  Char            buf [128];
  Char            ch;
  ObjectIdPtr     oip;
  CharPtr         str;
  long int        val;
  SeqFeatXrefPtr  xref;

  if (t == NULL || sfp == NULL) return;

  GetTitle (t, buf, sizeof (buf) - 1);
  if (StringHasNoText (buf)) {
    ClearFeatIDXrefs (sfp);
    return;
  }

  ClearFeatIDXrefs (sfp);

  oip = ObjectIdNew ();
  if (oip == NULL) return;

  str = buf;
  ch = *str;
  while (ch != '\0') {
    if (! IS_DIGIT (ch)) {
      all_digits = FALSE;
    }
    str++;
    ch = *str;
  }

  if (all_digits && sscanf (buf, "%ld", &val) == 1) {
    oip->id = (Int4) val;
  } else {
    oip->str = StringSave (buf);
  }

  xref = SeqFeatXrefNew ();
  if (xref != NULL) {
    xref->id.choice = 3;
    xref->id.value.ptrvalue = (Pointer) oip;
    xref->next = sfp->xref;
    sfp->xref = xref;
  }
}

static void GbqualsToVisStringDialog (SeqFeatPtr sfp, DialoG d, CharPtr qual)

{
  GBQualPtr   gbq;
  ValNodePtr  head = NULL;

  if (sfp == NULL || StringHasNoText (qual)) {
    PointerToDialog (d, NULL);
    return;
  }
  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (StringICmp (gbq->qual, qual) != 0) continue;
    ValNodeCopyStr (&head, 0, gbq->val);
  }
  PointerToDialog (d, head);
  ValNodeFreeData (head);
}

extern void VisStringDialogToGbquals (SeqFeatPtr sfp, DialoG d, CharPtr qual);
extern void VisStringDialogToGbquals (SeqFeatPtr sfp, DialoG d, CharPtr qual)

{
  GBQualPtr   gbq, gbqlast = NULL;
  ValNodePtr  head = NULL, vnp;
  CharPtr     str;

  if (sfp == NULL || StringHasNoText (qual)) return;
  head = DialogToPointer (d);
  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    gbqlast = gbq;
  }
  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    gbq = GBQualNew ();
    if (gbq == NULL) continue;
    gbq->qual = StringSave (qual);
    gbq->val = StringSave (str);
    if (gbqlast == NULL) {
      sfp->qual = gbq;
    } else {
      gbqlast->next = gbq;
    }
    gbqlast = gbq;
  }
  ValNodeFreeData (head);
}

extern void SeqFeatPtrToFieldPage (DialoG d, SeqFeatPtr sfp);

extern void SeqFeatPtrToCommon (FeatureFormPtr ffp, SeqFeatPtr sfp)

{
  GeneGatherList  ggl;
  GeneRefPtr      grp;
  GatherScope     gs;
  ProtRefPtr      prp;
  CharPtr         str;
  ValNodePtr      vnp;
  SeqFeatXrefPtr  xref;
  /*
  Char            ch;
  CharPtr         ptr;
  */

  if (ffp != NULL) {
    if (sfp != NULL) {
      str = StringSave (sfp->comment);
      /*
      ptr = str;
      if (ptr != NULL) {
        ch = *ptr;
        while (ch != '\0') {
          if (ch == '~') {
#ifdef WIN_MAC
            *ptr = '\015';
#else
            *ptr = '\n';
#endif
          }
          ptr++;
          ch = *ptr;
        }
      }
      */
      SetTitle (ffp->comment, str);
      SetTitle (ffp->title, sfp->title);
      SetValue (ffp->evidence, sfp->exp_ev + 1);
      if (GetAppProperty ("InternalNcbiSequin") == NULL) {
        if (sfp->exp_ev == 0) {
          SafeDisable (ffp->evidence);
        }
      }
      SetStatus (ffp->partial, sfp->partial);
      SetStatus (ffp->exception, sfp->excpt);
      SetStatus (ffp->pseudo, sfp->pseudo);
      SetTitle (ffp->exceptText, sfp->except_text);
      SetValue (ffp->useGeneXref, 1);
      SetTitle (ffp->geneSymbol, "");
      SetTitle (ffp->geneAllele, "");
      SetTitle (ffp->geneDesc, "");
      SetTitle (ffp->locusTag, "");
      ggl.ffp = ffp;
      ggl.omp = ObjMgrGet ();
      ggl.slp = sfp->location;
      ggl.genexref = NULL;
      grp = NULL;
      xref = sfp->xref;
      while (xref != NULL && xref->data.choice != SEQFEAT_PROT) {
        xref = xref->next;
      }
      if (xref != NULL) {
        prp = (ProtRefPtr) xref->data.value.ptrvalue;
        if (prp != NULL && ffp->protXrefName != NULL) {
          vnp = prp->name;
          if (vnp != NULL) {
            SetTitle (ffp->protXrefName, (CharPtr) vnp->data.ptrvalue);
          }
        }
      }
      xref = sfp->xref;
      while (xref != NULL && xref->data.choice != SEQFEAT_GENE) {
        xref = xref->next;
      }
      if (xref != NULL) {
        grp = (GeneRefPtr) xref->data.value.ptrvalue;
        ggl.genexref = grp;
      }
      ggl.xrefmatch = FALSE;
      ggl.idx = 2;
      ggl.val = 1;
      ggl.min = INT4_MAX;
      MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
      gs.seglevels = 1;
      gs.get_feats_location = TRUE;
	  MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
	  gs.ignore[OBJ_BIOSEQ] = FALSE;
	  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
	  gs.ignore[OBJ_SEQFEAT] = FALSE;
	  gs.ignore[OBJ_SEQANNOT] = FALSE;
      gs.scope = GetBestTopParentForItemID (ffp->input_entityID,
                                            ffp->input_itemID,
                                            ffp->input_itemtype);
      GatherEntity (ffp->input_entityID, (Pointer) &ggl, GeneMatchFunc, &gs);
      if (grp != NULL && StringHasNoText (grp->locus) && StringHasNoText (grp->allele) &&
            StringHasNoText (grp->desc) && StringHasNoText (grp->maploc) &&
            grp->db == NULL && grp->syn == NULL && grp->locus_tag == NULL) {
        SetValue (ffp->gene, 1);
        SetValue (ffp->useGeneXref, 3);
        SetTitle (ffp->geneSymbol, grp->locus);
        SetTitle (ffp->geneAllele, grp->allele);
        SetTitle (ffp->geneDesc, grp->desc);
        SetTitle (ffp->locusTag, grp->locus_tag);
        SafeHide (ffp->editGeneBtn);
        SafeHide (ffp->newGeneGrp);
      } else if (ggl.val == 1 && grp != NULL && (! ggl.xrefmatch)) {
        SetValue (ffp->gene, 2);
        SetValue (ffp->useGeneXref, 2);
        SetTitle (ffp->geneSymbol, grp->locus);
        SetTitle (ffp->geneAllele, grp->allele);
        SetTitle (ffp->geneDesc, grp->desc);
        SetTitle (ffp->locusTag, grp->locus_tag);
        SafeHide (ffp->editGeneBtn);
        SafeShow (ffp->newGeneGrp);
      } else if (grp != NULL && (! ggl.xrefmatch)) {
        SetValue (ffp->gene, 2);
        SetValue (ffp->useGeneXref, 2);
        SetTitle (ffp->geneSymbol, grp->locus);
        SetTitle (ffp->geneAllele, grp->allele);
        SetTitle (ffp->geneDesc, grp->desc);
        SetTitle (ffp->locusTag, grp->locus_tag);
        SafeHide (ffp->editGeneBtn);
        SafeShow (ffp->newGeneGrp);
      } else if (ggl.val > 2) {
        SetValue (ffp->gene, ggl.val);
        SafeShow (ffp->editGeneBtn);
        if (ffp->gene == ffp->geneList) {
          SetOffset (ffp->gene, 0, (Int2) (ggl.val - 1));
        }
        if (grp != NULL && (! ggl.xrefmatch)) {
          SetValue (ffp->useGeneXref, 2);
        }
      } else {
        SetValue (ffp->gene, ggl.val);
      }
      PointerToDialog (ffp->featcits, sfp->cit);
      PointerToDialog (ffp->dbxrefs, sfp->dbxref);
      PointerToDialog (ffp->gbquals, sfp->qual);
      SeqFeatPtrToFieldPage (ffp->gbquals, sfp);
      GbqualsToVisStringDialog (sfp, ffp->experiment, "experiment");
      GBQualsToInferenceDialog (ffp->inference, sfp);
      PointerToDialog (ffp->usrobjext, sfp->ext);
      VisitUserObjectsInUop (sfp->ext, (Pointer) ffp, SaveGoTermsInSfp);
      FeatIDtoText (ffp->featid, &(sfp->id));
      FeatXreftoText (ffp->fidxref, sfp);
    } else {
      SetTitle (ffp->comment, "");
      SetValue (ffp->evidence, 1);
      if (GetAppProperty ("InternalNcbiSequin") == NULL) {
        SafeDisable (ffp->evidence);
      }
      SetStatus (ffp->partial, FALSE);
      SetStatus (ffp->exception, FALSE);
      SetStatus (ffp->pseudo, FALSE);
      SetTitle (ffp->exceptText, "");
      SetValue (ffp->gene, 1);
      SafeHide (ffp->newGeneGrp);
      SafeHide (ffp->editGeneBtn);      
      SetValue (ffp->useGeneXref, 1);
      SetTitle (ffp->geneSymbol, "");
      SetTitle (ffp->geneAllele, "");
      SetTitle (ffp->geneDesc, "");
      SetTitle (ffp->locusTag, "");
      PointerToDialog (ffp->featcits, NULL);
      PointerToDialog (ffp->dbxrefs, NULL);
      PointerToDialog (ffp->gbquals, NULL);
      PointerToDialog (ffp->experiment, NULL);
      GBQualsToInferenceDialog (ffp->inference, NULL);
      PointerToDialog (ffp->usrobjext, NULL);
      ffp->goTermUserObj = NULL;
      FeatIDtoText (ffp->featid, NULL);
      FeatXreftoText (ffp->fidxref, NULL);
    }
  }
}

extern void CleanupEvidenceGBQuals (GBQualPtr PNTR prevgbq);
extern void CleanupEvidenceGBQuals (GBQualPtr PNTR prevgbq)

{
  GBQualPtr  gbq;
  GBQualPtr  next;
  Boolean    unlink;

  if (prevgbq == NULL) return;
  gbq = *prevgbq;
  while (gbq != NULL) {
    next = gbq->next;
    unlink = FALSE;
    if (StringICmp (gbq->qual, "experiment") == 0 || StringICmp (gbq->qual, "inference") == 0) {
      unlink = TRUE;
    }
    if (unlink) {
      *prevgbq = gbq->next;
      gbq->next = NULL;
      GBQualFree (gbq);
    } else {
      prevgbq = (GBQualPtr PNTR) &(gbq->next);
    }
    gbq = next;
  }
}

typedef struct replacesdata {
  FeatureFormPtr  ffp;
  SeqFeatPtr      sfp;
} ReplaceData, PNTR ReplaceDataPtr;

static Boolean ReplaceFeatureExtras (GatherContextPtr gcp)

{
  SeqFeatXrefPtr  curr;
  FeatureFormPtr  ffp;
  SeqFeatXrefPtr  PNTR last;
  SeqFeatXrefPtr  next;
  SeqFeatPtr      old;
  ReplaceDataPtr  rdp;
  SeqFeatPtr      sfp;

  rdp = (ReplaceDataPtr) gcp->userdata;
  if (rdp != NULL && rdp->sfp != NULL && rdp->ffp != NULL) {
    sfp = rdp->sfp;
    ffp = rdp->ffp;
    old = gcp->thisitem;
    if (old != NULL) {
      if (ffp->gbquals != NULL) {
        sfp->qual = DialogToPointer (ffp->gbquals);
      } else if (sfp->qual == NULL) {
        sfp->qual = old->qual;
        old->qual = NULL;
      }
      CleanupEvidenceGBQuals (&(sfp->qual));
      VisStringDialogToGbquals (sfp, ffp->experiment, "experiment");
      InferenceDialogToGBQuals (ffp->inference, sfp, TRUE);
      if (ffp->usrobjext != NULL) {
        sfp->ext = DialogToPointer (ffp->usrobjext);
      } else if (sfp->ext == NULL) {
        sfp->ext = old->ext;
        old->ext = NULL;
      }
      if (ffp->goTermUserObj != NULL) {
        sfp->ext = CombineUserObjects (sfp->ext, ffp->goTermUserObj);
      }
      if (ffp->featid != NULL) {
        TextToFeatID (ffp->featid, &(sfp->id));
      }
      /*
      if (sfp->cit == NULL) {
        sfp->cit = old->cit;
        old->cit = NULL;
      }
      */
      if (old->xref != NULL) {
        last = (SeqFeatXrefPtr PNTR) &(old->xref);
        curr = old->xref;
        while (curr != NULL) {
          next = curr->next;
          if (curr->data.choice == SEQFEAT_GENE) {
            *last = next;
            curr->next = NULL;
            SeqFeatXrefFree (curr);
          } else {
            last = &(curr->next);
          }
          curr = next;
        }
      }
      if (sfp->xref == NULL) {
        sfp->xref = old->xref;
      } else {
        curr = sfp->xref;
        while (curr->next != NULL) {
          curr = curr->next;
        }
        if (curr != NULL) {
          curr->next = old->xref;
        }
      }
      old->xref = NULL;
    }
  }
  return TRUE;
}

static Boolean GetOldFeatureLocation (GatherContextPtr gcp)
{
  SeqLocPtr PNTR  old_loc;
  SeqFeatPtr      old_feat;

  if (gcp == NULL || gcp->userdata == NULL)
  {
    return FALSE;
  }
  old_loc = (SeqLocPtr PNTR) gcp->userdata;
  if (*old_loc != NULL)
  {
    return TRUE;
  }
  old_feat = gcp->thisitem;
  if (old_feat != NULL)
  {
    *old_loc = (SeqLocPtr) AsnIoMemCopy (old_feat->location,
                                         (AsnReadFunc) SeqLocAsnRead,
                                         (AsnWriteFunc) SeqLocAsnWrite);
  }
  return TRUE;
}

static Boolean HasExceptionGBQual (SeqFeatPtr sfp)

{
  GBQualPtr  qual;
  Boolean    rsult;

  rsult = FALSE;
  if (sfp != NULL) {
    qual = sfp->qual;
    while (qual != NULL) {
      if (StringICmp (qual->qual, "exception") == 0) {
        if (! StringHasNoText (qual->val)) {
          rsult = TRUE;
        }
      }
      qual = qual->next;
    }
  }
  return rsult;
}

static void AddProtRefXref (SeqFeatPtr sfp, TexT protXrefName)

{
  ProtRefPtr      prp;
  Char            str [256];
  SeqFeatXrefPtr  xref;

  if (sfp == NULL || protXrefName == NULL) return;
  GetTitle (protXrefName, str, sizeof (str) - 1);
  if (StringHasNoText (str)) return;
  for (xref = sfp->xref; xref != NULL; xref = xref->next) {
    if (xref->data.choice == SEQFEAT_PROT) break;
  }
  if (xref == NULL) {
    xref = SeqFeatXrefNew ();
    if (xref != NULL) {
      prp = ProtRefNew ();
      xref->data.choice = SEQFEAT_PROT;
      xref->data.value.ptrvalue = (Pointer) prp;
      xref->next = sfp->xref;
      sfp->xref = xref;
    }
  }
  if (xref != NULL && xref->data.choice == SEQFEAT_PROT) {
    prp = (ProtRefPtr) xref->data.value.ptrvalue;
    xref->data.value.ptrvalue = ProtRefFree (prp);
    prp = CreateNewProtRef (str, NULL, NULL, NULL);
    xref->data.value.ptrvalue = (Pointer) prp;
  }
}

static Boolean IsPseudo (SeqFeatPtr sfp)
{
  Boolean   pseudo = FALSE;
  GBQualPtr qual;
  
  if (sfp != NULL)
  {
    if (sfp->pseudo)
    {
      pseudo = TRUE;
    }
    else
    {
      qual = sfp->qual;
      while (qual != NULL)
      {
        if (StringICmp (qual->qual, "pseudo") == 0) 
        {
          pseudo = TRUE;
        }
        qual = qual->next;
      }
    }
  }
  return pseudo;
}

static void RemoveProtXrefs (SeqFeatPtr sfp)
{
  SeqFeatXrefPtr  xref_prev = NULL, xref_next, xref;

  if (sfp == NULL || sfp->xref == NULL) return;

  for (xref = sfp->xref; xref != NULL; xref = xref_next)
  {
    xref_next = xref->next;
    if (xref->data.choice == SEQFEAT_PROT)
    {
      if (xref_prev == NULL)
      {
        sfp->xref = xref->next;
      }
      else
      {
        xref_prev->next = xref->next;
      }
      xref->next = NULL;
      xref->data.value.ptrvalue = ProtRefFree (xref->data.value.ptrvalue);
      SeqFeatXrefFree (xref);
    }
    else
    {
      xref_prev = xref;
    }
  }  
}

static CharPtr infDetails [] = {
  "unknown error",
  "empty inference string",
  "bad inference prefix",
  "bad inference body",
  "single inference field",
  "spaces in inference",
  "same species misused",
  "bad inference accession",
  "accession is missing version",
  "accession.version not public",
  NULL
};

extern Boolean TestInference (FeatureFormPtr ffp, CharPtr badInfQual, size_t len, CharPtr badInfMssg)

{
  GBQualPtr   gbq;
  Int2        inferenceCode;
  Boolean     rsult = TRUE;
  CharPtr     str;
  SeqFeatPtr  tmp;

  if (ffp == NULL) return FALSE;
  tmp = SeqFeatNew ();
  if (tmp == NULL) return FALSE;
  if (badInfQual != NULL) {
    *badInfQual = '\0';
  }
  if (badInfMssg != NULL) {
    *badInfMssg = '\0';
  }
  InferenceDialogToGBQuals (ffp->inference, tmp, FALSE);
  for (gbq = tmp->qual; gbq != NULL; gbq = gbq->next) {
    str = gbq->val;
    if (StringICmp (gbq->qual, "inference") != 0) continue;
    inferenceCode = ValidateInferenceQualifier (str, FALSE);
    if (inferenceCode != VALID_INFERENCE) {
  	  if (inferenceCode < VALID_INFERENCE || inferenceCode > ACC_VERSION_NOT_PUBLIC) {
  	    inferenceCode = VALID_INFERENCE;
  	  }
  	  if (badInfQual != NULL) {
        StringNCpy_0 (badInfQual, str, len);
  	  }
  	  if (badInfMssg != NULL) {
        StringCpy (badInfMssg, infDetails [(int) inferenceCode]);
  	  }
      rsult = FALSE;
      break;
    }
  }
  SeqFeatFree (tmp);
  return rsult;
}

static CharPtr infWarning1 =
"Bad inference qualifier!  You must conform to international nucleotide sequence database\nconventions for inference qualifiers!";

static CharPtr infWarning2 =
"Still has bad inference!  You must conform to international nucleotide sequence database\nconventions for inference qualifiers!";

static CharPtr infAccept =
"Do you want to accept changes with bad inference data?  The bad qualifier will be converted to a note!";

extern Boolean FeatFormReplaceWithoutUpdateProc (ForM f)

{
  Char            allele [128];
  MsgAnswer       ans;
  Int2            attempts = 3;
  Char            badInfMssg [32];
  Char            badInfQual [256];
  BioseqPtr       bsp;
  Char            ch;
  Char            desc [128];
  Int2            expev;
  SeqMgrFeatContext  fcontext;
  FeatureFormPtr  ffp;
  SeqFeatPtr      gene;
  GeneGatherList  ggl;
  GeneRefPtr      grp;
  GeneRefPtr      grpfeat;
  GatherScope     gs;
  SeqLocPtr       gslp;
  Boolean         hasNulls;
  Int2            i;
  Int4Ptr         intptr;
  Uint2           itemID;
  Char            locustag [128];
  Boolean         noLeft;
  Boolean         noRight;
  SeqEntryPtr     oldscope;
  OMProcControl   ompc;
  CharPtr         ptr;
  ReplaceData     rd;
  Boolean         rsult;
  SeqAnnotPtr     sap;
  SeqEntryPtr     sep;
  SeqFeatPtr      sfp;
  SeqLocPtr       slp;
  CharPtr         str;
  Char            symbol [128];
  Int2            usexref;
  Int2            val;
  ValNodePtr      vnp, err_list;
  SeqFeatXrefPtr  xref;
  SeqLocPtr       old_location = NULL; /* we need the old location of the feature 
                                        * if we're going to do a gene update    
                                        */
  Boolean         fix_interval_order;

  rsult = FALSE;
  ffp = (FeatureFormPtr) GetObjectExtra (f);
  if (ffp != NULL) {
    MemSet ((Pointer) &ompc, 0, sizeof (OMProcControl));
    ompc.input_entityID = ffp->input_entityID;
    ompc.input_itemID = ffp->input_itemID;
    ompc.input_itemtype = ffp->input_itemtype;
    ompc.output_itemtype = ffp->input_itemtype;
    sfp = SeqFeatNew ();
    if (sfp != NULL) {
      sep = GetTopSeqEntryForEntityID (ffp->input_entityID);
      oldscope = SeqEntrySetScope (sep);
      sfp->data.choice = FindFeatFromFeatDefType (ffp->this_subtype);
      switch (sfp->data.choice) {
        case SEQFEAT_BOND :
        case SEQFEAT_SITE :
        case SEQFEAT_PSEC_STR :
          intptr = (Int4Ptr) DialogToPointer (ffp->data);
          if (intptr != NULL) {
            sfp->data.value.intvalue = *intptr;
          }
          break;
        case SEQFEAT_COMMENT:
          sfp->data.value.ptrvalue = NULL;
          break;
        default :
          sfp->data.value.ptrvalue = DialogToPointer (ffp->data);
          break;
      }
      sfp->comment = SaveStringFromText (ffp->comment);
      ptr = sfp->comment;
      if (ptr != NULL) {
        ch = *ptr;
        while (ch != '\0') {
          if (ch < ' ' || ch > '~') {
            *ptr = '~';
          }
          ptr++;
          ch = *ptr;
        }
      }
      expev = GetValue (ffp->evidence);
      if (expev > 0 && expev <= 3) {
        sfp->exp_ev = expev - 1;
      } else {
        sfp->exp_ev = 0;
      }
      sfp->partial = GetStatus (ffp->partial);
      sfp->excpt = GetStatus (ffp->exception);
      sfp->pseudo = GetStatus (ffp->pseudo);
      sfp->except_text = SaveStringFromText (ffp->exceptText);
      sfp->title = NULL;
      sfp->product = DialogToPointer (ffp->product);
      sfp->location = DialogToPointer (ffp->location);
      if (sfp->location == NULL) {
        SeqEntrySetScope (oldscope);
        ErrPostEx (SEV_ERROR, 0, 0, "Feature must have a location!");
        err_list = TestDialog (ffp->location);
        DisplayErrorMessages ("Location Errors", err_list);
        err_list = ValNodeFree (err_list);    
        return FALSE;
      }
      if ((! ffp->acceptBadInf) && (! TestInference (ffp, badInfQual, sizeof (badInfQual), badInfMssg))) {
        (ffp->badInfAttempts)++;
        if (GetAppProperty ("InternalNcbiSequin") != NULL) {
          attempts = 2;
        }
        if (ffp->badInfAttempts < attempts) {
          if (ffp->badInfAttempts == 2) {
            Message (MSG_OK, "%s - Please fix the %s error in\n%s", infWarning2, badInfMssg, badInfQual);
          } else {
            Message (MSG_OK, "%s - Please fix the %s error in\n%s", infWarning1, badInfMssg, badInfQual);
          }
          return FALSE;
        } else {
          if (Message (MSG_YN, "%s", infAccept) == ANS_NO) return;
          ffp->acceptBadInf = TRUE;
        }
      }
      bsp = GetBioseqGivenSeqLoc (sfp->location, ffp->input_entityID);
      if (bsp != NULL) {
        if (SeqLocBadSortOrder (bsp, sfp->location)) {
          fix_interval_order = FALSE;
          if (StringICmp (sfp->except_text, "trans-splicing") == 0)
          {
            ans = Message (MSG_YN, "Your feature intervals are out of order, but this coding region has a trans-splicing exception. Continue without correcting interval order?");
            if (ans == ANS_NO)
            {
              fix_interval_order = TRUE;
            }
          }
          else
          {
            ans = Message (MSG_YN,
            "Feature location intervals are out of order.  Do you want them repaired?");
            if (ans == ANS_YES)
            {
              fix_interval_order = TRUE;
            }
          }
          if (fix_interval_order) {
            hasNulls = LocationHasNullsBetween (sfp->location);
            gslp = SeqLocMerge (bsp, sfp->location, NULL, FALSE, FALSE, hasNulls);
            if (gslp != NULL) {
              CheckSeqLocForPartial (sfp->location, &noLeft, &noRight);
              sfp->location = SeqLocFree (sfp->location);
              sfp->location = gslp;
              if (bsp->repr == Seq_repr_seg) {
                gslp = SegLocToParts (bsp, sfp->location);
                sfp->location = SeqLocFree (sfp->location);
                sfp->location = gslp;
              }
              FreeAllFuzz (sfp->location);
              SetSeqLocPartial (sfp->location, noLeft, noRight);
            }
          }
        }
      }
      if (CheckSeqLocForPartial (sfp->location, NULL, NULL)) {
        sfp->partial = TRUE;
      }
      sfp->cit = DialogToPointer (ffp->featcits);
      sfp->dbxref = DialogToPointer (ffp->dbxrefs);
      slp = AsnIoMemCopy (sfp->location, (AsnReadFunc) SeqLocAsnRead,
                          (AsnWriteFunc) SeqLocAsnWrite);
      usexref = GetValue (ffp->useGeneXref);
      if (ffp->gene != NULL) {
        val = GetValue (ffp->gene);
        if (usexref == 3 || (val > 1 && usexref == 2)) {
          grp = NULL;
          if (usexref == 3) {
            grp = GeneRefNew ();
          } else if (val == 2) {
            GetTitle (ffp->geneSymbol, symbol, sizeof (symbol));
            GetTitle (ffp->geneAllele, allele, sizeof (allele));
            GetTitle (ffp->geneDesc, desc, sizeof (desc));
            GetTitle (ffp->locusTag, locustag, sizeof (locustag));
            grp = CreateNewGeneRef (symbol, allele, desc, FALSE);
            if (! StringHasNoText (locustag)) {
              if (grp == NULL) {
                grp = GeneRefNew ();
              }
              grp->locus_tag = StringSave (locustag);
            }
          } else {
            vnp = ffp->geneNames;
            i = val - 3;
            while (i > 0 && vnp != NULL) {
              vnp = vnp->next;
              i--;
            }
            if (vnp != NULL) {
              if (vnp->choice == 1) {
                str = (CharPtr) vnp->data.ptrvalue;
                if (StringDoesHaveText (str)) {
                  grp = CreateNewGeneRef (str, NULL, NULL, FALSE);
                  gene = SeqMgrGetFeatureByLabel (bsp, str, SEQFEAT_GENE, 0, &fcontext);
                  if (gene != NULL && gene->data.choice == SEQFEAT_GENE) {
                    grpfeat = (GeneRefPtr) gene->data.value.ptrvalue;
                    if (grpfeat != NULL) {
                      grp->locus_tag = StringSaveNoNull (grpfeat->locus_tag);
                    }
                  }
                }
              } else if (vnp->choice == 2) {
                 grp = GeneRefNew ();
                if (grp != NULL) {
                  grp->desc = StringSave ((CharPtr) vnp->data.ptrvalue);
                }
             } else if (vnp->choice == 3) {
                grp = GeneRefNew ();
                if (grp != NULL) {
                  grp->locus_tag = StringSave ((CharPtr) vnp->data.ptrvalue);
                }
              }
            }
          }
          if (grp != NULL) {
            xref = SeqFeatXrefNew ();
            sfp->xref = xref;
            if (xref != NULL) {
              xref->data.choice = SEQFEAT_GENE;
              xref->data.value.ptrvalue = (Pointer) grp;
            }
          }
        }
      } else if (usexref == 3) {
        /* protein feature can now suppress gene on GenBank view */
        grp = GeneRefNew ();
        if (grp != NULL) {
          xref = SeqFeatXrefNew ();
          sfp->xref = xref;
          if (xref != NULL) {
            xref->data.choice = SEQFEAT_GENE;
            xref->data.value.ptrvalue = (Pointer) grp;
          }
        }
      }
      ompc.output_data = (Pointer) sfp;
      if (ompc.input_entityID == 0) {
        sfp->qual = DialogToPointer (ffp->gbquals);
        VisStringDialogToGbquals (sfp, ffp->experiment, "experiment");
        InferenceDialogToGBQuals (ffp->inference, sfp, TRUE);
        sfp->ext = DialogToPointer (ffp->usrobjext);
        if (ffp->goTermUserObj != NULL) {
          sfp->ext = CombineUserObjects (sfp->ext, ffp->goTermUserObj);
        }
        if (ffp->featid != NULL) {
          TextToFeatID (ffp->featid, &(sfp->id));
        }
        if (HasExceptionGBQual (sfp)) {
          sfp->excpt = TRUE;
        }
        if (ffp->fidxref != NULL) {
          TextToFeatXref (ffp->fidxref, sfp);
        }
        AddProtRefXref (sfp, ffp->protXrefName);
        if (! ObjMgrRegister (OBJ_SEQFEAT, (Pointer) sfp)) {
          Message (MSG_ERROR, "ObjMgrRegister failed");
        }
        SeqLocFree (slp);
        SeqEntrySetScope (oldscope);
        return TRUE;
      } else if (ompc.input_itemtype != OBJ_SEQFEAT) {
        sfp->qual = DialogToPointer (ffp->gbquals);
        VisStringDialogToGbquals (sfp, ffp->experiment, "experiment");
        InferenceDialogToGBQuals (ffp->inference, sfp, TRUE);
        sfp->ext = DialogToPointer (ffp->usrobjext);
        if (ffp->goTermUserObj != NULL) {
          sfp->ext = CombineUserObjects (sfp->ext, ffp->goTermUserObj);
        }
        if (ffp->featid != NULL) {
          TextToFeatID (ffp->featid, &(sfp->id));
        }
        if (HasExceptionGBQual (sfp)) {
          sfp->excpt = TRUE;
        }
        if (ffp->fidxref != NULL) {
          TextToFeatXref (ffp->fidxref, sfp);
        }
        AddProtRefXref (sfp, ffp->protXrefName);
        ompc.output_itemtype = OBJ_SEQFEAT;
        if (ompc.input_itemtype == OBJ_BIOSEQ) {
          bsp = GetBioseqGivenIDs (ompc.input_entityID, ompc.input_itemID, ompc.input_itemtype);
          if (bsp != NULL) {
            sap = bsp->annot;
            while (sap != NULL && (sap->name != NULL || sap->desc != NULL || sap->type != 1)) {
              sap = sap->next;
            }
            if (sap == NULL) {
              sap = SeqAnnotNew ();
              if (sap != NULL) {
                sap->type = 1;
                sap->next = bsp->annot;
                bsp->annot = sap;
              }
            }
            if (sap != NULL) {
              itemID = GetItemIDGivenPointer (ompc.input_entityID, OBJ_SEQANNOT, (Pointer) sap);
              if (itemID > 0) {
                ompc.input_itemID = itemID;
                ompc.input_itemtype = OBJ_SEQANNOT;
              }
            }
          }
        }
        if (! AttachDataForProc (&ompc, FALSE)) {
          Message (MSG_ERROR, "AttachDataForProc failed");
        }
        rsult = TRUE;
      } else {
        GatherItem (ompc.input_entityID, ompc.input_itemID, ompc.input_itemtype,
                    (Pointer) &old_location, GetOldFeatureLocation);
      
        rd.ffp = ffp;
        rd.sfp = sfp;
        GatherItem (ompc.input_entityID, ompc.input_itemID, ompc.input_itemtype,
                    (Pointer) &rd, ReplaceFeatureExtras);
        if (HasExceptionGBQual (sfp)) {
          sfp->excpt = TRUE;
        }
        if (ffp->fidxref != NULL) {
          TextToFeatXref (ffp->fidxref, sfp);
        }
        if (IsPseudo (sfp))
        {
          RemoveProtXrefs (sfp);
        }
        else
        {
          AddProtRefXref (sfp, ffp->protXrefName);
        }
        if (! ReplaceDataForProc (&ompc, FALSE)) {
          Message (MSG_ERROR, "ReplaceDataForProc failed");
        }
        rsult = TRUE;
      }
      if (ffp->gene != NULL && usexref == 1) {
        val = GetValue (ffp->gene);
        if (val == 2) {
          sep = GetBestTopParentForItemID (ffp->input_entityID,
                                           ffp->input_itemID,
                                           ffp->input_itemtype);
          /*
          sep = GetTopSeqEntryForEntityID (ffp->input_entityID);
          */
          if (sep != NULL) {
            sep = FindNucSeqEntry (sep);
          }
          if (sep != NULL && sep->data.ptrvalue != NULL) {
            GetTitle (ffp->geneSymbol, symbol, sizeof (symbol));
            GetTitle (ffp->geneAllele, allele, sizeof (allele));
            GetTitle (ffp->geneDesc, desc, sizeof (desc));
            GetTitle (ffp->locusTag, locustag, sizeof (locustag));
            grp = CreateNewGeneRef (symbol, allele, desc, FALSE);
            if (! StringHasNoText (locustag)) {
              if (grp == NULL) {
                grp = GeneRefNew ();
              }
              grp->locus_tag = StringSave (locustag);
            }
            if (grp != NULL) {
              sfp = CreateNewFeature (sep, NULL, SEQFEAT_GENE, NULL);
              if (sfp != NULL) {
                sfp->data.value.ptrvalue = (Pointer) grp;
                sfp->location = SeqLocFree (sfp->location);
                sfp->location = DialogToPointer (ffp->location);
                bsp = GetBioseqGivenSeqLoc (sfp->location, ffp->input_entityID);
                if (bsp != NULL) {
                  hasNulls = LocationHasNullsBetween (sfp->location);
                  gslp = SeqLocMerge (bsp, sfp->location, NULL, TRUE, FALSE, FALSE);
                  if (gslp != NULL) {
                    CheckSeqLocForPartial (sfp->location, &noLeft, &noRight);
                    sfp->location = SeqLocFree (sfp->location);
                    sfp->location = gslp;
                    if (bsp->repr == Seq_repr_seg) {
                      gslp = SegLocToPartsEx (bsp, sfp->location, TRUE);
                      sfp->location = SeqLocFree (sfp->location);
                      sfp->location = gslp;
                      hasNulls = LocationHasNullsBetween (sfp->location);
                      sfp->partial = (sfp->partial || hasNulls);
                    }
                    FreeAllFuzz (sfp->location);
                    SetSeqLocPartial (sfp->location, noLeft, noRight);
                    sfp->partial = (sfp->partial || noLeft || noRight);
                  }
                }
              }
            }
          }
        } else if (val > 2) {
          ggl.ffp = ffp;
          ggl.omp = ObjMgrGet ();
          ggl.slp = slp;
          ggl.genexref = NULL;
          ggl.xrefmatch = FALSE;
          ggl.idx = 2;
          ggl.val = val;
          ggl.min = INT4_MAX;
          ggl.old_feature_location = old_location;
          MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
          gs.seglevels = 1;
          gs.get_feats_location = TRUE;
	        MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
	        gs.ignore[OBJ_BIOSEQ] = FALSE;
	        gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
	        gs.ignore[OBJ_SEQFEAT] = FALSE;
	        gs.ignore[OBJ_SEQANNOT] = FALSE;
          gs.scope = GetBestTopParentForItemID (ffp->input_entityID,
                                                ffp->input_itemID,
                                                ffp->input_itemtype);
          GatherEntity (ffp->input_entityID, (Pointer) &ggl, GeneUpdateFunc, &gs);
        }
      }
      SeqLocFree (slp);
      SeqEntrySetScope (oldscope);
    }
  }
  old_location = SeqLocFree (old_location);
  return rsult;
}

extern void StdFeatFormActnProc (ForM f)

{
  FeatureFormPtr  ffp;

  if (FeatFormReplaceWithoutUpdateProc (f)) {
    ffp = (FeatureFormPtr) GetObjectExtra (f);
    if (ffp != NULL) {
      GetRidOfEmptyFeatsDescStrings (ffp->input_entityID, NULL);
      if (GetAppProperty ("InternalNcbiSequin") != NULL) {
        ExtendGeneFeatIfOnMRNA (ffp->input_entityID, NULL);
      }
      ObjMgrSetDirtyFlag (ffp->input_entityID, TRUE);
      ObjMgrSendMsg (OM_MSG_UPDATE, ffp->input_entityID,
                     ffp->input_itemID, ffp->input_itemtype);
    }
  }
}

extern void StdSeqFeatPtrToFeatFormProc (ForM f, Pointer data)

{
  FeatureFormPtr  ffp;
  SeqEntryPtr     oldsep;
  SeqEntryPtr     sep;
  SeqFeatPtr      sfp;
  Int4            val;

  ffp = (FeatureFormPtr) GetObjectExtra (f);
  if (ffp != NULL) {
    sep = GetTopSeqEntryForEntityID (ffp->input_entityID);
    oldsep = SeqEntrySetScope (sep);
    sfp = (SeqFeatPtr) data;
    if (sfp != NULL) {
      switch (sfp->data.choice) {
        case SEQFEAT_BOND :
        case SEQFEAT_SITE :
        case SEQFEAT_PSEC_STR :
          val = (Int4) sfp->data.value.intvalue;
          PointerToDialog (ffp->data, (Pointer) &(val));
          break;
        case SEQFEAT_COMMENT:
          break;
        default :
          PointerToDialog (ffp->data, sfp->data.value.ptrvalue);
          break;
      }
      SeqFeatPtrToCommon (ffp, sfp);
      PointerToDialog (ffp->location, sfp->location);
    }
    SeqEntrySetScope (oldsep);
  }
}

extern void StdInitFeatFormProc (ForM f)

{
  FeatureFormPtr  ffp;

  ffp = (FeatureFormPtr) GetObjectExtra (f);
  if (ffp != NULL) {
    PopulateGenePopup (ffp);
  }
}

extern void StdFeatFormAcceptButtonProc (ButtoN b)

{
  Int2            attempts = 3;
  Char            badInfMssg [32];
  Char            badInfQual [256];
  FeatureFormPtr  ffp;
  SeqLocPtr       slp;
  WindoW          w;
  ValNodePtr      err_list;
  SeqEntryPtr     oldscope, sep;

  if (b != NULL) {
    w = ParentWindow (b);
    ffp = (FeatureFormPtr) GetObjectExtra (b);
    if (ffp != NULL && ffp->form != NULL && ffp->actproc != NULL) {
      sep = GetTopSeqEntryForEntityID (ffp->input_entityID);
      oldscope = SeqEntrySetScope (sep);
      slp = DialogToPointer (ffp->location);
      if (slp == NULL) {
        ErrPostEx (SEV_ERROR, 0, 0, "Feature must have a location!");
        err_list = TestDialog (ffp->location);
        DisplayErrorMessages ("Location Errors", err_list);
        err_list = ValNodeFree (err_list);    
        SeqEntrySetScope (oldscope);
        return;
      }
      SeqLocFree (slp);
      if ((! ffp->acceptBadInf) && (! TestInference (ffp, badInfQual, sizeof (badInfQual), badInfMssg))) {
        (ffp->badInfAttempts)++;
        if (GetAppProperty ("InternalNcbiSequin") != NULL) {
          attempts = 2;
        }
        if (ffp->badInfAttempts < attempts) {
          if (ffp->badInfAttempts == 2) {
            Message (MSG_OK, "%s - Please fix the %s error in\n%s", infWarning2, badInfMssg, badInfQual);
          } else {
            Message (MSG_OK, "%s - Please fix the %s error in\n%s", infWarning1, badInfMssg, badInfQual);
          }
          return;
        } else {
          if (Message (MSG_YN, "%s", infAccept) == ANS_NO) return;
          ffp->acceptBadInf = TRUE;
        }
      }
      Hide (w);
      (ffp->actproc) (ffp->form);
      SeqEntrySetScope (oldscope);
    }
    Update ();
    if (ffp->leave_dlg_up == NULL || ! GetStatus (ffp->leave_dlg_up))
    {
      Remove (w);
    }
    else
    {
      /* set strand and sequence ID, but clear other location information */
      SendMessageToDialog (ffp->location, NUM_VIB_MSG + 1);
      SendMessageToDialog (ffp->location, NUM_VIB_MSG + 2);
      Show (w);
    }
  }
}

extern void StdFeatFormCleanupProc (GraphiC g, VoidPtr data)

{
  FeatureFormPtr  ffp;
  Uint2           userkey;

  ffp = (FeatureFormPtr) data;
  if (ffp != NULL) {
    ValNodeFreeData (ffp->geneNames);
    if (ffp->input_entityID > 0 && ffp->userkey > 0) {
      userkey = ffp->userkey;
      ffp->userkey = 0;
      ObjMgrFreeUserData (ffp->input_entityID, ffp->procid, ffp->proctype, userkey);
    }
  }
  StdCleanupFormProc (g, data);
}

extern ValNodePtr AddStringToValNodeChain (ValNodePtr head, CharPtr str, Uint1 choice)

{
  ValNodePtr  vnp;

  vnp = ValNodeNew (head);
  if (head == NULL) {
    head = vnp;
  }
  if (vnp != NULL) {
    vnp->choice = choice;
    vnp->data.ptrvalue = StringSave (str);
  }
  return head;
}

static void StripPeriods (CharPtr str)

{
  Char     ch;
  CharPtr  dst;

  if (str != NULL) {
    dst = str;
    ch = *str;
    while (ch != '\0') {
      if (ch != '.') {
        *dst = ch;
        dst++;
      }
      str++;
      ch = *str;
    }
    *dst = '\0';
  }
}

static void TrimLeadingSpaces (CharPtr str)

{
  Char     ch;
  CharPtr  dst;

  if (str != NULL && str [0] != '\0') {
    dst = str;
    ch = *str;
    while (ch != '\0' && ch <= ' ') {
      str++;
      ch = *str;
    }
    while (ch != '\0') {
      *dst = ch;
      dst++;
      str++;
      ch = *str;
    }
    *dst = '\0';
  }
}

extern CharPtr NameStdPtrToAuthorSpreadsheetString (NameStdPtr nsp);
extern CharPtr NameStdPtrToAuthorSpreadsheetString (NameStdPtr nsp)

{
  Char   first [256];
  Char   frstinits [64];
  Char   initials [64];
  Int2   j;
  Char   last [256];
  Char   middle [128];
  Char   str [512];
  Char   suffix [64];
  Char   suffixPosition[64];
  Int2   i;

  if (nsp == NULL) return NULL;
  str [0] = '\0';
  StringNCpy_0 (first, nsp->names [1], sizeof (first));
  TrimSpacesAroundString (first);
  StringNCpy_0 (initials, nsp->names [4], sizeof (initials));
  StripPeriods (initials);
  TrimLeadingSpaces (initials);
  StringNCpy_0 (last, nsp->names [0], sizeof (last));
  TrimLeadingSpaces (last);
  StringNCpy_0 (middle, nsp->names [2], sizeof (middle));
  TrimLeadingSpaces (middle);
  if (StringCmp (initials, "al") == 0 &&
      StringCmp (last, "et") == 0 &&
      first [0] == '\0') {
    initials [0] = '\0';
    StringCpy (last, "et al.");
  }
  if (first [0] == '\0') {
    StringNCpy_0 (first, initials, sizeof (first));
    if (IS_ALPHA (first [0])) {
      if (first [1] == '-') {
        first [3] = '\0';
      } else {
        first [1] = '\0';
      }
    } else {
      first [0] = '\0';
    }
  }
  FirstNameToInitials (first, frstinits, sizeof (frstinits) - 1);
  StripPeriods (first);
  TrimLeadingSpaces (first);
  if (first [0] != '\0') {
    StringCat (str, first);
  } else {
  }
  StringCat (str, "\t");
  j = 0;
  while (initials [j] != '\0' && TO_UPPER (initials [j]) == TO_UPPER (frstinits [j])) {
    j++;
  }
  if (initials [j] != '\0') {
    StringCat (str, initials + j);
  } else {
  }
  StringCat (str, "\t");
  StringCat (str, last);
  StringNCpy_0 (suffix, nsp->names [5], sizeof (suffix));
  for (i = 0; i <= NUMBER_OF_SUFFIXES; i++)
    if (StringICmp (suffix, name_suffix_labels [i]) == 0) {
      sprintf (suffixPosition, "%d", i);
      break;
    }
  if (i == NUMBER_OF_SUFFIXES)
    sprintf (suffixPosition, "%d", 0);
  StringCat (str, "\t");
  if (suffix [0] != '\0') {
    StringCat (str, suffixPosition);
  } else {
  }
  StringCat (str, "\t");
  StringCat (str, middle);
  StringCat (str, "\n");
  return StringSave (str);
}

static void StdAuthListPtrToAuthorDialog (DialoG d, Pointer data)

{
  AuthListPtr  alp;
  AuthorPtr    ap;
  ValNodePtr   head;
  Int2         j;
  ValNodePtr   names;
  NameStdPtr   nsp;
  PersonIdPtr  pid;
  TagListPtr   tlp;
  ValNodePtr   vnp;

  tlp = (TagListPtr) GetObjectExtra (d);
  alp = (AuthListPtr) data;
  if (tlp != NULL) {
    head = NULL;
    if (alp != NULL) {
      if (alp->choice == 1) {
        names = alp->names;
        while (names != NULL) {
          ap = names->data.ptrvalue;
          if (ap != NULL) {
            pid = ap->name;
            if (pid != NULL) {
              if (pid->choice == 2) {
                nsp = pid->data;
                if (nsp != NULL) {
                  vnp = ValNodeNew (head);
                  if (head == NULL) {
                    head = vnp;
                  }
                  if (vnp != NULL) {
                    vnp->data.ptrvalue = NameStdPtrToAuthorSpreadsheetString (nsp);
                  }
                }
              }
            }
          }
          names = names->next;
        }
      } else {
        Message (MSG_ERROR, "Unable to handle author type %d", (int) alp->choice);
      }
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

extern NameStdPtr AuthorSpreadsheetStringToNameStdPtr (CharPtr txt);
extern NameStdPtr AuthorSpreadsheetStringToNameStdPtr (CharPtr txt)

{
  Char        ch;
  CharPtr     first;
  Char        initials [64];
  Int2        j;
  Int2        k;
  Char        last;
  NameStdPtr  nsp;
  Char        periods [128];
  CharPtr     str;
  Char        str1 [64];
  CharPtr     suffix;
  Char        suffixVal [80];

  if (txt == NULL) return NULL;
  nsp = NameStdNew ();
  if (nsp == NULL) return NULL;
  nsp->names [0] = ExtractTagListColumn (txt, 2);
  TrimLeadingSpaces (nsp->names [0]);
  first = ExtractTagListColumn (txt, 0);
  StripPeriods (first);
  nsp->names [1] = StringSave (first);
  TrimLeadingSpaces (nsp->names [1]);
  FirstNameToInitials (first, str1, sizeof (str1) - 1);
  str = ExtractTagListColumn (txt, 1);
  StringNCat (str1, str, sizeof (str1) - 1);
  MemFree (str);
  j = 0;
  k = 0;
  ch = str1 [j];
  while (ch != '\0') {
    if (ch != ' ') {
      initials [k] = ch;
      k++;
    }
    j++;
    ch = str1 [j];
  }
  initials [k] = '\0';
  periods [0] = '\0';
          j = 0;
          ch = initials [j];
          while (ch != '\0') {
            if (ch == ',') {
              initials [j] = '.';
            }
            j++;
            ch = initials [j];
          }
          str = StringStr (initials, ".ST.");
          if (str != NULL) {
            *(str + 2) = 't';
          }
  j = 0;
  k = 0;
  ch = initials [j];
  while (ch != '\0') {
    if (ch == '-') {
      periods [k] = ch;
      k++;
      j++;
      ch = initials [j];
    } else if (ch == '.') {
      j++;
      ch = initials [j];
            } else if (ch == ' ') {
              j++;
              ch = initials [j];
    } else {
      periods [k] = ch;
              last = ch;
      k++;
      j++;
      ch = initials [j];
              if (ch == '\0') {
                if (! (IS_LOWER (last))) {
                  periods [k] = '.';
                  k++;
                }
              /* } else if (ch == '.' && initials [j + 1] == '\0') { */
              } else if (! (IS_LOWER (ch))) {
                periods [k] = '.';
                k++;
              }
    }
  }
  periods [k] = '\0';
  nsp->names [4] = StringSave (periods);
  TrimLeadingSpaces (nsp->names [4]);
  str = ExtractTagListColumn (txt, 3);
  StringNCpy_0 (str1, str, sizeof (str1));
  MemFree (str);
  j = 0;
  k = 0;
  ch = str1 [j];
  while (ch != '\0') {
    if (ch != ' ') {
      suffixVal [k] = ch;
      k++;
    }
    j++;
    ch = str1 [j];
  }
  suffixVal [k] = '\0';
  if (suffixVal [0] != '\0') {
    suffix = GetEnumName (atoi(suffixVal), name_suffix_alist);
    nsp->names [5] = StringSave (suffix);
    TrimLeadingSpaces (nsp->names [5]);
    if (StringHasNoText (nsp->names [5])) {
      nsp->names [5] = MemFree (nsp->names [5]);
    }
  }
  if (StringCmp (nsp->names [0], "et al") == 0) {
    nsp->names [0] = MemFree (nsp->names [0]);
    nsp->names [0] = StringSave ("et al.");
  }
  nsp->names [2] = ExtractTagListColumn (txt, 4);
  TrimLeadingSpaces (nsp->names [2]);
  if (StringHasNoText (nsp->names [2])) {
    nsp->names [2] = MemFree (nsp->names [2]);
  }
  return nsp;
}

static Pointer AuthorDialogToStdAuthListPtr (DialoG d)

{
  AuthListPtr  alp;
  AuthorPtr    ap;
  Char         ch;
  Int2         j;
  Int2         len;
  ValNodePtr   names;
  NameStdPtr   nsp;
  Boolean      okay;
  PersonIdPtr  pid;
  CharPtr      str;
  TagListPtr   tlp;
  ValNodePtr   vnp;

  alp = NULL;
  tlp = (TagListPtr) GetObjectExtra (d);
  if (tlp != NULL && tlp->vnp != NULL) {
    alp = AuthListNew ();
    if (alp != NULL) {
      alp->choice = 1;
      names = NULL;
      for (vnp = tlp->vnp; vnp != NULL; vnp = vnp->next) {
        str = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, 2);
        okay = FALSE;
        len = StringLen (str);
        for (j = 0; j < len; j++) {
          ch = str [j];
          if (ch != ' ' && ch != '\t' && ch != '\n') {
            okay = TRUE;
          }
        }
        MemFree (str);
        if (okay) {
          names = ValNodeNew (names);
          if (alp->names == NULL) {
            alp->names = names;
          }
          if (names != NULL) {
            ap = AuthorNew ();
            names->choice = 1;
            names->data.ptrvalue = ap;
            if (ap != NULL) {
              pid = PersonIdNew ();
              ap->name = pid;
              if (pid != NULL) {
                pid->choice = 2;
                nsp = AuthorSpreadsheetStringToNameStdPtr ((CharPtr) vnp->data.ptrvalue);
                pid->data = nsp;
              }
            }
          }
        }
      }
      if (alp->names == NULL) {
        alp = AuthListFree (alp);
      }
    }
  }
  return (Pointer) alp;
}

static void StrAuthListPtrToAuthorDialog (DialoG d, Pointer data)

{
  AuthListPtr  alp;
  ValNodePtr   head;
  Int2         j;
  ValNodePtr   names;
  Char         str [128];
  TagListPtr   tlp;
  ValNodePtr   vnp;

  tlp = (TagListPtr) GetObjectExtra (d);
  alp = (AuthListPtr) data;
  if (tlp != NULL) {
    head = NULL;
    if (alp != NULL) {
      if (alp->choice == 2 || alp->choice == 3) {
        names = alp->names;
        while (names != NULL) {
          StringNCpy_0 (str, names->data.ptrvalue, sizeof (str) - 2);
          StringCat (str, "\n");
          vnp = ValNodeNew (head);
          if (head == NULL) {
            head = vnp;
          }
          if (vnp != NULL) {
            vnp->data.ptrvalue = StringSave (str);
          }
          names = names->next;
        }
      } else {
        Message (MSG_ERROR, "Unable to handle author type %d", (int) alp->choice);
      }
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

static Pointer AuthorDialogToStrAuthListPtr (DialoG d)

{
  AuthListPtr  alp;
  Char         ch;
  Int2         j;
  Int2         len;
  ValNodePtr   names;
  Boolean      okay;
  CharPtr      str;
  TagListPtr   tlp;
  ValNodePtr   vnp;

  alp = NULL;
  tlp = (TagListPtr) GetObjectExtra (d);
  if (tlp != NULL && tlp->vnp != NULL) {
    alp = AuthListNew ();
    if (alp != NULL) {
      alp->choice = 2;
      names = NULL;
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
          names = ValNodeNew (names);
          if (alp->names == NULL) {
            alp->names = names;
          }
          if (names != NULL) {
            names->choice = 2;
            names->data.ptrvalue = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, 0);

          }
        }
      }
      if (alp->names == NULL) {
        alp = AuthListFree (alp);
      }
    }
  }
  return (Pointer) alp;
}

static void AuthListPtrToAuthorDialog (DialoG d, Pointer data)

{
  AuthorDialogPtr  adp;
  AuthListPtr      alp;

  adp = (AuthorDialogPtr) GetObjectExtra (d);
  if (adp != NULL) {
    alp = (AuthListPtr) data;
    if (alp != NULL) {
      adp->type = alp->choice;
    }
    if (adp->type == 1) {
      Hide (adp->strGrp);
      Show (adp->stdGrp);
      PointerToDialog (adp->stdAuthor, data);
    } else if (adp->type == 2 || adp->type == 3) {
      Hide (adp->stdGrp);
      Show (adp->strGrp);
      PointerToDialog (adp->strAuthor, data);
    }
  }
}

static Pointer AuthorDialogToAuthListPtr (DialoG d)

{
  AuthorDialogPtr  adp;
  AuthListPtr      alp;

  adp = (AuthorDialogPtr) GetObjectExtra (d);
  alp = NULL;
  if (adp != NULL) {
    if (adp->type == 1) {
      alp = (AuthListPtr) DialogToPointer (adp->stdAuthor);
    } else if (adp->type == 2 || adp->type == 3) {
      if ((alp = (AuthListPtr) DialogToPointer (adp->strAuthor)) != NULL)
        alp->choice = adp->type;
    }
  }
  return (Pointer) alp;
}

static void AuthorDialogMessage (DialoG d, Int2 mssg)

{
  AuthorDialogPtr  adp;

  adp = (AuthorDialogPtr) GetObjectExtra (d);
  if (adp != NULL) {
    switch (mssg) {
      case VIB_MSG_ENTER :
        if (adp->type == 1) {
          SendMessageToDialog (adp->stdAuthor, VIB_MSG_ENTER);
        } else if (adp->type == 2 || adp->type == 3) {
          SendMessageToDialog (adp->strAuthor, VIB_MSG_ENTER);
        }
        break;
      default :
        break;
    }
  }
}

static Boolean ReadAuthorDialog (DialoG d, CharPtr filename)

{
  AuthorDialogPtr  adp;
  AsnIoPtr         aip;
  AuthListPtr      alp;
  Char             path [PATH_MAX];

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  adp = (AuthorDialogPtr) GetObjectExtra (d);
  if (adp != NULL) {
    if (path [0] != '\0' || GetInputFileName (path, sizeof (path), "", "TEXT")) {
      aip = AsnIoOpen (path, "r");
      if (aip != NULL) {
        alp = AuthListAsnRead (aip, NULL);
        AsnIoClose (aip);
        if (alp != NULL) {
          PointerToDialog (adp->dialog, (Pointer) alp);
          alp = AuthListFree (alp);
          Update ();
          return TRUE;
        }
      }
    }
  }
  return FALSE;
}

static Boolean WriteAuthorDialog (DialoG d, CharPtr filename)

{
  AuthorDialogPtr  adp;
  AsnIoPtr         aip;
  AuthListPtr      alp;
  Char             path [PATH_MAX];
#ifdef WIN_MAC
  FILE            *f;
#endif

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  adp = (AuthorDialogPtr) GetObjectExtra (d);
  if (adp != NULL) {
    if (path [0] != '\0' || GetOutputFileName (path, sizeof (path), NULL)) {
#ifdef WIN_MAC
      f = FileOpen (path, "r");
      if (f != NULL) {
        FileClose (f);
      } else {
        FileCreate (path, "TEXT", "ttxt");
      }
#endif
      aip = AsnIoOpen (path, "w");
      if (aip != NULL) {
        alp = DialogToPointer (adp->dialog);
        AuthListAsnWrite (alp, aip, NULL);
        AsnIoClose (aip);
        alp = AuthListFree (alp);
        return TRUE;
      }
    }
  }
  return FALSE;
}

extern DialoG CreateAuthorDialog (GrouP prnt, Uint2 rows, Int2 spacing)

{
  AuthorDialogPtr  adp;
  GrouP            k;
  GrouP            p;
  PrompT           p1, p2, p3, p4, p5;
  TagListPtr       tlp;

  p = HiddenGroup (prnt, 0, 0, NULL);

  adp = (AuthorDialogPtr) MemNew (sizeof (AuthorDialog));
  if (adp != NULL) {

    SetObjectExtra (p, adp, StdCleanupExtraProc);
    adp->dialog = (DialoG) p;
    adp->todialog = AuthListPtrToAuthorDialog;
    adp->fromdialog = AuthorDialogToAuthListPtr;
    adp->testdialog = NULL;
    adp->importdialog = ReadAuthorDialog;
    adp->exportdialog = WriteAuthorDialog;
    adp->dialogmessage = AuthorDialogMessage;

    adp->strGrp = HiddenGroup (p, -1, 0, NULL);
    SetGroupSpacing (adp->strGrp, 3, 2);

    k = HiddenGroup (adp->strGrp, -4, 0, NULL);
    SetGroupSpacing (k, spacing, spacing);
    p1 = StaticPrompt (k, "Name", 0, 0, programFont, 'c');

    adp->strAuthor = CreateTagListDialog (adp->strGrp, rows, 1, spacing,
                                          author_types, str_author_widths, NULL,
                                          StrAuthListPtrToAuthorDialog,
                                          AuthorDialogToStrAuthListPtr);
    Hide (adp->strGrp);

    adp->stdGrp = HiddenGroup (p, -1, 0, NULL);
    SetGroupSpacing (adp->stdGrp, 3, 2);

    k = HiddenGroup (adp->stdGrp, -4, 0, NULL);
    SetGroupSpacing (k, spacing, spacing);
    p2 = StaticPrompt (k, "First Name", 0, 0, programFont, 'c');
    p3 = StaticPrompt (k, "M.I.", 0, 0, programFont, 'c');
    p4 = StaticPrompt (k, "Last Name", 0, 0, programFont, 'c');
    p5 = StaticPrompt (k, "Sfx", 0, 0, programFont, 'c');

    adp->stdAuthor = CreateTagListDialog (adp->stdGrp, rows, 4, spacing,
                                          author_types, std_author_widths,
                                          author_popups,
                                          StdAuthListPtrToAuthorDialog,
                                          AuthorDialogToStdAuthListPtr);
    adp->type = 1;

    tlp = (TagListPtr) GetObjectExtra (adp->strAuthor);
    if (tlp != NULL) {
      AlignObjects (ALIGN_JUSTIFY, (HANDLE) tlp->control [0], (HANDLE) p1, NULL);
    }
    tlp = (TagListPtr) GetObjectExtra (adp->stdAuthor);
    if (tlp != NULL) {
      AlignObjects (ALIGN_JUSTIFY, (HANDLE) tlp->control [0], (HANDLE) p2, NULL);
      AlignObjects (ALIGN_JUSTIFY, (HANDLE) tlp->control [1], (HANDLE) p3, NULL);
      AlignObjects (ALIGN_JUSTIFY, (HANDLE) tlp->control [2], (HANDLE) p4, NULL);
      AlignObjects (ALIGN_JUSTIFY, (HANDLE) tlp->control [3], (HANDLE) p5, NULL);
    }
  }

  return (DialoG) p;
}

typedef struct affildialog {
  DIALOG_MESSAGE_BLOCK
  TexT            affil;
  TexT            div;
  TexT            address;
  TexT            city;
  TexT            state;
  TexT            zip;
  TexT            country;
  TexT            phone;
  TexT            fax;
  TexT            email;
} AffilDialog, PNTR AffilDialogPtr;

static void AffilPtrToAffilDialog (DialoG d, Pointer data)

{
  AffilDialogPtr  adp;
  AffilPtr        ap;

  adp = (AffilDialogPtr) GetObjectExtra (d);
  ap = (AffilPtr) data;
  if (adp != NULL) {
    if (ap != NULL) {
      SafeSetTitle (adp->affil, ap->affil);
      SafeSetTitle (adp->div, ap->div);
      SafeSetTitle (adp->address, ap->street);
      SafeSetTitle (adp->city, ap->city);
      SafeSetTitle (adp->state, ap->sub);
      SafeSetTitle (adp->zip, ap->postal_code);
      SafeSetTitle (adp->country, ap->country);
      SafeSetTitle (adp->phone, ap->phone);
      SafeSetTitle (adp->fax, ap->fax);
      SafeSetTitle (adp->email, ap->email);
    } else {
      SafeSetTitle (adp->affil, "");
      SafeSetTitle (adp->div, "");
      SafeSetTitle (adp->address, "");
      SafeSetTitle (adp->city, "");
      SafeSetTitle (adp->state, "");
      SafeSetTitle (adp->zip, "");
      SafeSetTitle (adp->country, "");
      SafeSetTitle (adp->phone, "");
      SafeSetTitle (adp->fax, "");
      SafeSetTitle (adp->email, "");
    }
  }
}

static Pointer AffilDialogToAffilPtr (DialoG d)

{
  AffilDialogPtr  adp;
  AffilPtr        ap;

  ap = NULL;
  adp = (AffilDialogPtr) GetObjectExtra (d);
  if (adp != NULL) {
    ap = AffilNew ();
    if (ap != NULL) {
      ap->affil = SaveStringFromText (adp->affil);
      ap->div = SaveStringFromText (adp->div);
      ap->street = SaveStringFromText (adp->address);
      ap->city = SaveStringFromText (adp->city);
      ap->sub = SaveStringFromText (adp->state);
      ap->postal_code = SaveStringFromText (adp->zip);
      ap->country = SaveStringFromText (adp->country);
      ap->phone = SaveStringFromText (adp->phone);
      ap->fax = SaveStringFromText (adp->fax);
      ap->email = SaveStringFromText (adp->email);
      if (ap->div == NULL && ap->street == NULL && ap->city == NULL &&
           ap->sub == NULL && ap->postal_code == NULL && ap->country == NULL &&
           ap->phone == NULL && ap->fax == NULL && ap->email == NULL) {
        ap->choice = 1;
        if (ap->affil == NULL) {
          ap = AffilFree (ap);
        }
      } else {
        ap->choice = 2;
      }
    }
  }
  return (Pointer) ap;
}

static void AffilDialogMessage (DialoG d, Int2 mssg)

{
  AffilDialogPtr  adp;

  adp = (AffilDialogPtr) GetObjectExtra (d);
  if (adp != NULL) {
    switch (mssg) {
      case VIB_MSG_ENTER :
        Select (adp->affil);
        break;
      default :
        break;
    }
  }
}

static DialoG CreateAnAffilDialog (GrouP prnt, CharPtr title,
                                   Boolean publisher,
                                   Boolean split,
                                   Boolean proceedings,
                                   GrouP PNTR grp1, GrouP PNTR grp2)

{
  AffilDialogPtr  adp;
  GrouP           g;
  GrouP           g1, g2;
  GrouP           j;
  GrouP           m;
  GrouP           p;
  GrouP           q;
  GrouP           s;
#ifdef WIN_MAC
  Int2            wid = 20;
  Int2            zipw = 6;
  Int2            ewid = 8;
#else
  Int2            wid = 30;
  Int2            zipw = 10;
  Int2            ewid = 20;
#endif

  p = HiddenGroup (prnt, 0, 0, NULL);

  adp = (AffilDialogPtr) MemNew (sizeof (AffilDialog));
  if (adp != NULL) {

    SetObjectExtra (p, adp, StdCleanupExtraProc);
    adp->dialog = (DialoG) p;
    adp->todialog = AffilPtrToAffilDialog;
    adp->fromdialog = AffilDialogToAffilPtr;
    adp->testdialog = NULL;
    adp->dialogmessage = AffilDialogMessage;

    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    m = HiddenGroup (s, 0, 0, NULL);

    q = HiddenGroup (m, 2, 0, NULL);
    g1 = q;
    g2 = NULL;
    if (grp1 != NULL) {
      *grp1 = q;
    }
    if (grp2 != NULL) {
      *grp2 = NULL;
    }
    g = HiddenGroup (q, 0, 20, NULL);
    if (publisher) {
      StaticPrompt (g, "Publisher", 0, dialogTextHeight, programFont, 'l');
    } else if (proceedings) {
      StaticPrompt (g, "Location", 0, dialogTextHeight, programFont, 'l');
    } else {
      StaticPrompt (g, "Institution", 0, dialogTextHeight, programFont, 'l');
      StaticPrompt (g, "Department", 0, dialogTextHeight, programFont, 'l');
    }
    StaticPrompt (g, "Address", 0, dialogTextHeight, programFont, 'l');
    StaticPrompt (g, "City", 0, dialogTextHeight, programFont, 'l');
    StaticPrompt (g, "State/Province", 0, dialogTextHeight, programFont, 'l');
    StaticPrompt (g, "Country", 0, dialogTextHeight, programFont, 'l');
    if (! split) {
      if (! proceedings) {
        StaticPrompt (g, "", 0, stdLineHeight, programFont, 'l');
        StaticPrompt (g, "Phone", 0, dialogTextHeight, programFont, 'l');
        if (publisher) {
          StaticPrompt (g, "Internet Access", 0, dialogTextHeight, programFont, 'l');
        } else {
          StaticPrompt (g, "Email", 0, dialogTextHeight, programFont, 'l');
        }
      }
    }
    j = HiddenGroup (q, 0, 20, NULL);
    g = HiddenGroup (j, 0, 20, NULL);
    adp->affil = DialogText (g, "", wid, NULL);
    adp->div = NULL;
    if (! publisher && ! proceedings) {
      adp->div = DialogText (g, "", wid, NULL);
    }
    adp->address = DialogText (g, "", wid, NULL);
    adp->city = DialogText (g, "", wid, NULL);
    g = HiddenGroup (j, 3, 0, NULL);
    SetGroupSpacing (g, 20, 2);
    adp->state = DialogText (g, "", 10, NULL);
    if (! proceedings) {
      /* StaticPrompt (g, "Zip/Postal Code", 7 * stdCharWidth,
                    dialogTextHeight, programFont, 'l'); */
      StaticPrompt (g, "Zip/Postal Code", 0, dialogTextHeight, programFont, 'l');
      adp->zip = DialogText (g, "", zipw, NULL);
    }
    g = HiddenGroup (j, 0, 20, NULL);
    adp->country = DialogText (g, "", wid, NULL);
    if (split) {
      if (! proceedings) {
        q = HiddenGroup (m, 2, 0, NULL);
        g2 = q;
        if (grp2 != NULL) {
          *grp2 = q;
        }
        g = HiddenGroup (q, 0, 20, NULL);
        StaticPrompt (g, "", 0, stdLineHeight, programFont, 'l');
        StaticPrompt (g, "Phone", 0, dialogTextHeight, programFont, 'l');
        if (publisher) {
          StaticPrompt (g, "Internet Access", 0, dialogTextHeight, programFont, 'l');
        } else {
          StaticPrompt (g, "Email", 0, dialogTextHeight, programFont, 'l');
        }
        j = HiddenGroup (q, 0, 20, NULL);
      }
    }
    if (! proceedings) {
      StaticPrompt (j, "Please include country code for non-U.S. phone numbers.",
                    0, stdLineHeight, programFont, 'l');
      g = HiddenGroup (j, 3, 0, NULL);
      SetGroupSpacing (g, 20, 2);
      adp->phone = DialogText (g, "", 10, NULL);
      if (split) {
        StaticPrompt (g, "Fax", 0,
                      dialogTextHeight, programFont, 'l');
      } else {
        StaticPrompt (g, "Fax", /* 7 * stdCharWidth */ 0,
                      dialogTextHeight, programFont, 'l');
      }
      adp->fax = DialogText (g, "", 10, NULL);
      g = HiddenGroup (j, 0, 20, NULL);
      if (split) {
        adp->email = DialogText (g, "", ewid, NULL);
      } else {
        adp->email = DialogText (g, "", ewid, NULL);
      }
    }

    if (split) {
      AlignObjects (ALIGN_RIGHT, (HANDLE) adp->affil,
                    (HANDLE) adp->address, (HANDLE) adp->city,
                    (HANDLE) adp->zip, (HANDLE) adp->country,
                    (HANDLE) adp->div, NULL);
      AlignObjects (ALIGN_RIGHT, (HANDLE) adp->fax, (HANDLE) adp->email, NULL);
      AlignObjects (ALIGN_CENTER, (HANDLE) g1, (HANDLE) g2, NULL);
      Hide (g1);
      Hide (g2);
    } else {
      AlignObjects (ALIGN_RIGHT, (HANDLE) adp->affil,
                    (HANDLE) adp->address, (HANDLE) adp->city,
                    (HANDLE) adp->zip, (HANDLE) adp->country,
                    (HANDLE) adp->fax, (HANDLE) adp->email,
                    (HANDLE) adp->div, NULL);
    }
  }

  return (DialoG) p;
}

extern DialoG CreateAffilDialog (GrouP prnt, CharPtr title)

{
  return CreateAnAffilDialog (prnt, title, FALSE, FALSE, FALSE, NULL, NULL);
}

extern DialoG CreatePublisherAffilDialog (GrouP prnt, CharPtr title)

{
  return CreateAnAffilDialog (prnt, title, TRUE, FALSE, FALSE, NULL, NULL);
}

extern DialoG CreateProceedingsDialog (GrouP prnt, CharPtr title)

{
  return CreateAnAffilDialog (prnt, title, FALSE, FALSE, TRUE, NULL, NULL);
}

extern DialoG CreateExtAffilDialog (GrouP prnt, CharPtr title, GrouP PNTR grp1, GrouP PNTR grp2)

{
  return CreateAnAffilDialog (prnt, title, FALSE, TRUE, FALSE, grp1, grp2);
}

extern DialoG CreateExtPublisherAffilDialog (GrouP prnt, CharPtr title, GrouP PNTR grp1, GrouP PNTR grp2)

{
  return CreateAnAffilDialog (prnt, title, TRUE, TRUE, FALSE, grp1, grp2);
}

extern DialoG CreateExtProceedingsDialog (GrouP prnt, CharPtr title, GrouP PNTR grp1, GrouP PNTR grp2)

{
  return CreateAnAffilDialog (prnt, title, FALSE, TRUE, TRUE, grp1, grp2);
}

extern void DatePtrToVibrant (DatePtr dp, PopuP dateMonth, TexT dateDay, TexT dateYear)

{
  Int2  day;
  Char  str [32];
  Int2  year;

  if (dp != NULL) {
    if (dp->data [0] == 0) {
      DatePrint (dp, str);
    } else if (dp->data [0] == 1) {
      SetEnumPopup (dateMonth, months_alist, (UIEnum) dp->data [2]);
      day = (Int2) dp->data [3];
      if (day > 0 && day <= 31) {
        sprintf (str, "%d", (int) day);
        SafeSetTitle (dateDay, str);
      } else {
        SafeSetTitle (dateDay, "");
      }
      year = (Int2) dp->data [1];
      if (year > 0) {
        sprintf (str, "%d", (int) (year + 1900));
        SafeSetTitle (dateYear, str);
      } else {
        SafeSetTitle (dateYear, "");
      }
    } else {
      Message (MSG_ERROR, "Unknown date type");
    }
  } else {
    SafeSetValue (dateMonth, 1);
    SafeSetTitle (dateDay, "");
    SafeSetTitle (dateYear, "");
  }
}

extern DatePtr VibrantToDatePtr (PopuP dateMonth, TexT dateDay, TexT dateYear)

{
  Int2     day;
  Int2     dateType;
  DatePtr  dp;
  UIEnum   month;
  Char     str [32];
  Int2     year;

  dp = DateNew ();
  if (dp != NULL) {
    dateType = 1;
    dp->data [0] = (Uint1) dateType;
    if (dateType == 0) {
    } else if (dateType == 1) {
      GetTitle (dateYear, str, sizeof (str));
      if (! StringHasNoText (str)) {
        StrToInt (str, &year);
        if (year >= 1900) {
          dp->data [1] = (Uint1) (year - 1900);
        } else {
          /* dp->data [1] = 0; */
          dp = DateFree (dp);
          return dp;
        }
        if (GetEnumPopup (dateMonth, months_alist, &month)) {
          dp->data [2] = (Uint1) month;
        }
        GetTitle (dateDay, str, sizeof (str));
        StrToInt (str, &day);
        dp->data [3] = (Uint1) day;
      } else {
        dp = DateFree (dp);
      }
    } else {
      Message (MSG_ERROR, "Unknown date type");
    }
  }
  return dp;
}

static void GBQualPtrToQualsDialog (DialoG d, Pointer data)

{
  ValNodePtr  head;
  Int2        j;
  size_t      len;
  GBQualPtr   list;
  CharPtr     str;
  TagListPtr  tlp;
  ValNodePtr  vnp;

  tlp = (TagListPtr) GetObjectExtra (d);
  list = (GBQualPtr) data;
  if (tlp != NULL) {
    head = NULL;
    while (list != NULL) {
      vnp = ValNodeNew (head);
      if (head == NULL) {
        head = vnp;
      }
      if (vnp != NULL) {
        if (list->qual != NULL && list->val != NULL) {
          len = StringLen (list->qual) + StringLen (list->val);
          str = MemNew (len + 4);
          if (str != NULL) {
            StringCpy (str, list->qual);
            StringCat (str, "\t");
            StringCat (str, list->val);
            StringCat (str, "\n");
          }
          vnp->data.ptrvalue = str;
        }
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

static Pointer QualsDialogToGBQualPtr (DialoG d)

{
  Char         ch;
  GBQualPtr    gbq;
  GBQualPtr    gbqlast;
  GBQualPtr    head;
  Int2         j;
  Int2         len;
  Boolean      okay;
  CharPtr      str;
  TagListPtr   tlp;
  ValNodePtr   vnp;

  head = NULL;
  tlp = (TagListPtr) GetObjectExtra (d);
  if (tlp != NULL && tlp->vnp != NULL) {
    gbq = NULL;
    gbqlast = NULL;
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
        gbq = GBQualNew ();
        if (gbqlast == NULL) {
          head = gbq;
        } else {
          gbqlast->next = gbq;
        }
        gbqlast = gbq;
        if (gbq != NULL) {
          gbq->qual = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, 0);
          gbq->val = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, 1);
        }
      }
    }
  }
  return (Pointer) head;
}

Uint2 gbqual_types [] = {
  TAGLIST_TEXT, TAGLIST_TEXT
};

Uint2 gbqual_widths [] = {
  0, 0, 0, 0
};

extern DialoG CreateQualsDialog (GrouP h, Uint2 rows, Int2 spacing,
                                 Int2 width1, Int2 width2);
extern DialoG CreateQualsDialog (GrouP h, Uint2 rows, Int2 spacing,
                                 Int2 width1, Int2 width2)

{
  gbqual_widths [0] = width1;
  gbqual_widths [1] = width2;
  return CreateTagListDialog (h, rows, 2, spacing,
                              gbqual_types, gbqual_widths, NULL,
                              GBQualPtrToQualsDialog,
                              QualsDialogToGBQualPtr);
}

static void CreateSeqAlignLabel (SeqAlignPtr salp, CharPtr buf, Int4 buf_size)
{
  Int4     remaining_len, aln_pos, id_len;
  SeqIdPtr sip;
  CharPtr  buf_ptr;
  
  if (buf == NULL || buf_size < 1)
  {
    return;
  }
  MemSet (buf, 0, buf_size);
  if (salp == NULL)
  {
    return;
  }
  remaining_len = buf_size - 1;
  buf_ptr = buf;
  StringNCat (buf_ptr, "aln|", remaining_len);
  remaining_len -= 4;
  
  if (remaining_len <= 3)
  {
    return;
  }
  buf_ptr += 4;
  
  for (aln_pos = 1; aln_pos <= salp->dim && remaining_len > 0; aln_pos++)
  {
    sip = AlnMgr2GetNthSeqIdPtr(salp, aln_pos);
    SeqIdWrite (sip, buf_ptr, PRINTID_REPORT, remaining_len);
    id_len = StringLen (buf_ptr);
    remaining_len -= id_len;
    buf_ptr += id_len;
    /* put comma between IDs in list */
    if (aln_pos < salp->dim && remaining_len > 0)
    {
      StringCat (buf_ptr, ",");
      remaining_len -= 1;
      buf_ptr += 1;
    }
  }
  
  /* add ellipsis to indicate unshown IDs */
  if (remaining_len == 0 && aln_pos < salp->dim)
  {
    buf_ptr -= 3;
    StringCat (buf_ptr, "...");
  }
}

static void GetAlignmentsInSeqEntryCallback (SeqAnnotPtr sap, Pointer userdata)
{
  SeqAlignPtr PNTR salp_list;
  SeqAlignPtr salp, last_salp;
  
  if (sap == NULL || sap->type != 2 || userdata == NULL) 
  {
    return;
  }
  salp_list = (SeqAlignPtr PNTR) userdata;
  salp = (SeqAlignPtr) sap->data;
  if (salp == NULL) return;
  salp = SeqAlignListDup(salp);
  AlnMgr2IndexSeqAlign(salp);
  if (*salp_list == NULL)
  {
    *salp_list = salp; 
  }
  else
  {
    last_salp = *salp_list;
    while (last_salp->next != NULL)
    {
      last_salp = last_salp->next;
    }
    last_salp->next = salp;
  }
}

typedef struct intervalpage {
  DIALOG_MESSAGE_BLOCK
  DialoG             ivals;
  ButtoN             partial5;
  ButtoN             partial3;
  ButtoN             nullsBetween;
  Int2               count;
  SeqIdPtr           PNTR sip_list;
  EnumFieldAssoc     PNTR alist;
  EnumFieldAssocPtr  alists [5];
  Int4               PNTR lengths;
  Boolean            nucsOK;
  Boolean            protsOK;
  Boolean            showIdTags;
  FeatureFormPtr     ffp;
  IntEdPartialProc   proc;
  
  Int4               strand_col;  /* column for selecting strand.
                                   * 2 if nucsOK, -1 otherwise
                                   */
  Int4               seqid_col;   /* column for entering SeqIds,
                                   * which could be in a different place 
                                   * depending on nucsOK or -1 if not show_seqid
                                   */
  Int4               aln_col;     /* column for selecting alignments,
                                   * could be in a different place if not show_seqid
                                   * or -1 if not use_aln
                                   */
  Boolean            show_seqids; /* false when entering coordinates for an entire
                                   * alignments, true otherwise
                                   */
  
  /* for editing alignment intervals */
  SeqAlignPtr        PNTR salp_list; 
  TaglistCallback    PNTR callbacks;
  EnumFieldAssoc     PNTR aln_alist;
  Int2               aln_count;
  Int4               PNTR aln_lengths;
  Boolean            allow_nulls_in_list; /* some sequences may be all gap
                                           * in the specified interval -
                                           * if this value is FALSE, return
                                           * a NULL for the entire list,
                                           * otherwise allow NULLs to appear
                                           * in the list.
                                           */
  
} IntervalPage, PNTR IntervalPagePtr;

#define NUM_IVAL_ROWS  7
#define EXTRA_HEIGHT   2

static void ClearBspScratch (BioseqPtr bsp, Pointer userdata)

{
  if (bsp == NULL) return;
  bsp->idx.scratch = NULL;
}

static void AddToSipList (IntervalPagePtr ipp, BioseqPtr bsp)

{
  /*
  Int2      j;
  */
  SeqIdPtr  sip;

  if (ipp == NULL || bsp == NULL) return;
  if (bsp->idx.scratch != NULL) return;
  sip = SeqIdFindBest (bsp->id, 0);
  if (sip == NULL) return;
  /*
  for (j = 1; j <= ipp->count; j++) {
    if (SeqIdComp (sip, ipp->sip_list [j]) == SIC_YES) return;
  }
  */
  ipp->count++;
  ipp->sip_list [ipp->count] = SeqIdDup (sip);
  bsp->idx.scratch = (Pointer) bsp;
}

static void FillInProducts (SeqEntryPtr sep, Pointer mydata,
                            Int4 index, Int2 indent)

{
  BioseqPtr        bsp;
  IntervalPagePtr  ipp;
  SeqIdPtr         sip;
  SeqLocPtr        slp;
  ValNode          vn;

  if (sep != NULL && mydata != NULL && sep->choice == 1) {
    ipp = (IntervalPagePtr) mydata;
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp != NULL) {
      if ((ipp->nucsOK && ISA_na (bsp->mol)) ||
          (ipp->protsOK && ISA_aa (bsp->mol))) {
        AddToSipList (ipp, bsp);
      }
      if (bsp->repr == Seq_repr_seg && bsp->seq_ext != NULL) {
        vn.choice = SEQLOC_MIX;
        vn.next = NULL;
        vn.data.ptrvalue = bsp->seq_ext;
        slp = SeqLocFindNext (&vn, NULL);
        while (slp != NULL) {
          sip = SeqLocId (slp);
          if (sip != NULL) {
            bsp = BioseqFindCore (sip);
            if (bsp != NULL) {
              AddToSipList (ipp, bsp);
            } else {
              bsp = BioseqLockById (sip);
              if (bsp != NULL) {
                AddToSipList (ipp, bsp);
              }
              BioseqUnlockById (sip);
            }
          }
          slp = SeqLocFindNext (&vn, slp);
        }
      }
    }
  }
}

static Int4 SegmentedEntryList (SeqEntryPtr sep, Pointer mydata,
                                SeqEntryFunc mycallback,
                                Int4 index, Int2 indent)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  SeqLocPtr     slp;
  ValNode       vn;

  if (sep == NULL) return index;
  if (IS_Bioseq (sep)) {
    if (mycallback != NULL)
      (*mycallback) (sep, mydata, index, indent);
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp != NULL && bsp->repr == Seq_repr_seg && bsp->seq_ext != NULL) {
      vn.choice = SEQLOC_MIX;
      vn.next = NULL;
      vn.data.ptrvalue = bsp->seq_ext;
      slp = SeqLocFindNext (&vn, NULL);
      while (slp != NULL) {
        index++;
        slp = SeqLocFindNext (&vn, slp);
      }
    }
    return index + 1;
  }
  /*
  if (Bioseq_set_class (sep) == 4) return index;
  index++;
  */
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  sep = bssp->seq_set;
  indent++;
  while (sep != NULL) {
    index = SegmentedEntryList (sep, mydata, mycallback, index, indent);
    sep = sep->next;
  }
  return index;
}

#define SegmentedEntryCount( a )  SegmentedEntryList( a ,NULL,NULL,0,0);

static ENUM_ALIST(strand_alist)
{" ",             Seq_strand_unknown},  /* 0 */
{"Plus",          Seq_strand_plus},     /* 1 */
{"Minus",         Seq_strand_minus},    /* 2 */
{"Both",          Seq_strand_both},     /* 3 */
{"Reverse",       Seq_strand_both_rev}, /* 4 */
{"Other",         Seq_strand_other},    /* 255 */
END_ENUM_ALIST

static Boolean IsSeqIdInValNodeList (SeqIdPtr sip, ValNodePtr list)
{
  if (sip == NULL)
  {
    return FALSE;
  }
  while (list != NULL)
  {
    if (SeqIdComp (sip, list->data.ptrvalue) == SIC_YES)
    {
      return TRUE;
    }
    list = list->next;
  }
  return FALSE;
}

static Int4 FindLastStopForSeqId (SeqLocPtr slp, SeqIdPtr sip)
{
  Int4      last_stop = 0, new_stop;
  SeqLocPtr tmp_slp;
  SeqIdPtr  tmp_sip;
  
  tmp_slp = SeqLocFindNext (slp, NULL);
  while (tmp_slp != NULL)
  {
    tmp_sip = SeqLocId (tmp_slp);
    if (SeqIdComp (sip, tmp_sip) == SIC_YES)
    {
      new_stop = SeqLocStop (slp);
      if (new_stop > last_stop)
      {
        last_stop = new_stop;
      }
    }
    tmp_slp = SeqLocFindNext (slp, tmp_slp);
  }
  return last_stop;
}

extern void UpdateTagListPopupChoices (DialoG d, Int4 column);

/* We need to make sure that all IDs in the location are found in the enum list
 * for the interval editor ID Enum */
static void CorrectIntervalEditorSeqIdEnum (IntervalPagePtr ipp, SeqLocPtr slp)
{
  SeqLocPtr           tmp_slp;
  SeqIdPtr            sip;
  Int4                j;
  Boolean             found;
  ValNodePtr          missing_list = NULL, missing_vnp;
  Int4                new_count;
  SeqIdPtr PNTR       new_sip_list;
  EnumFieldAssoc PNTR new_alist;
  Int4 PNTR           new_lengths;
  BioseqPtr           bsp;
  Char                str [128];
  CharPtr             ptr;
  
  if (ipp == NULL || slp == NULL)
  {
    return;
  }
  
  tmp_slp = SeqLocFindNext (slp, NULL);
  while (tmp_slp != NULL)
  {
    sip = SeqLocId (tmp_slp);
    if (!IsSeqIdInValNodeList (sip, missing_list))
    {
      found = FALSE;
      for (j = 1; j <= ipp->count && !found; j++)
      {
        if (SeqIdComp (sip, ipp->sip_list [j]) == SIC_YES)
        {
          found = TRUE;
        }
      }
      /* this process takes longer, so don't combine it with the above
       * loop
       */
      if (!found)
      {
        for (j = 1; j <= ipp->count && !found; j++)
        {
          bsp = BioseqFindCore (ipp->sip_list [j]);
          if (bsp == NULL) {
            bsp = BioseqLockById (ipp->sip_list [j]);
            BioseqUnlockById (ipp->sip_list [j]);
          }
          if (bsp != NULL && SeqIdIn (sip, bsp->id))
          {
            found = TRUE;
          }
        }
      }
      if (!found)
      {
        ValNodeAddPointer (&missing_list, 0, sip);
      }
    }
    
    tmp_slp = SeqLocFindNext (slp, tmp_slp);
  }
  
  if (missing_list != NULL)
  {
    new_count = ipp->count + 4 + ValNodeLen (missing_list);
    new_sip_list = MemNew (sizeof (SeqIdPtr) * (size_t) new_count);
    new_alist = MemNew (sizeof (EnumFieldAssoc) * (size_t) new_count);
    new_lengths = MemNew (sizeof (Int4) * (size_t) new_count);
    /* first one is blank, remainder are actual data */
    for (j = 0; j < ipp->count + 1; j++)
    {
      new_sip_list [j] = ipp->sip_list [j];
      ipp->sip_list [j] = NULL;
      new_alist [j].name = ipp->alist [j].name;
      ipp->alist [j].name = NULL;
      new_alist [j].value = ipp->alist [j].value;
      new_lengths [j] = ipp->lengths [j];
    }
    
    missing_vnp = missing_list;
    while (j < new_count - 1 && missing_vnp != NULL)
    {
      new_sip_list [j] = SeqIdDup (missing_vnp->data.ptrvalue);
      new_alist [j].value = j;
      
      sip = new_sip_list [j];
      bsp = BioseqFindCore (sip);
      if (bsp == NULL) {
        bsp = BioseqLockById (sip);
        BioseqUnlockById (sip);
      }
      if (bsp != NULL)
      {
        new_lengths [j] = bsp->length;
        sip = SeqIdFindWorst (bsp->id);
      }
      else
      {
        new_lengths [j] = FindLastStopForSeqId (slp, sip);
      }
      SeqIdWrite (sip, str, PRINTID_REPORT, sizeof (str));
      ptr = StringChr (str, '|');
      if (ptr == NULL) {
        ptr = str;
      } else {
        ptr++;
      }
      new_alist [j].name = StringSave (ptr);
      missing_vnp = missing_vnp->next;
      j++;
    }
    /* set terminator for enum list */
    new_alist [j].name = NULL;
    new_alist [j].value = (UIEnum) 0;

    ipp->sip_list = MemFree (ipp->sip_list);
    ipp->sip_list = new_sip_list;
    ipp->alist = MemFree (ipp->alist);
    ipp->alist = new_alist;
    ipp->lengths = MemFree (ipp->lengths);
    ipp->lengths = new_lengths;
    ipp->count = j;
    
    ipp->alists [ipp->seqid_col] = ipp->alist;
    UpdateTagListPopupChoices (ipp->ivals, ipp->seqid_col);
  }
  missing_list = ValNodeFree (missing_list);
}

static void 
BuildIntervalString 
(IntervalPagePtr ipp,
 CharPtr         fuzz_from_ch,
 Int4            start,
 CharPtr         fuzz_to_ch,
 Int4            stop,
 Uint2           strand,
 Int4            seq,
 Int4            salp_num,
 Boolean         isInterval,
 Boolean         isPoint,
 CharPtr         buf,
 Int4            buf_len)
{
  CharPtr cp;
  Boolean need_tab = FALSE;
  
  if (ipp == NULL || fuzz_from_ch == NULL || fuzz_to_ch == NULL || buf == NULL || buf_len == 0)
  {
    return;
  }
  
  MemSet (buf, 0, buf_len);  
  
  if (isInterval)
  {
    sprintf (buf, "%s%ld\t%s%ld",
             fuzz_from_ch, (long) (start + 1),
             fuzz_to_ch, (long) (stop + 1));
    need_tab = TRUE;
  }
  else if (isPoint)
  {
    sprintf (buf, "%ld%s\t%ld%s",
             (long) (start + 1), fuzz_from_ch,
             (long) (stop + 1), fuzz_to_ch);
    need_tab = TRUE;
  }
  
  cp = buf + StringLen (buf);
  if (ipp->strand_col > -1)
  {
    if (need_tab)
    {
      StringCat (cp, "\t");
      cp++;
    }
    sprintf (cp, "%d", (int) strand);
    need_tab = TRUE;
    cp += StringLen (cp);
  }
  
  if (ipp->seqid_col > -1)
  {
    if (need_tab)
    {
      StringCat (cp, "\t");
      cp++;
    }
    sprintf (cp, "%d", (int) seq);
    need_tab = TRUE;
    cp += StringLen (cp);
  }
  
  if (ipp->aln_col > -1)
  {
    if (need_tab)
    {
      StringCat (cp, "\t");
      cp++;
    }
    sprintf (cp, "%d", (int) salp_num);
    need_tab = TRUE;
    cp += StringLen (cp);
  }
  StringCat (cp, "\n");
}

static void SeqLocPtrToIntervalPage (DialoG d, Pointer data)

{
  BioseqPtr        bsp;
  SeqLocPtr        firstSlp;
  Char             fuzz_from_ch [4];
  Char             fuzz_to_ch [4];
  ValNodePtr       head;
  SeqIdPtr         id;
  IntFuzzPtr       ifp;
  IntervalPagePtr  ipp;
  Boolean          isInterval;
  Boolean          isPoint;
  Int2             j;
  SeqLocPtr        lastSlp;
  SeqLocPtr        location;
  SeqLocPtr        next;
  SeqEntryPtr      oldscope;
  Boolean          partial5;
  Boolean          partial3;
  Int2             seq;
  SeqIntPtr        sip;
  SeqLocPtr        slp;
  SeqPntPtr        spp;
  Int4             start;
  Int4             stop;
  Char             str [255];
  Uint1            strand, aln_strand;
  TagListPtr       tlp;
  ValNodePtr       vnp;
  Int4             salp_num, salp_row;
  SeqIdPtr         tmp_sip, bsp_id;

  ipp = (IntervalPagePtr) GetObjectExtra (d);
  if (ipp == NULL) return;
  SafeSetStatus (ipp->partial5, FALSE);
  SafeSetStatus (ipp->partial3, FALSE);
  tlp = GetObjectExtra (ipp->ivals);
  if (tlp == NULL) return;

  location = (SeqLocPtr) data;
  partial5 = FALSE;
  partial3 = FALSE;
  head = NULL;

  if (location == NULL)
  {
  	if (ipp->count == 1)
  	{
       sprintf (str, "\t\t\t1\n");
       vnp = ValNodeNew (head);
       if (head == NULL) {
         head = vnp;
       }
       if (vnp != NULL) {
         vnp->data.ptrvalue = StringSave (str);
       }	
  	}
  }
  if (location != NULL) {
    CorrectIntervalEditorSeqIdEnum (ipp, location);
    firstSlp = NULL;
    lastSlp = NULL;
    slp = SeqLocFindNext (location, NULL);
    while (slp != NULL) {
      if (firstSlp == NULL) {
        firstSlp = slp;
      }
      lastSlp = slp;
      next = SeqLocFindNext (location, slp);
      if (slp->choice == SEQLOC_NULL) {
        SafeSetStatus (ipp->nullsBetween, TRUE);
        SafeShow (ipp->nullsBetween);
      } else {
        id = SeqLocId (slp);
        if (id != NULL) {
          bsp = BioseqFind (id);
          if (bsp == NULL) {
            oldscope = SeqEntrySetScope (NULL);
            if (oldscope != NULL) {
              bsp = BioseqFind (id);
              SeqEntrySetScope (oldscope);
            }
          }
          isInterval = TRUE;
          isPoint = FALSE;
          StringCpy (fuzz_from_ch, "");
          StringCpy (fuzz_to_ch, "");
          if (bsp == NULL) {
            start = SeqLocStart (slp);
            stop = SeqLocStop (slp);
          }
          else
          {
            start = GetOffsetInBioseq (slp, bsp, SEQLOC_START);
            stop = GetOffsetInBioseq (slp, bsp, SEQLOC_STOP);
          }
          if (start == stop && slp->choice == SEQLOC_PNT) {
            spp = (SeqPntPtr) slp->data.ptrvalue;
            if (spp != NULL) {
              ifp = spp->fuzz;
              if (ifp != NULL && ifp->choice == 4 && ifp->a ==  3) {
                isInterval = FALSE;
                isPoint = TRUE;
                StringCpy (fuzz_from_ch, "^");
                /* start--; */  /* compensate for other fix */
                stop++; /* compensate for other fix */
              }
            }
          }
          strand = SeqLocStrand (slp);
          if (strand > Seq_strand_both_rev && strand != Seq_strand_other) {
            strand = Seq_strand_unknown;
          }
          /*
          if (strand == Seq_strand_unknown) {
            strand = Seq_strand_plus;
          }
          */
          if (! ipp->nucsOK) {
            strand = 0;
          }
          seq = 0;
          if (ipp->sip_list != NULL && bsp != NULL) {
            for (j = 1; j <= ipp->count && seq == 0; j++) {
              if (SeqIdComp (SeqIdFindBest (bsp->id, 0), ipp->sip_list[j]) == SIC_YES) {
                seq = j;
              }
            }
          }
          
          salp_num = 0;
          salp_row = 0;
          if (seq > 0 && ipp->salp_list != NULL)
          {
            salp_num = 1;
            while (salp_num <= ipp->aln_count
                   && salp_row == 0)
            {
              tmp_sip = SeqIdPtrFromSeqAlign (ipp->salp_list [salp_num]);
              bsp = BioseqFind (ipp->sip_list [seq]);
              if (bsp != NULL)
              {
                for (bsp_id = bsp->id;
                     bsp_id != NULL && salp_row == 0; 
                     bsp_id = bsp_id->next)
                {
                  salp_row = SeqIdOrderInBioseqIdList(bsp_id, tmp_sip);
                }
              }
              if (salp_row < 1)
              {
                salp_num++;
              }
            }
            if (salp_row < 1)
            {
              salp_num = 0;
            }
            else
            {
              start = AlnMgr2MapBioseqToSeqAlign(ipp->salp_list [salp_num], start, salp_row);
              stop = AlnMgr2MapBioseqToSeqAlign(ipp->salp_list [salp_num], stop, salp_row);
              aln_strand = AlnMgr2GetNthStrand (ipp->salp_list [salp_num], salp_row);
              /* if reverse strand in alignment, reverse strand for location */
              if (aln_strand == Seq_strand_minus)
              {
                if (strand == Seq_strand_minus)
                {
                  strand = Seq_strand_plus;
                }
                else if (strand == Seq_strand_plus)
                {
                  strand = Seq_strand_minus;
                }
              }
            }
          }
          
          BuildIntervalString (ipp, fuzz_from_ch, start, fuzz_to_ch, stop,
                               strand, seq, salp_num, isInterval, isPoint,
                               str, sizeof (str));
          vnp = ValNodeNew (head);
          if (head == NULL) {
            head = vnp;
          }
          if (vnp != NULL) {
            vnp->data.ptrvalue = StringSave (str);
          }
        }
      }
      slp = next;
    }
    if (firstSlp != NULL) {
      if (firstSlp->choice == SEQLOC_INT && firstSlp->data.ptrvalue != NULL) {
        sip = (SeqIntPtr) firstSlp->data.ptrvalue;
        if (sip->strand == Seq_strand_minus || sip->strand == Seq_strand_both_rev) {
          ifp = sip->if_to;
          if (ifp != NULL && ifp->choice == 4 && ifp->a == 1) {
            partial5 = TRUE;
          }
        } else {
          ifp = sip->if_from;
          if (ifp != NULL && ifp->choice == 4 && ifp->a == 2) {
            partial5 = TRUE;
          }
        }
      } else if (firstSlp->choice == SEQLOC_PNT && firstSlp->data.ptrvalue != NULL) {
        spp = (SeqPntPtr) firstSlp->data.ptrvalue;
        if (spp->strand == Seq_strand_minus || spp->strand == Seq_strand_both_rev) {
          ifp = spp->fuzz;
          if (ifp != NULL && ifp->choice == 4 && ifp->a == 1) {
            partial5 = TRUE;
          }
        } else {
          ifp = spp->fuzz;
          if (ifp != NULL && ifp->choice == 4 && ifp->a == 2) {
            partial5 = TRUE;
          }
        }
      }
    }
    if (lastSlp != NULL) {
      if (lastSlp->choice == SEQLOC_INT && lastSlp->data.ptrvalue != NULL) {
        sip = (SeqIntPtr) lastSlp->data.ptrvalue;
        if (sip->strand == Seq_strand_minus || sip->strand == Seq_strand_both_rev) {
          ifp = sip->if_from;
          if (ifp != NULL && ifp->choice == 4 && ifp->a == 2) {
            partial3 = TRUE;
          }
        } else {
          ifp = sip->if_to;
          if (ifp != NULL && ifp->choice == 4 && ifp->a == 1) {
            partial3 = TRUE;
          }
        }
      } else if (lastSlp->choice == SEQLOC_PNT && lastSlp->data.ptrvalue != NULL) {
        spp = (SeqPntPtr) lastSlp->data.ptrvalue;
        if (spp->strand == Seq_strand_minus || spp->strand == Seq_strand_both_rev) {
          ifp = spp->fuzz;
          if (ifp != NULL && ifp->choice == 4 && ifp->a == 2) {
            partial3 = TRUE;
          }
        } else {
          ifp = spp->fuzz;
          if (ifp != NULL && ifp->choice == 4 && ifp->a == 1) {
            partial3 = TRUE;
          }
        }
      }
    }
  }
  SafeSetStatus (ipp->partial5, partial5);
  SafeSetStatus (ipp->partial3, partial3);

  SendMessageToDialog (tlp->dialog, VIB_MSG_RESET);
  tlp->vnp = head;
  SendMessageToDialog (tlp->dialog, VIB_MSG_REDRAW);
  for (j = 0, vnp = tlp->vnp; vnp != NULL; j++, vnp = vnp->next) {
  }
  tlp->max = MAX ((Int2) 0, (Int2) (j - tlp->rows + 1));
  CorrectBarMax (tlp->bar, tlp->max);
  CorrectBarPage (tlp->bar, tlp->rows - 1, tlp->rows - 1);
}

extern void SetSequenceAndStrandForIntervalPage (DialoG d)
{
  SendMessageToDialog (d, NUM_VIB_MSG + 1);
}

static Boolean 
ReadFromValueFromDialogLine 
(CharPtr line,
 Int4Ptr p_from_val,
 BoolPtr fuzz_after)
{
  CharPtr txt, ptr;
  Int4    len;
  Boolean okay;
  Int4    val;
  Int4    j;
  Char    ch;
  
  if (p_from_val == NULL || fuzz_after == NULL || StringHasNoText (line))
  {
    return FALSE;
  }

  txt = ExtractTagListColumn (line, 0);
  okay = FALSE;
  len = StringLen (txt);
  for (j = 0; j < len; j++) {
    ch = txt [j];
    if (ch != ' ' && ch != '\t' && ch != '\n') {
      okay = TRUE;
    }
  }
  if (okay) {
    *fuzz_after = FALSE;
    *p_from_val = 0;

    ptr = StringChr (txt, '<');
    if (ptr != NULL) {
      *ptr = ' ';
    }
    ptr = StringChr (txt, '>');
    if (ptr != NULL) {
      *ptr = ' ';
    }
    ptr = StringChr (txt, '^');
    if (ptr != NULL) {
      *fuzz_after = TRUE;
      *ptr = ' ';
    }
    TrimSpacesAroundString (txt);
    if (StrToLong (txt, &val)) {
      *p_from_val = val;
    } else {
      okay = FALSE;
    }
  } else {
    okay = FALSE;
  }
  MemFree (txt);
  return okay;
}

static Boolean 
ReadToValueFromDialogLine 
(CharPtr line,
 Int4Ptr p_to_val,
 Int4Ptr p_from_val,
 BoolPtr fuzz_after,
 BoolPtr isInterval,
 BoolPtr isPoint,
 Boolean partial5,
 Boolean partial3)
{
  CharPtr txt, ptr;
  Int4    val, tmp;
  Boolean okay = TRUE;
  
  if (p_to_val == NULL || p_from_val == NULL || fuzz_after == NULL
      || isInterval == NULL || isPoint == NULL)
  {
    return FALSE;
  }
  
  *p_to_val = 0;
  
  txt = ExtractTagListColumn (line, 1);
  if (! StringHasNoText (txt)) {
    ptr = StringChr (txt, '<');
    if (ptr != NULL) {
      *ptr = ' ';
    }
    ptr = StringChr (txt, '>');
    if (ptr != NULL) {
      *ptr = ' ';
    }
    ptr = StringChr (txt, '^');
    if (ptr != NULL) {
      *fuzz_after = TRUE;
      *ptr = ' ';
    }
    TrimSpacesAroundString (txt);
    if (StrToLong (txt, &val)) {
      *p_to_val = val;
      if (*fuzz_after && *p_to_val == *p_from_val + 1) {
        *isInterval = FALSE;
        *isPoint = TRUE;
        /* from++; */ /* this was causing point to be thrown off */
      } else if (*p_to_val == *p_from_val && (! partial5) && (! partial3)) {
        *isInterval = FALSE;
        *isPoint = TRUE;
      }
    } else {
      okay = FALSE;
    }
  } else {
    /*
    okay = FALSE;
    */
    *isInterval = FALSE;
    *isPoint = TRUE;
    *p_to_val = *p_from_val;
  }
  MemFree (txt);
  
  if (okay && *isInterval) {
    if (*p_from_val > *p_to_val) {
      tmp = *p_from_val;
      *p_from_val = *p_to_val;
      *p_to_val = tmp;
    }
  }
  

  return okay;
}

static Boolean 
ReadStrandFromDialogLine 
(CharPtr  line,
 Int4     strand_col,
 Uint2Ptr p_strand_val,
 Uint2Ptr p_prev_strand_val)
{
  CharPtr txt;
  Int2    val2;
  
  if (p_strand_val == NULL
      || p_prev_strand_val == NULL)
  {
    return FALSE;
  }
  
  *p_strand_val = Seq_strand_unknown;
  if (strand_col > -1) {
    txt = ExtractTagListColumn (line, strand_col);
    if (txt != NULL && StrToInt (txt, &val2)) {
      *p_strand_val = val2;
      if (*p_strand_val > Seq_strand_both_rev) {
        *p_strand_val = Seq_strand_other;
      }
      *p_prev_strand_val = *p_strand_val;
    } else {
      *p_strand_val = *p_prev_strand_val;
    }
    MemFree (txt);
  }
  if (*p_strand_val == Seq_strand_unknown) {
    *p_strand_val = Seq_strand_plus;
  }
  return TRUE;
}

static Boolean 
ReadSeqIdFromDialogLine 
(CharPtr       line, 
 Int4          seqid_col,
 SeqIdPtr PNTR p_sip_val,
 SeqIdPtr PNTR p_prev_sip_val,
 SeqIdPtr PNTR sip_list,
 Int4          sip_count)
{
  CharPtr txt;
  Int2    val2;
  Boolean okay = TRUE;
  
  if (p_sip_val == NULL || p_prev_sip_val == NULL
      || sip_count == 0 || sip_list == NULL
      || seqid_col < 0)
  {
    return FALSE;
  }
  *p_sip_val = NULL;
  txt = ExtractTagListColumn (line, seqid_col);
  if (txt != NULL) {
    if (! StrToInt (txt, &val2) || val2 <= 0)
    {
      if (*p_prev_sip_val != NULL)
      {
        *p_sip_val = *p_prev_sip_val;
      }
      else
      {
        okay = FALSE;
      }
    }
    else if (val2 <= sip_count)
    {
      *p_sip_val = sip_list [val2];
      *p_prev_sip_val = *p_sip_val;
    } else {
      okay = FALSE;
    }
  }
  else
  {
    okay = FALSE;
  }
  MemFree (txt);
  return okay;
}

static Boolean 
ReadSeqAlignFromDialogLine
(CharPtr          line,
 Int4             aln_col,
 SeqAlignPtr PNTR p_salp_val,
 SeqAlignPtr PNTR p_prev_salp_val,
 SeqAlignPtr PNTR salp_list,
 Int4             salp_count)
{
  CharPtr txt;
  Int2    val2;
  Boolean okay = TRUE;
  
  if (p_salp_val == NULL || p_prev_salp_val == NULL 
      || salp_count == 0 || salp_list == NULL
      || aln_col < 0)
  {
    return FALSE;
  }
  
  txt = ExtractTagListColumn (line, aln_col);
  if (txt != NULL) {
    if (! StrToInt (txt, &val2) || val2 <= 0)
    {
      if (*p_prev_salp_val != NULL)
      {
        *p_salp_val = *p_prev_salp_val;
      }
      else
      {
        okay = FALSE;
      }
    }
    else if (val2 <= salp_count)
    {
      *p_salp_val = salp_list [val2];
      *p_prev_salp_val = *p_salp_val;
    } else {
      okay = FALSE;
    }
  }
  else
  {
    if (*p_prev_salp_val == NULL)
    {
      okay = FALSE;
    }
    else
    {
      *p_salp_val = *p_prev_salp_val;
    }
 }
 MemFree (txt);
 return okay;
}

static Boolean 
GetBioseqAlignmentRow 
(SeqAlignPtr salp,
 BioseqPtr   bsp, 
 Int4Ptr     aln_row)
{
  SeqIdPtr  tmp_sip, bsp_id;
  if (salp == NULL || bsp == NULL || aln_row == NULL)
  {
    return FALSE;
  }
  
  *aln_row = 0;

  tmp_sip = SeqIdPtrFromSeqAlign (salp);
  
  for (bsp_id = bsp->id;
       bsp_id != NULL && *aln_row == 0; 
       bsp_id = bsp_id->next)
  {
    *aln_row = SeqIdOrderInBioseqIdList(bsp_id, tmp_sip);
  }
  if (*aln_row < 1)
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}

static void 
ListAlignmentsThatContainSequence 
(SeqAlignPtr PNTR salp_list, 
 Int4        salp_count,
 BioseqPtr   bsp,
 SeqIdPtr    sip,
 ValNodePtr  PNTR head)
{
  Int4       i;
  Char       str [34];
  Char       id_label [128];
  Int4       aln_row;
  ValNodePtr good_aln = NULL;
  CharPtr    none_found_fmt = "%s is not found in any alignments";
  CharPtr    one_found_fmt = "%s is found in %s";
  CharPtr    some_found_fmt = "%s is found in the following alignments";
  CharPtr    err_msg;
  CharPtr    cp;
  
  if (head == NULL || bsp == NULL || sip == NULL)
  {
    return;
  }
  SeqIdWrite (sip, id_label, PRINTID_REPORT, sizeof (id_label));
  
  /* indent the names of the alignments */
  str [0] = ' ';
  str [1] = ' ';
  str [2] = ' ';
  str [3] = ' ';
  for (i = 0; i < salp_count; i++)
  {
    if (GetBioseqAlignmentRow (salp_list [i], bsp, &aln_row))
    {
      CreateSeqAlignLabel (salp_list [i], str + 4, sizeof (str) - 4);
      ValNodeAddPointer (&good_aln, 0, StringSave (str));
    }           
  }
  if (good_aln == NULL)
  {
    err_msg = (CharPtr) MemNew ((StringLen (none_found_fmt) + StringLen (id_label))
                                * sizeof (Char));
    if (err_msg != NULL)
    {
      sprintf (err_msg, none_found_fmt, id_label);
      ValNodeAddPointer (head, 0, err_msg);
    }
  }
  else if (good_aln->next == NULL)
  {
    err_msg = (CharPtr) MemNew ((StringLen (one_found_fmt) 
                                 + StringLen (id_label)
                                 + StringLen (good_aln->data.ptrvalue))
                                * sizeof (Char));
    if (err_msg != NULL)
    {
      cp = (CharPtr) good_aln->data.ptrvalue;
      /* skip over indented space */
      cp += 4;
      sprintf (err_msg, one_found_fmt, id_label, cp);
      ValNodeAddPointer (head, 0, err_msg);
      good_aln = ValNodeFreeData (good_aln);
    }
  }
  else
  {
    err_msg = (CharPtr) MemNew ((StringLen (some_found_fmt) 
                                 + StringLen (id_label))
                                * sizeof (Char));
    if (err_msg != NULL)
    {
      sprintf (err_msg, some_found_fmt, id_label);
      ValNodeAddPointer (head, 0, err_msg);
      ValNodeLink (head, good_aln);
    }
  }
}

static void 
ListSequencesInAlignment 
(SeqIdPtr    PNTR sip_list,
 Int4        sip_count,
 SeqAlignPtr salp,
 ValNodePtr  PNTR head)
{
  SeqIdPtr   tmp_sip, aln_id;
  BioseqPtr  bsp;
  Int4       i;
  Boolean    found;
  Char       id_label [128];
  
  if (salp == NULL || head == NULL)
  {
    return;
  }

  tmp_sip = SeqIdPtrFromSeqAlign (salp);
  if (tmp_sip == NULL)
  {
    ValNodeAddPointer (head, 0, StringSave ("Selected alignment contains no sequences"));
    return;
  }
  /* indent list of sequence IDs */
  id_label [0] = ' ';
  id_label [1] = ' ';
  id_label [2] = ' ';
  id_label [3] = ' ';
  ValNodeAddPointer (head, 0, StringSave ("The selected alignment contains the following sequences:"));
  for (aln_id = tmp_sip; aln_id != NULL; aln_id = aln_id->next)
  {
    bsp = BioseqFind (aln_id);
    if (bsp == NULL)
    {
      continue;
    }
    found = FALSE;
    for (i = 0; i < sip_count && ! found; i++)
    {
      if (SeqIdIn (sip_list [i], bsp->id))
      {
        SeqIdWrite (sip_list [i], id_label + 4, PRINTID_REPORT, sizeof (id_label) - 4);
        ValNodeAddPointer (head, 0, StringSave (id_label));
        found = TRUE;
      }
    }
  }
}

extern Boolean AdjustFromForGap (Int4Ptr p_from, SeqAlignPtr salp, Int4 aln_len, Int4 aln_row)
{
  Int4 aln_from, aln_offset = 0;
  
  if (p_from == NULL || salp == NULL || aln_len == 0 || aln_row == 0
      || *p_from < 1)
  {
    return FALSE;
  }
  
  aln_from = AlnMgr2MapSeqAlignToBioseq(salp, (*p_from) - 1, aln_row);
  
  while (aln_from == -2 && (*p_from) + aln_offset < aln_len)
  {
    aln_offset ++;
    aln_from = AlnMgr2MapSeqAlignToBioseq(salp, (*p_from) + aln_offset - 1, aln_row);
  }
  if (aln_from < 0)
  {
    return FALSE;
  }
  else
  {
    *p_from = aln_from + 1;
    return TRUE;
  }
}

extern Boolean AdjustToForGap (Int4Ptr p_to, SeqAlignPtr salp, Int4 aln_row)
{
  Int4 aln_to, aln_offset = 0;
  
  if (p_to == NULL || salp == NULL || aln_row == 0
      || *p_to < 1)
  {
    return FALSE;
  }
  
  aln_to = AlnMgr2MapSeqAlignToBioseq(salp, (*p_to) - 1, aln_row);
  aln_offset = 0;
  while (aln_to == -2 && (*p_to) - 1 - aln_offset >= 0)
  {
    aln_offset ++;
    aln_to = AlnMgr2MapSeqAlignToBioseq(salp, (*p_to) - 1 - aln_offset, aln_row);
  }
  if (aln_to < 0)
  {
    return FALSE;
  }
  else
  {
    *p_to = aln_to + 1;
    return TRUE;
  }
}

static Boolean AdjustStrandForAlignment (Uint2Ptr p_strand, SeqAlignPtr salp, Int4 aln_row)
{
  Uint2 aln_strand;
  
  if (p_strand == NULL || salp == NULL || aln_row < 1)
  {
    return FALSE;
  }
  
  aln_strand = AlnMgr2GetNthStrand (salp, aln_row);
  if (aln_strand == Seq_strand_minus)
  {
    /* if alignment strand is minus, reverse strand of location */
    if (*p_strand == Seq_strand_minus)
    {
      *p_strand = Seq_strand_plus;
    }
    else if (*p_strand == Seq_strand_plus)
    {
      *p_strand = Seq_strand_minus;
    }
  }
  return TRUE;
}

static Boolean CoordinatesValidForBioseq (Int4 from, Int4 to, BioseqPtr bsp)
{
  if (bsp == NULL || from < 1 || to > bsp->length) 
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}

static ValNodePtr TestIntervalEditor (DialoG d)
{
  IntervalPagePtr  ipp;
  TagListPtr       tlp;
  ValNodePtr       head = NULL, vnp;
  Boolean          from_ok, to_ok, seqid_ok, salp_ok;
  Boolean          partial5, partial3;
  Int4             to, from, aln_row;
  Int4             interval_num;
  Uint2            strand, prev_strand = Seq_strand_unknown;
  SeqIdPtr         seqid, prev_seqid = NULL;
  SeqAlignPtr      salp, prev_salp = NULL;
  Char             err_msg[200];
  BioseqPtr        bsp = NULL;
  Int4             aln_len;
  Boolean          isInterval, isPoint, fuzz_before, fuzz_after;
  Int2             fuzz_from;
  Int2             fuzz_to;
  
  ipp = (IntervalPagePtr) GetObjectExtra (d);
  if (ipp == NULL) 
  {
    ValNodeAddPointer (&head, 0, StringSave ("No dialog data"));
    return head;
  }
  
  tlp = GetObjectExtra (ipp->ivals);
  if (tlp == NULL)
  {
    ValNodeAddPointer (&head, 0, StringSave ("No dialog data"));
    return head;
  }
  
  if (tlp->vnp == NULL)
  {
    ValNodeAddPointer (&head, 0, StringSave ("No location intervals listed!"));
  }

  partial5 = GetStatus (ipp->partial5);
  partial3 = GetStatus (ipp->partial3);

  for (vnp = tlp->vnp, interval_num = 1;
       vnp != NULL; 
       vnp = vnp->next, interval_num++) {
    isInterval = TRUE;
    isPoint = FALSE;
    fuzz_from = -1;
    fuzz_to = -1;
    fuzz_before = FALSE;
    fuzz_after = FALSE;
    from = 0;
    to = 0;

    from_ok = ReadFromValueFromDialogLine ((CharPtr) vnp->data.ptrvalue,
                                           &from, &fuzz_after);
    if (!from_ok)
    {
      /* we'll silently ignore lines that have no from value */
      continue;
    }
    
    to_ok = ReadToValueFromDialogLine ((CharPtr) vnp->data.ptrvalue, 
                                         &to, &from, &fuzz_after, &isInterval, &isPoint,
                                         partial5, partial3);
    if (!to_ok)
    {
      sprintf (err_msg, "Bad to value in interval %d", interval_num);
      ValNodeAddPointer (&head, 0, StringSave ("err_msg"));
    }
    strand = Seq_strand_unknown;
    ReadStrandFromDialogLine (vnp->data.ptrvalue, ipp->strand_col, 
                              &strand, &prev_strand);
    
    seqid_ok = ReadSeqIdFromDialogLine (vnp->data.ptrvalue, ipp->seqid_col,
                                         &seqid, &prev_seqid,
                                         ipp->sip_list, ipp->count);
    if (seqid_ok)
    {
      bsp = BioseqFind (seqid);
      if (bsp == NULL)
      {
        sprintf (err_msg, "Can't find bioseq for interval %d", interval_num);
        ValNodeAddPointer (&head, 0, StringSave (err_msg));
        seqid_ok = FALSE;
      }
    }
    else
    {
      sprintf (err_msg, "No sequence ID in interval %d", interval_num);
      ValNodeAddPointer (&head, 0, StringSave (err_msg));
    }

      
    /* get SeqAlign */
    salp = NULL;
    if (ipp->salp_list == NULL)
    {
      salp_ok = TRUE;
      if (bsp != NULL)
      {
        if (!CoordinatesValidForBioseq (from, to, bsp))
        {
          sprintf (err_msg, "Coordinates for interval %d are not in sequence (1-%d)",
                   interval_num, bsp->length);
          ValNodeAddPointer (&head, 0, StringSave (err_msg));
        }
      }
    }
    else
    {
      salp_ok = ReadSeqAlignFromDialogLine (vnp->data.ptrvalue, 
                                            ipp->aln_col,
                                            &salp, &prev_salp,
                                            ipp->salp_list,
                                            ipp->aln_count);
      if (!salp_ok)
      {
        sprintf (err_msg, "No alignment for interval %d", interval_num);
        ValNodeAddPointer (&head, 0, StringSave (err_msg));
      }
      
      if (salp_ok && seqid_ok)
      {
        aln_row = 0;
        if (GetBioseqAlignmentRow (salp, bsp, &aln_row))
        {
          aln_len = SeqAlignLength (salp);
          if (from < 1 || to > aln_len)
          {
            sprintf (err_msg, "Coordinates for interval %d are not in alignment interval (%d-%d)",
                     interval_num, 1, aln_len);
            ValNodeAddPointer (&head, 0, StringSave (err_msg));
          }
          else
          {
            /* check for locations in gaps */
            if (!AdjustFromForGap (&from, salp, aln_len, aln_row)
                || ! AdjustToForGap (&to, salp, aln_row))
            {
              sprintf (err_msg, "Interval %d is completely contained in a gap in the alignment.",
                       interval_num);
              ValNodeAddPointer (&head, 0, StringSave (err_msg));
            }
          }
        }
        else
        {
          sprintf (err_msg, "Sequence for interval %d not in selected alignment", interval_num);
          ValNodeAddPointer (&head, 0, StringSave (err_msg));
          ListAlignmentsThatContainSequence (ipp->salp_list, ipp->aln_count, 
                                             bsp, seqid, &head);
          ListSequencesInAlignment (ipp->sip_list, ipp->count, salp, &head);
        }
      }
    }
  }
  
  return head;
}


static Boolean 
AddLocToList 
(SeqIdPtr       seqid,
 Uint2          strand,
 Int4           from,
 Int4           to,
 Boolean        add_null,
 Int2           fuzz_from,
 Int2           fuzz_to,
 Boolean        fuzz_before,
 Boolean        fuzz_after,
 Boolean        partial5,
 Boolean        partial3,
 SeqLocPtr PNTR pslp)
{
  SeqLocPtr slp;
  SeqLocPtr tmploc1, tmploc2;
  Boolean   isInterval;
  Boolean   isPoint;
  Int4      tmp;
  
  if (pslp == NULL || seqid == NULL)
  {
    return FALSE;
  }
            
  if (add_null) {
    /* add NULL location between last location and this one */
    slp = ValNodeNew (NULL);
    if (slp != NULL) {
      slp->choice = SEQLOC_NULL;
      tmploc1 = *pslp;
      if (tmploc1 != NULL) {
        if (tmploc1->choice == SEQLOC_MIX) {
          tmploc2 = (ValNodePtr) (tmploc1->data.ptrvalue);
          if (tmploc2 != NULL) {
            while (tmploc2->next != NULL) {
              tmploc2 = tmploc2->next;
            }
            tmploc2->next = slp;
          }
        } else {
          tmploc2 = ValNodeNew (NULL);
          if (tmploc2 != NULL) {
            tmploc2->choice = SEQLOC_MIX;
            tmploc2->data.ptrvalue = (Pointer) tmploc1;
            tmploc1->next = slp;
            *pslp = tmploc2;
          }
        }
      }
    }
  }
    
  isInterval = TRUE;
  isPoint = FALSE;
  
  /* make sure from and to are in correct order */
  if (from > to)
  {
    tmp = from;
    from = to;
    to = tmp;
  }

  if (fuzz_after && to == from + 1) {
    isInterval = FALSE;
    isPoint = TRUE;
    /* from++; */ /* this was causing point to be thrown off */
  } else if (to == from && (! partial5) && (! partial3)) {
    isInterval = FALSE;
    isPoint = TRUE;
  }
          
  if (isInterval) {
    AddIntToSeqLoc (pslp, from - 1, to - 1, seqid,
                    fuzz_from, fuzz_to, strand);
  } else if (isPoint) {
    AddSeqLocPoint (pslp, seqid, from, fuzz_before, fuzz_after, strand);
  }
  return TRUE;
}

static void SetPartialsForOneLocation (SeqLocPtr master_slp, Boolean partial5, Boolean partial3)
{
  SeqLocPtr  firstSlp, lastSlp, tmp_slp;
  IntFuzzPtr ifp;
  SeqIntPtr  sip;
  SeqPntPtr  spp;
  
  if (master_slp == NULL)
  {
    return;
  }
  
  /* now set partials for location */
  firstSlp = NULL;
  lastSlp = NULL;
  tmp_slp = SeqLocFindNext (master_slp, NULL);
  while (tmp_slp != NULL) {
    if (firstSlp == NULL) {
      firstSlp = tmp_slp;
    }
    lastSlp = tmp_slp;
    tmp_slp = SeqLocFindNext (master_slp, tmp_slp);
  }
  if (firstSlp != NULL && partial5) {
    if (firstSlp->choice == SEQLOC_INT && firstSlp->data.ptrvalue != NULL) {
      sip = (SeqIntPtr) firstSlp->data.ptrvalue;
      ifp = IntFuzzNew ();
      if (ifp != NULL) {
        ifp->choice = 4;
        if (sip->strand == Seq_strand_minus || sip->strand == Seq_strand_both_rev) {
          sip->if_to = ifp;
          ifp->a = 1;
        } else {
          sip->if_from = ifp;
          ifp->a = 2;
        }
      }
    } else if (firstSlp->choice == SEQLOC_PNT && firstSlp->data.ptrvalue != NULL) {
      spp = (SeqPntPtr) firstSlp->data.ptrvalue;
      ifp = IntFuzzNew ();
      if (ifp != NULL) {
        ifp->choice = 4;
        if (spp->strand == Seq_strand_minus || spp->strand == Seq_strand_both_rev) {
          spp->fuzz = ifp;
          ifp->a = 1;
        } else {
          spp->fuzz = ifp;
          ifp->a = 2;
        }
      }
    }
  }
  if (lastSlp != NULL && partial3) {
    if (lastSlp->choice == SEQLOC_INT && lastSlp->data.ptrvalue != NULL) {
      sip = (SeqIntPtr) lastSlp->data.ptrvalue;
      ifp = IntFuzzNew ();
      if (ifp != NULL) {
        ifp->choice = 4;
        if (sip->strand == Seq_strand_minus || sip->strand == Seq_strand_both_rev) {
          sip->if_from = ifp;
          ifp->a = 2;
        } else {
          sip->if_to = ifp;
          ifp->a = 1;
        }
      }
    } else if (lastSlp->choice == SEQLOC_PNT && lastSlp->data.ptrvalue != NULL) {
      spp = (SeqPntPtr) lastSlp->data.ptrvalue;
      ifp = IntFuzzNew ();
      if (ifp != NULL) {
        ifp->choice = 4;
        if (spp->strand == Seq_strand_minus || spp->strand == Seq_strand_both_rev) {
          spp->fuzz = ifp;
          ifp->a = 2;
        } else {
          spp->fuzz = ifp;
          ifp->a = 1;
        }
      }
    }
  }
}

static SeqLocPtr 
ReadSingleSeqLoc 
(TagListPtr      tlp,
 IntervalPagePtr ipp,
 Boolean         partial5,
 Boolean         partial3, 
 Boolean         nullsBetween)
{
  Int4             from, to, aln_row, aln_len;
  Boolean          fuzz_after;
  Boolean          fuzz_before;
  Int2             fuzz_from;
  Int2             fuzz_to;
  Boolean          isInterval;
  Boolean          isPoint;
  Boolean          notFirst;
  Boolean          okay;
  SeqIdPtr         seqid, prev_sip;
  SeqLocPtr        master_slp;
  Uint2            strand, prev_strand;
  ValNodePtr       vnp;
  SeqAlignPtr      salp, prev_salp;
  BioseqPtr        bsp;

  if (tlp == NULL)
  {
    return NULL;
  }
  
  prev_sip = NULL;
  prev_salp = NULL;
  prev_strand = Seq_strand_unknown;
  master_slp = NULL;

  notFirst = FALSE;
  for (vnp = tlp->vnp; vnp != NULL; vnp = vnp->next) {
    if (StringHasNoText (vnp->data.ptrvalue))
    {
      continue;
    }
    fuzz_from = -1;
    fuzz_to = -1;
    fuzz_before = FALSE;
    fuzz_after = FALSE;
    from = 0;
    to = 0;
    isInterval = TRUE;
    isPoint = FALSE;
    strand = Seq_strand_unknown;
    seqid = NULL;
    salp = NULL;
    aln_row = 0;
    if (!ReadFromValueFromDialogLine ((CharPtr) vnp->data.ptrvalue,
                                           &from, &fuzz_after))
    {
      continue;
    }
    okay = ReadToValueFromDialogLine ((CharPtr) vnp->data.ptrvalue,
                                         &to, &from, &fuzz_after,
                                         &isInterval, &isPoint, 
                                         partial5, partial3)
           && ReadStrandFromDialogLine ((CharPtr) vnp->data.ptrvalue,
                                        ipp->strand_col,
                                        &strand, &prev_strand)
           && ReadSeqIdFromDialogLine ((CharPtr) vnp->data.ptrvalue, ipp->seqid_col,
                                       &seqid, &prev_sip,
                                       ipp->sip_list, ipp->count)
           && (ipp->salp_list == NULL 
               || ReadSeqAlignFromDialogLine ((CharPtr) vnp->data.ptrvalue, ipp->aln_col,
                                              &salp, &prev_salp, 
                                              ipp->salp_list, ipp->aln_count))
           && ((bsp = BioseqFind (seqid)) != NULL)
           && (salp != NULL || CoordinatesValidForBioseq (from, to, bsp))
           && (salp == NULL || 
                (GetBioseqAlignmentRow (salp, bsp, &aln_row)
                 && AdjustStrandForAlignment (&strand, salp, aln_row)
                 && (aln_len = SeqAlignLength (salp)) > 0
                 && AdjustFromForGap (&from, salp, aln_len, aln_row)
                 && AdjustToForGap (&to, salp, aln_row)))
           && AddLocToList (seqid, strand, from, to,
                            (nullsBetween  && notFirst),
                            fuzz_from, fuzz_to,
                            fuzz_before, fuzz_after,
                            partial5, partial3,
                            &(master_slp));
    if (okay)
    {
      notFirst = TRUE;
    }
    else
    {
      master_slp = SeqLocFree (master_slp);
      return NULL;
    }
  }
  
  /* now set partials for location */
  SetPartialsForOneLocation (master_slp, partial5, partial3);  
  return master_slp;
}

static SeqLocPtr PNTR FreeSeqLocArray (SeqLocPtr PNTR loc_list, Int4 num_loc)
{
  Int4 j;
  if (loc_list != NULL)
  {
    for (j = 0; j < num_loc; j++)
    {
      loc_list [j] = SeqLocFree (loc_list [j]);
    }
    loc_list = MemFree (loc_list);
  }
  return loc_list;
}

static SeqLocPtr 
ReadAlignedSeqLocList
(TagListPtr      tlp,
 IntervalPagePtr ipp,
 Boolean         partial5,
 Boolean         partial3, 
 Boolean         nullsBetween)
{
  Int4             from, to, aln_from, aln_to, aln_row, aln_len;
  Boolean          fuzz_after;
  Boolean          fuzz_before;
  Int2             fuzz_from;
  Int2             fuzz_to;
  Boolean          isInterval;
  Boolean          isPoint;
  Boolean          notFirst;
  Boolean          okay;
  SeqIdPtr         seqid;
  SeqLocPtr        master_slp, tmp_slp;
  Uint2            strand, prev_strand, aln_strand;
  ValNodePtr       vnp;
  SeqAlignPtr      salp, prev_salp;
  BioseqPtr        bsp;
  SeqLocPtr PNTR   loc_list;
  Int4             loc_num, max_locs;
  Boolean          asked_about_repair = FALSE;
  MsgAnswer        ans = ANS_YES;

  if (tlp == NULL || ipp == NULL || ipp->salp_list == NULL)
  {
    return NULL;
  }
  
  prev_salp = NULL;
  prev_strand = Seq_strand_unknown;
  master_slp = NULL;
  
  max_locs = 0;
  for (vnp = tlp->vnp; vnp != NULL; vnp = vnp->next)
  {
    if (!ReadSeqAlignFromDialogLine ((CharPtr) vnp->data.ptrvalue, ipp->aln_col,
                                     &salp, &prev_salp, 
                                     ipp->salp_list, ipp->aln_count))
    {
      return NULL;
    }
    else
    {
      max_locs = MAX (max_locs, salp->dim);
    }
  }
  if (max_locs == 0)
  {
    return NULL;
  }
  
  loc_list = (SeqLocPtr PNTR) MemNew (max_locs * sizeof (SeqLocPtr));
  if (loc_list == NULL)
  {
    return NULL;
  }

  prev_salp = NULL;

  okay = TRUE;
  notFirst = FALSE;
  for (vnp = tlp->vnp; vnp != NULL && okay; vnp = vnp->next) {
    if (StringHasNoText (vnp->data.ptrvalue))
    {
      continue;
    }
    fuzz_from = -1;
    fuzz_to = -1;
    fuzz_before = FALSE;
    fuzz_after = FALSE;
    from = 0;
    to = 0;
    isInterval = TRUE;
    isPoint = FALSE;
    strand = Seq_strand_unknown;
    seqid = NULL;
    salp = NULL;
    aln_row = 0;
    if (! ReadFromValueFromDialogLine ((CharPtr) vnp->data.ptrvalue,
                                           &from, &fuzz_after))
    {
      continue;
    }
    okay = ReadToValueFromDialogLine ((CharPtr) vnp->data.ptrvalue,
                                       &to, &from, &fuzz_after,
                                       &isInterval, &isPoint, 
                                       partial5, partial3)
           && ReadStrandFromDialogLine ((CharPtr) vnp->data.ptrvalue,
                                        ipp->strand_col,
                                        &strand, &prev_strand)
           && ReadSeqAlignFromDialogLine ((CharPtr) vnp->data.ptrvalue, ipp->aln_col,
                                              &salp, &prev_salp, 
                                              ipp->salp_list, ipp->aln_count);
    if (okay)
    {
      for (loc_num = 0; loc_num < salp->dim; loc_num++)
      {
        aln_row = loc_num + 1;
        aln_from = from;
        aln_to = to;
        aln_strand = strand;
        seqid = AlnMgr2GetNthSeqIdPtr (salp, aln_row);
        okay = ((bsp = BioseqFind (seqid)) != NULL)
                 && AdjustStrandForAlignment (&aln_strand, salp, aln_row)
                 && (aln_len = SeqAlignLength (salp)) > 0;
        if (okay)
        {
          if (ipp->allow_nulls_in_list)
          {
            if (! AdjustFromForGap (&aln_from, salp, aln_len, aln_row)
                || ! AdjustToForGap (&aln_to, salp, aln_row)
                || ! AddLocToList (seqid, aln_strand, aln_from, aln_to,
                                   (nullsBetween  && notFirst),
                                   fuzz_from, fuzz_to,
                                   fuzz_before, fuzz_after,
                                   partial5, partial3,
                                   &(loc_list [loc_num])))
            {
              loc_list [loc_num] = ValNodeNew (NULL);
              loc_list [loc_num]->choice = SEQLOC_NULL;
            }
          }
          else
          {
            okay = AdjustFromForGap (&aln_from, salp, aln_len, aln_row)
                   && AdjustToForGap (&aln_to, salp, aln_row)
                   && AddLocToList (seqid, aln_strand, aln_from, aln_to,
                              (nullsBetween  && notFirst),
                              fuzz_from, fuzz_to,
                              fuzz_before, fuzz_after,
                              partial5, partial3,
                              &(loc_list [loc_num]));
          }
        }
      }
    }
      
    notFirst = TRUE;
  }
  
  if (!okay)
  {
    loc_list = FreeSeqLocArray (loc_list, max_locs);
    return NULL;
  }
  
  /* now fix intervals that are out of order and set partials for locations */
  for (loc_num = 0; loc_num < max_locs && loc_list [loc_num] != NULL; loc_num++)
  {
    bsp = BioseqFindFromSeqLoc (loc_list [loc_num]);
    if (bsp != NULL && SeqLocBadSortOrder (bsp, loc_list [loc_num])) 
    {
      if (!asked_about_repair)
      {
        ans = Message (MSG_YN,
            "Feature location intervals are out of order.  Do you want them repaired?");
        asked_about_repair = TRUE;
      }
      if (ans == ANS_YES) {
        tmp_slp = SeqLocMerge (bsp, loc_list [loc_num], NULL, FALSE, FALSE, nullsBetween);
        loc_list [loc_num] = SeqLocFree (loc_list [loc_num]);
        loc_list [loc_num] = tmp_slp;
        if (bsp->repr == Seq_repr_seg) {
          tmp_slp = SegLocToParts (bsp, loc_list [loc_num]);
          loc_list [loc_num] = SeqLocFree (loc_list [loc_num]);
          loc_list [loc_num] = tmp_slp;
        }
        FreeAllFuzz (loc_list [loc_num]);
      }
    }
    SetSeqLocPartial (loc_list [loc_num], partial5, partial3);
  }
  
  /* now make chain */
  master_slp = loc_list [0];
  tmp_slp = loc_list [0];
  for (tmp_slp = loc_list [0], loc_num = 1; 
       tmp_slp != NULL && loc_num < max_locs;
       loc_num++)
  {
    tmp_slp->next = loc_list [loc_num];
    tmp_slp = tmp_slp->next;
  }
  
  loc_list = MemFree (loc_list);
  
  return master_slp;
}

static Pointer IntervalPageToSeqLocPtr (DialoG d)

{
  IntervalPagePtr  ipp;
  Boolean          nullsBetween;
  Boolean          partial5;
  Boolean          partial3;
  SeqLocPtr        master_slp;
  TagListPtr       tlp;

  ipp = (IntervalPagePtr) GetObjectExtra (d);
  if (ipp == NULL) return NULL;
  tlp = GetObjectExtra (ipp->ivals);
  if (tlp == NULL) return NULL;
  

  nullsBetween = GetStatus (ipp->nullsBetween);
  partial5 = GetStatus (ipp->partial5);
  partial3 = GetStatus (ipp->partial3);
  
  if (ipp->seqid_col > -1)
  {
    master_slp = ReadSingleSeqLoc (tlp, ipp, partial5, partial3, nullsBetween);
  }
  else
  {
    master_slp = ReadAlignedSeqLocList (tlp, ipp, partial5, partial3, nullsBetween);
  }
  return (Pointer) master_slp;
}

static void CleanupIntervalPage (GraphiC g, VoidPtr data)

{
  IntervalPagePtr  ipp;
  Int2             j;

  ipp = (IntervalPagePtr) data;
  if (ipp != NULL) {
    /* free seq ID list */
    if (ipp->sip_list != NULL)
    {
      for (j = 0; j <= ipp->count + 1; j++) {
        ipp->sip_list [j] = SeqIdFree (ipp->sip_list [j]);
      }
    }
    MemFree (ipp->sip_list);
    if (ipp->alist != NULL) {
      for (j = 0; j <= ipp->count + 1; j++) {
        MemFree (ipp->alist [j].name);
      }
    }
    MemFree (ipp->alist);
    MemFree (ipp->lengths);

    /* free list of alignments */
    if (ipp->salp_list != NULL)
    {
      for (j = 0; j <= ipp->aln_count + 1; j++) {
        ipp->salp_list [j] = SeqAlignFree (ipp->salp_list [j]); 
      }
    }
    MemFree (ipp->salp_list);
    /* free alignment tags */
    if (ipp->aln_alist != NULL) {
      for (j = 0; j <= ipp->aln_count + 1; j++) {
        MemFree (ipp->aln_alist [j].name);
      }
    }
    MemFree (ipp->aln_alist);
    /* free alignment lengths */
    MemFree (ipp->aln_lengths);

    
    /* free callback list */
    MemFree (ipp->callbacks);
  }
  MemFree (data);
}

Uint2 interval_types [] = {
  TAGLIST_TEXT, TAGLIST_TEXT, TAGLIST_POPUP, TAGLIST_POPUP, TAGLIST_POPUP
};

Uint2 interval_widths [] = {
  5, 5, 0, 0, 0
};

static Boolean ReadSeqLocDialog (DialoG d, CharPtr filename)

{
  AsnIoPtr         aip;
  IntervalPagePtr  ipp;
  SeqLocPtr        slp;
  Char             path [PATH_MAX];

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  ipp = (IntervalPagePtr) GetObjectExtra (d);
  if (ipp != NULL) {
    if (path [0] != '\0' || GetInputFileName (path, sizeof (path), "", "TEXT")) {
      aip = AsnIoOpen (path, "r");
      if (aip != NULL) {
        slp = SeqLocAsnRead (aip, NULL);
        AsnIoClose (aip);
        if (slp != NULL) {
          SeqLocPtrToIntervalPage (ipp->dialog, slp);
          slp = SeqLocFree (slp);
          Update ();
          return TRUE;
        }
      }
    }
  }
  return FALSE;
}

static Boolean WriteSeqLocDialog (DialoG d, CharPtr filename)

{
  AsnIoPtr         aip;
  IntervalPagePtr  ipp;
  SeqLocPtr        slp;
  Char             path [PATH_MAX];
#ifdef WIN_MAC
  FILE            *f;
#endif

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  ipp = (IntervalPagePtr) GetObjectExtra (d);
  if (ipp != NULL) {
    if (path [0] != '\0' || GetOutputFileName (path, sizeof (path), NULL)) {
#ifdef WIN_MAC
      f = FileOpen (path, "r");
      if (f != NULL) {
        FileClose (f);
      } else {
        FileCreate (path, "TEXT", "ttxt");
      }
#endif
      aip = AsnIoOpen (path, "w");
      if (aip != NULL) {
        slp = IntervalPageToSeqLocPtr (ipp->dialog);
        SeqLocAsnWrite (slp, aip, NULL);
        AsnIoClose (aip);
        slp = SeqLocFree (slp);
        return TRUE;
      }
    }
  }
  return FALSE;
}

static void SetOnlySequenceAndStrand (IntervalPagePtr ipp)
{
  TagListPtr       tlp;
  ValNodePtr       vnp;
  CharPtr          cp;
  CharPtr          tabptr;
  ValNodePtr       saved_list;
  
  if (ipp == NULL) return;
  tlp = GetObjectExtra (ipp->ivals);
  if (tlp == NULL) return;
  saved_list = tlp->vnp;
  tlp->vnp = NULL;
  SendMessageToDialog (tlp->dialog, VIB_MSG_RESET);
  tlp->vnp = saved_list;
  for (vnp = tlp->vnp; vnp != NULL; vnp = vnp->next) 
  {
    cp = vnp->data.ptrvalue;
    if (cp != NULL)
    {
      tabptr = StringChr (cp, '\t');
      if (tabptr != NULL)
      {
      	tabptr = StringChr (tabptr + 1, '\t');
      }
      if (tabptr != NULL)
      {
      	cp[0] = '\t';
      	cp++;
      	while (*tabptr != 0)
      	{
      	  *cp = *tabptr;
      	  cp++;
      	  tabptr++;
      	}
      	*cp = 0;
      }
    }
  }

  SendMessageToDialog (tlp->dialog, VIB_MSG_REDRAW);  
}

static void ClearLocationPartialCheckboxes (IntervalPagePtr ipp)
{
  if (ipp != NULL)
  {
    SetStatus (ipp->partial5, FALSE);
    SetStatus (ipp->partial3, FALSE); 
    SetStatus (ipp->nullsBetween, FALSE); 
    if (ipp->proc != NULL) {
      ipp->proc (ipp->ffp, FALSE, FALSE, FALSE);
    }
  }
}

static void IntervalEditorMessage (DialoG d, Int2 mssg)

{
  IntervalPagePtr  ipp;

  ipp = (IntervalPagePtr) GetObjectExtra (d);
  if (ipp != NULL) {
    if (mssg == VIB_MSG_INIT) {
      SeqLocPtrToIntervalPage (d, NULL);
    } else if (mssg == VIB_MSG_ENTER) {
      SendMessageToDialog (ipp->ivals, VIB_MSG_ENTER);
    } else if (mssg == VIB_MSG_RESET) {
    }
    else if (mssg == NUM_VIB_MSG + 1)
    {
      SetOnlySequenceAndStrand (ipp);
    }
    else if (mssg == NUM_VIB_MSG + 2)
    {
      ClearLocationPartialCheckboxes (ipp);
    }
  }
}

extern void StdFeatIntEdPartialCallback (FeatureFormPtr ffp, Boolean partial5, Boolean partial3, Boolean order)

{
  if (ffp == NULL) return;
  SafeSetStatus (ffp->partial, (partial5 || partial3 || order));
  Update ();
}

static void ChangedPartialProc (ButtoN b)

{
  IntervalPagePtr  ipp;

  ipp = (IntervalPagePtr) GetObjectExtra (b);
  if (ipp == NULL) return;
  if (ipp->proc != NULL) {
    ipp->proc (ipp->ffp, GetStatus (ipp->partial5),
               GetStatus (ipp->partial3), GetStatus (ipp->nullsBetween));
  }
}

extern DialoG CreateTagListDialogEx (GrouP h, Uint2 rows, Uint2 cols,
                                     Int2 spacing, Uint2Ptr types,
                                     Uint2Ptr textWidths, EnumFieldAssocPtr PNTR alists,
                                     Boolean useBar, Boolean noExtend,
                                     ToDialogFunc tofunc, FromDialogFunc fromfunc);

extern DialoG CreateIntervalEditorDialogExEx (GrouP h, CharPtr title, Uint2 rows,
                                              Int2 spacing, SeqEntryPtr sep,
                                              Boolean nucsOK, Boolean protsOK,
                                              Boolean useBar, Boolean showPartials,
                                              Boolean allowGaps, FeatureFormPtr ffp,
                                              IntEdPartialProc proc, 
                                              Boolean use_aln, Boolean show_seqid,
                                              TaglistCallback tlp_callback, 
                                              Pointer callback_data,
                                              Boolean allow_nulls_in_list)

{
  BioseqPtr        bsp;
  Int4             count;
  GrouP            f;
  IntervalPagePtr  ipp;
  Int2             j;
  GrouP            m;
  GrouP            p;
  PrompT           p1;
  PrompT           p2;
  PrompT           p3;
  CharPtr          ptr;
  GrouP            q;
  GrouP            s;
  Boolean          showIdTags;
  SeqIdPtr         sip;
  Char             str [128];
  TagListPtr       tlp;
  SeqAlignPtr      salp_list = NULL, salp;
  Int4             num_cols;

  p = HiddenGroup (h, 1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  ipp = (IntervalPagePtr) MemNew (sizeof (IntervalPage));
  if (ipp != NULL) {

    SetObjectExtra (p, ipp, CleanupIntervalPage);
    ipp->dialog = (DialoG) p;
    ipp->todialog = SeqLocPtrToIntervalPage;
    ipp->fromdialog = IntervalPageToSeqLocPtr;
    ipp->dialogmessage = IntervalEditorMessage;
    ipp->testdialog = TestIntervalEditor;
    ipp->importdialog = ReadSeqLocDialog;
    ipp->exportdialog = WriteSeqLocDialog;

    ipp->ffp = ffp;
    ipp->proc = proc;
     
    ipp->allow_nulls_in_list = allow_nulls_in_list;

    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    m = HiddenGroup (s, -1, 0, NULL);
    /*
    SetGroupSpacing (m, 10, 10);
    */

    ipp->nucsOK = nucsOK;
    ipp->protsOK = protsOK;
    ipp->showIdTags = FALSE;
    ipp->count = 0;
        
    if (ipp->nucsOK)
    {
      ipp->strand_col = 2;
    }
    else
    {
      ipp->strand_col = -1;
    }
    
    ipp->aln_col = -1;
    if (show_seqid)
    {
      if (ipp->nucsOK)
      {
        ipp->seqid_col = 3;
      }
      else
      {
        ipp->seqid_col = 2;
      }
      if (use_aln)
      {
        ipp->aln_col = ipp->seqid_col + 1;
      }
    }
    else
    {
      ipp->seqid_col = -1;
      if (use_aln)
      {
        if (ipp->nucsOK)
        {
          ipp->aln_col = 3;
        }
        else
        {
          ipp->aln_col = 2;
        }
      }
    }
    
    /* set up callbacks */
    ipp->callbacks = (TaglistCallback PNTR) MemNew (5 * sizeof (TaglistCallback));
    if (ipp->callbacks != NULL)
    {
      ipp->callbacks [0] = NULL;
      ipp->callbacks [1] = NULL;
      ipp->callbacks [2] = NULL;
      ipp->callbacks [3] = NULL;
      ipp->callbacks [4] = NULL;
  
      if (ipp->seqid_col > -1)
      {
        ipp->callbacks [ipp->seqid_col] = tlp_callback;
      }
      if (ipp->aln_col > -1)
      {
        ipp->callbacks [ipp->aln_col] = tlp_callback;
      }
    }

    if (sep != NULL) {
      if (use_aln)
      {
        VisitAnnotsInSep (sep, &salp_list, GetAlignmentsInSeqEntryCallback);
        count = 4;
        salp = salp_list;
        while (salp != NULL)
        {
          count++;
          salp = salp->next;
        }
        
        ipp->salp_list = (SeqAlignPtr PNTR) MemNew (sizeof (SeqAlignPtr) * count);
        ipp->aln_alist = MemNew (sizeof (EnumFieldAssoc) * (size_t) count);
        ipp->aln_lengths = MemNew (sizeof (Int4) * (size_t) count);
        ipp->aln_count = 0;
        if (ipp->salp_list != NULL && ipp->aln_alist != NULL && ipp->aln_lengths != NULL)
        {
          /* first one is NULL */
          ipp->salp_list [0] = NULL;
          ipp->aln_lengths [0] = 0;
          ipp->aln_alist [0].name = StringSave ("     ");
          ipp->aln_alist [0].value = (UIEnum) 0;
          ipp->aln_count ++;
          salp = salp_list;
          while (salp != NULL)
          {
            ipp->salp_list [ipp->aln_count] = salp;
            ipp->aln_lengths [ipp->aln_count] = SeqAlignLength (ipp->salp_list [ipp->aln_count]);
            CreateSeqAlignLabel (ipp->salp_list [ipp->aln_count], str, 30);
            ipp->aln_alist [ipp->aln_count].name = StringSave (str);
            ipp->aln_alist [ipp->aln_count].value = (UIEnum) ipp->aln_count;
                        
            salp = salp->next;
            /* sever chain */
            ipp->salp_list [ipp->aln_count]->next = NULL;
            ipp->aln_count++;
          }
        }
        /* add end of list marker */
        ipp->aln_alist [ipp->aln_count].name = NULL;
        ipp->aln_alist [ipp->aln_count].value = (UIEnum) 0;

      }
      else
      {
        ipp->salp_list = NULL;
        ipp->aln_alist = NULL;
        ipp->aln_lengths = NULL;
        ipp->aln_count = 0;
      }
      count = SegmentedEntryCount (sep);
      count += 4;
      ipp->sip_list = MemNew (sizeof (SeqIdPtr) * (size_t) count);
      ipp->alist = MemNew (sizeof (EnumFieldAssoc) * (size_t) count);
      ipp->lengths = MemNew (sizeof (Int4) * (size_t) count);
      ipp->count = 0;

      if (ipp->sip_list != NULL && ipp->alist != NULL && ipp->lengths != NULL) {
        VisitBioseqsInSep (sep, NULL, ClearBspScratch);
        SeqEntryExplore (sep, (Pointer) ipp, FillInProducts);
        VisitBioseqsInSep (sep, NULL, ClearBspScratch);
        j = 0;
        ipp->alist [j].name = StringSave ("     ");
        ipp->alist [j].value = (UIEnum) 0;
        for (j = 1; j <= ipp->count; j++) {
          sip = ipp->sip_list [j];
          if (sip != NULL) {
            bsp = BioseqFindCore (sip);
            if (bsp == NULL) {
              bsp = BioseqLockById (sip);
              BioseqUnlockById (sip);
            }
            if (bsp != NULL)
            {
              ipp->lengths [j] = bsp->length;
              sip = SeqIdFindWorst (bsp->id);
            }
            SeqIdWrite (sip, str, PRINTID_REPORT, sizeof (str));
            ptr = StringChr (str, '|');
            showIdTags = FALSE;
            if (ptr == NULL) {
              ptr = str;
            } else if (showIdTags) {
              ptr = str;
            } else {
              ptr++;
            }
            ipp->alist [j].name = StringSave (ptr);
            ipp->alist [j].value = (UIEnum) j;
          }
        }
        /* add end of list marker */
        j = ipp->count + 1;
        ipp->alist [j].name = NULL;
        ipp->alist [j].value = (UIEnum) 0;
      }
#ifdef WIN_MOTIF
      if (ipp->count > 31) {
        if (ipp->seqid_col > -1)
        {
          interval_types [ipp->seqid_col] = TAGLIST_LIST;
        }
        if (ipp->aln_col > -1)
        {
          interval_types [ipp->aln_col] = TAGLIST_LIST;
        }
      } else {
        interval_types [2] = TAGLIST_POPUP;
        interval_types [3] = TAGLIST_POPUP;
        interval_types [4] = TAGLIST_POPUP;
      }
#endif

    } else {
      ipp->alist = MemNew (sizeof (EnumFieldAssoc) * (size_t) 4);
      if (ipp->alist != NULL) {
        j = 0;
        ipp->alist [j].name = StringSave ("     ");
        ipp->alist [j].value = (UIEnum) 0;
        j = 1;
        ipp->alist [j].name = NULL;
        ipp->alist [j].value = (UIEnum) 0;

      }
    }
    
    ipp->alists [0] = NULL;
    ipp->alists [1] = NULL;
    ipp->alists [2] = NULL;
    ipp->alists [3] = NULL;
    ipp->alists [4] = NULL;
    if (ipp->strand_col > -1)
    {
      ipp->alists [ipp->strand_col] = strand_alist;
    }
    if (ipp->seqid_col > -1)
    {
      ipp->alists [ipp->seqid_col] = ipp->alist;
    }
    if (ipp->aln_col > -1)
    {
      ipp->alists [ipp->aln_col] = ipp->aln_alist;
    }

    q = NULL;
    if (showPartials) {
      q = HiddenGroup (m, 4, 0, NULL);
      SetGroupSpacing (q, 20, 2);
      if (nucsOK) {
        ipp->partial5 = CheckBox (q, "5' Partial", ChangedPartialProc);
      } else {
        ipp->partial5 = CheckBox (q, "NH2 Partial", ChangedPartialProc);
      }
      SetObjectExtra (ipp->partial5, ipp, NULL);
      if (nucsOK) {
        ipp->partial3 = CheckBox (q, "3' Partial", ChangedPartialProc);
      } else {
        ipp->partial3 = CheckBox (q, "CO2H Partial", ChangedPartialProc);
      }
      SetObjectExtra (ipp->partial3, ipp, NULL);
    }

    f = HiddenGroup (m, 5, 0, NULL);
    StaticPrompt (f, "From", 5 * stdCharWidth, 0, programFont, 'c');
    StaticPrompt (f, "To", 5 * stdCharWidth, 0, programFont, 'c');
    p1 = NULL;
    p2 = NULL;
    p3 = NULL;
    if (ipp->strand_col > -1)
    {
      p1 = StaticPrompt (f, "Strand", 0, 0, programFont, 'c');
    }
    if (ipp->seqid_col > -1)
    {
      p2 = StaticPrompt (f, "SeqID", 0, 0, programFont, 'c');
    }
    if (ipp->aln_col > -1)
    {
      p3 = StaticPrompt (f, "Alignment", 0, 0, programFont, 'c');
    }

    f = HiddenGroup (m, 0, 4, NULL);
    SetGroupSpacing (f, 0, 0);

    num_cols = 2;
    if (ipp->strand_col > -1)
    {
      num_cols ++;
    }
    if (ipp->seqid_col > -1)
    {
      num_cols ++;
    }
    if (ipp->aln_col > -1)
    {
      num_cols ++;
    }
    
    ipp->ivals = CreateTagListDialogExEx (f, rows, num_cols, spacing,
                                        interval_types, interval_widths, ipp->alists,
                                        useBar, FALSE, NULL, NULL,
                                        ipp->callbacks, callback_data, FALSE);

    /* put back static interval_types values that may have been changed */
    interval_types [2] = TAGLIST_POPUP;
    interval_types [3] = TAGLIST_POPUP;
    interval_types [4] = TAGLIST_POPUP;
    ipp->nullsBetween = CheckBox (m, "'order' (intersperse intervals with gaps)", ChangedPartialProc);
    SetObjectExtra (ipp->nullsBetween, ipp, NULL);

    AlignObjects (ALIGN_CENTER, (HANDLE) ipp->ivals,
                  (HANDLE) q, (HANDLE) ipp->nullsBetween, NULL);
    tlp = (TagListPtr) GetObjectExtra (ipp->ivals);
    if (tlp != NULL) {
      if (ipp->strand_col > -1)
      {
        AlignObjects (ALIGN_JUSTIFY, (HANDLE) tlp->control [ipp->strand_col], (HANDLE) p1, NULL);
      }
      if (ipp->seqid_col > -1)
      {
        AlignObjects (ALIGN_JUSTIFY, (HANDLE) tlp->control [ipp->seqid_col], (HANDLE) p2, NULL);
      }
      if (ipp->aln_col > -1)
      {
        AlignObjects (ALIGN_JUSTIFY, (HANDLE) tlp->control [ipp->aln_col], (HANDLE) p3, NULL);
      }
    }
  }

  return (DialoG) p;
}

typedef struct intervalchoice 
{
  DIALOG_MESSAGE_BLOCK
  GrouP  seq_or_aln;
  DialoG aln_dlg;
  DialoG seq_dlg;  
} IntervalChoiceData, PNTR IntervalChoicePtr;

static void IntervalChoiceCallback (Pointer userdata)
{
  IntervalChoicePtr dlg;
  SeqLocPtr         slp, slp_tmp;
  Boolean           ok_for_aln = FALSE;
  BioseqPtr         bsp;
  SeqAlignPtr       salp, salp_next;

  if (userdata == NULL)
  {
    return;
  }

  dlg = (IntervalChoicePtr) userdata;
  
  if (GetValue (dlg->seq_or_aln) == 1)
  {
    slp = DialogToPointer (dlg->seq_dlg);
    if (slp != NULL)
    {
      ok_for_aln = TRUE;
      slp_tmp = SeqLocFindNext (slp, NULL);
      while (slp_tmp != NULL && ok_for_aln)
      {
        bsp = BioseqFindFromSeqLoc (slp_tmp);
        salp = FindAlignmentsForBioseq (bsp);
        if (salp == NULL)
        {
          ok_for_aln = FALSE;
        }
        else
        {
          while (salp != NULL)
          {
            salp_next = salp->next;
            salp->next = NULL;
            salp = SeqAlignFree (salp);
            salp = salp_next;
          }
        }
        slp_tmp = SeqLocFindNext (slp, slp_tmp);
      }
    }
    if (ok_for_aln)  
    {
      Enable (dlg->seq_or_aln);
    }
    else
    {
      Disable (dlg->seq_or_aln);
    }
  }
}

static void SeqLocToIntervalChoiceEditor (DialoG d, Pointer userdata)
{
  IntervalChoicePtr dlg;
  
  dlg = (IntervalChoicePtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  PointerToDialog (dlg->seq_dlg, userdata);
  PointerToDialog (dlg->aln_dlg, userdata); 
  IntervalChoiceCallback (dlg);
}

static Pointer IntervalChoiceEditorToSeqLoc (DialoG d)
{
  IntervalChoicePtr dlg;
  
  dlg = (IntervalChoicePtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  if (GetValue (dlg->seq_or_aln) == 1)
  {
    return DialogToPointer (dlg->seq_dlg);
  }
  else
  {
    return DialogToPointer (dlg->aln_dlg);
  }
}


static Boolean WriteIntervalChoiceDialog (DialoG d, CharPtr filename)

{
  AsnIoPtr          aip;
  IntervalChoicePtr dlg;
  SeqLocPtr         slp;
  Char              path [PATH_MAX];
#ifdef WIN_MAC
  FILE            *f;
#endif

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  dlg = (IntervalChoicePtr) GetObjectExtra (d);
  if (dlg != NULL) {
    if (path [0] != '\0' || GetOutputFileName (path, sizeof (path), NULL)) {
#ifdef WIN_MAC
      f = FileOpen (path, "r");
      if (f != NULL) {
        FileClose (f);
      } else {
        FileCreate (path, "TEXT", "ttxt");
      }
#endif
      aip = AsnIoOpen (path, "w");
      if (aip != NULL) {
        slp = IntervalChoiceEditorToSeqLoc (dlg->dialog);
        SeqLocAsnWrite (slp, aip, NULL);
        AsnIoClose (aip);
        slp = SeqLocFree (slp);
        return TRUE;
      }
    }
  }
  return FALSE;
}


static Boolean ReadIntervalChoiceDialog (DialoG d, CharPtr filename)

{
  AsnIoPtr          aip;
  IntervalChoicePtr dlg;
  SeqLocPtr         slp;
  Char              path [PATH_MAX];

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  dlg = (IntervalChoicePtr) GetObjectExtra (d);
  if (dlg != NULL) {
    if (path [0] != '\0' || GetInputFileName (path, sizeof (path), "", "TEXT")) {
      aip = AsnIoOpen (path, "r");
      if (aip != NULL) {
        slp = SeqLocAsnRead (aip, NULL);
        AsnIoClose (aip);
        if (slp != NULL) {
          SeqLocToIntervalChoiceEditor (dlg->dialog, slp);
          slp = SeqLocFree (slp);
          Update ();
          return TRUE;
        }
      }
    }
  }
  return FALSE;
}


static void IntervalChoiceEditorMessage (DialoG d, Int2 mssg)

{
  IntervalChoicePtr dlg;

  dlg = (IntervalChoicePtr) GetObjectExtra (d);
  if (dlg != NULL) {  
    if (GetValue (dlg->seq_or_aln) == 1)
    {
      SendMessageToDialog (dlg->seq_dlg, mssg);
    }
    else
    {
      SendMessageToDialog (dlg->aln_dlg, mssg);
    }
  }
}

static ValNodePtr TestIntervalChoiceEditor (DialoG d)
{
  IntervalChoicePtr dlg;

  dlg = (IntervalChoicePtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  else if (GetValue (dlg->seq_or_aln) == 1)
  {
    return TestDialog (dlg->seq_dlg);
  }
  else
  {
    return TestDialog (dlg->aln_dlg);
  }
}

static void DisplayErrorMessagesOk (ButtoN b)
{
  BoolPtr pdone;
  
  pdone = (BoolPtr) GetObjectExtra (b);
  if (pdone != NULL)
  {
    *pdone = TRUE;
  }
}

extern void DisplayErrorMessages (CharPtr title, ValNodePtr err_list)
{
  WindoW     w;
  GrouP      h;
  DoC        doc;
  ButtoN     b;
  ValNodePtr vnp;
  Boolean    done = FALSE;
  
  w = MovableModalWindow(-20, -13, -10, -10, title, NULL);
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  doc = DocumentPanel (h, stdCharWidth * 27, stdLineHeight * 8);
  SetDocAutoAdjust (doc, TRUE);
  for (vnp = err_list; vnp != NULL; vnp = vnp->next)
  {
    if (!StringHasNoText (vnp->data.ptrvalue))
    {
      AppendText (doc, vnp->data.ptrvalue, NULL, NULL, programFont);
    }
  }
  b = PushButton (h, "OK", DisplayErrorMessagesOk);
  SetObjectExtra (b, &done, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) doc,
                              (HANDLE) b, 
                              NULL);
  Show (w);
  Select (w);
  while (!done)
  {
    ProcessExternalEvent ();
    Update ();
  }
  ProcessAnEvent ();
  Remove (w);
}

static void ShowIntervalChoice (IntervalChoicePtr dlg)
{
  if (dlg == NULL)
  {
    return;
  }
  if (GetValue (dlg->seq_or_aln) == 1)
  {
    Show (dlg->seq_dlg);
    Hide (dlg->aln_dlg);
  }
  else
  {
    Show (dlg->aln_dlg);
    Hide (dlg->seq_dlg);
  }
}

static void ChangeIntervalChoice (GrouP g)
{
  IntervalChoicePtr dlg;
  SeqLocPtr         slp;
  ValNodePtr        err_list = NULL;
  CharPtr           title = NULL;

  dlg = (IntervalChoicePtr) GetObjectExtra (g);
  if (dlg == NULL)
  {
    return;
  }
  if (GetValue (dlg->seq_or_aln) == 1)
  {
    title = "Unable to translate to sequence coordinates";
    err_list = TestDialog (dlg->aln_dlg);
    if (err_list == NULL)
    {
      slp = DialogToPointer (dlg->aln_dlg);
      PointerToDialog (dlg->seq_dlg, slp);
      slp = SeqLocFree (slp);
      err_list = TestDialog (dlg->seq_dlg);
    }    
    if (err_list != NULL)
    {
      SetValue (dlg->seq_or_aln, 2);       
    }
  }
  else
  {
    title = "Unable to translate to alignment coordinates";
    err_list = TestDialog (dlg->seq_dlg);
    if (err_list == NULL)
    {
      slp = DialogToPointer (dlg->seq_dlg);
      PointerToDialog (dlg->aln_dlg, slp);
      slp = SeqLocFree (slp);
      err_list = TestDialog (dlg->aln_dlg);
    }
    if (err_list != NULL)
    {
      SetValue (dlg->seq_or_aln, 1);
    }
  }
  ShowIntervalChoice (dlg);
  if (err_list != NULL)
  {
    DisplayErrorMessages (title, err_list);
    err_list = ValNodeFreeData (err_list);
  }
}

static DialoG CreateIntervalEditorDialogAlnChoice (GrouP h, CharPtr title, Uint2 rows,
                                            Int2 spacing, SeqEntryPtr sep,
                                            Boolean nucsOK, Boolean protsOK,
                                            Boolean useBar, Boolean showPartials,
                                            Boolean allowGaps, FeatureFormPtr ffp,
                                            IntEdPartialProc proc)
{
  SeqAlignPtr       salp_list = NULL;
  GrouP             p, g;
  IntervalChoicePtr dlg;
  
  VisitAnnotsInSep (sep, &salp_list, GetAlignmentsInSeqEntryCallback);
  if (salp_list == NULL)
  {
    return CreateIntervalEditorDialogExEx (h, title, rows, spacing, sep,
                                         nucsOK, protsOK, useBar, showPartials,
                                         allowGaps, ffp, proc, FALSE, TRUE,
                                         NULL, NULL, FALSE);
  }
  else
  {
    dlg = (IntervalChoicePtr) MemNew (sizeof (IntervalChoiceData));
    if (dlg == NULL)
    {
      return NULL;
    }
  
    p = HiddenGroup (h, -1, 0, NULL);
    SetGroupSpacing (p, 10, 10);
    SetObjectExtra (p, dlg, StdCleanupExtraProc);
    
    dlg->dialog = (DialoG) p;
    dlg->todialog = SeqLocToIntervalChoiceEditor;
    dlg->fromdialog = IntervalChoiceEditorToSeqLoc;
    dlg->dialogmessage = IntervalChoiceEditorMessage;
    dlg->testdialog = TestIntervalChoiceEditor;
    dlg->exportdialog = WriteIntervalChoiceDialog;
    dlg->importdialog = ReadIntervalChoiceDialog;

    g = HiddenGroup (p, 0, 0, NULL);
    dlg->seq_dlg = CreateIntervalEditorDialogExEx (g, title, rows, spacing, sep,
                                         nucsOK, protsOK, useBar, showPartials,
                                         allowGaps, ffp, proc, FALSE, TRUE,
                                         IntervalChoiceCallback, dlg, FALSE);
                                         
    dlg->aln_dlg = CreateIntervalEditorDialogExEx (g, title, rows, spacing, sep,
                                         nucsOK, protsOK, useBar, showPartials,
                                         allowGaps, ffp, proc, TRUE, TRUE,
                                         IntervalChoiceCallback, dlg, FALSE);
    AlignObjects (ALIGN_CENTER, (HANDLE)dlg->seq_dlg, (HANDLE) dlg->aln_dlg, NULL);                                         
    dlg->seq_or_aln = HiddenGroup (p, 2, 0, ChangeIntervalChoice);
    SetObjectExtra (dlg->seq_or_aln, dlg, NULL);
    RadioButton (dlg->seq_or_aln, "Sequence Coordinates");
    RadioButton (dlg->seq_or_aln, "Alignment Coordinates");
    SetValue (dlg->seq_or_aln, 1);
    
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) dlg->seq_or_aln, NULL);
    ShowIntervalChoice (dlg);
    IntervalChoiceCallback (dlg);
    return (DialoG) p; 
  }
}

extern DialoG CreateIntervalEditorDialogEx (GrouP h, CharPtr title, Uint2 rows,
                                            Int2 spacing, SeqEntryPtr sep,
                                            Boolean nucsOK, Boolean protsOK,
                                            Boolean useBar, Boolean showPartials,
                                            Boolean allowGaps, FeatureFormPtr ffp,
                                            IntEdPartialProc proc)
{
  return CreateIntervalEditorDialogAlnChoice (h, title, rows, spacing, sep,
                                         nucsOK, protsOK, useBar, showPartials,
                                         allowGaps, ffp, proc);
}

extern DialoG CreateIntervalEditorDialog (GrouP h, CharPtr title, Uint2 rows,
                                          Int2 spacing, SeqEntryPtr sep,
                                          Boolean nucsOK, Boolean protsOK)

{
  return CreateIntervalEditorDialogEx (h, title, rows, spacing, sep,
                                       nucsOK, protsOK, TRUE, TRUE, FALSE,
                                       NULL, NULL);
}

static void ValNodePtrToVisStringDialog (DialoG d, Pointer data)

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

static Pointer VisStringDialogToValNodePtr (DialoG d)

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
          list->data.ptrvalue = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, 0);
        }
      }
    }
  }
  return (Pointer) head;
}

static void DbtagPtrToDbtagDialog (DialoG d, Pointer data, Boolean readOnly)

{
  DbtagPtr     dp;
  ValNodePtr   head;
  Int2         j;
  size_t       len;
  ValNodePtr   list;
  ObjectIdPtr  oid;
  CharPtr      ptr;
  CharPtr      str;
  TagListPtr   tlp;
  Char         tmp [16];
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
        dp = (DbtagPtr) list->data.ptrvalue;
        if (dp != NULL && dp->db != NULL && dp->tag != NULL) {
          oid = dp->tag;
          ptr = NULL;
          if (oid->str != NULL) {
            ptr = oid->str;
          } else {
            sprintf (tmp, "%ld", (long) oid->id);
            ptr = tmp;
          }
          len = StringLen (dp->db) + StringLen (ptr);
          str = MemNew (len + 4);
          if (str != NULL) {
            StringCpy (str, dp->db);
            StringCat (str, "\t");
            StringCat (str, ptr);
            StringCat (str, "\n");
          }
          vnp->data.ptrvalue = str;
        }
      }
      list = list->next;
    }
    SendMessageToDialog (tlp->dialog, VIB_MSG_RESET);
    tlp->vnp = head;
    SendMessageToDialog (tlp->dialog, VIB_MSG_REDRAW);
    for (j = 0, vnp = tlp->vnp; vnp != NULL; j++, vnp = vnp->next) {
    }
    if (readOnly) {
      tlp->max = MAX ((Int2) 0, (Int2) (j - tlp->rows));
      CorrectBarMax (tlp->bar, tlp->max);
      CorrectBarPage (tlp->bar, tlp->rows - 1, tlp->rows - 1);
      if (tlp->max > 0) {
        SafeShow (tlp->bar);
      } else {
        SafeHide (tlp->bar);
      }
    } else {
      tlp->max = MAX ((Int2) 0, (Int2) (j - tlp->rows + 1));
      CorrectBarMax (tlp->bar, tlp->max);
      CorrectBarPage (tlp->bar, tlp->rows - 1, tlp->rows - 1);
    }
  }
}

static void DbtagPtrToRODbtagDialog (DialoG d, Pointer data)

{
  DbtagPtrToDbtagDialog (d, data, TRUE);
}

static void DbtagPtrToRWDbtagDialog (DialoG d, Pointer data)

{
  DbtagPtrToDbtagDialog (d, data, FALSE);
}

static Pointer DbtagDialogToDbtagPtr (DialoG d)

{
  Boolean      alldigits;
  Char         ch;
  DbtagPtr     dp;
  ValNodePtr   head;
  Int2         j;
  Boolean      leadingzero;
  Int2         len;
  ValNodePtr   list;
  Boolean      notallzero;
  ObjectIdPtr  oid;
  Boolean      okay;
  CharPtr      ptr;
  CharPtr      str;
  TagListPtr   tlp;
  CharPtr      tmp;
  long int     val;
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
          dp = DbtagNew ();
          list->data.ptrvalue = (Pointer) dp;
          if (dp != NULL) {
            dp->db = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, 0);
            oid = ObjectIdNew ();
            dp->tag = oid;
            if (oid != NULL) {
              tmp = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, 1);
              TrimSpacesAroundString (tmp);
              if (tmp != NULL) {
                leadingzero = FALSE;
                notallzero = FALSE;
                alldigits = TRUE;
                ptr = tmp;
                ch = *ptr;
                if (ch == '0') {
                  leadingzero = TRUE;
                }
                while (ch != '\0') {
                  if (ch == ' ' || ch == '+' || ch == '-') {
                  } else if (! (IS_DIGIT (ch))) {
                    alldigits = FALSE;
                  } else if ('1'<= ch && ch <='9') {
                    notallzero = TRUE;
                  }
                  ptr++;
                  ch = *ptr;
                }
                if (alldigits && (! (leadingzero && notallzero)) && sscanf (tmp, "%ld", &val) == 1) {
                  oid->id = (Int4) val;
                  MemFree (tmp);
                } else {
                  oid->str = tmp;
                }
              } else {
                oid->str = StringSave ("?");
              }
            }
          }
        }
      }
    }
  }
  return (Pointer) head;
}

Uint2 visstring_types [] = {
  TAGLIST_TEXT, TAGLIST_TEXT
};

Uint2 readonlystring_types [] = {
  TAGLIST_PROMPT, TAGLIST_PROMPT
};

Uint2 visstring_widths [] = {
  0, 0, 0, 0
};

extern DialoG CreateVisibleStringDialog (GrouP h, Uint2 rows,
                                         Int2 spacing, Int2 width)

{
  visstring_widths [0] = width;
  return CreateTagListDialog (h, rows, 1, spacing,
                              visstring_types, visstring_widths, NULL,
                              ValNodePtrToVisStringDialog,
                              VisStringDialogToValNodePtr);
}

extern DialoG CreateDbtagDialog (GrouP h, Uint2 rows, Int2 spacing,
                                 Int2 width1, Int2 width2)

{
  DialoG      d;
  TagListPtr  tlp;

  visstring_widths [0] = width1;
  visstring_widths [1] = width2;
  if (GetAppProperty ("ReadOnlyDbTags") == NULL) {
    return CreateTagListDialog (h, rows, 2, spacing,
                                visstring_types, visstring_widths, NULL,
                                DbtagPtrToRWDbtagDialog,
                                DbtagDialogToDbtagPtr);
  } else {
    d = CreateTagListDialog (h, rows, 2, spacing,
                             readonlystring_types, visstring_widths, NULL,
                             DbtagPtrToRODbtagDialog,
                             DbtagDialogToDbtagPtr);
    tlp = (TagListPtr) GetObjectExtra (d);
    if (tlp != NULL) {
      Hide (tlp->bar);
    }
    return d;
  }
}

