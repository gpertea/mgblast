/*   nps2gps.c
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
* File Name:  nps2gps.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   5/12/05
*
* $Revision: 1.7 $
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
#include <seqport.h>
#include <explore.h>
#include <subutil.h>
#include <toasn3.h>
#include <pmfapi.h>

#define NPS2GPSAPP_VER "1.5"

CharPtr NPS2GPSAPPLICATION = NPS2GPSAPP_VER;

typedef struct n2gdata {
  CharPtr  results;
  CharPtr  outfile;
  Boolean  failure;
  Boolean  lock;
} N2GData, PNTR N2GPtr;

typedef struct npsseqs {
  BioseqPtr  nuc;
  BioseqPtr  prot;
} NpsSeqs, PNTR NpsSeqsPtr;

static void FindNucProtSeqs (
  BioseqPtr bsp,
  Pointer userdata
)

{
  NpsSeqsPtr  nsp;

  if (bsp == NULL) return;
  nsp = (NpsSeqsPtr) userdata;
  if (nsp == NULL) return;

  if (ISA_na (bsp->mol)) {
    nsp->nuc = bsp;
  } else if (ISA_aa (bsp->mol)) {
    nsp->prot = bsp;
  }
}

static void LclMakeNucProtCDS (
  BioseqSetPtr bssp,
  Pointer userdata
)

{
  CodeBreakPtr  cbp;
  SeqFeatPtr    cds;
  CdRegionPtr   crp;
  SeqFeatPtr    mrna;
  NpsSeqs       ns;
  Boolean       partial5, partial3;
  SeqFeatPtr    sfp;
  SeqIdPtr      sip;
  SeqLocPtr     slp;
  Int4          start, stop;
  SeqFeatPtr    temp;

  ns.nuc = NULL;
  ns.prot = NULL;
  if (VisitBioseqsInSet (bssp, (Pointer) &ns, FindNucProtSeqs) != 2) return;
  if (ns.nuc == NULL || ns.prot == NULL) return;

  cds = SeqMgrGetCDSgivenProduct (ns.prot, NULL);
  mrna = SeqMgrGetRNAgivenProduct (ns.nuc, NULL);
  if (cds == NULL || mrna == NULL) return;

  CheckSeqLocForPartial (cds->location, &partial5, &partial3);

  start = GetOffsetInLoc (cds->location, mrna->location, SEQLOC_START);
  stop = GetOffsetInLoc (cds->location, mrna->location, SEQLOC_STOP);

  if (start < 0 || start >= ns.nuc->length ||
      stop < 0 || stop >= ns.nuc->length) return;

  sip = SeqIdFindBest (ns.nuc->id, 0);
  if (sip == NULL) return;

  /* copy cds feature fields to paste into new cds feature */
  temp = AsnIoMemCopy (cds,
                       (AsnReadFunc) SeqFeatAsnRead,
                       (AsnWriteFunc) SeqFeatAsnWrite);
  if (temp == NULL) return;

  sfp = CreateNewFeatureOnBioseq (ns.nuc, SEQFEAT_CDREGION, NULL);
  if (sfp == NULL) return;

  sfp->location = SeqLocFree (sfp->location);
  if (StringISearch (cds->except_text, "ribosomal slippage") == NULL &&
      StringISearch (cds->except_text, "ribosome slippage") == NULL &&
      StringISearch (cds->except_text, "trans splicing") == NULL &&
      StringISearch (cds->except_text, "trans-splicing") == NULL &&
      StringISearch (cds->except_text, "artificial frameshift") == NULL) {
    sfp->location = AddIntervalToLocation (NULL, sip, start, stop, partial5, partial3);
  } else {
    slp = SeqLocFindNext (cds->location, NULL);
    while (slp != NULL) {
      start = GetOffsetInLoc (slp, mrna->location, SEQLOC_START);
      stop = GetOffsetInLoc (slp, mrna->location, SEQLOC_STOP);
      sfp->location = AddIntervalToLocation (sfp->location, sip, start,
                                             stop, partial5, partial3);
      slp = SeqLocFindNext (cds->location, slp);
    }
    sfp->location = SeqLocMergeEx (ns.nuc, sfp->location, NULL, FALSE, TRUE, FALSE, FALSE);
  }
  SetSeqFeatProduct (sfp, ns.prot);

  /* paste fields from temp copy of original cds */
  crp = (CdRegionPtr) temp->data.value.ptrvalue;
  sfp->data.value.ptrvalue = (Pointer) crp;

  sfp->partial = temp->partial;
  sfp->excpt = temp->excpt;
  sfp->comment = temp->comment;
  sfp->qual = temp->qual;
  sfp->title = temp->title;
  sfp->ext = temp->ext;
  sfp->cit = temp->cit;
  sfp->exp_ev = temp->exp_ev;
  sfp->xref = temp->xref;
  sfp->dbxref = temp->dbxref;
  sfp->pseudo = temp->pseudo;
  sfp->except_text = temp->except_text;

  MemFree (temp); /* do not SeqFeatFree */

  /* update code break locations */
  for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
    CheckSeqLocForPartial (cbp->loc, &partial5, &partial3);
    start = GetOffsetInLoc (cbp->loc, mrna->location, SEQLOC_START);
    stop = GetOffsetInLoc (cbp->loc, mrna->location, SEQLOC_STOP);
    if (start < 0 || start >= ns.nuc->length ||
        stop < 0 || stop >= ns.nuc->length) continue;
	cbp->loc = SeqLocFree (cbp->loc);
	cbp->loc = AddIntervalToLocation (NULL, sip, start, stop, partial5, partial3);;
  }
}

/* copy gene from contig onto nuc-prot, single interval on cdna bioseq */

static void LclCopyGene (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  BioseqPtr          bsp;
  SeqMgrFeatContext  gcontext;
  SeqFeatPtr         gene, copy, temp;
  GeneRefPtr         grp, xref;
  Boolean            partial5, partial3;

  /* input mrna features are multi-interval on contig */

  if (sfp->data.choice != SEQFEAT_RNA) return;

  /* find cdna product of mrna */

  bsp = BioseqFindFromSeqLoc (sfp->product);
  if (bsp == NULL) return;

  /* check for gene xref */

  xref = SeqMgrGetGeneXref (sfp);
  if (xref != NULL) {
    if (SeqMgrGeneIsSuppressed (xref)) return;

    /* copy gene xref for new gene feature */

    grp = AsnIoMemCopy (xref,
                        (AsnReadFunc) GeneRefAsnRead,
                        (AsnWriteFunc) GeneRefAsnWrite);
    if (grp == NULL) return;

    /* make new gene feature on full-length of cdna */

    copy = CreateNewFeatureOnBioseq (bsp, SEQFEAT_GENE, NULL);
    if (copy == NULL) return;

    copy->data.value.ptrvalue = grp;
    return;
  }

  /* overlapping gene should be single interval on contig */

  gene = SeqMgrGetOverlappingGene (sfp->location, &gcontext);
  if (gene == NULL) return;

  CheckSeqLocForPartial (gene->location, &partial5, &partial3);

  /* copy gene feature fields to paste into new gene feature */

  temp = AsnIoMemCopy (gene,
                       (AsnReadFunc) SeqFeatAsnRead,
                       (AsnWriteFunc) SeqFeatAsnWrite);
  if (temp == NULL) return;

  /* make new gene feature on full-length of cdna */

  copy = CreateNewFeatureOnBioseq (bsp, SEQFEAT_GENE, NULL);
  if (copy == NULL) {
    SeqFeatFree (temp);
    return;
  }

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

  SetSeqLocPartial (copy->location, partial5, partial3);

  SeqLocFree (temp->location);
  MemFree (temp); /* do not SeqFeatFree */
}

static void LclAddMrnaTitles (
  SeqLocPtr slp,
  Pointer userdata
)

{
  BioseqPtr          bsp;
  SeqMgrFeatContext  ccontext;
  CharPtr            cdslabel = NULL;
  SeqMgrFeatContext  gcontext;
  CharPtr            genelabel = NULL;
  size_t             len;
  CharPtr            organism;
  SeqFeatPtr         sfp;
  CharPtr            str;

  if (slp == NULL) return;
  bsp = BioseqFindFromSeqLoc (slp);
  if (bsp == NULL) return;
  if (! ISA_na (bsp->mol)) return;
  organism = (CharPtr) userdata;
  if (BioseqGetTitle (bsp) != NULL) return;
  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, 0, &gcontext);
  if (sfp != NULL) {
    genelabel = gcontext.label;
    if (StringHasNoText (genelabel)) {
      genelabel = NULL;
    }
  }
  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &ccontext);
  if (sfp != NULL) {
    cdslabel = ccontext.label;
    if (StringHasNoText (cdslabel)) {
      cdslabel = NULL;
    }
  }
  len = StringLen (organism) + StringLen (genelabel) + StringLen (cdslabel) +
        StringLen (" mRNA, complete cds.") + 10;
  str = (CharPtr) MemNew (len * sizeof (Char));
  if (str == NULL) return;
  str [0] = '\0';

  if (! StringHasNoText (organism)) {
    StringCat (str, organism);
  }
  if (cdslabel != NULL) {
    StringCat (str, " ");
    StringCat (str, cdslabel);
  }
  if (genelabel != NULL) {
      StringCat (str, " (");
      StringCat (str, genelabel);
      StringCat (str, ")");
  }
  if (cdslabel != NULL && genelabel != NULL) {
    if (ccontext.partialL || ccontext.partialR) {
      StringCat (str, " mRNA, partial cds.");
    } else {
      StringCat (str, " mRNA, complete cds.");
    }
  } else if (genelabel != NULL) {
    StringCat (str, " mRNA.");
  }
  SeqDescrAddPointer (&(bsp->descr), Seq_descr_title, (Pointer) str);
}

static SeqIdPtr MakeIdFromLocusTag (SeqFeatPtr mrna)

{
  Char        buf [64], suffix [8];
  SeqFeatPtr  gene;
  GeneRefPtr  grp;

  if (mrna == NULL) return NULL;
  suffix [0] = '\0';
  grp = SeqMgrGetGeneXref (mrna);
  if (grp != NULL) {
    if (SeqMgrGeneIsSuppressed (grp)) return NULL;
  }
  if (grp == NULL) {
    gene = SeqMgrGetOverlappingGene (mrna->location, NULL);
    if (gene != NULL) {
      grp = (GeneRefPtr) gene->data.value.ptrvalue;
      if (gene->id.choice == 1) {
        (gene->id.value.intvalue)++;
        sprintf (suffix, "_%ld", (long) gene->id.value.intvalue);
      }
    }
  }
  if (grp != NULL) {
    if (StringDoesHaveText (grp->locus_tag)) {
      /* StringCpy (buf, "lcl|"); */
      StringCpy (buf, "gnl|MTRACK|");
      StringCat (buf, grp->locus_tag);
      StringCat (buf, suffix);
      return MakeSeqID (buf);
    }
  }
  return NULL;
}

static void InstantiateMrnaIntoProt (SeqFeatPtr cds, SeqFeatPtr mrna, Int2Ptr ctrp, Boolean ambig)

{
  ByteStorePtr  bs;
  Int4          len;
  BioseqPtr     mbsp, pbsp;
  MolInfoPtr    mip;
  SeqEntryPtr   msep, psep;
  Boolean       partial5, partial3;
  CharPtr       rnaseq;
  ValNodePtr    vnp;

  if (cds == NULL || mrna == NULL) return;
  if (cds->product == NULL || mrna->product != NULL) return;

  pbsp = BioseqFindFromSeqLoc (cds->product);
  if (pbsp  == NULL) return;
  psep = SeqMgrGetSeqEntryForData (pbsp);
  if (psep == NULL) return;

  rnaseq = GetSequenceByFeature (mrna);
  if (rnaseq == NULL) return;
  len = (Int4) StringLen (rnaseq);

  bs = BSNew (len + 2);
  if (bs == NULL) return;
  BSWrite (bs, (VoidPtr) rnaseq, len);
  MemFree (rnaseq);

  mbsp = BioseqNew ();
  if (mbsp == NULL) return;
  mbsp->repr = Seq_repr_raw;
  mbsp->mol = Seq_mol_rna;
  mbsp->seq_data_type = Seq_code_iupacna;
  mbsp->seq_data = bs;
  mbsp->length = BSLen (bs);
  BioseqPack (mbsp);

  /*
  if (! ambig) {
    mbsp->id = MakeIdFromLocusTag (mrna);
  }
  */
  /* now adds _# suffix to general Seq-id if ambiguous */
  mbsp->id = MakeIdFromLocusTag (mrna);
  if (mbsp->id == NULL) {
    mbsp->id = MakeNewProteinSeqIdEx (mrna->location, NULL, NULL, ctrp);
  }
  CheckSeqLocForPartial (mrna->location, &partial5, &partial3);
  SeqMgrAddToBioseqIndex (mbsp);

  msep = SeqEntryNew ();
  if (msep == NULL) return;
  msep->choice = 1;
  msep->data.ptrvalue = (Pointer) mbsp;
  SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) mbsp, msep);

  mip = MolInfoNew ();
  if (mip == NULL) return;
  mip->biomol = MOLECULE_TYPE_MRNA;
  if (partial5 && partial3) {
    mip->completeness = 5;
  } else if (partial5) {
    mip->completeness = 3;
  } else if (partial3) {
    mip->completeness = 4;
  }
  vnp = CreateNewDescriptor (msep, Seq_descr_molinfo);
  if (vnp != NULL) {
    vnp->data.ptrvalue = (Pointer) mip;
  }

  SetSeqFeatProduct (mrna, mbsp);

  /* add cDNA to protein, promoting to nuc-prot set component of genomic product set */

  AddSeqEntryToSeqEntry (psep, msep, FALSE);
}

static void RemoveFeatIDs (SeqFeatPtr sfp, Pointer userdata)

{
  if (sfp == NULL) return;
  if (sfp->id.choice == 1) {
    sfp->id.choice = 0;
    sfp->id.value.intvalue = 0;
  }
}

typedef struct loopdata {
  Boolean     ambig;
  Int2        count;
  SeqFeatPtr  mrna;
} LoopData, PNTR LoopDataPtr;

static Boolean LIBCALLBACK FindSingleMrnaProc (
  SeqFeatPtr sfp,
  SeqMgrFeatContextPtr context
)

{
  LoopDataPtr  ldp;

  ldp = (LoopDataPtr) context->userdata;

  if (sfp->id.choice == 0) {
    (ldp->count)++;
    ldp->mrna = sfp;
  } else {
    ldp->ambig = TRUE;
  }

  return TRUE;
}

static void LoopThroughCDSs (BioseqPtr bsp, Int2Ptr ctrp, N2GPtr ngp)

{
  Int2               count;
  SeqMgrFeatContext  fcontext;
  Boolean            goOn;
  LoopData           ld;
  SeqFeatPtr         sfp;
  CharPtr            str;

  /* loop through CDS features, finding single unused mRNA partner */

  goOn = TRUE;
  while (goOn) {
    goOn = FALSE;
    sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
    while (sfp != NULL) {
      if (sfp->id.choice == 0) {
        ld.ambig = FALSE;
        ld.count = 0;
        ld.mrna = NULL;
        count = SeqMgrGetAllOverlappingFeatures (sfp->location, FEATDEF_mRNA, NULL, 0,
                                                 CHECK_INTERVALS, (Pointer) &ld, FindSingleMrnaProc);
        if (ld.count == 1 && ld.mrna != NULL) {
          InstantiateMrnaIntoProt (sfp, ld.mrna, ctrp, ld.ambig);
          sfp->id.choice = 1;
          ld.mrna->id.choice = 1;
          goOn = TRUE;
        } else {
          str = fcontext.label;
          if (StringHasNoText (str)) {
            str = "?";
          }
          Message (MSG_POSTERR, " CDS '%s' has %d overlapping mRNA", str, (int) ld.count);
        }
      }
      sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext);
    }
  }

  /* if any CDSs were not promoted, record failure */

  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
  while (sfp != NULL) {
    if (sfp->id.choice == 0) {
      if (ngp != NULL) {
        ngp->failure = TRUE;
      }
      str = fcontext.label;
      if (StringHasNoText (str)) {
        str = "?";
      }
      Message (MSG_POSTERR, "Failed to match CDS '%s' to mRNA", str);
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext);
  }
}

static void NPStoGPS (Uint2 entityID, CharPtr filename, N2GPtr ngp, SeqDescrPtr descr)

{
  BioSourcePtr       biop;
  BioseqPtr          bsp;
  BioseqSetPtr       bssp;
  Int2               ctr = 1;
  SeqMgrFeatContext  gcontext, mcontext;
  SeqFeatPtr         gene, mrna, sfp, lastsfp;
  GeneRefPtr         grp;
  SeqEntryPtr        old, top;
  CharPtr            organism;
  OrgRefPtr          orp;
  Uint2              parenttype;
  Pointer            parentptr;
  SeqAnnotPtr        sap, annot, lastsap;
  SeqDescrPtr        sdp;

  top = GetTopSeqEntryForEntityID (entityID);
  if (top == NULL || top->choice != 2) return;
  bssp = (BioseqSetPtr) top->data.ptrvalue;
  if (bssp == NULL || bssp->_class != BioseqseqSet_class_nuc_prot) return;

  if (StringHasNoText (filename)) {
    filename = "?";
  }

  bsp = FindNucBioseq (top);
  if (bsp == NULL) {
    Message (MSG_OK, "Unable to find nucleotide Bioseq in %s", filename);
    if (ngp != NULL) {
      ngp->failure = TRUE;
    }
    return;
  }

  GetSeqEntryParent (top, &parentptr, &parenttype);

  bssp->_class = BioseqseqSet_class_gen_prod_set;

  /* move feature table from nuc-prot set to nucleotide bioseq */

  sap = bssp->annot;
  bssp->annot = NULL;
  annot = bsp->annot;
  if (annot == NULL) {
    bsp->annot = sap;
  } else if (sap->next == NULL && annot->next == NULL &&
             sap->type == 1 && annot->type == 1 &&
             sap->data != NULL && annot->data != NULL) {
    lastsfp = (SeqFeatPtr) annot->data;
    while (lastsfp->next != NULL) {
      lastsfp = lastsfp->next;
    }
    sfp = (SeqFeatPtr) sap->data;
    lastsfp->next = sfp;
    sap->data = NULL;
    SeqAnnotFree (sap);
  } else {
    lastsap = bsp->annot;
    while (lastsap->next != NULL) {
      lastsap = lastsap->next;
    }
    lastsap->next = sap;
  }

  VisitFeaturesInSep (top, NULL, RemoveFeatIDs);

  old = SeqEntrySetScope (top);

  ctr = (Int2) VisitBioseqsInSep (top, NULL, NULL) + 1;

  /* count mRNAs per overlapping gene */

  mrna = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_mRNA, &mcontext);
  while (mrna != NULL) {
    grp = SeqMgrGetGeneXref (mrna);
    if (grp == NULL) {
      gene = SeqMgrGetOverlappingGene (mrna->location, NULL);
      if (gene != NULL) {
        if (gene->id.choice == 0) {
          gene->id.choice = 1;
          gene->id.value.intvalue = 1;
        } else if (gene->id.choice == 1) {
          (gene->id.value.intvalue)++;
        }
      }
    }
    mrna = SeqMgrGetNextFeature (bsp, mrna, 0, FEATDEF_mRNA, &mcontext);
  }

  /* only leave genes with multiple mRNAs marked */

  gene = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, 0, &gcontext);
  while (gene != NULL) {
    if (gene->id.choice == 1) {
      if (gene->id.value.intvalue == 1) {
        /* if count was 1, clear id */
        gene->id.value.intvalue = 0;
        gene->id.choice = 0;
      } else {
        /* if count was > 1, just reset count, do not clear id */
        gene->id.value.intvalue = 0;
      }
    }
    gene = SeqMgrGetNextFeature (bsp, gene, SEQFEAT_GENE, 0, &gcontext);
  }

  /*
  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
  while (sfp != NULL) {
    mrna = SeqMgrGetOverlappingFeature (sfp->location, FEATDEF_mRNA, NULL, 0,
                                        NULL, CHECK_INTERVALS, &mcontext);
    if (mrna != NULL) {
      InstantiateMrnaIntoProt (sfp, mrna, &ctr);
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext);
  }
  */

  /* make cDNA bioseq from mRNA feature, package with protein product */

  LoopThroughCDSs (bsp, &ctr, ngp);

  SeqMgrLinkSeqEntry (top, parenttype, parentptr);

  SeqEntrySetScope (old);

  SeqMgrClearFeatureIndexes (bssp->idx.entityID, NULL);

  VisitFeaturesInSep (top, NULL, RemoveFeatIDs);

  /* need to reindex to get mRNA and CDS features from cDNA and protein */
  SeqMgrIndexFeatures (bssp->idx.entityID, NULL);
  VisitSetsInSet (bssp, NULL, LclMakeNucProtCDS);

  /* need to reindex before copying genes, instantiating protein titles */
  SeqMgrIndexFeatures (bssp->idx.entityID, NULL);

  VisitFeaturesInSep (top, NULL, LclCopyGene);

  /* need to reindex before instantiating mRNA titles */
  SeqMgrIndexFeatures (bssp->idx.entityID, NULL);

  organism = NULL;
  for (sdp = descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice != Seq_descr_source) continue;
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    if (biop == NULL) continue;
    orp = biop->org;
    if (orp == NULL) continue;
    if (! StringHasNoText (orp->taxname)) {
      organism = orp->taxname;
    }
  }

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_mRNA, &mcontext);
  while (sfp != NULL) {
    LclAddMrnaTitles (sfp->product, organism);
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_mRNA, &mcontext);
  }

  SeqMgrClearFeatureIndexes (bssp->idx.entityID, NULL);

  move_cds (top);
}

static void ProcessOneRecord (
  CharPtr filename,
  Pointer userdata
)

{
  AsnIoPtr      aip;
  BioseqPtr     bsp;
  ValNodePtr    bsplist;
  BioseqSetPtr  bssp;
  Pointer       dataptr = NULL;
  Uint2         datatype, entityID = 0;
  SeqDescrPtr   descr1, descr2;
  Char          file [FILENAME_MAX], id [42], path [PATH_MAX];
  FILE          *fp;
  N2GPtr        ngp;
  SeqEntryPtr   nsep, sep;
  CharPtr       ptr, str;
  SeqIdPtr      sip;

  if (StringHasNoText (filename)) return;
  ngp = (N2GPtr) userdata;
  if (ngp == NULL) return;

  fp = FileOpen (filename, "r");
  if (fp == NULL) {
    Message (MSG_POSTERR, "Failed to open '%s'", filename);
    return;
  }

  dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, FALSE, FALSE, FALSE);

  FileClose (fp);

  entityID = ObjMgrRegister (datatype, dataptr);

  if (entityID < 1 || dataptr == NULL) {
    Message (MSG_POSTERR, "Data read failed for input file '%s'", filename);
    return;
  }

  if (datatype == OBJ_SEQSUB || datatype == OBJ_SEQENTRY ||
        datatype == OBJ_BIOSEQ || datatype == OBJ_BIOSEQSET) {

    sep = GetTopSeqEntryForEntityID (entityID);

    if (sep == NULL) {
      sep = SeqEntryNew ();
      if (sep != NULL) {
        if (datatype == OBJ_BIOSEQ) {
          bsp = (BioseqPtr) dataptr;
          sep->choice = 1;
          sep->data.ptrvalue = bsp;
          SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);
        } else if (datatype == OBJ_BIOSEQSET) {
          bssp = (BioseqSetPtr) dataptr;
          sep->choice = 2;
          sep->data.ptrvalue = bssp;
          SeqMgrSeqEntry (SM_BIOSEQSET, (Pointer) bssp, sep);
        } else {
          sep = SeqEntryFree (sep);
        }
      }
      sep = GetTopSeqEntryForEntityID (entityID);
    }

    if (sep != NULL) {
      bsplist = NULL;
      if (ngp->lock) {
        bsplist = LockFarComponents (sep);
      }
    
      str = StringRChr (filename, DIRDELIMCHR);
      if (str != NULL) {
        str++;
      } else {
        str = filename;
      }
      StringCpy (file, str);

      bsp = FindNucBioseq (sep);
      id [0] = '\0';
      if (bsp != NULL) {
        sip = SeqIdFindWorst (bsp->id);
        SeqIdWrite (sip, id, PRINTID_REPORT, sizeof (id));
        Message (MSG_POST, "Processing %s (%s)", id, file);
      }

      SeqMgrIndexFeatures (entityID, NULL);

      /* remove pubs and biosources from nuc-prot set and nucleotide bioseq */

      nsep = FindNucSeqEntry (sep);
      descr1 = ExtractBioSourceAndPubs (sep);
      descr2 = ExtractBioSourceAndPubs (nsep);
      if (descr1 == NULL) {
        descr1 = descr2;
      } else if (descr2 != NULL) {
        ValNodeLink (&descr1, descr2);
      }

      /* convert from nuc-prot set to genomic-product set */

      NPStoGPS (entityID, file, ngp, descr1);

      /* put pubs and biosources onto genomic-product set */

      ReplaceBioSourceAndPubs (sep, descr1);

      if (ngp->failure) {
        Message (MSG_POST, "Not saving %s", file);
      } else if (StringDoesHaveText (ngp->outfile)) {
        aip = AsnIoOpen (ngp->outfile, "w");
        if (aip != NULL) {
          SeqEntryAsnWrite (sep, aip, NULL);
          AsnIoClose (aip);
        }
      } else {
        str = StringRChr (filename, DIRDELIMCHR);
        if (str != NULL) {
          *str = '\0';
          str++;
        } else {
          str = filename;
        }
        if (StringDoesHaveText (ngp->results)) {
          StringCpy (path, ngp->results);
        } else {
          StringCpy (path, filename);
        }
        if (str != NULL) {
          ptr = StringRChr (str, '.');
          if (ptr != NULL) {
            *ptr = '\0';
          }
          StringCpy (file, str);
          StringCat (file, ".gps");
          FileBuildPath (path, NULL, file);
          aip = AsnIoOpen (path, "w");
          if (aip != NULL) {
            SeqEntryAsnWrite (sep, aip, NULL);
            AsnIoClose (aip);
          }
        }
      }

      bsplist = UnlockFarComponents (bsplist);
    }

  } else {

    Message (MSG_POSTERR, "Datatype %d not recognized", (int) datatype);
  }

  ObjMgrFree (datatype, dataptr);
}

/* Args structure contains command-line arguments */

#define p_argInputPath   0
#define r_argOutputPath  1
#define i_argInputFile   2
#define o_argOutputFile  3
#define f_argFilter      4
#define x_argSuffix      5
#define R_argRemote      6
#define L_argLockFar     7


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
  {"File Selection Suffix", ".ent", NULL, NULL,
    TRUE, 'x', ARG_STRING, 0.0, 0, NULL},
  {"Remote Fetching from ID", "F", NULL, NULL,
    TRUE, 'R', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Lock Components in Advance", "F", NULL, NULL,
    TRUE, 'L', ARG_BOOLEAN, 0.0, 0, NULL},
};

Int2 Main (void)

{
  Char     app [64];
  CharPtr  directory, filter, infile, outfile, results, suffix;
  N2GData  ngd;
  N2GPtr   ngp;
  Boolean  remote;

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

  sprintf (app, "nps2gps %s", NPS2GPSAPPLICATION);
  if (! GetArgs (app, sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  MemSet ((Pointer) &ngd, 0, sizeof (N2GData));
  ngp = &ngd;
  ngd.failure = FALSE;
  ngd.lock = (Boolean) myargs [L_argLockFar].intvalue;

  directory = (CharPtr) myargs [p_argInputPath].strvalue;
  results = (CharPtr) myargs [r_argOutputPath].strvalue;
  if (StringHasNoText (results)) {
    results = NULL;
  }
  infile = (CharPtr) myargs [i_argInputFile].strvalue;
  outfile = (CharPtr) myargs [o_argOutputFile].strvalue;
  filter = (CharPtr) myargs [f_argFilter].strvalue;
  suffix = (CharPtr) myargs [x_argSuffix].strvalue;

  remote = (Boolean) myargs [R_argRemote].intvalue;

  /* register fetch function */

  if (remote) {
    PubSeqFetchEnable ();
  }

  /* process input file */

  if (StringDoesHaveText (directory)) {

    ngd.results = results;

    DirExplore (directory, filter, suffix, TRUE, ProcessOneRecord, (Pointer) &ngd);

  } else if (StringDoesHaveText (infile) && StringDoesHaveText (outfile)) {

    ngd.outfile = outfile;

    ProcessOneRecord (infile, (Pointer) &ngd);
  }

  /* close fetch function */

  if (remote) {
    PubSeqFetchDisable ();
  }

  if (ngd.failure) return 1;

  return 0;
}

