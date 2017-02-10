/*   asn2gb.c
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
* File Name:  asn2gb.c
*
* Author:  Karl Sirotkin, Tom Madden, Tatiana Tatusov, Jonathan Kans
*
* Version Creation Date:   10/21/98
*
* $Revision: 6.88 $
*
* File Description:  New GenBank flatfile generator application
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
#include <objgbseq.h>
#include <objtseq.h>
#include <sequtil.h>
#include <sqnutils.h>
#include <explore.h>
#include <asn2gnbp.h>

#define ASN2GB_APP_VER "3.2"

CharPtr ASN2GB_APPLICATION = ASN2GB_APP_VER;

static void SaveSeqEntry (
  SeqEntryPtr sep,
  CharPtr filename
)

{
  AsnIoPtr  aop;

  if (sep == NULL) return;
  aop = AsnIoOpen (filename, "w");
  if (aop != NULL) {
    SeqEntryAsnWrite (sep, aop, NULL);
  }
  AsnIoClose (aop);
}

static void SaveAsn2gnbk (
  SeqEntryPtr sep,
  CharPtr filename,
  FmtType format,
  ModType mode,
  StlType style,
  FlgType flags,
  LckType locks,
  CstType custom
)

{
  FILE  *fp;

  if (sep == NULL) return;
  fp = FileOpen (filename, "w");
  if (fp != NULL) {
    SeqEntryToGnbk (sep, NULL, format, mode, style, flags, locks, custom, NULL, fp);
  }
  FileClose (fp);
}

static void GetFirstGoodBioseq (
  BioseqPtr bsp,
  Pointer userdata
)

{
  BioseqPtr PNTR bspp;

  bspp = (BioseqPtr PNTR) userdata;
  if (*bspp != NULL) return;
  *bspp = bsp;
}

NLM_EXTERN void AsnPrintNewLine PROTO((AsnIoPtr aip));

static void SaveTinySeqs (
  BioseqPtr bsp,
  Pointer userdata
)

{
  AsnIoPtr  aip;
  TSeqPtr   tsp;

  if (bsp == NULL) return;
  aip = (AsnIoPtr) userdata;

  tsp = BioseqToTSeq (bsp);
  if (tsp == NULL) return;

  TSeqAsnWrite (tsp, aip, NULL);
  AsnPrintNewLine (aip);
  AsnIoFlush (aip);

  TSeqFree (tsp);
}

static void SaveTinyStreams (
  BioseqPtr bsp,
  Pointer userdata
)

{
  AsnIoPtr  aip;

  if (bsp == NULL) return;
  aip = (AsnIoPtr) userdata;

  BioseqAsnWriteAsTSeq (bsp, aip, NULL);
  AsnPrintNewLine (aip);
  AsnIoFlush (aip);
}

static Int2 HandleSingleRecord (
  CharPtr inputFile,
  CharPtr outputFile,
  FmtType format,
  FmtType altformat,
  ModType mode,
  StlType style,
  FlgType flags,
  LckType locks,
  CstType custom,
  XtraPtr extra,
  Int2 type,
  Boolean binary,
  Boolean compressed,
  Int4 from,
  Int4 to,
  Uint1 strand,
  Int4 itemID,
  Boolean do_tiny_seq,
  Boolean do_fasta_stream
)

{
  AsnIoPtr      aip;
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Pointer       dataptr = NULL;
  Uint2         datatype = 0;
  Uint2         entityID;
  FILE          *fp;
  FILE          *ofp = NULL;
  ObjMgrPtr     omp;
  SeqEntryPtr   sep;
  SeqFeatPtr    sfp;
  SeqInt        sint;
  SeqLocPtr     slp = NULL;
  ValNode       vn;

  if (type == 1) {
    fp = FileOpen (inputFile, "r");
    if (fp == NULL) {
      Message (MSG_POSTERR, "FileOpen failed for input file '%s'", inputFile);
      return 1;
    }

    dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, FALSE, FALSE, FALSE);

    FileClose (fp);

    entityID = ObjMgrRegister (datatype, dataptr);

  } else if (type >= 2 && type <= 5) {
    aip = AsnIoOpen (inputFile, binary? "rb" : "r");
    if (aip == NULL) {
      Message (MSG_POSTERR, "AsnIoOpen failed for input file '%s'", inputFile);
      return 1;
    }

    switch (type) {
      case 2 :
        dataptr = (Pointer) SeqEntryAsnRead (aip, NULL);
        datatype = OBJ_SEQENTRY;
        break;
      case 3 :
        dataptr = (Pointer) BioseqAsnRead (aip, NULL);
        datatype = OBJ_BIOSEQ;
        break;
      case 4 :
        dataptr = (Pointer) BioseqSetAsnRead (aip, NULL);
        datatype = OBJ_BIOSEQSET;
        break;
      case 5 :
        dataptr = (Pointer) SeqSubmitAsnRead (aip, NULL);
        datatype = OBJ_SEQSUB;
        break;
      default :
        break;
    }

    AsnIoClose (aip);

    entityID = ObjMgrRegister (datatype, dataptr);

  } else {
    Message (MSG_POSTERR, "Input format type '%d' unrecognized", (int) type);
    return 1;
  }

  if (dataptr == NULL) {
    Message (MSG_POSTERR, "Data read failed for input file '%s'", inputFile);
    return 1;
  }

  if (datatype == OBJ_SEQSUB || datatype == OBJ_SEQENTRY ||
        datatype == OBJ_BIOSEQ || datatype == OBJ_BIOSEQSET) {

/*
#ifdef WIN_MAC
#if __profile__
    ProfilerSetStatus (TRUE);
#endif
#endif
*/

    entityID = SeqMgrIndexFeatures (entityID, NULL);

/*
#ifdef WIN_MAC
#if __profile__
    ProfilerSetStatus (FALSE);
#endif
#endif
*/

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
      if (extra == NULL) {
        FileRemove (outputFile);
#ifdef WIN_MAC
        FileCreate (outputFile, "TEXT", "ttxt");
#endif
        ofp = FileOpen (outputFile, "w");
      }

      if ((from > 0 && to > 0) || strand == Seq_strand_minus) {
        bsp = NULL;
        if (format == GENPEPT_FMT) {
          VisitSequencesInSep (sep, (Pointer) &bsp, VISIT_PROTS, GetFirstGoodBioseq);
        } else {
          VisitSequencesInSep (sep, (Pointer) &bsp, VISIT_NUCS, GetFirstGoodBioseq);
        }
        if (bsp != NULL) {
          if (strand == Seq_strand_minus && from == 0 && to == 0) {
            from = 1;
            to = bsp->length;
          }
          MemSet ((Pointer) &vn, 0, sizeof (ValNode));
          MemSet ((Pointer) &sint, 0, sizeof (SeqInt));
          sint.from = from - 1;
          sint.to = to - 1;
          sint.strand = strand;
          sint.id = SeqIdFindBest (bsp->id, 0);
          vn.choice = SEQLOC_INT;
          vn.data.ptrvalue = (Pointer) &sint;
          slp = &vn;
        }
      } else if (itemID > 0) {
        sfp = SeqMgrGetDesiredFeature (entityID, 0, itemID, 0, NULL, NULL);
        if (sfp != NULL) {
          slp = sfp->location;
        }
      }

      if (do_tiny_seq) {
        aip = AsnIoNew (ASNIO_TEXT_OUT | ASNIO_XML, ofp, NULL, NULL, NULL);
        VisitBioseqsInSep (sep, (Pointer) aip, SaveTinySeqs);
        AsnIoFree (aip, FALSE);
      } else if (do_fasta_stream) {
        aip = AsnIoNew (ASNIO_TEXT_OUT | ASNIO_XML, ofp, NULL, NULL, NULL);
        VisitBioseqsInSep (sep, (Pointer) aip, SaveTinyStreams);
        AsnIoFree (aip, FALSE);
      } else {
        SeqEntryToGnbk (sep, slp, format, mode, style, flags, locks, custom, extra, ofp);
        if (altformat != 0) {
          SeqEntryToGnbk (sep, slp, altformat, mode, style, flags, locks, custom, extra, ofp);
        }
      }
      if (ofp != NULL) {
        FileClose (ofp);
      }
    }
  } else {
    Message (MSG_POSTERR, "Datatype %d not recognized", (int) datatype);
  }

  omp = ObjMgrGet ();
  ObjMgrReapOne (omp);
  ObjMgrFreeCache (0);
  FreeSeqIdGiCache ();

  ObjMgrFree (datatype, dataptr);

  return 0;
}

typedef struct hasgidata {
  Int4     gi;
  CharPtr  accn;
  Boolean  found;
} HasGiData, PNTR HasGiPtr;

static void LookForGi (
  SeqEntryPtr sep,
  Pointer mydata,
  Int4 index,
  Int2 indent
)

{
  BioseqPtr     bsp;
  HasGiPtr      hgp;
  SeqIdPtr      sip;
  TextSeqIdPtr  tsip;

  if (sep == NULL) return;
  if (! IS_Bioseq (sep)) return;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL) return;
  hgp = (HasGiPtr) mydata;
  if (hgp == NULL) return;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    switch (sip->choice) {
      case SEQID_GI :
        if (sip->data.intvalue == hgp->gi) {
          hgp->found = TRUE;
          return;
        }
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
        tsip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tsip != NULL && hgp->accn!= NULL &&
            StringICmp (tsip->accession, hgp->accn) == 0) {
          hgp->found = TRUE;
          return;
        }
        break;
      default :
        break;
    }
  }
}

static Boolean SeqEntryHasGi (
  SeqEntryPtr sep,
  CharPtr accn
)

{
  HasGiData  hgd;
  long int   val;

  if (sep == NULL || StringHasNoText (accn)) return FALSE;
  MemSet ((Pointer) &hgd, 0, sizeof (HasGiData));
  if (sscanf (accn, "%ld", &val) == 1) {
    hgd.gi = (Int4) val;
  } else {
    hgd.accn = accn;
  }
  hgd.found = FALSE;
  SeqEntryExplore (sep, (Pointer) (&hgd), LookForGi);
  return hgd.found;
}

static void FreeUnpubAffil (
  PubdescPtr pdp,
  Pointer userdata
)

{
  AuthListPtr  alp;
  CitGenPtr    cgp;
  ValNodePtr   vnp;

  if (pdp == NULL) return;
  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice != PUB_Gen) continue;
    cgp = (CitGenPtr) vnp->data.ptrvalue;
    if (cgp == NULL) continue;
    if (cgp->cit != NULL) {
      if (StringNICmp (cgp->cit, "submitted", 8) == 0 ||
                       StringNICmp (cgp->cit, "unpublished", 11) == 0 ||
                       StringNICmp (cgp->cit, "in press", 8) == 0 ||
                       StringNICmp (cgp->cit, "to be published", 15) == 0) {
        cgp->cit = MemFree (cgp->cit);
        cgp->cit = StringSave ("Unpublished");
      }
    }
    alp = cgp->authors;
    if (alp == NULL) continue;
    alp->affil = AffilFree (alp->affil);
  }
}

static void LookForRefSeq (
  BioseqPtr bsp,
  Pointer userdata
)

{
  BoolPtr   hasRefseqP;
  SeqIdPtr  sip;

  hasRefseqP = (BoolPtr) userdata;
  if (*hasRefseqP) return;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_OTHER) {
      *hasRefseqP = TRUE;
      return;
    }
  }
}

static CharPtr fffmt [] = {
  "",
  "genbank",
  "embl",
  "genbank",
  "embl",
  "ftable",
  NULL
};

static CharPtr ffmod [] = {
  "",
  "release",
  "entrez",
  "gbench",
  "dump",
  NULL
};

static CharPtr ffstl [] = {
  "",
  "normal",
  "segment",
  "master",
  "contig",
  NULL
};

static CharPtr ffvew [] = {
  "",
  "nuc",
  "nuc",
  "prot",
  "prot",
  "nuc",
  NULL
};

static void CompareFlatFiles (
  CharPtr path1,
  CharPtr path2,
  CharPtr path3,
  SeqEntryPtr sep,
  FILE* fp,
  FmtType format,
  FmtType altformat,
  ModType mode,
  StlType style,
  FlgType flags,
  LckType locks,
  CstType custom,
  XtraPtr extra,
  Int2 batch,
  CharPtr ffdiff,
  CharPtr asn2flat,
  Boolean useFfdiff
)

{
#ifdef OS_UNIX
  AsnIoPtr     aip;
  Char         arguments [128];
  BioseqPtr    bsp;
  Char         buf [256];
  Char         cmmd [256];
  size_t       ct;
  int          diff;
  FILE         *fpo;
  SeqEntryPtr  fsep;

  if (sep == NULL) return;

  if (batch == 1) {

    SeqEntryToGnbk (sep, NULL, format, mode, style, flags, locks, custom, extra, fp);
    if (altformat != 0) {
      SeqEntryToGnbk (sep, NULL, altformat, mode, style, flags, locks, custom, extra, fp);
    }
    return; /* just make report, nothing to diff */

  } else if (batch == 2) {

#ifdef ASN2GNBK_SUPPRESS_UNPUB_AFFIL
    VisitPubdescsInSep (sep, NULL, FreeUnpubAffil);
#endif

    SaveAsn2gnbk (sep, path1, format, SEQUIN_MODE, style, flags, locks, custom);
    SaveAsn2gnbk (sep, path2, format, RELEASE_MODE, style, flags, locks, custom);

    if (useFfdiff) {
      sprintf (cmmd, "%s -o %s -n %s -d reports", ffdiff, path1, path2);
      system (cmmd);

      sprintf (cmmd, "rm %s; rm %s", path1, path2);
      system (cmmd);
    } else {
      sprintf (cmmd, "sort %s | uniq -c > %s.suc; rm %s", path1, path1, path1);
      system (cmmd);

      sprintf (cmmd, "sort %s | uniq -c > %s.suc; rm %s", path2, path2, path2);
      system (cmmd);

      sprintf (cmmd, "diff %s.suc %s.suc > %s", path1, path2, path3);
      system (cmmd);

      sprintf (cmmd, "cat %s", path3);
      fpo = popen (cmmd, "r");
      if (fpo != NULL) {
        while ((ct = fread (buf, 1, sizeof (buf), fpo)) > 0) {
          fwrite (buf, 1, ct, fp);
          fflush (fp);
        }
        pclose (fpo);
      }

      sprintf (cmmd, "rm %s.suc; rm %s.suc", path1, path2);
      system (cmmd);
    }

  } else if (batch == 3) {

    SaveAsn2gnbk (sep, path1, format, mode, style, flags, locks, custom);

    aip = AsnIoOpen (path3, "w");
    if (aip == NULL) return;

    SeqEntryAsnWrite (sep, aip, NULL);
    AsnIoClose (aip);

    fsep = FindNthBioseq (sep, 1);
    if (fsep == NULL || fsep->choice != 1) return;
    bsp = (BioseqPtr) fsep->data.ptrvalue;
    if (bsp == NULL) return;
    SeqIdWrite (bsp->id, buf, PRINTID_FASTA_LONG, sizeof (buf));

    arguments [0] = '\0';
    sprintf (arguments, "-format %s -mode %s -style %s -view %s",
             fffmt [(int) format], ffmod [(int) mode], ffstl [(int) style], ffvew [(int) format]);

    sprintf (cmmd, "%s %s -i %s -o %s", asn2flat, arguments, path3, path2);
    system (cmmd);

    if (useFfdiff) {
      sprintf (cmmd, "%s -o %s -n %s -d reports", ffdiff, path1, path2);
      system (cmmd);

      sprintf (cmmd, "rm %s; rm %s", path1, path2);
      system (cmmd);
    } else {
      sprintf (cmmd, "sort %s | uniq -c > %s.suc; rm %s", path1, path1, path1);
      system (cmmd);

      sprintf (cmmd, "sort %s | uniq -c > %s.suc; rm %s", path2, path2, path2);
      system (cmmd);

      sprintf (cmmd, "diff %s.suc %s.suc > %s", path1, path2, path3);
      system (cmmd);

      sprintf (cmmd, "cat %s", path3);
      fpo = popen (cmmd, "r");
      if (fpo != NULL) {
        while ((ct = fread (buf, 1, sizeof (buf), fpo)) > 0) {
          fwrite (buf, 1, ct, fp);
          fflush (fp);
        }
        pclose (fpo);
      }

      sprintf (cmmd, "rm %s.suc; rm %s.suc", path1, path2);
      system (cmmd);
    }

  } else if (batch == 4) {

    aip = AsnIoOpen (path3, "w");
    if (aip == NULL) return;

    SeqEntryAsnWrite (sep, aip, NULL);
    AsnIoClose (aip);

    fsep = FindNthBioseq (sep, 1);
    if (fsep == NULL || fsep->choice != 1) return;
    bsp = (BioseqPtr) fsep->data.ptrvalue;
    if (bsp == NULL) return;
    SeqIdWrite (bsp->id, buf, PRINTID_FASTA_LONG, sizeof (buf));

    arguments [0] = '\0';
    sprintf (arguments, "-format %s -mode %s -style %s -view %s -nocleanup",
             fffmt [(int) format], ffmod [(int) mode], ffstl [(int) style], ffvew [(int) format]);

    sprintf (cmmd, "%s %s -i %s -o %s", asn2flat, arguments, path3, path1);
    system (cmmd);

    arguments [0] = '\0';
    sprintf (arguments, "-format %s -mode %s -style %s -view %s",
             fffmt [(int) format], ffmod [(int) mode], ffstl [(int) style], ffvew [(int) format]);

    sprintf (cmmd, "%s %s -i %s -o %s", asn2flat, arguments, path3, path2);
    system (cmmd);

    if (useFfdiff) {
      sprintf (cmmd, "%s -o %s -n %s -d reports", ffdiff, path1, path2);
      system (cmmd);

      sprintf (cmmd, "rm %s; rm %s", path1, path2);
      system (cmmd);
    } else {
      sprintf (cmmd, "sort %s | uniq -c > %s.suc; rm %s", path1, path1, path1);
      system (cmmd);

      sprintf (cmmd, "sort %s | uniq -c > %s.suc; rm %s", path2, path2, path2);
      system (cmmd);

      sprintf (cmmd, "diff %s.suc %s.suc > %s", path1, path2, path3);
      system (cmmd);

      sprintf (cmmd, "cat %s", path3);
      fpo = popen (cmmd, "r");
      if (fpo != NULL) {
        while ((ct = fread (buf, 1, sizeof (buf), fpo)) > 0) {
          fwrite (buf, 1, ct, fp);
          fflush (fp);
        }
        pclose (fpo);
      }

      sprintf (cmmd, "rm %s.suc; rm %s.suc", path1, path2);
      system (cmmd);
    }

  } else if (batch == 5) {

    aip = AsnIoOpen (path3, "w");
    if (aip == NULL) return;

    SeqEntryAsnWrite (sep, aip, NULL);
    AsnIoClose (aip);

    if (FindNucBioseq (sep) != NULL) {

      sprintf (cmmd, "./oldasn2gb -i %s -o %s", path3, path1);
      system (cmmd);

      sprintf (cmmd, "./newasn2gb -i %s -o %s", path3, path2);
      system (cmmd);

    } else {

      sprintf (cmmd, "./oldasn2gb -f p -i %s -o %s", path3, path1);
      system (cmmd);

      sprintf (cmmd, "./newasn2gb -f p -i %s -o %s", path3, path2);
      system (cmmd);

    }

    sprintf (cmmd, "diff -b %s %s > %s", path1, path2, path3);
    diff = system (cmmd);

    if (diff > 0) {
      sprintf (cmmd, "cat %s", path3);
      fpo = popen (cmmd, "r");
      if (fpo != NULL) {
        fprintf (fp, "\nasn2gb difference in %s\n", buf);
        fflush (fp);
        while ((ct = fread (buf, 1, sizeof (buf), fpo)) > 0) {
          fwrite (buf, 1, ct, fp);
          fflush (fp);
        }
        pclose (fpo);
      }
    }
  }

#else

  SeqEntryToGnbk (sep, NULL, format, mode, style, flags, locks, custom, extra, fp);
  if (altformat != 0) {
    SeqEntryToGnbk (sep, NULL, altformat, mode, style, flags, locks, custom, extra, fp);
  }
#endif
}

static void CheckOrder (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
#ifdef ASN2GNBK_IGNORE_OUT_OF_ORDER
  BoolPtr    bp;
  BioseqPtr  bsp;
#endif
#ifdef ASN2GNBK_REPAIR_OUT_OF_ORDER
  BioseqPtr  bsp;
  SeqLocPtr  gslp;
  Boolean    hasNulls;
  Boolean    noLeft;
  Boolean    noRight;
#endif

  /* ignore order of bonds in heterogen features from PDB */

  if (sfp->data.choice == SEQFEAT_HET) return;

#ifdef ASN2GNBK_IGNORE_OUT_OF_ORDER
  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp != NULL && SeqLocBadSortOrder (bsp, sfp->location)) {
    bp = (BoolPtr) userdata;
    *bp = TRUE;
  }
#endif
#ifdef ASN2GNBK_REPAIR_OUT_OF_ORDER
  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp != NULL && SeqLocBadSortOrder (bsp, sfp->location)) {
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
#endif
}

static Int2 HandleMultipleRecords (
  CharPtr inputFile,
  CharPtr outputFile,
  FmtType format,
  FmtType altformat,
  ModType mode,
  StlType style,
  FlgType flags,
  LckType locks,
  CstType custom,
  XtraPtr extra,
  Int2 type,
  Int2 batch,
  Boolean binary,
  Boolean compressed,
  Boolean propOK,
  CharPtr ffdiff,
  CharPtr asn2flat,
  CharPtr accn,
  FILE *logfp
)

{
  AsnIoPtr        aip;
  AsnModulePtr    amp;
  AsnTypePtr      atp, atp_bss, atp_desc, atp_sbp, atp_se, atp_ssp;
  Boolean         atp_se_seen = FALSE;
  BioseqPtr       bsp;
  BioseqSetPtr    bssp;
  Char            buf [41];
  Char            cmmd [256];
  CitSubPtr       csp = NULL;
  SeqDescrPtr     descr = NULL;
  FILE            *fp;
  SeqEntryPtr     fsep;
  Boolean         hasgi;
  Boolean         hasRefSeq;
  Char            longest [41];
  Int4            numrecords = 0;
  FILE            *ofp = NULL;
  ObjMgrPtr       omp;
  Boolean         outOfOrder;
  ObjValNode      ovn;
  Char            path1 [PATH_MAX];
  Char            path2 [PATH_MAX];
  Char            path3 [PATH_MAX];
  Pubdesc         pd;
  SubmitBlockPtr  sbp = NULL;
  SeqEntryPtr     sep;
  time_t          starttime, stoptime, worsttime;
  SeqDescrPtr     subcit = NULL;
  FILE            *tfp;
  Boolean         useFfdiff;
  ValNode         vn;
#ifdef OS_UNIX
  CharPtr         gzcatprog;
  int             ret;
  Boolean         usedPopen = FALSE;
#endif

  if (StringHasNoText (inputFile)) return 1;

#ifndef OS_UNIX
  if (compressed) {
    Message (MSG_POSTERR, "Can only decompress on-the-fly on UNIX machines");
    return 1;
  }
#endif

  amp = AsnAllModPtr ();
  if (amp == NULL) {
    Message (MSG_POSTERR, "Unable to load AsnAllModPtr");
    return 1;
  }

  atp_ssp = AsnFind ("Seq-submit");
  if (atp_ssp == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-submit");
    return 1;
  }

  atp_sbp = AsnFind ("Seq-submit.sub");
  if (atp_sbp == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-submit.sub");
    return 1;
  }

  atp_bss = AsnFind ("Bioseq-set");
  if (atp_bss == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq-set");
    return 1;
  }

  atp_desc = AsnFind ("Bioseq-set.descr");
  if (atp_desc == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq-set.descr");
    return 1;
  }

  atp_se = AsnFind ("Bioseq-set.seq-set.E");
  if (atp_se == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq-set.seq-set.E");
    return 1;
  }

#ifdef OS_UNIX
  if (compressed) {
    gzcatprog = getenv ("NCBI_UNCOMPRESS_BINARY");
    if (gzcatprog != NULL) {
      sprintf (cmmd, "%s %s", gzcatprog, inputFile);
    } else {
      ret = system ("gzcat -h >/dev/null 2>&1");
      if (ret == 0) {
        sprintf (cmmd, "gzcat %s", inputFile);
      } else if (ret == -1) {
        Message (MSG_POSTERR, "Unable to fork or exec gzcat in ScanBioseqSetRelease");
        return 1;
      } else {
        ret = system ("zcat -h >/dev/null 2>&1");
        if (ret == 0) {
          sprintf (cmmd, "zcat %s", inputFile);
        } else if (ret == -1) {
          Message (MSG_POSTERR, "Unable to fork or exec zcat in ScanBioseqSetRelease");
          return 1;
        } else {
          Message (MSG_POSTERR, "Unable to find zcat or gzcat in ScanBioseqSetRelease - please edit your PATH environment variable");
          return 1;
        }
      }
    }
    fp = popen (cmmd, /* binary? "rb" : */ "r");
    usedPopen = TRUE;
  } else {
    fp = FileOpen (inputFile, binary? "rb" : "r");
  }
#else
  fp = FileOpen (inputFile, binary? "rb" : "r");
#endif
  if (fp == NULL) {
    Message (MSG_POSTERR, "FileOpen failed for input file '%s'", inputFile);
    return 1;
  }

  aip = AsnIoNew (binary? ASNIO_BIN_IN : ASNIO_TEXT_IN, fp, NULL, NULL, NULL);
  if (aip == NULL) {
    Message (MSG_ERROR, "AsnIoNew failed for input file '%s'", inputFile);
    return 1;
  }

  if ((batch == 1 || batch == 3 || batch == 4 || batch == 5 || format != GENBANK_FMT) && extra == NULL) {
    ofp = FileOpen (outputFile, "w");
    if (ofp == NULL) {
      AsnIoClose (aip);
      Message (MSG_POSTERR, "FileOpen failed for output file '%s'", outputFile);
      return 1;
    }
  }

  TmpNam (path1);
  tfp = FileOpen (path1, "w");
  fprintf (tfp, "\n");
  FileClose (tfp);

  TmpNam (path2);
  tfp = FileOpen (path2, "w");
  fprintf (tfp, "\n");
  FileClose (tfp);

  TmpNam (path3);
  tfp = FileOpen (path3, "w");
  fprintf (tfp, "\n");
  FileClose (tfp);

  if (type == 4) {
    atp = atp_bss;
  } else if (type == 5) {
    atp = atp_ssp;
  } else {
    Message (MSG_ERROR, "Batch processing type not set properly");
    return 1;
  }

  longest [0] = '\0';
  worsttime = 0;

  while ((atp = AsnReadId (aip, amp, atp)) != NULL) {
    if (atp == atp_se) {
      atp_se_seen = TRUE;
      sep = SeqEntryAsnRead (aip, atp);

      /* propagate descriptors from the top-level set */

      if (propOK && descr != NULL && sep != NULL && sep->data.ptrvalue != NULL) {
        if (sep->choice == 1) {
          bsp = (BioseqPtr) sep->data.ptrvalue;
          ValNodeLink (&(bsp->descr),
                       AsnIoMemCopy ((Pointer) descr,
                                     (AsnReadFunc) SeqDescrAsnRead,
                                     (AsnWriteFunc) SeqDescrAsnWrite));
        } else if (sep->choice == 2) {
          bssp = (BioseqSetPtr) sep->data.ptrvalue;
          ValNodeLink (&(bssp->descr),
                       AsnIoMemCopy ((Pointer) descr,
                                     (AsnReadFunc) SeqDescrAsnRead,
                                     (AsnWriteFunc) SeqDescrAsnWrite));
        }
      }

      /* propagate submission citation as descriptor onto each Seq-entry */

      if (subcit != NULL && sep != NULL && sep->data.ptrvalue != NULL) {
        if (sep->choice == 1) {
          bsp = (BioseqPtr) sep->data.ptrvalue;
          ValNodeLink (&(bsp->descr),
                       AsnIoMemCopy ((Pointer) subcit,
                                     (AsnReadFunc) SeqDescrAsnRead,
                                     (AsnWriteFunc) SeqDescrAsnWrite));
        } else if (sep->choice == 2) {
          bssp = (BioseqSetPtr) sep->data.ptrvalue;
          ValNodeLink (&(bssp->descr),
                       AsnIoMemCopy ((Pointer) subcit,
                                     (AsnReadFunc) SeqDescrAsnRead,
                                     (AsnWriteFunc) SeqDescrAsnWrite));
        }
      }

      fsep = FindNthBioseq (sep, 1);
      if (fsep != NULL && fsep->choice == 1) {
        bsp = (BioseqPtr) fsep->data.ptrvalue;
        if (bsp != NULL) {
          SeqIdWrite (bsp->id, buf, PRINTID_FASTA_LONG, sizeof (buf));
#ifdef OS_UNIX
          if (batch != 1) {
            printf ("%s\n", buf);
            fflush (stdout);
            if (batch != 3 && batch != 4 && batch != 5) {
              if (ofp != NULL) {
                fprintf (ofp, "%s\n", buf);
                fflush (ofp);
              }
            }
          }
#endif
          if (logfp != NULL) {
            fprintf (logfp, "%s\n", buf);
            fflush (logfp);
          }
        }
      }

      hasgi = SeqEntryHasGi (sep, accn);
      if (hasgi) {
        sprintf (buf, "%s.before", accn);
        SaveSeqEntry (sep, buf);
        sprintf (buf, "%s.gbff.before", accn);
        SaveAsn2gnbk (sep, buf, format, SEQUIN_MODE, NORMAL_STYLE, 0, 0, 0);
        if (ofp != NULL) {
          FileClose (ofp);
        }
        AsnIoClose (aip);
        return 0;
      }
      outOfOrder = FALSE;
#ifdef ASN2GNBK_IGNORE_OUT_OF_ORDER
      VisitFeaturesInSep (sep, (Pointer) &outOfOrder, CheckOrder);
#endif
#ifdef ASN2GNBK_REPAIR_OUT_OF_ORDER
      VisitFeaturesInSep (sep, (Pointer) &outOfOrder, CheckOrder);
#endif
      if ((! outOfOrder) && StringHasNoText (accn)) {
        if ((format != GENPEPT_FMT && SeqEntryHasNucs (sep)) ||
           (format == GENPEPT_FMT && SeqEntryHasProts (sep))) {

          hasRefSeq = FALSE;
          VisitBioseqsInSep (sep, (Pointer) &hasRefSeq, LookForRefSeq);
          if (hasRefSeq) {
            if (batch != 1 && format == GENBANK_FMT && ofp == NULL && extra == NULL) {
              ofp = FileOpen (outputFile, "w");
              if (ofp == NULL) {
                ofp = stdout;
              }
            }
          }

          starttime = GetSecs ();
          useFfdiff = (Boolean) (format == GENBANK_FMT && (! hasRefSeq));
          CompareFlatFiles (path1, path2, path3, sep, ofp,
                            format, altformat, mode, style, flags, locks,
                            custom, extra, batch, ffdiff, asn2flat, useFfdiff);
          stoptime = GetSecs ();
          if (stoptime - starttime > worsttime) {
            worsttime = stoptime - starttime;
            StringCpy (longest, buf);
          }
          numrecords++;
        }
      }
      SeqEntryFree (sep);
      omp = ObjMgrGet ();
      ObjMgrReapOne (omp);
      ObjMgrFreeCache (0);
      FreeSeqIdGiCache ();
    } else if (atp == atp_desc && (! atp_se_seen)) {
      descr = SeqDescrAsnRead (aip, atp);
    } else if (atp == atp_sbp) {
      sbp = SubmitBlockAsnRead (aip, atp);
      if (sbp != NULL) {
        csp = sbp->cit;
        if (csp != NULL) {
          MemSet ((Pointer) &ovn, 0, sizeof (ObjValNode));
          MemSet ((Pointer) &pd, 0, sizeof (Pubdesc));
          MemSet ((Pointer) &vn, 0, sizeof (ValNode));
          vn.choice = PUB_Sub;
          vn.data.ptrvalue = (Pointer) csp;
          vn.next = NULL;
          pd.pub = &vn;
          ovn.vn.choice = Seq_descr_pub;
          ovn.vn.data.ptrvalue = (Pointer) &pd;
          ovn.vn.next = NULL;
          ovn.vn.extended = 1;
          subcit = (SeqDescrPtr) &ovn;
        }
      }
    } else {
      AsnReadVal (aip, atp, NULL);
    }
  }

  if (ofp != NULL) {
    FileClose (ofp);
  }

  AsnIoFree (aip, FALSE);

  SeqDescrFree (descr);
  SubmitBlockFree (sbp);

#ifdef OS_UNIX
  if (usedPopen) {
    pclose (fp);
  } else {
    FileClose (fp);
  }
#else
  FileClose (fp);
#endif

  if (logfp != NULL && (! StringHasNoText (longest))) {
    fprintf (logfp, "Longest processing time %ld seconds on %s\n",
             (long) worsttime, longest);
    fprintf (logfp, "Total number of records %ld\n", (long) numrecords);
    fflush (logfp);
  }

  sprintf (cmmd, "rm %s; rm %s; rm %s", path1, path2, path3);
  system (cmmd);

  return 0;
}

#include <lsqfetch.h>
#include <pmfapi.h>
#ifdef INTERNAL_NCBI_ASN2GB
#include <accpubseq.h>
#endif

static void ProcessOneSeqEntry (
  SeqEntryPtr sep,
  CharPtr outputFile,
  FmtType format,
  FmtType altformat,
  ModType mode,
  StlType style,
  FlgType flags,
  LckType locks,
  CstType custom,
  XtraPtr extra,
  Boolean do_tiny_seq,
  Boolean do_fasta_stream
)


{
  AsnIoPtr  aip;
  FILE      *ofp = NULL;

  if (sep == NULL) return;

  if (extra == NULL) {
    FileRemove (outputFile);
#ifdef WIN_MAC
    FileCreate (outputFile, "TEXT", "ttxt");
#endif
    ofp = FileOpen (outputFile, "w");
  }

  if (do_tiny_seq) {
    aip = AsnIoNew (ASNIO_TEXT_OUT | ASNIO_XML, ofp, NULL, NULL, NULL);
    VisitBioseqsInSep (sep, (Pointer) aip, SaveTinySeqs);
    AsnIoFree (aip, FALSE);
  } else if (do_fasta_stream) {
    aip = AsnIoNew (ASNIO_TEXT_OUT | ASNIO_XML, ofp, NULL, NULL, NULL);
    VisitBioseqsInSep (sep, (Pointer) aip, SaveTinyStreams);
    AsnIoFree (aip, FALSE);
  } else {
    SeqEntryToGnbk (sep, NULL, format, mode, style, flags, locks, custom, extra, ofp);
    if (altformat != 0) {
      SeqEntryToGnbk (sep, NULL, altformat, mode, style, flags, locks, custom, extra, ofp);
    }
  }
  if (ofp != NULL) {
    FileClose (ofp);
  }
}

static SeqEntryPtr SeqEntryFromAccnOrGi (
  CharPtr accn
)

{
  Boolean      alldigits;
  BioseqPtr    bsp;
  Char         ch;
  CharPtr      ptr;
  SeqEntryPtr  sep = NULL;
  SeqIdPtr     sip;
  Int4         uid = 0;
  long int     val;
  ValNode      vn;

  if (StringHasNoText (accn)) return NULL;

  TrimSpacesAroundString (accn);

  alldigits = TRUE;
  ptr = accn;
  ch = *ptr;
  while (ch != '\0') {
    if (! IS_DIGIT (ch)) {
      alldigits = FALSE;
    }
    ptr++;
    ch = *ptr;
  }

  if (alldigits) {
    if (sscanf (accn, "%ld", &val) == 1) {
      uid = (Int4) val;
    }
  } else {
    sip = SeqIdFromAccessionDotVersion (accn);
    if (sip != NULL) {
      uid = GetGIForSeqId (sip);
      SeqIdFree (sip);
    }
  }

  if (uid > 0) {
    sep = PubSeqSynchronousQuery (uid, 0, -1);
    if (sep != NULL) {
      MemSet ((Pointer) &vn, 0, sizeof (ValNode));
      vn.choice = SEQID_GI;
      vn.data.intvalue = uid;
      bsp = BioseqFind (&vn);
      if (bsp != NULL) {
        sep = SeqMgrGetSeqEntryForData ((Pointer) bsp);
      }
    }
  }

  return sep;
}

/* Args structure contains command-line arguments */

#define i_argInputFile    0
#define o_argOutputFile   1
#define f_argFormat       2
#define m_argMode         3
#define s_argStyle        4
#define g_argFlags        5
#define h_argLock         6
#define u_argCustom       7
#define a_argType         8
#define t_argBatch        9
#define b_argBinary      10
#define c_argCompressed  11
#define p_argPropagate   12
#define l_argLogFile     13
#define r_argRemote      14
#define A_argAccession   15
#ifdef OS_UNIX
#define q_argFfDiff      16
#define n_argAsn2Flat    17
#define j_argFrom        18
#define k_argTo          19
#define d_argStrand      20
#define y_argItemID      21
#ifdef ENABLE_ARG_X
#define x_argAccnToSave  22
#endif
#endif

Args myargs [] = {
  {"Input File Name", "stdin", NULL, NULL,
    FALSE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"Output File Name", "stdout", NULL, NULL,
    FALSE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Format (b GenBank, e EMBL, p GenPept, t Feature Table, x INSDSet)", "b", NULL, NULL,
    FALSE, 'f', ARG_STRING, 0.0, 0, NULL},
  {"Mode (r Release, e Entrez, s Sequin, d Dump)", "s", NULL, NULL,
    FALSE, 'm', ARG_STRING, 0.0, 0, NULL},
  {"Style (n Normal, s Segment, m Master, c Contig)", "n", NULL, NULL,
    FALSE, 's', ARG_STRING, 0.0, 0, NULL},
  {"Bit Flags (1 HTML, 2 XML, 4 ContigFeats, 8 ContigSrcs, 16 FarTransl)", "0", NULL, NULL,
    FALSE, 'g', ARG_INT, 0.0, 0, NULL},
  {"Lock/Lookup Flags (8 LockProd, 16 LookupComp, 64 LookupProd)", "0", NULL, NULL,
    FALSE, 'h', ARG_INT, 0.0, 0, NULL},
  {"Custom Flags (2 HideMostImpFeats, 4 HideSnpFeats)", "0", NULL, NULL,
    FALSE, 'u', ARG_INT, 0.0, 0, NULL},
  {"ASN.1 Type (a Any, e Seq-entry, b Bioseq, s Bioseq-set, m Seq-submit, t Batch Bioseq-set, u Batch Seq-submit)", "a", NULL, NULL,
    TRUE, 'a', ARG_STRING, 0.0, 0, NULL},
  {"Batch (1 Report, 2 Sequin/Release, 3 asn2gb/asn2flat, 4 asn2flat BSEC/nocleanup, 5 oldasn2gb/newasn2gb)", "0", "0", "5",
    FALSE, 't', ARG_INT, 0.0, 0, NULL},
  {"Input File is Binary", "F", NULL, NULL,
    TRUE, 'b', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Batch File is Compressed", "F", NULL, NULL,
    TRUE, 'c', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Propagate Top Descriptors", "F", NULL, NULL,
    TRUE, 'p', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Log file", NULL, NULL, NULL,
    TRUE, 'l', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Remote Fetching", "F", NULL, NULL,
    TRUE, 'r', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Accession to Fetch", NULL, NULL, NULL,
    TRUE, 'A', ARG_STRING, 0.0, 0, NULL},
#ifdef OS_UNIX
#ifdef PROC_I80X86
  {"Ffdiff Executable", "ffdiff", NULL, NULL,
    TRUE, 'q', ARG_FILE_IN, 0.0, 0, NULL},
  {"Asn2Flat Executable", "asn2flat", NULL, NULL,
    TRUE, 'n', ARG_FILE_IN, 0.0, 0, NULL},
#else
  {"Ffdiff Executable", "/netopt/genbank/subtool/bin/ffdiff", NULL, NULL,
    TRUE, 'q', ARG_FILE_IN, 0.0, 0, NULL},
  {"Asn2Flat Executable", "asn2flat", NULL, NULL,
    TRUE, 'n', ARG_FILE_IN, 0.0, 0, NULL},
#endif
  {"SeqLoc From", "0", NULL, NULL,
    TRUE, 'j', ARG_INT, 0.0, 0, NULL},
  {"SeqLoc To", "0", NULL, NULL,
    TRUE, 'k', ARG_INT, 0.0, 0, NULL},
  {"SeqLoc Minus Strand", "F", NULL, NULL,
    TRUE, 'd', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Feature itemID", "0", NULL, NULL,
    TRUE, 'y', ARG_INT, 0.0, 0, NULL},
#ifdef ENABLE_ARG_X
  {"Accession to extract", NULL, NULL, NULL,
    TRUE, 'x', ARG_STRING, 0.0, 0, NULL},
#endif
#endif
};


#define HTML_XML_ASN_MASK (CREATE_HTML_FLATFILE | CREATE_XML_GBSEQ_FILE | CREATE_ASN_GBSEQ_FILE)

Int2 Main (
  void
)

{
  CharPtr      accn = NULL;
  CharPtr      accntofetch = NULL;
  AsnIoPtr     aip = NULL;
  FmtType      altformat = (FmtType) 0;
  Char         app [64];
  CharPtr      asn2flat = NULL;
  AsnTypePtr   atp = NULL;
  Int2         batch = 0;
  Boolean      binary = FALSE;
  Boolean      compressed = FALSE;
  CstType      custom;
  Boolean      do_gbseq = FALSE;
  Boolean      do_insdseq = FALSE;
  Boolean      do_tiny_seq = FALSE;
  Boolean      do_fasta_stream = FALSE;
  XtraPtr      extra = NULL;
  CharPtr      ffdiff = NULL;
  FlgType      flags;
  FmtType      format = GENBANK_FMT;
  Int4         from = 0;
  GBSeq        gbsq;
  GBSet        gbst;
  Int4         itemID = 0;
  LckType      locks;
  CharPtr      logfile = NULL;
  FILE         *logfp = NULL;
  ModType      mode = SEQUIN_MODE;
  Boolean      propOK = FALSE;
  Boolean      remote = FALSE;
  Int2         rsult = 0;
  time_t       runtime, starttime, stoptime;
  SeqEntryPtr  sep;
  CharPtr      str;
  Uint1        strand = Seq_strand_plus;
  StlType      style = NORMAL_STYLE;
  Int4         to = 0;
  Int2         type = 0;
  Char         xmlbuf [128];
  XtraBlock    xtra;

  /* standard setup */

  ErrSetFatalLevel (SEV_MAX);
  ErrClearOptFlags (EO_SHOW_USERSTR);
  ErrSetLogfile ("stderr", ELOG_APPEND);
  UseLocalAsnloadDataAndErrMsg ();
  ErrPathReset ();

  if (! AllObjLoad ()) {
    Message (MSG_POSTERR, "AllObjLoad failed");
    return 1;
  }
  if (! SubmitAsnLoad ()) {
    Message (MSG_POSTERR, "SubmitAsnLoad failed");
    return 1;
  }
  if (! FeatDefSetLoad ()) {
    Message (MSG_POSTERR, "FeatDefSetLoad failed");
    return 1;
  }
  if (! SeqCodeSetLoad ()) {
    Message (MSG_POSTERR, "SeqCodeSetLoad failed");
    return 1;
  }
  if (! GeneticCodeTableLoad ()) {
    Message (MSG_POSTERR, "GeneticCodeTableLoad failed");
    return 1;
  }

  /* process command line arguments */

  sprintf (app, "asn2gb %s", ASN2GB_APPLICATION);
  if (! GetArgs (app, sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  if (myargs [b_argBinary].intvalue) {
    binary = TRUE;
  } else {
    binary = FALSE;
  }

  if (myargs [c_argCompressed].intvalue) {
    compressed = TRUE;
  } else {
    compressed = FALSE;
  }

  if (myargs [p_argPropagate].intvalue) {
    propOK = TRUE;
  } else {
    propOK = FALSE;
  }

  str = myargs [f_argFormat].strvalue;
  if (StringICmp (str, "bp") == 0 || StringICmp (str, "pb") == 0) {
    format = GENBANK_FMT;
    altformat = GENPEPT_FMT;

  } else if (StringICmp (str, "b") == 0) {
    format = GENBANK_FMT;
  } else if (StringICmp (str, "e") == 0) {
    format = EMBL_FMT;
  } else if (StringICmp (str, "p") == 0) {
    format = GENPEPT_FMT;
  } else if (StringICmp (str, "t") == 0) {
    format = FTABLE_FMT;
  
  } else if (StringICmp (str, "q") == 0) {
    do_gbseq = TRUE;
    format = GENBANK_FMT;
  } else if (StringICmp (str, "r") == 0) {
    do_gbseq = TRUE;
    format = GENPEPT_FMT;

  } else if (StringICmp (str, "x") == 0) {
    do_gbseq = TRUE;
    do_insdseq = TRUE;
    format = GENBANK_FMT;
  } else if (StringCmp (str, "y") == 0) {
    do_tiny_seq = TRUE;
    format = GENBANK_FMT;
  } else if (StringCmp (str, "Y") == 0) {
    do_fasta_stream = TRUE;
    format = GENBANK_FMT;
  } else if (StringICmp (str, "z") == 0) {
    do_gbseq = TRUE;
    do_insdseq = TRUE;
    format = GENPEPT_FMT;
  } else {
    format = GENBANK_FMT;
  }

  str = myargs [m_argMode].strvalue;
  if (StringICmp (str, "r") == 0) {
    mode = RELEASE_MODE;
  } else if (StringICmp (str, "e") == 0) {
    mode = ENTREZ_MODE;
  } else if (StringICmp (str, "s") == 0) {
    mode = SEQUIN_MODE;
  } else if (StringICmp (str, "d") == 0) {
    mode = DUMP_MODE;
  } else {
    mode = SEQUIN_MODE;
  }

  str = myargs [s_argStyle].strvalue;
  if (StringICmp (str, "n") == 0) {
    style = NORMAL_STYLE;
  } else if (StringICmp (str, "s") == 0) {
    style = SEGMENT_STYLE;
  } else if (StringICmp (str, "m") == 0) {
    style = MASTER_STYLE;
  } else if (StringICmp (str, "c") == 0) {
    style = CONTIG_STYLE;
  } else {
    style = NORMAL_STYLE;
  }

  flags = (FlgType) myargs [g_argFlags].intvalue;

  locks = (LckType) myargs [h_argLock].intvalue;

  custom = (CstType) myargs [u_argCustom].intvalue;

  str = myargs [a_argType].strvalue;
  if (StringICmp (str, "a") == 0) {
    type = 1;
  } else if (StringICmp (str, "e") == 0) {
    type = 2;
  } else if (StringICmp (str, "b") == 0) {
    type = 3;
  } else if (StringICmp (str, "s") == 0) {
    type = 4;
  } else if (StringICmp (str, "m") == 0) {
    type = 5;
  } else if (StringICmp (str, "t") == 0) {
    batch = 1;
    type = 4;
  } else if (StringICmp (str, "u") == 0) {
    batch = 1;
    type = 5;
  } else {
    type = 1;
  }

  if (myargs [t_argBatch].intvalue > 0) {
    batch = (Int2) myargs [t_argBatch].intvalue;
  }

  if ((binary || compressed) && batch == 0) {
    if (type == 1) {
      Message (MSG_FATAL, "-b or -c cannot be used without -t or -a");
      return 1;
    }
  }

  remote = (Boolean) myargs [r_argRemote].intvalue;

  accntofetch = (CharPtr) myargs [A_argAccession].strvalue;
  if (StringDoesHaveText (accntofetch)) {
    remote = TRUE;
  }

  if (remote) {
#ifdef INTERNAL_NCBI_ASN2GB
    if (! PUBSEQBioseqFetchEnable ("asn2gb", FALSE)) {
      Message (MSG_POSTERR, "PUBSEQBioseqFetchEnable failed");
      return 1;
    }
#else
    PubSeqFetchEnable ();
#endif
    PubMedFetchEnable ();
    LocalSeqFetchInit (FALSE);
  }

  logfile = (CharPtr) myargs [l_argLogFile].strvalue;
  if (! StringHasNoText (logfile)) {
    logfp = FileOpen (logfile, "w");
  }

#ifdef OS_UNIX
  ffdiff = myargs [q_argFfDiff].strvalue;
  asn2flat = myargs [n_argAsn2Flat].strvalue;

  from = myargs [j_argFrom].intvalue;
  to = myargs [k_argTo].intvalue;
  if (myargs [d_argStrand].intvalue) {
    strand = Seq_strand_minus;
  } else {
    strand = Seq_strand_plus;
  }
  itemID = myargs [y_argItemID].intvalue;

#ifdef ENABLE_ARG_X
  if (! StringHasNoText (myargs [x_argAccnToSave].strvalue)) {
    accn = myargs [x_argAccnToSave].strvalue;
  }
#endif
#endif

  if (GetAppParam ("NCBI", "SETTINGS", "XMLPREFIX", NULL, xmlbuf, sizeof (xmlbuf))) {
    AsnSetXMLmodulePrefix (StringSave (xmlbuf));
  }

  if (do_gbseq) {
    if (! objgbseqAsnLoad ()) {
      Message (MSG_POSTERR, "objgbseqAsnLoad failed");
      return 1;
    }
    if (! objinsdseqAsnLoad ()) {
      Message (MSG_POSTERR, "objinsdseqAsnLoad failed");
      return 1;
    }
    MemSet ((Pointer) &xtra, 0, sizeof (XtraBlock));
    MemSet ((Pointer) &gbsq, 0, sizeof (GBSeq));
    xtra.gbseq = &gbsq;
    if ((flags & HTML_XML_ASN_MASK) == CREATE_ASN_GBSEQ_FILE) {
      aip = AsnIoOpen (myargs [o_argOutputFile].strvalue, "w");
    } else {
      aip = AsnIoOpen (myargs [o_argOutputFile].strvalue, "wx");
    }
    if (aip == NULL) {
      Message (MSG_POSTERR, "AsnIoOpen failed");
      return 1;
    }
    xtra.aip = aip;
    if ((Boolean) ((flags & PRODUCE_OLD_GBSEQ) != 0)) {
      do_insdseq = FALSE;
    }
    if (do_insdseq) {
      atp = AsnLinkType (NULL, AsnFind ("INSDSet"));
      xtra.atp = AsnLinkType (NULL, AsnFind ("INSDSet.E"));
    } else {
      atp = AsnLinkType (NULL, AsnFind ("GBSet"));
      xtra.atp = AsnLinkType (NULL, AsnFind ("GBSet.E"));
      flags |= PRODUCE_OLD_GBSEQ;
    }
    if (atp == NULL || xtra.atp == NULL) {
      Message (MSG_POSTERR, "AsnLinkType or AsnFind failed");
      return 1;
    }
    extra = &xtra;
    MemSet ((Pointer) &gbst, 0, sizeof (GBSet));
    AsnOpenStruct (aip, atp, (Pointer) &gbst);
  }

  starttime = GetSecs ();

  if (StringDoesHaveText (accntofetch)) {

    if (remote) {
      sep = SeqEntryFromAccnOrGi (accntofetch);
      if (sep != NULL) {
        ProcessOneSeqEntry (sep, myargs [o_argOutputFile].strvalue,
                            format, altformat, mode, style, flags, locks,
                            custom, extra, do_tiny_seq, do_fasta_stream);
        SeqEntryFree (sep);
      }
    }

  } else if (batch != 0 || accn != NULL) {

    rsult = HandleMultipleRecords (myargs [i_argInputFile].strvalue,
                                   myargs [o_argOutputFile].strvalue,
                                   format, altformat, mode, style, flags, locks,
                                   custom, extra, type, batch, binary, compressed,
                                   propOK, ffdiff, asn2flat, accn, logfp);
  } else {

    rsult = HandleSingleRecord (myargs [i_argInputFile].strvalue,
                                myargs [o_argOutputFile].strvalue,
                                format, altformat, mode, style, flags, locks,
                                custom, extra, type, binary, compressed,
                                from, to, strand, itemID, do_tiny_seq, do_fasta_stream);
  }

  if (aip != NULL) {
    AsnCloseStruct (aip, atp, NULL);
    AsnPrintNewLine (aip);
    AsnIoClose (aip);
  }

  stoptime = GetSecs ();
  runtime = stoptime - starttime;
  if (logfp != NULL) {
    fprintf (logfp, "Finished in %ld seconds\n", (long) runtime);
    FileClose (logfp);
  }

  if (remote) {
    LocalSeqFetchDisable ();
    PubMedFetchDisable ();
#ifdef INTERNAL_NCBI_ASN2GB
    PUBSEQBioseqFetchDisable ();
#else
    PubSeqFetchDisable ();
#endif
  }

  return rsult;
}

