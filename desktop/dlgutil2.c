/*   dlgutil2.c
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
* File Name:  dlgutil2.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   1/22/95
*
* $Revision: 6.108 $
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
#include <utilpub.h>
#include <objfeat.h>
#include <objseq.h>
#include <toasn3.h>
#ifdef WIN_MOTIF
#include <netscape.h>
#endif

typedef struct datepage {
  DIALOG_MESSAGE_BLOCK
  TexT          year;
  PopuP         month;
  TexT          day;
} DatePage, PNTR DatePagePtr;

extern CharPtr SaveStringFromTextAndStripNewlines (TexT t)

{
  Char     ch;
  size_t   len;
  CharPtr  ptr;
  CharPtr  str;

  len = TextLength (t);
  if (len > 0) {
    str = MemNew (len + 1);
    if (str != NULL) {
      GetTitle (t, str, len + 1);
      ptr = str;
      ch = *ptr;
      while (ch != '\0') {
        if (ch < ' ') {
          *ptr = ' ';
        }
        ptr++;
        ch = *ptr;
      }
      TrimSpacesAroundString (str);
      if (StringHasNoText (str)) {
        str = MemFree (str);
      }
      return str;
    } else {
      return NULL;
    }
  } else {
    return NULL;
  }
}

static void DatePtrToDatePage (DialoG d, Pointer data)

{
  DatePtr      dp;
  DatePagePtr  dpp;

  dpp = (DatePagePtr) GetObjectExtra (d);
  dp = (DatePtr) data;
  if (dpp != NULL) {
    DatePtrToVibrant (dp, dpp->month, dpp->day, dpp->year);
  }
}

static Pointer DatePageToDatePtr (DialoG d)

{
  DatePtr      dp;
  DatePagePtr  dpp;

  dp = NULL;
  dpp = (DatePagePtr) GetObjectExtra (d);
  if (dpp != NULL) {
    dp = VibrantToDatePtr (dpp->month, dpp->day, dpp->year);
  }
  return (Pointer) dp;
}

extern DialoG CreateDateDialog (GrouP prnt, CharPtr title)

{
  DatePagePtr  dpp;
  GrouP        f;
  GrouP        m;
  GrouP        p;
  GrouP        s;

  p = HiddenGroup (prnt, 1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  dpp = (DatePagePtr) MemNew (sizeof (DatePage));
  if (dpp) {

    SetObjectExtra (p, dpp, StdCleanupExtraProc);
    dpp->dialog = (DialoG) p;
    dpp->todialog = DatePtrToDatePage;
    dpp->fromdialog = DatePageToDatePtr;
    dpp->testdialog = NULL;

    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    m = HiddenGroup (s, -1, 0, NULL);
    /*
    SetGroupSpacing (m, 10, 10);
    */

    f = HiddenGroup (m, -6, 0, NULL);
    StaticPrompt (f, "Month", 0, popupMenuHeight, programFont, 'l');
    dpp->month = PopupList (f, TRUE, NULL);
    InitEnumPopup (dpp->month, months_alist, NULL);
    SetValue (dpp->month, 1);
    StaticPrompt (f, "Day", 0, dialogTextHeight, programFont, 'l');
    dpp->day = DialogText (f, "", 4, NULL);
    StaticPrompt (f, "Year", 0, dialogTextHeight, programFont, 'l');
    dpp->year = DialogText (f, "", 6, NULL);
  }

  return (DialoG) p;
}

typedef struct featcit {
  DIALOG_MESSAGE_BLOCK
  DoC         citdoc;
  ValNodePtr  pubset;
} FeatCitPage, PNTR FeatCitPagePtr;

static ParData cofParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static ColData cofColFmt = {0, 0, 0, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, TRUE};

static ParData cofeParFmt = {TRUE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static ColData cofeColFmt = {0, 0, 0, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, TRUE};

static Uint1 diamondSym [] = {
  0x00, 0x10, 0x38, 0x7C, 0x38, 0x10, 0x00, 0x00
};

static void DrawCitOnFeat (DoC d, RectPtr r, Int2 item, Int2 firstLine)

{
  Int2  lineHeight;
  RecT  rct;

  if (d != NULL && r != NULL && item > 0 && firstLine == 0) {
    if (item == 1) return; /* citonfeattxt or citonfeathdr explanatory text */
    GetItemParams (d, item, NULL, NULL, NULL, &lineHeight, NULL);
    rct = *r;
    rct.left += 1;
    rct.right = rct.left + 7;
    rct.top += (lineHeight - 7) / 2;
    rct.bottom = rct.top + 7;
    CopyBits (&rct, diamondSym);
    /*
    x = r->left + 1;
    y = r->top + lineHeight / 2;
    MoveTo (x, y);
    LineTo (x + 5, y);
    */
  }
}

static int LIBCALLBACK CompareStrings (VoidPtr ptr1, VoidPtr ptr2)

{
  CharPtr  str1;
  CharPtr  str2;

  if (ptr1 != NULL && ptr2 != NULL) {
    str1 = *((CharPtr PNTR) ptr1);
    str2 = *((CharPtr PNTR) ptr2);
    if (str1 != NULL && str2 != NULL) {
      return StringICmp (str1, str2);
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}

static CharPtr citonfeattxt =
"Press 'Edit Citations' to attach publications to this feature. \
Publications must first be added to the record. For biological \
justification, and not to credit the sequencer, create publication \
with the 'Cites a feature on the sequence' scope.\n";

static CharPtr citonfeathdr =
"Press 'Edit Citations' to change publications.\n\n";

static void PubsetPtrToFeatCitPage (DialoG d, Pointer data)

{
  Int2            count;
  FeatCitPagePtr  fpp;
  Int2            i;
  Char            label [128];
  ObjMgrPtr       omp;
  ObjMgrTypePtr   omtp;
  ValNodePtr      ppr;
  ValNodePtr      psp;
  RecT            r;
  CharPtr         PNTR strs;

  fpp = (FeatCitPagePtr) GetObjectExtra (d);
  psp = (ValNodePtr) data;
  if (fpp != NULL) {
    Reset (fpp->citdoc);
    fpp->pubset = PubSetFree (fpp->pubset);
    fpp->pubset = AsnIoMemCopy (data,
                                (AsnReadFunc) PubSetAsnRead,
                                (AsnWriteFunc) PubSetAsnWrite);
    SetDocAutoAdjust (fpp->citdoc, FALSE);
    ObjectRect (fpp->citdoc, &r);
    InsetRect (&r, 4, 4);
    cofColFmt.pixWidth = r.right - r.left;
    cofColFmt.pixInset = 10;
    omp = ObjMgrGet ();
    if (omp == NULL) return;
    omtp = ObjMgrTypeFind (omp, OBJ_SEQFEAT_CIT, NULL, NULL);
    if (omtp == NULL || omtp->labelfunc == NULL) return;
    if (psp != NULL && psp->data.ptrvalue != NULL) {
      count = 0;
      for (ppr = psp->data.ptrvalue; ppr != NULL; ppr = ppr->next) {
        count++;
      }
      if (count > 0) {
        strs = MemNew (sizeof (CharPtr) * (size_t) (count + 1));
        if (strs != NULL) {
          i = 0;
          for (ppr = psp->data.ptrvalue; ppr != NULL; ppr = ppr->next) {
            (*(omtp->labelfunc)) (ppr, label, 127, OM_LABEL_CONTENT);
            strs [i] = StringSave (label);
            i++;
          }
          HeapSort (strs, count, sizeof (CharPtr), CompareStrings);
          AppendText (fpp->citdoc, citonfeathdr, &cofParFmt, NULL, programFont);
          for (i = 0; i < count; i++) {
            AppendText (fpp->citdoc, strs [i],
                        &cofParFmt, &cofColFmt, programFont);
          }
          for (i = 0; i < count; i++) {
            strs [i] = MemFree (strs [i]);
          }
          MemFree (strs);
        } else {
          AppendText (fpp->citdoc, citonfeattxt, &cofParFmt, NULL, programFont);
        }
      }
    } else {
      AppendText (fpp->citdoc, citonfeattxt, &cofParFmt, NULL, programFont);
    }
    SetDocShade (fpp->citdoc, DrawCitOnFeat, NULL, NULL, NULL);
    UpdateDocument (fpp->citdoc, 0, 0);
  }
}

static Pointer FeatCitPageToPubsetPtr (DialoG d)

{
  FeatCitPagePtr  fpp;
  ValNodePtr      psp;

  psp = NULL;
  fpp = (FeatCitPagePtr) GetObjectExtra (d);
  if (fpp != NULL) {
    psp = AsnIoMemCopy (fpp->pubset,
                        (AsnReadFunc) PubSetAsnRead,
                        (AsnWriteFunc) PubSetAsnWrite);
  }
  return (Pointer) psp;
}

static void CleanupCitOnFeatProc (GraphiC g, VoidPtr data)

{
  FeatCitPagePtr  fpp;

  fpp = (FeatCitPagePtr) data;
  if (fpp != NULL) {
    PubSetFree (fpp->pubset);
  }
  MemFree (data);
}

static DialoG CreateCitOnFeatDialog (GrouP h, CharPtr title)

{
  FeatCitPagePtr  fpp;
  Int2            lineHeight;
  GrouP           m;
  GrouP           p;
  GrouP           s;

  p = HiddenGroup (h, 1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  fpp = (FeatCitPagePtr) MemNew (sizeof (FeatCitPage));
  if (fpp != NULL) {

    SetObjectExtra (p, fpp, CleanupCitOnFeatProc);
    fpp->dialog = (DialoG) p;
    fpp->todialog = PubsetPtrToFeatCitPage;
    fpp->fromdialog = FeatCitPageToPubsetPtr;
    fpp->testdialog = NULL;

    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    m = HiddenGroup (s, -1, 0, NULL);
    /*
    SetGroupSpacing (m, 10, 10);
    */

    SelectFont (programFont);
    lineHeight = LineHeight ();
    SelectFont (systemFont);
    cofParFmt.minHeight = lineHeight + 2;
    StaticPrompt (m, "Citations on Feature", 25 * stdCharWidth, 0, programFont, 'c');
    fpp->citdoc = DocumentPanel (m, 25 * stdCharWidth, 5 * cofParFmt.minHeight);
    SetObjectExtra (fpp->citdoc, fpp, NULL);
    SetDocAutoAdjust (fpp->citdoc, FALSE);
    fpp->pubset = NULL;
  }

  return (DialoG) p;
}

typedef struct citlist {
  Uint2    entityID;
  Uint2    itemID;
  Uint2    itemtype;
  CharPtr  label;
} CitListData, PNTR CitListDataPtr;

typedef struct featcitedit {
  DIALOG_MESSAGE_BLOCK
  DoC             allcitdoc;
  Int2            clickedItem;
  Int2            clickedRow;
  Int2            numitems;
  Int2            lineheight;
  Int2            index;
  BoolPtr         chosen;
  ValNodePtr      citlist;
  ValNodePtr      psp;
  ObjMgrPtr       omp;
  Int2            entityID;
} FeatCitEdit, PNTR FeatCitEditPtr;

static void CleanupFeatCitForm (GraphiC g, VoidPtr data)

{
  FeatCitEditPtr  fcep;

  fcep = (FeatCitEditPtr) data;
  if (fcep != NULL) {
    MemFree (fcep->chosen);
  }
  MemFree (data);
}

static void DrawFeatCit (DoC d, RectPtr r, Int2 item, Int2 firstLine)

{
  FeatCitEditPtr  fcep;
  RecT            rct;

  fcep = (FeatCitEditPtr) GetObjectExtra (d);
  if (fcep != NULL && r != NULL && item > 0 && firstLine == 0) {
    rct = *r;
    rct.left += 1;
    rct.right = rct.left + fcep->lineheight;
    rct.bottom = rct.top + (rct.right - rct.left);
    FrameRect (&rct);
    if (item > 0 && item <= fcep->numitems) {
      if (fcep->chosen != NULL && fcep->chosen [item - 1]) {
        MoveTo (rct.left, rct.top);
        LineTo (rct.right - 1, rct.bottom - 1);
        MoveTo (rct.left, rct.bottom - 1);
        LineTo (rct.right - 1, rct.top);
      }
    }
  }
}

static void ClickFeatCit (DoC d, PoinT pt)

{
}

static void ReleaseFeatCit (DoC d, PoinT pt)

{
  Int2            col;
  FeatCitEditPtr  fcep;
  Int2            item;
  RecT            rct;
  Int2            row;

  fcep = (FeatCitEditPtr) GetObjectExtra (d);
  if (fcep != NULL && fcep->chosen != NULL) {
    MapDocPoint (d, pt, &item, &row, &col, &rct);
    rct.left += 1;
    rct.right = rct.left + fcep->lineheight;
    rct.bottom = rct.top + (rct.right - rct.left);
    if (row == 1 && col == 1 && item > 0 && item <= fcep->numitems && PtInRect (pt, &rct)) {
      if (fcep->chosen [item - 1]) {
        fcep->chosen [item - 1] = FALSE;
      } else {
        fcep->chosen [item - 1] = TRUE;
      }
      InsetRect (&rct, -1, -1);
      InvalRect (&rct);
      Update ();
    }
  }
}

static Boolean GatherAllCits (GatherContextPtr gcp)

{
  CitListDataPtr  cldp;
  FeatCitEditPtr  fcep;
  Char            label [128];
  ObjMgrTypePtr   omtp;
  PubdescPtr      pdp;
  ValNodePtr      sdp;
  SeqFeatPtr      sfp;
  ValNodePtr      vnp;

  if (gcp == NULL) return TRUE;
  fcep = (FeatCitEditPtr) gcp->userdata;
  if (fcep == NULL) return TRUE;
  label [0] = '\0';
  pdp = NULL;
  if (gcp->thistype == OBJ_SEQDESC) {
    sdp = (ValNodePtr) gcp->thisitem;
    if (sdp == NULL || sdp->choice != Seq_descr_pub) return TRUE;
    pdp = sdp->data.ptrvalue;
  } else if (gcp->thistype == OBJ_SEQFEAT) {
    sfp = (SeqFeatPtr) gcp->thisitem;
    if (sfp == NULL || sfp->data.choice != SEQFEAT_PUB) return TRUE;
    pdp = sfp->data.value.ptrvalue;
  } else return TRUE;
  if (pdp == NULL) return TRUE;
  omtp = ObjMgrTypeFind (fcep->omp, gcp->thistype, NULL, NULL);
  if (omtp == NULL || omtp->labelfunc == NULL) return TRUE;
  (*(omtp->labelfunc)) (gcp->thisitem, label, 127, OM_LABEL_CONTENT);
  cldp = (CitListDataPtr) MemNew (sizeof (CitListData));
  if (cldp != NULL) {
    vnp = ValNodeNew (fcep->citlist);
    if (fcep->citlist == NULL) {
      fcep->citlist = vnp;
    }
    if (vnp != NULL) {
      vnp->data.ptrvalue = cldp;
      cldp->entityID = gcp->entityID;
      cldp->itemID = gcp->itemID;
      cldp->itemtype = gcp->thistype;
      cldp->label = StringSave (label);
      (fcep->numitems)++;
    } else {
      MemFree (cldp);
      return TRUE;
    }
  }
  return TRUE;
}

static Boolean PrintCitOnFeatItem (GatherContextPtr gcp)

{
  CharPtr  PNTR  rsultp;
  ValNodePtr     sdp;
  SeqFeatPtr     sfp;
  Boolean        success;

  rsultp = (CharPtr PNTR) gcp->userdata;
  if (rsultp != NULL) {
    success = FALSE;
    switch (gcp->thistype) {
      case OBJ_SEQDESC :
        sdp = (ValNodePtr) gcp->thisitem;
        if (sdp->data.ptrvalue != NULL) {
          success = StdFormatPrint ((Pointer) sdp, (AsnWriteFunc) SeqDescAsnWrite,
                                    "StdSeqDesc", spop);
        } else {
          *rsultp = StringSave ("Empty Descriptor\n");
          success = TRUE;
        }
        break;
      case OBJ_SEQFEAT :
        sfp = (SeqFeatPtr) gcp->thisitem;
        if (sfp != NULL && (sfp->data.choice == 10 || sfp->data.value.ptrvalue != NULL)) {
          success = StdFormatPrint ((Pointer) sfp, (AsnWriteFunc) SeqFeatAsnWrite,
                                    "StdSeqFeat", spop);
        } else {
          *rsultp = StringSave ("Empty Feature\n");
          success = TRUE;
        }
        break;
      default :
        break;
    }
    if (success) {
      if (spop->ptr != NULL && *((CharPtr) (spop->ptr)) != '\0') {
        *rsultp = spop->ptr;
        spop->ptr = NULL;
      } else {
        *rsultp = StringSave ("Empty Data\n");
      }
    } else {
      *rsultp = StringSave ("Data Failure\n");
    }
  }
  return TRUE;
}

static CharPtr CitOnFeatPrintProc (DoC doc, Int2 item, Pointer data)

{
  unsigned int  entityID;
  unsigned int  itemID;
  unsigned int  itemtype;
  CharPtr       rsult;
  CharPtr       str;

  rsult = NULL;
  if (data != NULL) {
    str = (CharPtr) data;
    if (sscanf (str, "%u %u %u", &entityID, &itemID, &itemtype) == 3) {
      GatherItem ((Uint2) entityID, (Uint2) itemID, (Uint2) itemtype,
                  (Pointer) &rsult, PrintCitOnFeatItem);
    }
  } else {
    rsult = StringSave ("Null Data\n");
  }
  return rsult;
}

static int LIBCALLBACK CompareCitList (VoidPtr ptr1, VoidPtr ptr2)

{
  CitListDataPtr  cldp1;
  CitListDataPtr  cldp2;
  CharPtr         str1;
  CharPtr         str2;

  if (ptr1 != NULL && ptr2 != NULL) {
    cldp1 = (CitListDataPtr) ptr1;
    cldp2 = (CitListDataPtr) ptr2;
    if (cldp1 != NULL && cldp2 != NULL) {
      str1 = cldp1->label;
      str2 = cldp2->label;
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

static Boolean MakeMinimalCitOnFeatItem (GatherContextPtr gcp)

{
  PubdescPtr  pdp;
  ValNodePtr  ppr;
  ValNodePtr  sdp;
  SeqFeatPtr  sfp;
  ValNodePtr  vnp;
  ValNodePtr  PNTR  vnpp;

  vnpp = (ValNodePtr PNTR) gcp->userdata;
  if (vnpp != NULL) {
    pdp = NULL;
    switch (gcp->thistype) {
      case OBJ_SEQDESC :
        sdp = (ValNodePtr) gcp->thisitem;
        if (sdp->data.ptrvalue != NULL) {
          pdp = (PubdescPtr) sdp->data.ptrvalue;
        }
        break;
      case OBJ_SEQFEAT :
        sfp = (SeqFeatPtr) gcp->thisitem;
        if (sfp != NULL && (sfp->data.choice == 10 || sfp->data.value.ptrvalue != NULL)) {
          pdp = (PubdescPtr) sfp->data.value.ptrvalue;
        }
        break;
      default :
        break;
    }
    if (pdp != NULL) {
      ppr = NULL;
      vnp = ValNodeNew (NULL);
      if (vnp != NULL) {
        vnp->choice = PUB_Equiv;
        vnp->data.ptrvalue = pdp->pub;
        ppr = MinimizePub (vnp);
        ValNodeFree (vnp);
      }
      vnp = ValNodeNew (*vnpp);
      if (*vnpp == NULL) {
        *vnpp = vnp;
      }
      if (vnp != NULL) {
        vnp->choice = PUB_Equiv;
        vnp->data.ptrvalue = ppr;
      }
    }
  }
  return TRUE;
}

static Boolean MatchMinimalCits (GatherContextPtr gcp)

{
  FeatCitEditPtr  fcep;
  PubdescPtr      pdp;
  ValNodePtr      ppr;
  ValNodePtr      sdp;
  SeqFeatPtr      sfp;
  ValNodePtr      vnp;

  fcep = (FeatCitEditPtr) gcp->userdata;
  if (fcep != NULL) {
    pdp = NULL;
    switch (gcp->thistype) {
      case OBJ_SEQDESC :
        sdp = (ValNodePtr) gcp->thisitem;
        if (sdp->data.ptrvalue != NULL) {
          pdp = (PubdescPtr) sdp->data.ptrvalue;
        }
        break;
      case OBJ_SEQFEAT :
        sfp = (SeqFeatPtr) gcp->thisitem;
        if (sfp != NULL && (sfp->data.choice == 10 || sfp->data.value.ptrvalue != NULL)) {
          pdp = (PubdescPtr) sfp->data.value.ptrvalue;
        }
        break;
      default :
        break;
    }
    if (pdp != NULL && fcep->psp != NULL) {
      vnp = ValNodeNew (NULL);
      if (vnp != NULL) {
        vnp->choice = PUB_Equiv;
        vnp->data.ptrvalue = pdp->pub;
        for (ppr = fcep->psp->data.ptrvalue; ppr != NULL; ppr = ppr->next) {
          if (PubLabelMatch (vnp, ppr) == 0) {
            fcep->chosen [fcep->index] = TRUE;
          }
        }
        ValNodeFree (vnp);
      }
    }
  }
  return TRUE;
}

static void CitListToDialog (DialoG d, Pointer userdata)
{
  FeatCitEditPtr  fcep;
  CitListDataPtr  cldp;
  CitListDataPtr  cldpp;
  Int2            count;
  GatherScope     gs;
  Int2            i;
  Int2            j;
  Char            last [128];
  CharPtr         ptr;
  RecT            r;
  Char            str [34];
  ValNodePtr      vnp;

  fcep = (FeatCitEditPtr) GetObjectExtra (d);
  if (fcep == NULL)
  {
    return;
  }
  Reset (fcep->allcitdoc);
  
  fcep->citlist = ValNodeFree (fcep->citlist);
  fcep->numitems = 0;
  fcep->chosen = MemFree (fcep->chosen);
  
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  MemSet ((Pointer) (gs.ignore), (int) (TRUE), (size_t) (OBJ_MAX * sizeof (Boolean)));
  gs.ignore [OBJ_SEQDESC] = FALSE;
  gs.ignore [OBJ_SEQANNOT] = FALSE;
  gs.ignore [OBJ_SEQFEAT] = FALSE;
  gs.seglevels = 1;
  GatherEntity (fcep->entityID, (Pointer) fcep, GatherAllCits, &gs);
  count = fcep->numitems;
  cldpp = (CitListDataPtr) MemNew (sizeof (CitListData) * (size_t) (count + 1));
  if (cldpp != NULL) {
    for (i = 0, vnp = fcep->citlist; i < count && vnp != NULL; i++, vnp = vnp->next) {
      cldp = vnp->data.ptrvalue;
      if (cldp != NULL) {
        cldpp [i] = *cldp;
        cldp->label = NULL;
      }
    }
    fcep->citlist = ValNodeFreeData (fcep->citlist);
    HeapSort (cldpp, count, sizeof (CitListData), CompareCitList);
    last [0] = '\0';
    ObjectRect (fcep->allcitdoc, &r);
    InsetRect (&r, 4, 4);
    cofeColFmt.pixWidth = r.right - r.left;
    cofeColFmt.pixInset = 20;
    fcep->numitems = 0;
    fcep->chosen = MemNew (sizeof (Boolean) * (size_t) (count + 1));

    fcep->psp = (ValNodePtr) userdata;
    for (i = 0, j = 0; i < count; i++) {
      cldp = &(cldpp [i]);
      if (cldp != NULL) {
        if (last == NULL || StringCmp (cldp->label, last) != 0) {
          sprintf (str, "%d %d %d", (int) cldp->entityID,
                   (int) cldp->itemID, (int) cldp->itemtype);
          ptr = StringSave (str);
          (fcep->numitems)++;
          AppendItem (fcep->allcitdoc, CitOnFeatPrintProc, (Pointer) ptr,
                      TRUE, 5, &cofeParFmt, &cofeColFmt, programFont);
          if (fcep->chosen != NULL) {
            fcep->index = j;
            GatherItem ((Uint2) cldp->entityID, (Uint2) cldp->itemID,
                        (Uint2) cldp->itemtype, (Pointer) fcep, MatchMinimalCits);
          }
          j++;
        }
        StringNCpy_0 (last, cldp->label, sizeof (last));
        cldpp [i].label = MemFree (cldpp [i].label);
      }
    }
    fcep->psp = NULL;
  }
}

static Pointer DialogToMinimizedCitList (DialoG d)
{
  FeatCitEditPtr  fcep;
  Pointer         data;
  unsigned int    entityID;
  Int2            i;
  unsigned int    itemID;
  unsigned int    itemtype;
  Int2            numItems;
  ValNodePtr      ppr = NULL, psp = NULL;
  CharPtr         str;
  
  fcep = (FeatCitEditPtr) GetObjectExtra (d);
  if (fcep == NULL)
  {
    return NULL;
  }
  
  if (fcep->chosen != NULL) {
    GetDocParams (fcep->allcitdoc, &numItems, NULL);
    for (i = 1; i <= numItems; i++) {
      if (fcep->chosen [i - 1]) {
        GetItemParams (fcep->allcitdoc, i, NULL, NULL, NULL, NULL, &data);
        if (data != NULL) {
          str = (CharPtr) data;
          if (sscanf (str, "%u %u %u", &entityID, &itemID, &itemtype) == 3) {
            GatherItem ((Uint2) entityID, (Uint2) itemID, (Uint2) itemtype,
                        (Pointer) &ppr, MakeMinimalCitOnFeatItem);
          }
        }
      }
    }
  }
 
  if (ppr != NULL)
  {
    psp = ValNodeNew (NULL);
    if (psp != NULL) {
      psp->choice = 1;
      psp->data.ptrvalue = ppr;
    }
  }
  return psp;
}

extern DialoG FeatCitEditDialog (GrouP parent, Uint2 entityID)
{
  FeatCitEditPtr  fcep;
  GrouP           p;
  
  fcep = (FeatCitEditPtr) MemNew (sizeof (FeatCitEdit));
  if (fcep == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (parent, -1, 0, NULL);
  SetObjectExtra (p, fcep, CleanupFeatCitForm);

  fcep->dialog = (DialoG) p;
  fcep->todialog = CitListToDialog;
  fcep->fromdialog = DialogToMinimizedCitList;
  fcep->dialogmessage = NULL;
  fcep->testdialog = NULL;

  fcep->entityID = entityID;
  fcep->omp = ObjMgrGet ();
  fcep->numitems = 0;

  fcep->allcitdoc = DocumentPanel (p, 25 * stdCharWidth, 15 * stdLineHeight);
  SetObjectExtra (fcep->allcitdoc, fcep, NULL);
  SetDocProcs (fcep->allcitdoc, ClickFeatCit, NULL, ReleaseFeatCit, NULL);
  SetDocShade (fcep->allcitdoc, DrawFeatCit, NULL, NULL, NULL);
  
  SelectFont (programFont);
  fcep->lineheight = LineHeight ();
  SelectFont (systemFont);

  return (DialoG) p;
}

typedef struct featcitationform
{
  FeatureFormPtr ffp;
  DialoG         citation_list;
   
} FeatCitationFormData, PNTR FeatCitationFormPtr;

static void AcceptFeatCit (ButtoN b)

{
  FeatCitationFormPtr  fcfp;
  FeatureFormPtr       ffp;
  ValNodePtr           psp = NULL;

  fcfp = (FeatCitationFormPtr) GetObjectExtra (b);
  if (fcfp == NULL)
  {
    return;
  }
  
  psp = DialogToPointer (fcfp->citation_list);
  ffp = (FeatureFormPtr) fcfp->ffp;
  if (ffp != NULL) {
    PointerToDialog (ffp->featcits, (Pointer) psp);
    PubSetFree (psp);
  }
  Remove (ParentWindow (b));
}

static void EditFeatCitsProc (ButtoN b)

{
  ButtoN          btn;
  GrouP           c;
  FeatCitationFormPtr  fcfp;
  FeatureFormPtr  ffp;
  WindoW          w;
  ValNodePtr      psp;

  ffp = (FeatureFormPtr) GetObjectExtra (b);
  if (ffp == NULL)
  {
    return;
  }
  fcfp = (FeatCitationFormPtr) MemNew (sizeof (FeatCitationFormData));
  if (fcfp == NULL) 
  {
    return;
  }

  WatchCursor ();
  Update ();
  w = MovableModalWindow (-50, -33, -10, -10, "Citations", NULL);
  SetObjectExtra (w, fcfp, StdCleanupExtraProc);
  
  fcfp->ffp = ffp;
  fcfp->citation_list = FeatCitEditDialog (w, ffp->input_entityID);

  c = HiddenGroup (w, 4, 0, NULL);
  SetGroupSpacing (c, 10, 3);
  btn = PushButton (c, "Accept", AcceptFeatCit);
  SetObjectExtra (btn, fcfp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  AlignObjects (ALIGN_CENTER, (HANDLE) fcfp->citation_list, (HANDLE) c, NULL);
  RealizeWindow (w);
  
  psp = DialogToPointer (ffp->featcits);
  PointerToDialog (fcfp->citation_list, psp);
  PubSetFree (psp);

  Show (w);
  Select (w);
  ArrowCursor ();
  Update ();
}

static void ChangeGenePopupOrList (Handle gene)

{
  FeatureFormPtr  ffp;
  Int2            val;

  ffp = (FeatureFormPtr) GetObjectExtra (gene);
  if (ffp != NULL) {
    val = GetValue (ffp->gene);
    if (val == 1) {
      SafeHide (ffp->newGeneGrp);
      SafeHide (ffp->editGeneBtn);
    } else if (val == 2) {
      SafeHide (ffp->editGeneBtn);
      SafeShow (ffp->newGeneGrp);
    } else {
      SafeHide (ffp->newGeneGrp);
      SafeShow (ffp->editGeneBtn);
    }
  }
}

static void ChangeSubGroup (VoidPtr data, Int2 newval, Int2 oldval)

{
  FeatureFormPtr  ffp;

  ffp = (FeatureFormPtr) data;
  if (ffp != NULL) {
    if (ffp->commonPage >= 0 && ffp->commonPage <= 6) {
      SafeHide (ffp->commonSubGrp [ffp->commonPage]);
    }
    ffp->commonPage = newval;
    if (ffp->commonPage >= 0 && ffp->commonPage <= 6) {
      SafeShow (ffp->commonSubGrp [ffp->commonPage]);
    }
  }
}

typedef struct fieldpage {
  DIALOG_MESSAGE_BLOCK
  Int2               numfields;
  CharPtr            PNTR fields;
  TexT               PNTR values;
  ButtoN             PNTR boxes;
} FieldPage, PNTR FieldPagePtr;

static Boolean ShouldBeAGBQual (SeqFeatPtr sfp, Int2 qual, Boolean allowProductGBQual)

{
  if (qual < 0) return FALSE;
  if (allowProductGBQual && qual == GBQUAL_product) return TRUE;
  if (qual == GBQUAL_citation ||
      qual == GBQUAL_db_xref ||
      qual == GBQUAL_evidence ||
      qual == GBQUAL_exception ||
      qual == GBQUAL_gene ||
      qual == GBQUAL_label ||
      qual == GBQUAL_locus_tag ||
      qual == GBQUAL_note ||
      qual == GBQUAL_partial ||
      qual == GBQUAL_product ||
      qual == GBQUAL_pseudo ||
      qual == GBQUAL_rpt_unit ||
      qual == GBQUAL_experiment ||
      qual == GBQUAL_inference) {
    return FALSE;
  }
  if (qual == GBQUAL_map && (sfp == NULL ||
      (sfp->idx.subtype != FEATDEF_repeat_region && sfp->idx.subtype != FEATDEF_gap))) return FALSE;
  if (qual == GBQUAL_operon && (sfp == NULL || sfp->idx.subtype != FEATDEF_operon)) return FALSE;
  if (Nlm_GetAppProperty ("SequinUseEMBLFeatures") == NULL) {
    if (qual == GBQUAL_usedin) {
      return FALSE;
    }
  }
  return TRUE;
}

static CharPtr TrimParenthesesAndCommasAroundGBString (CharPtr str)

{
  Uchar    ch;	/* to use 8bit characters in multibyte languages */
  CharPtr  dst;
  CharPtr  ptr;

  if (str != NULL && str [0] != '\0') {
    dst = str;
    ptr = str;
    ch = *ptr;
    while (ch != '\0' && (ch < ' ' || ch == '(' || ch == ',')) {
      ptr++;
      ch = *ptr;
    }
    while (ch != '\0') {
      *dst = ch;
      dst++;
      ptr++;
      ch = *ptr;
    }
    *dst = '\0';
    dst = NULL;
    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      if (ch != ')' && ch != ',') {
        dst = NULL;
      } else if (dst == NULL) {
        dst = ptr;
      }
      ptr++;
      ch = *ptr;
    }
    if (dst != NULL) {
      *dst = '\0';
    }
  }
  return str;
}

static CharPtr CombineSplitGBQual (CharPtr origval, CharPtr newval)

{
  size_t   len;
  CharPtr  str = NULL;

  if (StringStr (origval, newval) != NULL) return origval;
  len = StringLen (origval) + StringLen (newval) + 5;
  str = MemNew (sizeof (Char) * len);
  if (str == NULL) return origval;
  TrimParenthesesAndCommasAroundGBString (origval);
  TrimParenthesesAndCommasAroundGBString (newval);
  StringCpy (str, "(");
  StringCat (str, origval);
  StringCat (str, ",");
  StringCat (str, newval);
  StringCat (str, ")");
  /* free original string, knowing return value will replace it */
  MemFree (origval);
  return str;
}

static void CombineTitle (TexT t, CharPtr val)

{
  CharPtr  str;

  if (t == NULL || StringHasNoText (val)) return;
  str = SaveStringFromText (t);
  if (StringDoesHaveText (str)) {
    str = CombineSplitGBQual (str, val);
    SetTitle (t, str);
  } else {
    SetTitle (t, val);
  }
  MemFree (str);
}

static void GBQualPtrToFieldPage (DialoG d, Pointer data)

{
  Char          empty [4];
  FieldPagePtr  fpf;
  Int2          i;
  GBQualPtr     list;

  fpf = (FieldPagePtr) GetObjectExtra (d);
  list = (GBQualPtr) data;
  if (fpf != NULL) {
    while (list != NULL) {
      if (list->qual != NULL && list->val != NULL) {
        for (i = 0; i < fpf->numfields; i++) {
          if (StringICmp (list->qual, fpf->fields [i]) == 0) {
            if (fpf->values [i] != NULL) {
              if (*(list->val) == '\0') {
                empty [0] = '"';
                empty [1] = '"';
                empty [2] = '\0';
                SetTitle (fpf->values [i], empty);
              } else if (StringICmp (list->qual, "rpt_type") == 0 ||
                         StringICmp (list->qual, "rpt_unit") == 0 ||
                         StringICmp (list->qual, "rpt_unit_range") == 0 ||
                         StringICmp (list->qual, "rpt_unit_seq") == 0 ||
                         StringICmp (list->qual, "replace") == 0 ||
                         StringICmp (list->qual, "compare") == 0 ||
                         StringICmp (list->qual, "old_locus_tag") == 0 ||
                         StringICmp (list->qual, "usedin") == 0) {
                CombineTitle (fpf->values [i], list->val);
              } else {
                SetTitle (fpf->values [i], list->val);
              }
            } else if (fpf->boxes [i] != NULL) {
              SetStatus (fpf->boxes [i], TRUE);
            }
          }
        }
      }
      list = list->next;
    }
  }
}

extern void SeqFeatPtrToFieldPage (DialoG d, SeqFeatPtr sfp);
extern void SeqFeatPtrToFieldPage (DialoG d, SeqFeatPtr sfp)

{
  /*
  FieldPagePtr  fpf;
  Int2          i;

  fpf = (FieldPagePtr) GetObjectExtra (d);
  if (fpf != NULL && sfp != NULL) {
        for (i = 0; i < fpf->numfields; i++) {
          if (StringICmp ("exception", fpf->fields [i]) == 0) {
            if (fpf->values [i] != NULL) {
              if (sfp->except_text == NULL || *(sfp->except_text) == '\0') {
                SetTitle (fpf->values [i], "");
              } else {
                SetTitle (fpf->values [i], sfp->except_text);
              }
            }
          } else if (StringICmp ("pseudo", fpf->fields [i]) == 0) {
            if (fpf->boxes [i] != NULL && (! GetStatus (fpf->boxes [i]))) {
              SetStatus (fpf->boxes [i], sfp->pseudo);
            }
          }
        }
  }
  */

}

static Pointer FieldPageToGBQualPtr (DialoG d)

{
  FieldPagePtr  fpf;
  GBQualPtr     gbq;
  GBQualPtr     gbqlast;
  GBQualPtr     head;
  Int2          i;

  head = NULL;
  fpf = (FieldPagePtr) GetObjectExtra (d);
  if (fpf != NULL) {
    gbq = NULL;
    gbqlast = NULL;
    for (i = 0; i < fpf->numfields; i++) {
      if (fpf->fields [i] == NULL ||
          StringHasNoText (fpf->fields [i]) ||
          (fpf->values [i] == NULL && fpf->boxes [i] ==NULL) ||
          (fpf->values [i] != NULL && TextHasNoText (fpf->values [i]))) {
      } else if (fpf->boxes [i] != NULL && (! GetStatus (fpf->boxes [i]))) {
      } else {
        gbq = GBQualNew ();
        if (gbqlast == NULL) {
          head = gbq;
        } else {
          gbqlast->next = gbq;
        }
        gbqlast = gbq;
        if (gbq != NULL) {
          gbq->qual = StringSave (fpf->fields [i]);
          if (fpf->values [i] != NULL) {
            gbq->val = SaveStringFromText (fpf->values [i]);
          } else {
            gbq->val = StringSave ("");
          }
        }
      }
    }
  }
  return (Pointer) head;
}

static void CleanupFieldsPage (GraphiC g, VoidPtr data)

{
  FieldPagePtr  fpf;
  Int2          i;

  fpf = (FieldPagePtr) data;
  if (fpf != NULL) {
    if (fpf->fields != NULL) {
      for (i = 0; i < fpf->numfields; i++) {
        MemFree (fpf->fields [i]);
      }
    }
    MemFree (fpf->fields);
    MemFree (fpf->values);
    MemFree (fpf->boxes);
  }
  MemFree (data);
}

#define LEGAL_FEATURE    1
#define ILLEGAL_FEATURE  2

extern DialoG CreateImportFields (GrouP h, CharPtr name, SeqFeatPtr sfp, Boolean allowProductGBQual)

{
  FieldPagePtr    fpf;
  GrouP           g;
  GBQualPtr       gbq;
  Boolean         hasillegal;
  Int2            i;
  ImpFeatPtr      imp;
  Int2            index;
  Boolean         is_gap = FALSE;
  Int2            j;
  Int2            max;
  Int2            num;
  GrouP           p;
  PrompT          ppt;
  Int2            qual;
  Int2            seen [ParFlat_TOTAL_GBQUAL];
  SematicFeatPtr  sefp;
  Char            str [64];
  Int2            wid;

  p = HiddenGroup (h, 1, 0, NULL);
  fpf = (FieldPagePtr) MemNew (sizeof (FieldPage));
  if (fpf != NULL) {

    SetObjectExtra (p, fpf, CleanupFieldsPage);
    fpf->dialog = (DialoG) p;
    fpf->todialog = GBQualPtrToFieldPage;
    fpf->fromdialog = FieldPageToGBQualPtr;
    fpf->testdialog = NULL;

    if (sfp != NULL && sfp->data.choice == SEQFEAT_IMP) {
      imp = (ImpFeatPtr) sfp->data.value.ptrvalue;
      if (imp != NULL && StringICmp (imp->key, "gap") == 0) {
        is_gap = TRUE;
      }
    }

    gbq = NULL;
    if (sfp != NULL) {
      gbq = sfp->qual;
    }

    j = 0;
    hasillegal = FALSE;
    for (i = 0; i < ParFlat_TOTAL_GBQUAL; i++) {
      seen [i] = 0;
    }
    num = 0;

    if (name != NULL) {
      index = GBFeatKeyNameValid (&name, FALSE);
      if (index >= 0) {

        sefp = &(ParFlat_GBFeat [index]);

        for (i = 0; i < sefp->mand_num; i++) {
          qual = sefp->mand_qual [i];
          if (qual > -1 && ShouldBeAGBQual (sfp, qual, allowProductGBQual)) {
            seen [qual] = LEGAL_FEATURE;
          }
        }
        if (StringCmp (name, "repeat_region") == 0) {
          seen [GBQUAL_map] = TRUE;
        }
        if (StringCmp (name, "operon") == 0) {
          seen [GBQUAL_operon] = TRUE;
        }
        for (i = 0; i < sefp->opt_num; i++) {
          qual = sefp->opt_qual [i];
          if (qual > -1 && ShouldBeAGBQual (sfp, qual, allowProductGBQual)) {
            seen [qual] = LEGAL_FEATURE;
          }
        }
      }
    } else if (sfp != NULL && sfp->data.choice == SEQFEAT_CDREGION) {
      /*
      seen [GBQUAL_exception] = LEGAL_FEATURE;
      seen [GBQUAL_pseudo] = LEGAL_FEATURE;
      */
    }

    while (gbq != NULL) {
      qual = GBQualNameValid (gbq->qual);
      if (qual > -1) {
        if (seen [qual] == 0 && qual != GBQUAL_experiment && qual != GBQUAL_inference) {
          seen [qual] = ILLEGAL_FEATURE;
          hasillegal = TRUE;
        }
      }
      gbq = gbq->next;
    }

    SelectFont (programFont);
    max = 0;
    for (i = 0; i < ParFlat_TOTAL_GBQUAL; i++) {
      if (seen [i] > 0) {
        num++;
        if (seen [i] == LEGAL_FEATURE) {
          StringNCpy_0 (str, ParFlat_GBQual_names [i].name, sizeof (str));
          wid = StringWidth (str) + 2;
        } else {
          str [0] = '\0';
          if (name != NULL) {
            StringCpy (str, "*");
          }
          StringNCat (str, ParFlat_GBQual_names [i].name, sizeof (str) - 2);
          wid = StringWidth (str) + 2;
        }
        if (wid > max) {
          max = wid;
        }
      }
    }
    SelectFont (systemFont);

    fpf->numfields = num;
    fpf->fields = MemNew (sizeof (CharPtr) * (num + 1));
    fpf->values = MemNew (sizeof (TexT) * (num + 1));
    fpf->boxes = MemNew (sizeof (ButtoN) * (num + 1));

    g = HiddenGroup (p, 2, 0, NULL);
    j = 0;
    for (i = 0; i < ParFlat_TOTAL_GBQUAL; i++) {
      if (seen [i] == LEGAL_FEATURE) {
        fpf->fields [j] = StringSave (ParFlat_GBQual_names [i].name);
        ppt = StaticPrompt (g, fpf->fields [j], max, dialogTextHeight, programFont, 'l');
        if (ParFlat_GBQual_names [i].gbclass != Class_none) {
          fpf->values [j] = DialogText (g, "", 20, NULL);
          /*
          if (i == GBQUAL_estimated_length && is_gap) {
            Disable (fpf->values [j]);
          }
          */
        } else {
          fpf->boxes [j] = CheckBox (g, "", NULL);
          AlignObjects (ALIGN_MIDDLE, (HANDLE) ppt, (HANDLE) fpf->boxes [j], NULL);
        }
        j++;
      }
    }
    if (hasillegal && name != NULL) {
      StaticPrompt (p, "Illegal Qualifiers", 0, 0, programFont, 'c');
    }
    g = HiddenGroup (p, 2, 0, NULL);
    for (i = 0; i < ParFlat_TOTAL_GBQUAL; i++) {
      if (seen [i] == ILLEGAL_FEATURE) {
        fpf->fields [j] = StringSave (ParFlat_GBQual_names [i].name);
        str [0] = '\0';
        if (name != NULL) {
          StringCpy (str, "*");
        }
        StringNCat (str, fpf->fields [j], sizeof (str) - 2);
        ppt = StaticPrompt (g, str, max, dialogTextHeight, programFont, 'l');
        if (ParFlat_GBQual_names [i].gbclass != Class_none) {
          fpf->values [j] = DialogText (g, "", 20, NULL);
        } else {
          fpf->boxes [j] = CheckBox (g, "", NULL);
          AlignObjects (ALIGN_MIDDLE, (HANDLE) ppt, (HANDLE) fpf->boxes [j], NULL);
        }
        j++;
      }
    }

    if (j == 0) {
      StaticPrompt (p, "See Attributes page to set legal qualifiers for this feature.",
                    0, 0, programFont, 'c');
    }
  }
  return (DialoG) p;
}

static void ChangeCannedMessage (PopuP p)

{
  FeatureFormPtr  ffp;
  Int2            val;

  ffp = (FeatureFormPtr) GetObjectExtra (p);
  if (ffp == NULL) return;
  val = GetValue (p);
  switch (val) {
    case 1 :
      if (Message (MSG_YN, "Clear the explanation field?") == ANS_YES) {
        SetTitle (ffp->exceptText, "");
        if (Message (MSG_YN, "Clear the exception flag?") == ANS_YES) {
          SetStatus (ffp->exception, FALSE);
        }
      }
      break;
    case 2 :
      SetTitle (ffp->exceptText, "RNA editing");
      SetStatus (ffp->exception, TRUE);
      break;
    case 3 :
      SetTitle (ffp->exceptText, "reasons given in citation");
      SetStatus (ffp->exception, TRUE);
      break;
    case 4 :
      SetTitle (ffp->exceptText, "ribosomal slippage");
      SetStatus (ffp->exception, TRUE);
      break;
    case 5 :
      SetTitle (ffp->exceptText, "trans-splicing");
      SetStatus (ffp->exception, TRUE);
      break;
    case 6 :
      SetTitle (ffp->exceptText, "artificial frameshift");
      SetStatus (ffp->exception, TRUE);
      break;
    case 7 :
      SetTitle (ffp->exceptText, "nonconsensus splice site");
      SetStatus (ffp->exception, TRUE);
      break;
    case 8 :
      SetTitle (ffp->exceptText, "rearrangement required for product");
      SetStatus (ffp->exception, TRUE);
      break;
    case 9 :
      SetTitle (ffp->exceptText, "alternative start codon");
      SetStatus (ffp->exception, TRUE);
      break;
    default :
      break;
  }
}

static CharPtr crossRefWarn =
"A gene mapped by cross-reference does not extend the range\n\
of the indicated gene feature.  Overlap is the usual case, and\n\
it does extend the selected gene.";

static CharPtr suppressWarn =
"This will suppress display of a gene qualifier even though\n\
there is a gene feature that overlaps this one.";

static void GeneXrefWarn (GrouP g)

{
  Int2  val;
  Boolean indexerVersion;

  indexerVersion = (Boolean) (GetAppProperty ("InternalNcbiSequin") != NULL);
  if (indexerVersion) return;

  val = GetValue (g);
  if (val == 2) {
    Message (MSG_OK, "%s", crossRefWarn);
  } else if (val == 3) {
    Message (MSG_OK, "%s", suppressWarn);
  }
}

static CharPtr  commonRadioFormTabs [] = {
  "General", "Comment", "Citations", "Cross-Refs", "Evidence", "Identifiers", NULL, NULL
};

static CharPtr  commonNoCitFormTabs [] = {
  "General", "Comment", "Cross-Refs", "Evidence", "Identifiers", NULL, NULL
};

static DialoG CreateInferenceDialog (GrouP h, Uint2 rows, Int2 spacing, Int2 width);
static DialoG NewCreateInferenceDialog (GrouP prnt);
extern void Nlm_LaunchGeneFeatEd (ButtoN b);

extern GrouP CreateCommonFeatureGroupEx (GrouP h, FeatureFormPtr ffp,
                                         SeqFeatPtr sfp, Boolean hasGeneControl,
                                         Boolean hasCitationTab, Boolean hasGeneSuppress)

{
  ButtoN     b;
  GrouP      c;
  PopuP      canned;
  Boolean    cdsQuals;
  GrouP      f;
  GrouP      g;
  GBQualPtr  gbq;
  Boolean    hasQuals;
  Boolean    indexerVersion;
  Char       just;
  GrouP      k;
  GrouP      m;
  GrouP      p;
  Int2       page;
  PrompT     ppt1, ppt2;
  ButtoN     pseudo = NULL;
  GrouP      q;
  GrouP      r;
  GrouP      t;
  GrouP      v;
  GrouP      x;
  GrouP      y;

  c = NULL;
  if (ffp != NULL) {
    hasQuals = FALSE;
    cdsQuals = FALSE;
    if (ffp->gbquals == NULL && sfp != NULL && sfp->qual != NULL) {
      for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
        if (StringCmp (gbq->qual, "experiment") != 0 && StringCmp (gbq->qual, "inference") != 0) {
          hasQuals = TRUE;
        }
      }
      /*
      if (GetAppProperty ("InternalNcbiSequin") != NULL) {
        if (sfp->data.choice == SEQFEAT_CDREGION) {
          cdsQuals = TRUE;
        }
      }
      */
    }
    m = HiddenGroup (h, -1, 0, NULL);
    SetGroupSpacing (m, 10, 10);
    ffp->commonPage = 0;
    if (cdsQuals) {
    } else if (hasQuals) {
      commonRadioFormTabs [6] = "Qualifiers";
      commonNoCitFormTabs [5] = "Qualifiers";
    }
    if (hasCitationTab) {
      ffp->commonRadio = CreateFolderTabs (m, commonRadioFormTabs, ffp->commonPage,
                                           0, 0, PROGRAM_FOLDER_TAB,
                                           ChangeSubGroup, (Pointer) ffp);
    } else {
      ffp->commonRadio = CreateFolderTabs (m, commonNoCitFormTabs, ffp->commonPage,
                                           0, 0, PROGRAM_FOLDER_TAB,
                                           ChangeSubGroup, (Pointer) ffp);
    }
    commonRadioFormTabs [6] = NULL;
    commonNoCitFormTabs [5] = NULL;
    p = HiddenGroup (m, 0, 0, NULL);

    c = HiddenGroup (p, -1, 0, NULL);
    SetGroupSpacing (c, 10, 10);
    f = HiddenGroup (c, -1, 0, NULL);
    r = HiddenGroup (f, 7, 0, NULL);
    StaticPrompt (r, "Flags", 0, popupMenuHeight, programFont, 'l');
    ffp->partial = CheckBox (r, "Partial", NULL);
    indexerVersion = (Boolean) (GetAppProperty ("InternalNcbiSequin") != NULL);
    if (! indexerVersion) {
      Disable (ffp->partial);
    }
    if (ffp->pseudo == NULL) {
      ffp->pseudo = CheckBox (r, "Pseudo", NULL);
      pseudo = ffp->pseudo; /* allows pseudo control on earlier feature-specific page */
    }
    StaticPrompt (r, "Evidence", 0, popupMenuHeight, programFont, 'l');
    ffp->evidence = PopupList (r, TRUE, NULL);
    PopupItem (ffp->evidence, " ");
    PopupItem (ffp->evidence, "Experimental");
    PopupItem (ffp->evidence, "Non-Experimental");
    AlignObjects (ALIGN_MIDDLE, (HANDLE) ffp->partial,
                  (HANDLE) ffp->evidence, (HANDLE) pseudo, NULL);
    r = HiddenGroup (f, -3, 0, NULL);
    ffp->exception = CheckBox (r, "Exception", NULL);
    StaticPrompt (r, "Explanation", 0, dialogTextHeight, programFont, 'l');
    ffp->exceptText = DialogText (r, "", 12, NULL);
    AlignObjects (ALIGN_MIDDLE, (HANDLE) ffp->exception,
                 (HANDLE) ffp->exceptText, NULL);
    if (ffp->this_subtype == FEATDEF_CDS) {
      StaticPrompt (r, "Standard explanation", 0, popupMenuHeight, programFont, 'l');
      canned = PopupList (r, TRUE, ChangeCannedMessage);
      SetObjectExtra (canned, (Pointer) ffp, NULL);
      PopupItem (canned, " ");
      PopupItem (canned, "RNA editing");
      PopupItem (canned, "reasons given in citation");
      PopupItem (canned, "ribosomal slippage");
      PopupItem (canned, "trans-splicing");
      PopupItem (canned, "artificial frameshift");
      PopupItem (canned, "nonconsensus splice site");
      PopupItem (canned, "rearrangement required");
      PopupItem (canned, "alternative start codon");
      if (sfp != NULL && sfp->excpt) {
        if (StringICmp (sfp->except_text, "RNA editing") == 0) {
          SetValue (canned, 2);
        } else if (StringICmp (sfp->except_text, "reasons given in citation") == 0 ||
                   StringICmp (sfp->except_text, "reasons cited in publication") == 0) {
          SetValue (canned, 3);
        } else if (StringICmp (sfp->except_text, "ribosomal slippage") == 0 ||
                   StringICmp (sfp->except_text, "ribosome slippage") == 0) {
          SetValue (canned, 4);
        } else if (StringICmp (sfp->except_text, "trans-splicing") == 0 ||
                   StringICmp (sfp->except_text, "trans splicing") == 0) {
          SetValue (canned, 5);
        } else if (StringICmp (sfp->except_text, "artificial frameshift") == 0) {
          SetValue (canned, 6);
        } else if (StringICmp (sfp->except_text, "non-consensus splice site") == 0 ||
                   StringICmp (sfp->except_text, "nonconsensus splice site") == 0) {
          SetValue (canned, 7);
        } else if (StringICmp (sfp->except_text, "rearrangement required for product") == 0) {
          SetValue (canned, 8);
        } else if (StringICmp (sfp->except_text, "alternative start codon") == 0) {
          SetValue (canned, 9);
        }
      } else {
        SetValue (canned, 1);
      }
    }

    if (cdsQuals) {
      /*
      ffp->gbquals = CreateImportFields (c, NULL, sfp, FALSE);
      */
    }

    g = NULL;
    k = NULL;
    ffp->gene = NULL;
    ffp->genePopup = NULL;
    ffp->geneList = NULL;
    ffp->geneNames = NULL;
    ffp->useGeneXref = NULL;
    ffp->newGeneGrp = NULL;
    ffp->geneSymbol = NULL;
    ffp->geneAllele = NULL;
    ffp->geneDesc = NULL;
    ffp->locusTag = NULL;
    for (page = 0; page < 8; page++) {
      ffp->commonSubGrp [page] = NULL;
    }
    page = 0;

    if (hasGeneControl) {
      g = HiddenGroup (c, -2, 0, NULL);
      StaticPrompt (g, "Gene", 0, popupMenuHeight, programFont, 'l');
      v = HiddenGroup (g, 0, 0, NULL);
      ffp->genePopup = PopupList (v, TRUE, (PupActnProc) ChangeGenePopupOrList);
      SetObjectExtra (ffp->genePopup, (Pointer) ffp, NULL);
      PopupItem (ffp->genePopup, " ");
      PopupItem (ffp->genePopup, "New:");
      SetValue (ffp->genePopup, 1);
      Hide (ffp->genePopup);
      ffp->geneList = SingleList (v, 6, 3, (LstActnProc) ChangeGenePopupOrList);
      SetObjectExtra (ffp->geneList, (Pointer) ffp, NULL);
      ListItem (ffp->geneList, " ");
      ListItem (ffp->geneList, "New:");
      SetValue (ffp->geneList, 1);
      Hide (ffp->geneList);
      k = HiddenGroup (c, 3, 0, NULL);
      StaticPrompt (k, "Map by", 0, stdLineHeight, programFont, 'l');
      ffp->useGeneXref = HiddenGroup (k, 3, 0, GeneXrefWarn);
      SetObjectExtra (ffp->useGeneXref, ffp, NULL);
      RadioButton (ffp->useGeneXref, "Overlap");
      RadioButton (ffp->useGeneXref, "Cross-reference");
      RadioButton (ffp->useGeneXref, "Suppress");
      SetValue (ffp->useGeneXref, 1);
      y = HiddenGroup (c, 0, 0, NULL);
      ffp->newGeneGrp = HiddenGroup (y, 2, 0, NULL);
      StaticPrompt (ffp->newGeneGrp, "Gene Symbol", 0, dialogTextHeight, programFont, 'l');
      ffp->geneSymbol = DialogText (ffp->newGeneGrp, "", 20, NULL);
      StaticPrompt (ffp->newGeneGrp, "Allele", 0, dialogTextHeight, programFont, 'l');
      ffp->geneAllele = DialogText (ffp->newGeneGrp, "", 20, NULL);
      StaticPrompt (ffp->newGeneGrp, "Description", 0, dialogTextHeight, programFont, 'l');
      ffp->geneDesc = DialogText (ffp->newGeneGrp, "", 20, NULL);
      StaticPrompt (ffp->newGeneGrp, "Locus Tag", 0, dialogTextHeight, programFont, 'l');
      ffp->locusTag = DialogText (ffp->newGeneGrp, "", 20, NULL);
      Hide (ffp->newGeneGrp);
      ffp->editGeneBtn = PushButton (y, "Edit Gene Feature", Nlm_LaunchGeneFeatEd);
      SetObjectExtra (ffp->editGeneBtn, ffp, NULL);
      Hide (ffp->editGeneBtn);
    } else if (hasGeneSuppress) {
      k = HiddenGroup (c, 3, 0, NULL);
      StaticPrompt (k, "Map by", 0, stdLineHeight, programFont, 'l');
      ffp->useGeneXref = HiddenGroup (k, 3, 0, GeneXrefWarn);
      SetObjectExtra (ffp->useGeneXref, ffp, NULL);
      RadioButton (ffp->useGeneXref, "Overlap");
      b = RadioButton (ffp->useGeneXref, "Cross-reference");
      Disable (b);
      RadioButton (ffp->useGeneXref, "Suppress");
      SetValue (ffp->useGeneXref, 1);
      y = HiddenGroup (c, 0, 0, NULL);
    }
    ffp->commonSubGrp [page] = c;
    page++;
    AlignObjects (ALIGN_CENTER, (HANDLE) f, (HANDLE) g,
                  (HANDLE) k, (HANDLE) ffp->newGeneGrp,
                  (HANDLE) ffp->editGeneBtn, NULL);

    c = HiddenGroup (p, -1, 0, NULL);
    SetGroupSpacing (c, 10, 10);
    q = HiddenGroup (c, 0, 2, NULL);
    StaticPrompt (q, "Comment", 25 * stdCharWidth, 0, programFont, 'c');
    if (GetAppProperty ("InternalNcbiSequin") != NULL) {
      ffp->comment = ScrollText (q, 25, 10, programFont, TRUE, NULL);
    } else {
      ffp->comment = ScrollText (q, 25, 5, programFont, TRUE, NULL);
    }
    ffp->commonSubGrp [page] = c;
    Hide (ffp->commonSubGrp [page]);
    page++;

    if (hasCitationTab) {
      c = HiddenGroup (p, -1, 0, NULL);
      SetGroupSpacing (c, 10, 10);
      ffp->featcits = CreateCitOnFeatDialog (c, NULL);
      b = PushButton (c, "Edit Citations", EditFeatCitsProc);
      SetObjectExtra (b, ffp, NULL);
      ffp->commonSubGrp [page] = c;
      AlignObjects (ALIGN_CENTER, (HANDLE) ffp->featcits, (HANDLE) b, NULL);
      Hide (ffp->commonSubGrp [page]);
      page++;
    }

    c = HiddenGroup (p, -1, 0, NULL);
    SetGroupSpacing (c, 10, 10);
    if (GetAppProperty ("ReadOnlyDbTags") == NULL) {
      just = 'c';
    } else {
      just = 'l';
      StaticPrompt (c, "This page is read-only", 15 * stdCharWidth, 0, programFont, 'c');
    }
    t = HiddenGroup (c, 2, 0, NULL);
    StaticPrompt (t, "Database", 7 * stdCharWidth, 0, programFont, just);
    StaticPrompt (t, "Object ID", 8 * stdCharWidth, 0, programFont, just);
    ffp->dbxrefs = CreateDbtagDialog (c, 3, -1, 7, 8);
    ffp->commonSubGrp [page] = c;
    Hide (ffp->commonSubGrp [page]);
    page++;

    c = HiddenGroup (p, -1, 0, NULL);
    SetGroupSpacing (c, 10, 10);
    q = HiddenGroup (c, 0, -6, NULL);
    ppt1 = StaticPrompt (q, "Experiment", 0, 0, programFont, 'c');
    ffp->experiment = CreateVisibleStringDialog (q, 3, -1, 15);
    ppt2 = StaticPrompt (q, "Inference", 0, 0, programFont, 'c');
    /*
    ffp->inference = CreateInferenceDialog (q, 3, 2, 15);
    */
    ffp->inference = NewCreateInferenceDialog (q);
    AlignObjects (ALIGN_CENTER, (HANDLE) ppt1, (HANDLE) ffp->experiment,
                  (HANDLE) ppt2, (HANDLE) ffp->inference, NULL);
    ffp->commonSubGrp [page] = c;
    Hide (ffp->commonSubGrp [page]);
    page++;

    c = HiddenGroup (p, -1, 0, NULL);
    SetGroupSpacing (c, 10, 10);
    t = HiddenGroup (c, 2, 0, NULL);
    SetGroupSpacing (t, 10, 30);
    StaticPrompt (t, "Feature ID for this feature", 0, dialogTextHeight, programFont, 'l');
    ffp->featid = DialogText (t, "", 8, NULL);
    StaticPrompt (t, "ID Xref to associated feature", 0, dialogTextHeight, programFont, 'l');
    ffp->fidxref = DialogText (t, "", 8, NULL);
    if (! indexerVersion) {
      Disable (ffp->featid);
      Disable (ffp->fidxref);
    }
    ffp->commonSubGrp [page] = c;
    Hide (ffp->commonSubGrp [page]);
    page++;

    c = HiddenGroup (p, -1, 0, NULL);
    SetGroupSpacing (c, 10, 10);
    x = NULL;
    /* z = NULL; */
    if (hasQuals && ffp->gbquals == NULL) {
      x = HiddenGroup (c, -1, 0, NULL);
      /*
      z = HiddenGroup (x, -4, 0, NULL);
      SetGroupSpacing (z, -1, 0);
      StaticPrompt (z, "Qualifier", 7 * stdCharWidth, 0, programFont, 'c');
      StaticPrompt (z, "Value", 10 * stdCharWidth, 0, programFont, 'c');
      ffp->gbquals = CreateQualsDialog (x, 5, -1, 7, 10);
      */
      ffp->gbquals = CreateImportFields (x, NULL, sfp, FALSE);
    }
    ffp->commonSubGrp [page] = c;
    Hide (ffp->commonSubGrp [page]);
    page++;

    AlignObjects (ALIGN_CENTER, (HANDLE) ffp->commonRadio, (HANDLE) ffp->commonSubGrp [0],
                  (HANDLE) ffp->commonSubGrp [1], (HANDLE) ffp->commonSubGrp [2],
                  (HANDLE) ffp->commonSubGrp [3], (HANDLE) ffp->commonSubGrp [4],
                  (HANDLE) ffp->commonSubGrp [5], (HANDLE) ffp->commonSubGrp [6],
                  (HANDLE) ffp->gbquals, NULL);
  }
  return c;
}

extern GrouP CreateCommonFeatureGroup (GrouP h, FeatureFormPtr ffp,
                                       SeqFeatPtr sfp, Boolean hasGeneControl,
                                       Boolean hasCitationTab)

{
  return CreateCommonFeatureGroupEx (h, ffp, sfp, hasGeneControl, hasCitationTab, FALSE);
}

static Boolean DlgutilFindBspItem (GatherContextPtr gcp)

{
  BioseqPtr  PNTR bspp;

  bspp = (BioseqPtr PNTR) gcp->userdata;
  if (bspp != NULL && gcp->thistype == OBJ_BIOSEQ) {
    *bspp = (BioseqPtr) gcp->thisitem;
  }
  return TRUE;
}

extern void SetNewFeatureDefaultInterval (FeatureFormPtr ffp)

{
  BioseqPtr     bsp;
  SelStructPtr  sel;
  SeqIntPtr     sip;
  SeqLocPtr     slp;

  if (ffp == NULL) return;
  bsp = NULL;
  GatherItem (ffp->input_entityID, ffp->input_itemID, ffp->input_itemtype,
              (Pointer) (&bsp), DlgutilFindBspItem);
  if (bsp == NULL) return;
  slp = NULL;
  sel = ObjMgrGetSelected ();
  if (sel != NULL && sel->next == NULL && sel->entityID == ffp->input_entityID &&
      sel->itemID == ffp->input_itemID && sel->itemtype == ffp->input_itemtype) {
    if (sel->regiontype == 1 && sel->region != NULL) {
      if (GetBioseqGivenSeqLoc ((SeqLocPtr) sel->region, ffp->input_entityID) == bsp) {
        slp = AsnIoMemCopy (sel->region,
                            (AsnReadFunc) SeqLocAsnRead,
                            (AsnWriteFunc) SeqLocAsnWrite);
      }
    }
  }
  if (slp == NULL) {
    slp = ValNodeNew (NULL);
    if (slp == NULL) return;
    sip = SeqIntNew ();
    if (sip == NULL) {
      slp = SeqLocFree (slp);
      return;
    }
    slp->choice = SEQLOC_INT;
    slp->data.ptrvalue = (Pointer) sip;
    sip->from = 0;
    sip->to = bsp->length - 1;
    sip->strand = Seq_strand_plus;
    sip->id = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (bsp->id, 0)));
  }
  if (slp != NULL) {
    PointerToDialog (ffp->location, (Pointer) slp);
    SeqLocFree (slp);
  }
}

extern Boolean FileToScrollText (TexT t, CharPtr path)

{
  FILE     *fp;
  Int4     len;
  Int4     read_len;
  Int4     max;
  CharPtr  str;
#ifdef WIN_MAC
  CharPtr  p;
#endif
#if (defined(OS_DOS) || defined (OS_NT))
  CharPtr  p;
  CharPtr  q;
#endif

  if (t != NULL && path != NULL && *path != '\0') {
    len = FileLength (path);
    max = (Int4) INT2_MAX;
#ifdef WIN_MOTIF
    max = INT4_MAX;
#endif
#ifdef WIN_MSWIN
    max = INT4_MAX;
#endif
#ifdef WIN_MAC
#ifdef OS_UNIX_DARWIN
    max = INT4_MAX;
#endif
#endif
    if (len > 0 && len < max - 4) {
      str = MemNew (sizeof (char) * (len + 3));
      if (str != NULL) {
        fp = FileOpen (path, "r");
        if (fp != NULL) {
          read_len = FileRead (str, sizeof (char), (size_t) len, fp);
          str [ read_len ] = 0;
#if (defined(OS_DOS) || defined (OS_NT))
          p = str;
          q = str;
          while (*p) {
            if (*p == '\r') {
              p++;
            } else {
              *q = *p;
              p++;
              q++;
            }
          }
          *q = '\0';
#endif
#ifdef WIN_MAC
          p = str;
          while (*p) {
            if (*p == '\r' || *p == '\n') {
              *p = '\015';
            }
            p++;
          }
#endif
          FileClose (fp);
          SetTitle (t, str);
        }
        MemFree (str);
        return TRUE;
      }
    }
  }
  return FALSE;
}

extern void ScrollTextToFile (TexT t, CharPtr path)

{
  FILE     *fp;
  size_t   len;
  CharPtr  str;
#ifdef WIN_MAC
  CharPtr  p;
#endif

  if (t != NULL && path != NULL && *path != '\0') {
    len = TextLength (t);
    if (len > 0) {
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
        str = MemNew (sizeof (char) * (len + 3));
        if (str != NULL) {
          GetTitle (t, str, len + 1);
#ifdef WIN_MAC
          p = str;
          while (*p) {
            if (*p == '\r' || *p == '\n') {
              *p = '\n';
            }
            p++;
          }
#endif
          FileWrite (str, sizeof (char), len, fp);
          MemFree (str);
        }
        FileClose (fp);
      }
    }
  }
}

extern void FileToClipboard (CharPtr path)

{
  FILE     *fp;
  Int4     len;
  Int4     max;
  CharPtr  str;
#ifdef WIN_MAC
  CharPtr  p;
#endif
#if (defined(OS_DOS) || defined (OS_NT))
  CharPtr  p;
  CharPtr  q;
#endif

  if (path != NULL && *path != '\0') {
    len = FileLength (path);
#ifdef WIN_MOTIF
    max = INT4_MAX;
#else
    max = (Int4) INT2_MAX;
#endif
    if (len > 0 && len < max - 4) {
      str = MemNew (sizeof (char) * (len + 3));
      if (str != NULL) {
        fp = FileOpen (path, "r");
        if (fp != NULL) {
          FileRead (str, sizeof (char), (size_t) len, fp);
#if (defined(OS_DOS) || defined (OS_NT))
          p = str;
          q = str;
          while (*p) {
            if (*p == '\r') {
              p++;
            } else {
              *q = *p;
              p++;
              q++;
            }
          }
          *q = '\0';
#endif
#ifdef WIN_MAC
          p = str;
          while (*p) {
            if (*p == '\r' || *p == '\n') {
              *p = '\015';
            }
            p++;
          }
#endif
          FileClose (fp);
          Nlm_StringToClipboard (str);
        }
        MemFree (str);
      }
    }
  }
}

typedef struct textviewform {
  FORM_MESSAGE_BLOCK

  DoC              doc;
  TexT             text;

  TexT             find;
} TextViewForm, PNTR TextViewFormPtr;

static ParData txtParFmt = {FALSE, FALSE, FALSE, FALSE, TRUE, 0, 0};
static ColData txtColFmt = {0, 0, 80, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, TRUE};

static void ResizeTextViewer (WindoW w)

{
  Int2             height;
  RecT             r;
  RecT             s;
  TextViewFormPtr  tfp;
  Int2             width;

  tfp = (TextViewFormPtr) GetObjectExtra (w);
  if (tfp != NULL) {
    WatchCursor ();
    ObjectRect (w, &r);
    width = r.right - r.left;
    height = r.bottom - r.top;
    if (tfp->doc != NULL) {
      GetPosition (tfp->doc, &s);
      s.right = width - s.left;
      s.bottom = height - s.left;
      SetPosition (tfp->doc, &s);
      AdjustPrnt (tfp->doc, &s, FALSE);
      txtColFmt.pixWidth = screenRect.right - screenRect.left;
      txtColFmt.pixInset = 8;
      if (Visible (tfp->doc) && AllParentsVisible (tfp->doc)) {
        UpdateDocument (tfp->doc, 0, 0);
      }
    }
    if (tfp->text != NULL) {
      GetPosition (tfp->text, &s);
      s.right = width - s.left;
      s.bottom = height - s.left;
      SetPosition (tfp->text, &s);
      AdjustPrnt (tfp->text, &s, FALSE);
    }
    ArrowCursor ();
    Update ();
  }
}

static Boolean ExportTextViewForm (ForM f, CharPtr filename)

{
  FILE             *fp;
  Char             path [PATH_MAX];
  TextViewFormPtr  tfp;

  tfp = (TextViewFormPtr) GetObjectExtra (f);
  if (tfp != NULL && (tfp->doc != NULL || tfp->text != NULL)) {
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
      if (tfp->doc != NULL) {
        fp = FileOpen (path, "w");
        if (fp != NULL) {
          SaveDocument (tfp->doc, fp);
          FileClose (fp);
          return TRUE;
        }
      } else if (tfp->text != NULL) {
        ScrollTextToFile (tfp->text, path);
        return TRUE;
      }
    }
  }
  return FALSE;
}

static void LIBCALL PrintTextViewForm (Pointer formDataPtr)

{
  TextViewFormPtr  tfp;

  tfp = (TextViewFormPtr) formDataPtr;
  if (tfp != NULL && tfp->doc != NULL) {
    PrintDocument (tfp->doc);
  }
}

static void TextViewFormMessage (ForM f, Int2 mssg)

{
  TextViewFormPtr  tfp;

  tfp = (TextViewFormPtr) GetObjectExtra (f);
  if (tfp != NULL) {
    switch (mssg) {
      case VIB_MSG_EXPORT :
        ExportTextViewForm (f, NULL);
        break;
      case VIB_MSG_CLOSE :
        Remove (f);
        break;
      case VIB_MSG_PRINT :
        PrintTextViewForm (tfp);
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
        if (tfp->appmessage != NULL) {
          tfp->appmessage (f, mssg);
        }
        break;
    }
  }
}


static void FindInGeneralText (ButtoN b)

{
  Int2             actual;
  Char             buf [1030];
  Char             ch;
  Int2             cnt;
  Int4             cntr;
  Int2             first;
  FILE             *fp;
  Char             lastch;
  Int4             line;
  ValNodePtr       matches;
  Int2             next;
  Int2             offset = 0;
  Char             path [PATH_MAX];
  CharPtr          ptr;
  Int2             state;
  CharPtr          str;
  TextFsaPtr       tbl;
  TextViewFormPtr  tfp;
  Int4             max;

  tfp = (TextViewFormPtr) GetObjectExtra (b);
  if (tfp == NULL) return;
  if (tfp->doc != NULL) {
    GetOffset (tfp->doc, NULL, &offset);
  } else if (tfp->text != NULL) {
    GetOffset (tfp->text, NULL, &offset);
  }
  first = -1;
  next = -1;
  max = INT2_MAX;

  str = SaveStringFromText (tfp->find);

  if (StringDoesHaveText (str)) {
    TmpNam (path);
    if (ExportForm (tfp->form, path)) {
      tbl = TextFsaNew ();
      if (tbl != NULL) {
        TextFsaAdd (tbl, str);
        fp = FileOpen (path, "r");
        if (fp != NULL) {
          line = 0;
          state = 0;
          cntr = FileLength (path);
          cnt = (Int2) MIN (cntr, 1024L);
          lastch = '\0';
          while (cnt > 0 && cntr > 0 && line <= max && 
                 ((next == -1 && offset > first) || first == -1)) {
            actual = (Int2) FileRead (buf, 1, cnt, fp);
            if (actual > 0) {
              cnt = actual;
              buf [cnt] = '\0';
              ptr = buf;
              ch = *ptr;
              while (ch != '\0') {
                if (ch == '\n' || ch == '\r') {
                  if (ch == '\n' && lastch == '\r') {
                    /* do not increment line */
                  } else if (ch == '\r' && lastch == '\n') {
                    /* do not increment line */
                  } else {
                    line++;
                  }
                }
                state = TextFsaNext (tbl, state, ch, &matches);
                if (matches != NULL) {
                  if (first == -1) {
                    first = line;
                  }
                  if (next == -1 && line > offset) {
                    next = line;
                  }
                }
                lastch = ch;
                ptr++;
                ch = *ptr;
              }
              cntr -= cnt;
              cnt = (Int2) MIN (cntr, 1024L);
            } else {
              cnt = 0;
              cntr = 0;
            }
          }
        }
        FileClose (fp);
      }
      TextFsaFree (tbl);
    }
    FileRemove (path);
  }
  MemFree (str);
  if (line > max) {
    Message (MSG_ERROR, "Too many lines for search");
  }

  if (next >= 0) {
    offset = next;
  } else if (first >= 0) {
    offset = first;
  } else return;
  if (tfp->doc != NULL) {
    SetOffset (tfp->doc, 0, offset);
    Update ();
  } else if (tfp->text != NULL) {
    SetOffset (tfp->text, 0, offset);
    Update ();
  }
}

typedef struct repopulateviewer
{
  Nlm_RepopulateViewer   repopulate_func; 
  Pointer                repopulate_data;
  Nlm_RepopulateDataFree free_data_func;
  TextViewFormPtr        tfp;
  FonT                   fnt;
} RepopulateViewerData, PNTR RepopulateViewerPtr;

static void CleanupRepopulateViewer (Nlm_GraphiC g, Nlm_VoidPtr data)
{
  RepopulateViewerPtr rp;
  
  rp = (RepopulateViewerPtr) data;
  if (rp != NULL && rp->free_data_func != NULL)
  {
    (rp->free_data_func)(rp->repopulate_data);
  }
  
  StdCleanupExtraProc (g, data);
}

static void RepopulateViewer (ButtoN b)
{
  RepopulateViewerPtr rp;
  CharPtr             new_path;
  
  rp = (RepopulateViewerPtr) GetObjectExtra (b);
  
  if (rp == NULL || rp->repopulate_func == NULL)
  {
    return;
  }
  
  new_path = (rp->repopulate_func) (rp->repopulate_data);
  if (new_path == NULL)
  {
    return;
  }
  
  if (rp->tfp->text != NULL)
  {
    FileToScrollText (rp->tfp->text, new_path);
  }
  else if (rp->tfp->doc != NULL)
  {
    txtColFmt.pixWidth = screenRect.right - screenRect.left;
    txtColFmt.pixInset = 8;
    DisplayFancy (rp->tfp->doc, new_path, &txtParFmt, &txtColFmt, rp->fnt, 0);
  }
  FileRemove (new_path);
  new_path = MemFree (new_path);
}

static void 
LaunchGeneralTextViewerEx 
(CharPtr path,
 CharPtr title, 
 Boolean useScrollText,
 Nlm_RepopulateViewer   repopulate_func, 
 Pointer                repopulate_data,
 Nlm_RepopulateDataFree free_data_func)

{
  ButtoN              b;
  FonT                fnt;
  GrouP               g;
  Int2                pixheight;
  Int2                pixwidth;
  StdEditorProcsPtr   sepp;
  TextViewFormPtr     tfp;
  TextViewProcsPtr    tvpp;
  RepopulateViewerPtr rp;
  WindoW              w;
#ifndef WIN_MAC
  MenU                m;
#endif

  tfp = (TextViewFormPtr) MemNew (sizeof (TextViewForm));
  if (tfp == NULL) return;

  w = DocumentWindow (-50, -33, -10, -10, title, StdCloseWindowProc, ResizeTextViewer);
  SetObjectExtra (w, tfp, StdCleanupFormProc);
  tfp->form = (ForM) w;
  tfp->exportform = ExportTextViewForm;
  tfp->formmessage = TextViewFormMessage;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    tfp->appmessage = sepp->handleMessages;
  }

  fnt = programFont;
  pixwidth = 35 * stdCharWidth + 17;
  pixheight = 20 * stdLineHeight;

  tvpp = (TextViewProcsPtr) GetAppProperty ("TextDisplayForm");
  if (tvpp != NULL) {
    pixwidth = MAX (pixwidth, tvpp->minPixelWidth);
    pixheight = MAX (pixheight, tvpp->minPixelHeight);
    if (tvpp->displayFont != NULL) {
      fnt = tvpp->displayFont;
    }
    if (tvpp->activateForm != NULL) {
      SetActivate (w, tvpp->activateForm);
    }
  }

#ifndef WIN_MAC
  m = PulldownMenu (w, "File");
  FormCommandItem (m, "Export...", (BaseFormPtr) tfp, VIB_MSG_EXPORT);
  SeparatorItem (m);
  FormCommandItem (m, "Close", (BaseFormPtr) tfp, VIB_MSG_CLOSE);
  if (tvpp != NULL && tvpp->useScrollText) {
    m = PulldownMenu (w, "Edit");
    FormCommandItem (m, CUT_MENU_ITEM, (BaseFormPtr) tfp, VIB_MSG_CUT);
    FormCommandItem (m, COPY_MENU_ITEM, (BaseFormPtr) tfp, VIB_MSG_COPY);
    FormCommandItem (m, PASTE_MENU_ITEM, (BaseFormPtr) tfp, VIB_MSG_PASTE);
    FormCommandItem (m, CLEAR_MENU_ITEM, (BaseFormPtr) tfp, VIB_MSG_DELETE);
  }
#endif

  /* right now Find button is only in indexer Sequin */
  if (useScrollText) {
    g = HiddenGroup (w, 5, 0, NULL);
    b = PushButton (g, "Find", FindInGeneralText);
    SetObjectExtra (b, tfp, NULL);
    tfp->find = DialogText (g, "", 10, NULL);
    if (repopulate_func != NULL)
    {
      rp = (RepopulateViewerPtr) MemNew (sizeof (RepopulateViewerData));
      if (rp != NULL)
      {
        rp->repopulate_func = repopulate_func;
        rp->repopulate_data = repopulate_data;
        rp->free_data_func = free_data_func;
        rp->tfp = tfp;
        rp->fnt = fnt;
        b = PushButton (g, "Repopulate", RepopulateViewer);
        SetObjectExtra (b, rp, CleanupRepopulateViewer);
      }
    }
  }

  if (useScrollText) {
    tfp->text = ScrollText (w, (pixwidth + stdCharWidth - 1) / stdCharWidth,
                            (pixheight + stdLineHeight - 1) / stdLineHeight,
                            fnt, FALSE, NULL);
    SetObjectExtra (tfp->text, tfp, NULL);
    RealizeWindow (w);
    if (! FileToScrollText (tfp->text, path)) {
      /* SetTitle (tfp->text, "(Text is too large to be displayed in this control.)"); */
      Remove (w);
      LaunchGeneralTextViewerEx (path, title, FALSE,
                                 repopulate_func, repopulate_data, free_data_func);
      return;
    }
  } else {
    tfp->doc = DocumentPanel (w, pixwidth, pixheight);
    SetObjectExtra (tfp->doc, tfp, NULL);
    RealizeWindow (w);
    txtColFmt.pixWidth = screenRect.right - screenRect.left;
    txtColFmt.pixInset = 8;
    DisplayFancy (tfp->doc, path, &txtParFmt, &txtColFmt, fnt, 0);
    /* document.c: SaveTableItem does not strip preceeding tabs if tabCount is 0 */
  }
  Show (w);
  Select (w);
  Update ();
}

extern void LaunchGeneralTextViewerWithRepopulate 
(CharPtr                path,
 CharPtr                title, 
 Nlm_RepopulateViewer   repopulate_func,
 Pointer                repopulate_data,
 Nlm_RepopulateDataFree free_data_func)
{
  TextViewProcsPtr tvpp;

  tvpp = (TextViewProcsPtr) GetAppProperty ("TextDisplayForm");
  if (tvpp != NULL && tvpp->useScrollText) {
    LaunchGeneralTextViewerEx (path, title, TRUE, 
                               repopulate_func, repopulate_data, free_data_func);
  } else {
    LaunchGeneralTextViewerEx (path, title, FALSE, 
                               repopulate_func, repopulate_data, free_data_func);
  }
}

extern void LaunchGeneralTextViewer (CharPtr path, CharPtr title)

{
  TextViewProcsPtr tvpp;

  tvpp = (TextViewProcsPtr) GetAppProperty ("TextDisplayForm");
  if (tvpp != NULL && tvpp->useScrollText) {
    LaunchGeneralTextViewerEx (path, title, TRUE, NULL, NULL, NULL);
  } else {
    LaunchGeneralTextViewerEx (path, title, FALSE, NULL, NULL, NULL);
  }
}

extern void LaunchAsnTextViewer (Pointer from, AsnWriteFunc writefunc, CharPtr title)

{
  AsnIoPtr  aip;
  Char      path [PATH_MAX];

  if (from == NULL || writefunc == NULL) return;
  if (StringHasNoText (title)) {
    title = "General ASN.1 Text Viewer";
  }

  TmpNam (path);
  aip = AsnIoOpen (path, "w");
  if (aip != NULL) {
    (*writefunc) (from, aip, NULL);
    AsnIoClose (aip);
    LaunchGeneralTextViewer (path, title);
  }
  FileRemove (path);
}

#ifndef WIN_MAC
extern void CreateStdEditorFormMenus (WindoW w)

{
  BaseFormPtr   bfp;
  MenU          m;

  bfp = (BaseFormPtr) GetObjectExtra (w);
  if (bfp != NULL) {
    m = PulldownMenu (w, "File");
    if (bfp->importform != NULL || bfp->exportform != NULL) {
      if (bfp->importform != NULL) {
        FormCommandItem (m, "Import...", bfp, VIB_MSG_IMPORT);
      }
      if (bfp->exportform != NULL) {
        FormCommandItem (m, "Export...", bfp, VIB_MSG_EXPORT);
      }
      SeparatorItem (m);
    }
    FormCommandItem (m, "Close", bfp, VIB_MSG_CLOSE);
    m = PulldownMenu (w, "Edit");
    FormCommandItem (m, CUT_MENU_ITEM, bfp, VIB_MSG_CUT);
    FormCommandItem (m, COPY_MENU_ITEM, bfp, VIB_MSG_COPY);
    FormCommandItem (m, PASTE_MENU_ITEM, bfp, VIB_MSG_PASTE);
    FormCommandItem (m, CLEAR_MENU_ITEM, bfp, VIB_MSG_DELETE);
  }
}
#endif

static Boolean DlgutilGetLowestStackSeqEntry (GatherContextPtr gcp)

{
  BaseFormPtr  bfp;
  Int2         i;

  if (gcp == NULL) return TRUE;
  bfp = (BaseFormPtr) gcp->userdata;
  if (bfp == NULL) return TRUE;
  if (gcp->gatherstack != NULL && gcp->numstack > 0) {
    for (i = 0; i < gcp->numstack; i++) {
      if (gcp->gatherstack [i].itemtype == OBJ_BIOSEQ ||
          gcp->gatherstack [i].itemtype == OBJ_BIOSEQSET) {
        bfp->input_itemID = gcp->gatherstack [i].itemID;
        bfp->input_itemtype = gcp->gatherstack [i].itemtype;
      }
    }
  }
  return FALSE;
}

extern Boolean SetClosestParentIfDuplicating (BaseFormPtr bfp)

{
  Uint2              itemID;
  Uint2              itemtype;
  StdEditorProcsPtr  sepp;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (bfp == NULL || sepp == NULL || (! sepp->duplicateExisting)) return FALSE;
  itemID = bfp->input_itemID;
  itemtype = bfp->input_itemtype;
  GatherItem (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype,
              (Pointer) bfp, DlgutilGetLowestStackSeqEntry);
  if (itemID == bfp->input_itemID && itemtype == bfp->input_itemtype) {
    return FALSE;
  }
  return TRUE;
}

/*****************************************************************************
*
*   Bond and Point SeqLoc dialogs
*
*****************************************************************************/

typedef struct pointpage {
  DIALOG_MESSAGE_BLOCK
  TexT               point;
  Int2               count;
  SeqEntryPtr        PNTR bsptr;
  EnumFieldAssoc     PNTR alist;
  PopuP              strand;
  PopuP              seqIdx;
  Boolean            nucsOK;
  Boolean            protsOK;
  Boolean            showIdTags;
} PointPage, PNTR PointPagePtr;

typedef struct bondpage {
  DIALOG_MESSAGE_BLOCK
  DialoG             pointA;
  DialoG             pointB;
} BondPage, PNTR BondPagePtr;

static void FillInPointProducts (SeqEntryPtr sep, Pointer mydata,
                                 Int4 index, Int2 indent)

{
  BioseqPtr     bsp;
  PointPagePtr  ppp;

  if (sep != NULL && mydata != NULL && sep->choice == 1) {
    ppp = (PointPagePtr) mydata;
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp != NULL) {
      if ((ppp->nucsOK && ISA_na (bsp->mol)) ||
          (ppp->protsOK && ISA_aa (bsp->mol))) {
        ppp->count++;
        ppp->bsptr [ppp->count] = sep;
      }
    }
  }
}

static ENUM_ALIST(strand_alist)
{" ",             Seq_strand_unknown},  /* 0 */
{"Plus",          Seq_strand_plus},     /* 1 */
{"Minus",         Seq_strand_minus},    /* 2 */
{"Both",          Seq_strand_both},     /* 3 */
{"Reverse",       Seq_strand_both_rev}, /* 4 */
{"Other",         Seq_strand_other},    /* 255 */
END_ENUM_ALIST

static void SeqPntPtrToPointPage (DialoG d, Pointer data)


{
  BioseqPtr     bsp;
  Int2          j;
  PointPagePtr  ppp;
  SeqEntryPtr   sep;
  Int2          seq;
  SeqPntPtr     spp;
  Char          str [16];
  Uint1         strand;

  ppp = (PointPagePtr) GetObjectExtra (d);
  if (ppp == NULL) return;
  spp = (SeqPntPtr) data;
  if (spp != NULL) {
    sprintf (str, "%ld", (long) (spp->point + 1));
    SafeSetTitle (ppp->point, str);
    seq = 0;
    strand = 0;
    bsp = BioseqFind (spp->id);
    if (bsp != NULL) {
      strand = spp->strand;
      if (strand > Seq_strand_both_rev && strand != Seq_strand_other) {
        strand = Seq_strand_unknown;
      }
      if (ppp->bsptr != NULL) {
        for (j = 1; j <= ppp->count && seq == 0; j++) {
          sep = ppp->bsptr [j];
          if (sep != NULL && sep->choice == 1) {
            if (bsp == (BioseqPtr) sep->data.ptrvalue) {
              seq = j;
            }
          }
        }
      }
    }
    SetEnumPopup (ppp->strand, strand_alist, (UIEnum) strand);
    SetEnumPopup (ppp->seqIdx, ppp->alist, (UIEnum) seq);
  } else {
    SafeSetTitle (ppp->point, "");
    SafeSetValue (ppp->strand, 0);
    SafeSetValue (ppp->seqIdx, 0);
  }
}

static Pointer PointPageToSeqPntPtr (DialoG d)

{
  BioseqPtr     bsp;
  PointPagePtr  ppp;
  SeqEntryPtr   sep;
  UIEnum        seq;
  SeqPntPtr     spp;
  Char          str [16];
  UIEnum        strand;
  Int4          val;

  ppp = (PointPagePtr) GetObjectExtra (d);
  if (ppp == NULL) return NULL;
  spp = NULL;
  GetTitle (ppp->point, str, sizeof (str) - 1);
  if (StrToLong (str, &val) && val > 0) {
    if (GetEnumPopup (ppp->seqIdx, ppp->alist, &seq) &&
        seq > 0 && seq <= ppp->count) {
      spp = SeqPntNew ();
      if (spp != NULL) {
        spp->point = val - 1;
        if (GetEnumPopup (ppp->strand, strand_alist, &strand)) {
          spp->strand = (Uint1) strand;
        }
        sep = ppp->bsptr [(Int2) seq];
        if (sep != NULL && sep->choice == 1) {
          bsp = (BioseqPtr) sep->data.ptrvalue;
          if (bsp != NULL) {
            spp->id = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (bsp->id, 0)));
          }
        }
      }
    }
  }
  return (Pointer) spp;
}

static void PointEditorMessage (DialoG d, Int2 mssg)

{
  PointPagePtr  ppp;

  ppp = (PointPagePtr) GetObjectExtra (d);
  if (ppp != NULL) {
    if (mssg == VIB_MSG_INIT) {
      SeqPntPtrToPointPage (d, NULL);
    } else if (mssg == VIB_MSG_ENTER) {
      Select (ppp->point);
    } else if (mssg == VIB_MSG_RESET) {
    }
  }
}

static void CleanupPointPage (GraphiC g, VoidPtr data)

{
  Int2          j;
  PointPagePtr  ppp;

  ppp = (PointPagePtr) data;
  if (ppp != NULL) {
    MemFree (ppp->bsptr);
    if (ppp->alist != NULL) {
      for (j = 0; j <= ppp->count + 1; j++) {
        MemFree (ppp->alist [j].name);
      }
    }
    MemFree (ppp->alist);
  }
  MemFree (data);
}

static DialoG CreatePointEditorDialog (GrouP h, CharPtr title, SeqEntryPtr sep,
                                       Boolean nucsOK, Boolean protsOK)

{
  BioseqPtr     bsp;
  Int4          count;
  GrouP         f;
  Int2          j;
  GrouP         m;
  GrouP         p;
  PointPagePtr  ppp;
  CharPtr       ptr;
  GrouP         s;
  Boolean       showIdTags;
  SeqIdPtr      sip;
  Char          str [128];

  p = HiddenGroup (h, 1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  ppp = (PointPagePtr) MemNew (sizeof (PointPage));
  if (ppp != NULL) {

    SetObjectExtra (p, ppp, CleanupPointPage);
    ppp->dialog = (DialoG) p;
    ppp->todialog = SeqPntPtrToPointPage;
    ppp->fromdialog = PointPageToSeqPntPtr;
    ppp->dialogmessage = PointEditorMessage;
    ppp->testdialog = NULL;
    ppp->importdialog = NULL;
    ppp->exportdialog = NULL;

    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    m = HiddenGroup (s, -1, 0, NULL);
    /*
    SetGroupSpacing (m, 10, 10);
    */

    ppp->nucsOK = nucsOK;
    ppp->protsOK = protsOK;
    ppp->showIdTags = FALSE;
    ppp->count = 0;

    if (sep != NULL) {
      count = SeqEntryCount (sep);
      count += 4;
      ppp->bsptr = MemNew (sizeof (BioseqPtr) * (size_t) count);
      ppp->alist = MemNew (sizeof (EnumFieldAssoc) * (size_t) count);
      ppp->count = 0;

      if (ppp->bsptr != NULL && ppp->alist != NULL) {
        SeqEntryExplore (sep, (Pointer) ppp, FillInPointProducts);
        j = 0;
        ppp->alist [j].name = StringSave ("     ");
        ppp->alist [j].value = (UIEnum) 0;
        for (j = 1; j <= ppp->count; j++) {
          sep = ppp->bsptr [j];
          if (sep != NULL && sep->choice == 1 && sep->data.ptrvalue != NULL) {
            bsp = (BioseqPtr) sep->data.ptrvalue;
            sip = SeqIdFindWorst (bsp->id);
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
            ppp->alist [j].name = StringSave (ptr);
            ppp->alist [j].value = (UIEnum) j;
          }
        }
        j = ppp->count + 1;
        ppp->alist [j].name = NULL;
        ppp->alist [j].value = (UIEnum) 0;
      }

    } else {
      ppp->alist = MemNew (sizeof (EnumFieldAssoc) * (size_t) 4);
      if (ppp->alist != NULL) {
        j = 0;
        ppp->alist [j].name = StringSave ("     ");
        ppp->alist [j].value = (UIEnum) 0;
        j = 1;
        ppp->alist [j].name = NULL;
        ppp->alist [j].value = (UIEnum) 0;

      }
    }

    f = HiddenGroup (m, 6, 0, NULL);
    /*StaticPrompt (f, "Point", 0, dialogTextHeight, programFont, 'l');*/
    ppp->point = DialogText (f, "", 5, NULL);
    if (nucsOK) {
      /*StaticPrompt (f, "Strand", 0, popupMenuHeight, programFont, 'c');*/
      ppp->strand = PopupList (f, TRUE, NULL);
      InitEnumPopup (ppp->strand, strand_alist, NULL);
    }
    /*StaticPrompt (f, "SeqID", 0, popupMenuHeight, programFont, 'l');*/
    ppp->seqIdx = PopupList (f, TRUE, NULL);
    InitEnumPopup (ppp->seqIdx, ppp->alist, NULL);
  }

  return (DialoG) p;
}

static void SeqLocPtrToBondPage (DialoG d, Pointer data)


{
  BondPagePtr  bpp;
  SeqBondPtr   sbp;
  SeqIdPtr     sip;
  SeqLocPtr    slp;
  SeqPnt       sqp;

  bpp = (BondPagePtr) GetObjectExtra (d);
  if (bpp == NULL) return;
  slp = (SeqLocPtr) data;
  if (slp != NULL) {
    if (slp->choice == SEQLOC_BOND) {
      sbp = (SeqBondPtr) slp->data.ptrvalue;
      if (sbp != NULL) {
        PointerToDialog (bpp->pointA, (Pointer) sbp->a);
        PointerToDialog (bpp->pointB, (Pointer) sbp->b);
      }
    } else {
      sip = SeqLocId (slp);
      if (sip != NULL) {
        sqp.strand = SeqLocStrand (slp);
        sqp.id = sip;
        sqp.point = SeqLocStart (slp);
        PointerToDialog (bpp->pointA, (Pointer) &sqp);
        sqp.point = SeqLocStop (slp);
        PointerToDialog (bpp->pointB, (Pointer) &sqp);
      }
    }
  }
}

static Pointer BondPageToSeqLocPtr (DialoG d)

{
  BondPagePtr  bpp;
  SeqBondPtr   sbp;
  SeqLocPtr    slp;

  bpp = (BondPagePtr) GetObjectExtra (d);
  if (bpp == NULL) return NULL;
  slp = NULL;
  sbp = SeqBondNew ();
  if (sbp != NULL) {
    slp = ValNodeNew (NULL);
    if (slp != NULL) {
      slp->choice = SEQLOC_BOND;
      slp->data.ptrvalue = (Pointer) sbp;
      sbp->a = DialogToPointer (bpp->pointA);
      sbp->b = DialogToPointer (bpp->pointB);
    } else {
      SeqBondFree (sbp);
    }
  }
  return (Pointer) slp;
}

static void BondEditorMessage (DialoG d, Int2 mssg)

{
  BondPagePtr  bpp;

  bpp = (BondPagePtr) GetObjectExtra (d);
  if (bpp != NULL) {
    if (mssg == VIB_MSG_INIT) {
      SeqLocPtrToBondPage (d, NULL);
    } else if (mssg == VIB_MSG_ENTER) {
      SendMessageToDialog (bpp->pointA, VIB_MSG_ENTER);
    } else if (mssg == VIB_MSG_RESET) {
    }
  }
}

extern DialoG CreateBondEditorDialog (GrouP h, CharPtr title, SeqEntryPtr sep);

extern DialoG CreateBondEditorDialog (GrouP h, CharPtr title, SeqEntryPtr sep)

{
  BondPagePtr  bpp;
  GrouP        f;
  GrouP        m;
  GrouP        p;
  GrouP        s;

  p = HiddenGroup (h, 1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  bpp = (BondPagePtr) MemNew (sizeof (BondPage));
  if (bpp != NULL) {

    SetObjectExtra (p, bpp, StdCleanupExtraProc);
    bpp->dialog = (DialoG) p;
    bpp->todialog = SeqLocPtrToBondPage;
    bpp->fromdialog = BondPageToSeqLocPtr;
    bpp->dialogmessage = BondEditorMessage;
    bpp->testdialog = NULL;
    bpp->importdialog = NULL;
    bpp->exportdialog = NULL;

    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    m = HiddenGroup (s, -1, 0, NULL);
    /*
    SetGroupSpacing (m, 10, 10);
    */

    f = HiddenGroup (m, 2, 0, NULL);
    StaticPrompt (f, "From", 0, popupMenuHeight, programFont, 'l');
    bpp->pointA = CreatePointEditorDialog (f, NULL, sep, FALSE, TRUE);
    StaticPrompt (f, "(To)", 0, popupMenuHeight, programFont, 'l');
    bpp->pointB = CreatePointEditorDialog (f, NULL, sep, FALSE, TRUE);
  }

  return (DialoG) p;
}

void GetRidOfEmptyFeatsDescStrings (Uint2 entityID, SeqEntryPtr sep)

{
  if (entityID < 1 && sep == NULL) return;
  if (entityID > 0 && sep == NULL) {
    sep = GetTopSeqEntryForEntityID (entityID);
  }
  if (sep == NULL) return;
  SeqEntryExplore (sep, NULL, GetRidOfEmptyFeatsDescCallback);
}

extern Int2 LIBCALLBACK StdVibrantEditorMsgFunc (OMMsgStructPtr ommsp)

{
  BaseFormPtr    bfp;
  OMUserDataPtr  omudp;

  omudp = (OMUserDataPtr)(ommsp->omuserdata);
  if (omudp == NULL) return OM_MSG_RET_ERROR;
  bfp = (BaseFormPtr) omudp->userdata.ptrvalue;
  if (bfp == NULL) return OM_MSG_RET_ERROR;
  switch (ommsp->message) {
    case OM_MSG_DEL:
      Remove (bfp->form);
      Update ();
      break;
    default :
      break;
  }
  return OM_MSG_RET_OK;
}

/* launch url section */

NLM_EXTERN void LaunchEntrezURL (CharPtr database, Int4 uid, CharPtr format)

{
#ifdef WIN_MOTIF
  NS_Window  window = NULL;
#endif

  Char  url [256];

  if (uid < 1 || StringHasNoText (database) || StringHasNoText (format)) return;
  sprintf (url,
           "http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=%s&list_uids=%ld&dopt=%s",
            database, (long) uid, format);

#ifdef WIN_MAC
  Nlm_SendURLAppleEvent (url, "MOSS", NULL);
#endif
#ifdef WIN_MSWIN
  if (! Nlm_MSWin_OpenDocument (url)) {
    Message (MSG_POST, "Unable to launch browser");
  }
#endif
#ifdef WIN_MOTIF
  if (! NS_OpenURL (&window, url, NULL, TRUE)) {
    Message (MSG_POST, "Unable to launch browser");
  }
  NS_WindowFree (window);
#endif
}

extern void ModalAcceptButton (ButtoN b)
{
  ModalAcceptCancelPtr acp;
  
  acp = (ModalAcceptCancelPtr) GetObjectExtra (b);
  if (acp != NULL)
  {
    acp->accepted = TRUE;
  }
}

extern void ModalCancelButton (ButtoN b)
{
  ModalAcceptCancelPtr acp;
  
  acp = (ModalAcceptCancelPtr) GetObjectExtra (b);
  if (acp != NULL)
  {
    acp->cancelled = TRUE;
  }
}

extern void ModalThirdOptionButton (ButtoN b)
{
  ModalAcceptCancelPtr acp;
  
  acp = (ModalAcceptCancelPtr) GetObjectExtra (b);
  if (acp != NULL)
  {
    acp->third_option = TRUE;
  }
}

typedef struct tabledisplay 
{
  DIALOG_MESSAGE_BLOCK
  PaneL panel;
  ValNodePtr row_list;
  Int4 frozen_header;
  Int4 frozen_left;
  Int4 table_inset;
  Int4 char_width;
  Int4 descent;
  FonT display_font;
  TableDisplayDblClick dbl_click;
  Pointer dbl_click_data;
  TableDisplayLeftInRed left_in_red;
  Pointer left_in_red_data;
} TableDisplayData, PNTR TableDisplayPtr;

extern ValNodePtr FreeTableDisplayRowList (ValNodePtr row_list)
{
  ValNodePtr row_vnp, column_list;
  
  if (row_list != NULL)
  {
    /* free table text */
    for (row_vnp = row_list; row_vnp != NULL; row_vnp = row_vnp->next)
    {
      column_list = (ValNodePtr) row_vnp->data.ptrvalue;
      row_vnp->data.ptrvalue = ValNodeFreeData (column_list);
    }
    row_list = ValNodeFree (row_list);
  }
  return row_list;
}

extern void PrintTableDisplayRowListToFile (ValNodePtr row_list, FILE *fp)
{
  ValNodePtr row_vnp, col_vnp, column_list;
  CharPtr    txt_val;
  
  if (row_list == NULL || fp == NULL)
  {
    return;
  }
  
  for (row_vnp = row_list; row_vnp != NULL; row_vnp = row_vnp->next)
  {
    column_list = (ValNodePtr) row_vnp->data.ptrvalue;
    for (col_vnp = column_list; col_vnp != NULL; col_vnp = col_vnp->next)
    {
      txt_val = (CharPtr) col_vnp->data.ptrvalue;
      if (!StringHasNoText (txt_val))
      {
        fprintf (fp, "%s", txt_val);
      }
      if (col_vnp->next == NULL)
      {
        fprintf (fp, "\n");
      }
      else
      {
        fprintf (fp, "\t");
      }
    }
  }
}

static ValNodePtr ValNodeStringListCopy (ValNodePtr orig_list)
{
  ValNodePtr new_list = NULL;
  
  if (orig_list == NULL)
  {
    return NULL;
  }
  
  new_list = ValNodeNew (NULL);
  new_list->choice = orig_list->choice;
  new_list->data.ptrvalue = StringSave (orig_list->data.ptrvalue);
  new_list->next = ValNodeStringListCopy (orig_list->next);
  return new_list;
}

extern ValNodePtr CopyTableDisplayRowList (ValNodePtr row_list)
{
  ValNodePtr new_row_list = NULL;
  
  if (row_list == NULL)
  {
    return NULL; 
  }
  
  new_row_list = ValNodeNew (NULL);
  new_row_list->choice = row_list->choice;
  new_row_list->data.ptrvalue = ValNodeStringListCopy (row_list->data.ptrvalue);
  new_row_list->next = CopyTableDisplayRowList (row_list->next);
  return new_row_list;
}

static void CleanupTableDisplayDialog (GraphiC g, VoidPtr data)
{
  TableDisplayPtr dlg;

  dlg = (TableDisplayPtr) data;
  if (dlg != NULL) {
    dlg->row_list = FreeTableDisplayRowList (dlg->row_list);
  }
  StdCleanupExtraProc (g, data);
}

static void UpdateTableDisplayDialogScrollBars (TableDisplayPtr dlg)
{
  BaR  sb_vert;
  BaR  sb_horiz;
  Int4 start_row, start_col;
  Int4 num_rows, num_columns, visible_rows;
  Int4 new_vmax, new_hmax, old_vmax, old_hmax;
  RecT r;
  Int4 x, y;
  
  if (dlg == NULL)
  {
    return;
  }
  
  sb_vert  = GetSlateVScrollBar ((SlatE) dlg->panel);
  sb_horiz = GetSlateHScrollBar ((SlatE) dlg->panel);
  
  start_row = GetBarValue (sb_vert) + dlg->frozen_header;
  start_col = GetBarValue (sb_horiz) + dlg->frozen_left;
    
  if (dlg->row_list == NULL)
  {
    num_rows = 0;
    num_columns = 0;
  }
  else
  {
    num_rows = ValNodeLen (dlg->row_list);
    num_columns = ValNodeLen (dlg->row_list->data.ptrvalue);
  }

  ObjectRect (dlg->panel, &r);
  InsetRect (&r, dlg->table_inset, dlg->table_inset);
  x = r.left + 1;
  y = r.top + stdLineHeight;
    
  visible_rows = (r.bottom - r.top - 2 * dlg->table_inset) / stdLineHeight - dlg->frozen_header;
  new_vmax = num_rows - visible_rows - 1;
  new_hmax = num_columns - dlg->frozen_left - 1;
  if (new_vmax < 0)
  {
    new_vmax = 0;
  }
  if (new_hmax < 0)
  {
    new_hmax = 0;
  }
  old_vmax = GetBarMax (sb_vert);
  old_hmax = GetBarMax (sb_horiz);
  
  if (old_vmax != new_vmax)
  {
    CorrectBarMax (sb_vert, new_vmax);
    if (start_row > new_vmax + dlg->frozen_header)
    {
      start_row = new_vmax + dlg->frozen_header;
    }
    CorrectBarValue (sb_vert, start_row - dlg->frozen_header);
    CorrectBarPage (sb_vert, 1, 1);
  }
  
  if (old_hmax != new_hmax)
  {
    CorrectBarMax (sb_horiz, new_hmax);
    if (start_col > new_hmax + dlg->frozen_left)
    {
      start_col = new_hmax + dlg->frozen_left;
    }
    CorrectBarValue (sb_horiz, start_col - dlg->frozen_left);
    CorrectBarPage (sb_horiz, 1, 1);
  }  
}

static void RowsToTableDisplayDialog (DialoG d, Pointer userdata)
{
  TableDisplayPtr dlg;
  RecT            r;
	WindoW          temport;

  dlg = (TableDisplayPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  dlg->row_list = FreeTableDisplayRowList (dlg->row_list);
  dlg->row_list = CopyTableDisplayRowList (userdata);
  UpdateTableDisplayDialogScrollBars (dlg);
  temport = SavePort (dlg->panel);
  Select (dlg->panel);
  ObjectRect (dlg->panel, &r);
  InvalRect (&r);  
  Update ();
  RestorePort (temport);
}

static Pointer TableDisplayDialogToRows (DialoG d)
{
  TableDisplayPtr dlg;

  dlg = (TableDisplayPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  
  return CopyTableDisplayRowList (dlg->row_list);
}

static Int4 
PrepareTableDisplayTextBuffer 
(CharPtr buf,
 Int4    remaining_chars_in_row,
 Int4    col_width,
 CharPtr data)
{
  Int4 chars_to_paint = 0;
  if (buf == NULL)
  {
    return 0;
  }
  
  if (remaining_chars_in_row < col_width)
  {
    chars_to_paint = remaining_chars_in_row;
    StringNCpy (buf, data, chars_to_paint);
    buf [chars_to_paint] = 0;
  }
  else
  {
    chars_to_paint = col_width;
    StringNCpy (buf, data, chars_to_paint);
    buf [chars_to_paint] = 0;
    if (StringLen (data) > chars_to_paint && chars_to_paint > 2)
    {
      buf [chars_to_paint - 1] = '.';
      buf [chars_to_paint - 2] = '.';
      buf [chars_to_paint - 3] = '.';
    }
  }
  return chars_to_paint;
}

static void DrawTableDisplayLine (Int4 x, Int4 y, 
 ValNodePtr header_row,
 ValNodePtr data_row,
 CharPtr    buf,
 Int4       row_length,
 Int4       frozen_left,
 Int4       start_col,
 Int4       char_width,
 Int4       descent,
 Boolean    left_in_red)
{
  ValNodePtr header_vnp, data_vnp;
  Int4       x_offset, chars_to_paint, col_num;
  PoinT      pt1, pt2;
  RecT       rct;
  
  /* draw left margin */
  
  for (header_vnp = header_row, data_vnp = data_row, x_offset = 0, col_num = 0;
       header_vnp != NULL && data_vnp != NULL && x_offset < row_length && col_num < frozen_left;
       header_vnp = header_vnp->next, data_vnp = data_vnp->next, col_num++)
  {
    Gray ();
    InvertColors ();
    if (left_in_red)
    {
      Red ();
    }
    else
    {
      White ();
    }
    chars_to_paint = PrepareTableDisplayTextBuffer (buf, 
                                   (row_length - x_offset) / char_width,
                                   header_vnp->choice,
                                   data_vnp->data.ptrvalue);
        
    LoadRect (&rct, x + x_offset, y + descent,
              x + x_offset + (chars_to_paint + 2) * char_width, 
              y - stdLineHeight + descent);
    EraseRect (&rct);

    PaintStringEx ( (CharPtr)buf, x + x_offset, y);
    x_offset += (chars_to_paint + 2) * char_width;
    InvertColors ();
    Black ();
  }
  
  if (frozen_left > 0)
  {
    pt1.x = x + x_offset - 1;
    pt1.y = y;
    pt2.x = x + x_offset - 1;
    pt2.y = y - stdLineHeight;
    DrawLine (pt1, pt2);
  }

  
  while (col_num < start_col && header_vnp != NULL && data_vnp != NULL)
  {
    col_num++;
    header_vnp = header_vnp->next;
    data_vnp = data_vnp->next;
  }
  
  /* draw unfrozen columns */
  while (header_vnp != NULL && data_vnp != NULL && x_offset < row_length)
  {
    chars_to_paint = MIN (header_vnp->choice, (row_length - x_offset)/char_width);
    StringNCpy (buf, data_vnp->data.ptrvalue, chars_to_paint);
    buf [chars_to_paint] = 0;
    chars_to_paint = PrepareTableDisplayTextBuffer (buf, 
                                   (row_length - x_offset) / char_width,
                                   header_vnp->choice,
                                   data_vnp->data.ptrvalue);

    PaintStringEx ( (CharPtr)buf, x + x_offset, y);
    x_offset += (chars_to_paint + 2) * char_width;
    header_vnp = header_vnp->next;
    data_vnp = data_vnp->next;
  }
}

static void OnDrawTableDisplay (PaneL p)
{
  TableDisplayPtr dlg;
  BaR             sb_vert, sb_horiz;
  Int4            start_row, start_col;
  RecT            r;
  Int4            x, y, row, row_length;
  CharPtr         row_buffer;
  ValNodePtr      row_vnp;
  PoinT           pt1, pt2;
  Boolean         left_in_red;

  dlg = (TableDisplayPtr) GetObjectExtra (p);
  if (dlg == NULL)
  {
    return;
  }
  
  sb_vert  = GetSlateVScrollBar ((SlatE) p);
  sb_horiz = GetSlateHScrollBar ((SlatE) p);
  
  start_row = GetBarValue (sb_vert) + dlg->frozen_header;
  start_col = GetBarValue (sb_horiz) + dlg->frozen_left;

  ObjectRect (p, &r);
  InsetRect (&r, dlg->table_inset, dlg->table_inset);
  x = r.left + 1;
  y = r.top + stdLineHeight;

  SelectFont (programFont); 
  
  row_length = r.right - r.left - 2;
  row_buffer = (CharPtr) MemNew (((row_length / dlg->char_width) + 1) * sizeof (Char));
  
  for (row = 0, row_vnp = dlg->row_list;
       row < dlg->frozen_header && y <= r.bottom - 2 * dlg->table_inset && row_vnp != NULL;
       row++, row_vnp = row_vnp->next)
  {
    DrawTableDisplayLine (x, y, dlg->row_list->data.ptrvalue, row_vnp->data.ptrvalue,
                          row_buffer, row_length, dlg->frozen_left, start_col, 
                          dlg->char_width, dlg->descent, FALSE);
    y += stdLineHeight;
  }
  
  while (row < start_row && row_vnp != NULL)
  {
    row++;
    row_vnp = row_vnp->next;
  }
  
  while (row_vnp != NULL && y <= r.bottom - 2 * dlg->table_inset)
  {
    left_in_red = FALSE;
    if (dlg->left_in_red != NULL)
    {
      left_in_red = (dlg->left_in_red) (row, dlg->row_list, dlg->left_in_red_data);
    }
    DrawTableDisplayLine (x, y, dlg->row_list->data.ptrvalue, row_vnp->data.ptrvalue,
                          row_buffer, row_length, dlg->frozen_left, start_col,
                          dlg->char_width, dlg->descent, left_in_red);
    row_vnp = row_vnp->next;
    y += stdLineHeight;
    row++;
  }
  
  /* draw line to separate header from remaining lines */
  if (dlg->frozen_header > 0)
  {
    Black ();
    pt1.x = x;
    pt1.y = r.top + stdLineHeight + dlg->descent;
    pt2.x = x + row_length;
    pt2.y = r.top + stdLineHeight + dlg->descent;
    DrawLine (pt1, pt2);
  }
  

}

static void OnVScrollTableDisplay (BaR sb, SlatE s, Int4 newval, Int4 oldval)
{
  RecT   r;
  WindoW temport;

  temport = SavePort (s);
  Select (s);
  ObjectRect (s, &r);
  InvalRect (&r);  
  RestorePort (temport);
  Update ();
}

static void OnHScrollTableDisplay (BaR sb, SlatE s, Int4 newval, Int4 oldval)
{
  RecT   r;
  WindoW temport;

  temport = SavePort (s);
  Select (s);
  ObjectRect (s, &r);
  InvalRect (&r);  
  RestorePort (temport);
  Update ();
}

static PoinT GetTableDisplayCell (TableDisplayPtr dlg, PoinT pt)
{
  BaR sb_horiz;
  BaR sb_vert;
  Int4 start_row, start_col;
  RecT r;
  PoinT cell_coord;
  Int4  x, y;
  ValNodePtr header_vnp;
  Int4  col_width;
  
  cell_coord.x = 0;
  cell_coord.y = 0;
  
  if (dlg == NULL || dlg->row_list == NULL)
  {
    return cell_coord;
  }
  
  sb_vert  = GetSlateVScrollBar ((SlatE) dlg->panel);
  sb_horiz = GetSlateHScrollBar ((SlatE) dlg->panel);
  
  start_row = GetBarValue (sb_vert) + dlg->frozen_header;
  start_col = GetBarValue (sb_horiz) + dlg->frozen_left;
  
  ObjectRect (dlg->panel, &r);
  InsetRect (&r, dlg->table_inset, dlg->table_inset);
  x = pt.x - r.left;
  y = pt.y - r.top;
  
  cell_coord.y = y / stdLineHeight;
  
  if (cell_coord.y >= dlg->frozen_header)
  {
    cell_coord.y += GetBarValue (sb_vert);
  }

  header_vnp = dlg->row_list->data.ptrvalue;
  
  col_width = 0;
  while (header_vnp != NULL && col_width + (header_vnp->choice + 2) * dlg->char_width < x
         && cell_coord.x < dlg->frozen_left)
  {
    cell_coord.x++;
    col_width += (header_vnp->choice + 2) * dlg->char_width;
    header_vnp = header_vnp->next;
  }
  
  if (cell_coord.x >= dlg->frozen_left)
  {
    /* skip over unfrozen columns not currently displayed */
    while (header_vnp != NULL && cell_coord.x < start_col)
    {
      header_vnp = header_vnp->next;
      cell_coord.x++;
    }
  
    while (header_vnp != NULL && col_width + (header_vnp->choice + 2) * dlg->char_width < x)
    {
      cell_coord.x++;
      col_width += (header_vnp->choice + 2) * dlg->char_width;
      header_vnp = header_vnp->next;
    }
  }
  return cell_coord;
}

extern CharPtr GetRowListCellText (ValNodePtr row_list, Int4 row, Int4 column)
{
  ValNodePtr row_vnp, col_vnp;
  Int4       row_num, col_num;
  
  if (row_list == NULL || row < 0 || column < 0)
  {
    return NULL;
  }
  
  for (row_vnp = row_list, row_num = 0;
       row_vnp != NULL && row_num < row;
       row_vnp = row_vnp->next, row_num++)
  {
  }
  if (row_num != row || row_vnp == NULL)
  {
    return NULL;
  }
  for (col_vnp = row_vnp->data.ptrvalue, col_num = 0;
       col_vnp != NULL && col_num < column;
       col_vnp = col_vnp->next, col_num++)
  {
  }
  if (col_num != column || col_vnp == NULL)
  {
    return NULL;
  }
  else
  {
    return StringSave (col_vnp->data.ptrvalue);
  }  
}

static CharPtr TableDisplayGetTextForCell (TableDisplayPtr dlg, PoinT pt)
{
  if (dlg == NULL || dlg->row_list == NULL || pt.x < 0 || pt.y < 0)
  {
    return NULL;
  }
  
  return GetRowListCellText (dlg->row_list, pt.y, pt.x);
}

static void TableDisplayOnClick (PaneL p, PoinT pt)
{
  TableDisplayPtr dlg;
  Boolean         dbl_click;
  PoinT           cell_coord;
  PoinT           header_coord;
  CharPtr         cell_text;
  CharPtr         header_text;
  
  dlg = (TableDisplayPtr) GetObjectExtra (p);
  if (dlg == NULL)
  {
    return;
  }
  
  dbl_click = dblClick;
  if (dbl_click && dlg->dbl_click != NULL)
  {
    cell_coord = GetTableDisplayCell (dlg, pt);
    cell_text = TableDisplayGetTextForCell (dlg, cell_coord);
    header_coord.x = cell_coord.x;
    header_coord.y = 0;
    header_text = TableDisplayGetTextForCell (dlg, header_coord);
    (dlg->dbl_click) (cell_coord, header_text, cell_text, dlg->dbl_click_data);
    MemFree (cell_text);
    MemFree (header_text);
  }
}

extern FonT GetTableDisplayDefaultFont (void)
{
  FonT display_font = NULL;
  
#ifdef WIN_MAC
  display_font = ParseFont ("Monaco, 9");
#endif
#ifdef WIN_MSWIN
  display_font = ParseFont ("Courier, 9");
#endif
#ifdef WIN_MOTIF
  display_font = ParseFont ("fixed, 12");
#endif  
  return display_font;
}

extern DialoG TableDisplayDialog (GrouP parent, Int4 width, Int4 height,
                                  Int4 frozen_header, Int4 frozen_left,
                                  TableDisplayDblClick dbl_click,
                                  Pointer dbl_click_data,
                                  TableDisplayLeftInRed left_in_red,
                                  Pointer left_in_red_data)
{
  TableDisplayPtr dlg;
  GrouP           p;
  
  dlg = (TableDisplayPtr) MemNew (sizeof (TableDisplayData));
  if (dlg == NULL)
  {
    return NULL;
  }
  p = HiddenGroup (parent, -1, 0, NULL);
  SetObjectExtra (p, dlg, CleanupTableDisplayDialog);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = RowsToTableDisplayDialog;
  dlg->fromdialog = TableDisplayDialogToRows;
  dlg->dialogmessage = NULL;
  dlg->testdialog = NULL;
  
  dlg->row_list = NULL;
  dlg->frozen_header = frozen_header;
  dlg->frozen_left = frozen_left;
  dlg->table_inset = 4;
  dlg->dbl_click = dbl_click;
  dlg->dbl_click_data = dbl_click_data;
  dlg->left_in_red = left_in_red;
  dlg->left_in_red_data = left_in_red_data;
  
  dlg->display_font = GetTableDisplayDefaultFont ();

  SelectFont (dlg->display_font);
  dlg->char_width  = CharWidth ('0');
  dlg->descent = Descent ();
  
  dlg->panel = AutonomousPanel4 (p, width, height, OnDrawTableDisplay,
                               OnVScrollTableDisplay, OnHScrollTableDisplay,
                               sizeof (TableDisplayData), NULL, NULL); 
  SetObjectExtra (dlg->panel, dlg, NULL);
  SetPanelClick(dlg->panel, TableDisplayOnClick, NULL, NULL, NULL);
  
  return (DialoG) p;  
}

typedef struct multiselectdialog
{
  DIALOG_MESSAGE_BLOCK
  DoC                  doc;
  ValNodePtr           selected_list;
  ParData              listPar;
  ColData              listCol;
  Int4                 num_choices;
  Nlm_ChangeNotifyProc change_notify;
  Pointer              change_userdata;    
} MultiSelectDialogData, PNTR MultiSelectDialogPtr;

static void DataToMultiSelectionDialog (DialoG d, Pointer userdata)
{
  MultiSelectDialogPtr dlg;
  ValNodePtr           vnp;
  Boolean              all_selected = FALSE;
  
  dlg = (MultiSelectDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  
  dlg->selected_list = ValNodeFree (dlg->selected_list);
  for (vnp = (ValNodePtr) userdata; vnp != NULL && ! all_selected; vnp = vnp->next)
  {
    if (vnp->data.intvalue == 1)
    {
      all_selected = TRUE;
    }
    else
    {
      ValNodeAddInt (&(dlg->selected_list), vnp->data.intvalue, vnp->data.intvalue);
    }
  }
  if (all_selected)
  {
    dlg->selected_list = ValNodeFree (dlg->selected_list);
    ValNodeAddInt (&(dlg->selected_list), 1, 1);
  }
  InvalDocRows (dlg->doc, 0, 0, 0);
}

static Pointer MultiSelectionDialogToData (DialoG d)
{
  MultiSelectDialogPtr dlg;
  ValNodePtr           output_list = NULL, vnp;
 
  dlg = (MultiSelectDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;
  
  for (vnp = dlg->selected_list; vnp != NULL; vnp = vnp->next)
  {
    ValNodeAddInt (&output_list, vnp->choice, vnp->choice);
  }
  return output_list;
}

static void 
AddChoiceToSelection 
(Int2                 choice_num,
 Boolean              remove_other_choices,
 MultiSelectDialogPtr dlg)
{
  ValNodePtr already_sel, vnp;
  
  if (dlg == NULL || dlg->doc == NULL || choice_num < 1)
  {
    return;
  }
  
  if (choice_num == 1)
  {
    remove_other_choices = TRUE;
  }
  else 
  {
    /* if we have added a choice other than "All", remove the "All"
     * selection if it was present.
     */
    ValNodeFree (ValNodeExtractList (&(dlg->selected_list), 1));
    InvalDocRows (dlg->doc, 1, 1, 1);
  }
  
  already_sel = ValNodeExtractList (&(dlg->selected_list), choice_num);
  if (already_sel == NULL)
  {
    if (remove_other_choices)
    {
      /* delete old selections */
      for (vnp = dlg->selected_list; vnp != NULL; vnp = vnp->next)
      {
        InvalDocRows (dlg->doc, vnp->choice, 1, 1);
      }
      dlg->selected_list = ValNodeFree (dlg->selected_list);      
      ValNodeAddInt (&dlg->selected_list, choice_num, choice_num);
      InvalDocRows (dlg->doc, choice_num, 1, 1);
    }
    else
    {
      /* add new selection */
      ValNodeAddInt (&(dlg->selected_list), choice_num, choice_num);
      InvalDocRows (dlg->doc, choice_num, 1, 1);
    }
  }
  else
  {
    already_sel = ValNodeFree (already_sel);
    InvalDocRows (dlg->doc, choice_num, 1, 1); 
  }
  if (dlg->change_notify != NULL)
  {
    (dlg->change_notify) (dlg->change_userdata);
  }
}

static void SelectChoice (DoC d, PoinT pt)
{
  Int2      item, row;
  Boolean   remove_other_choices;
  MultiSelectDialogPtr dlg;
  
  remove_other_choices = ! ctrlKey;
  
  dlg = (MultiSelectDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  
  MapDocPoint (d, pt, &item, &row, NULL, NULL);
  AddChoiceToSelection (item, remove_other_choices, dlg);
}

static Boolean ChoiceHighlight (DoC doc, Int2 item, Int2 row, Int2 col)
{
  MultiSelectDialogPtr dlg;
  ValNodePtr           vnp;
  
  dlg = (MultiSelectDialogPtr) GetObjectExtra (doc);
  if (dlg == NULL) return FALSE;
  
  for (vnp = dlg->selected_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == item)
    {
      return TRUE;
    }
  }
  return FALSE;
}

static void ChoiceOnKey (SlatE s, Char ch)
{
  MultiSelectDialogPtr dlg;
  CharPtr              str;
  Int2                 start_pos;
  Boolean              found = FALSE;
  ValNodePtr           vnp;
  
  dlg = (MultiSelectDialogPtr) GetObjectExtra (s);
  if (dlg == NULL) return;

  if ( (int) ch == 0 ) return;
  
  if (isalpha (ch))
  {
    /* find the position of the last choice we added */
    for (vnp = dlg->selected_list; vnp != NULL && vnp->next != NULL; vnp = vnp->next)
    {}
    
    if (vnp == NULL)
    {
      GetOffset (dlg->doc, NULL, &start_pos);
      /* start pos is one less than document row */
      start_pos ++;
      /* want to start at row after top row */
      start_pos ++;
    }
    else
    {
      /* want to start at row after currently selected row */
      start_pos = vnp->choice + 1;
    }
    
    while (!found && start_pos <= dlg->num_choices)
    {
      str = GetDocText (dlg->doc, start_pos, 1, 1);
      if (tolower (str [0]) == tolower (ch))
      {
        SetOffset (dlg->doc, 0, start_pos - 1);
        AddChoiceToSelection (start_pos, TRUE, dlg);
        found = TRUE;
      }
      str = MemFree (str);
      start_pos ++;
    }
    if (!found)
    {
      /* start searching at the top of the list */
      start_pos = 1;
      while (!found && start_pos <= dlg->num_choices)
      {
        str = GetDocText (dlg->doc, start_pos, 1, 1);
        if (tolower (str [0]) == tolower (ch))
        {
          SetOffset (dlg->doc, 0, start_pos - 1);
          AddChoiceToSelection (start_pos, TRUE, dlg);
          found = TRUE;
        }
        str = MemFree (str);
        start_pos ++;
      }
    }
  }
  else if (ch == NLM_DOWN)
  {
    /* down key */
    if (dlg->selected_list == NULL)
    {
      GetOffset (dlg->doc, NULL, &start_pos);
      start_pos ++;
    }
    else
    {
      start_pos = dlg->selected_list->choice;
    }
    start_pos ++;
    
    if (start_pos <= dlg->num_choices)
    {
      SetOffset (dlg->doc, 0, start_pos - 1);
      AddChoiceToSelection (start_pos, TRUE, dlg); 
      InvalDocRows (dlg->doc, 0, 0, 0);
    }
  }
  else if (ch == NLM_UP)
  {
    /* up key */
    if (dlg->selected_list == NULL)
    {
      GetOffset (dlg->doc, NULL, &start_pos);
      start_pos ++;
    }
    else
    {
      start_pos = dlg->selected_list->choice;
    }
    start_pos --;
    if (start_pos > 0)
    {
      SetOffset (dlg->doc, 0, start_pos - 1);
      AddChoiceToSelection (start_pos, TRUE, dlg); 
      InvalDocRows (dlg->doc, 0, 0, 0);
    }
  }
}

static DialoG 
MultiSelectDialog 
(GrouP      parent,
 ValNodePtr choice_list,
 Int4       list_height,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  MultiSelectDialogPtr dlg;
  GrouP                p;
  Int4                 height;
  Int4                 width;
  RecT                 r;
  ValNodePtr           vnp;
  
  if (choice_list == NULL)
  {
    return NULL;
  }
  dlg = (MultiSelectDialogPtr) MemNew (sizeof (MultiSelectDialogData));
  if (dlg == NULL)
  {
    return NULL;
  }
  
  p = HiddenGroup (parent, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = DataToMultiSelectionDialog;
  dlg->fromdialog = MultiSelectionDialogToData;
  dlg->dialogmessage = NULL;
  dlg->testdialog = NULL;  
  
  dlg->selected_list = NULL;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  
  SelectFont (systemFont);
  width = 0;
  for (vnp = choice_list; vnp != NULL; vnp = vnp->next)
  {
    width = MAX (width, StringWidth (vnp->data.ptrvalue));
  }
  
  dlg->num_choices = ValNodeLen (choice_list) + 1;
  
  height = LineHeight ();
  dlg->doc = DocumentPanel (p, width, height * list_height);
  SetObjectExtra (dlg->doc, dlg, NULL);
  
  dlg->listPar.openSpace = FALSE;
  dlg->listPar.keepWithNext = FALSE;
  dlg->listPar.keepTogether = FALSE;
  dlg->listPar.newPage = FALSE;
  dlg->listPar.tabStops = FALSE;
  dlg->listPar.minLines = 0;
  dlg->listPar.minHeight = 0;

  ObjectRect (dlg->doc, &r);
  InsetRect (&r, 4, 4);
  dlg->listCol.pixWidth = r.right - r.left;
  dlg->listCol.pixInset = 0;
  dlg->listCol.charWidth = 160;
  dlg->listCol.charInset = 0;
  dlg->listCol.font = systemFont;
  dlg->listCol.just = 'l';
  dlg->listCol.wrap = FALSE;
  dlg->listCol.bar = FALSE;
  dlg->listCol.underline = FALSE;
  dlg->listCol.left = FALSE;
  dlg->listCol.last = TRUE;  
  
	AppendText (dlg->doc, "All", &(dlg->listPar), &(dlg->listCol), programFont);
  for (vnp = choice_list; vnp != NULL; vnp = vnp->next)
  {
	  AppendText (dlg->doc, vnp->data.ptrvalue, &(dlg->listPar), &(dlg->listCol), programFont);
  }
  SetDocAutoAdjust (dlg->doc, FALSE);
  SetDocProcs (dlg->doc, SelectChoice, NULL, NULL, NULL);
  SetDocShade (dlg->doc, NULL, NULL, ChoiceHighlight, NULL);
  SetSlateChar ((SlatE) dlg->doc, ChoiceOnKey);
  InvalDocument (dlg->doc);

  return (DialoG) p;
}

typedef struct selectiondialog
{
  DIALOG_MESSAGE_BLOCK
  DialoG     multi_select_dlg;
  LisT       list_ctrl;
  PopuP      popup_ctrl;
  Int4       num_choices;
  CharPtr    err_msg;
  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} SelectionDialogData, PNTR SelectionDialogPtr;

static void ResetSelectionDialog (SelectionDialogPtr dlg)
{  
  if (dlg != NULL)
  {
    if (dlg->multi_select_dlg != NULL)
    {
      PointerToDialog (dlg->multi_select_dlg, NULL);
    }
    else if (dlg->list_ctrl != NULL)
    {
      SetValue (dlg->list_ctrl, 0);
    }
    else if (dlg->popup_ctrl != NULL)
    {
      SetValue (dlg->popup_ctrl, 0);
    }
  } 
}

static void SelectionDialogChanged (LisT l)
{
  SelectionDialogPtr dlg;

  dlg = (SelectionDialogPtr) GetObjectExtra (l);
  if (dlg == NULL)
  {
    return;
  }
    
  if (dlg->change_notify != NULL)
  {
    (dlg->change_notify)(dlg->change_userdata);
  } 
}

static void SelectionDialogPopupChanged (PopuP p)
{
  SelectionDialogPtr dlg;
 
  dlg = (SelectionDialogPtr) GetObjectExtra (p);
  if (dlg == NULL)
  {
    return;
  }
    
  if (dlg->change_notify != NULL)
  {
    (dlg->change_notify)(dlg->change_userdata);
  } 
}

static void SelectionListToSelectionDialog (DialoG d, Pointer userdata)
{
  SelectionDialogPtr dlg;
  ValNodePtr         selected_list;

  dlg = (SelectionDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  ResetSelectionDialog (dlg);  
  selected_list = (ValNodePtr) userdata;
  if (dlg->multi_select_dlg != NULL)
  {
    PointerToDialog (dlg->multi_select_dlg, selected_list);
  }
  else if (dlg->list_ctrl != NULL)
  {
    if (selected_list == NULL)
    {
      SetValue (dlg->list_ctrl, 0);
    }
    else
    {
      SetValue (dlg->list_ctrl, selected_list->data.intvalue);
    }
  }
  else if (dlg->popup_ctrl)
  {
    if (selected_list == NULL)
    {
      SetValue (dlg->popup_ctrl, 0);
    }
    else
    {
      SetValue (dlg->popup_ctrl, selected_list->data.intvalue);
    }
  }
}


static Pointer SelectionDialogToSelectionList (DialoG d)
{
  SelectionDialogPtr dlg;
  ValNodePtr         sel_list = NULL, vnp;
  Int4               i = 0;

  dlg = (SelectionDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  
  if (dlg->multi_select_dlg != NULL)
  {
    sel_list = (ValNodePtr) DialogToPointer (dlg->multi_select_dlg);
    if (sel_list != NULL && sel_list->choice == 1)
    {
      sel_list = ValNodeFree (sel_list);
      for (i = 2; i <= dlg->num_choices; i++)
      {
        ValNodeAddInt (&sel_list, 0, i - 1);
      }
    }
    else
    {
      for (vnp = sel_list; vnp != NULL; vnp = vnp->next)
      {
        vnp->choice = 0;
        vnp->data.intvalue = vnp->data.intvalue - 1;
      }
    }
  }
  else
  {
    if (dlg->list_ctrl != NULL)
    {
      i = GetValue (dlg->list_ctrl);
    }
    else if (dlg->popup_ctrl != NULL)
    {
      i = GetValue (dlg->popup_ctrl);
    }
    if (i > 0)
    {
      ValNodeAddInt (&sel_list, 0, i);
    }
  }
  return (Pointer) sel_list;
}

static void SelectionDialogMessage (DialoG d, Int2 mssg)

{
  SelectionDialogPtr dlg;
  ValNode            vn;

  dlg = (SelectionDialogPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        /* reset list */
        ResetSelectionDialog (dlg);
        break;
      case VIB_MSG_ENTER :
        if (dlg->multi_select_dlg != NULL)
        {
          Select (dlg->multi_select_dlg);
        }
        else if (dlg->list_ctrl != NULL)
        {
          Select (dlg->list_ctrl);
        }
        else if (dlg->popup_ctrl != NULL)
        {
          Select (dlg->popup_ctrl);
        }
        break;
      case NUM_VIB_MSG + 1:
        if (dlg->multi_select_dlg != NULL)
        {
          vn.next = NULL;
          vn.choice = 1;
          vn.data.intvalue = 1;
          PointerToDialog (dlg->multi_select_dlg, &vn);
        }
        else if (dlg->list_ctrl != NULL)
        {
          SetItemStatus (dlg->list_ctrl, 1, TRUE);
        }
        else if (dlg->popup_ctrl != NULL)
        {
          SetValue (dlg->popup_ctrl, 1);
        }
        SelectionDialogChanged (dlg->list_ctrl);
        break;
      default :
        break;
    }
  }
}

static ValNodePtr TestSelectionDialog (DialoG d)

{
  SelectionDialogPtr dlg;
  ValNodePtr         head = NULL, vnp;
  Boolean            any_selected = FALSE;

  dlg = (SelectionDialogPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    if (dlg->multi_select_dlg != NULL)
    {
      vnp = DialogToPointer (dlg->multi_select_dlg);
      if (vnp != NULL)
      {
        any_selected = TRUE;
        vnp = ValNodeFree (vnp);
      }
    }
    else if (dlg->list_ctrl != NULL)
    {
      if (GetValue (dlg->list_ctrl) > 0)
      {
        any_selected = TRUE;
      }
    }
    else if (dlg->popup_ctrl != NULL)
    {
      if (GetValue (dlg->popup_ctrl) > 0)
      {
        any_selected = TRUE;
      }
    }
    if (!any_selected)
    {
      head = AddStringToValNodeChain (head, dlg->err_msg, 1);
    }
  }
  return head;
}

/* err_msg is the message to put in the results from TestDialog if nothing is selected */
/* choice_list should be a valnode list of strings to use for the names of the choices. */
/* All is automatically included as a choice if allow_multi is true. */
/* The ValNodeList returned is a list of integers indicating the position of the item
 * in the list - 1 is the first item, 2 is the second item, etc. */
extern DialoG SelectionDialogEx 
(GrouP h,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Boolean                  allow_multi,
 CharPtr                  err_msg,
 ValNodePtr               choice_list,
 Int2                     list_height,
 Boolean                  force_list)

{
  SelectionDialogPtr  dlg;
  GrouP               p;
  ValNodePtr          vnp;
  Int4                num_choices;
  
  if (choice_list == NULL)
  {
    return NULL;
  }
  
  dlg = (SelectionDialogPtr) MemNew (sizeof (SelectionDialogData));
  if (dlg == NULL)
  {
    return NULL;
  }
  
  p = HiddenGroup (h, 0, 2, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = SelectionListToSelectionDialog;
  dlg->fromdialog = SelectionDialogToSelectionList;
  dlg->dialogmessage = SelectionDialogMessage;
  dlg->testdialog = TestSelectionDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;
  dlg->err_msg = err_msg;
  
  num_choices = ValNodeLen (choice_list);
  
  if (allow_multi)
  {
    dlg->multi_select_dlg = MultiSelectDialog (p, choice_list, list_height,
                                               change_notify, change_userdata);
    dlg->num_choices = num_choices + 1;                                               
  }
  else
  {
    if (num_choices < 20 && ! force_list)
    {
      dlg->popup_ctrl = PopupList (p, TRUE, SelectionDialogPopupChanged);
      SetObjectExtra (dlg->popup_ctrl, dlg, NULL);
      for (vnp = choice_list; vnp != NULL; vnp = vnp->next) {
        PopupItem (dlg->popup_ctrl, vnp->data.ptrvalue);
      }
    }
    else
    {
      dlg->list_ctrl = SingleList (p, 8, list_height, SelectionDialogChanged);
      SetObjectExtra (dlg->list_ctrl, dlg, NULL);
      for (vnp = choice_list; vnp != NULL; vnp = vnp->next) {
        ListItem (dlg->list_ctrl, vnp->data.ptrvalue);
      }      
    }
    dlg->num_choices = num_choices;
  }
  
  return (DialoG) p;
}

extern DialoG SelectionDialog 
(GrouP h,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Boolean                  allow_multi,
 CharPtr                  err_msg,
 ValNodePtr               choice_list,
 Int2                     list_height)
{
  return SelectionDialogEx (h, change_notify, change_userdata, allow_multi,
                          err_msg, choice_list, list_height, FALSE);
}

typedef struct valnodeselection
{
  DIALOG_MESSAGE_BLOCK
  DialoG           list_dlg;
  ValNodePtr       choice_list;
  
  FreeValNodeProc     free_vn_proc;
  CopyValNodeDataProc copy_vn_proc;
  MatchValNodeProc    match_vn_proc;
  RemapValNodeProc    remap_vn_proc;
  
} ValNodeSelectionData, PNTR ValNodeSelectionPtr;

static void ValNodeSelectionListToDialog (DialoG d, Pointer userdata)
{
  ValNodeSelectionPtr dlg;
  ValNodePtr          item_list, vnp_list, vnp_sel, pos_list = NULL;
  Int4                i;
  Boolean             found;

  dlg = (ValNodeSelectionPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  /* reset list control */
  PointerToDialog (dlg->list_dlg, NULL);  
  
  item_list = (ValNodePtr) userdata;
  for (vnp_list = item_list; vnp_list != NULL; vnp_list = vnp_list->next)
  {
    found = FALSE;
    for (vnp_sel = dlg->choice_list, i = 1;
         vnp_sel != NULL && !found;
         vnp_sel = vnp_sel->next, i++)
    {
      if ((dlg->match_vn_proc)(vnp_sel, vnp_list))
      {
        found = TRUE;
        ValNodeAddInt (&pos_list, 0, i);
      }
    }
  }
  PointerToDialog (dlg->list_dlg, pos_list);
  ValNodeFree (pos_list);  
}

static Pointer ValNodeSelectionDialogToList (DialoG d)
{
  ValNodeSelectionPtr dlg;
  ValNodePtr          item_list = NULL, vnp_list, pos_list, vnp_pos;
  ValNodePtr          vnp_copy, vnp_last = NULL, vnp_test;
  Int4                i;
  Boolean             found;

  dlg = (ValNodeSelectionPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  
  pos_list = DialogToPointer (dlg->list_dlg);
  for (vnp_pos = pos_list; vnp_pos != NULL; vnp_pos = vnp_pos->next)
  {
    for (i = 1, vnp_list = dlg->choice_list;
         i < vnp_pos->data.intvalue && vnp_list != NULL;
         i++, vnp_list = vnp_list->next)
    {
    }
    if (i == vnp_pos->data.intvalue && vnp_list != NULL)
    {
      /* make sure we don't already have this value in the list */
      for (vnp_test = item_list, found = FALSE;
           vnp_test != NULL && !found;
           vnp_test = vnp_test->next)
      {
        found = (dlg->match_vn_proc) (vnp_list, vnp_test);
      }
      
      if (found)
      {
        continue;
      }
      vnp_copy = (dlg->copy_vn_proc) (vnp_list);
      if (vnp_last == NULL)
      {
        item_list = vnp_copy;
      }
      else
      {
        vnp_last->next = vnp_copy;
      }
      vnp_last = vnp_copy;
    }
  }
  if (dlg->remap_vn_proc != NULL)
  {
    item_list = (dlg->remap_vn_proc) (item_list);
  }
  return item_list;  
}

static void CleanupValNodeSelectionDialogForm (GraphiC g, VoidPtr data)

{
  ValNodeSelectionPtr dlg;
  ValNodePtr          vnp;

  dlg = (ValNodeSelectionPtr) data;
  if (dlg != NULL) {
    if (dlg->free_vn_proc != NULL)
    {
      for (vnp = dlg->choice_list; vnp != NULL; vnp = vnp->next)
      {
        (dlg->free_vn_proc) (vnp);
      }
    }
    dlg->choice_list = ValNodeFree (dlg->choice_list);
  }
  StdCleanupExtraProc (g, data);
}

static void ValNodeSelectionDialogMessage (DialoG d, Int2 mssg)

{
  ValNodeSelectionPtr dlg;

  dlg = (ValNodeSelectionPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        /* reset list */
        PointerToDialog (dlg->list_dlg, NULL);
        break;
      case VIB_MSG_SELECT:
        Select (dlg->list_dlg);
        break;
      case VIB_MSG_ENTER :
        Select (dlg->list_dlg);
        break;
      case NUM_VIB_MSG + 1:
        SendMessageToDialog (dlg->list_dlg, NUM_VIB_MSG + 1);
        break;
      default :
        break;
    }
  }
}

static ValNodePtr TestValNodeSelectionDialog (DialoG d)

{
  ValNodeSelectionPtr  dlg;
  ValNodePtr           head = NULL;

  dlg = (ValNodeSelectionPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    head = TestDialog (dlg->list_dlg);
  }
  return head;
}

extern DialoG ValNodeSelectionDialogEx
(GrouP h,
 ValNodePtr               choice_list,
 Int2                     list_height,
 NameFromValNodeProc      name_proc,
 FreeValNodeProc          free_vn_proc,
 CopyValNodeDataProc      copy_vn_proc,
 MatchValNodeProc         match_vn_proc,
 CharPtr                  err_name,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Boolean                  allow_multi,
 Boolean                  force_list,
 RemapValNodeProc         remap_vn_proc)
{
  ValNodeSelectionPtr  dlg;
  GrouP                p;
  ValNodePtr           choice_name_list = NULL, vnp;

  if (choice_list == NULL || name_proc == NULL
      || copy_vn_proc == NULL || match_vn_proc == NULL)
  {
    return NULL;
  }
  
  dlg = (ValNodeSelectionPtr) MemNew (sizeof (ValNodeSelectionData));
  if (dlg == NULL)
  {
    return NULL;
  }
  
  p = HiddenGroup (h, 1, 0, NULL);
  SetObjectExtra (p, dlg, CleanupValNodeSelectionDialogForm);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = ValNodeSelectionListToDialog;
  dlg->fromdialog = ValNodeSelectionDialogToList;
  dlg->dialogmessage = ValNodeSelectionDialogMessage;
  dlg->testdialog = TestValNodeSelectionDialog;
  
  dlg->choice_list = choice_list;
  dlg->free_vn_proc = free_vn_proc;
  dlg->copy_vn_proc = copy_vn_proc;
  dlg->match_vn_proc = match_vn_proc;
  dlg->remap_vn_proc = remap_vn_proc;

  for (vnp = choice_list; vnp != NULL; vnp = vnp->next)
  {
    ValNodeAddPointer (&choice_name_list, 0, (name_proc) (vnp));
  }

  dlg->list_dlg = SelectionDialogEx (p, change_notify, change_userdata,
                                   allow_multi, err_name, choice_name_list, 
                                   list_height, force_list);
  ValNodeFreeData (choice_name_list);  
  
  return (DialoG) p;
}

extern DialoG ValNodeSelectionDialog
(GrouP h,
 ValNodePtr               choice_list,
 Int2                     list_height,
 NameFromValNodeProc      name_proc,
 FreeValNodeProc          free_vn_proc,
 CopyValNodeDataProc      copy_vn_proc,
 MatchValNodeProc         match_vn_proc,
 CharPtr                  err_name,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Boolean                  allow_multi)
{
  return ValNodeSelectionDialogEx (h, choice_list, list_height,
                                   name_proc, free_vn_proc,
                                   copy_vn_proc, match_vn_proc, err_name,
                                   change_notify, change_userdata,
                                   allow_multi, FALSE, NULL);
}

extern DialoG EnumAssocSelectionDialog 
(GrouP                 h,
 Nlm_EnumFieldAssocPtr eap,
 CharPtr               err_name,
 Boolean               allow_multi,
 Nlm_ChangeNotifyProc  change_notify,
 Pointer               change_userdata)

{
  DialoG     dlg;
  ValNodePtr choice_list = NULL;
  
  if (eap == NULL)
  {
    return NULL;
  }

  while (eap->name != NULL)
  {
    if (!StringHasNoText (eap->name))
    {
      ValNodeAddPointer (&choice_list, eap->value, StringSave (eap->name));
    }
    eap++;
  }
  
  /* note - the ValNodeSelectionDialog will free the qual_choice_list when done */                                            
  dlg = ValNodeSelectionDialog (h, choice_list, TALL_SELECTION_LIST, ValNodeStringName,
                                ValNodeSimpleDataFree, ValNodeStringCopy,
                                ValNodeChoiceMatch, err_name, 
                                change_notify, change_userdata, allow_multi);

  return dlg;
}

extern CharPtr ValNodeStringName (ValNodePtr vnp)
{
  if (vnp == NULL || vnp->data.ptrvalue == NULL)
  {
    return NULL;
  }
  else
  {
    return StringSave (vnp->data.ptrvalue);
  }
}

extern void ValNodeSimpleDataFree (ValNodePtr vnp)
{
  if (vnp != NULL && vnp->data.ptrvalue != NULL)
  {
    vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
  }
}

extern ValNodePtr ValNodeStringCopy (ValNodePtr vnp)
{
  ValNodePtr vnp_copy = NULL;
  if (vnp != NULL)
  {
    ValNodeAddPointer (&vnp_copy, vnp->choice, StringSave (vnp->data.ptrvalue));
  }
  return vnp_copy;
}

extern Boolean ValNodeChoiceMatch (ValNodePtr vnp1, ValNodePtr vnp2)
{
  if (vnp1 == NULL || vnp2 == NULL)
  {
    return FALSE;
  }
  if (vnp1->choice == vnp2->choice)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

extern Boolean ValNodeStringMatch (ValNodePtr vnp1, ValNodePtr vnp2)
{
  if (vnp1 == NULL || vnp2 == NULL)
  {
    return FALSE;
  }
  if (StringCmp (vnp1->data.ptrvalue, vnp2->data.ptrvalue) == 0)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

typedef struct sequenceselection
{
  DIALOG_MESSAGE_BLOCK
  DialoG     sequence_list_dlg;
  ValNodePtr sequence_choice_list;
} SequenceSelectionData, PNTR SequenceSelectionPtr;

static void CleanupSequenceSelectionDialogForm (GraphiC g, VoidPtr data)

{
  SequenceSelectionPtr dlg;

  dlg = (SequenceSelectionPtr) data;
  if (dlg != NULL) {
    dlg->sequence_choice_list = ValNodeFree (dlg->sequence_choice_list);
  }
  StdCleanupExtraProc (g, data);
}

static void ResetSequenceSelectionDialog (SequenceSelectionPtr dlg)
{  
  if (dlg != NULL)
  {
    PointerToDialog (dlg->sequence_list_dlg, NULL);
  }
}

static void SequenceSelectionListToSequenceSelectionDialog (DialoG d, Pointer userdata)
{
  SequenceSelectionPtr dlg;
  ValNodePtr           sequence_list, vnp_list, vnp_sel, pos_list = NULL;
  Int4                 i;
  SeqIdPtr             sip;
  Boolean              found;

  dlg = (SequenceSelectionPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  ResetSequenceSelectionDialog (dlg);  
  sequence_list = (ValNodePtr) userdata;
  for (vnp_list = sequence_list; vnp_list != NULL; vnp_list = vnp_list->next)
  {
    sip = (SeqIdPtr) vnp_list->data.ptrvalue;
    found = FALSE;
    while (sip != NULL && ! found)
    {
      for (vnp_sel = dlg->sequence_choice_list, i = 1;
           vnp_sel != NULL && !found;
           vnp_sel = vnp_sel->next, i++)
      {
        found = SeqIdIn (sip, vnp_sel->data.ptrvalue);
        if (found)
        {
          ValNodeAddInt (&pos_list, 0, i);
        }
      }
      sip = sip->next;
    }
  }
  PointerToDialog (dlg->sequence_list_dlg, pos_list);
  ValNodeFree (pos_list);
}

static Pointer SequenceSelectionDialogToSequenceSelectionList (DialoG d)
{
  SequenceSelectionPtr dlg;
  ValNodePtr           sequence_list = NULL, vnp_list, pos_list, vnp_pos;
  Int4                 i;

  dlg = (SequenceSelectionPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  
  pos_list = DialogToPointer (dlg->sequence_list_dlg);
  for (vnp_pos = pos_list; vnp_pos != NULL; vnp_pos = vnp_pos->next)
  {
    for (i = 1, vnp_list = dlg->sequence_choice_list;
         i < vnp_pos->data.intvalue && vnp_list != NULL;
         i++, vnp_list = vnp_list->next)
    {
    }
    if (i == vnp_pos->data.intvalue && vnp_list != NULL)
    {
      ValNodeAddPointer (&sequence_list, 0, vnp_list->data.ptrvalue);
    }
  }
  return sequence_list;
}

static void 
GetSequenceChoiceList 
(SeqEntryPtr sep,
 ValNodePtr PNTR list, 
 Boolean show_nucs, 
 Boolean show_prots)
{
  BioseqPtr                bsp;
  BioseqSetPtr             bssp;
  
  if (sep == NULL) return;
  
  if (IS_Bioseq (sep))
  {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return;
    if (!show_nucs && ISA_na (bsp->mol))
    {
      return;
    }
    if (!show_prots && ISA_aa (bsp->mol))
    {
      return;
    }
    ValNodeAddPointer (list, 0, bsp->id);
  }
  else
  {
  	bssp = (BioseqSetPtr) sep->data.ptrvalue;
    for (sep = bssp->seq_set; sep != NULL; sep = sep->next) 
    {
      GetSequenceChoiceList (sep, list, show_nucs, show_prots);
    }
  }
}

static void SequenceSelectionDialogMessage (DialoG d, Int2 mssg)

{
  SequenceSelectionPtr dlg;

  dlg = (SequenceSelectionPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        /* reset list */
        ResetSequenceSelectionDialog (dlg);
        break;
      case VIB_MSG_SELECT:
        Select (dlg->sequence_list_dlg);
        break;
      case VIB_MSG_ENTER :
        Select (dlg->sequence_list_dlg);
        break;
      case NUM_VIB_MSG + 1:
        SendMessageToDialog (dlg->sequence_list_dlg, NUM_VIB_MSG + 1);
        break;
      default :
        break;
    }
  }
}

static ValNodePtr TestSequenceSelectionDialog (DialoG d)

{
  SequenceSelectionPtr dlg;
  ValNodePtr           head = NULL;

  dlg = (SequenceSelectionPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    head = TestDialog (dlg->sequence_list_dlg);
  }
  return head;
}


extern DialoG SequenceSelectionDialog 
(GrouP h,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Boolean                  allow_multi,
 Boolean                  show_nucs,
 Boolean                  show_prots,
 Uint2                    entityID)

{
  SequenceSelectionPtr  dlg;
  GrouP                 p;
  ValNodePtr                vnp;
  SeqEntryPtr               sep;
  SeqIdPtr                  sip;
  Char                      tmp[128];
  ValNodePtr            choice_name_list = NULL;
  
  if (!show_nucs && ! show_prots)
  {
    return NULL;
  }
  
  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL)
  {
    return NULL;
  }

  dlg = (SequenceSelectionPtr) MemNew (sizeof (SequenceSelectionData));
  if (dlg == NULL)
  {
    return NULL;
  }
  
  p = HiddenGroup (h, 1, 0, NULL);
  SetObjectExtra (p, dlg, CleanupSequenceSelectionDialogForm);

  dlg->dialog = (DialoG) p;
  dlg->todialog = SequenceSelectionListToSequenceSelectionDialog;
  dlg->fromdialog = SequenceSelectionDialogToSequenceSelectionList;
  dlg->dialogmessage = SequenceSelectionDialogMessage;
  dlg->testdialog = TestSequenceSelectionDialog;

  dlg->sequence_choice_list = NULL;
  GetSequenceChoiceList (sep, &dlg->sequence_choice_list, show_nucs, show_prots);
  
  
  for (vnp = dlg->sequence_choice_list; vnp != NULL; vnp = vnp->next) {
    sip = SeqIdFindWorst ((SeqIdPtr) vnp->data.ptrvalue);
    SeqIdWrite (sip, tmp, PRINTID_REPORT, sizeof (tmp));
    ValNodeAddPointer (&choice_name_list, 0, StringSave (tmp));
  }

  dlg->sequence_list_dlg = SelectionDialog (p, change_notify, change_userdata,
                                            allow_multi, "sequence",
                                            choice_name_list, TALL_SELECTION_LIST);
  ValNodeFreeData (choice_name_list); 
  return (DialoG) p;
}

/*
static CharPtr inferencePrefix [] = {
  "",
  "similar to sequence",
  "similar to AA sequence",
  "similar to DNA sequence",
  "similar to RNA sequence",
  "similar to RNA sequence, mRNA",
  "similar to RNA sequence, EST",
  "similar to RNA sequence, other RNA",
  "profile",
  "nucleotide motif",
  "protein motif",
  "ab initio prediction",
  NULL
};

ENUM_ALIST(inference_alist)
  { " ",                     0 },
  { "similar to sequence",   1 },
  { "similar to protein",    2 },
  { "similar to DNA",        3 },
  { "similar to RNA",        4 },
  { "similar to mRNA",       5 },
  { "similar to EST",        6 },
  { "similar to other RNA",  7 },
  { "profile",               8 },
  { "nucleotide motif",      9 },
  { "protein motif",        10 },
  { "ab initio prediction", 11 },
END_ENUM_ALIST

Uint2 inference_types [] = {
  TAGLIST_POPUP, TAGLIST_TEXT
};

Uint2 inference_widths [] = {
  0, 0
};

static EnumFieldAssocPtr inference_popups [] = {
  inference_alist, NULL
};

extern void GBQualsToInferenceDialog (DialoG d, SeqFeatPtr sfp)

{
  Int2               best;
  Char               ch;
  GBQualPtr          gbq;
  ValNodePtr         head = NULL;
  Int2               j;
  ValNodePtr         last = NULL;
  size_t             len;
  CharPtr            rest;
  CharPtr            str;
  TagListPtr         tlp;
  Char               tmp [32];
  ValNodePtr         vnp;

  tlp = (TagListPtr) GetObjectExtra (d);
  if (tlp == NULL) return;

  if (sfp != NULL) {
    for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
      if (StringICmp (gbq->qual, "inference") != 0) continue;
      if (StringHasNoText (gbq->val)) continue;
      vnp = ValNodeNew (last);
      if (vnp == NULL) continue;
      if (head == NULL) {
        head = vnp;
      }
      last = vnp;

      rest = NULL;
      best = -1;
      for (j = 0; inferencePrefix [j] != NULL; j++) {
        len = StringLen (inferencePrefix [j]);
        if (StringNICmp (gbq->val, inferencePrefix [j], len) != 0) continue;
        rest = gbq->val + len;
        best = j;
      }
      if (best >= 0 && inferencePrefix [best] != NULL) {
        if (rest != NULL) {
          ch = *rest;
          while (IS_WHITESP (ch) || ch == ':') {
            rest++;
            ch = *rest;
          }
        }
        len = StringLen (rest);
        str = MemNew (len + 16);
        if (str != NULL) {
          sprintf (tmp, "%d", (int) best);
          StringCpy (str, tmp);
          StringCat (str, "\t");
          StringCat (str, rest);
          StringCat (str, "\n");
        }
        vnp->data.ptrvalue = str;
      } else {
        len + StringLen (gbq->val);
        str = MemNew (len + 8);
        if (str != NULL) {
          StringCpy (str, "0");
          StringCat (str, "\t");
          StringCat (str, gbq->val);
          StringCat (str, "\n");
        }
        vnp->data.ptrvalue = str;
      }
    }
  }

  SendMessageToDialog (tlp->dialog, VIB_MSG_RESET);
  tlp->vnp = head;
  SendMessageToDialog (tlp->dialog, VIB_MSG_REDRAW);
  for (j = 0, vnp = tlp->vnp; vnp != NULL; j++, vnp = vnp->next) continue;
  tlp->max = MAX ((Int2) 0, (Int2) (j - tlp->rows + 1));
  CorrectBarMax (tlp->bar, tlp->max);
  CorrectBarPage (tlp->bar, tlp->rows - 1, tlp->rows - 1);
}

static void VisStringDialogToGbquals (SeqFeatPtr sfp, DialoG d, CharPtr qual)

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

extern void InferenceDialogToGBQuals (DialoG d, SeqFeatPtr sfp)

{
  GBQualPtr   gbq;
  GBQualPtr   gbqlast = NULL;
  Int2        j;
  size_t      len;
  CharPtr     prefix;
  CharPtr     ptr;
  CharPtr     rest;
  CharPtr     str;
  TagListPtr  tlp;
  Int2        val;
  ValNodePtr  vnp;

  tlp = (TagListPtr) GetObjectExtra (d);
  if (tlp == NULL || sfp == NULL) return;

  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    gbqlast = gbq;
  }

  for (vnp = tlp->vnp; vnp != NULL; vnp = vnp->next) {
    if (StringHasNoText ((CharPtr) vnp->data.ptrvalue)) continue;
    ptr = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, 0);
    TrimSpacesAroundString (ptr);
    prefix = NULL;
    if (StrToInt (ptr, &val)) {
      for (j = 0; inferencePrefix [j] != NULL; j++) {
        if (j == val) {
          prefix = inferencePrefix [j];
        }
      }
    }
    MemFree (ptr);
    rest = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, 1);
    TrimSpacesAroundString (rest);
    if (StringDoesHaveText (prefix)) {
      len = StringLen (prefix) + StringLen (rest);
      str = (CharPtr) MemNew (len + 8);
      if (str != NULL) {
        if (StringDoesHaveText (prefix)) {
          StringCpy (str, prefix);
          if (StringDoesHaveText (rest)) {
            if (StringNICmp (rest, "(same species)", 14) != 0) {
              StringCat (str, ":");
            } else {
              StringCat (str, " ");
            }
          }
        }
        if (StringDoesHaveText (rest)) {
          StringCat (str, rest);
        }
        gbq = GBQualNew ();
        if (gbq != NULL) {
          gbq->qual = StringSave ("inference");
          gbq->val = str;
          if (gbqlast == NULL) {
            sfp->qual = gbq;
          } else {
            gbqlast->next = gbq;
          }
          gbqlast = gbq;
        }
      }
    }
    MemFree (rest);
  }
}

static DialoG CreateInferenceDialog (GrouP h, Uint2 rows, Int2 spacing, Int2 width)

{
  inference_widths [1] = width;
  return CreateTagListDialog (h, rows, 2, spacing,
                              inference_types, inference_widths,
                              inference_popups, NULL, NULL);
}
*/

/* ************************ */

/* inference dialog controls, utility functions */

typedef struct inferevid {
  CharPtr  prefix;     /* from inferencePrefix     */
  Boolean  species;    /* optional (same species)  */
  CharPtr  database;   /* INSD, RefSeq, etc.       */
  CharPtr  db_other;   /* other database           */
  CharPtr  accession;  /* accession.version        */
  CharPtr  program;    /* common analysis program  */
  CharPtr  pr_other;   /* other program            */
  CharPtr  version;    /* program version          */
  CharPtr  basis1;     /* profile or motif         */
  CharPtr  basis2;     /*  evidence_basis texts    */
} InferEvid, PNTR InferEvidPtr;

typedef struct inferdialog {
  DIALOG_MESSAGE_BLOCK

  DoC           inferdoc;
  Int2          currItem;

  PopuP         prefix;
  ButtoN        species;
  PopuP         database;
  TexT          db_other;
  TexT          accession;
  PopuP         program;
  TexT          pr_other;
  TexT          version;
  TexT          basis1;
  TexT          basis2;

  GrouP         inf_accn_group;
  GrouP         other_db_group;
  GrouP         inf_prog_group;
  GrouP         other_pr_group;
  GrouP         inf_free_group;

  Int2          numInf;
  InferEvidPtr  evidence [128];

} InferDialog, PNTR InferDialogPtr;

static InferEvidPtr InferEvidNew (
  void
)

{
  InferEvidPtr  iep;

  iep = MemNew (sizeof (InferEvid));
  if (iep == NULL) return NULL;

  return iep;
}

static InferEvidPtr InferEvidFree (
  InferEvidPtr iep
)

{
  if (iep == NULL) return NULL;

  MemFree (iep->prefix);
  MemFree (iep->database);
  MemFree (iep->accession);
  MemFree (iep->program);
  MemFree (iep->version);
  MemFree (iep->basis1);
  MemFree (iep->basis2);

  return MemFree (iep);
}

static InferEvidPtr GetInferEvid (
  InferDialogPtr idp,
  Int2 item
)

{
  InferEvidPtr  iep;

  if (idp == NULL || item < 0 || item > 127) return NULL;
  iep = idp->evidence [item];
  if (iep != NULL) return iep;

  iep = InferEvidNew ();
  if (iep != NULL) {
    /*
    iep->prefix = StringSave (" ");
    iep->database = StringSave (" ");
    iep->db_other = StringSave ("");
    iep->accession = StringSave ("");
    iep->program = StringSave (" ");
    iep->pr_other = StringSave ("");
    iep->version = StringSave ("");
    iep->basis1 = StringSave ("");
    iep->basis2 = StringSave ("");
    */
  }
  idp->evidence [item] = iep;
  return iep;
}

/* inference DoC object tables */

#define NUM_INFERENCE_LINES 3

static ParData  inferParFmt = { FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0 };

static ColData  inferColFmt [] = {
  {0, 5, 25, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, FALSE}, /* class     */
  {0, 5, 25, 2, NULL, 'l', FALSE, TRUE,  FALSE, FALSE, TRUE}   /* specifics */
};

static CharPtr inferencePrefix [] = {
  "",
  "similar to sequence",
  "similar to AA sequence",
  "similar to DNA sequence",
  "similar to RNA sequence",
  "similar to RNA sequence, mRNA",
  "similar to RNA sequence, EST",
  "similar to RNA sequence, other RNA",
  "profile",
  "nucleotide motif",
  "protein motif",
  "ab initio prediction",
  NULL
};

ENUM_ALIST(inference_alist)
  { " ",                     0 },
  { "similar to sequence",   1 },
  { "similar to protein",    2 },
  { "similar to DNA",        3 },
  { "similar to RNA",        4 },
  { "similar to mRNA",       5 },
  { "similar to EST",        6 },
  { "similar to other RNA",  7 },
  { "profile",               8 },
  { "nucleotide motif",      9 },
  { "protein motif",        10 },
  { "ab initio prediction", 11 },
END_ENUM_ALIST

static CharPtr accnTypePrefix [] = {
  "",
  "GenBank",
  "EMBL",
  "DDBJ",
  "INSD",
  "RefSeq",
  "UniProt",
  "?",
  NULL
};

ENUM_ALIST(accn_type_alist)
  { " ",       0 },
  { "GenBank", 1 },
  { "EMBL",    2 },
  { "DDBJ",    3 },
  { "INSD",    4 },
  { "RefSeq",  5 },
  { "UniProt", 6 },
  { "Other",   7 },
END_ENUM_ALIST

static CharPtr programPrefix [] = {
  "",
  "tRNAscan",
  "Genscan",
  "?",
  NULL
};

ENUM_ALIST(program_alist)
  { " ",        0 },
  { "tRNAscan", 1 },
  { "Genscan",  2 },
  { "Other",    3 },
END_ENUM_ALIST

static CharPtr PrintInferTable (
  DoC d,
  Int2 item,
  Pointer data
)

{
  Char            buf [256];
  InferDialogPtr  idp;
  InferEvidPtr    iep;

  idp = (InferDialogPtr) GetObjectExtra (d);
  if (idp == NULL || item < 1 || item > 127) return NULL;
  iep = GetInferEvid (idp, item);
  if (iep == NULL) return NULL;

  buf [0] = '\0';

  if (StringHasNoText (iep->prefix)) {
    StringCat (buf, " \t \n");
    return StringSave (buf);
  }

  StringCat (buf, iep->prefix);

  StringCat (buf, "\t");

  if (StringNICmp (iep->prefix, "similar to ", 11) == 0) {
    if (StringDoesHaveText (iep->accession)) {
      if (StringCmp (iep->database, "Other") == 0) {
        if (StringDoesHaveText (iep->db_other)) {
          StringCat (buf, iep->db_other);
          StringCat (buf, ":");
        }
      } else if (StringDoesHaveText (iep->database)) {
        StringCat (buf, iep->database);
        StringCat (buf, ":");
      }
      StringCat (buf, iep->accession);
    }
  } else if (StringNICmp (iep->prefix, "ab initio ", 10) == 0) {
    if (StringCmp (iep->program, "Other") == 0) {
      if (StringDoesHaveText (iep->pr_other)) {
        StringCat (buf, iep->pr_other);
        if (StringDoesHaveText (iep->version)) {
          StringCat (buf, ":");
          StringCat (buf, iep->version);
        }
      }
    } else if (StringDoesHaveText (iep->program)) {
      StringCat (buf, iep->program);
      if (StringDoesHaveText (iep->version)) {
        StringCat (buf, ":");
        StringCat (buf, iep->version);
      }
    }
  } else if (StringDoesHaveText (iep->basis1)) {
    StringCat (buf, iep->basis1);
    if (StringDoesHaveText (iep->basis2)) {
      StringCat (buf, ":");
      StringCat (buf, iep->basis2);
    }
  } else {
    StringCat (buf, " ");
  }

  StringCat (buf, "\n");
  return StringSave (buf);
}

static void ShowInferenceGroup (
  InferDialogPtr idp
)

{
  CharPtr  str;
  UIEnum   val;

  if (idp == NULL) return;
  if (GetEnumPopup (idp->prefix, inference_alist, &val)) {
    if (val >= 1 && val <= 7) {
      SafeHide (idp->inf_prog_group);
      SafeHide (idp->inf_free_group);
      SafeShow (idp->inf_accn_group);
      SafeShow (idp->species);
      str = GetEnumPopupByName (idp->database, accn_type_alist);
      if (StringCmp (str, "Other") == 0) {
        SafeShow (idp->other_db_group);
      } else {
        SafeHide (idp->other_db_group);
      }
      MemFree (str);
    } else if (val >= 8 && val <= 10) {
      SafeHide (idp->inf_accn_group);
      SafeHide (idp->species);
      SafeHide (idp->inf_prog_group);
      SafeShow (idp->inf_free_group);
    } else if (val == 11) {
      SafeHide (idp->inf_accn_group);
      SafeHide (idp->species);
      SafeHide (idp->inf_free_group);
      SafeShow (idp->inf_prog_group);
      str = GetEnumPopupByName (idp->program, program_alist);
      if (StringCmp (str, "Other") == 0) {
        SafeShow (idp->other_pr_group);
      } else {
        SafeHide (idp->other_pr_group);
      }
      MemFree (str);
    } else {
      SafeHide (idp->inf_accn_group);
      SafeHide (idp->species);
      SafeHide (idp->inf_prog_group);
      SafeHide (idp->inf_free_group);
    }
  } else {
    SafeHide (idp->inf_accn_group);
    SafeHide (idp->species);
    SafeHide (idp->inf_prog_group);
    SafeHide (idp->inf_free_group);
  }
  Update ();
}

static void SafeSetEnumPopupByName (PopuP lst, EnumFieldAssocPtr al, CharPtr name)

{
  if (StringDoesHaveText (name)) {
    SetEnumPopupByName (lst, al, name);
  } else {
    SetEnumPopupByName (lst, al, " ");
  }
}

static void ChangeInferTableSelect (
  DoC d,
  Int2 item,
  Int2 row,
  Int2 col,
  Boolean dblClck
)

{
  InferDialogPtr  idp;
  InferEvidPtr    iep;
  Int2            itemOld1, itemOld2;

  idp = (InferDialogPtr) GetObjectExtra (d);
  if (idp == NULL) return;
  if (item == 0 || row == 0 || col == 0) return;

  GetDocHighlight (d, &itemOld1, &itemOld2);
  SetDocHighlight (d, item, item);
  UpdateDocument (d, itemOld1, itemOld2);
  UpdateDocument (d, item, item);
  idp->currItem = item;

  iep = GetInferEvid (idp, item);
  if (iep != NULL) {
    ResetClip ();
    SafeSetEnumPopupByName (idp->prefix, inference_alist, iep->prefix);

    SafeSetStatus (idp->species, iep->species);
    SafeSetEnumPopupByName (idp->database, accn_type_alist, iep->database);
    SafeSetTitle (idp->db_other, iep->db_other);
    SafeSetTitle (idp->accession, iep->accession);

    SafeSetEnumPopupByName (idp->program, program_alist, iep->program);
    SafeSetTitle (idp->pr_other, iep->pr_other);
    SafeSetTitle (idp->version, iep->version);

    SafeSetTitle (idp->basis1, iep->basis1);
    SafeSetTitle (idp->basis2, iep->basis2);

    ShowInferenceGroup (idp);
  }

  Update ();
}

static void CheckExtendInferTable (
  InferDialogPtr idp
)

{
  Int2  numItems;

  if (idp == NULL) return;

  GetDocParams (idp->inferdoc, &numItems, NULL);
  if (idp->currItem == numItems) {
    AppendItem (idp->inferdoc, PrintInferTable, idp, FALSE, 1,
                &inferParFmt, inferColFmt, systemFont);
  }

  Update ();
}

static void ChangeInferPrefix (
  PopuP p
)

{
  AlistDialogPtr  adp;
  InferDialogPtr  idp;
  InferEvidPtr    iep;
  CharPtr         str;

  adp = (AlistDialogPtr) GetObjectExtra (p);
  if (adp == NULL) return;
  idp = (InferDialogPtr) adp->userdata;
  if (idp == NULL) return;
  iep = GetInferEvid (idp, idp->currItem);
  if (iep == NULL) return;

  str = GetEnumPopupByName (idp->prefix, inference_alist);
  iep->prefix = MemFree (iep->prefix);
  iep->prefix = str; /* allocated by GetEnumPopupByName */

  ShowInferenceGroup (idp);

  UpdateDocument (idp->inferdoc, idp->currItem, idp->currItem);
  Update ();

  CheckExtendInferTable (idp);
}

static void ChangeSameSpecies (
  ButtoN b
)

{
  InferDialogPtr  idp;
  InferEvidPtr    iep;

  idp = (InferDialogPtr) GetObjectExtra (b);
  if (idp == NULL) return;
  iep = GetInferEvid (idp, idp->currItem);
  if (iep == NULL) return;

  iep->species = (Boolean) (GetStatus (b));

  ShowInferenceGroup (idp);

  UpdateDocument (idp->inferdoc, idp->currItem, idp->currItem);
  Update ();

  CheckExtendInferTable (idp);
}

static CharPtr insdmessage =
"GenBank, EMBL, and DDBJ records are part of the International Nucleotide " \
"Sequence Database collaboration.\nThe database prefix for the /inference " \
"qualifier in these cases is INSD by collaboration policy.";

static void ChangeInferDatabase (
  PopuP p
)


{
  AlistDialogPtr  adp;
  InferDialogPtr  idp;
  InferEvidPtr    iep;
  CharPtr         str;

  adp = (AlistDialogPtr) GetObjectExtra (p);
  if (adp == NULL) return;
  idp = (InferDialogPtr) adp->userdata;
  if (idp == NULL) return;
  iep = GetInferEvid (idp, idp->currItem);
  if (iep == NULL) return;

  str = GetEnumPopupByName (idp->database, accn_type_alist);
  if (StringCmp (str, "GenBank") == 0 ||
      StringCmp (str, "EMBL") == 0 ||
      StringCmp (str, "DDBJ") == 0) {
    if (GetAppProperty ("InternalNcbiSequin") == NULL) {
      Message (MSG_OK, "%s", insdmessage);
    }
    SetEnumPopupByName (idp->database, accn_type_alist, "INSD");
    str = MemFree (str);
    str = StringSave ("INSD");
  }
  iep->database = MemFree (iep->database);
  iep->database = str; /* allocated by GetEnumPopupByName */

  ShowInferenceGroup (idp);

  UpdateDocument (idp->inferdoc, idp->currItem, idp->currItem);
  Update ();

  CheckExtendInferTable (idp);
}

static void ChangeInferDbOther (
  TexT t
)

{
  InferDialogPtr  idp;
  InferEvidPtr    iep;

  idp = (InferDialogPtr) GetObjectExtra (t);
  if (idp == NULL) return;
  iep = GetInferEvid (idp, idp->currItem);
  if (iep == NULL) return;

  iep->db_other = MemFree (iep->db_other);
  iep->db_other = SaveStringFromText (t);

  ShowInferenceGroup (idp);

  UpdateDocument (idp->inferdoc, idp->currItem, idp->currItem);
  Update ();

  CheckExtendInferTable (idp);
}

static void ChangeInferAccession (
  TexT t
)

{
  InferDialogPtr  idp;
  InferEvidPtr    iep;

  idp = (InferDialogPtr) GetObjectExtra (t);
  if (idp == NULL) return;
  iep = GetInferEvid (idp, idp->currItem);
  if (iep == NULL) return;

  iep->accession = MemFree (iep->accession);
  iep->accession = SaveStringFromText (t);

  ShowInferenceGroup (idp);

  UpdateDocument (idp->inferdoc, idp->currItem, idp->currItem);
  Update ();

  CheckExtendInferTable (idp);
}

static void ChangeInferProgram (
  PopuP p
)


{
  AlistDialogPtr  adp;
  InferDialogPtr  idp;
  InferEvidPtr    iep;
  CharPtr         str;

  adp = (AlistDialogPtr) GetObjectExtra (p);
  if (adp == NULL) return;
  idp = (InferDialogPtr) adp->userdata;
  if (idp == NULL) return;
  iep = GetInferEvid (idp, idp->currItem);
  if (iep == NULL) return;

  str = GetEnumPopupByName (idp->program, program_alist);
  iep->program = MemFree (iep->program);
  iep->program = str; /* allocated by GetEnumPopupByName */

  ShowInferenceGroup (idp);

  UpdateDocument (idp->inferdoc, idp->currItem, idp->currItem);
  Update ();

  CheckExtendInferTable (idp);
}

static void ChangeInferPrOther (
  TexT t
)

{
  InferDialogPtr  idp;
  InferEvidPtr    iep;

  idp = (InferDialogPtr) GetObjectExtra (t);
  if (idp == NULL) return;
  iep = GetInferEvid (idp, idp->currItem);
  if (iep == NULL) return;

  iep->pr_other = MemFree (iep->pr_other);
  iep->pr_other = SaveStringFromText (t);

  ShowInferenceGroup (idp);

  UpdateDocument (idp->inferdoc, idp->currItem, idp->currItem);
  Update ();

  CheckExtendInferTable (idp);
}

static void ChangeInferVersion (
  TexT t
)

{
  InferDialogPtr  idp;
  InferEvidPtr    iep;

  idp = (InferDialogPtr) GetObjectExtra (t);
  if (idp == NULL) return;
  iep = GetInferEvid (idp, idp->currItem);
  if (iep == NULL) return;

  iep->version = MemFree (iep->version);
  iep->version = SaveStringFromText (t);

  ShowInferenceGroup (idp);

  UpdateDocument (idp->inferdoc, idp->currItem, idp->currItem);
  Update ();

  CheckExtendInferTable (idp);
}

static void ChangeInferBasis1 (
  TexT t
)

{
  InferDialogPtr  idp;
  InferEvidPtr    iep;

  idp = (InferDialogPtr) GetObjectExtra (t);
  if (idp == NULL) return;
  iep = GetInferEvid (idp, idp->currItem);
  if (iep == NULL) return;

  iep->basis1 = MemFree (iep->basis1);
  iep->basis1 = SaveStringFromText (t);

  ShowInferenceGroup (idp);

  UpdateDocument (idp->inferdoc, idp->currItem, idp->currItem);
  Update ();

  CheckExtendInferTable (idp);
}

static void ChangeInferBasis2 (
  TexT t
)

{
  InferDialogPtr  idp;
  InferEvidPtr    iep;

  idp = (InferDialogPtr) GetObjectExtra (t);
  if (idp == NULL) return;
  iep = GetInferEvid (idp, idp->currItem);
  if (iep == NULL) return;

  iep->basis2 = MemFree (iep->basis2);
  iep->basis2 = SaveStringFromText (t);

  ShowInferenceGroup (idp);

  UpdateDocument (idp->inferdoc, idp->currItem, idp->currItem);
  Update ();

  CheckExtendInferTable (idp);
}

static Boolean StringInList (CharPtr str, CharPtr PNTR list)

{
  Int2  i;

  if (str == NULL || list == NULL) return FALSE;

  for (i = 0; list [i] != NULL; i++) {
    if (StringICmp (str, list[i]) == 0) return TRUE;
  }

  return FALSE;
}

extern void GBQualsToInferenceDialog (DialoG d, SeqFeatPtr sfp)

{
  Int2            best;
  Char            ch;
  GBQualPtr       gbq;
  Int2            i, j, k;
  InferDialogPtr  idp;
  InferEvidPtr    iep;
  size_t          len;
  CharPtr         rest;
  CharPtr         str;
  CharPtr         tmp;

  idp = (InferDialogPtr) GetObjectExtra (d);
  if (idp == NULL) return;

  if (sfp == NULL || sfp->qual == NULL) {
    Reset (idp->inferdoc);
    SetValue (idp->prefix, 0);
    SetStatus (idp->species, FALSE);
    SetValue (idp->database, 0);
    SetTitle (idp->db_other, "");
    SetTitle (idp->accession, "");
    SetValue (idp->program, 0);
    SetTitle (idp->pr_other, "");
    SetTitle (idp->version, "");
    SetTitle (idp->basis1, "");
    SetTitle (idp->basis2, "");
    SafeHide (idp->inf_accn_group);
    SafeHide (idp->inf_prog_group);
    SafeHide (idp->inf_free_group);
    idp->numInf = 0;
    idp->currItem = 1;
    for (i = 0; i < NUM_INFERENCE_LINES; i++) {
      AppendItem (idp->inferdoc, PrintInferTable, idp, FALSE, 1,
                  &inferParFmt, inferColFmt, systemFont);
    }
    SetDocHighlight (idp->inferdoc, 1, 1);
    return;
  }

  idp->numInf = 0;
  idp->currItem = 1;
  Reset (idp->inferdoc);

  for (k = 0; k < 128; k++) {
    iep = idp->evidence [k];
    InferEvidFree (iep);
    idp->evidence [k] = NULL;
  }

  for (gbq = sfp->qual, k = 0; gbq != NULL; gbq = gbq->next) {
    if (StringICmp (gbq->qual, "inference") != 0) continue;
    if (StringHasNoText (gbq->val)) continue;

    rest = NULL;
    best = -1;
    for (j = 0; inferencePrefix [j] != NULL; j++) {
      len = StringLen (inferencePrefix [j]);
      if (StringNICmp (gbq->val, inferencePrefix [j], len) != 0) continue;
      rest = gbq->val + len;
      best = j;
    }

    k++;
    iep = GetInferEvid (idp, k);
    if (iep == NULL) continue;

    str = NULL;
    if (best > 0 && inferencePrefix [best] != NULL) {
      iep->prefix = MemFree (iep->prefix);
      iep->prefix = StringSave(GetEnumName ((UIEnum) best, inference_alist));

      if (rest != NULL) {
        ch = *rest;
        while (IS_WHITESP (ch)) {
          rest++;
          ch = *rest;
        }
        if (StringNICmp (rest, "(same species)", 14) == 0) {
          iep->species = TRUE;
          rest += 14;
        }
        ch = *rest;
        while (IS_WHITESP (ch) || ch == ':') {
          rest++;
          ch = *rest;
        }
      }
      if (StringDoesHaveText (rest)) {
        str = StringSave (rest);
      }
      tmp = StringChr (str, ':');
      if (tmp != NULL) {
        *tmp = '\0';
        tmp++;
        TrimSpacesAroundString (str);
        TrimSpacesAroundString (tmp);
      } else {
        TrimSpacesAroundString (str);
      }
      if (StringNICmp (iep->prefix, "similar to ", 11) == 0) {
        if (StringInList (str, accnTypePrefix)) {
          iep->database = MemFree (iep->database);
          iep->database = StringSaveNoNull (str);
          iep->accession = MemFree (iep->accession);
          iep->accession = StringSaveNoNull (tmp);
        } else if (tmp != NULL) {
          iep->database = MemFree (iep->database);
          iep->database = StringSaveNoNull ("Other");
          iep->db_other = MemFree (iep->db_other);
          iep->db_other = StringSaveNoNull (str);
          iep->accession = MemFree (iep->accession);
          iep->accession = StringSaveNoNull (tmp);
        } else {
          iep->database = MemFree (iep->database);
          iep->database = StringSaveNoNull (" ");
          iep->db_other = MemFree (iep->db_other);
          iep->accession = MemFree (iep->accession);
          iep->accession = StringSaveNoNull (str);
        }
      } else if (StringNICmp (iep->prefix, "ab initio ", 10) == 0) {
        if (StringInList (str, programPrefix)) {
          iep->program = MemFree (iep->program);
          iep->program = StringSaveNoNull (str);
        } else {
          iep->program = MemFree (iep->program);
          iep->program = StringSaveNoNull ("Other");
          iep->pr_other = MemFree (iep->pr_other);
          iep->pr_other = StringSaveNoNull (str);
        }
        iep->version = MemFree (iep->version);
        iep->version = StringSaveNoNull (tmp);
      } else {
        iep->basis1 = MemFree (iep->basis1);
        iep->basis1 = StringSaveNoNull (str);
        iep->basis2 = MemFree (iep->basis2);
        iep->basis2 = StringSaveNoNull (tmp);
      }

    } else {
      iep->prefix = StringSave ("???");
      str = StringSave (gbq->val);
      tmp = StringChr (str, ':');
      if (tmp != NULL) {
        *tmp = '\0';
        tmp++;
        TrimSpacesAroundString (str);
        TrimSpacesAroundString (tmp);
      } else {
        TrimSpacesAroundString (str);
      }
      iep->basis1 = MemFree (iep->basis1);
      iep->basis1 = StringSaveNoNull (str);
      iep->basis2 = MemFree (iep->basis2);
      iep->basis2 = StringSaveNoNull (tmp);
    }

    MemFree (str);

    AppendItem (idp->inferdoc, PrintInferTable, idp, FALSE, 1,
                &inferParFmt, inferColFmt, systemFont);

    (idp->numInf)++;
  }

  AppendItem (idp->inferdoc, PrintInferTable, idp, FALSE, 1,
              &inferParFmt, inferColFmt, systemFont);
  k++;

  while (k < NUM_INFERENCE_LINES) {
    AppendItem (idp->inferdoc, PrintInferTable, idp, FALSE, 1,
                &inferParFmt, inferColFmt, systemFont);
    k++;
  }

  ShowInferenceGroup (idp);

  UpdateDocument (idp->inferdoc, 0, 0);

  ChangeInferTableSelect (idp->inferdoc, 1, 1, 1, FALSE);

  Update ();
}

extern void InferenceDialogToGBQuals (DialoG d, SeqFeatPtr sfp, Boolean convertBadToNote)

{
  CharPtr         first = NULL, second = NULL;
  GBQualPtr       gbq, lastgbq;
  InferDialogPtr  idp;
  InferEvidPtr    iep;
  Int2            k, numItems;
  size_t          len;
  CharPtr         prefix = NULL;
  CharPtr         speciesies = NULL;
  CharPtr         str;
  UIEnum          val;

  idp = (InferDialogPtr) GetObjectExtra (d);
  if (idp == NULL || sfp == NULL) return;

  lastgbq = NULL;
  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    lastgbq = gbq;
  }

  GetDocParams (idp->inferdoc, &numItems, NULL);
  for (k = 1; k <= numItems; k++) {
    iep = GetInferEvid (idp, k);
    if (iep == NULL) continue;

    if (StringHasNoText (iep->prefix)) continue;
    gbq = GBQualNew ();
    if (gbq == NULL) continue;

    gbq->qual = StringSave ("inference");

    if (WhereInEnumPopup (inference_alist, iep->prefix, &val)) {
      if (val > 0 && val <= 11) {
        prefix = inferencePrefix [(int) val];
      }
    }
    if (StringNICmp (iep->prefix, "similar to ", 11) == 0) {
      if (iep->species) {
        speciesies = " (same species)";
      }
      if (StringDoesHaveText (iep->accession)) {
        if (StringCmp (iep->database, "Other") == 0) {
          if (StringDoesHaveText (iep->db_other)) {
            first = iep->db_other;
          }
        } else if (StringDoesHaveText (iep->database)) {
          first = iep->database;
        }
        second = iep->accession;
      }
    } else if (StringNICmp (iep->prefix, "ab initio ", 10) == 0) {
      if (StringCmp (iep->program, "Other") == 0) {
        if (StringDoesHaveText (iep->pr_other)) {
          first = iep->pr_other;
        }
      } else if (StringDoesHaveText (iep->program)) {
        first = iep->program;
      }
      second = iep->version;
    } else {
      if (StringDoesHaveText (iep->basis1)) {
        first = iep->basis1;
        second = iep->basis2;
      }
    }

    len = StringLen (prefix) + StringLen (speciesies) + StringLen (first) + StringLen (second);
    str = MemNew (len + 5);
    if (str != NULL) {
      StringCpy (str, prefix);
      StringCat (str, speciesies);
      StringCat (str, ":");
      StringCat (str, first);
      StringCat (str, ":");
      StringCat (str, second);
      gbq->val = StringSave (str);
      MemFree (str);
    } else {
      gbq->val = StringSave ("?");
    }

    /* do not allow saving of bad qualifier */
    if (convertBadToNote &&
        ValidateInferenceQualifier (gbq->val, FALSE) != VALID_INFERENCE) {
      if (StringNICmp (gbq->val, "similar to ", 11) == 0) {
        len = StringLen ("similar to ") + StringLen (first) + StringLen (second);
        str = MemNew (len + 5);
        if (str != NULL) {
          StringCpy (str, "similar to ");
          if (StringDoesHaveText (first)) {
            StringCat (str, first);
            if (StringDoesHaveText (second)) {
              StringCat (str, ":");
              StringCat (str, second);
            }
          } else if (StringDoesHaveText (second)) {
            StringCat (str, second);
          }
          gbq->val = MemFree (gbq->val);
          gbq->val = StringSave (str);
          MemFree (str);
        }
      }
      gbq->qual = MemFree (gbq->qual);
      gbq->qual = StringSave ("note");
    }

    if (sfp->qual == NULL) {
      sfp->qual = gbq;
    }
    if (lastgbq != NULL) {
      lastgbq->next = gbq;
    }
    lastgbq = gbq;
  }
}

static void CleanupInferProc (GraphiC g, VoidPtr data)

{
  InferDialogPtr  idp;
  InferEvidPtr    iep;
  Int2            k;

  idp = (InferDialogPtr) data;
  if (idp != NULL) {
    for (k = 0; k < 128; k++) {
      iep = idp->evidence [k];
      InferEvidFree (iep);
      idp->evidence [k] = NULL;
    }
  }
  StdCleanupExtraProc (g, data);
}

static DialoG NewCreateInferenceDialog (
  GrouP prnt
)

{
  GrouP           cts, tbl, g0, g1, g2, g3, g4, g5, p;
  FonT            fnt;
  Int2            i, hgt, wid;
  InferDialogPtr  idp;

  idp = (InferDialogPtr) MemNew (sizeof (InferDialog));
  if (idp == NULL) return NULL;

  p = HiddenGroup (prnt, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  SetObjectExtra (p, idp, CleanupInferProc);
  idp->dialog = (DialoG) p;
  /*
  idp->todialog = GBQualToInferTable;
  idp->fromdialog = InferTableToGBQual;
  */

  SelectFont (systemFont);
  hgt = LineHeight ();
  inferColFmt [0].pixWidth = MaxAlistWidths (inference_alist) + 5;
  inferColFmt [1].pixWidth = 25 * StringWidth ("X") + 5;
  SelectFont (systemFont);

  wid = 0;
  for (i = 0; i < 2; i++) {
    wid += inferColFmt [i].pixWidth;
  }

  tbl = HiddenGroup (p, -1, 0, NULL);
  SetGroupSpacing (tbl, 10, 5);
  SetGroupMargins (tbl, 5, 5);

  g0 = HiddenGroup (tbl, 15, 0, NULL);
  SetGroupSpacing (g0, 0, 3);
#ifdef WIN_MSWIN
  fnt = systemFont;
#else
  fnt = programFont;
#endif
  /*
  StaticPrompt (g0, "Category", inferColFmt [0].pixWidth, 0, fnt, 'c');
  StaticPrompt (g0, "Explanation", inferColFmt [1].pixWidth, 0, fnt, 'c');
  */

  idp->inferdoc = DocumentPanel (tbl, wid + 2, NUM_INFERENCE_LINES * hgt + 2);
  SetObjectExtra (idp->inferdoc, idp, NULL);
  SetDocCache (idp->inferdoc, NULL, NULL, NULL);
  SetDocNotify (idp->inferdoc, ChangeInferTableSelect);
  idp->numInf = 0;

  for (i = 0; i < NUM_INFERENCE_LINES; i++) {
    AppendItem (idp->inferdoc, PrintInferTable, idp, FALSE, 1,
                &inferParFmt, inferColFmt, systemFont);
  }

  cts = HiddenGroup (p, -1, 0, NULL);
  SetGroupSpacing (cts, 10, 10);
  SetGroupMargins (cts, 5, 5);

  g1 = HiddenGroup (cts, -10, 0, NULL);
  SetGroupSpacing (g1, 5, 5);

  StaticPrompt (g1, "Category", 0, popupMenuHeight, programFont, 'l');
  idp->prefix = CreateEnumPopupDialog (g1, TRUE, ChangeInferPrefix, inference_alist, (UIEnum) 0, idp);

  idp->species = CheckBox (g1, "(same species)", ChangeSameSpecies);
  SetObjectExtra (idp->species, idp, NULL);
  Hide (idp->species);

  g2 = HiddenGroup (cts, 0, 0, NULL);
  SetGroupSpacing (g2, 5, 5);

  g3 = HiddenGroup (g2, -3, 0, NULL);
  SetGroupSpacing (g3, 5, 5);

  StaticPrompt (g3, "Database", 0, dialogTextHeight, programFont, 'l');
  idp->database = CreateEnumPopupDialog (g3, TRUE, ChangeInferDatabase, accn_type_alist, (UIEnum) 0, idp);
  idp->other_db_group = HiddenGroup (g3, -4, 0, NULL);
  StaticPrompt (idp->other_db_group, ":", 0, dialogTextHeight, programFont, 'l');
  idp->db_other = DialogText (idp->other_db_group, "", 8, ChangeInferDbOther);
  SetObjectExtra (idp->db_other, idp, NULL);
  Hide (idp->other_db_group);

  StaticPrompt (g3, "Accession", 0, dialogTextHeight, programFont, 'l');
  idp->accession = DialogText (g3, "", 10, ChangeInferAccession);
  SetObjectExtra (idp->accession, idp, NULL);

  idp->inf_accn_group = g3;
  Hide (idp->inf_accn_group);

  g4 = HiddenGroup (g2, -3, 0, NULL);
  SetGroupSpacing (g4, 5, 5);

  StaticPrompt (g4, "Program", 0, dialogTextHeight, programFont, 'l');
  idp->program = CreateEnumPopupDialog (g4, TRUE, ChangeInferProgram, program_alist, (UIEnum) 0, idp);
  idp->other_pr_group = HiddenGroup (g4, -4, 0, NULL);
  StaticPrompt (idp->other_pr_group, ":", 0, dialogTextHeight, programFont, 'l');
  idp->pr_other = DialogText (idp->other_pr_group, "", 8, ChangeInferPrOther);
  SetObjectExtra (idp->pr_other, idp, NULL);
  Hide (idp->other_pr_group);

  StaticPrompt (g4, "Program Version", 0, dialogTextHeight, programFont, 'l');
  idp->version = DialogText (g4, "", 3, ChangeInferVersion);
  SetObjectExtra (idp->version, idp, NULL);

  idp->inf_prog_group = g4;
  Hide (idp->inf_prog_group);

  g5 = HiddenGroup (g2, 2, 0, NULL);
  SetGroupSpacing (g5, 5, 5);

  StaticPrompt (g5, "Program or Database", 0, dialogTextHeight, programFont, 'l');
  idp->basis1 = DialogText (g5, "", 10, ChangeInferBasis1);
  SetObjectExtra (idp->basis1, idp, NULL);

  StaticPrompt (g5, "Version or Accession", 0, dialogTextHeight, programFont, 'l');
  idp->basis2 = DialogText (g5, "", 10, ChangeInferBasis2);
  SetObjectExtra (idp->basis2, idp, NULL);

  idp->inf_free_group = g5;
  Hide (idp->inf_free_group);

  AlignObjects (ALIGN_CENTER, (HANDLE) tbl, (HANDLE) cts, NULL);

  idp->numInf = 0;
  idp->currItem = 1;
  SetDocHighlight (idp->inferdoc, 1, 1);

  return (DialoG) p;
}

