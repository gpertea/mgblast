/*   dlogutil.h
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
* File Name:  dlogutil.h
*
* Author:  Jonathan Kans
*
* Version Creation Date:   1/22/95
*
* $Revision: 6.51 $
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

#ifndef _DLOGUTIL_
#define _DLOGUTIL_

#include <vibrant.h>
#include <sqnutils.h>
#include <prtutil.h>
#include <asn.h>
#include <objall.h>
#include <objfeat.h>
#include <objfdef.h>
#include <sequtil.h>


#ifdef __cplusplus
extern "C" {
#endif

#define CHANGE_VIEW_NOTABS        0
#define CHANGE_VIEW_FOLDERTABS    1
#define CHANGE_VIEW_TEXTTABS      2
#define CHANGE_VIEW_RADIOBUTTONS  3
#define CHANGE_VIEW_POPUP         4

extern StdPrintOptionsPtr  spop;

extern Boolean SetupPrintOptions (void);
extern void FreePrintOptions (void);

extern EnumFieldAssoc  months_alist [];

extern void DatePtrToVibrant (DatePtr dp, PopuP dateMonth, TexT dateDay, TexT dateYear);
extern DatePtr VibrantToDatePtr (PopuP dateMonth, TexT dateDay, TexT dateYear);

extern ValNodePtr AddStringToValNodeChain (ValNodePtr head, CharPtr str, Uint1 choice);

#define DESCRIPTOR_FORM_BLOCK       \
  FORM_MESSAGE_BLOCK                \
  DialoG          data;


#define FEATURE_FORM_BLOCK          \
  DESCRIPTOR_FORM_BLOCK             \
  DialoG          commonRadio;      \
  GrouP           commonSubGrp [8]; \
  Int2            commonPage;       \
  ButtoN          partial;          \
  ButtoN          exception;        \
  ButtoN          pseudo;           \
  DialoG          experiment;       \
  DialoG          inference;        \
  TexT            comment;          \
  TexT            title;            \
  TexT            exceptText;       \
  PopuP           evidence;         \
  Handle          gene;             \
  PopuP           genePopup;        \
  LisT            geneList;         \
  ValNodePtr      geneNames;        \
  GrouP           useGeneXref;      \
  GrouP           newGeneGrp;       \
  TexT            geneSymbol;       \
  TexT            geneAllele;       \
  TexT            geneDesc;         \
  TexT            locusTag;         \
  ButtoN          editGeneBtn;      \
  TexT            protXrefName;     \
  DialoG          location;         \
  DialoG          product;          \
  DialoG          featcits;         \
  DialoG          dbxrefs;          \
  DialoG          gbquals;          \
  DialoG          usrobjext;        \
  TexT            featid;           \
  TexT            fidxref;          \
  ButtoN          leave_dlg_up;     \
  UserObjectPtr   goTermUserObj;    \
  Int2            badInfAttempts;   \
  Boolean         acceptBadInf;

typedef struct descform {
  DESCRIPTOR_FORM_BLOCK
} DescriptorForm, PNTR DescriptorFormPtr;

typedef struct featform {
  FEATURE_FORM_BLOCK
} FeatureForm, PNTR FeatureFormPtr;

extern void SetDescriptorPropagate (BioseqSetPtr bssp);
extern Int2 LIBCALLBACK DescriptorPropagate (Pointer data);
extern Boolean DescFormReplaceWithoutUpdateProc (ForM f);
extern void StdDescFormActnProc (ForM f);

extern void StdDescFormCleanupProc (GraphiC g, VoidPtr data);

extern Boolean FeatFormReplaceWithoutUpdateProc (ForM f);
extern void StdFeatFormActnProc (ForM f);

extern void StdSeqFeatPtrToFeatFormProc (ForM f, Pointer data);
extern void StdInitFeatFormProc (ForM f);
extern void StdFeatFormCleanupProc (GraphiC g, VoidPtr data);
extern void StdFeatFormAcceptButtonProc (ButtoN b);

extern void InferenceDialogToGBQuals (DialoG d, SeqFeatPtr sfp, Boolean convertBadToNote);
extern void GBQualsToInferenceDialog (DialoG d, SeqFeatPtr sfp);

/*
extern void ExtendGeneFeatIfOnMRNA (Uint2 entityID, SeqEntryPtr sep);
*/

extern OMUserDataPtr ItemAlreadyHasEditor (Uint2 entityID, Uint2 itemID, Uint2 itemtype, Uint2 procid);
extern Int2 LIBCALLBACK StdVibrantEditorMsgFunc (OMMsgStructPtr ommsp);

extern Boolean TestInference (FeatureFormPtr ffp, CharPtr badInfQual, size_t len, CharPtr badInfMssg);

/*****************************************************************************
*
*   Various specific spreadsheet or editor dialogs.
*
*****************************************************************************/

extern DialoG CreateAuthorDialog (GrouP prnt, Uint2 rows, Int2 spacing);

extern DialoG CreateDateDialog (GrouP prnt, CharPtr title);

extern DialoG CreateAffilDialog (GrouP prnt, CharPtr title);
extern DialoG CreatePublisherAffilDialog (GrouP prnt, CharPtr title);
extern DialoG CreateProceedingsDialog (GrouP prnt, CharPtr title);

/* The Ext versions split the affil dialogs into two superimposed pages,
   both of which are initially hidden, that can go in separate folder tabs */

extern DialoG CreateExtAffilDialog (GrouP prnt, CharPtr title,
                                    GrouP PNTR grp1, GrouP PNTR grp2);
extern DialoG CreateExtPublisherAffilDialog (GrouP prnt, CharPtr title,
                                             GrouP PNTR grp1, GrouP PNTR grp2);
extern DialoG CreateExtProceedingsDialog (GrouP prnt, CharPtr title,
                                          GrouP PNTR grp1, GrouP PNTR grp2);

extern GrouP CreateCommonFeatureGroup (GrouP h, FeatureFormPtr ffp,
                                       SeqFeatPtr sfp, Boolean hasGeneControl,
                                       Boolean hasCitationTab);
extern GrouP CreateCommonFeatureGroupEx (GrouP h, FeatureFormPtr ffp,
                                         SeqFeatPtr sfp, Boolean hasGeneControl,
                                         Boolean hasCitationTab, Boolean hasGeneSuppress);
extern void UpdateGeneLocation (SeqFeatPtr gene, SeqLocPtr old_feat_loc,
                                SeqLocPtr new_feat_loc, Uint2 entityID);
extern void PopulateGenePopup (FeatureFormPtr ffp);
extern void SeqFeatPtrToCommon (FeatureFormPtr ffp, SeqFeatPtr sfp);
extern void SetNewFeatureDefaultInterval (FeatureFormPtr ffp);
extern DialoG CreateImportFields (GrouP h, CharPtr name, SeqFeatPtr sfp, Boolean allowProductGBQual);

extern DialoG CreateIntervalEditorDialog (GrouP h, CharPtr title, Uint2 rows,
                                          Int2 spacing, SeqEntryPtr sep,
                                          Boolean nucsOK, Boolean protsOK);

typedef void (*IntEdPartialProc) (FeatureFormPtr ffp, Boolean partial5, Boolean partial3, Boolean order);

extern DialoG CreateIntervalEditorDialogEx (GrouP h, CharPtr title, Uint2 rows,
                                            Int2 spacing, SeqEntryPtr sep,
                                            Boolean nucsOK, Boolean protsOK,
                                            Boolean useBar, Boolean showPartials,
                                            Boolean allowGaps, FeatureFormPtr ffp,
                                            IntEdPartialProc proc);
                                            
extern DialoG CreateIntervalEditorDialogExEx (GrouP h, CharPtr title, Uint2 rows,
                                              Int2 spacing, SeqEntryPtr sep,
                                              Boolean nucsOK, Boolean protsOK,
                                              Boolean useBar, Boolean showPartials,
                                              Boolean allowGaps, FeatureFormPtr ffp,
                                              IntEdPartialProc proc, 
                                              Boolean use_aln, Boolean show_seqid,
                                              TaglistCallback tlp_callback, 
                                              Pointer callback_data,
                                              Boolean allow_nulls_in_list);

extern void DisplayErrorMessages (CharPtr title, ValNodePtr err_list);

extern Boolean AdjustFromForGap (Int4Ptr p_from, SeqAlignPtr salp, Int4 aln_len, Int4 aln_row);
extern Boolean AdjustToForGap (Int4Ptr p_to, SeqAlignPtr salp, Int4 aln_row);
                                            
extern void SetSequenceAndStrandForIntervalPage (DialoG d);                                            
extern void StdFeatIntEdPartialCallback (FeatureFormPtr ffp, Boolean partial5, Boolean partial3, Boolean order);

extern DialoG CreateVisibleStringDialog (GrouP h, Uint2 rows,
                                         Int2 spacing, Int2 width);

extern DialoG CreateDbtagDialog (GrouP h, Uint2 rows, Int2 spacing,
                                 Int2 width1, Int2 width2);

/*****************************************************************************
*
*   Miscellaneous functions.
*
*****************************************************************************/

/*
*  The TextViewProcsPtr may be registered with a call to SetAppProperty
*  e.g., SetAppProperty ("TextDisplayForm", &viewprocs), where viewprocs
*  is a persistent structure filled with callback function pointers specific
*  for a given application.
*/

typedef struct textviewprocs {
  WndActnProc      activateForm;
  WndActnProc      closeForm;
  WndActnProc      createMenus;

  FormMessageFunc  handleMessages;

  Int2             minPixelWidth;
  Int2             minPixelHeight;
  Boolean          useScrollText;

  FonT             displayFont;

  Char             screenMode;

} TextViewProcs, PNTR TextViewProcsPtr;

typedef CharPtr (*Nlm_RepopulateViewer) PROTO ((Pointer userdata));
typedef Pointer (*Nlm_RepopulateDataFree) PROTO ((Pointer userdata));

extern void LaunchGeneralTextViewer (CharPtr path, CharPtr title);
extern void LaunchAsnTextViewer (Pointer from, AsnWriteFunc writefunc, CharPtr title);
extern void LaunchGeneralTextViewerWithRepopulate 
         (CharPtr path, CharPtr title, 
         Nlm_RepopulateViewer repopulate_func, Pointer repopulate_data,
         Nlm_RepopulateDataFree free_data_func);

extern Boolean FileToScrollText (TexT t, CharPtr path);
extern void ScrollTextToFile (TexT t, CharPtr path);

extern void FileToClipboard (CharPtr path);

/*
*  The StdEditorProcsPtr may be registered with a call to SetAppProperty
*  e.g., SetAppProperty ("StdEditorForm", &viewprocs), where viewprocs
*  is a persistent structure filled with callback function pointers specific
*  for a given application.
*/

typedef struct stdeditorprocs {
  WndActnProc      activateForm;
  FormMessageFunc  handleMessages;

  Boolean          duplicateExisting;

  Char             screenMode;

} StdEditorProcs, PNTR StdEditorProcsPtr;

#ifndef WIN_MAC
extern void CreateStdEditorFormMenus (WindoW w);
#endif

/*
*  When duplicating instead of editing, the editors call the following
*  function to change the baseFormPtr itemID and itemtype to the closest
*  Bioseq or BioseqSet, thus adding the new item to that parent.
*/

extern Boolean SetClosestParentIfDuplicating (BaseFormPtr bfp);

/*
*  The HelpMessageFunc may be registered with a call to SetAppProperty
*  e.g., SetAppProperty ("HelpMessageProc", helpproc), where helpproc
*  is an application-supplied callback function that takes character
*  strings for the heading and section.
*/

typedef void (*HelpMessageFunc) (CharPtr, CharPtr);

/*****************************************************************************
*
*   Miscellaneous functions.
*
*****************************************************************************/

extern CharPtr SaveStringFromTextAndStripNewlines (TexT t);

/* GetRidOfEmptyFeatsDescStrings takes either an entityID or a SeqEntryPtr */

extern void GetRidOfEmptyFeatsDescStrings (Uint2 entityID, SeqEntryPtr sep);

/*****************************************************************************
*
*   LaunchEntrezURL constructs a web-Entrez URL query.  Parameter choices are:
*
*     Database      Report Formats
*
*      PubMed        DocSum, Brief, Abstract, Citation, MEDLINE, ASN.1, ExternalLink
*      Nucleotide    DocSum, Brief, GenBank, ASN.1, FASTA, ExternalLink
*      Protein       DocSum, Brief, GenPept, ASN.1, FASTA, ExternalLink
*      Genome        DocSum, Brief, ASN.1, ExternalLink
*      Popset        DocSum, Brief, ASN.1, FASTA, ExternalLink
*      Structure     DocSum, Brief
*      OMIM          DocSum, Text, Synopsis, Variants, MiniMIM, ASN.1, ExternalLink
*      Taxonomy      DocSum, Brief, TxInfo, TxTree, ExternalLink
*
*****************************************************************************/

NLM_EXTERN void LaunchEntrezURL (CharPtr database, Int4 uid, CharPtr format);

/* this function creates a dialog for selecting publications in a record */
extern DialoG FeatCitEditDialog (GrouP parent, Uint2 entityID);

/* This structure and the functions after it are used for
 * modal dialogs with two choices.
 */
typedef struct modalacceptcancel
{
  Boolean accepted;
  Boolean cancelled; 
  Boolean third_option; 
} ModalAcceptCancelData, PNTR ModalAcceptCancelPtr;

extern void ModalAcceptButton (ButtoN b);
extern void ModalCancelButton (ButtoN b);
extern void ModalThirdOptionButton (ButtoN b);

typedef void (*TableDisplayDblClick) PROTO((PoinT, CharPtr, CharPtr, Pointer));
typedef Boolean (*TableDisplayLeftInRed) PROTO ((Int4, ValNodePtr, Pointer));
extern DialoG TableDisplayDialog (GrouP parent, Int4 width, Int4 height,
                                  Int4 frozen_header, Int4 frozen_left,
                                  TableDisplayDblClick dbl_click,
                                  Pointer dbl_click_data,
                                  TableDisplayLeftInRed,
                                  Pointer left_in_red_data);
extern ValNodePtr FreeTableDisplayRowList (ValNodePtr row_list);
extern void PrintTableDisplayRowListToFile (ValNodePtr row_list, FILE *fp);
extern ValNodePtr CopyTableDisplayRowList (ValNodePtr row_list);
extern CharPtr GetRowListCellText (ValNodePtr row_list, Int4 row, Int4 column);
extern FonT GetTableDisplayDefaultFont (void);

typedef  void  (*Nlm_ChangeNotifyProc) PROTO ((Pointer));

#define TALL_SELECTION_LIST 8
#define SHORT_SELECTION_LIST 4
/* err_msg is the message to put in the results from TestDialog if nothing is selected */
/* choice_list should be a valnode list of strings to use for the names of the choices. */
/* All is automatically included as a choice if allow_multi is true. */
/* The ValNodeList returned is a list of integers indicating the position of the item
 * in the list: 1 is the first item, 2 is the second item, etc. */
extern DialoG SelectionDialog 
(GrouP h,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Boolean                  allow_multi,
 CharPtr                  err_msg,
 ValNodePtr               choice_list,
 Int2                     list_height);

extern DialoG SelectionDialogEx 
(GrouP h,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Boolean                  allow_multi,
 CharPtr                  err_msg,
 ValNodePtr               choice_list,
 Int2                     list_height,
 Boolean                  force_list);

/* This function should free just the data associated with the ValNode */
typedef void (*FreeValNodeProc) PROTO ((ValNodePtr));
/* This function should copy just this ValNode, not the chain */
typedef ValNodePtr (*CopyValNodeDataProc) PROTO ((ValNodePtr));
typedef ValNodePtr (*RemapValNodeProc) PROTO ((ValNodePtr));
typedef Boolean (*MatchValNodeProc) PROTO ((ValNodePtr, ValNodePtr));
typedef CharPtr (*NameFromValNodeProc) PROTO ((ValNodePtr));

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
 RemapValNodeProc         remap_vn_proc);
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
 Boolean                  allow_multi);

extern DialoG EnumAssocSelectionDialog 
(GrouP                 h,
 Nlm_EnumFieldAssocPtr eap,
 CharPtr               err_name,
 Boolean               allow_multi,
 Nlm_ChangeNotifyProc  change_notify,
 Pointer               change_userdata);

extern CharPtr ValNodeStringName (ValNodePtr vnp);
extern void ValNodeSimpleDataFree (ValNodePtr vnp);
extern ValNodePtr ValNodeStringCopy (ValNodePtr vnp);
extern Boolean ValNodeChoiceMatch (ValNodePtr vnp1, ValNodePtr vnp2);
extern Boolean ValNodeStringMatch (ValNodePtr vnp1, ValNodePtr vnp2);

extern DialoG SequenceSelectionDialog 
(GrouP h,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Boolean                  allow_multi,
 Boolean                  show_nucs,
 Boolean                  show_prots,
 Uint2                    entityID);


#ifdef __cplusplus
}
#endif

#endif /* ndef _DLOGUTIL_ */

