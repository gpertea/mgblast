/*   seqsub.c
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
* File Name:  seqsub.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   1/22/95
*
* $Revision: 6.33 $
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

#include <seqsub.h>
#include <gather.h>
#include <biosrc.h>
#include <pubdesc.h>
#include <asn2gnbp.h>
#include <document.h>
#include <explore.h>
#include <subutil.h>

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

typedef struct contactnamedialog 
{
  DIALOG_MESSAGE_BLOCK
  TexT            firstname;
  TexT            middleinit;
  TexT            lastname;
  PopuP           suffix;  
} ContactNameDialogData, PNTR ContactNameDialogPtr;

extern CharPtr NameStdPtrToAuthorSpreadsheetString (NameStdPtr nsp);

static void AuthorToContactNameDialog (DialoG d, Pointer userdata)
{
  ContactNameDialogPtr dlg;
  AuthorPtr            ap;
  PersonIdPtr          pid;
  NameStdPtr           nsp;
  CharPtr              str;
  CharPtr              txt;
  
  dlg = (ContactNameDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  ap = (AuthorPtr) userdata;
  if (ap == NULL)
  {
    SafeSetTitle (dlg->firstname, "");
    SafeSetTitle (dlg->middleinit, "");  
    SafeSetTitle (dlg->lastname, "");
    SafeSetValue (dlg->suffix, 0);  
  }
  else
  {
    pid = ap->name;
    if (pid != NULL && pid->choice == 2) {
      nsp = pid->data;
      if (nsp != NULL) {
        str = NameStdPtrToAuthorSpreadsheetString (nsp);
        if (str != NULL) {
          txt = ExtractTagListColumn (str, 0);
          SafeSetTitle (dlg->firstname, txt);
          MemFree (txt);
          txt = ExtractTagListColumn (str, 1);
          SafeSetTitle (dlg->middleinit, txt);
          MemFree (txt);
          txt = ExtractTagListColumn (str, 2);
          SafeSetTitle (dlg->lastname, txt);
          MemFree (txt);
          txt = ExtractTagListColumn (str, 3);
          SafeSetValue (dlg->suffix, atoi(txt)+1);
          MemFree (str);
        }
      }
    }
  }
}

extern NameStdPtr AuthorSpreadsheetStringToNameStdPtr (CharPtr txt);

static Pointer ContactNameDialogToAuthor (DialoG d)
{
  ContactNameDialogPtr dlg;
  NameStdPtr           nsp;
  AuthorPtr            ap;
  PersonIdPtr          pid;
  Char                 str [128];
  CharPtr              txt;
  Char                 suffix [32];
  Uint2                suffixVal;
  
  dlg = (ContactNameDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  
  nsp = NULL;
  ap = AuthorNew ();
  if (ap != NULL) {
    pid = PersonIdNew ();
    ap->name = pid;
    if (pid != NULL) {
      pid->choice = 2;
      str [0] = '\0';
      txt = SaveStringFromText (dlg->firstname);
      StringCat (str, txt);
      StringCat (str, "\t");
      MemFree (txt);
      txt = SaveStringFromText (dlg->middleinit);
      StringCat (str, txt);
      StringCat (str, "\t");
      MemFree (txt);
      txt = SaveStringFromText (dlg->lastname);
      StringCat (str, txt);
      StringCat (str, "\t");
      MemFree (txt);
      suffixVal = GetValue (dlg->suffix);
      if (suffixVal < 1)
      {
        suffixVal = 1;
      }
      sprintf (suffix, "%d", (int) (suffixVal - 1));
      StringCat (str, suffix);
      StringCat (str, "\n");
      txt = StringSave (str);
      nsp = AuthorSpreadsheetStringToNameStdPtr (txt);
      MemFree (txt);
      pid->data = nsp;
    } 
  }
  return (Pointer) ap;
}

static ValNodePtr TestContactNameDialog (DialoG d)
{
  ContactNameDialogPtr dlg;
  ValNodePtr           head;
  
  dlg = (ContactNameDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }


  head = NULL;
  dlg = (ContactNameDialogPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    if (TextHasNoText (dlg->firstname)) {
      head = AddStringToValNodeChain (head, "You must specify a first name for the contact!", 0);
    }
    if (TextHasNoText (dlg->lastname)) {
      head = AddStringToValNodeChain (head, "You must specify a last name for the contact!", 0);
    }
  }
  return head;
  
}

static DialoG ContactNameDialog (GrouP parent)
{
  ContactNameDialogPtr dlg;
  GrouP                g;
  
  dlg = (ContactNameDialogPtr) MemNew (sizeof (ContactNameDialogData));
  if (dlg == NULL)
  {
    return NULL;
  }
  
  g = HiddenGroup (parent, 4, 0, NULL);
  SetGroupSpacing (g, -1, 2);
  
  SetObjectExtra (g, dlg, StdCleanupExtraProc);
  dlg->dialog = (DialoG) g;
  dlg->todialog = AuthorToContactNameDialog;
  dlg->fromdialog = ContactNameDialogToAuthor;
  dlg->testdialog = TestContactNameDialog;
  
  StaticPrompt (g, "First Name", 0, 0, programFont, 'c');
  StaticPrompt (g, "M.I.", 0, 0, programFont, 'c');
  StaticPrompt (g, "Last Name", 0, 0, programFont, 'c');
  StaticPrompt (g, "Sfx", 0, 0, programFont, 'c');
  dlg->firstname = DialogText (g, "", 8, NULL);
  dlg->middleinit = DialogText (g, "", 4, NULL);
  dlg->lastname = DialogText (g, "", 9, NULL);
  dlg->suffix = PopupList (g, TRUE, NULL);
  InitEnumPopup (dlg->suffix, name_suffix_alist, NULL);
  SetEnumPopup (dlg->suffix, name_suffix_alist, 0);

  return (DialoG) g;
}

typedef struct contactpage {
  DIALOG_MESSAGE_BLOCK
  DialoG          contact_name_dlg;
  TexT            firstname;
  TexT            middleinit;
  TexT            lastname;
  PopuP           suffix;
  DialoG          affil;
  GrouP           contactGrp [3];
  DialoG          citsub_dlg;
  ButtoN          copy_btn;
  Int4            currentPage;
  DialoG          tbs;           /* folder tabs dialog */
} ContactPage, PNTR ContactPagePtr;

static void ContactInfoPtrToContactPage (DialoG d, Pointer data)

{
  ContactInfoPtr  cip;
  ContactPagePtr  cpp;

  cpp = (ContactPagePtr) GetObjectExtra (d);
  cip = (ContactInfoPtr) data;
  if (cpp != NULL) {
    PointerToDialog (cpp->contact_name_dlg, NULL);
    PointerToDialog (cpp->affil, NULL);
    if (cip != NULL && cip->contact != NULL) {
      PointerToDialog (cpp->contact_name_dlg, cip->contact);  
      PointerToDialog (cpp->affil, cip->contact->affil);
    }
  }
}

static Pointer ContactPageToContactInfoPtr (DialoG d)

{
  ContactInfoPtr  cip;
  ContactPagePtr  cpp;


  cip = NULL;
  cpp = (ContactPagePtr) GetObjectExtra (d);
  if (cpp != NULL) {
    cip = ContactInfoNew ();
    if (cip != NULL) {
      cip->contact = DialogToPointer (cpp->contact_name_dlg);
      if (cip->contact == NULL)
      {
        cip->contact = AuthorNew ();
      }
      if (cip->contact != NULL)
      {
        cip->contact->affil = DialogToPointer (cpp->affil);
      }
      if (cip->contact->affil == NULL
          || cip->contact->name == NULL
          || cip->contact->name->data == NULL)
      {
        cip = ContactInfoFree (cip);
      }
    }
  }
  return (Pointer) cip;
}

static Boolean ReadContactDialog (DialoG d, CharPtr filename)

{
  AsnIoPtr        aip;
  ContactInfoPtr  cip;
  ContactPagePtr  cpp;
  Char            path [PATH_MAX];

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  cpp = (ContactPagePtr) GetObjectExtra (d);
  if (cpp != NULL) {
    if (path [0] != '\0' || GetInputFileName (path, sizeof (path), "", "TEXT")) {
      aip = AsnIoOpen (path, "r");
      if (aip != NULL) {
        cip = ContactInfoAsnRead (aip, NULL);
        AsnIoClose (aip);
        if (cip != NULL) {
          ContactInfoPtrToContactPage (cpp->dialog, cip);
          cip = ContactInfoFree (cip);
          Update ();
          Select (cpp->firstname);
          return TRUE;
        }
      }
    }
  }
  return FALSE;
}

static Boolean WriteContactDialog (DialoG d, CharPtr filename)

{
  AsnIoPtr        aip;
  ContactInfoPtr  cip;
  ContactPagePtr  cpp;
  Char            path [PATH_MAX];
#ifdef WIN_MAC
  FILE            *f;
#endif

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  cpp = (ContactPagePtr) GetObjectExtra (d);
  if (cpp != NULL) {
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
        cip = ContactPageToContactInfoPtr (cpp->dialog);
        ContactInfoAsnWrite (cip, aip, NULL);
        AsnIoClose (aip);
        cip = ContactInfoFree (cip);
        return TRUE;
      }
    }
  }
  return FALSE;
}

static ValNodePtr TestContactDialog (DialoG d)

{
  ContactPagePtr  cpp;
  ValNodePtr      head;
  AffilPtr        affil;

  head = NULL;
  cpp = (ContactPagePtr) GetObjectExtra (d);
  if (cpp != NULL) {
    head = TestDialog (cpp->contact_name_dlg);
  }
  affil = (AffilPtr) DialogToPointer (cpp->affil);
  if (affil == NULL)
  {
    ValNodeAddPointer (&head, 1, StringSave ("You must supply an affiliation for the contact"));
  }
  
  return head;
}

static CharPtr contactTabs [] = {
  "Name", "Affiliation", "Contact", NULL
};

static void ChangeContactSubPage (VoidPtr data, Int2 newval, Int2 oldval)

{
  ContactPagePtr  cpp;

  cpp = (ContactPagePtr) data;
  if (cpp != NULL) {
    cpp->currentPage = newval;
    if (oldval >= 0 && oldval <= 2) {
      SafeHide (cpp->contactGrp [oldval]);
    }
    if (newval >= 0 && newval <= 2) {
      SafeShow (cpp->contactGrp [newval]);
    }
    
    if (cpp->currentPage == 0)
    {
      SafeSetTitle (cpp->copy_btn, "Copy Contact Name to Citation");
      SafeShow (cpp->copy_btn);
    }
    else if (cpp->currentPage == 1)
    {
      SafeSetTitle (cpp->copy_btn, "Copy Affiliation to Citation");
      SafeShow (cpp->copy_btn);
    }
    else
    {
      SafeHide (cpp->copy_btn);
    }
    
    Update ();
  }
}

static void CopyToCitation (ButtoN b)
{
  ContactPagePtr  cpp;
  CitSubPtr       csp;
  AuthorPtr       ap;

  cpp = (ContactPagePtr) GetObjectExtra (b);
  if (cpp == NULL || cpp->citsub_dlg == NULL)
  {
    return;
  }
  if (cpp->currentPage != 0 && cpp->currentPage != 1)
  {
    return;
  }
  
  csp = (CitSubPtr) DialogToPointer (cpp->citsub_dlg);
  if (csp == NULL)
  {
    csp = CitSubNew ();
    if (csp == NULL)
    {
      return;
    }
  }
  
  if (csp->authors == NULL)
  {
    csp->authors = AuthListNew ();
    if (csp->authors == NULL)
    {
      return;
    }
    else
    {
      csp->authors->choice = 1;
    }
  }
  
  if (cpp->currentPage == 0)
  {
    /* copy name */
    ap = (AuthorPtr) DialogToPointer (cpp->contact_name_dlg);
    if (ap != NULL)
    {
      ValNodeAddPointer (&(csp->authors->names), 1, ap);
    }
  }
  else if (cpp->currentPage == 1)
  {
    /* copy affiliation */
    csp->authors->affil = AffilFree (csp->authors->affil);
    csp->authors->affil = DialogToPointer (cpp->affil);
  }
  
  PointerToDialog (cpp->citsub_dlg, csp);
  csp = CitSubFree (csp);
}

static void ContactDialogMessage (DialoG d, Int2 mssg)

{
  ContactPagePtr dlg;
  Int2           pageval;

  dlg = (ContactPagePtr) GetObjectExtra (d);
  if (dlg != NULL) {
    if (mssg == VIB_MSG_ENTER)
    {
      SetValue (dlg->tbs, 0);
    }
    else if (mssg > NUM_VIB_MSG)
    {
      pageval = mssg - NUM_VIB_MSG - 1;
      if (pageval < 3)
      {
        SetValue (dlg->tbs, pageval);
      }
    }
  }
}


extern DialoG CreateContactDialog (GrouP h, CharPtr title, DialoG citsub_dlg)

{
  ContactPagePtr  cpp;
  GrouP           k;
  GrouP           m;
  GrouP           p;
  GrouP           s;

  p = HiddenGroup (h, 1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  cpp = (ContactPagePtr) MemNew (sizeof (ContactPage));
  if (cpp != NULL) {

    SetObjectExtra (p, cpp, StdCleanupExtraProc);
    cpp->dialog = (DialoG) p;
    cpp->todialog = ContactInfoPtrToContactPage;
    cpp->fromdialog = ContactPageToContactInfoPtr;
    cpp->testdialog = TestContactDialog;
    cpp->importdialog = ReadContactDialog;
    cpp->exportdialog = WriteContactDialog;
    cpp->dialogmessage = ContactDialogMessage;
    
    cpp->citsub_dlg = citsub_dlg;
    cpp->currentPage = 0;

    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    m = HiddenGroup (s, -1, 0, NULL);
    SetGroupSpacing (m, 10, 10);

    cpp->tbs = CreateFolderTabs (m, contactTabs, 0, 0, 0,
                            PROGRAM_FOLDER_TAB,
                            ChangeContactSubPage, (Pointer) cpp);
    k = HiddenGroup (m, 0, 0, NULL);

    cpp->contactGrp [0] = HiddenGroup (k, -1, 0, NULL);
    cpp->contact_name_dlg = ContactNameDialog (cpp->contactGrp [0]);

    cpp->affil = CreateExtAffilDialog (k, NULL,
                                       &(cpp->contactGrp [1]),
                                       &(cpp->contactGrp [2]));

    if (cpp->citsub_dlg != NULL)
    {
      cpp->copy_btn = PushButton (m, "Copy Contact Name to Citation", CopyToCitation);
      SetObjectExtra (cpp->copy_btn, cpp, NULL);
    }
    else
    {
      cpp->copy_btn = NULL;
    }

    AlignObjects (ALIGN_CENTER, (HANDLE) cpp->tbs,
                                (HANDLE) cpp->contact_name_dlg, 
                                (HANDLE) cpp->affil, 
                                (HANDLE) cpp->copy_btn, 
                                NULL);
  }

  return (DialoG) p;
}

typedef struct citsubpage {
  DIALOG_MESSAGE_BLOCK
  DialoG          authors;
  TexT            consortium;
  DialoG          date;
  DialoG          affil;
  TexT            descr;
  GrouP           citsubGrp [5];
  DialoG          tbs;           /* folder tabs dialog */
} CitsubPage, PNTR CitsubPagePtr;

static AuthListPtr AddConsortiumToAuthList (AuthListPtr alp, TexT consortium)

{
  AuthorPtr    ap;
  ValNodePtr   names;
  PersonIdPtr  pid;

  if (TextHasNoText (consortium)) return alp;
  if (alp == NULL) {
    alp = AuthListNew ();
    alp->choice = 1;
  }
  pid = PersonIdNew ();
  if (pid == NULL) return NULL;
  pid->choice = 5;
  pid->data = SaveStringFromText (consortium);
  ap = AuthorNew ();
  if (ap == NULL) return NULL;
  ap->name = pid;
  names = ValNodeAdd (&(alp->names));
  names->choice = 1;
  names->data.ptrvalue = ap;
  return alp;
}

static void AuthListToConsortium (AuthListPtr alp, TexT consortium)

{
  AuthorPtr    ap;
  ValNodePtr   names;
  PersonIdPtr  pid;
  CharPtr      str;

  if (alp == NULL || consortium == NULL) return;
  if (alp->choice != 1) return;
  for (names = alp->names; names != NULL; names = names->next) {
    ap = names->data.ptrvalue;
    if (ap == NULL) continue;
    pid = ap->name;
    if (pid == NULL || pid->choice != 5) continue;
    str = (CharPtr) pid->data;
    SafeSetTitle (consortium, str);
  }
}

static void CitSubPtrToCitsubPage (DialoG d, Pointer data)

{
  AuthListPtr    alp;
  CitsubPagePtr  cpp;
  CitSubPtr      csp;
  WindoW         tempPort;

  cpp = (CitsubPagePtr) GetObjectExtra (d);
  csp = (CitSubPtr) data;
  if (cpp != NULL) {
    PointerToDialog (cpp->authors, NULL);
    PointerToDialog (cpp->date, NULL);
    PointerToDialog (cpp->affil, NULL);
    SafeSetTitle (cpp->descr, "");
    tempPort = SavePort (cpp->authors);
    if (csp != NULL) {
      alp = csp->authors;
      PointerToDialog (cpp->authors, (Pointer) alp);
      AuthListToConsortium (alp, cpp->consortium);
      if (alp != NULL) {
        PointerToDialog (cpp->affil, alp->affil);
      }
      if (csp->date != NULL) {
        PointerToDialog (cpp->date, csp->date);
      } else if (csp->imp != NULL) {
        PointerToDialog (cpp->date, csp->imp->date);
      }
      SafeSetTitle (cpp->descr, csp->descr);
    }
    RestorePort (tempPort);
  }
}

static Pointer CitsubPageToCitSubPtr (DialoG d)

{
  AuthListPtr    alp;
  CitsubPagePtr  cpp;
  CitSubPtr      csp;

  csp = NULL;
  cpp = (CitsubPagePtr) GetObjectExtra (d);
  if (cpp != NULL) {
    csp = CitSubNew ();
    if (csp != NULL) {
      alp = (AuthListPtr) DialogToPointer (cpp->authors);
      alp = AddConsortiumToAuthList (alp, cpp->consortium);
      if (alp != NULL) {
        alp->affil = DialogToPointer (cpp->affil);
      }
      csp->authors = alp;
      csp->date = DialogToPointer (cpp->date);
      csp->descr = SaveStringFromTextAndStripNewlines (cpp->descr);
      if (alp == NULL) {
        csp = CitSubFree (csp);
      }
    }
  }
  return (Pointer) csp;
}

static Boolean ReadCitsubDialog (DialoG d, CharPtr filename)

{
  AsnIoPtr       aip;
  CitsubPagePtr  cpp;
  CitSubPtr      csp;
  Char           path [PATH_MAX];

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  cpp = (CitsubPagePtr) GetObjectExtra (d);
  if (cpp != NULL) {
    if (path [0] != '\0' || GetInputFileName (path, sizeof (path), "", "TEXT")) {
      aip = AsnIoOpen (path, "r");
      if (aip != NULL) {
        csp = CitSubAsnRead (aip, NULL);
        AsnIoClose (aip);
        if (csp != NULL) {
          CitSubPtrToCitsubPage (cpp->dialog, csp);
          csp = CitSubFree (csp);
          Update ();
          return TRUE;
        }
      }
    }
  }
  return FALSE;
}

static Boolean WriteCitsubDialog (DialoG d, CharPtr filename)

{
  AsnIoPtr       aip;
  CitsubPagePtr  cpp;
  CitSubPtr      csp;
  Char           path [PATH_MAX];
#ifdef WIN_MAC
  FILE            *f;
#endif

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  cpp = (CitsubPagePtr) GetObjectExtra (d);
  if (cpp != NULL) {
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
        csp = CitsubPageToCitSubPtr (cpp->dialog);
        CitSubAsnWrite (csp, aip, NULL);
        AsnIoClose (aip);
        csp = CitSubFree (csp);
        return TRUE;
      }
    }
  }
  return FALSE;
}

static ValNodePtr TestCitsubDialog (DialoG d)

{
  CitsubPagePtr  cpp;
  ValNodePtr     head;
  AuthListPtr    alp = NULL;
  AffilPtr       affil = NULL;

  head = NULL;
  cpp = (CitsubPagePtr) GetObjectExtra (d);
  if (cpp != NULL) {
    alp = (AuthListPtr) DialogToPointer (cpp->authors);
    alp = AddConsortiumToAuthList (alp, cpp->consortium);
    if (alp == NULL)
    {
      ValNodeAddPointer (&head, 0, StringSave ("You must supply authors for the citation!"));
    }
    alp = AuthListFree (alp);
    
    affil = (AffilPtr) DialogToPointer (cpp->affil);
    if (affil == NULL)
    {
      ValNodeAddPointer (&head, 1, StringSave ("You must supply an affiliation for the citation!"));
    }
    affil = AffilFree (affil);
  }
  return head;
}

static CharPtr citsubTabs [] = {
  "Names", "Affiliation", "Contact", "Description", NULL, NULL
};

static void ChangeCitsubSubPage (VoidPtr data, Int2 newval, Int2 oldval)

{
  CitsubPagePtr  cpp;

  cpp = (CitsubPagePtr) data;
  if (cpp != NULL) {
    if (oldval >= 0 && oldval <= 4) {
      SafeHide (cpp->citsubGrp [oldval]);
    }
    if (newval >= 0 && newval <= 4) {
      SafeShow (cpp->citsubGrp [newval]);
    }
    Update ();
  }
}

static void CitSubDialogMessage (DialoG d, Int2 mssg)
{
  CitsubPagePtr dlg;
  Int2          pageval;

  dlg = (CitsubPagePtr) GetObjectExtra (d);
  if (dlg != NULL) {
    if (mssg == VIB_MSG_ENTER)
    {
      SetValue (dlg->tbs, 0);
    }
    else if (mssg > NUM_VIB_MSG)
    {
      pageval = mssg - NUM_VIB_MSG - 1;
      if (pageval < 3)
      {
        SetValue (dlg->tbs, pageval);
      }
    }
  }
}

extern DialoG CreateCitSubDialog (GrouP h, CharPtr title, CitSubPtr csp)

{
  CitsubPagePtr  cpp;
  GrouP          k;
  GrouP          m;
  GrouP          p;
  GrouP          q;
  GrouP          s;

  p = HiddenGroup (h, 1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  cpp = (CitsubPagePtr) MemNew (sizeof (CitsubPage));
  if (cpp != NULL) {

    SetObjectExtra (p, cpp, StdCleanupExtraProc);
    cpp->dialog = (DialoG) p;
    cpp->todialog = CitSubPtrToCitsubPage;
    cpp->fromdialog = CitsubPageToCitSubPtr;
    cpp->testdialog = TestCitsubDialog;
    cpp->importdialog = ReadCitsubDialog;
    cpp->exportdialog = WriteCitsubDialog;
    cpp->dialogmessage = CitSubDialogMessage;

    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    m = HiddenGroup (s, -1, 0, NULL);
    SetGroupSpacing (m, 10, 10);

    citsubTabs [4] = NULL;
    if (csp != NULL) {
      if (csp->date != NULL) {
        citsubTabs [4] = "Date";
      } else if (csp->imp != NULL && csp->imp->date != NULL) {
        citsubTabs [4] = "Date";
      } else if (GetAppProperty ("InternalNcbiSequin") != NULL) {
        citsubTabs [4] = "Date";
      }
    }
    else if (GetAppProperty ("InternalNcbiSequin") != NULL) {
      citsubTabs [4] = "Date";
    }
    cpp->tbs = CreateFolderTabs (m, citsubTabs, 0, 0, 0,
                            PROGRAM_FOLDER_TAB,
                            ChangeCitsubSubPage, (Pointer) cpp);
    citsubTabs [4] = NULL;
    k = HiddenGroup (m, 0, 0, NULL);

    cpp->citsubGrp [0] = HiddenGroup (k, -1, 0, NULL);
    cpp->authors = CreateAuthorDialog (cpp->citsubGrp [0], 3, -1);
    q = HiddenGroup (cpp->citsubGrp [0], 2, 0, NULL);
    StaticPrompt (q, "Consortium", 0, stdLineHeight, programFont, 'l');
    cpp->consortium = DialogText (q, "", 16, NULL);

    cpp->affil = CreateExtAffilDialog (k, NULL,
                                       &(cpp->citsubGrp [1]),
                                       &(cpp->citsubGrp [2]));

    cpp->citsubGrp [3] = HiddenGroup (k, 0, 2, NULL);
    StaticPrompt (cpp->citsubGrp [3], "Description", 0, 0, programFont, 'c');
    cpp->descr = ScrollText (cpp->citsubGrp [3], 30, 4, programFont, TRUE, NULL);
    Hide (cpp->citsubGrp [3]);

    cpp->citsubGrp [4] = HiddenGroup (k, -1, 0, NULL);
    cpp->date = CreateDateDialog (cpp->citsubGrp [4], NULL);
    Hide (cpp->citsubGrp [4]);

    AlignObjects (ALIGN_CENTER, (HANDLE) cpp->tbs, (HANDLE) cpp->authors,
                  (HANDLE) q, (HANDLE) cpp->citsubGrp [3],
                  (HANDLE) cpp->date, (HANDLE) cpp->affil, NULL);
  }

  return (DialoG) p;
}

typedef struct submitpage {
  DIALOG_MESSAGE_BLOCK
  GrouP           subtype;
  GrouP           releaseGrp;
  GrouP           hup;
  GrouP           releaseDate;
  DialoG          date;
  TexT            userID;
  TexT            comment;
  Char            tool [64];
} SubmitPage, PNTR SubmitPagePtr;

static void ChangeSubtype (GrouP g)

{
  Boolean        update;
  SubmitPagePtr  spp;

  spp = (SubmitPagePtr) GetObjectExtra (g);
  if (spp != NULL && spp->subtype != NULL) {
    update = (Boolean) (GetValue (spp->subtype) == 2);
    if (update) {
      SafeHide (spp->releaseGrp);
    } else {
      SafeShow (spp->releaseGrp);
    }
  }
}

static void ChangeHup (GrouP g)

{
  Boolean        hup;
  SubmitPagePtr  spp;

  spp = (SubmitPagePtr) GetObjectExtra (g);
  if (spp != NULL) {
    hup = (Boolean) (GetValue (spp->hup) == 2);
    if (hup) {
      SafeShow (spp->releaseDate);
    } else {
      SafeHide (spp->releaseDate);
    }
  }
}

static void SubmitBlockPtrToSubmitPage (DialoG d, Pointer data)

{
  Char            ch;
  Int2            i;
  CharPtr         ptr;
  SubmitBlockPtr  sbp;
  SubmitPagePtr   spp;
  CharPtr         str;
  WindoW          tempPort;

  spp = (SubmitPagePtr) GetObjectExtra (d);
  sbp = (SubmitBlockPtr) data;
  if (spp != NULL) {
    tempPort = SavePort (spp->comment);
    SafeSetValue (spp->hup, 2);
    SafeSetValue (spp->subtype, 1);
    /*
    SafeSetTitle (spp->userID, "");
    */
    SafeSetTitle (spp->comment, "");
    PointerToDialog (spp->date, NULL);
    if (sbp != NULL) {
      if (sbp->hup) {
        SafeSetValue (spp->hup, 2);
      } else {
        SafeSetValue (spp->hup, 1);
      }
      i = sbp->subtype;
      if (i > 3 || i < 0) {
        i = 4;
      }
      SafeSetValue (spp->subtype, i);
      /*
      SafeSetTitle (spp->userID, sbp->user_tag);
      */
      str = StringSave (sbp->comment);
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
      SafeSetTitle (spp->comment, str);
      MemFree (str);
      PointerToDialog (spp->date, sbp->reldate);
      spp->tool [0] = '\0';
      if (! StringHasNoText (sbp->tool)) {
        StringNCpy_0 (spp->tool, sbp->tool, sizeof (spp->tool));
      }
    }
    ChangeSubtype (spp->subtype);
    ChangeHup (spp->hup);
    RestorePort (tempPort);
  }
}

static Pointer SubmitPageToSubmitBlockPtr (DialoG d)

{
  Char            ch;
  Int2            i;
  CharPtr         os;
  CharPtr         ptr;
  SubmitBlockPtr  sbp;
  CharPtr         sequin_app_version;
  SubmitPagePtr   spp;
  Char            tmp [64];

  sbp = NULL;
  spp = (SubmitPagePtr) GetObjectExtra (d);
  if (spp != NULL) {
    sbp = SubmitBlockNew ();
    if (sbp != NULL) {
      sbp->hup = (Boolean) (GetValue (spp->hup) == 2);
      if (spp->subtype != NULL) {
        i = GetValue (spp->subtype);
        if (i > 3 || i < 0) {
          i = 255;
        }
        sbp->subtype = i;
      } else {
        sbp->subtype = 1;
      }
      /*
      sbp->user_tag = SaveStringFromText (spp->userID);
      */
      sbp->comment = SaveStringFromTextAndStripNewlines (spp->comment);
      ptr = sbp->comment;
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
      if (StringHasNoText (spp->tool)) {
        sequin_app_version = (CharPtr) GetAppProperty ("SequinAppVersion");
        if (sequin_app_version != NULL) {
          os = GetOpSysString ();
          if (os != NULL) {
            sprintf (tmp, "Sequin %s - %s", sequin_app_version, os);
          } else {
            sprintf (tmp, "Sequin %s", sequin_app_version);
          }
        } else {
          StringCpy (tmp, "Sequin");
        }
        sbp->tool = StringSave (tmp);
      } else {
        sbp->tool = StringSave (spp->tool);
      }
      sbp->reldate = DialogToPointer (spp->date);
    }
  }
  return (Pointer) sbp;
}

static ValNodePtr TestSubmitPage (DialoG d)

{
  DatePtr        dp;
  SubmitPagePtr  spp;
  ValNodePtr     head;

  head = NULL;
  spp = (SubmitPagePtr) GetObjectExtra (d);
  if (spp != NULL) {
    switch (GetValue (spp->subtype)) {
      case 1 :
        switch (GetValue (spp->hup)) {
          case 0 :
            head = AddStringToValNodeChain (head,
              "You must state whether your data can be released immediately or not", 0);
            break;
          case 1 :
            break;
          case 2 :
            dp = DialogToPointer (spp->date);
            if (dp == NULL) {
              head = AddStringToValNodeChain (head, "You must specify a date for release", 0);
            }
            DateFree (dp);
            break;
          default :
            break;
        }
        break;
      case 2 :
        break;
      default :
        break;
    }
  }
  return head;
}

extern DialoG CreateSubmitDataDialog (GrouP h, CharPtr title, Boolean newOnly, Boolean defaultAsUpdate)

{
  GrouP          g;
  GrouP          k;
  GrouP          n;
  GrouP          p;
  GrouP          s;
  GrouP          t;
  SubmitPagePtr  spp;

  p = HiddenGroup (h, 1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  spp = (SubmitPagePtr) MemNew (sizeof (SubmitPage));
  if (spp != NULL) {
    SetObjectExtra (p, spp, StdCleanupExtraProc);
    spp->dialog = (DialoG) p;
    spp->todialog = SubmitBlockPtrToSubmitPage;
    spp->fromdialog = SubmitPageToSubmitBlockPtr;
    spp->testdialog = TestSubmitPage;

    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    n = HiddenGroup (s, -1, 0, NULL);
    SetGroupSpacing (n, 3, 15);

    t = NULL;
    if (! newOnly) {
      t = HiddenGroup (n, 3, 0, NULL);
      SetGroupSpacing (t, 3, 10);
      StaticPrompt (t, "Submission Type", 0, stdLineHeight, programFont, 'l');
      spp->subtype = HiddenGroup (t, 4, 0, ChangeSubtype);
      SetObjectExtra (spp->subtype, spp, NULL);
      RadioButton (spp->subtype, "New");
      RadioButton (spp->subtype, "Update");
      if (defaultAsUpdate) {
        SetValue (spp->subtype, 2);
      } else {
        SetValue (spp->subtype, 1);
      }
    }

    k = HiddenGroup (n, 0, 0, NULL);

    spp->releaseGrp = HiddenGroup (k, -1, 0, NULL);
    g = HiddenGroup (spp->releaseGrp, 3, 0, NULL);
    SetGroupSpacing (g, 3, 10);
    StaticPrompt (g, "May we release this record before publication?",
                  0, stdLineHeight, programFont, 'l');
    spp->hup = HiddenGroup (g, 3, 0, ChangeHup);
    SetObjectExtra (spp->hup, spp, NULL);
    RadioButton (spp->hup, "Yes");
    RadioButton (spp->hup, "No");
    SetValue (spp->hup, 1);
    spp->releaseDate = HiddenGroup (spp->releaseGrp, 10, 0, NULL);
    StaticPrompt (spp->releaseDate, "Release Date: ", 0, popupMenuHeight, programFont, 'l');
    spp->date = CreateDateDialog (spp->releaseDate, NULL);
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) spp->releaseDate, NULL);
    Hide (spp->releaseDate);
    if (GetValue (spp->subtype) == 2) {
      SafeHide (spp->releaseGrp);
    }

    g = NULL;
    if (! newOnly) {
      g = HiddenGroup (n, 0, 2, NULL);
      StaticPrompt (g, "Special Instructions to Database Staff", 0, 0, programFont, 'c');
      spp->comment = ScrollText (g, 30, 4, programFont, TRUE, NULL);
    }

    AlignObjects (ALIGN_CENTER, (HANDLE) spp->releaseGrp,
                  (HANDLE) t, (HANDLE) g, NULL);
  }

  return (DialoG) p;
}

#define SUBMISSION_PAGE   0
#define CONTACT_PAGE      1
#define CITATION_PAGE     2

typedef struct submitform {
  FORM_MESSAGE_BLOCK
  GrouP           pages [3];
  Boolean         visited [3];

  Int2            currentPage;
  Boolean         firstTime;
  Boolean         contactSeen;

  DialoG          submit;
  DialoG          contact;
  DialoG          citsub;

  ButtoN          process;
} SubmitForm, PNTR SubmitFormPtr;

static void SubmitBlockPtrToSubmitForm (ForM f, Pointer data)

{
  SubmitFormPtr   sbfp;
  SubmitBlockPtr  sbp;

  sbfp = (SubmitFormPtr) GetObjectExtra (f);
  sbp = (SubmitBlockPtr) data;
  if (sbfp != NULL) {
    sbfp->firstTime = FALSE;
    if (sbp) {
      PointerToDialog (sbfp->submit, (Pointer) sbp);
      PointerToDialog (sbfp->contact, (Pointer) sbp->contact);
      PointerToDialog (sbfp->citsub, (Pointer) sbp->cit);
    } else {
      PointerToDialog (sbfp->submit, NULL);
      PointerToDialog (sbfp->contact, NULL);
      PointerToDialog (sbfp->citsub, NULL);
    }
  }
}

static Pointer SubmitFormToSubmitBlockPtr (ForM f)

{
  CitSubPtr       csp;
  SubmitFormPtr   sbfp;
  SubmitBlockPtr  sbp;

  sbp = NULL;
  sbfp = (SubmitFormPtr) GetObjectExtra (f);
  if (sbfp != NULL) {
    sbp = (SubmitBlockPtr) DialogToPointer (sbfp->submit);
    if (sbp != NULL) {
      sbp->contact = (ContactInfoPtr) DialogToPointer (sbfp->contact);
      sbp->cit = (CitSubPtr) DialogToPointer (sbfp->citsub);
      if (sbp->contact == NULL || sbp->cit == NULL) {
        sbp = SubmitBlockFree (sbp);
      }
      csp = sbp->cit;
      if (csp->date == NULL) {
        csp->date = DateCurr ();
      }
    }
  }
  return (Pointer) sbp;
}

static ValNodePtr TestSubmitForm (ForM f)

{
  ValNodePtr     head;
  SubmitFormPtr  sbfp;
  ValNodePtr     vnp;

  head = NULL;

  sbfp = (SubmitFormPtr) GetObjectExtra (f);
  if (sbfp != NULL) {

    vnp = TestDialog (sbfp->contact);
    if (vnp != NULL) {
      head = AddStringToValNodeChain (head, "Contact error messages:", 0);
      ValNodeLink (&head, vnp);
    }

    vnp = TestDialog (sbfp->submit);
    if (vnp != NULL) {
      head = AddStringToValNodeChain (head, "Submit error messages:", 0);
      ValNodeLink (&head, vnp);
    }

    vnp = TestDialog (sbfp->citsub);
    if (vnp != NULL) {
      head = AddStringToValNodeChain (head, "Citation error messages:", 0);
      ValNodeLink (&head, vnp);
    }

  }
  return head;
}

static Boolean ReadSubmitBlock (ForM f, CharPtr filename)

{
  Pointer         dataptr;
  Uint2           datatype;
  Uint2           entityID;
  SubmitFormPtr   sbfp;
  SubmitBlockPtr  sbp;
  SeqSubmitPtr    ssp;
  Char            path [PATH_MAX];

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  sbfp = (SubmitFormPtr) GetObjectExtra (f);
  if (sbfp != NULL) {
    if (path [0] != '\0' || GetInputFileName (path, sizeof (path), "", "TEXT")) {
      dataptr = ObjMgrGenericAsnTextFileRead (path, &datatype, &entityID);
      if (dataptr != NULL && entityID > 0) {
        sbp = NULL;
        switch (datatype) {
          case OBJ_SUBMIT_BLOCK :
            sbp = (SubmitBlockPtr) dataptr;
            break;
          case OBJ_SEQSUB :
            ssp = (SeqSubmitPtr) dataptr;
            if (ssp != NULL) {
              sbp = ssp->sub;
            }
            break;
          default :
            break;
        }
        if (sbp != NULL) {
          SubmitBlockPtrToSubmitForm (f, sbp);
        }
        ObjMgrDelete (datatype, dataptr);
        if (sbp != NULL) {
          Update ();
          return TRUE;
        }
      }
    }
  }
  return FALSE;
}

static Boolean WriteSubmitBlock (ForM f, CharPtr filename)

{
  AsnIoPtr        aip;
  SubmitFormPtr   sbfp;
  SubmitBlockPtr  sbp;
  Char            path [PATH_MAX];
#ifdef WIN_MAC
  FILE            *fp;
#endif

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  sbfp = (SubmitFormPtr) GetObjectExtra (f);
  if (sbfp != NULL) {
    if (path [0] != '\0' || GetOutputFileName (path, sizeof (path), NULL)) {
#ifdef WIN_MAC
      fp = FileOpen (path, "r");
      if (fp != NULL) {
        FileClose (fp);
      } else {
        FileCreate (path, "TEXT", "ttxt");
      }
#endif
      aip = AsnIoOpen (path, "w");
      if (aip != NULL) {
        sbp = SubmitFormToSubmitBlockPtr (f);
        SubmitBlockAsnWrite (sbp, aip, NULL);
        AsnIoClose (aip);
        sbp = SubmitBlockFree (sbp);
        return TRUE;
      }
    }
  }
  return FALSE;
}

static Boolean ImportSubmitForm (ForM f, CharPtr filename)

{
  SubmitFormPtr  sbfp;

  sbfp = (SubmitFormPtr) GetObjectExtra (f);
  if (sbfp != NULL) {
    switch (sbfp->currentPage) {
      case SUBMISSION_PAGE :
        return ReadSubmitBlock (f, filename);
      case CONTACT_PAGE :
        return ImportDialog (sbfp->contact, filename);
      case CITATION_PAGE :
        return ImportDialog (sbfp->citsub, filename);
      default :
        break;
    }
  }
  return FALSE;
}

static Boolean ExportSubmitForm (ForM f, CharPtr filename)

{
  SubmitFormPtr  sbfp;

  sbfp = (SubmitFormPtr) GetObjectExtra (f);
  if (sbfp != NULL) {
    switch (sbfp->currentPage) {
      case SUBMISSION_PAGE :
        return WriteSubmitBlock (f, filename);
      case CONTACT_PAGE :
        return ExportDialog (sbfp->contact, filename);
      case CITATION_PAGE :
        return ExportDialog (sbfp->citsub, filename);
      default :
        break;
    }
  }
  return FALSE;
}

static void CopyContactToCitAuthors (SubmitFormPtr sbfp)

{
  AuthListPtr     alp;
  AuthorPtr       ap;
  ContactInfoPtr  cip;
  CitsubPagePtr   cpp;
  ValNodePtr      names;

  if (sbfp == NULL) return;
  cpp = (CitsubPagePtr) GetObjectExtra (sbfp->citsub);
  if (cpp == NULL) return;
  cip = (ContactInfoPtr) DialogToPointer (sbfp->contact);
  if (cip == NULL) return;
  ap = NULL;
  if (cip->contact != NULL) {
    ap = (AuthorPtr) AsnIoMemCopy (cip->contact,
                                   (AsnReadFunc) AuthorAsnRead,
                                   (AsnWriteFunc) AuthorAsnWrite);
  }
  ContactInfoFree (cip);
  if (ap == NULL) return;
  alp = AuthListNew ();
  if (alp != NULL) {
    alp->choice = 1;
    names = ValNodeNew (NULL);
    alp->choice = 1;
    alp->names = names;
    if (names != NULL) {
      names->choice = 1;
      names->data.ptrvalue = ap;
    }
    if (ap == NULL) {
      alp = AuthListFree (alp);
    }
    if (alp != NULL) {
       PointerToDialog (cpp->authors, (Pointer) alp);
    }
  }
  alp = AuthListFree (alp);
}

static void SetSubmitBlockImportExportItems (SubmitFormPtr sbfp)

{
  CitSubPtr  csp;
  IteM       exportItm;
  IteM       importItm;

  if (sbfp != NULL) {
    importItm = FindFormMenuItem ((BaseFormPtr) sbfp, VIB_MSG_IMPORT);
    exportItm = FindFormMenuItem ((BaseFormPtr) sbfp, VIB_MSG_EXPORT);
    switch (sbfp->currentPage) {
      case SUBMISSION_PAGE :
        SafeSetTitle (importItm, "Import Submitter Info...");
        SafeSetTitle (exportItm, "Export Submitter Info...");
        SafeEnable (importItm);
        SafeEnable (exportItm);
        break;
      case CONTACT_PAGE :
        SafeSetTitle (importItm, "Import Contact...");
        SafeSetTitle (exportItm, "Export Contact...");
        SafeEnable (importItm);
        SafeEnable (exportItm);
        break;
      case CITATION_PAGE :
        csp = (CitSubPtr) DialogToPointer (sbfp->citsub);
        if (csp == NULL) {
          CopyContactToCitAuthors (sbfp);
        }
        CitSubFree (csp);
        SafeSetTitle (importItm, "Import Citation...");
        SafeSetTitle (exportItm, "Export Citation...");
        SafeEnable (importItm);
        SafeEnable (exportItm);
        break;
      default :
        break;
    }
  }
}

static void ChangeSubmitBlockPage (VoidPtr data, Int2 newval, Int2 oldval)

{
  /*
  ContactInfoPtr  cip;
  CitSubPtr       csp;
  */
  Int2            i;
  SubmitFormPtr   sbfp;
  Int2            sum;

  sbfp = (SubmitFormPtr) data;
  if (sbfp != NULL) {
    sbfp->currentPage = newval;
    SafeHide (sbfp->pages [oldval]);
    Update ();
    if (sbfp->firstTime && sbfp->contactSeen && newval == CITATION_PAGE) {
      /*
      cip = (ContactInfoPtr) DialogToPointer (sbfp->contact);
      if (cip != NULL) {
        csp = CitSubFromContactInfo (cip);
        if (csp != NULL) {
          PointerToDialog (sbfp->citsub, (Pointer) csp);
          CitSubFree (csp);
        }
        ContactInfoFree (cip);
      }
      */
      sbfp->firstTime = FALSE;
    }
    SafeShow (sbfp->pages [newval]);
    sbfp->visited [sbfp->currentPage] = TRUE;
    if (newval == CONTACT_PAGE) {
      sbfp->contactSeen = TRUE;
    }
    sum = 0;
    for (i = 0; i < 3; i++) {
      if (sbfp->visited [i]) {
        sum++;
      }
    }
    if (sum >= 3) {
      SafeEnable (sbfp->process);
    }
    SetSubmitBlockImportExportItems (sbfp);
    Update ();
  }
}

static CharPtr  submitFormTabs [] = {
  "Submission", "Contact", "Citation", NULL
};

static void SubmitBlockFormMessage (ForM f, Int2 mssg)

{
  SubmitFormPtr  sbfp;

  sbfp = (SubmitFormPtr) GetObjectExtra (f);
  if (sbfp != NULL) {
    switch (mssg) {
      case VIB_MSG_IMPORT :
        ImportSubmitForm (f, NULL);
        break;
      case VIB_MSG_EXPORT :
        ExportSubmitForm (f, NULL);
        break;
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
        if (sbfp->appmessage != NULL) {
          sbfp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static void SubmitBlockFormActivate (WindoW w)

{
  SubmitFormPtr  sbfp;

  sbfp = (SubmitFormPtr) GetObjectExtra (w);
  if (sbfp != NULL) {
    if (sbfp->activate != NULL) {
      sbfp->activate (w);
    }
    SetSubmitBlockImportExportItems (sbfp);
  }
}

static void AcceptSubmitBlockFormButtonProc (ButtoN b)

{
  Boolean         failed;
  SubmitFormPtr   sbfp;
  SubmitBlockPtr  sbp;

  sbp = NULL;
  failed = FALSE;
  sbfp = (SubmitFormPtr) GetObjectExtra (b);
  if (sbfp == NULL) return;
  sbp = (SubmitBlockPtr) DialogToPointer (sbfp->submit);
  if (sbp != NULL) {
    sbp->contact = (ContactInfoPtr) DialogToPointer (sbfp->contact);
    sbp->cit = (CitSubPtr) DialogToPointer (sbfp->citsub);
    if (sbp->contact == NULL && sbp->cit == NULL) {
      failed = TRUE;
      Message (MSG_OK, "Requires contact and citation information");
    } else if (sbp->contact == NULL) {
      failed = TRUE;
      Message (MSG_OK, "Requires contact information");
    } else if (sbp->cit == NULL) {
      failed = TRUE;
      Message (MSG_OK, "Requires citation information");
    } else if (sbp->hup && sbp->reldate == NULL) {
      failed = TRUE;
      Message (MSG_OK, "Requires release date");
    }
  }
  sbp = SubmitBlockFree (sbp);
  if (failed) return;
  Hide (sbfp->form);
  if (b != NULL) {
    if (sbfp != NULL && sbfp->form != NULL && sbfp->actproc != NULL) {
      (sbfp->actproc) (sbfp->form);
    }
  }
  Update ();
  Remove (sbfp->form);
}

static void SubmitBlockFormCleanupProc (GraphiC g, VoidPtr data)

{
  SubmitFormPtr  sbfp;
  Uint2          userkey;

  sbfp = (SubmitFormPtr) data;
  if (sbfp != NULL) {
    if (sbfp->input_entityID > 0 && sbfp->userkey > 0) {
      userkey = sbfp->userkey;
      sbfp->userkey = 0;
      ObjMgrFreeUserData (sbfp->input_entityID, sbfp->procid, sbfp->proctype, userkey);
    }
  }
  StdCleanupExtraProc (g, data);
}

extern ForM CreateSubmitBlockForm (Int2 left, Int2 top, CharPtr title,
                                   Boolean newOnly, Boolean defaultAsUpdate, SubmitBlockPtr sbp,
                                   BtnActnProc cnclproc, FormActnFunc actproc)

{
  ButtoN             b;
  GrouP              c;
  CitSubPtr          csp;
  GrouP              h;
  Int2               i;
  GrouP              j;
  GrouP              q;
  SubmitFormPtr      sbfp;
  StdEditorProcsPtr  sepp;
  DialoG             tbs;
  WindoW             w;

  w = NULL;
  sbfp = MemNew (sizeof (SubmitForm));
  if (sbfp != NULL) {
    if (cnclproc == NULL) {
      w = FixedWindow (left, top, -10, -10, title, StdCloseWindowProc);
    } else {
      w = FixedWindow (left, top, -10, -10, title, NULL);
    }
    SetObjectExtra (w, sbfp, SubmitBlockFormCleanupProc);
    sbfp->form = (ForM) w;
    sbfp->actproc = actproc;
    sbfp->toform = SubmitBlockPtrToSubmitForm;
    sbfp->fromform = SubmitFormToSubmitBlockPtr;
    sbfp->testform = TestSubmitForm;
    sbfp->importform = ImportSubmitForm;
    sbfp->exportform = ExportSubmitForm;
    sbfp->formmessage = SubmitBlockFormMessage;

#ifndef WIN_MAC
    CreateStdEditorFormMenus (w);
#endif

    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      sbfp->activate = sepp->activateForm;
      sbfp->appmessage = sepp->handleMessages;
    }
    SetActivate (w, SubmitBlockFormActivate);

    j = HiddenGroup (w, -1, 0, NULL);
    SetGroupSpacing (j, 10, 10);

    tbs = CreateFolderTabs (j, submitFormTabs, 0, 0, 0,
                            SYSTEM_FOLDER_TAB,
                            ChangeSubmitBlockPage, (Pointer) sbfp);
    sbfp->currentPage = SUBMISSION_PAGE;
    sbfp->firstTime = TRUE;
    sbfp->contactSeen = FALSE;
    for (i = 0; i < 3; i++) {
      sbfp->visited [i] = FALSE;
    }
    sbfp->visited [sbfp->currentPage] = TRUE;

    h = HiddenGroup (w, 0, 0, NULL);

    q = HiddenGroup (h, -1, 0, NULL);
    SetGroupSpacing (q, 10, 20);
    sbfp->submit = CreateSubmitDataDialog (q, "", newOnly, defaultAsUpdate);
    sbfp->pages [SUBMISSION_PAGE] = q;
    Hide (sbfp->pages [SUBMISSION_PAGE]);

    q = HiddenGroup (h, -1, 0, NULL);
    SetGroupSpacing (q, 10, 20);
    sbfp->contact = CreateContactDialog (q, "", NULL);
    sbfp->pages [CONTACT_PAGE] = q;
    Hide (sbfp->pages [CONTACT_PAGE]);

    q = HiddenGroup (h, -1, 0, NULL);
    SetGroupSpacing (q, 10, 20);
    csp = NULL;
    if (sbp != NULL) {
      csp = sbp->cit;
    }
    sbfp->citsub = CreateCitSubDialog (q, "", csp);
    sbfp->pages [CITATION_PAGE] = q;
    Hide (sbfp->pages [CITATION_PAGE]);

    c = HiddenGroup (w, 2, 0, NULL);
    SetGroupSpacing (c, 10, 2);
    sbfp->process = PushButton (c, "Accept", AcceptSubmitBlockFormButtonProc);
    SetObjectExtra (sbfp->process, sbfp, NULL);
    if (cnclproc == NULL) {
      cnclproc = StdCancelButtonProc;
    }
    b = PushButton (c, "Cancel", cnclproc);
    SetObjectExtra (b, sbfp, NULL);

    AlignObjects (ALIGN_CENTER, (HANDLE) sbfp->pages [SUBMISSION_PAGE],
                  (HANDLE) sbfp->pages [CONTACT_PAGE],
                  (HANDLE) sbfp->pages [CITATION_PAGE],
                  (HANDLE) j, (HANDLE) c, NULL);

    RealizeWindow (w);

    Show (sbfp->pages [sbfp->currentPage]);
  }
  return (ForM) w;
}

typedef struct submitblockdlg
{
  DIALOG_MESSAGE_BLOCK
  GrouP           pages [3];
  Boolean         visited [3];

  Int2            currentPage;
  Boolean         firstTime;
  Boolean         contactSeen;

  DialoG          submit;
  DialoG          contact;
  DialoG          citsub;
  DialoG          tbs;
  
  ButtoN          copy_btn;
  ButtoN          import_btn;
  ButtoN          export_btn;

} SubmitBlockDialogData, PNTR SubmitBlockDialogPtr;

static void ChangeSubmitBlockDialogPage (VoidPtr data, Int2 newval, Int2 oldval)

{
  /*
  ContactInfoPtr  cip;
  CitSubPtr       csp;
  */
  Int2            i;
  SubmitBlockDialogPtr   sbfp;
  Int2            sum;

  sbfp = (SubmitBlockDialogPtr) data;
  if (sbfp != NULL) {
    sbfp->currentPage = newval;
    SafeHide (sbfp->pages [oldval]);
    Update ();
    if (sbfp->firstTime && sbfp->contactSeen && newval == CITATION_PAGE) {
      /*
      cip = (ContactInfoPtr) DialogToPointer (sbfp->contact);
      if (cip != NULL) {
        csp = CitSubFromContactInfo (cip);
        if (csp != NULL) {
          PointerToDialog (sbfp->citsub, (Pointer) csp);
          CitSubFree (csp);
        }
        ContactInfoFree (cip);
      }
      */
      sbfp->firstTime = FALSE;
    }
    SafeShow (sbfp->pages [newval]);
    sbfp->visited [sbfp->currentPage] = TRUE;
    if (newval == CONTACT_PAGE) {
      sbfp->contactSeen = TRUE;
    }
    sum = 0;
    for (i = 0; i < 3; i++) {
      if (sbfp->visited [i]) {
        sum++;
      }
    }
    
    /* set title for import button */
    switch (newval)
    {
      case CITATION_PAGE:
        SetTitle (sbfp->import_btn, "Import Citation");
        Show (sbfp->import_btn);
        SetTitle (sbfp->export_btn, "Export Citation");
        Show (sbfp->export_btn);
        break;
      case CONTACT_PAGE:
        SetTitle (sbfp->import_btn, "Import Contact");
        Show (sbfp->import_btn);
        SetTitle (sbfp->export_btn, "Export Contact");
        Show (sbfp->export_btn);
        break;
      case SUBMISSION_PAGE:
        SetTitle (sbfp->import_btn, "Import Submitter Info");
        Show (sbfp->import_btn);
        SetTitle (sbfp->export_btn, "Export Submitter Info");
        Show (sbfp->export_btn);
        break;
      default:
        Hide (sbfp->import_btn);
        Hide (sbfp->export_btn);
        break;
    }
    
    Update ();
  }
}

static void SubmitBlockPointerToDialog (DialoG d, Pointer userdata)
{
  SubmitBlockDialogPtr dlg;
  SubmitBlockPtr       sbp;
  
  dlg = (SubmitBlockDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  sbp = (SubmitBlockPtr) userdata;
  
  if (sbp == NULL)
  {
    PointerToDialog (dlg->submit, NULL);
    PointerToDialog (dlg->contact, NULL);
    PointerToDialog (dlg->citsub, NULL);
  }
  else
  {
    PointerToDialog (dlg->submit, sbp);
    PointerToDialog (dlg->contact, sbp->contact);
    PointerToDialog (dlg->citsub, sbp->cit);
  }

}

static Pointer SubmitBlockDialogToPointer (DialoG d)
{
  SubmitBlockDialogPtr dlg;
  SubmitBlockPtr       sbp;
  
  dlg = (SubmitBlockDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  
  sbp = DialogToPointer (dlg->submit);
  if (sbp != NULL)
  {
    sbp->contact = DialogToPointer (dlg->contact);
    sbp->cit = DialogToPointer (dlg->citsub);
  }
  return sbp;
}

#define MAX_SUBMIT_BLOCK_TABS 5

static ValNodePtr TestSubmitBlockDialog (DialoG d)
{
  SubmitBlockDialogPtr dlg;
  ValNodePtr           errlist = NULL, new_errlist, tmp;
  
  dlg = (SubmitBlockDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  
  new_errlist = TestDialog (dlg->submit);
  tmp = new_errlist;
  while (tmp != NULL)
  {
    tmp->choice *= MAX_SUBMIT_BLOCK_TABS;
    tmp->choice += SUBMISSION_PAGE;
    tmp = tmp->next;
  }
  ValNodeLink (&errlist, new_errlist);
  new_errlist = TestDialog (dlg->contact);
  tmp = new_errlist;
  while (tmp != NULL)
  {
    tmp->choice *= MAX_SUBMIT_BLOCK_TABS;
    tmp->choice += CONTACT_PAGE;
    tmp = tmp->next;
  }
  ValNodeLink (&errlist, new_errlist);
  new_errlist = TestDialog (dlg->citsub);
  tmp = new_errlist;
  while (tmp != NULL)
  {
    tmp->choice *= MAX_SUBMIT_BLOCK_TABS;
    tmp->choice += CITATION_PAGE;
    tmp = tmp->next;
  }
  ValNodeLink (&errlist, new_errlist);
  
  return errlist;
}

static void SubmitBlockDialogMessage (DialoG d, Int2 mssg)

{
  SubmitBlockDialogPtr dlg;
  Int2                 this_page, subpage, pageval;

  dlg = (SubmitBlockDialogPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    if (mssg == VIB_MSG_ENTER)
    {
      SetValue (dlg->tbs, 0);
    }
    else if (mssg > NUM_VIB_MSG)
    {
      pageval = mssg - NUM_VIB_MSG - 1;
      
      this_page = pageval % MAX_SUBMIT_BLOCK_TABS;
      subpage = pageval / MAX_SUBMIT_BLOCK_TABS;
      if (this_page < 3)
      {
        SetValue (dlg->tbs, this_page);
        switch (this_page)
        {
          case SUBMISSION_PAGE:
            SendMessageToDialog (dlg->submit, NUM_VIB_MSG + 1 + subpage);
            break;
          case CONTACT_PAGE:
            SendMessageToDialog (dlg->contact, NUM_VIB_MSG + 1 + subpage);
            break;
          case CITATION_PAGE:
            SendMessageToDialog (dlg->citsub, NUM_VIB_MSG + 1 + subpage);
            break;
        }
      }
    }
  }
}

static Boolean FileToSubmitDialog (DialoG d, CharPtr filename)
{
  SubmitBlockDialogPtr dlg;
  Pointer              dataptr;
  Uint2                datatype;
  Uint2                entityID;
  SubmitBlockPtr       sbp;
  SeqSubmitPtr         ssp;
  Char                 path [PATH_MAX];
  Boolean              rval = FALSE;

  dlg = (SubmitBlockDialogPtr) GetObjectExtra (d);
  
  if (dlg == NULL)
  {
    return FALSE;
  }
  
  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  if (path [0] != '\0' || GetInputFileName (path, sizeof (path), "", "TEXT")) 
  {
    dataptr = ObjMgrGenericAsnTextFileRead (path, &datatype, &entityID);
    if (dataptr == NULL)
    {
      Message (MSG_ERROR, "Unable to read %s", path);
      return FALSE;
    }

    if (datatype == OBJ_SUBMIT_BLOCK)
    {
      sbp = (SubmitBlockPtr) dataptr;
      PointerToDialog (dlg->dialog, sbp);
      rval = TRUE;
    }
    else if (datatype == OBJ_SEQSUB)
    {
      ssp = (SeqSubmitPtr) dataptr;
      if (ssp != NULL)
      {
        PointerToDialog (dlg->dialog, ssp->sub);
      }
      rval = TRUE;
    }
    else
    {
      Message (MSG_ERROR, "Wrong data type in file!");
    }
    ObjMgrDelete (datatype, dataptr);
  }
  
  return rval;
}

static Boolean ReadSubmitBlockDialog (DialoG d, CharPtr filename)

{
  SubmitBlockDialogPtr dlg;
  Boolean              rval = FALSE;
  
  dlg = (SubmitBlockDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return FALSE;
  }
  
  switch (dlg->currentPage) {
    case SUBMISSION_PAGE :
      rval = FileToSubmitDialog (dlg->dialog, filename);
      break;
    case CONTACT_PAGE:
      rval = ReadContactDialog (dlg->contact, filename);
      break;
    case CITATION_PAGE :
      rval = ReadCitsubDialog (dlg->citsub, filename);
      break;
  }
  return rval;
}

static void ImportSubmitBlockBtn (ButtoN b)
{
  SubmitBlockDialogPtr dlg;

  dlg = (SubmitBlockDialogPtr) GetObjectExtra (b);
  if (dlg == NULL)
  {
    return;
  }
  
  ReadSubmitBlockDialog (dlg->dialog, NULL);
}

static Boolean SubmitBlockDialogToFile (DialoG d, CharPtr filename)
{
  SubmitBlockDialogPtr dlg;
  Char                 path [PATH_MAX];
  SubmitBlockPtr       sbp;
  AsnIoPtr             aip;
#ifdef WIN_MAC
  FILE                 *f;
#endif

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));

  dlg = (SubmitBlockDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return FALSE;
  }

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
      sbp = SubmitBlockDialogToPointer (d);
      SubmitBlockAsnWrite (sbp, aip, NULL);
      AsnIoClose (aip);
      sbp = SubmitBlockFree (sbp);
      return TRUE;
    }
  }
  return FALSE;
}

static Boolean WriteSubmitBlockDialog (DialoG d, CharPtr filename)

{
  SubmitBlockDialogPtr dlg;
  Boolean              rval = FALSE;

  dlg = (SubmitBlockDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return FALSE;
  }
  
  switch (dlg->currentPage) {
    case SUBMISSION_PAGE :
      rval = SubmitBlockDialogToFile (dlg->dialog, filename);
      break;
    case CONTACT_PAGE:
      rval = WriteContactDialog (dlg->contact, filename);
      break;
    case CITATION_PAGE :
      rval = WriteCitsubDialog (dlg->citsub, filename);
      break;
  }
  return rval;

}

static void ExportSubmitBlockBtn (ButtoN b)
{
  SubmitBlockDialogPtr dlg;

  dlg = (SubmitBlockDialogPtr) GetObjectExtra (b);
  if (dlg == NULL)
  {
    return;
  }
  
  WriteSubmitBlockDialog (dlg->dialog, NULL);
}

extern DialoG SubmitBlockDialog (GrouP parent, Boolean newOnly, Boolean defaultAsUpdate)
{
  SubmitBlockDialogPtr dlg;
  GrouP                p, h, q, g;
  Int4                 i;
  
  dlg = (SubmitBlockDialogPtr) MemNew (sizeof (SubmitBlockDialogData));
  if (dlg == NULL)
  {
    return NULL;
  }
  
  p = HiddenGroup (parent, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);

  dlg->dialog = (DialoG) p;
  dlg->todialog = SubmitBlockPointerToDialog;
  dlg->fromdialog = SubmitBlockDialogToPointer;
  dlg->dialogmessage = SubmitBlockDialogMessage;
  dlg->testdialog = TestSubmitBlockDialog;
  dlg->importdialog = ReadSubmitBlockDialog;

  
  dlg->tbs = CreateFolderTabs (p, submitFormTabs, 0, 0, 0,
                          SYSTEM_FOLDER_TAB,
                          ChangeSubmitBlockDialogPage, (Pointer) dlg);
  dlg->currentPage = SUBMISSION_PAGE;
  dlg->firstTime = TRUE;
  dlg->contactSeen = FALSE;
  for (i = 0; i < 3; i++) {
    dlg->visited [i] = FALSE;
  }
  dlg->visited [dlg->currentPage] = TRUE;

  h = HiddenGroup (p, 0, 0, NULL);

  q = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (q, 10, 20);
  dlg->submit = CreateSubmitDataDialog (q, "", newOnly, defaultAsUpdate);
  dlg->pages [SUBMISSION_PAGE] = q;
  Hide (dlg->pages [SUBMISSION_PAGE]);

  q = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (q, 10, 20);
  dlg->citsub = CreateCitSubDialog (q, "", NULL);
  dlg->pages [CITATION_PAGE] = q;
  Hide (dlg->pages [CITATION_PAGE]);

  q = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (q, 10, 20);
  dlg->contact = CreateContactDialog (q, "", dlg->citsub);
  dlg->pages [CONTACT_PAGE] = q;
  Hide (dlg->pages [CONTACT_PAGE]);
  
  g = HiddenGroup (p, 2, 0, NULL);
  SetGroupSpacing (g, 10, 10);
  dlg->import_btn = PushButton (g, "Import Submitter Info", ImportSubmitBlockBtn);
  SetObjectExtra  (dlg->import_btn, dlg, NULL);
  dlg->export_btn = PushButton (g, "Export Submitter Info", ExportSubmitBlockBtn);
  SetObjectExtra  (dlg->export_btn, dlg, NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->tbs,
                              (HANDLE) dlg->pages [SUBMISSION_PAGE],
                              (HANDLE) dlg->pages [CONTACT_PAGE],
                              (HANDLE) dlg->pages [CITATION_PAGE],
                              (HANDLE) g,
                              NULL);


  Show (dlg->pages [dlg->currentPage]);
  
  return (DialoG) p;
}

static void SubmitBlockFormActnProc (ForM f)

{
  OMProcControl   ompc;
  SubmitFormPtr   sbfp;
  SubmitBlockPtr  sbp;

  sbfp = (SubmitFormPtr) GetObjectExtra (f);
  if (sbfp != NULL) {
    MemSet ((Pointer) &ompc, 0, sizeof (OMProcControl));
    ompc.input_entityID = sbfp->input_entityID;
    ompc.input_itemID = sbfp->input_itemID;
    ompc.input_itemtype = sbfp->input_itemtype;
    ompc.output_itemtype = sbfp->input_itemtype;
    sbp = (SubmitBlockPtr) FormToPointer (sbfp->form);
    if (sbp != NULL) {
      ompc.output_data = (Pointer) sbp;
      if (ompc.input_entityID == 0) {
        if (! ObjMgrRegister (OBJ_SUBMIT_BLOCK, (Pointer) sbp)) {
          Message (MSG_ERROR, "ObjMgrRegister failed");
        }
      } else if (ompc.input_itemtype != OBJ_SUBMIT_BLOCK) {
        ompc.output_itemtype = OBJ_SUBMIT_BLOCK;
        if (! AttachDataForProc (&ompc, FALSE)) {
          Message (MSG_ERROR, "AttachDataForProc failed");
        }
        ObjMgrSendMsg (OM_MSG_UPDATE, sbfp->input_entityID,
                       sbfp->input_itemID, sbfp->input_itemtype);
      } else {
        if (! ReplaceDataForProc (&ompc, FALSE)) {
          Message (MSG_ERROR, "ReplaceDataForProc failed");
        }
        ObjMgrSendMsg (OM_MSG_UPDATE, sbfp->input_entityID,
                       sbfp->input_itemID, sbfp->input_itemtype);
      }
    }
  }
}

extern Int2 LIBCALLBACK SubmitBlockGenFunc (Pointer data)

{
  ObjMgrDataPtr     omdp;
  OMProcControlPtr  ompcp;
  OMUserDataPtr     omudp;
  SubmitFormPtr     sbfp;
  SubmitBlockPtr    sbp = NULL;
  SeqSubmitPtr      ssp;
  WindoW            w;

  ompcp = (OMProcControlPtr) data;
  sbp = NULL;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_SUBMIT_BLOCK :
      sbp = (SubmitBlockPtr) ompcp->input_data;
      break;
    case OBJ_SEQSUB :
      break;
    case OBJ_SEQSUB_CIT :
      omdp = ObjMgrGetData (ompcp->input_entityID);
      if (omdp != NULL && omdp->datatype == OBJ_SEQSUB) {
        ssp = (SeqSubmitPtr) omdp->dataptr;
        if (ssp != NULL && ssp->datatype == 1) {
          sbp = ssp->sub;
        }
       }
      break;
    case 0 :
      break;
    default :
      return OM_MSG_RET_ERROR;
  }
  omudp = ItemAlreadyHasEditor (ompcp->input_entityID, ompcp->input_itemID,
                                ompcp->input_itemtype, ompcp->proc->procid);
  if (omudp != NULL) {
    sbfp = (SubmitFormPtr) omudp->userdata.ptrvalue;
    if (sbfp != NULL) {
      Select (sbfp->form);
    }
    return OM_MSG_RET_DONE;
  }
  w = (WindoW) CreateSubmitBlockForm (-50, -33,
                                      "Submission Instructions",
                                      FALSE, FALSE, sbp, NULL,
                                      SubmitBlockFormActnProc);
  sbfp = (SubmitFormPtr) GetObjectExtra (w);
  if (sbfp != NULL) {
    sbfp->input_entityID = ompcp->input_entityID;
    sbfp->input_itemID = ompcp->input_itemID;
    sbfp->input_itemtype = ompcp->input_itemtype;
    sbfp->this_itemtype = OBJ_SUBMIT_BLOCK;
    sbfp->this_subtype = 0;
    if (ompcp->input_itemtype == OBJ_SEQSUB_CIT && ompcp->input_itemID == 1) {
      sbfp->input_itemtype = OBJ_SUBMIT_BLOCK;
    }
    sbfp->procid = ompcp->proc->procid;
    sbfp->proctype = ompcp->proc->proctype;
    sbfp->userkey = OMGetNextUserKey ();
    omudp = ObjMgrAddUserData (ompcp->input_entityID, ompcp->proc->procid,
	                           OMPROC_EDIT, sbfp->userkey);

    if (omudp != NULL) {
      omudp->userdata.ptrvalue = (Pointer) sbfp;
      omudp->messagefunc = StdVibrantEditorMsgFunc;
    }
    if (sbp != NULL) {
      PointerToForm (sbfp->form, (Pointer) sbp);
    }
  }
  Show (w);
  Select (w);
  return OM_MSG_RET_DONE;
}

static AuthListPtr NameToAuthList (CharPtr first, CharPtr last)

{
  AuthListPtr  alp;
  AuthorPtr    ap;
  ValNodePtr   names;
  NameStdPtr   nsp;
  PersonIdPtr  pid;

  alp = NULL;
  if (last != NULL && last [0] != '\0') {
    alp = AuthListNew ();
    if (alp != NULL) {
      alp->choice = 1;
      names = ValNodeNew (NULL);
      alp->names = names;
      if (names != NULL) {
        ap = AuthorNew ();
        names->choice = 1;
        names->data.ptrvalue = ap;
        if (ap != NULL) {
          pid = PersonIdNew ();
          ap->name = pid;
          if (pid != NULL) {
            pid->choice = 2;
            nsp = NameStdNew ();
            pid->data = nsp;
            if (nsp != NULL) {
              nsp->names [0] = StringSave (last);
              nsp->names [1] = StringSave (first);
            }
          }
        }
      }
    }
  }
  return alp;
}

extern CitSubPtr CitSubFromContactInfo (ContactInfoPtr cip)

{
  AuthListPtr  alp;
  AuthorPtr    ap;
  CitSubPtr    csp;
  NameStdPtr   nsp;
  PersonIdPtr  pid;

  csp = NULL;
  if (cip != NULL) {
    csp = CitSubNew ();
    if (csp != NULL) {
      alp = NULL;
      ap = cip->contact;
      if (ap != NULL) {
        pid = ap->name;
        if (pid != NULL && pid->choice == 2) {
          nsp = pid->data;
          if (nsp != NULL) {
            alp = NameToAuthList (nsp->names [1], nsp->names [0]);
            if (alp != NULL) {
              alp->affil = AsnIoMemCopy (ap->affil,
                                         (AsnReadFunc) AffilAsnRead,
                                         (AsnWriteFunc) AffilAsnWrite);
            }
          }
        }
      }
      csp->authors = alp;
      csp->date = DateCurr ();
    }
  }
  return csp;
}

typedef struct commentdlg
{
  DIALOG_MESSAGE_BLOCK
  TexT comment_txt;
} CommentDlgData, PNTR CommentDlgPtr;

static void CommentDlgFromPtr (DialoG d, Pointer userdata)
{
  CommentDlgPtr dlg;
  
  dlg = (CommentDlgPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  if (userdata == NULL)
  {
    SetTitle (dlg->comment_txt, "");
  }
  else
  {
    SetTitle (dlg->comment_txt, (CharPtr) userdata);
  }
}

static Pointer CommentFromDlg (DialoG d)
{
  CommentDlgPtr dlg;
  
  dlg = (CommentDlgPtr) GetObjectExtra (d);
  if (dlg == NULL || TextHasNoText (dlg->comment_txt))
  {
    return NULL;
  }
  else
  {
    return SaveStringFromText (dlg->comment_txt);
  }
}

static DialoG CommentDialog (GrouP parent)
{
  CommentDlgPtr dlg;
  GrouP         p;
  PrompT        s;
  
  dlg = (CommentDlgPtr) MemNew (sizeof (CommentDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (parent, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = CommentDlgFromPtr;
  dlg->fromdialog = CommentFromDlg;
  dlg->dialogmessage = NULL;
  dlg->testdialog = NULL;
  
  s = StaticPrompt (p, "Comment", 25 * stdCharWidth, 0, programFont, 'c');
  dlg->comment_txt = ScrollText (p, 25, 10, programFont, TRUE, NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) s, (HANDLE) dlg->comment_txt, NULL);
  
  return (DialoG) p;
}

static ENUM_ALIST(molinfo_biomol_nucX_alist)
  {" ",                      0},
  {"Genomic",                1},
  {"Precursor RNA",          2},
  {"mRNA [cDNA]",            3},
  {"Ribosomal RNA",          4},
  {"Transfer RNA",           5},
  {"Small nuclear RNA",      6},
  {"Small cytoplasmic RNA",  7},
  {"Other-Genetic",          9},
  {"Genomic-mRNA",          10},
  {"cRNA",                  11},
  {"Small nucleolar RNA",   12},
  {"Transcribed RNA",       13},
  {"Other",                255},
END_ENUM_ALIST

typedef struct molinfodlg
{
  DIALOG_MESSAGE_BLOCK
  PopuP moltype;
} SeqSubMolInfoDlgData, PNTR SeqSubMolInfoDlgPtr;

static void SeqSubMolInfoDlgFromPtr (DialoG d, Pointer userdata)
{
  SeqSubMolInfoDlgPtr dlg;
  MolInfoPtr          mip;
  
  dlg = (SeqSubMolInfoDlgPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  mip = (MolInfoPtr) userdata;
  if (mip == NULL)
  {
    SetEnumPopup (dlg->moltype, molinfo_biomol_nucX_alist, 0);
  }
  else
  {
    SetEnumPopup (dlg->moltype, molinfo_biomol_nucX_alist, mip->biomol);
  }  
}

static Pointer MolInfoFromSeqSubMolInfoDlgDlg (DialoG d)
{
  SeqSubMolInfoDlgPtr dlg;
  MolInfoPtr          mip = NULL;
  UIEnum              val;
  
  dlg = (SeqSubMolInfoDlgPtr) GetObjectExtra (d);
  if (dlg != NULL)
  {
    if (GetEnumPopup (dlg->moltype, molinfo_biomol_nucX_alist, &val)
        && val > 0)
    {
      mip = MolInfoNew ();
      if (mip != NULL)
      {
        mip->biomol = val;
      }
    }
  }
  return mip;
}

static DialoG SeqSubMolInfoDlg (GrouP parent)
{
  SeqSubMolInfoDlgPtr dlg;
  GrouP         p;
  PrompT        s;
  
  dlg = (SeqSubMolInfoDlgPtr) MemNew (sizeof (SeqSubMolInfoDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }
  
  p = HiddenGroup (parent, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = SeqSubMolInfoDlgFromPtr;
  dlg->fromdialog = MolInfoFromSeqSubMolInfoDlgDlg;
  dlg->dialogmessage = NULL;
  dlg->testdialog = NULL;
  
  s = StaticPrompt (p, "Molecule Type", 25 * stdCharWidth, 0, programFont, 'c');
  dlg->moltype = PopupList (p, TRUE, NULL);
  InitEnumPopup (dlg->moltype, molinfo_biomol_nucX_alist, NULL);
  SetEnumPopup (dlg->moltype, molinfo_biomol_nucX_alist, 1);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) s, (HANDLE) dlg->moltype, NULL);
  
  return (DialoG) p;
}

typedef void (*TemplatePreviewCallbackProc) PROTO ((Uint2, Uint2, Uint4,
                                                    BlockType, Pointer, 
                                                    Boolean, Int4));

typedef struct templatepreviewdialog
{
  DIALOG_MESSAGE_BLOCK
  DoC            preview_doc;

  Asn2gbJobPtr   ajp;
  SeqSubmitPtr   ssp;
  TemplatePreviewCallbackProc click_callback;
  Pointer        callback_data;
   
} TemplatePreviewDialogData, PNTR TemplatePreviewDialogPtr;

static void CleanUpTemplatePreviewDialog (GraphiC g, VoidPtr data)
{
  TemplatePreviewDialogPtr dlg;
  
  dlg = (TemplatePreviewDialogPtr) data;
  if (dlg != NULL)
  {
    asn2gnbk_cleanup (dlg->ajp);
    dlg->ssp = SeqSubmitFree (dlg->ssp);
  }
  MemFree (dlg);
}

static Int4 GetRefNumForPub (SeqSubmitPtr ssp, Int4 itemID)
{
  SeqEntryPtr  sep;
  BioseqPtr    bsp;
  BioseqSetPtr bssp;
  SeqDescrPtr  sdp = NULL;
  Int4         ref_num = 0;
  ObjValNodePtr ovp;
  
  if (ssp == NULL || ssp->datatype != OBJ_SEQENTRY || ssp->data == NULL)
  {
    return -1;
  }
  
  sep = (SeqEntryPtr) ssp->data;
  if (IS_Bioseq (sep))
  {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp != NULL)
    {
      sdp = bsp->descr;
    }
  }
  else if (IS_Bioseq_set (sep))
  {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL)
    {
      sdp = bssp->descr;
    }
  }
  while (sdp != NULL)
  {
    if (sdp->choice == Seq_descr_pub)
    {
      if (sdp->extended != 0) 
      {
        ovp = (ObjValNodePtr) sdp;
        if (ovp->idx.itemID == itemID)
        {
          return ref_num;
        }
      }
      ref_num++;
    }
    sdp = sdp->next;
  }
  return -1;
}

static void ClickPreviewItem (DoC d, PoinT pt)

{
  TemplatePreviewDialogPtr dlg;
  Int2                     itemClicked = 0;
  BaseBlockPtr             bbp;
  Int4                     par_num, item_num, ref_num = 0;

  dlg = (TemplatePreviewDialogPtr) GetObjectExtra (d);
  if (dlg == NULL || dlg->ajp == NULL || dlg->click_callback == NULL) return;
  MapDocPoint (d, pt, &itemClicked, NULL, NULL, NULL);
  if (itemClicked < 0 || itemClicked > dlg->ajp->numParagraphs)
  {
    return;
  }
  else
  {
    par_num = 0;
    item_num = 1;
    while (par_num < dlg->ajp->numParagraphs
           && dlg->ajp->paragraphArray[par_num]->blocktype != SOURCE_BLOCK
           && dlg->ajp->paragraphArray[par_num]->blocktype !=  REFERENCE_BLOCK
           && dlg->ajp->paragraphArray[par_num]->blocktype != COMMENT_BLOCK
           && dlg->ajp->paragraphArray[par_num]->blocktype !=  SOURCEFEAT_BLOCK
           && dlg->ajp->paragraphArray[par_num]->blocktype != LOCUS_BLOCK)
    {
      par_num++;
    }
    while (item_num < itemClicked && par_num < dlg->ajp->numParagraphs)
    {
      item_num++;
      par_num++;
      while (par_num < dlg->ajp->numParagraphs
             && dlg->ajp->paragraphArray[par_num]->blocktype != SOURCE_BLOCK
             && dlg->ajp->paragraphArray[par_num]->blocktype !=  REFERENCE_BLOCK
             && dlg->ajp->paragraphArray[par_num]->blocktype != COMMENT_BLOCK
             && dlg->ajp->paragraphArray[par_num]->blocktype !=  SOURCEFEAT_BLOCK
             && dlg->ajp->paragraphArray[par_num]->blocktype != LOCUS_BLOCK)
      {
        par_num++;
      }
    }
    if (par_num < dlg->ajp->numParagraphs)
    {
      bbp = dlg->ajp->paragraphArray [par_num];
      if (bbp->blocktype == REFERENCE_BLOCK
          && bbp->itemtype != OBJ_SEQSUB_CIT)
      {
        ref_num = GetRefNumForPub (dlg->ssp, bbp->itemID);
      }
      (dlg->click_callback) (bbp->entityID, bbp->itemtype, bbp->itemID, 
                             bbp->blocktype, dlg->callback_data, dblClick, ref_num);
    }
  }
}

static void ReformatPreviewLocusLine (CharPtr locus_line, SeqDescrPtr mol_sdp)
{
  CharPtr cp, dst;
  Int4    copy_len;
  
  if (StringHasNoText (locus_line) || StringNCmp (locus_line, "LOCUS", 5) != 0)
  {
    return;
  }
  
  cp = locus_line + 5;
  if (mol_sdp == NULL)
  {
    *cp = 0;
    return;
  }
  
  /* skip over whitespace after "LOCUS" */
  cp += StringSpn (cp, " \t");
  dst = cp;
  
  /* remove ID */
  cp += StringCSpn (cp, " \t");
  cp += StringSpn (cp, " \t");
  
  /* remove length */
  cp += StringCSpn (cp, " \t");
  cp += StringSpn (cp, " \t");
  cp += StringCSpn (cp, " \t");
  cp += StringSpn (cp, " \t");
  
  /* copy moltype and topology */
  copy_len = StringCSpn (cp, " \t");
  copy_len += StringSpn (cp + copy_len, " \t");
  copy_len += StringCSpn (cp + copy_len, " \t");
  StringNCpy (dst, cp, copy_len);
  dst [copy_len] = 0;
}

static void SeqSubmitToTemplatePreviewDialog (DialoG d, Pointer userdata)
{
  TemplatePreviewDialogPtr dlg;
  ErrSev                   level;
  Int4                     index;
  BaseBlockPtr             bbp;
  CharPtr                  string;
  BioseqPtr                bsp = NULL;
  SeqEntryPtr              sep;
  SeqSubmitPtr             ssp;
  SeqDescrPtr              mol_sdp = NULL;
  
  dlg = (TemplatePreviewDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }

  Reset (dlg->preview_doc);
  dlg->ajp = asn2gnbk_cleanup (dlg->ajp);

  dlg->ssp = SeqSubmitFree (dlg->ssp);

  ssp = (SeqSubmitPtr) userdata;
  if (ssp == NULL || ssp->datatype != OBJ_SEQENTRY || ssp->data == NULL)
  {
    return;
  }
  
  if (ssp->sub == NULL || ssp->sub->cit == NULL || ssp->sub->contact == NULL)
  {
    AppendText (dlg->preview_doc, "Need citation and contact before preview can be generated!",
                NULL, NULL, programFont);
    UpdateDocument (dlg->preview_doc, 0, 0);
    Update ();
    return;
  }
  
  
  dlg->ssp = (SeqSubmitPtr) AsnIoMemCopy (ssp, (AsnReadFunc) SeqSubmitAsnRead,
                                               (AsnWriteFunc) SeqSubmitAsnWrite);

  sep = (SeqEntryPtr) dlg->ssp->data;
  if (IS_Bioseq (sep))
  {
    bsp = sep->data.ptrvalue;
    if (bsp != NULL)
    {
      mol_sdp = bsp->descr;
      while (mol_sdp != NULL && mol_sdp->choice != Seq_descr_molinfo)
      {
        mol_sdp = mol_sdp->next;
      }
    }
  }
  
  if (bsp != NULL)
  {
    level = ErrSetMessageLevel (SEV_MAX);

    dlg->ajp = asn2gnbk_setup (bsp, NULL, NULL, (FmtType)GENBANK_FMT, SEQUIN_MODE, NORMAL_STYLE, 0, 0, 0, NULL);
    if (dlg->ajp != NULL) {
      for (index = 0; index < dlg->ajp->numParagraphs; index++) {
        bbp = dlg->ajp->paragraphArray [index];
        if (bbp->blocktype == SOURCE_BLOCK
            || bbp->blocktype == REFERENCE_BLOCK
            || bbp->blocktype == COMMENT_BLOCK
            || bbp->blocktype == SOURCEFEAT_BLOCK)
        {
          string = asn2gnbk_format (dlg->ajp, (Int4) index);
          if (string != NULL && *string != '\0') {
            AppendText (dlg->preview_doc, string, NULL, NULL, programFont);
          }
          MemFree (string);
        }
        else if (bbp->blocktype == LOCUS_BLOCK)
        {
          string = asn2gnbk_format (dlg->ajp, (Int4) index);
          if (string != NULL && *string != '\0') {
            ReformatPreviewLocusLine (string, mol_sdp);
            AppendText (dlg->preview_doc, string, NULL, NULL, programFont);
          }
          MemFree (string);
        }
      }
    }

    ErrSetMessageLevel (level);    
  }

  UpdateDocument (dlg->preview_doc, 0, 0);
  Update ();
  
}

static DialoG 
TemplatePreviewDialog 
(GrouP parent, 
 TemplatePreviewCallbackProc click_callback,
 Pointer callback_data)
{
  TemplatePreviewDialogPtr dlg;
  GrouP                    p;
  PrompT                   s;
  
  dlg = (TemplatePreviewDialogPtr) MemNew (sizeof (TemplatePreviewDialogData));
  if (dlg == NULL)
  {
    return NULL;
  }
  
  p = HiddenGroup (parent, -1, 0, NULL);
  SetObjectExtra (p, dlg, CleanUpTemplatePreviewDialog);
  SetGroupSpacing (p, 10, 10);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = SeqSubmitToTemplatePreviewDialog;
  dlg->fromdialog = NULL;
  dlg->dialogmessage = NULL;
  dlg->testdialog = NULL;

  dlg->click_callback = click_callback;
  dlg->callback_data = callback_data;
  
  s = StaticPrompt (p, "Flatfile Preview", 0, 0, programFont, 'c');
  dlg->preview_doc = DocumentPanel (p, stdCharWidth * 35, stdLineHeight * 10);
  SetDocAutoAdjust (dlg->preview_doc, TRUE);
  SetObjectExtra (dlg->preview_doc, dlg, NULL);
  SetDocProcs (dlg->preview_doc, ClickPreviewItem, NULL, NULL, NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) s, (HANDLE) dlg->preview_doc, NULL);
  
  return (DialoG) p;
}

/*****************************************************************************
*
*   PubdescMatch(pub1, pub2)
*
*****************************************************************************/
static Boolean LIBCALL PubdescMatch (PubdescPtr pdp1, PubdescPtr pdp2)
{
    if (pdp1 == NULL && pdp2 == NULL)
    {
        return TRUE;
    }
    else if (pdp1 == NULL || pdp2 == NULL)
    {
        return FALSE;
    }
    else if (pdp1->reftype != pdp2->reftype)
    {
        return FALSE;
    }
    else if (StringCmp (pdp1->name, pdp2->name) != 0
             || StringCmp (pdp1->fig, pdp2->fig) != 0
             || StringCmp (pdp1->maploc, pdp2->maploc) != 0
             || StringCmp (pdp1->seq_raw, pdp2->seq_raw) != 0
             || StringCmp (pdp1->comment, pdp2->comment) != 0)
    {
        return FALSE;
    }
    else if (PubMatch (pdp1->pub, pdp2->pub) != 0)
    {
        return FALSE;
    }
    else
    {
        return TRUE;
    }
}

static Boolean MolInfoMatch (MolInfoPtr mip1, MolInfoPtr mip2)
{
  if (mip1 == NULL && mip2 == NULL)
  {
    return TRUE;
  }
  else if (mip1 == NULL || mip2 == NULL)
  {
    return FALSE;
  }
  else if (mip1->biomol != mip2->biomol
           || mip1->tech != mip2->tech
           || mip1->completeness != mip2->completeness
           || StringCmp (mip1->techexp, mip2->techexp) != 0)
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}

static Boolean SeqDescrListMatch (SeqDescrPtr sdp_list1, SeqDescrPtr sdp_list2)
{
  while (sdp_list1 != NULL && sdp_list2 != NULL)
  {
    if (sdp_list1->choice != sdp_list2->choice)
    {
      return FALSE;
    }
    if (sdp_list1->choice == Seq_descr_comment)
    {
      /* compare comments */
      if (StringCmp (sdp_list1->data.ptrvalue, sdp_list2->data.ptrvalue) != 0)
      {
        return FALSE;
      }
    }
    else if (sdp_list1->choice == Seq_descr_pub)
    {
      /* compare publications */
      if (! PubdescMatch (sdp_list1->data.ptrvalue, sdp_list2->data.ptrvalue))
      {
        return FALSE;
      }
    }
    else if (sdp_list1->choice == Seq_descr_source)
    {
      /* compare sources */
      if (! BioSourceMatch (sdp_list1->data.ptrvalue, sdp_list2->data.ptrvalue))
      {
        return FALSE;
      }
    }
    else if (sdp_list1->choice == Seq_descr_molinfo)
    {
      if (! MolInfoMatch (sdp_list1->data.ptrvalue, sdp_list2->data.ptrvalue))
      {
        return FALSE;
      }
    }
    else
    {
      /* won't compare other kinds of descriptors */
      return FALSE;
    }
    sdp_list1 = sdp_list1->next;
    sdp_list2 = sdp_list2->next;
  }
  if (sdp_list1 == NULL && sdp_list2 == NULL)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

static Boolean SeqSubmitMatch (SeqSubmitPtr ssp1, SeqSubmitPtr ssp2)
{
  BioseqPtr bsp1 = NULL, bsp2 = NULL;
  BioseqSetPtr bssp1 = NULL, bssp2 = NULL;
  SeqEntryPtr  sep1 = NULL, sep2 = NULL;
  SeqDescrPtr  sdp_list1 = NULL, sdp_list2 = NULL;
  
  if (ssp1 == NULL && ssp2 == NULL)
  {
    return TRUE;
  }
  else if (ssp1 == NULL || ssp2 == NULL)
  {
    return FALSE;
  }
  else if (!SubmitBlockMatch (ssp1->sub, ssp2->sub))
  {
    return FALSE;
  }
  else if (ssp1->datatype != ssp2->datatype)
  {
    return FALSE;
  }
  else if (ssp1->data == NULL && ssp2->data == NULL)
  {
    return TRUE;
  }
  else if (ssp1->data == NULL || ssp2->data == NULL)
  {
    return FALSE;
  }
  
  if (ssp1->datatype == OBJ_SEQENTRY)
  {
    sep1 = (SeqEntryPtr) ssp1->data;
    if (IS_Bioseq (sep1) && sep1->data.ptrvalue != NULL)
    {
      bsp1 = (BioseqPtr) sep1->data.ptrvalue;
    }
    else if (IS_Bioseq_set (sep1))
    {
      bssp1 = (BioseqSetPtr) sep1->data.ptrvalue;
    }
    sep2 = (SeqEntryPtr) ssp2->data;
    if (IS_Bioseq (sep2) && sep2->data.ptrvalue != NULL)
    {
      bsp2 = (BioseqPtr) sep2->data.ptrvalue;
    }
    else if (IS_Bioseq_set (sep2))
    {
      bssp2 = (BioseqSetPtr) sep2->data.ptrvalue;
    }
  }
  else if (ssp1->datatype == OBJ_BIOSEQ)
  {
    bsp1 = (BioseqPtr) ssp1->data;
    bsp2 = (BioseqPtr) ssp2->data;
  }
  else if (ssp1->datatype == OBJ_BIOSEQSET)
  {
    bssp1 = (BioseqSetPtr) ssp1->data;
    bssp2 = (BioseqSetPtr) ssp2->data;
  }
  else
  {
    /* won't compare other types */
    return FALSE;
  }
  
  if (bsp1 != NULL)
  {
    sdp_list1 = bsp1->descr;
  }
  else if (bssp1 != NULL)
  {
    sdp_list1 = bssp1->descr;
  }
  
  if (bsp2 != NULL)
  {
    sdp_list2 = bsp2->descr;
  }
  else if (bssp2 != NULL)
  {
    sdp_list2 = bssp2->descr;
  }
  
  return SeqDescrListMatch (sdp_list1, sdp_list2);
  
}

#define TEMPLATE_SUBMITBLOCK_PAGE  0
#define TEMPLATE_ORGANISM_PAGE     1
#define TEMPLATE_MOLINFO_PAGE      2
#define TEMPLATE_COMMENT_PAGE      3
#define TEMPLATE_REFERENCES_PAGE   4
#define NUM_TEMPLATE_PAGES         5

typedef struct submittemplateeditor
{
  FORM_MESSAGE_BLOCK
  DialoG sbt_dlg;
  DialoG src_dlg;
  DialoG comment_dlg;
  DialoG reference_dlg;
  DialoG molinfo_dlg;
  DialoG tbs;
  GrouP  pages [NUM_TEMPLATE_PAGES];
  ButtoN clear_page_btn;
  
  Int2   currentPage;
  Int2   tagFromPage [NUM_TEMPLATE_PAGES];
  
  DialoG preview_dlg;
  
  CloseSubmitTemplateEditorFunc close_proc;
  Pointer      closedata;
  SeqSubmitPtr last_saved_ssp;
  CharPtr      last_loaded_filename;
} SubmitTemplateEditorData, PNTR SubmitTemplateEditorPtr;

static CharPtr  submitTemplateFormTabs [] = {
  "Submit-Block", "Organism", "Molecule", "Comment", "References", NULL
};

static void FixBioSourceLineage (BioSourcePtr biop)
{
  OrgModPtr mod, lastmod = NULL;
  CharPtr   new_str;
  
  if (biop == NULL || biop->org == NULL || biop->org->orgname == NULL
      || StringHasNoText (biop->org->orgname->lineage))
  {
    return;
  }
  
  mod = biop->org->orgname->mod;
  while (mod != NULL && mod->subtype != ORGMOD_old_lineage)
  {
    lastmod = mod;
    mod = mod->next;
  }
  
  if (mod == NULL)
  {
    mod = OrgModNew ();
    if (mod != NULL)
    {
      mod->subtype = ORGMOD_old_lineage;
      mod->subname = StringSave (biop->org->orgname->lineage);
      if (lastmod == NULL)
      {
        biop->org->orgname->mod = mod;
      }
      else
      {
        lastmod->next = mod;
      }
    }
    else if (StringCmp (mod->subname, biop->org->orgname->lineage) == 0)
    {
      /* do nothing, old_lineage already has copy */
    }
    else if (StringHasNoText (mod->subname))
    {
      mod->subname = MemFree (mod->subname);
      mod->subname = StringSave (biop->org->orgname->lineage);
    }
    else
    {
      new_str = (CharPtr) MemNew ((StringLen (mod->subname) 
                                   + StringLen (biop->org->orgname->lineage)
                                   + 3)
                                  * sizeof (Char));
      if (new_str != NULL)
      {
        sprintf (new_str, "%s; %s", mod->subname, biop->org->orgname->lineage);
        mod->subname = MemFree (mod->subname);
        mod->subname = new_str;
      }
    }
  }
}

static Pointer SeqSubmitFromTemplateEditor (ForM f)
{
  SubmitTemplateEditorPtr dlg;
  BioseqPtr        bsp;
  SeqEntryPtr      sep;
  SeqSubmitPtr     ssp;
  SubmitBlockPtr   sbp;
  SeqDescrPtr      sdp, pub_sdp;
  BioSourcePtr     biop;
  CharPtr          comment;
  Uint2            entityID;
  ObjMgrDataPtr    omdp;
  MolInfoPtr       mip;

  dlg = (SubmitTemplateEditorPtr) GetObjectExtra (f);
  if (dlg == NULL)
  {
    return NULL;
  }

  sbp = (SubmitBlockPtr) DialogToPointer (dlg->sbt_dlg);
  
  
  if (sbp == NULL)
  {
    Message (MSG_ERROR, "No submit block!");
    return NULL;
  }
  
  ssp = SeqSubmitNew ();
  if (ssp == NULL)
  {
    Message (MSG_ERROR, "Unable to allocate memory for seq-submit");
    return NULL;
  }

  bsp = BioseqNew ();
  bsp->mol = Seq_mol_dna;
  bsp->repr = Seq_repr_raw;
  bsp->id = MakeSeqID ("lcl|tmp_id");
  sep = SeqEntryNew ();
  sep->choice = 1;
  sep->data.ptrvalue = bsp;
    
  biop = DialogToPointer (dlg->src_dlg);
  if (biop != NULL)
  {
    FixBioSourceLineage (biop);
    sdp = CreateNewDescriptor (sep, Seq_descr_source);
    sdp->data.ptrvalue = biop;
  }
    
  sdp = NULL;   
  comment = DialogToPointer (dlg->comment_dlg);
  if (comment != NULL)
  {
    sdp = CreateNewDescriptor (sep, Seq_descr_comment);
    sdp->data.ptrvalue = comment;
  }
  
  mip = (MolInfoPtr) DialogToPointer (dlg->molinfo_dlg);
  if (mip != NULL)
  {
    sdp = CreateNewDescriptor (sep, Seq_descr_molinfo);
    sdp->data.ptrvalue = mip;
  }
  
  pub_sdp = (SeqDescrPtr) DialogToPointer (dlg->reference_dlg);
  if (pub_sdp != NULL)
  {
    ValNodeLink (&(bsp->descr), pub_sdp);
  }

  ssp->datatype = OBJ_SEQENTRY;
  ssp->data = (Pointer) sep;
  
  SeqMgrLinkSeqEntry (sep, OBJ_SEQSUB, ssp);

  ssp->sub = sbp;

  entityID = ObjMgrGetEntityIDForPointer (ssp);
  omdp = ObjMgrGetData (entityID);
  SeqMgrIndexFeatures (entityID, NULL);
  
  return (Pointer) ssp;  
}

static void DrawTemplatePreview (SubmitTemplateEditorPtr dlg)
{
  SeqSubmitPtr     ssp;

  if (dlg == NULL)
  {
    return;
  }

  ssp = FormToPointer (dlg->form);
  PointerToDialog (dlg->preview_dlg, ssp);
  ssp = SeqSubmitFree (ssp);  
}


static void ChangeSubmitTemplatePage (VoidPtr data, Int2 newval, Int2 oldval)

{
  SubmitTemplateEditorPtr dlg;

  dlg = (SubmitTemplateEditorPtr) data;
  if (dlg != NULL) {
    dlg->currentPage = newval;
    SafeHide (dlg->pages [oldval]);
    Update ();
    switch (dlg->tagFromPage [newval]) {
      case TEMPLATE_SUBMITBLOCK_PAGE :
        SendMessageToDialog (dlg->sbt_dlg, VIB_MSG_ENTER);
        SetTitle (dlg->clear_page_btn, "Clear Submit Block");
        break;
      case TEMPLATE_ORGANISM_PAGE :
        SendMessageToDialog (dlg->src_dlg, VIB_MSG_ENTER);
        SetTitle (dlg->clear_page_btn, "Clear Organism");
        break;
      case TEMPLATE_COMMENT_PAGE :
        SendMessageToDialog (dlg->comment_dlg, VIB_MSG_ENTER);
        SetTitle (dlg->clear_page_btn, "Clear Comment");
        break;
      case TEMPLATE_REFERENCES_PAGE :
        SendMessageToDialog (dlg->reference_dlg, VIB_MSG_ENTER);
        SetTitle (dlg->clear_page_btn, "Clear References");
        break;
      case TEMPLATE_MOLINFO_PAGE :
        SendMessageToDialog (dlg->molinfo_dlg, VIB_MSG_ENTER);
        SetTitle (dlg->clear_page_btn, "Clear Molecule Type");
        break;
      default :
        break;
    }
    SafeShow (dlg->pages [newval]);
    DrawTemplatePreview (dlg);
    Update ();
  }
}
static CharPtr submit_template_editor_quit_warning =
"You have not saved your template since your last change - are you sure you want to lose your changes?";

static void TemplatePreviewClick 
(Uint2      entityID,
 Uint2      itemtype,
 Uint4      itemID,
 BlockType  blocktype,
 Pointer    userdata,
 Boolean    was_double,
 Int4       ref_num)
{
  SubmitTemplateEditorPtr dlg;

  dlg = (SubmitTemplateEditorPtr) userdata;
  if (userdata == NULL)
  {
    return;
  }
  
  switch (blocktype)
  {
    case SOURCE_BLOCK:
      SetValue (dlg->tbs, TEMPLATE_ORGANISM_PAGE);
      break;
    case REFERENCE_BLOCK:
      if (itemtype == OBJ_SEQSUB_CIT)
      {
        SetValue (dlg->tbs, TEMPLATE_SUBMITBLOCK_PAGE);
        SendMessageToDialog (dlg->sbt_dlg, NUM_VIB_MSG + 3);
      }
      else
      {
        SetValue (dlg->tbs, TEMPLATE_REFERENCES_PAGE);
        if (was_double)
        {
          EditPublicationInDialog (dlg->reference_dlg, ref_num);
        }
      }
      break;
    case COMMENT_BLOCK:
      SetValue (dlg->tbs, TEMPLATE_COMMENT_PAGE);
      break;
    case SOURCEFEAT_BLOCK:
      SetValue (dlg->tbs, TEMPLATE_ORGANISM_PAGE);
      break;
    case LOCUS_BLOCK:
      SetValue (dlg->tbs, TEMPLATE_MOLINFO_PAGE);
      break;
    default:
      break;
  }
}

static Boolean OkToQuitSBTEditor (SubmitTemplateEditorPtr dlg)
{
  SeqSubmitPtr this_ssp;
  MsgAnswer    ans = ANS_OK;
  
  if (dlg == NULL)
  {
    return TRUE;
  }
  
  this_ssp = (SeqSubmitPtr) FormToPointer (dlg->form);
  if (!SeqSubmitMatch (this_ssp, dlg->last_saved_ssp))
  {
    ans = Message (MSG_OKC, submit_template_editor_quit_warning);
  }
  
  this_ssp = SeqSubmitFree (this_ssp);
  
  if (ans == ANS_CANCEL)
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}

static SeqDescrPtr ExtractPubCommentMolinfoOrgDescriptors (SeqDescrPtr PNTR sdp_list)
{
  SeqDescrPtr sdp, sdp_prev = NULL, sdp_next;
  SeqDescrPtr pubs_and_comments = NULL;
  
  if (sdp_list == NULL)
  {
    return NULL;
  }
  
  sdp = *sdp_list;
  while (sdp != NULL)
  {
    sdp_next = sdp->next;
    if (sdp->choice == Seq_descr_pub 
        || sdp->choice == Seq_descr_comment
        || sdp->choice == Seq_descr_molinfo
        || sdp->choice == Seq_descr_source)
    {
      if (sdp_prev == NULL)
      {
        *sdp_list = sdp->next;
      }
      else
      {
        sdp_prev->next = sdp_next;
      }
      sdp->next = NULL;
      ValNodeLink (&pubs_and_comments, (ValNodePtr) sdp);
    }
    else
    {
      sdp_prev = sdp;
    }
   
    
    sdp = sdp_next;
  }
  return pubs_and_comments;
}

static Boolean ExportSeqSubmitForm (ForM f, CharPtr filename)
{
  SubmitTemplateEditorPtr dlg;
  SeqSubmitPtr            ssp;
  Char                    path [PATH_MAX];
  AsnIoPtr                aip;
  FILE                    *fp;
  SeqDescrPtr             write_sdp;
  BioseqSetPtr            bssp;
  BioseqPtr               bsp;
  SeqEntryPtr             sep;
  ValNodePtr              err_list;
  SeqDescrPtr             comment_pub_list = NULL;

  dlg = (SubmitTemplateEditorPtr) GetObjectExtra (f);
  if (dlg == NULL)
  {
    return FALSE;
  }
  
  err_list = TestDialog (dlg->sbt_dlg);
  
  if (err_list != NULL)
  {
    SendMessageToDialog (dlg->sbt_dlg, NUM_VIB_MSG + 1 + err_list->choice);
    DisplayErrorMessages ("Submit-Block Errors", err_list);
    ValNodeFreeData (err_list);
    return FALSE;
  }
  
  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  
  if (!GetOutputFileName (path, sizeof (path), NULL))
  {
    return FALSE;
  }
  
  fp = FileOpen (path, "w");
  if (fp == NULL)
  {
    Message (MSG_ERROR, "Unable to open %s", path);
    return FALSE;
  }

  ssp = (SeqSubmitPtr) FormToPointer (dlg->form);
  if (ssp == NULL)
  {
    return FALSE;
  }
  
  if (ssp->datatype == OBJ_SEQENTRY && ssp->data != NULL)
  {
    sep = ssp->data;
    if (IS_Bioseq (sep) && sep->data.ptrvalue != NULL)
    {
      bsp = (BioseqPtr) sep->data.ptrvalue;
      comment_pub_list = ExtractPubCommentMolinfoOrgDescriptors (&(bsp->descr));
    }
    else if (IS_Bioseq_set (sep) && sep->data.ptrvalue != NULL)
    {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      comment_pub_list = ExtractPubCommentMolinfoOrgDescriptors (&(bssp->descr));
    }
  }
  
  aip = AsnIoNew(ASNIO_TEXT_OUT, fp, NULL, NULL, NULL);
  SubmitBlockAsnWrite (ssp->sub, aip, NULL);
  
	AsnIoFlush(aip);
  AsnIoReset(aip);

  /* write out comment and pub descriptors */
  write_sdp = comment_pub_list;
  while (write_sdp != NULL)
  {
    SeqDescAsnWrite (write_sdp, aip, NULL);  
	  AsnIoFlush(aip);
    AsnIoReset(aip);
    write_sdp = write_sdp->next;
  }

  AsnIoClose (aip);
  
  ssp = SeqSubmitFree (ssp);
  comment_pub_list = SeqDescrFree (comment_pub_list);

  dlg->last_saved_ssp = SeqSubmitFree (dlg->last_saved_ssp);
  dlg->last_saved_ssp = (SeqSubmitPtr) FormToPointer (dlg->form);

  dlg->last_loaded_filename = MemFree (dlg->last_loaded_filename);
  dlg->last_loaded_filename = StringSave (path);    

  return TRUE;
}

static void ExportSubmissionTemplate (IteM i)
{
  SubmitTemplateEditorPtr dlg;

  dlg = (SubmitTemplateEditorPtr) GetObjectExtra (i);
  if (dlg == NULL)
  {
    return;
  }
  
  ExportSeqSubmitForm (dlg->form, dlg->last_loaded_filename);
}

static void SaveSubmissionTemplate (ButtoN b)
{
  SubmitTemplateEditorPtr dlg;

  dlg = (SubmitTemplateEditorPtr) GetObjectExtra (b);
  if (dlg == NULL)
  {
    return;
  }
  
  ExportSeqSubmitForm (dlg->form, dlg->last_loaded_filename);
}

static void ClearSubmissionTemplate (SubmitTemplateEditorPtr dlg)
{
  if (dlg == NULL)
  {
    return;
  }
  
  PointerToDialog (dlg->sbt_dlg, NULL);
  PointerToDialog (dlg->comment_dlg, NULL);
  PointerToDialog (dlg->src_dlg, NULL);
  PointerToDialog (dlg->reference_dlg, NULL);
  PointerToDialog (dlg->molinfo_dlg, NULL);
  PointerToDialog (dlg->preview_dlg, NULL);
}

static void ClearSubmissionTemplateItem (IteM i)
{
  SubmitTemplateEditorPtr dlg;

  dlg = (SubmitTemplateEditorPtr) GetObjectExtra (i);
  if (dlg == NULL)
  {
    return;
  }
  
  ClearSubmissionTemplate (dlg);
}

static void ClearSubmissionTemplateButton (ButtoN b)
{
  SubmitTemplateEditorPtr dlg;

  dlg = (SubmitTemplateEditorPtr) GetObjectExtra (b);
  if (dlg == NULL)
  {
    return;
  }
  
  ClearSubmissionTemplate (dlg);
}

static void ClearTemplatePage (ButtoN b)
{
  SubmitTemplateEditorPtr dlg;

  dlg = (SubmitTemplateEditorPtr) GetObjectExtra (b);
  if (dlg == NULL)
  {
    return;
  }
  
  switch (dlg->currentPage)
  {
    case TEMPLATE_SUBMITBLOCK_PAGE :
      PointerToDialog (dlg->sbt_dlg, NULL);
      break;
    case TEMPLATE_ORGANISM_PAGE :
      PointerToDialog (dlg->src_dlg, NULL);    
      break;
    case TEMPLATE_COMMENT_PAGE :
      PointerToDialog (dlg->comment_dlg, NULL);    
      break;
    case TEMPLATE_REFERENCES_PAGE :
      PointerToDialog (dlg->reference_dlg, NULL);    
      break;
    case TEMPLATE_MOLINFO_PAGE:
      PointerToDialog (dlg->molinfo_dlg, NULL);
    default :
      break;
  }
  DrawTemplatePreview (dlg);
  Update ();  
}

static Boolean GetBioSourceFromSeqEntry (SeqEntryPtr sep, BioSourcePtr PNTR pbiop)
{
  BioseqPtr    bsp = NULL;
  BioseqSetPtr bssp = NULL;
  Boolean      rval = TRUE;
  SeqDescrPtr  sdp;
  
  if (pbiop == NULL)
  {
    return FALSE;
  }
  
  if (sep == NULL)
  {
    return TRUE;
  }
  
  if (sep->data.ptrvalue != NULL)
  {
    if (IS_Bioseq (sep))
    {
      bsp = (BioseqPtr) sep->data.ptrvalue;
      sdp = bsp->descr;
    }
    else if (IS_Bioseq_set (sep))
    {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      rval = GetBioSourceFromSeqEntry (bssp->seq_set, pbiop);
      sdp = bssp->descr;
    }
    
    while (rval && sdp != NULL)
    {
      if (sdp->choice == Seq_descr_source && sdp->data.ptrvalue != NULL)
      {
        if (*pbiop != NULL)
        {
          Message (MSG_ERROR, "Found more than one biosource in SeqEntry!");
          rval = FALSE;
        }
        else
        {
          *pbiop = sdp->data.ptrvalue;
        }
      }
      sdp = sdp->next;
    }
    
  }
  if (rval)
  {
    rval = GetBioSourceFromSeqEntry (sep->next, pbiop);
  }
  return rval;
}

static Boolean GetBioSourceFromSeqSubmit (SeqSubmitPtr ssp, BioSourcePtr PNTR pbiop)
{
  SeqEntryPtr sep;
  Boolean     rval;
  
  if (pbiop == NULL || ssp == NULL || ssp->data == NULL || ssp->datatype != OBJ_SEQENTRY)
  {
    return FALSE;
  }
  *pbiop = NULL;
    
  sep = (SeqEntryPtr) ssp->data;  
    
  rval = GetBioSourceFromSeqEntry (sep, pbiop);
  return rval;
}

static Boolean ImportSeqSubmitForm (ForM f, CharPtr filename)
{
  SubmitTemplateEditorPtr dlg;
  SeqSubmitPtr            ssp = NULL;
  Char                    path [PATH_MAX];
  Boolean                 bad_data_found = FALSE;
  Pointer                 dataptr;
  Uint2                   datatype;
  Uint2                   entityID;
  FILE                    *fp;
  SubmitBlockPtr          sbp = NULL;
  SeqDescrPtr             sdp_comment, sdp_pub = NULL;
  CharPtr                 comment = NULL;
  BioSourcePtr            biop = NULL;
  MolInfoPtr              mip = NULL;
  
  dlg = (SubmitTemplateEditorPtr) GetObjectExtra (f);
  if (dlg == NULL)
  {
    return FALSE;
  }
  
  if (!OkToQuitSBTEditor (dlg))
  {
    return FALSE;
  }
  
  
  if (!GetInputFileName (path, sizeof (path), NULL, NULL))
  {
    return FALSE;
  }
  
  fp = FileOpen (path, "r");
  if (fp == NULL)
  {
    Message (MSG_ERROR, "Unable to open %s", path);
    return FALSE;
  }
  
  dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, &entityID,
                                    FALSE, FALSE, FALSE, FALSE);
  
  while (dataptr != NULL && ! bad_data_found)
  {
    if (datatype == OBJ_SEQSUB)
    {
      if (ssp != NULL)
      {
        Message (MSG_ERROR, "Found more than one SeqSubmit!  Bad file!");
        bad_data_found = TRUE;
      }
      else if (sbp != NULL)
      {
        Message (MSG_ERROR, "Found Separate SeqSubmit and SubmitBlock!  Bad file!");
        bad_data_found = TRUE;
      }
      else
      {
        ssp = (SeqSubmitPtr) dataptr;
        if (ssp->data == NULL || ssp->datatype != OBJ_SEQENTRY)
        {
          Message (MSG_ERROR, "SeqSubmit contained data other than a SeqEntry! Bad file!");
          bad_data_found = TRUE;
        }
        if (!GetBioSourceFromSeqSubmit (ssp, &biop))
        {
          bad_data_found = TRUE;
        }
      }
    }
    else if (datatype == OBJ_SUBMIT_BLOCK)
    {
      if (sbp != NULL)
      {
        Message (MSG_ERROR, "Found more than one SubmitBlock!  Bad file!");
        bad_data_found = TRUE;
      }
      else if (ssp != NULL)
      {
        Message (MSG_ERROR, "Found Separate SeqSubmit and SubmitBlock!  Bad file!");
        bad_data_found = TRUE;
      }
      else
      {
        sbp = (SubmitBlockPtr) dataptr;
      }
    }
    else if (datatype == OBJ_SEQDESC)
    {
      sdp_comment = (SeqDescrPtr) dataptr;
      if (sdp_comment->choice == Seq_descr_pub)
      {
        ValNodeLink (&sdp_pub, sdp_comment);
        sdp_comment = NULL;
      }
      else if (sdp_comment->choice == Seq_descr_molinfo)
      {
        if (mip != NULL)
        {
          Message (MSG_ERROR, "Found more than one molinfo descriptor!  Bad file!");
          bad_data_found = TRUE;
        }
        else
        {
          mip = (MolInfoPtr) sdp_comment->data.ptrvalue;
        }
        sdp_comment = NULL;
      }
      else if (sdp_comment->choice == Seq_descr_source)
      {
        if (biop != NULL)
        {
          Message (MSG_ERROR, "Found more than one biosource!  Bad file!");
          bad_data_found = TRUE;
        }
        else
        {
          biop = sdp_comment->data.ptrvalue;
          sdp_comment->data.ptrvalue = NULL;
        }
      }
      else if (sdp_comment->choice != Seq_descr_comment)
      {
        Message (MSG_ERROR, "Found descriptor other than comment or pub!  Bad file!");
        bad_data_found = TRUE;
      }
      else if (comment != NULL)
      {
        Message (MSG_ERROR, "Found more than one comment!  Cannot edit this file!");
        bad_data_found = TRUE;
      }
      else
      {
        comment = StringSave (sdp_comment->data.ptrvalue);
      }
      sdp_comment = SeqDescrFree (sdp_comment);
    }
    else
    {
      Message (MSG_ERROR, "Read unrecognized data!  Bad file!");
    }
    
    dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, &entityID,
                                    FALSE, FALSE, FALSE, FALSE);    
  }
  
  FileClose (fp);
  
  if (!bad_data_found)
  {
    PointerToDialog (dlg->comment_dlg, comment);
    PointerToDialog (dlg->src_dlg, biop);
    if (ssp == NULL)
    {
      PointerToDialog (dlg->sbt_dlg, sbp);
    }
    else
    {
      PointerToDialog (dlg->sbt_dlg, ssp->sub);
    }
    
    PointerToDialog (dlg->reference_dlg, sdp_pub);
    
    if (mip == NULL)
    {
      /* necessary so that we don't carry over values from the last file */
      mip = MolInfoNew ();
    }
    PointerToDialog (dlg->molinfo_dlg, mip);
    
    DrawTemplatePreview (dlg);
    dlg->last_loaded_filename = MemFree (dlg->last_loaded_filename);
    dlg->last_loaded_filename = StringSave (path);
    dlg->last_saved_ssp = SeqSubmitFree (dlg->last_saved_ssp);
    dlg->last_saved_ssp = FormToPointer (dlg->form);
  }
    
  ssp = SeqSubmitFree (ssp);
  sbp = SubmitBlockFree (sbp);
  comment = MemFree (comment);
  sdp_pub = SeqDescrFree (sdp_pub);
  mip = MolInfoFree (mip);
  
  return bad_data_found;  
}

static void ImportSubmissionTemplateFileMenu (IteM i)
{
  SubmitTemplateEditorPtr dlg;
  
  dlg = (SubmitTemplateEditorPtr) GetObjectExtra (i);
  if (dlg == NULL)
  {
    return;
  }
    
  ImportSeqSubmitForm (dlg->form, NULL);
}

static void ImportSubmissionTemplateFileButton (ButtoN b)
{
  SubmitTemplateEditorPtr dlg;
  
  dlg = (SubmitTemplateEditorPtr) GetObjectExtra (b);
  if (dlg == NULL)
  {
    return;
  }
  
  ImportSeqSubmitForm (dlg->form, NULL);
}

static void CleanUpSubmitTemplateEditor (GraphiC g, VoidPtr data)
{
  SubmitTemplateEditorPtr dlg;
  
  dlg = (SubmitTemplateEditorPtr) data;
  if (dlg != NULL)
  {
    dlg->last_saved_ssp = SeqSubmitFree (dlg->last_saved_ssp);
    dlg->last_loaded_filename = MemFree (dlg->last_loaded_filename);
  }
  MemFree (dlg);
}

static void CloseSBTWindowProc (WindoW w)
{
  SubmitTemplateEditorPtr dlg;
  
  dlg = (SubmitTemplateEditorPtr) GetObjectExtra (w);
  if (OkToQuitSBTEditor (dlg))
  {
    ObjMgrFreeUserData(dlg->input_entityID, 0, 0, 0);
    dlg->close_proc (dlg->closedata, w);
  }
}

static void QuitSbt (IteM i)
{
  SubmitTemplateEditorPtr dlg;
  
  dlg = (SubmitTemplateEditorPtr) GetObjectExtra (i);
  if (OkToQuitSBTEditor (dlg))
  {
    ObjMgrFreeUserData(dlg->input_entityID, 0, 0, 0);
    dlg->close_proc (dlg->closedata, (WindoW) dlg->form);
  }
}

static void SetupSubmissionTemplateMenus (WindoW w, SubmitTemplateEditorPtr dlg)
{
  MenU edit_menu;
  IteM localItem;
  
  if (w == NULL)
  {
    return;
  }
  
  edit_menu = PulldownMenu (w, "File");
  
  localItem = CommandItem (edit_menu, "Load Submission Template", ImportSubmissionTemplateFileMenu);
  SetObjectExtra (localItem, dlg, NULL);
  localItem = CommandItem (edit_menu, "Save Submission Template", ExportSubmissionTemplate);
  SetObjectExtra (localItem, dlg, NULL);
  localItem = CommandItem (edit_menu, "Clear Submission Template", ClearSubmissionTemplateItem);
  SetObjectExtra (localItem, dlg, NULL);
  
  localItem = CommandItem (edit_menu, "Quit", QuitSbt);
  SetObjectExtra (localItem, dlg, NULL);
}

static void DrawTemplatePreviewBtn (ButtoN b)
{
  SubmitTemplateEditorPtr dlg;
  
  dlg = (SubmitTemplateEditorPtr) GetObjectExtra (b);
  DrawTemplatePreview (dlg);
}

static Int2 LIBCALLBACK SeqSubFormMsgFunc (OMMsgStructPtr ommsp)
{
  WindoW             currentport,
                     temport;
  OMUserDataPtr      omudp;
  SubmitTemplateEditorPtr dlg = NULL;
  
  omudp = (OMUserDataPtr)(ommsp->omuserdata);
  if (omudp == NULL) return OM_MSG_RET_ERROR;
  dlg = (SubmitTemplateEditorPtr) omudp->userdata.ptrvalue;
  if (dlg == NULL) return OM_MSG_RET_ERROR;

  currentport = (WindoW)dlg->form;
  temport = SavePort (currentport);
  UseWindow (currentport);
  Select (currentport);
  switch (ommsp->message) 
  {
      case OM_MSG_UPDATE:
          SendMessageToDialog (dlg->reference_dlg, VIB_MSG_REDRAW);
          DrawTemplatePreview (dlg);
          break;

      default:
          break;
  }
  RestorePort (temport);
  UseWindow (temport);
  return OM_MSG_RET_OK;
}

extern ForM 
CreateSubmitTemplateEditorForm 
(Int2 left, Int2 top, CharPtr title,
 CloseSubmitTemplateEditorFunc close_proc, Pointer closedata)
{

  WindoW w;
  SubmitTemplateEditorPtr dlg;
  GrouP                   h, tab_page_grp, k, preview_grp, op_grp;
  Int4                    page_num = 0;
  ButtoN                  b, b2;
  OMUserDataPtr           omudp;
 
  dlg = (SubmitTemplateEditorPtr) MemNew (sizeof (SubmitTemplateEditorData));
  if (dlg == NULL)
  {
    return NULL;
  }
  
  w = FixedWindow (left, top, -10, -10, "Submission Template Editor", CloseSBTWindowProc);
  SetObjectExtra (w, dlg, CleanUpSubmitTemplateEditor);
  dlg->form = (ForM) w;
  dlg->fromform = SeqSubmitFromTemplateEditor;
  dlg->importform = ImportSeqSubmitForm;
  dlg->exportform = ExportSeqSubmitForm;
  
  dlg->close_proc = close_proc;
  dlg->closedata = closedata;
  
  dlg->last_saved_ssp = NULL;
  dlg->last_loaded_filename = NULL;
  
  SetupSubmissionTemplateMenus (w, dlg);
  
  k = HiddenGroup (w, -1, 0, NULL);
  h = HiddenGroup (k, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);


  dlg->tbs = CreateFolderTabs (h, submitTemplateFormTabs, TEMPLATE_SUBMITBLOCK_PAGE,
                                    0, 0, SYSTEM_FOLDER_TAB,
                                    ChangeSubmitTemplatePage, (Pointer) dlg);
  
  tab_page_grp = HiddenGroup (h, 0, 0, NULL);
  dlg->currentPage = TEMPLATE_SUBMITBLOCK_PAGE;
  page_num = 0;

  /* submit block page */
  dlg->sbt_dlg = SubmitBlockDialog (tab_page_grp, TRUE, FALSE);

  dlg->pages [page_num] = (GrouP)dlg->sbt_dlg;
  dlg->tagFromPage [page_num] = TEMPLATE_SUBMITBLOCK_PAGE;
  Hide (dlg->pages [page_num]);
  page_num++;

  /* organism page */
  dlg->src_dlg = BioSourceDialog (tab_page_grp);
  dlg->pages [page_num] = (GrouP)dlg->src_dlg;
  dlg->tagFromPage [page_num] = TEMPLATE_ORGANISM_PAGE;
  Hide (dlg->pages [page_num]);
  page_num++;

  /* molinfo page */
  dlg->molinfo_dlg = SeqSubMolInfoDlg (tab_page_grp);
  dlg->pages [page_num] = (GrouP)dlg->molinfo_dlg;
  dlg->tagFromPage [page_num] = TEMPLATE_MOLINFO_PAGE;
  Hide (dlg->pages [page_num]);
  page_num++;

  /* comment page */
  dlg->comment_dlg = CommentDialog (tab_page_grp);
  dlg->pages [page_num] = (GrouP)dlg->comment_dlg;
  dlg->tagFromPage [page_num] = TEMPLATE_COMMENT_PAGE;
  Hide (dlg->pages [page_num]);
  page_num++;
  
  /* reference page */
  dlg->reference_dlg = PublicationListDialog (tab_page_grp);
  dlg->pages [page_num] = (GrouP)dlg->reference_dlg;
  dlg->tagFromPage [page_num] = TEMPLATE_REFERENCES_PAGE;
  Hide (dlg->pages [page_num]);
  page_num++;
  
  Show (dlg->pages [dlg->currentPage]);
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->sbt_dlg,
                              (HANDLE) dlg->src_dlg,
                              (HANDLE) dlg->comment_dlg,
                              (HANDLE) dlg->reference_dlg,
                              (HANDLE) dlg->molinfo_dlg,
                              NULL);
  
  dlg->clear_page_btn = PushButton (h, "Clear Submit Block", ClearTemplatePage);
  SetObjectExtra (dlg->clear_page_btn, dlg, NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->tbs,
                              (HANDLE) tab_page_grp,
                              (HANDLE) dlg->clear_page_btn,
                              NULL);

  preview_grp = HiddenGroup (k, -1, 0, NULL);
  dlg->preview_dlg = TemplatePreviewDialog (preview_grp, TemplatePreviewClick, dlg);
  b = PushButton (preview_grp, "Update Flatfile Preview", DrawTemplatePreviewBtn);
  SetObjectExtra (b, dlg, NULL);
  op_grp = HiddenGroup (preview_grp, 4, 0, NULL);
  SetGroupSpacing (op_grp, 10, 10);
  b2 = PushButton (op_grp, "Load Submission Template File", ImportSubmissionTemplateFileButton);
  SetObjectExtra (b2, dlg, NULL);
  b2 = PushButton (op_grp, "Save Submission Template File", SaveSubmissionTemplate);
  SetObjectExtra (b2, dlg, NULL);
  b2 = PushButton (op_grp, "Clear Submission Template", ClearSubmissionTemplateButton);
  SetObjectExtra (b2, dlg, NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->preview_dlg, 
                              (HANDLE) b, 
                              (HANDLE) op_grp, 
                              NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) h, (HANDLE) preview_grp, NULL);
                              
  /* register to receive update messages */
  dlg->userkey = OMGetNextUserKey ();
  dlg->procid = 0;
  dlg->proctype = OMPROC_EDIT;
  omudp = ObjMgrAddUserData (dlg->input_entityID, dlg->procid, dlg->proctype, dlg->userkey);
  if (omudp != NULL) {
    omudp->userdata.ptrvalue = (Pointer) dlg;
    omudp->messagefunc = SeqSubFormMsgFunc;
  }
  
  dlg->last_saved_ssp = (SeqSubmitPtr) FormToPointer (dlg->form);
                              
  return (ForM) w;  
}

