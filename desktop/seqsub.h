/*   seqsub.h
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
* File Name:  seqsub.h
*
* Author:  Jonathan Kans
*
* Version Creation Date:   1/22/95
*
* $Revision: 6.4 $
*
* File Description: 
*
* Modifications:  
* --------------------------------------------------------------------------
* $Log: seqsub.h,v $
* Revision 6.4  2005/05/24 13:28:20  bollin
* added optional button to Contact Info dialog name and affiliation tabs to
* allow user to copy name or affiliation into the supplied citation dialog
*
* Revision 6.3  2005/05/19 20:25:35  bollin
* added dialog for editing Submission Template files.  Still in progress.
*
* Revision 6.2  2002/12/20 18:31:08  kans
* added REGISTER_SEQSUBCIT_EDIT, support for OBJ_SEQSUB_CIT in SubmitBlockGenFunc
*
* Revision 6.1  2000/07/08 20:44:01  vakatov
* Get all "#include" out of the 'extern "C" { }' scope;  other cleanup...
*
* ==========================================================================
*/

#ifndef _SEQSUB_
#define _SEQSUB_

#include <dlogutil.h>
#include <objsub.h>


#ifdef __cplusplus
extern "C" {
#endif

extern DialoG CreateContactDialog (GrouP h, CharPtr title, DialoG citsub_dlg);
extern DialoG CreateCitSubDialog (GrouP h, CharPtr title, CitSubPtr csp);
extern DialoG CreateSubmitDataDialog (GrouP h, CharPtr title, Boolean newOnly, Boolean defaultAsUpdate);

#define REGISTER_SUBMITBLOCK_EDIT ObjMgrProcLoad(OMPROC_EDIT,"Edit SubmitBlock","SubmitBlock",OBJ_SUBMIT_BLOCK,0,OBJ_SUBMIT_BLOCK,0,NULL,SubmitBlockGenFunc,PROC_PRIORITY_DEFAULT)
#define REGISTER_SEQSUBCIT_EDIT ObjMgrProcLoad(OMPROC_EDIT,"Edit SeqSubCit","SeqSubCit",OBJ_SEQSUB_CIT,0,OBJ_SEQSUB_CIT,0,NULL,SubmitBlockGenFunc,PROC_PRIORITY_DEFAULT)

extern ForM CreateSubmitBlockForm (Int2 left, Int2 top, CharPtr title,
                                   Boolean newOnly, Boolean defaultAsUpdate, SubmitBlockPtr sbp,
                                   BtnActnProc cnclproc, FormActnFunc actproc);
extern Int2 LIBCALLBACK SubmitBlockGenFunc (Pointer data);

extern CitSubPtr CitSubFromContactInfo (ContactInfoPtr cip);

extern DialoG SubmitBlockDialog (GrouP parent, Boolean newOnly, Boolean defaultAsUpdate);

typedef void (*CloseSubmitTemplateEditorFunc) PROTO ((Pointer, WindoW));

extern ForM CreateSubmitTemplateEditorForm (Int2 left, Int2 top, CharPtr title,
                                            CloseSubmitTemplateEditorFunc close_proc,
                                            Pointer closedata);


#ifdef __cplusplus
}
#endif

#endif /* ndef _SEQSUB_ */

