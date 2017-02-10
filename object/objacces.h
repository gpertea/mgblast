/*  objacces.h
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE                          
*               National Center for Biotechnology Information
*                                                                          
*  This software/database is a "United States Government Work" under the   
*  terms of the United States Copyright Act.  It was written as part of    
*  the author's official duties as a United States Government employee and 
*  thus cannot be copyrighted.  This software/database is freely available 
*  to the public for use. The National Library of Medicine and the U.S.    
*  Government have not placed any restriction on its use or reproduction.  
*                                                                          
*  Although all reasonable efforts have been taken to ensure the accuracy  
*  and reliability of the software and data, the NLM and the U.S.          
*  Government do not and cannot warrant the performance or results that    
*  may be obtained by using this software or data. The NLM and the U.S.    
*  Government disclaim all warranties, express or implied, including       
*  warranties of performance, merchantability or fitness for any particular
*  purpose.                                                                
*                                                                          
*  Please cite the author in any work or product based on this material.   
*
* ===========================================================================
*
* File Name:  objacces.h
*
* Author:  James Ostell
*   
* Version Creation Date: 1/1/91
*
* $Revision: 6.0 $
*
* File Description:  Object manager interface for module NCBI-Access
*
* Modifications:  
* --------------------------------------------------------------------------
* Date	   Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* $Log: objacces.h,v $
* Revision 6.0  1997/08/25 18:49:07  madden
* Revision changed to 6.0
*
* Revision 4.1  1997/06/19 18:40:35  vakatov
* [WIN32,MSVC++]  Adopted for the "NCBIOBJ.LIB" DLL'ization
*
* Revision 4.0  1995/07/26 13:48:06  ostell
* force revision to 4.0
*
 * Revision 3.1  1995/05/15  21:22:00  ostell
 * added Log line
 *
*
*
* ==========================================================================
*/

#ifndef _NCBI_Access_
#define _NCBI_Access_

#ifndef _ASNTOOL_
#include <asn.h>
#endif

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*****************************************************************************
*
*   loader
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL AccessAsnLoad PROTO((void));

/*****************************************************************************
*
*   internal structures for NCBI-Access objects
*
*****************************************************************************/

/*****************************************************************************
*
*   Link-set
*
*****************************************************************************/
                                /* NOTE: in the uids array, there are num+1
                                   uids, with the last being 0 */
typedef struct linkset 
{
    Int4 num;                    /* number of uids */
    Int4Ptr uids;                /* uids of type "target" */
    Int4Ptr weights;             /* weights, if neighbors */
} LinkSet, PNTR LinkSetPtr;

NLM_EXTERN LinkSetPtr LIBCALL LinkSetNew PROTO((void));
NLM_EXTERN LinkSetPtr LIBCALL LinkSetFree PROTO((LinkSetPtr ufp));
NLM_EXTERN LinkSetPtr LIBCALL LinkSetAsnRead PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN Boolean    LIBCALL LinkSetAsnWrite PROTO((LinkSetPtr ufp, AsnIoPtr aip, AsnTypePtr atp));

#ifdef __cplusplus
}
#endif

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif

