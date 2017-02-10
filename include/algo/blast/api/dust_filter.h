/* $Id: dust_filter.h,v 1.1 2005/07/19 14:16:54 madden Exp $
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
* Author: Tom Madden
*
*/

/** @file dust_filter.h
* Dust filtering API for the C version of rewritten BLAST engine.
*/

#ifndef _DUST_FILTER_ 
#define _DUST_FILTER_ 

#ifdef __cplusplus
extern "C" {
#endif

#include <objloc.h>
#include <algo/blast/api/blast_options_api.h>

/** @addtogroup CToolkitAlgoBlast
 *
 * @{
 */

/** Finds dust filtering locations.
 * Does nothing if dust filtering is not required.
 * @param query_seqloc List of query locations [in]
 * @param options object containing dust filter options [in]
 * @param mask_loc List of dust mask locations [out]
 * @return Status.
 */
Int2
Blast_FindDustSeqLoc(SeqLoc* query_seqloc,
                             const SBlastOptions* options,
                             SeqLoc* *mask_loc);

/* @} */

#ifdef __cplusplus
}
#endif

#endif
