/* $Id: blast_message_api.h,v 1.1 2006/04/26 12:44:09 madden Exp $
***************************************************************************
*                                                                         *
*                             COPYRIGHT NOTICE                            *
*                                                                         *
* This software/database is categorized as "United States Government      *
* Work" under the terms of the United States Copyright Act.  It was       *
* produced as part of the author's official duties as a Government        *
* employee and thus can not be copyrighted.  This software/database is    *
* freely available to the public for use without a copyright notice.      *
* Restrictions can not be placed on its present or future use.            *
*                                                                         *
* Although all reasonable efforts have been taken to ensure the accuracy  *
* and reliability of the software and data, the National Library of       *
* Medicine (NLM) and the U.S. Government do not and can not warrant the   *
* performance or results that may be obtained by using this software,     *
* data, or derivative works thereof.  The NLM and the U.S. Government     *
* disclaim any and all warranties, expressed or implied, as to the        *
* performance, merchantability or fitness for any particular purpose or   *
* use.                                                                    *
*                                                                         *
* In any work or product derived from this material, proper attribution   *
* of the author(s) as the source of the software or data would be         *
* appreciated.                                                            *
***************************************************************************
* Author: Tom Madden
*
*/

/** @file blast_messge_api.h
 * API for C Toolkit BLAST messages
 */

#ifndef _BLAST_MESSAGE_API_H_
#define _BLAST_MESSAGE_API_H_

#ifdef __cplusplus
extern "C" {
#endif

#ifndef NCBI_C_TOOLKIT
#define NCBI_C_TOOLKIT
#endif

#include <ncbierr.h>
#include <objloc.h>
#include <algo/blast/core/blast_message.h>
#include <algo/blast/core/blast_query_info.h>

/** @addtogroup CToolkitAlgoBlast
 *
 * @{
 */

/** C Toolkit specific API for messages from BLAST. 
*/
typedef struct SBlastMessage {
    struct SBlastMessage *next; /**< next message. */
    ErrSev sev; /**< one of SEV_NONE=0, SEV_INFO, SEV_WARNING, SEV_ERROR, SEV_REJECT, SEV_FATAL, SEV_MAX */
    char* message; /**< error string. Should be set, but may be NULL if unknown error */
    SeqId* query_id;  /**< used to label query that produced above message, may be NULL. */
    Boolean believe_query;  /**< specifies whether Query ID was parsed. */
} SBlastMessage;


/** Appends an SBlastMessage to the chain of messages, or creates new one
 * if NULL.
 * @param blast_message linked list to be appended to. [in|out]
 * @param sev severity level [in]
 * @param message text to be printed for user. [in]
 * @param seqid describes the query having problem, may be NULL if message applies to setup [in]
 * @param believe_query whether Query ID was parsed [in]
 */
void SBlastMessageWrite(SBlastMessage** blast_message, ErrSev sev, const char* message, 
   SeqId* seqid, Boolean believe_query);

/** Frees all memory associated with SBlastMessage*
 * @param message object to be freed [in]
 * @return NULL
 */
SBlastMessage* SBlastMessageFree(SBlastMessage* message);


/** Duplicates entire chain of old messages.
 * @param old object(s) to be duplicated. [in]
 * @return chain of new messages
 */
SBlastMessage* SBlastMessageDup(SBlastMessage* old);

/** Copies messages on (the internal) Blast_Message format 
 * to (C API) SBlastMessage format.
 * @param old object in Blast_Message format to be copied.
 * @param query_slp pointers to queries [in]
 * @param query_info maps context to query number [in]
 * @param believe_query whether Query ID was parsed [in]
 * @return object in SBlastMessage format.
 */
SBlastMessage* Blast_MessageToSBlastMessage(const Blast_Message* old, SeqLoc* query_slp, 
     const BlastQueryInfo* query_info, Boolean believe_query);

/** Posts message using ErrPostEx after prefixing message identifying query.
 * Returns the highest ErrSev status found.
 * @param message Object to be ErrPostEx'd.
 * @return highest severity encountered. 
 */
ErrSev SBlastMessageErrPost(const SBlastMessage* message);

/* @} */

#ifdef __cplusplus
}
#endif

#endif  /* !_BLAST_MESSAGE_API_H_ */
