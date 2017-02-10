/* $Id: blast_message_api.c,v 1.2 2006/04/26 14:11:33 kans Exp $
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
*                                                                         *
*  Author: Tom Madden                                                     *
**************************************************************************/

/** @file blast_message_api.c
 * Produces messages for the C API, using the internal structure from core. 
 */

#include <sequtil.h>
#include <blfmtutl.h>
#include <algo/blast/api/blast_message_api.h>
#include <algo/blast/core/blast_message.h>

/** @addtogroup CToolkitAlgoBlast
 *
 * @{
 */

SBlastMessage* SBlastMessageFree(SBlastMessage* message)
{
    SBlastMessage* var = message;
    SBlastMessage* next = NULL;

    while (var)
    {
         next = var->next;
         sfree((var->message));
         SeqIdSetFree(var->query_id);
         var = next;
    }

    return NULL;
}

void SBlastMessageWrite(SBlastMessage** blast_message, ErrSev sev, const char* message, SeqId* seqid, Boolean believe_query)
{
     SBlastMessage* new_message = (SBlastMessage*) calloc(1, sizeof(SBlastMessage));
     ASSERT(blast_message);

     new_message->sev = sev;
     new_message->message = strdup(message);
     if (seqid)
     {
        new_message->query_id = SeqIdDup(seqid);
        new_message->believe_query = believe_query;
     }

     if (*blast_message)
     {
     	SBlastMessage* var = *blast_message;
     	while (var->next)
     	   var = var->next;
        var->next = new_message;
     }
     else
     {
         *blast_message = new_message;
     }

     return;
}

SBlastMessage* SBlastMessageDup(SBlastMessage* old)
{
    SBlastMessage* retval = NULL;
    
    while (old)
    {
	SBlastMessageWrite(&retval, old->sev, old->message, old->query_id, old->believe_query);
        old = old->next;
    }
    return retval;
}


ErrSev s_EBlastSeverity2ErrSev(EBlastSeverity severity)
{
    ErrSev retval = SEV_NONE;

    switch (severity)
    {
         case eBlastSevInfo:
           retval = SEV_INFO;
           break;
         case eBlastSevWarning:
           retval = SEV_WARNING;
           break;
         case eBlastSevError:
           retval = SEV_ERROR;
           break;
         case eBlastSevFatal:
           retval = SEV_FATAL;
           break;
    }

    return retval;
   
}

SeqLoc* s_Context2SeqLoc(int context, SeqLoc* query_slp, const BlastQueryInfo* query_info)
{
      SeqLoc* slp = NULL;

      if (context != -1 && query_slp != NULL && query_info != NULL)
      {
          int query_num = query_info->contexts[context].query_index;
          int index = 0;
          while (query_slp)
          {
               if (index == query_num)
               {
                  slp = query_slp;
                  break;
               }
               query_slp = query_slp->next;
               index++;
          }

         
      }
      return slp;
}

/** Compares two SBlastMessage elements and returns TRUE if they are identical.
 * FALSE is returned if they are not or one or both is NULL.
 * @param message_1 first object to compare [in]
 * @param message_2 second object to compare [in]
 * @return TRUE if identical
 */
Boolean s_SBlastMessageCompare(SBlastMessage* message_1, SBlastMessage* message_2)
{
    if (message_1 == NULL || message_2 == NULL)
        return FALSE;

    if (message_1->sev != message_2->sev)
       return FALSE;

    if (message_1->message == NULL || message_2->message == NULL)
       return FALSE;

    if (StringCmp(message_1->message, message_2->message))
       return FALSE;

    if ((message_1->query_id == NULL && message_2->query_id) ||
        (message_1->query_id && message_2->query_id == NULL))
       return FALSE;

    if (message_1->query_id && message_2->query_id && 
        SeqIdComp(message_1->query_id, message_2->query_id) != SIC_YES)
       return FALSE;

    if (message_1->believe_query != message_2->believe_query)
       return FALSE;

    return TRUE;
}

SBlastMessage* Blast_MessageToSBlastMessage(const Blast_Message* old, SeqLoc* query_slp, const BlastQueryInfo* query_info, Boolean believe_query)
{
    SBlastMessage* retval = NULL;
    SBlastMessage* last = NULL;
    SBlastMessage* var = NULL;

    while (old)
    {
         SeqId* sip = NULL;
         SeqLoc* slp = s_Context2SeqLoc(old->context, query_slp, query_info);

         if (slp)
           sip = SeqLocId(slp);
             
         SBlastMessageWrite(&retval, s_EBlastSeverity2ErrSev(old->severity), old->message, sip, believe_query);
         old = old->next;
    }

    /* 
    This loop removes a redundant message only if it is the next message.  It removes the second one so that
    the head of the list is not changed. 
    */
    var = retval;
    while (var)
    {
       if (s_SBlastMessageCompare(var, last) == TRUE)
       {
          last->next = var->next;
          var->next = NULL;
          SBlastMessageFree(var);
          var = last->next;
       }
       else
       {
          last = var;
          var = var->next;
       }
    }
    
    return retval;
}

ErrSev SBlastMessageErrPost(const SBlastMessage* message)
{
        ErrSev max_sev = SEV_NONE;
	Uint1 err_id = 0;
        
        while (message)
        {
           err_id = BlastSetUserErrorString(NULL, message->query_id, message->believe_query);
           ErrPostEx(message->sev, 0, 0, message->message); 
           BlastDeleteUserErrorString(err_id);
           err_id = 0;
           max_sev = MAX(max_sev, message->sev);
           message = message->next;
        }
        return max_sev;
}



/* @} */

