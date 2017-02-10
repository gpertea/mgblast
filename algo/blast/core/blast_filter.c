/* $Id: blast_filter.c,v 1.82 2006/04/20 19:27:47 madden Exp $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
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
 */

/** @file blast_filter.c
 * All code related to query sequence masking/filtering for BLAST
 */

#ifndef SKIP_DOXYGEN_PROCESSING
static char const rcsid[] = 
    "$Id: blast_filter.c,v 1.82 2006/04/20 19:27:47 madden Exp $";
#endif /* SKIP_DOXYGEN_PROCESSING */

#include <algo/blast/core/blast_def.h>
#include <algo/blast/core/blast_util.h>
#include <algo/blast/core/blast_filter.h>
#include <algo/blast/core/blast_dust.h>
#include <algo/blast/core/blast_seg.h>
#include "blast_inline.h"

/****************************************************************************/
/* Constants */
const Uint1 kNuclMask = 14;     /* N in BLASTNA */
const Uint1 kProtMask = 21;     /* X in NCBISTDAA */


/** Allowed length of the filtering options string. */
#define BLASTOPTIONS_BUFFER_SIZE 128


/** Copies filtering commands for one filtering algorithm from "instructions" to
 * "buffer". 
 * ";" is a delimiter for the commands for different algorithms, so it stops
 * copying when a ";" is found. 
 * Example filtering string: "m L; R -d rodents.lib"
 * @param instructions filtering commands [in] 
 * @param buffer filled with filtering commands for one algorithm. [out]
*/
static const char *
s_LoadOptionsToBuffer(const char *instructions, char* buffer)
{
	Boolean not_started=TRUE;
	char* buffer_ptr;
	const char *ptr;
	Int4 index;

	ptr = instructions;
	buffer_ptr = buffer;
	for (index=0; index<BLASTOPTIONS_BUFFER_SIZE && *ptr != NULLB; index++)
	{
		if (*ptr == ';')
		{	/* ";" is a delimiter for different filtering algorithms. */
			ptr++;
			break;
		}
		/* Remove blanks at the beginning. */
		if (not_started && *ptr == ' ')
		{
			ptr++;
		}
		else
		{
			not_started = FALSE;
			*buffer_ptr = *ptr;
			buffer_ptr++; ptr++;
		}
	}

	*buffer_ptr = NULLB;

	if (not_started == FALSE)
	{	/* Remove trailing blanks. */
		buffer_ptr--;
		while (*buffer_ptr == ' ' && buffer_ptr > buffer)
		{
			*buffer_ptr = NULLB;
			buffer_ptr--;
		}
	}

	return ptr;
}

/** Parses repeat filtering options string.
 * @param repeat_options Input character string [in]
 * @param dbname Database name for repeats filtering [out]
 */
static Int2  
s_ParseRepeatOptions(const char* repeat_options, char** dbname)
{
    char* ptr;

    ASSERT(dbname);
    *dbname = NULL;

    if (!repeat_options)
        return 0;
    
    ptr = strstr(repeat_options, "-d");
    if (ptr) {
        ptr += 2;
        while (*ptr == ' ' || *ptr == '\t')
            ++ptr;
        *dbname = strdup(ptr);
    }
    return 0;
}

/** Parses options used for dust.
 * @param ptr buffer containing instructions. [in]
 * @param level sets level for dust. [out]
 * @param window sets window for dust [out]
 * @param linker sets linker for dust. [out]
*/
static Int2
s_ParseDustOptions(const char *ptr, int* level, int* window, int* linker)

{
	char buffer[BLASTOPTIONS_BUFFER_SIZE];
	int arg, index, index1, window_pri=-1, linker_pri=-1, level_pri=-1;

	arg = 0;
	index1 = 0;
	for (index=0; index<BLASTOPTIONS_BUFFER_SIZE; index++)
	{
		if (*ptr == ' ' || *ptr == NULLB)
		{
	                long tmplong;
			buffer[index1] = NULLB;
			index1 = 0;
			switch(arg) {
				case 0:
					sscanf(buffer, "%ld", &tmplong);
					level_pri = tmplong;
					break;
				case 1:
					sscanf(buffer, "%ld", &tmplong);
					window_pri = tmplong;
					break;
				case 2:
					sscanf(buffer, "%ld", &tmplong);
					linker_pri = tmplong;
					break;
				default:
					break;
			}

			arg++;
			while (*ptr == ' ')
				ptr++;

			/* end of the buffer. */
			if (*ptr == NULLB)
				break;
		}
		else
		{
			buffer[index1] = *ptr; ptr++;
			index1++;
		}
	}
        if (arg != 0 && arg != 3)
           return 1;

	*level = level_pri; 
	*window = window_pri; 
	*linker = linker_pri; 

	return 0;
}

/** parses a string to set three seg options. 
 * @param ptr buffer containing instructions [in]
 * @param window returns "window" for seg algorithm. [out]
 * @param locut returns "locut" for seg. [out]
 * @param hicut returns "hicut" for seg. [out]
*/
static Int2
s_ParseSegOptions(const char *ptr, Int4* window, double* locut, double* hicut)

{
	char buffer[BLASTOPTIONS_BUFFER_SIZE];
	Int4 arg, index, index1; 

	arg = 0;
	index1 = 0;
	for (index=0; index<BLASTOPTIONS_BUFFER_SIZE; index++)
	{
		if (*ptr == ' ' || *ptr == NULLB)
		{
	                long tmplong;
	                double tmpdouble;
			buffer[index1] = NULLB;
			index1 = 0;
			switch(arg) {
				case 0:
					sscanf(buffer, "%ld", &tmplong);
					*window = tmplong;
					break;
				case 1:
					sscanf(buffer, "%le", &tmpdouble);
					*locut = tmpdouble;
					break;
				case 2:
					sscanf(buffer, "%le", &tmpdouble);
					*hicut = tmpdouble;
					break;
				default:
					break;
			}

			arg++;
			while (*ptr == ' ')
				ptr++;

			/* end of the buffer. */
			if (*ptr == NULLB)
				break;
		}
		else
		{
			buffer[index1] = *ptr; ptr++;
			index1++;
		}
	}
        if (arg != 0 && arg != 3)
           return 1;

	return 0;
}



Int2
BlastFilteringOptionsFromString(EBlastProgramType program_number, const char* instructions, 
       SBlastFilterOptions* *filtering_options, Blast_Message* *blast_message)
{
        Boolean mask_at_hash = FALSE; /* the default. */
        char* buffer;
        const char* ptr = instructions;
        char error_buffer[1024];
        Int2 status = 0;
        SSegOptions* segOptions = NULL;
        SDustOptions* dustOptions = NULL;
        SRepeatFilterOptions* repeatOptions = NULL;
   
        *filtering_options = NULL;
        if (blast_message)
            *blast_message = NULL;

        if (instructions == NULL || strcasecmp(instructions, "F") == 0)
        {
             SBlastFilterOptionsNew(filtering_options, eEmpty);
             return status;
        }

        buffer = (char*) calloc(strlen(instructions), sizeof(char));
	/* allow old-style filters when m cannot be followed by the ';' */
	if (ptr[0] == 'm' && ptr[1] == ' ')
	{
		mask_at_hash = TRUE;
		ptr += 2;
	}

	while (*ptr != NULLB)
	{
		if (*ptr == 'S')
		{
                        SSegOptionsNew(&segOptions);
			ptr = s_LoadOptionsToBuffer(ptr+1, buffer);
			if (buffer[0] != NULLB)
			{
                                int window;
                                double locut, hicut; 
				status = s_ParseSegOptions(buffer, &window, &locut, &hicut);
                                if (status)
                                {
                                     segOptions = SSegOptionsFree(segOptions);
                                     sprintf(error_buffer, "Error parsing filter string: %s", buffer);
                                     if (blast_message)
                                       Blast_MessageWrite(blast_message, eBlastSevError, kBlastMessageNoContext,
                                            error_buffer);
                                     sfree(buffer);
                                     return status;
                                }
                                segOptions->window = window;
                                segOptions->locut = locut;
                                segOptions->hicut = hicut;
			}
		}
		else if (*ptr == 'D')
		{
                        SDustOptionsNew(&dustOptions);
			ptr = s_LoadOptionsToBuffer(ptr+1, buffer);
			if (buffer[0] != NULLB)
                        {
                                int window, level, linker;
				status = s_ParseDustOptions(buffer, &level, &window, &linker);
                                if (status)
                                {
                                     dustOptions = SDustOptionsFree(dustOptions);
                                     sprintf(error_buffer, "Error parsing filter string: %s", buffer);
                                     if (blast_message)
                                       Blast_MessageWrite(blast_message, eBlastSevError, kBlastMessageNoContext,
                                            error_buffer);
                                     sfree(buffer);
                                     return status;
                                }
                                dustOptions->level = level;
                                dustOptions->window = window;
                                dustOptions->linker = linker;
                        }
		}
                else if (*ptr == 'R')
                {
                        SRepeatFilterOptionsNew(&repeatOptions);
			ptr = s_LoadOptionsToBuffer(ptr+1, buffer);
			if (buffer[0] != NULLB)
                        {
                             char* dbname = NULL;
                             status = s_ParseRepeatOptions(buffer, &dbname); 
                             if (status)
                             {
                                  repeatOptions = SRepeatFilterOptionsFree(repeatOptions);
                                  sprintf(error_buffer, "Error parsing filter string: %s", buffer);
                                  if (blast_message)
                                     Blast_MessageWrite(blast_message, eBlastSevError, kBlastMessageNoContext,
                                            error_buffer);
                                   sfree(buffer);
                                   return status;
                             }
                             if (dbname)
                             {
                                 sfree(repeatOptions->database);
                                 repeatOptions->database = dbname;
                             }
                        }
                }
		else if (*ptr == 'L' || *ptr == 'T')
		{ /* do low-complexity filtering; dust for blastn, otherwise seg.*/
                        if (program_number == eBlastTypeBlastn)
                            SDustOptionsNew(&dustOptions);
                        else
                            SSegOptionsNew(&segOptions);
			ptr++;
		}
		else if (*ptr == 'm')
		{
		        mask_at_hash = TRUE;
                        ptr++;
                }
                else
                {    /* Nothing applied */
                         ptr++;
                }
	}
        sfree(buffer);

        status = SBlastFilterOptionsNew(filtering_options, eEmpty);
        if (status)
            return status;

        (*filtering_options)->dustOptions = dustOptions;
        (*filtering_options)->segOptions = segOptions;
        (*filtering_options)->repeatFilterOptions = repeatOptions;
        (*filtering_options)->mask_at_hash = mask_at_hash;

        return status;
}


BlastSeqLoc* BlastSeqLocNew(BlastSeqLoc** head, Int4 from, Int4 to)
{
   BlastSeqLoc* loc = (BlastSeqLoc*) calloc(1, sizeof(BlastSeqLoc));
   if ( !loc ) {
       return NULL;
   }
   loc->ssr = (SSeqRange*) calloc(1, sizeof(SSeqRange));
   loc->ssr->left = from;
   loc->ssr->right = to;

   return BlastSeqLocAppend(head, loc);
}

BlastSeqLoc* BlastSeqLocAppend(BlastSeqLoc** head, BlastSeqLoc* node)
{
    if ( !node ) {
        return NULL;
    }

    if (head)
    {
        if (*head)
        {
            BlastSeqLoc* tmp = *head;
            while (tmp->next)
               tmp = tmp->next;
            tmp->next = node;
        }
        else
        {
            *head = node;
        }
    }
        
    return node;
}

/** Makes a copy of the BlastSeqLoc and also a copy of the 
 * SSRange element.  Does not copy BlastSeqLoc that is pointed
 * to by "next".
 * @param source the object to be copied [in]
 * @return another BlastSeqLoc*
 */
static BlastSeqLoc* s_BlastSeqLocNodeDup(BlastSeqLoc* source)
{
    if ( !source ) {
        return NULL;
    }
    ASSERT(source->ssr);
    return BlastSeqLocNew(NULL, source->ssr->left, source->ssr->right);
}

/** Prepend node to the head of the list and return the new head of the list */
static BlastSeqLoc* s_BlastSeqLocPrepend(BlastSeqLoc* head, BlastSeqLoc* node)
{
    if ( !node ) {
        return NULL;
    }
    node->next = head;
    return node;
}

/** Reverse elements in the list 
 * @param head pointer to pointer to the head of the list. After this call,
 * this is set to NULL [in|out]
 * @return the new head of the list or NULL if argument is NULL
 */
static BlastSeqLoc* s_BlastSeqLocListReverse(BlastSeqLoc** head)
{
    BlastSeqLoc* retval = NULL;     /* return value */
    BlastSeqLoc* itr = NULL;        /* iterator */

    if ( !head ) {
        return NULL;
    }

    for (itr = *head; itr; itr = itr->next) {
        retval = s_BlastSeqLocPrepend(retval, s_BlastSeqLocNodeDup(itr));
    }
    *head = BlastSeqLocFree(*head);
    return retval;
}

BlastSeqLoc* BlastSeqLocNodeFree(BlastSeqLoc* loc)
{
    if ( !loc ) {
        return NULL;
    }
    sfree(loc->ssr);
    sfree(loc);
    return NULL;
}

BlastSeqLoc* BlastSeqLocFree(BlastSeqLoc* loc)
{
    while (loc) {
        BlastSeqLoc* next_loc = loc->next;
        loc = BlastSeqLocNodeFree(loc);
        loc = next_loc;
    }
    return NULL;
}

BlastSeqLoc* BlastSeqLocListDup(BlastSeqLoc* head)
{
    BlastSeqLoc* retval = NULL;
    BlastSeqLoc* retval_tail = NULL;

    for (; head; head = head->next) {
        retval_tail = BlastSeqLocAppend(retval_tail ? &retval_tail : &retval, 
                                        s_BlastSeqLocNodeDup(head));
    }

    return retval;
}

BlastMaskLoc* BlastMaskLocNew(Int4 total)
{
    BlastMaskLoc* retval = (BlastMaskLoc *) calloc(1, sizeof(BlastMaskLoc));
    retval->total_size = total;
    if (total > 0)
        retval->seqloc_array = (BlastSeqLoc **) calloc(total, 
                                                       sizeof(BlastSeqLoc *));
    return retval;
}

BlastMaskLoc* BlastMaskLocFree(BlastMaskLoc* mask_loc)
{
   Int4 index;

   if (mask_loc == NULL)
      return NULL;

   for (index=0; index<mask_loc->total_size; index++)
   {
      if (mask_loc->seqloc_array != NULL)
         BlastSeqLocFree(mask_loc->seqloc_array[index]);
   }
   sfree(mask_loc->seqloc_array);
   sfree(mask_loc);
   return NULL;
}

Int2 BlastMaskLocDNAToProtein(BlastMaskLoc* mask_loc, 
                              const BlastQueryInfo* query_info)
{
    Uint4 seq_index;

    if (!mask_loc)
        return 0;

    /* Check that the array size in BlastMaskLoc corresponds to the number
       of contexts in BlastQueryInfo. */
    ASSERT(mask_loc->total_size == query_info->last_context + 1);

    /* Loop over multiple DNA sequences */
    for (seq_index = 0; seq_index < (Uint4)query_info->num_queries; 
         ++seq_index) { 
        const Uint4 kCtxIndex = NUM_FRAMES * seq_index;
        BlastSeqLoc* dna_seqloc = mask_loc->seqloc_array[kCtxIndex];
        BlastSeqLoc** prot_seqloc = &(mask_loc->seqloc_array[kCtxIndex]);
        Int4 dna_length = BlastQueryInfoGetQueryLength(query_info,
                                                       eBlastTypeBlastx,
                                                       seq_index);
        Int4 context;

        /* Reproduce this mask for all 6 frames, with translated coordinates */
        for (context = 0; context < NUM_FRAMES; ++context) {
            BlastSeqLoc* prot_head=NULL;
            BlastSeqLoc* seqloc_var;
            Int2 frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);

            for (seqloc_var = dna_seqloc; seqloc_var; 
                 seqloc_var = seqloc_var->next) {
                Int4 from, to;
                SSeqRange* seq_range = seqloc_var->ssr;
                if (frame < 0) {
                    from = (dna_length + frame - seq_range->right)/CODON_LENGTH;
                    to = (dna_length + frame - seq_range->left)/CODON_LENGTH;
                } else {
                    from = (seq_range->left - frame + 1)/CODON_LENGTH;
                    to = (seq_range->right - frame + 1)/CODON_LENGTH;
                }
                BlastSeqLocNew(&prot_head, from, to);
            }
            prot_seqloc[context] = prot_head;
        }
        BlastSeqLocFree(dna_seqloc);
    }

    return 0;
}


Int2 BlastMaskLocProteinToDNA(BlastMaskLoc* mask_loc, 
                              const BlastQueryInfo* query_info)
{
   Int2 status = 0;
   Int4 index;

   /* If there is not mask, there is nothing to convert to DNA coordinates,
      hence just return. */
   if (!mask_loc) 
      return 0;

   /* Check that the array size in BlastMaskLoc corresponds to the number
      of contexts in BlastQueryInfo. */
   ASSERT(mask_loc->total_size == query_info->last_context + 1);

   /* Loop over all DNA sequences */
   for (index=0; index < query_info->num_queries; ++index)
   {
       Int4 frame_start = index*NUM_FRAMES;
       Int4 frame_index;
       Int4 dna_length = BlastQueryInfoGetQueryLength(query_info,
                                                      eBlastTypeBlastx,
                                                      index);
       /* Loop over all frames of one DNA sequence */
       for (frame_index=frame_start; frame_index<(frame_start+NUM_FRAMES); 
            frame_index++) {
           BlastSeqLoc* loc;
           Int2 frame = 
               BLAST_ContextToFrame(eBlastTypeBlastx, frame_index % NUM_FRAMES);
           /* Loop over all mask locations for a given frame */
           for (loc = mask_loc->seqloc_array[frame_index]; loc; loc = loc->next) {
               Int4 from=0, to=0;
               SSeqRange* seq_range = loc->ssr;
               if (frame < 0) {
                   to = dna_length - CODON_LENGTH*seq_range->left + frame;
                   from = dna_length - CODON_LENGTH*seq_range->right + frame + 1;
               } else {
                   from = CODON_LENGTH*seq_range->left + frame - 1;
                   to = CODON_LENGTH*seq_range->right + frame - 1;
               }
               seq_range->left = from;
               seq_range->right = to;
           }
       }
   }
   return status;
}

/** Used for qsort, compares two SeqLoc's by starting position. */
static int s_SeqRangeSortByStartPosition(const void *vp1, const void *vp2)
{
   BlastSeqLoc* v1 = *((BlastSeqLoc**) vp1);
   BlastSeqLoc* v2 = *((BlastSeqLoc**) vp2);
   SSeqRange* loc1 = (SSeqRange*) v1->ssr;
   SSeqRange* loc2 = (SSeqRange*) v2->ssr;
   
   if (loc1->left < loc2->left)
      return -1;
   else if (loc1->left > loc2->left)
      return 1;
   else
      return 0;
}

/** Calculates number of links in a chain of BlastSeqLoc's.
 * @param var Chain of BlastSeqLoc structures [in]
 * @return Number of links in the chain.
 */
static Int4 s_BlastSeqLocLen(BlastSeqLoc* var)
{
     Int4 count=0;
   
     while (var)
     {
         count++;
         var = var->next;
     }

     return count;
}

/** Sort a list of BlastSeqLoc structures
 * @param list List of BlastSeqLoc's [in]
 * @param compar Comparison function to use [in]
 * @return Sorted list.  
 */
static BlastSeqLoc* 
s_BlastSeqLocSort (BlastSeqLoc* list,
                 int (*compar )(const void *, const void *))
{
        BlastSeqLoc* tmp,** head;
        Int4 count, i;

        if (list == NULL) 
           return NULL;

        count = s_BlastSeqLocLen (list);
        head = (BlastSeqLoc* *) calloc (((size_t) count + 1), sizeof (BlastSeqLoc*));
        for (tmp = list, i = 0; tmp != NULL && i < count; i++) {
                head [i] = tmp;
                tmp = tmp->next;
        }

        qsort (head, (size_t) count, sizeof (BlastSeqLoc*), compar);
        for (i = 0; i < count; i++) {
                tmp = head [i];
                tmp->next = head [i + 1];
        }
        list = head [0];
        sfree (head);

        return list;
}

BlastSeqLoc*
BlastSeqLocCombine(BlastSeqLoc* mask_loc, Int4 link_value)
{
   BlastSeqLoc* retval = NULL;
   BlastSeqLoc* retval_tail = NULL;
   Int4 start, stop;	/* Used to merge overlapping SeqLoc's. */
   BlastSeqLoc* loc_head=NULL,* loc_var=NULL;
   
   if ( !mask_loc ) {
      return NULL;
   }

   /* Copy the BlastSeqLoc-s and sort them by starting position. */
   loc_head = loc_var = s_BlastSeqLocSort(BlastSeqLocListDup(mask_loc),
                                          s_SeqRangeSortByStartPosition);
   if ( !loc_head ) {
       return NULL;
   }
   start = loc_head->ssr->left;
   stop = loc_head->ssr->right;

   for (; loc_var; loc_var = loc_var->next) {
       SSeqRange* ssr = loc_var->next ? loc_var->next->ssr : NULL;

       if (ssr && ((stop + link_value) > ssr->left)) {
          stop = MAX(stop, ssr->right);
       } else {
          /* Cache the tail of the list to avoid the overhead of traversing the
           * list when appending to it */
          retval_tail = BlastSeqLocNew((retval_tail ? &retval_tail : &retval), 
                                       start, stop);
          if (loc_var->next) {
              start = ssr->left;
              stop = ssr->right;
          }
       }
   }

   /* Free memory allocated for the temporary list of BlastSeqLoc-s */
   BlastSeqLocFree(loc_head);

   return retval;
}

Int2 
BLAST_ComplementMaskLocations(EBlastProgramType program_number, 
   const BlastQueryInfo* query_info, 
   const BlastMaskLoc* mask_loc, BlastSeqLoc* *complement_mask) 
{
   Int4 context;
   const Boolean kIsNucl = (program_number == eBlastTypeBlastn);
   BlastSeqLoc* tail = NULL;    /* Pointer to the tail of the complement_mask
                                   linked list */

   if (complement_mask == NULL)
	return -1;

   *complement_mask = NULL;

   for (context = query_info->first_context; 
        context <= query_info->last_context; ++context) {

      Boolean first = TRUE;	/* Specifies beginning of query. */
      Boolean last_interval_open=TRUE; /* if TRUE last interval needs to be closed. */
      Int4 start_offset, end_offset, filter_start, filter_end;
      Int4 left=0, right; /* Used for left/right extent of a region. */
      BlastSeqLoc* loc = NULL;

      if (query_info->contexts[context].is_valid == FALSE) {
          continue;
      }

      start_offset = query_info->contexts[context].query_offset;
      end_offset = query_info->contexts[context].query_length 
          + start_offset - 1;
      ASSERT(start_offset <= end_offset);

      /* mask_loc NULL is simply the case that NULL was passed in, which we 
         take to mean that nothing on query is masked. */
      if (mask_loc == NULL || mask_loc->seqloc_array[context] == NULL) {
         /* Cache the tail of the list to avoid the overhead of traversing the
          * list when appending to it */
         tail = BlastSeqLocNew(tail ? &tail : complement_mask, 
                               start_offset, end_offset);
         continue;
      }
      
      if (BlastIsReverseStrand(kIsNucl, context)) {
         mask_loc->seqloc_array[context] = 
             s_BlastSeqLocListReverse(&mask_loc->seqloc_array[context]);
      }
      loc = mask_loc->seqloc_array[context];

      first = TRUE;
      for ( ; loc; loc = loc->next) {
         SSeqRange* seq_range = loc->ssr;
         if (BlastIsReverseStrand(kIsNucl, context)) {
            filter_start = end_offset - seq_range->right;
            filter_end = end_offset - seq_range->left;
         } else {
            filter_start = start_offset + seq_range->left;
            filter_end = start_offset + seq_range->right;
         }
         /* The canonical "state" at the top of this 
            while loop is that both "left" and "right" have
            been initialized to their correct values.  
            The first time this loop is entered in a call to 
            the function this is not true and the following "if" 
            statement moves everything to the canonical state. */
         if (first) {
            last_interval_open = TRUE;
            first = FALSE;
            
            if (filter_start > start_offset) {
               /* beginning of sequence not filtered */
               left = start_offset;
            } else {
               /* beginning of sequence filtered */
               left = filter_end + 1;
               continue;
            }
         }

         right = filter_start - 1;

         /* Cache the tail of the list to avoid the overhead of traversing the
          * list when appending to it */
         tail = BlastSeqLocNew((tail ? &tail : complement_mask), left, right);
         if (filter_end >= end_offset) {
            /* last masked region at end of sequence */
            last_interval_open = FALSE;
            break;
         } else {
            left = filter_end + 1;
         }
      }

      if (last_interval_open) {
         /* Need to finish SSeqRange* for last interval. */
         right = end_offset;
         /* Cache the tail of the list to avoid the overhead of traversing the
          * list when appending to it */
         tail = BlastSeqLocNew((tail ? &tail : complement_mask), left, right);
      }
   }
   return 0;
}


Int2
BlastSetUp_Filter(EBlastProgramType program_number, 
                  Uint1* sequence, 
                  Int4 length, 
                  Int4 offset, 
                  const SBlastFilterOptions* filter_options, 
                  BlastSeqLoc** seqloc_retval, 
                  Blast_Message* *blast_message)
{
	Int2 status=0;		/* return value. */

    ASSERT(filter_options);
    ASSERT(seqloc_retval);

    *seqloc_retval = NULL;

    status = SBlastFilterOptionsValidate(program_number, filter_options, 
                                         blast_message);
    if (status)
       return status;

	if (filter_options->segOptions)
	{
        SSegOptions* seg_options = filter_options->segOptions;
        SegParameters* sparamsp=NULL;

        sparamsp = SegParametersNewAa();
        sparamsp->overlaps = TRUE;
        if (seg_options->window > 0)
            sparamsp->window = seg_options->window;
        if (seg_options->locut > 0.0)
            sparamsp->locut = seg_options->locut;
        if (seg_options->hicut > 0.0)
            sparamsp->hicut = seg_options->hicut;

		status = SeqBufferSeg(sequence, length, offset, sparamsp, 
                              seqloc_retval);
		SegParametersFree(sparamsp);
		sparamsp = NULL;
	}

	return status;
}

BlastSeqLoc*
BlastSeqLocReverse(const BlastSeqLoc* filter_in, Int4 query_length)
{
   BlastSeqLoc* filter_out = NULL;   /* Return variable. */

   while (filter_in)
   {
        Int4 start, stop;
        SSeqRange* loc = filter_in->ssr;
        
        start = query_length - 1 - loc->right;
        stop = query_length - 1 - loc->left;
        BlastSeqLocNew(&filter_out, start, stop);
        filter_in = filter_in->next;
   }

   return filter_out;
}

/** Calculates the mask locations one context at a time.
 * @param query_blk sequence [in]
 * @param query_info information about sequences [in]
 * @param context which context is this?  [in]
 * @param program_number program (blastn, blastp, etc.) [in]
 * @param filter_options instructions for producing mask [in]
 * @param filter_out results of filtering operations [out]
 * @param blast_message any error or warning messages [out]
 * @return zero on success
 */
static Int2
s_GetFilteringLocationsForOneContext(BLAST_SequenceBlk* query_blk, 
                                     const BlastQueryInfo* query_info, 
                                     Int4 context, 
                                     EBlastProgramType program_number, 
                                     const SBlastFilterOptions* filter_options, 
                                     BlastSeqLoc* *filter_out, 
                                     Blast_Message* *blast_message)
{
    Int2 status = 0;
    Int4 query_length = 0;      /* Length of query described by SeqLocPtr. */
    Int4 context_offset;
    BlastSeqLoc *filter_slp = NULL;     /* SeqLocPtr computed for filtering. */
    Uint1 *buffer;              /* holds sequence for plus strand or protein. */

    const Boolean kIsNucl = (program_number == eBlastTypeBlastn);

    context_offset = query_info->contexts[context].query_offset;
    buffer = &query_blk->sequence[context_offset];

    if (query_info->contexts[context].is_valid == FALSE) {
          return 0;
    }

    query_length = query_info->contexts[context].query_length;

    status = BlastSetUp_Filter(program_number, 
                               buffer, 
                               query_length, 
                               0, 
                               filter_options, 
                               &filter_slp, 
                               blast_message);
    if (status)
         return status;

    if (BlastIsReverseStrand(kIsNucl, context) == TRUE)
    {  /* Reverse this as it's on minus strand. */
          BlastSeqLoc* tmp = BlastSeqLocReverse(filter_slp, query_length);
          filter_slp = BlastSeqLocFree(filter_slp);
          filter_slp = tmp;
    }

    /* Extract the mask locations corresponding to this query 
       (frame, strand), detach it from other masks.
       NB: for translated search the mask locations are expected in 
       protein coordinates. The nucleotide locations must be converted
       to protein coordinates prior to the call to BLAST_MainSetUp.
    */
    {
        /* Auxiliary locations for lower-case masking or any other masking
         * which occurred outside of CORE BLAST */
        BlastSeqLoc *lcase_mask_slp = NULL; 
        if (query_blk->lcase_mask && query_blk->lcase_mask->seqloc_array)
        {
            ASSERT(context < query_blk->lcase_mask->total_size);
            lcase_mask_slp = query_blk->lcase_mask->seqloc_array[context];
            /* Set location list to NULL, to allow safe memory deallocation, 
              ownership transferred to filter_slp below. */
            query_blk->lcase_mask->seqloc_array[context] = NULL;
        }

        /* Attach the lower case mask locations to the filter locations and 
           combine them */
        BlastSeqLocAppend(&filter_slp, lcase_mask_slp);
    }

    *filter_out = BlastSeqLocCombine(filter_slp, 0);
    filter_slp = BlastSeqLocFree(filter_slp);

	return 0;
}

Int2
BlastSetUp_GetFilteringLocations(BLAST_SequenceBlk* query_blk, 
                                 const BlastQueryInfo* query_info, 
                                 EBlastProgramType program_number, 
                                 const SBlastFilterOptions* filter_options, 
                                 BlastMaskLoc** filter_maskloc, 
                                 Blast_Message** blast_message)
{
    Int2 status = 0;
    Int4 context = 0; /* loop variable. */
    const int kNumContexts = query_info->last_context + 1;

    ASSERT(query_info && query_blk && filter_maskloc);

    ASSERT(blast_message);
    ASSERT(kNumContexts == 
           query_info->num_queries*BLAST_GetNumberOfContexts(program_number));
    *filter_maskloc = BlastMaskLocNew(kNumContexts);

    for (context = query_info->first_context;
         context <= query_info->last_context; ++context) {
  
        BlastSeqLoc *filter_per_context = NULL;
        status = s_GetFilteringLocationsForOneContext(query_blk, 
                                                      query_info, 
                                                      context, 
                                                      program_number, 
                                                      filter_options, 
                                                      &filter_per_context, 
                                                      blast_message);
        if (status) {
            Blast_MessageWrite(blast_message, eBlastSevError, context,
                                   "Failure at filtering");
            return status;
        }

    /* NB: for translated searches filter locations are returned in 
           protein coordinates, because the DNA lengths of sequences are 
           not available here. The caller must take care of converting 
           them back to nucleotide coordinates. */
         (*filter_maskloc)->seqloc_array[context] = filter_per_context;
    }
    return 0;
}

void
Blast_MaskTheResidues(Uint1 * buffer, Int4 length, Boolean is_na,
                      const BlastSeqLoc* mask_loc, Boolean reverse, Int4 offset)
{
    ASSERT(buffer);
    for (; mask_loc; mask_loc = mask_loc->next) {

        Int4 index, start, stop;
        const Uint1 kMaskingLetter = is_na ? kNuclMask : kProtMask;
        
        if (reverse) {
            start = length - 1 - mask_loc->ssr->right;
            stop = length - 1 - mask_loc->ssr->left;
        } else {
            start = mask_loc->ssr->left;
            stop = mask_loc->ssr->right;
        }
        
        start -= offset;
        stop -= offset;
        
        ASSERT(start < length);
        ASSERT(stop < length);
        
        for (index = start; index <= stop; index++)
            buffer[index] = kMaskingLetter;
    }
}

void
BlastSetUp_MaskQuery(BLAST_SequenceBlk* query_blk, 
                     const BlastQueryInfo* query_info, 
                     const BlastMaskLoc *filter_maskloc, 
                     EBlastProgramType program_number)
{
    const Boolean kIsNucl = (program_number == eBlastTypeBlastn);
    Int4 context; /* loop variable. */

    ASSERT(query_blk);
    ASSERT(query_info);
    ASSERT(filter_maskloc);

    for (context = query_info->first_context;
         context <= query_info->last_context; ++context) {
      
        Int4 query_length = 0;
        Int4 context_offset = 0;
        Uint1 *buffer = NULL;              /* holds sequence */

        if (query_info->contexts[context].is_valid == FALSE) {
          continue;
        }

        query_length = query_info->contexts[context].query_length;

        context_offset = query_info->contexts[context].query_offset;
        buffer = &query_blk->sequence[context_offset];
        ASSERT(buffer);

        Blast_MaskTheResidues(buffer, query_length, kIsNucl, 
                              filter_maskloc->seqloc_array[context],
                              BlastIsReverseStrand(kIsNucl, context), 0);
    }
}
