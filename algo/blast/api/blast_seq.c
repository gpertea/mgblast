static char const rcsid[] = "$Id: blast_seq.c,v 1.82 2006/04/20 19:31:44 madden Exp $";
/*
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's offical duties as a United States Government employee and
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
*  Author: Ilya Dondoshansky
* ===========================================================================*/

/** @file blast_seq.c
 * Functions converting between SeqLocs and structures used in BLAST.
 */

#include <seqport.h>
#include <sequtil.h>
#include <objloc.h>
#include <readdb.h>
#include <algo/blast/api/blast_seq.h>
#include <algo/blast/core/blast_filter.h>
#include <algo/blast/core/blast_util.h>
#include <algo/blast/core/blast_encoding.h>
#include <algo/blast/core/blast_setup.h> /* For BlastSeqLoc_RestrictToInterval */
#include <algo/blast/core/blast_inline.h>

/** @addtogroup CToolkitAlgoBlast
 *
 * @{
 */


/** Structure used for hash-based comparison of sequence IDs */
typedef struct SeqIdHash {
    SeqId *id;           /**< The ID of this entry */
    Int4 query_index;    /**< index of query with this ID */
    Int4 next_id;        /**< Offset of the next hash entry in a chain */
} SeqIdHash;

Boolean
BlastSeqlocsHaveDuplicateIDs(SeqLoc* query_seqlocs)
{
   Boolean retval = FALSE;
   const Int4 kNumSeqs = ValNodeLen(query_seqlocs);
   const Int4 kLog2HashSize = 11;
   SeqIdHash *id_entries;
   Uint4 *hashtable;
   Int4 curr_id_num;
   SeqLocPtr slp;

   if (kNumSeqs == 1)
      return FALSE;

   /* allocate hashtable */
   hashtable = (Uint4 *)calloc(1 << kLog2HashSize, sizeof(Uint4));
   id_entries = (SeqIdHash *)malloc((kNumSeqs + 1) * sizeof(SeqIdHash));

   for (slp = query_seqlocs, curr_id_num = 1; slp; slp = slp->next) {

       Uint4 hashval;
       SeqIdPtr id = SeqLocId(slp);
       Char buffer[64];

       /* hash the ID of the next query sequence */
       SeqIdLabel(id, buffer, sizeof(buffer), OM_LABEL_CONTENT);
       hashval = readdb_sequence_hash(buffer, strlen(buffer));
       hashval = hashval >> (32 - kLog2HashSize);
       if (hashtable[hashval] != 0) {
          Int4 offset = hashtable[hashval];
          /* check the chain of ID's in the resulting hash
             entry for a duplicate */
          while (offset != 0) {
             SeqIdHash *curr_entry = id_entries + offset;
             if (SeqIdMatch(id, curr_entry->id)) {
                 retval = TRUE;
                 goto clean_up;
             }
             offset = curr_entry->next_id;
          }
       }

       /* no duplicate found; add the ID to the hashtable */
       id_entries[curr_id_num].id = id;
       id_entries[curr_id_num].next_id = hashtable[hashval];
       hashtable[hashval] = curr_id_num++;
   }

clean_up:
   sfree(hashtable);
   sfree(id_entries);
   return retval;
}

/** Converts a SeqLocPtr to a BlastSeqLoc, used for formatting.
 * @param slp SeqLocPtr to be converted [in]
 * @param head_loc BlastSeqLoc returned from last call [in]
 * @return pointer to BlastSeqLoc
 */
static BlastSeqLoc* 
s_BlastSeqLocFromSeqLoc(SeqLocPtr slp, BlastSeqLoc* head_loc)
{
   BlastSeqLoc* last_loc = head_loc;

   if (slp == NULL)
      return NULL;

   if (slp->choice == SEQLOC_PACKED_INT)
      slp = (SeqLocPtr) slp->data.ptrvalue;

   for ( ; slp; slp = slp->next) {
      SeqIntPtr si = (SeqIntPtr) slp->data.ptrvalue;
      if (!head_loc) {
         last_loc = head_loc = BlastSeqLocNew(&last_loc, si->from, si->to);
      } else {
         last_loc = BlastSeqLocNew(&last_loc, si->from, si->to);
      }
   }
   return head_loc;
}

BlastMaskLoc* 
BlastMaskLocFromSeqLoc(SeqLoc* mask_seqlocs, SeqLoc* query_seqlocs, 
                       EBlastProgramType program_number)
{
    const Int4 kNumSeqs = ValNodeLen(query_seqlocs);
    BlastMaskLoc* retval = NULL;
    Int4 query_index = 0;
    const unsigned int kNumContexts = BLAST_GetNumberOfContexts(program_number);
    const Int4 kLog2HashSize = 11;
    SeqIdHash *id_entries;
    Uint4 *hashtable;
    Int4 curr_id_num = 1;
    SeqLocPtr query_slp, mask_slp;

    if (!mask_seqlocs)
        return NULL;

    retval = BlastMaskLocNew(kNumSeqs*kNumContexts);

    /* create hashtable for query IDs */
    hashtable = (Uint4 *)calloc(1 << kLog2HashSize, sizeof(Uint4));
    id_entries = (SeqIdHash *)malloc((kNumSeqs + 1) * sizeof(SeqIdHash));
 
    /* add the ID of each query sequence to the hashtable */
    for (query_slp = query_seqlocs; query_slp; query_slp = query_slp->next) {
        Uint4 hashval;
        SeqIdPtr seq_id = SeqLocId(query_slp);
        Char buffer[64];
 
        SeqIdLabel(seq_id, buffer, sizeof(buffer), OM_LABEL_CONTENT);
        hashval = readdb_sequence_hash(buffer, strlen(buffer));
        hashval = hashval >> (32 - kLog2HashSize);
 
        id_entries[curr_id_num].id = seq_id;
        id_entries[curr_id_num].query_index = query_index++;
        id_entries[curr_id_num].next_id = hashtable[hashval];
        hashtable[hashval] = curr_id_num++;
    }

    /* for each mask location, find the query sequence containing
       that mask and add to the list of filter locations for
       that query. Note that this assumes IDs for all query
       sequences are unique */

    for (mask_slp = mask_seqlocs; mask_slp; mask_slp = mask_slp->next) {
       SeqLocPtr current_mask = (SeqLocPtr) mask_slp->data.ptrvalue;
       Uint4 hashval;
       SeqIdPtr mask_id;
       Char buffer[64];

       if (current_mask == NULL)
           continue;

       mask_id = SeqLocId(current_mask);
       SeqIdLabel(mask_id, buffer, sizeof(buffer), OM_LABEL_CONTENT);
       hashval = readdb_sequence_hash(buffer, strlen(buffer));
       hashval = hashval >> (32 - kLog2HashSize);

       /* examine only the query IDs that hash to the same value */
       if (hashtable[hashval] != 0) {
          Int4 offset = hashtable[hashval];
          while (offset != 0) {

             SeqIdHash *q_entry = id_entries + offset;

             if (SeqIdMatch(mask_id, q_entry->id)) {
                Int4 context_idx = kNumContexts * q_entry->query_index;
                retval->seqloc_array[context_idx] = 
                              s_BlastSeqLocFromSeqLoc(current_mask,
                                         retval->seqloc_array[context_idx]);
                break;
             }
             offset = q_entry->next_id;
          }
       }
    }

    sfree(hashtable);
    sfree(id_entries);

    /* iterate through the query sequences and compute
       the complement of the filtering locations for each */

    for (query_slp = query_seqlocs, query_index = 0; 
         query_slp; 
         query_slp = query_slp->next, query_index++) {

        const int kCtxIndex = kNumContexts * query_index; /* context index */
         
        if (retval->seqloc_array[kCtxIndex])
        {
            const Boolean kIsNa = Blast_QueryIsNucleotide(program_number) &&
                !Blast_QueryIsTranslated(program_number) &&
                !Blast_ProgramIsPhiBlast(program_number);
            BlastSeqLoc_RestrictToInterval(&retval->seqloc_array[kCtxIndex], 
                                           SeqLocStart(query_slp), 
                                           SeqLocStop(query_slp));
            if (kIsNa) {
                /* N.B.: Unlike in the C++ APIs, this logic is only applied to
                 * non-translated nucleotide queries. See comment for
                 * BlastMaskLocDNAToProtein */
                Uint1 strand = SeqLocStrand(query_slp);
                if (strand == Seq_strand_minus) {
                    retval->seqloc_array[kCtxIndex+1] = 
                        retval->seqloc_array[kCtxIndex];
                    retval->seqloc_array[kCtxIndex] = NULL;
                } else if (strand == Seq_strand_plus) {
                    retval->seqloc_array[kCtxIndex+1] = NULL;
                } else {
                    retval->seqloc_array[kCtxIndex+1] = 
                        BlastSeqLocListDup(retval->seqloc_array[kCtxIndex]);
                }
            }
        }
    }
    
    return retval;
}

SeqLoc*
Blast_ValNodeMaskListFree(SeqLoc* mask_loc)
{
    ValNode* mask_var;
    for (mask_var = mask_loc; mask_var; mask_var = mask_var->next) 
        SeqLocSetFree((SeqLoc*)mask_var->data.ptrvalue);
    mask_loc = ValNodeFree(mask_loc);
    return mask_loc;
}

SeqLocPtr BlastMaskLocToSeqLoc(EBlastProgramType program_number, 
                               const BlastMaskLoc* mask_loc, 
                               SeqLoc* query_loc)
{
   SeqLocPtr retval = NULL, retval_tail = NULL;
   Int4 index;
   const Boolean k_translate = Blast_QueryIsTranslated(program_number);
   const Uint1 k_num_frames = BLAST_GetNumberOfContexts(program_number);
   const Boolean kIsNucl = (program_number == eBlastTypeBlastn);
   SeqLoc* slp;
   Boolean all_minus = TRUE;

   if (mask_loc == NULL || mask_loc->seqloc_array == NULL)
      return NULL;

   for (slp = query_loc; slp; slp = slp->next)
   {
        Uint1 strand = SeqLocStrand(slp);
        if (strand != Seq_strand_minus)
        {
           all_minus = FALSE;
           break;
        }
   }

   for (index=0, slp = query_loc; slp; ++index, slp = slp->next)
   {
      const int kCtxIndex = k_num_frames * index; /* context index */
      Int4 tmp_index;
      Int4 slp_from = SeqLocStart(slp);
      SeqIdPtr seqid = SeqLocId(slp);
      for (tmp_index=kCtxIndex; tmp_index<(kCtxIndex+k_num_frames); tmp_index++)
      {
         BlastSeqLoc* loc = NULL;
         SeqLocPtr mask_slp_head = NULL, mask_slp_tail = NULL;
         if (all_minus || BlastIsReverseStrand(kIsNucl , tmp_index) == FALSE)
         {
            for (loc = mask_loc->seqloc_array[tmp_index]; loc; loc = loc->next)
            {
               SeqIntPtr si = SeqIntNew();
               si->from = loc->ssr->left + slp_from;
               si->to = loc->ssr->right + slp_from;
               si->id = SeqIdDup(seqid);
               /* Append the pointer, but also keep track of the tail of the list
                * so that appending to the list is a constant operation */
               mask_slp_tail = ValNodeAddPointer
                   ( (mask_slp_tail ? &mask_slp_tail : &mask_slp_head), 
                     SEQLOC_INT, si);
            }
         }

         if (mask_slp_head) {
            SeqLocPtr new_mask_slp = ValNodeAddPointer(NULL, SEQLOC_PACKED_INT, 
                                             mask_slp_head);
            Uint1 tmp_choice = 0;
            /* The 'choice' of the SeqLoc in masks should show the frame,
               with values 1..6 when queries are translated; otherwise
               it does not matter. */
            if (k_translate)  
                tmp_choice = (tmp_index % NUM_FRAMES) + 1;
            else
                tmp_choice = 0;

            /* Append the pointer, but also keep track of the tail of the list
             * so that appending to the list is a constant operation */
            retval_tail = ValNodeAddPointer
                ( (retval_tail ? &retval_tail : &retval), 
                  tmp_choice, new_mask_slp);
        }
      }
   }
   return retval;
}

/** Set field values for one element of the context array of a
 * concatenated query.  All previous contexts should have already been
 * assigned correct values.
 * @param qinfo  Query info structure containing contexts. [in/out]
 * @param index  Index of the context to fill. [in]
 * @param length Length of this context. [in]
 */
static void
s_QueryInfoSetContextInfo(BlastQueryInfo*   qinfo,
                          Uint4             index,
                          Uint4             length)
{
    if (index) {
        Uint4 prev_loc = qinfo->contexts[index-1].query_offset;
        Uint4 prev_len = qinfo->contexts[index-1].query_length;
        
        Uint4 shift = prev_len ? prev_len + 1 : 0;
        
        qinfo->contexts[index].query_offset = prev_loc + shift;
        qinfo->contexts[index].query_length = length;
        if (length == 0)
            qinfo->contexts[index].is_valid = FALSE;

    } else {
        /* First context */
        qinfo->contexts[0].query_offset = 0;
        qinfo->contexts[0].query_length = length;
        if (length == 0)
            qinfo->contexts[0].is_valid = FALSE;
    }
}

/** Sets up the query information structure with all contexts' data.
 * @param slp List of query Seq-loc's [in]
 * @param program Type of BLAST program [in]
 * @param query_info_ptr Pointer to the structure to populate. [out]
 */
static Int4 
s_QueryInfoSetUp(SeqLocPtr slp, EBlastProgramType program, 
                 BlastQueryInfo** query_info_ptr)
{
   Uint4 length, protein_length;
   Boolean translate = 
      (program == eBlastTypeBlastx || program == eBlastTypeTblastx ||
       program == eBlastTypeRpsTblastn);
   Boolean is_na = (program == eBlastTypeBlastn || 
                    program == eBlastTypePhiBlastn);
   Int2 num_frames, frame;
   Uint1 strand;
   BlastQueryInfo* query_info;
   Int4 index;
   Uint4 max_length = 0;

   if (translate)
      num_frames = NUM_FRAMES;
   else if (is_na)
      num_frames = 2;
   else
      num_frames = 1;

   if ((query_info = BlastQueryInfoNew(program, ValNodeLen(slp))) == NULL) 
      return -1;

   if ((strand = SeqLocStrand(slp)) == Seq_strand_minus) {
      if (translate)
         query_info->first_context = 3;
      else
         query_info->first_context = 1;
   }
   
   /* Fill the context offsets */
   for (index = 0; slp; slp = slp->next, index += num_frames) {
      length = SeqLocLen(slp);  /* FIXME: could return -1 */
      strand = SeqLocStrand(slp);
      if (translate) {
         Int2 first_frame, last_frame;
         if (strand == Seq_strand_plus) {
            first_frame = 0;
            last_frame = 2;
         } else if (strand == Seq_strand_minus) {
            first_frame = 3;
            last_frame = 5;
         } else {
            first_frame = 0;
            last_frame = 5;
         }

         /* Set the unused initial contexts if any */
         for (frame = 0; frame < first_frame; ++frame) {
             s_QueryInfoSetContextInfo(query_info, index+frame, 0);
         }
         
         for (frame = first_frame; frame <= last_frame; ++frame) {
            protein_length = (length - frame%CODON_LENGTH)/CODON_LENGTH;
            max_length = MAX(max_length, protein_length);

            s_QueryInfoSetContextInfo(query_info,
                                     index+frame,
                                     protein_length);
         }

         /* Set the unused trailing contexts if any */
         for (frame = last_frame + 1; frame < num_frames; ++frame) {
             s_QueryInfoSetContextInfo(query_info, index+frame, 0);
         }
      } else {
         max_length = MAX(max_length, length);
         
         if (is_na) {
            if (strand == Seq_strand_plus) {
                s_QueryInfoSetContextInfo(query_info, index,   length);
                s_QueryInfoSetContextInfo(query_info, index+1, 0);
            } else if (strand == Seq_strand_minus) {
                s_QueryInfoSetContextInfo(query_info, index,   0);
                s_QueryInfoSetContextInfo(query_info, index+1, length);
            } else {
                s_QueryInfoSetContextInfo(query_info, index,   length);
                s_QueryInfoSetContextInfo(query_info, index+1, length);
            }
         } else {
             s_QueryInfoSetContextInfo(query_info, index, length);
         }
      }
   }
   query_info->max_length = max_length;

   *query_info_ptr = query_info;
   return 0;
}

/** Given a SeqLoc, fills a preallocated sequence buffer in the correct 
 * encoding.
 * @param slp SeqLoc structure to get data from. [in]
 * @param encoding What encoding to fill sequence buffer in? [in]
 * @param buffer Buffer to fill.
 */
static Int2 
s_SeqLocReadSequence(SeqLocPtr slp, EBlastEncoding encoding, Uint1** buffer)
{
   Uint1* buffer_var = *buffer;
   Int4 size, index;

   if (!buffer_var || !slp)
      return -1;

   size = SeqPortStreamLoc(slp, STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL,
                          buffer_var, NULL);

   switch (encoding) {
   case eBlastEncodingProtein: 
      for (index = 0; index < size; index++) {
         Uint1 letter = buffer_var[index];
         if (letter == 'U' || letter == 'O' || letter == 'J')
            buffer_var[index] = AMINOACID_TO_NCBISTDAA['X'];
         else
            buffer_var[index] = AMINOACID_TO_NCBISTDAA[letter];
      }
      break;
   case eBlastEncodingNcbi4na:
      for (index = 0; index < size; index++) 
         buffer_var[index] = IUPACNA_TO_NCBI4NA[buffer_var[index]];
      break;
   case eBlastEncodingNucleotide:
      for (index = 0; index < size; index++) 
         buffer_var[index] = IUPACNA_TO_BLASTNA[buffer_var[index]];
      break;
   default:
       /* This function should not be called for any other encodings - 
          return an error status. */
       return -1;
   }

   *buffer = buffer_var + size;
   return 0;
}

/** Fills sequence buffer for a single SeqLoc; fills both strands if necessary.
 * @param slp Sequence location [in]
 * @param encoding Encoding to use for the sequence buffer. [in]
 * @param add_sentinel_bytes Should sentinel bytes be added at the ends of the 
 *                           buffer? [in]
 * @param both_strands Should buffer include both strands for a nucleotide 
 *                     sequence? [in]
 * @param buffer Buffer to populate. Must be already allocated. [in] [out]
 */
static Int2 
s_SeqLocFillSequenceBuffer(SeqLocPtr slp, EBlastEncoding encoding, 
    Boolean add_sentinel_bytes, Boolean both_strands, Uint1* buffer)
{
   Uint1* buffer_var;
   Uint1 sentinel = (encoding == eBlastEncodingNucleotide 
                     ? NCBI4NA_TO_BLASTNA[NULLB] 
                     : NULLB);
   Uint1 seq_code, strand;

   buffer_var = buffer;

   if (add_sentinel_bytes) {
      *buffer_var = sentinel;
      ++buffer_var;
   }

   if (encoding == eBlastEncodingProtein) {
      seq_code = Seq_code_ncbistdaa;
      strand = Seq_strand_unknown;
   } else {
      seq_code = Seq_code_ncbi4na;
      strand = SeqLocStrand(slp);
   }

   s_SeqLocReadSequence(slp, encoding, &buffer_var);

   if (add_sentinel_bytes)
      *buffer_var = sentinel;

   if (both_strands && strand == Seq_strand_both) {
      SeqLocPtr tmp_slp=NULL;

      ++buffer_var;

      tmp_slp = SeqLocIntNew(SeqLocStart(slp), SeqLocStop(slp),
                             Seq_strand_minus, SeqLocId(slp));
            
      s_SeqLocReadSequence(tmp_slp, encoding, &buffer_var);
      if (add_sentinel_bytes)
         *buffer_var = sentinel;

      SeqLocFree(tmp_slp);
   }

   return 0;
}

Int2 BLAST_GeneticCodeFind(Int4 gc, Uint1** genetic_code)
{
   ValNodePtr vnp;
   GeneticCodePtr gcp;
   char* gen_code_eaa = NULL;
   Uint1* gen_code_stdaa = NULL;
   Int4 gen_code_length = 0, index;
   SeqMapTablePtr smtp;

   gcp = GeneticCodeFind(gc, NULL);
   for (vnp = (ValNodePtr)gcp->data.ptrvalue; vnp != NULL; 
        vnp = vnp->next) {
      if (vnp->choice == 3) {  /* ncbieaa */
         gen_code_eaa = (char*)vnp->data.ptrvalue;
         break;
      }
   }

   if (!gen_code_eaa)
      return -1;
   smtp = SeqMapTableFind(Seq_code_ncbistdaa, Seq_code_ncbieaa);
   gen_code_length = StrLen(gen_code_eaa);
   *genetic_code = gen_code_stdaa = (Uint1*) calloc(gen_code_length+1, 1);

   if (!gen_code_stdaa)
      return -2;

   for (index = 0; index < gen_code_length; ++index) {
      gen_code_stdaa[index] = 
         SeqMapTableConvert(smtp, gen_code_eaa[index]);
   }
   
   return 0;
}

/** s_GetSequence
 * Purpose:     Get the sequence for the BLAST engine, put in a Uint1 buffer
 * @param slp SeqLoc to extract sequence for [in]
 * @param query_info The query information structure, pre-initialized,
 *                   but filled here [in]
 * @param query_options Query setup options, containing the genetic code for
 *                      translation. N.B.: its strand_option field is ignored [in]
 * @param num_frames How many frames to get for this sequence? [in]
 * @param encoding In what encoding to retrieve the sequence? [in]
 * @param buffer_out Buffer to hold plus strand or protein [out]
 * @param buffer_length Length of buffer allocated [out]
 */
static Int2 
s_GetSequence(SeqLocPtr slp, BlastQueryInfo* query_info, 
   const QuerySetUpOptions* query_options, Uint1 num_frames, 
   EBlastEncoding encoding, Uint1* *buffer_out, Int4 *buffer_length)
{
   Int2		status=0; /* return value. */
   Int4 total_length; /* Total length of all queries/frames/strands */
   Int4		index; /* Loop counter */
   SeqLocPtr	slp_var; /* loop variable */
   Uint1*	buffer; /* buffer to fill. */
   Boolean add_sentinel_bytes = TRUE;
   Uint1* genetic_code=NULL;
   Boolean translate = FALSE;
   Int4 offset = 0;

   if (query_info) {
       *buffer_length = total_length = QueryInfo_GetSeqBufLen(query_info);
   } else {
      /* Subject sequence in 2 sequences comparison */
      *buffer_length = SeqLocLen(slp);
      /* allow two extra bytes for sentinels or a trailing
         null appended by the low-level sequence conversion */
      total_length = (*buffer_length) + 2;
      if (encoding == eBlastEncodingNcbi4na) {
         /* Searches with translated subjects (tblastn, tblastx) */
         add_sentinel_bytes = FALSE;
      }
   }

   if (num_frames == NUM_FRAMES) {
      /* Sequence must be translated in 6 frames. This can only happen
         for query - subject sequences are translated later. */
      Int4 gc;
      
      translate = TRUE;
      gc = (query_options ? query_options->genetic_code : 1);

      if ((status = BLAST_GeneticCodeFind(gc, &genetic_code)) != 0)
         return status;
   }

   *buffer_out = buffer = (Uint1 *) malloc((total_length)*sizeof(Uint1));
   
   for (index = 0, slp_var = slp; slp_var; 
        slp_var = slp_var->next, index += num_frames)
   {
      if (translate) {
         Uint1* na_buffer, *buffer_rev = NULL;
         Int4 context, context_start, context_end;
         Int4 na_length;
         Uint1 strand;
         

         na_length = SeqLocLen(slp_var);
         strand = SeqLocStrand(slp_var);
         /* Retrieve nucleotide sequence in an auxiliary buffer; 
            then translate into the appropriate place in the 
            preallocated buffer */
         if (strand == Seq_strand_plus) {
            na_buffer = (Uint1 *) malloc(na_length + 2);
            context_start = 0;
            context_end = 2;
         } else if (strand == Seq_strand_minus) {
            na_buffer = (Uint1 *) malloc(na_length + 2);
            context_start = 3;
            context_end = 5;
         } else {
            na_buffer = (Uint1*) malloc(2*na_length + 3);
            context_start = 0;
            context_end = 5;
         }
         s_SeqLocFillSequenceBuffer(slp_var, encoding, TRUE, TRUE, na_buffer);
         if (strand == Seq_strand_both)
            buffer_rev = na_buffer + na_length + 1;
	 else if (strand == Seq_strand_minus)
	    buffer_rev = na_buffer;

         for (context = context_start; context <= context_end; context++) {
             offset = query_info->contexts[index+context].query_offset;
             
             BLAST_GetTranslation(na_buffer+1, buffer_rev, na_length,
                                  BLAST_ContextToFrame(eBlastTypeBlastx, context),
                                  &buffer[offset], genetic_code);
         }
         sfree(na_buffer);
      } else {
         /* This can happen both for query and subject, so query_info 
            might not be initialized here. */
          if (query_info)
              offset = query_info->contexts[index].query_offset;
          
          s_SeqLocFillSequenceBuffer(slp_var, encoding, add_sentinel_bytes, 
                                     (Boolean)(num_frames == 2), &buffer[offset]);
      }
      /* For subjects, do only one SeqLoc at a time */
      if (!query_info)
         break;
   }

   sfree(genetic_code);

   return status;
}

Int2 BLAST_SetUpQuery(EBlastProgramType program_number, 
        SeqLocPtr query_slp, const QuerySetUpOptions* query_options, 
        SeqLoc* masking_locs, BlastQueryInfo** query_info, 
        BLAST_SequenceBlk* *query_blk)
{
   Uint1* buffer;	/* holds sequence for plus strand or protein. */
   Int4 buffer_length;
   Int2 status;
   Uint1 num_frames;
   EBlastEncoding encoding;

   if (query_slp == NULL || query_options == NULL ||
       query_info == NULL || query_blk == NULL)
      return -1;

   if ((status = s_QueryInfoSetUp(query_slp, program_number, query_info)))
      return status;

   if (program_number == eBlastTypeBlastn || 
       program_number == eBlastTypePhiBlastn) {
      encoding = eBlastEncodingNucleotide;
      num_frames = 2;
   } else if (program_number == eBlastTypeBlastp ||
              program_number == eBlastTypeRpsBlast ||
              program_number == eBlastTypeTblastn ||
              program_number == eBlastTypePhiBlastp) {
      encoding = eBlastEncodingProtein;
      num_frames = 1;
   } else { /* blastx or rpstblastn, which is also essentially blastx */
      encoding = eBlastEncodingNcbi4na;
      num_frames = NUM_FRAMES;
   }

   if ((status=s_GetSequence(query_slp, *query_info, query_options,
                  num_frames, encoding, &buffer, &buffer_length)))
      return status; 
        
   /* Do not count the first and last sentinel bytes in the 
      query length */
   if ((status=BlastSetUp_SeqBlkNew(buffer, buffer_length-2, 
                                    0, query_blk, TRUE)))
      return status;

   if (masking_locs) {
       BlastMaskLoc* lcase_mask = BlastMaskLocFromSeqLoc(masking_locs, 
                                                         query_slp,
                                                         program_number);
       if (Blast_QueryIsTranslated(program_number))
           BlastMaskLocDNAToProtein(lcase_mask, *query_info);
       (*query_blk)->lcase_mask = lcase_mask;
       (*query_blk)->lcase_mask_allocated = TRUE;
   }

   return 0;
}

Int2 BLAST_SetUpSubject(EBlastProgramType program_number, 
        SeqLocPtr subject_slp, BLAST_SequenceBlk** subject)
{
   Int2 status = 0;
   Uint1* subject_buffer = NULL; /* Buffer for the compressed subject 
                                      sequence in two sequences case */
   Int4 buffer_length=0; /* Length of subject sequence for two sequences 
                            case */
   EBlastEncoding encoding;
   const Boolean kNucleotide = (program_number == eBlastTypeBlastn || 
                                program_number == eBlastTypePhiBlastn);
   const Boolean kTranslated = (program_number == eBlastTypeTblastn ||
                                program_number == eBlastTypeTblastx);

   if (kNucleotide)
      encoding = eBlastEncodingNucleotide;
   else if (kTranslated) {
      encoding = eBlastEncodingNcbi4na;
   } else {
      encoding = eBlastEncodingProtein;
   }

   if ((status = s_GetSequence(subject_slp, NULL, NULL, 1, encoding,
                                   &subject_buffer, &buffer_length)))
      return status;
   
   /* Initialize the sequence block, saving the sequence buffer in 
      'sequence_start'. */
   if ((status=BlastSetUp_SeqBlkNew(subject_buffer, buffer_length,
                                    0, subject, TRUE)))
      return status;

   /* If subject sequence is nucleotide, create compressed sequence buffer
      and save it in 'sequence'. For blastn, the sentinel bytes should not 
      be included in the packed sequence. */
   if (kNucleotide)
      ++subject_buffer;

   if (kNucleotide || kTranslated) {
      BLAST_PackDNA(subject_buffer, buffer_length, encoding, 
                    &((*subject)->sequence));
      (*subject)->sequence_allocated = TRUE;
   }

   return 0;
}
/* @} */

