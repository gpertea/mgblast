/* $Id: blast_options_api.c,v 1.18 2006/04/26 12:45:28 madden Exp $
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
*  Author: Ilya Dondoshansky                                              *
**************************************************************************/

/** @file blast_options_api.c
 * Functions for C toolkit applications to perform a BLAST search, using the 
 * core engine, shared between C and C++ toolkits.
 */

#include <algo/blast/api/blast_options_api.h>
#include <algo/blast/core/blast_util.h>
#include <algo/blast/core/blast_filter.h>
#include <algo/blast/api/blast_seq.h>

/** @addtogroup CToolkitAlgoBlast
 *
 * @{
 */

Int2 SBlastOptionsNew(const char* program_name, SBlastOptions** options_out,
                      Blast_SummaryReturn* extra_returns)
{
   QuerySetUpOptions* query_options=NULL;	
   LookupTableOptions* lookup_options=NULL;
   BlastInitialWordOptions* word_options=NULL;
   BlastScoringOptions* score_options=NULL;
   BlastExtensionOptions* ext_options=NULL;
   BlastHitSavingOptions* hit_options=NULL;
   BlastEffectiveLengthsOptions* eff_len_options=NULL;
   PSIBlastOptions* psi_options = NULL;
   BlastDatabaseOptions* db_options = NULL;
   SBlastOptions* options;
   EBlastProgramType program = eBlastTypeUndefined;
   Int2 status = 0;

   if (!options_out || !extra_returns)
       return -1;

   BlastProgram2Number(program_name, &program);
   if (program == eBlastTypeUndefined) {
       char message[256];

       sprintf(message, 
               "Program name %s is not supported. The supported programs "
               "are blastn, blastp, blastx, tblastn, tblastx, rpsblast, "
               "rpstblastn\n", program_name);
       SBlastMessageWrite(&extra_returns->error, SEV_ERROR, message, NULL, FALSE);
       return -1;
   }

   status = 
       BLAST_InitDefaultOptions(program, &lookup_options, &query_options, 
           &word_options, &ext_options, &hit_options, &score_options, 
           &eff_len_options, &psi_options, &db_options);
   
   if (status) {
       *options_out = NULL;
       SBlastMessageWrite(&extra_returns->error, SEV_ERROR, "Failed to initialize default options\n", NULL, FALSE);
       return status;
   }

   if (program == eBlastTypeTblastn || program == eBlastTypeRpsTblastn ||
       program == eBlastTypeTblastx) {
       if ((status = BLAST_GeneticCodeFind(db_options->genetic_code, 
                                           &db_options->gen_code_string)))
           return status;
   }
   
   *options_out = options = (SBlastOptions*) calloc(1, sizeof(SBlastOptions));
   options->program = program;
   options->query_options = query_options;
   options->lookup_options = lookup_options;
   options->word_options = word_options;
   options->ext_options = ext_options;
   options->score_options = score_options;
   options->hit_options = hit_options;
   options->eff_len_options = eff_len_options;
   options->psi_options = psi_options;
   options->db_options = db_options;
   options->num_cpus = 1;
   options->believe_query = FALSE;

   /* Set default filter string to low complexity filtering. */
   SBlastOptionsSetFilterString(options, "T");

   return status;
}

SBlastOptions* SBlastOptionsFree(SBlastOptions* options)
{
    if (options) {
        LookupTableOptionsFree(options->lookup_options);
        BlastQuerySetUpOptionsFree(options->query_options);
        BlastExtensionOptionsFree(options->ext_options);
        BlastHitSavingOptionsFree(options->hit_options);
        BlastInitialWordOptionsFree(options->word_options);
        BlastScoringOptionsFree(options->score_options);
        BlastEffectiveLengthsOptionsFree(options->eff_len_options);
        PSIBlastOptionsFree(options->psi_options);
        BlastDatabaseOptionsFree(options->db_options);
        sfree(options);
    }
    return NULL;
}

Int2 SBlastOptionsSetEvalue(SBlastOptions* options, double evalue)
{
    if (options && options->hit_options) {
        options->hit_options->expect_value = evalue;
        return 0;
    } else
        return -1;
}

Int2 SBlastOptionsSetWordSize(SBlastOptions* options, Int4 word_size)
{
    if (options && options->lookup_options) {
        options->lookup_options->word_size = word_size;
        return 0;
    } else
        return -1;
}

Int2 SBlastOptionsSetThreshold(SBlastOptions* options, Int4 threshold)
{

    if (!options || !options->lookup_options || !options->score_options)
        return -1;

    if (threshold < 0)
       return -2;

    if (Blast_QueryIsNucleotide(options->program) == TRUE && Blast_QueryIsTranslated(options->program) == FALSE)
        return 0;

   if (threshold == 0)
   {
     Int2 status=0;
     if ((status=BLAST_GetSuggestedThreshold(options->program, options->score_options->matrix, &threshold)) != 0)
         return status;
   }

   options->lookup_options->threshold = threshold;

   return 0;
}

Int2 SBlastOptionsSetWindowSize(SBlastOptions* options, Int4 window_size)
{

   if (!options || !options->score_options || !options->word_options)
       return -1;

   if (window_size < 0)
       return -2;

   if (Blast_QueryIsNucleotide(options->program) == TRUE && Blast_QueryIsTranslated(options->program) == FALSE)
        return 0;

   if (window_size == 0)
   {
     Int2 status=0;
     if ((status=BLAST_GetSuggestedWindowSize(options->program, options->score_options->matrix, &window_size)) != 0)
         return status;
   }

   options->word_options->window_size = window_size;
}

Int2 SBlastOptionsSetDiscMbParams(SBlastOptions* options, Int4 template_length,
                                 Int4 template_type)
{
    if (!options || !options->lookup_options)
        return -1;

    options->lookup_options->lut_type = MB_LOOKUP_TABLE;
    options->lookup_options->mb_template_length = template_length;
    options->lookup_options->mb_template_type = template_type;
    
    return 0;
}

Int2 SBlastOptionsSetMatrixAndGapCosts(SBlastOptions* options, 
                                       const char* matrix_name, 
                                       Int4 gap_open, Int4 gap_extend)
{
    Int2 status = 0;

    if (!matrix_name || !options || !options->score_options)
        return -1;
    
    /* Reward penalty do not apply to blastn. */
    if (options->program == eBlastTypeBlastn)
        return 0;

    status = BLAST_FillScoringOptions(options->score_options, options->program,
                        FALSE, -1, -1, matrix_name, gap_open, gap_extend);
    if (status != 0)
        return status;
    
    if (gap_open < 0 || gap_extend < 0)
    {
        Int4 gap_open_priv = 0;
        Int4 gap_extend_priv = 0;

        BLAST_GetProteinGapExistenceExtendParams(matrix_name, &gap_open_priv, &gap_extend_priv);
        if (gap_open < 0)
            gap_open = gap_open_priv;
        if (gap_extend < 0)
            gap_extend = gap_extend_priv;
    }

    options->score_options->gap_open = gap_open;
    options->score_options->gap_extend = gap_extend;

    return status;
}

Int2 SBlastOptionsSetRewardPenaltyAndGapCosts(SBlastOptions* options, 
                                       Int4 reward, Int4 penalty,
                                       Int4 gap_open, Int4 gap_extend,
                                       Boolean greedy)
{
    Int2 status = 0;

    if (reward <= 0 || penalty >= 0 || !options || !options->score_options)
        return -1;

    /* Reward penalty only apply to blastn. */
    if (options->program != eBlastTypeBlastn)
        return 0;
    
    status = BLAST_FillScoringOptions(options->score_options, options->program,
                        greedy, penalty, reward, NULL, gap_open, gap_extend);
    if (status != 0)
        return status;
    
    if (gap_open < 0 || gap_extend < 0)
    {
        Int4 gap_open_priv;
        Int4 gap_extend_priv;
      
        if (greedy)
        {
             gap_open_priv = BLAST_GAP_OPEN_MEGABLAST;
             gap_extend_priv = BLAST_GAP_EXTN_MEGABLAST;
        }
        else
        {
             gap_open_priv = BLAST_GAP_OPEN_NUCL;
             gap_extend_priv = BLAST_GAP_EXTN_NUCL;
        }

        status = BLAST_GetNucleotideGapExistenceExtendParams(reward, penalty, &gap_open_priv, &gap_extend_priv);
        if (status)
           return status;

        if (gap_open < 0)
            gap_open = gap_open_priv;
        if (gap_extend < 0)
            gap_extend = gap_extend_priv;
    }

    options->score_options->gap_open = gap_open;
    options->score_options->gap_extend = gap_extend;

    return status;
}

Int2 SBlastOptionsSetFilterString(SBlastOptions* options, const char* str)
{
    Int2 status = 0;
    if (!options || !options->query_options)
        return -1;

    /* Reset filtering options */
    sfree(options->query_options->filter_string);
    options->query_options->filter_string = strdup(str);
    options->query_options->filtering_options =
         SBlastFilterOptionsFree(options->query_options->filtering_options);
    status = BlastFilteringOptionsFromString(options->program, 
        options->query_options->filter_string, 
        &options->query_options->filtering_options, NULL);
    return status;
}

Int2 SBlastOptionsSetDbGeneticCode(SBlastOptions* options, Int4 gc)
{
    Int2 status = 0;

    if (gc == 0)
        return 0;

    if (!options || !options->db_options)
        return -1;

    /* If previously set genetic code is the same as the new one, there is no
       need to do anything. */
    if (options->db_options->genetic_code != gc) {
        options->db_options->genetic_code = gc;
        /* Free old genetic code string. */
        sfree(options->db_options->gen_code_string);
        status = BLAST_GeneticCodeFind(gc, &options->db_options->gen_code_string);
    }

    return status;
    
}

Boolean SBlastOptionsGetMaskAtHash(const SBlastOptions* options)
{
    ASSERT(options && options->query_options &&
           options->query_options->filtering_options);

    return options->query_options->filtering_options->mask_at_hash;
}

Int2 SBlastOptionsSetBelieveQuery(SBlastOptions* options, Boolean believe_query)
{
    Int2 status = 0;

    if (!options)
       return -1;

    options->believe_query = believe_query;

    return status;
}

Boolean SBlastOptionsGetBelieveQuery(const SBlastOptions* options)
{

    ASSERT(options);

    return options->believe_query;
}

/* @} */

