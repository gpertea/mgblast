/* $Id: gapinfo.c,v 1.18 2006/02/15 15:06:00 madden Exp $
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
 * Author: Ilya Dondoshansky
 *
 */

/** @file gapinfo.c
 * Initialization and freeing of structures containing traceback information.
 */


#ifndef SKIP_DOXYGEN_PROCESSING
static char const rcsid[] = 
    "$Id: gapinfo.c,v 1.18 2006/02/15 15:06:00 madden Exp $";
#endif /* SKIP_DOXYGEN_PROCESSING */

#include <algo/blast/core/gapinfo.h>

GapStateArrayStruct* 
GapStateFree(GapStateArrayStruct* state_struct)

{
    GapStateArrayStruct* next;

    while (state_struct) {
        next = state_struct->next;
        sfree(state_struct->state_array);
        sfree(state_struct);
        state_struct = next;
    }

    return NULL;
}

/* see gapinfo.h for description */
GapEditScript* 
GapEditScriptNew(Int4 size)

{
    GapEditScript* new;

    if (size <= 0) 
       return NULL;

    new = (GapEditScript*) calloc(1, sizeof(GapEditScript));
    if (new)
    {
         new->size = size;
         new->op_type = (EGapAlignOpType*) calloc(size, sizeof(EGapAlignOpType));
         new->num = (Int4*) calloc(size, sizeof(Int4));
    }
    return new;
}

GapEditScript*
GapEditScriptDelete(GapEditScript* old)
{

    if (old)
    {
       sfree(old->op_type);
       sfree(old->num);
       sfree(old);
    }
    return old;
}

/* see gapinfo.h for description */
GapEditScript* 
GapEditScriptDup(const GapEditScript* old)
{
    GapEditScript* new = NULL;
    if (old)
    {
        new = GapEditScriptNew(old->size);
        if (new)
        {
            memcpy(new->op_type, old->op_type, (old->size)*sizeof(EGapAlignOpType));
            memcpy(new->num, old->num, (old->size)*sizeof(Int4));
        }
    }
    return new;
}

Int2
GapEditScriptPartialCopy(GapEditScript* new, int offset, const GapEditScript* old, int start, int stop)
{
     int size = 1+stop-start;
     int new_index = 0;
     int old_index = 0;

     if (old == NULL || new == NULL || stop-start+1 > new->size)
        return -1;

     old_index = start;
     for (new_index=offset; new_index<size+offset; new_index++)
     {
         new->num[new_index] = old->num[old_index];
         new->op_type[new_index] = old->op_type[old_index];
         old_index++;
     }

     return 0;
}

/** Ensures that a preliminary edit script has enough memory allocated
 *  to hold a given number of edit operations
 *
 *  @param edit_block The script to examine [in/modified]
 *  @param total_ops Number of operations the script must support [in]
 *  @return 0 if successful, nonzero otherwise
 */
static Int2 
s_GapPrelimEditBlockRealloc(GapPrelimEditBlock *edit_block, Int4 total_ops)
{
    if (edit_block->num_ops_allocated <= total_ops) {
        Int4 new_size = total_ops * 2;
        GapPrelimEditScript *new_ops;
    
        new_ops = realloc(edit_block->edit_ops, new_size * 
                                sizeof(GapPrelimEditScript));
        if (new_ops == NULL)
            return -1;

        edit_block->edit_ops = new_ops;
        edit_block->num_ops_allocated = new_size;
    }
    return 0;
}

/** Add an edit operation to an edit script

  @param edit_block The edit script to update [in/modified]
  @param op_type The edit operation to add [in]
  @param num_ops The number of operations of the specified type [in]
  @return 0 on success, nonzero otherwise
*/
static Int2 
s_GapPrelimEditBlockAddNew(GapPrelimEditBlock *edit_block, 
                           EGapAlignOpType op_type, Int4 num_ops)
{
    if (s_GapPrelimEditBlockRealloc(edit_block, edit_block->num_ops + 2) != 0)
        return -1;

    ASSERT(op_type != eGapAlignInvalid);

    edit_block->last_op = op_type;
    edit_block->edit_ops[edit_block->num_ops].op_type = op_type;
    edit_block->edit_ops[edit_block->num_ops].num = num_ops;
    edit_block->num_ops++;

    return 0;
}

void
GapPrelimEditBlockAdd(GapPrelimEditBlock *edit_block, 
                 EGapAlignOpType op_type, Int4 num_ops)
{
    if (num_ops == 0)
        return;

    if (edit_block->last_op == op_type)
        edit_block->edit_ops[edit_block->num_ops-1].num += num_ops;
    else
        s_GapPrelimEditBlockAddNew(edit_block, op_type, num_ops);
}

GapPrelimEditBlock *
GapPrelimEditBlockNew(void)
{
    GapPrelimEditBlock *edit_block = malloc(sizeof(GapPrelimEditBlock));
    if (edit_block != NULL) {
        edit_block->edit_ops = NULL;
        edit_block->num_ops_allocated = 0;
        edit_block->num_ops = 0;
        edit_block->last_op = eGapAlignInvalid;
        s_GapPrelimEditBlockRealloc(edit_block, 100);
    }
    return edit_block;
}

GapPrelimEditBlock *
GapPrelimEditBlockFree(GapPrelimEditBlock *edit_block)
{
    if (edit_block == NULL)
        return NULL;

    sfree(edit_block->edit_ops);
    sfree(edit_block);
    return NULL;
}

void
GapPrelimEditBlockReset(GapPrelimEditBlock *edit_block)
{
    if (edit_block) {
        edit_block->num_ops = 0;
        edit_block->last_op = eGapAlignInvalid;
    }
}

void
GapPrelimEditBlockAppend(GapPrelimEditBlock *edit_block1,
                         GapPrelimEditBlock *edit_block2)
{
    Int4 index;
    GapPrelimEditScript *op;

    for (index = 0, op = edit_block2->edit_ops; index < edit_block2->num_ops; 
         ++index, ++op) {
        GapPrelimEditBlockAdd(edit_block1, op->op_type, op->num);
    }
}
