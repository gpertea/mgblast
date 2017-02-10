/*  $Id: showalignwrap.h,v 1.1 2006/01/05 17:28:57 jianye Exp $
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
 * Author:  Jian Ye
 *
 * File Description:
 *   c adaptor to showalign.hpp
 *
 */

#ifndef SHOWALIGN_WRAP_H
#define SHOWALIGN_WRAP_H

#include <ncbi.h>
#include <objseq.h>
#include <objsset.h>
#include <sequtil.h>
#include <seqport.h>
#include <tofasta.h>
#include <blast.h>
#include <blastpri.h>

#include <gapxdrop.h>
#include <salsap.h>
#include <objalign.h>
#include <accentr.h>
#include <objblst3.h>



#ifdef __cplusplus 
extern "C" {
#endif
void DisplayAlign(SeqAlignPtr align, Int4 line_len, BioseqPtr query, 
                  BioseqPtr subject, CharPtr program, 
                  Int4Ptr PNTR matrix,
                  ValNodePtr mask,  Int4 mask_char, Int4 mask_color,
                  Boolean cds_translation, Int4 view, FILE *fp);


#ifdef __cplusplus
}
#endif

#endif


/* 
*============================================================
*$Log: showalignwrap.h,v $
*Revision 1.1  2006/01/05 17:28:57  jianye
*for wblast2 to use new formatter
*

*/
