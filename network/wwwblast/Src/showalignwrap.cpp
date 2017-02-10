/*  $Id: showalignwrap.cpp,v 1.2 2006/03/07 17:06:30 jianye Exp $
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
 *   c adaptor to showalign.cpp
 *
 */


#include <corelib/ncbiapp.hpp>
#include <corelib/ncbiargs.hpp>
#include <corelib/ncbienv.hpp>
#include <corelib/ncbistre.hpp>
#include <serial/iterator.hpp>
#include <serial/objistr.hpp>
#include <serial/objostr.hpp>
#include <serial/serial.hpp>
#include <objtools/blast_format/showalign.hpp>
#include <showalignwrap.h>
#include <ctools/asn_converter.hpp>
#include <objtools/data_loaders/blastdb/bdbloader.hpp>
#include <objmgr/object_manager.hpp>
#include <objmgr/gbloader.hpp>
#include <algo/blast/api/blast_aux.hpp>



USING_SCOPE(ncbi);
USING_SCOPE(objects);


static void  
s_BlastFormat_ConvertToCppLoc(list <CRef<blast::CSeqLocInfo> >& target_mask,
                              ValNodePtr source_mask)
{
    DECLARE_ASN_CONVERTER(CSeq_loc, SeqLoc, convertor_feat); 
    ValNodePtr  mask_temp = source_mask;
    while(mask_temp){
        int frame(0);
        switch(mask_temp->choice) {
        case Blast_mask_frame_minus1: frame = -1; break;
        case Blast_mask_frame_minus2: frame = -2; break;
        case Blast_mask_frame_minus3: frame = -3; break;
        case Blast_mask_frame_plus1:  frame =  1; break;
        case Blast_mask_frame_plus2:  frame =  2; break;
        case Blast_mask_frame_plus3:  frame =  3; break;
        default:                      frame =  0; break;
        }
      

        SeqLocPtr slp = (SeqLocPtr) mask_temp->data.ptrvalue;
        
        if (slp->choice == SEQLOC_PACKED_INT){
            for (ValNodePtr vnpTemp2 = (ValNodePtr) slp->data.ptrvalue;  vnpTemp2; vnpTemp2 = vnpTemp2->next){
                CRef<CSeq_loc> seqloc(new CSeq_loc());
                convertor_feat.FromC((SeqLocPtr)vnpTemp2, &*seqloc); 
                CSeq_interval* interval = &seqloc->SetInt();
                CRef<blast::CSeqLocInfo> seqlocinfo
                    (new blast::CSeqLocInfo(interval, frame));
                
                target_mask.push_back(seqlocinfo);
            }
        } else if (slp->choice == SEQLOC_INT) {
            CRef<CSeq_loc> seqloc(new CSeq_loc());
            convertor_feat.FromC(slp, &*seqloc); 
            CRef<blast::CSeqLocInfo> seqlocinfo
                (new blast::CSeqLocInfo(&seqloc->SetInt(), frame)); 
            target_mask.push_back(seqlocinfo);
        }
        mask_temp=mask_temp->next;
        
    }
} 

void DisplayAlign(SeqAlignPtr align, Int4 line_len, BioseqPtr query,
                  BioseqPtr subject, CharPtr program, Int4Ptr PNTR matrix,
                  ValNodePtr mask, Int4 mask_char, Int4 mask_color,
                  Boolean cds_translation, Int4 view, FILE *fp)
{
   
    CSeq_align_set seqalignSet;
    DECLARE_ASN_CONVERTER(CSeq_align, SeqAlign, convertor_seqalign);
    
    for (SeqAlignPtr salpTemp = align; salpTemp; salpTemp = salpTemp->next){
        CRef<CSeq_align> csal(new CSeq_align());
        convertor_seqalign.FromC(salpTemp, &*csal);
        seqalignSet.Set().push_back(csal);
    }
    
    CRef<CObjectManager> obj;
    obj = CObjectManager::GetInstance();
    CGBDataLoader::RegisterInObjectManager(*obj);       
    CRef<CScope> scope (new CScope(*obj));
    scope->AddDefaults();
    
  
    DECLARE_ASN_CONVERTER(CBioseq, Bioseq, convertor2);  
  
    CRef<CBioseq> cbsp(new CBioseq());
    convertor2.FromC(subject, &*cbsp); 
    CRef<CSeq_entry> entry(new CSeq_entry());
    entry->SetSeq(*cbsp);
    scope->AddTopLevelSeqEntry(*entry);
    
    DECLARE_ASN_CONVERTER(CBioseq, Bioseq, convertor3); 
    CRef<CBioseq> cbsp2(new CBioseq());
    convertor3.FromC(query, &*cbsp2); 
    CRef<CSeq_entry> entry2(new CSeq_entry());
    entry2->SetSeq(*cbsp2);
    scope->AddTopLevelSeqEntry(*entry2);
    //auto_ptr<CObjectOStream> out2(CObjectOStream::Open(eSerial_AsnText, cout));
    //*out2 << seqalignSet;
    int option = 0;

    list <CRef<blast::CSeqLocInfo> > mask_loc;

    s_BlastFormat_ConvertToCppLoc(mask_loc, mask);
  
    int the_matrix[CDisplaySeqalign::ePMatrixSize][CDisplaySeqalign::ePMatrixSize];
    if(matrix){
        for(int i = 0; i < CDisplaySeqalign::ePMatrixSize; i ++){
            for(int j = 0; j < CDisplaySeqalign::ePMatrixSize; j ++){
                the_matrix[i][j] = matrix[i][j];
            }
        }
    }
    CDisplaySeqalign cds(seqalignSet,  *scope, &mask_loc, NULL, 
                         matrix ? the_matrix : NULL);
    CNcbiOstrstream oss;
    if(strcmp(program, "tblastx") == 0){
        option += CDisplaySeqalign::eTranslateNucToNucAlignment;
    }
    option += CDisplaySeqalign::eShowBlastStyleId;
    if (cds_translation) {
        option += CDisplaySeqalign::eShowCdsFeature;
    }
    option += CDisplaySeqalign::eShowBlastInfo;
    option += CDisplaySeqalign::eShowNoDeflineInfo;
    option += CDisplaySeqalign::eHtml;
    if (view == 2){
        option += CDisplaySeqalign::eColorDifferentBases;
        option += CDisplaySeqalign::eShowIdentity;
    } else {
        option += CDisplaySeqalign::eShowMiddleLine;
    }
    
    if (mask_char == 2){
        cds.SetSeqLocChar(CDisplaySeqalign::eLowerCase);
    } else if (mask_char == 0){
        if(strcmp(program, "blastn") == 0){
            cds.SetSeqLocChar (CDisplaySeqalign::eN);
        }else {
            cds.SetSeqLocChar (CDisplaySeqalign::eX);
        }
    }
    
    if (mask_color == 2){
        cds.SetSeqLocColor(CDisplaySeqalign::eGrey);
    } else if (mask_color == 3){
        cds.SetSeqLocColor(CDisplaySeqalign::eRed);
    }
    
    if(strcmp(program, "blastn") == 0){//blastn
        cds.SetMiddleLineStyle (CDisplaySeqalign::eBar);
        cds.SetAlignType(CDisplaySeqalign::eNuc);
    } else {
        cds.SetMiddleLineStyle (CDisplaySeqalign::eChar);
        cds.SetAlignType(CDisplaySeqalign::eProt);
    }
   
    cds.SetAlignOption(option);
    cds.DisplaySeqalign(oss);
    string output = (string)(CNcbiOstrstreamToString(oss));
    fprintf(fp, "%s", output.c_str());
   
}


/* 
*============================================================
*$Log: showalignwrap.cpp,v $
*Revision 1.2  2006/03/07 17:06:30  jianye
*correct mask char
*
*Revision 1.1  2006/01/05 17:28:57  jianye
*for wblast2 to use new formatter
*

*/
