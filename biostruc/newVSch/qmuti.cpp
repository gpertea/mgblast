/*
 * $Id: qmuti.cpp,v 1.1 2005/07/26 17:11:46 chenj Exp $
 *
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *            National Center for Biotechnology Information (NCBI)
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government do not place any restriction on its use or reproduction.
 *  We would, however, appreciate having the NCBI and the author cited in
 *  any work or product based on this material
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 * ===========================================================================
 *
 * Author: Jie Chen
 * Source file for VSMmdb.cgi
 *
 * $Log: qmuti.cpp,v $
 * Revision 1.1  2005/07/26 17:11:46  chenj
 * Making linux VSMmdb.cgi
 *
 *
 *
 */



#include "qmuti.hpp"


using namespace ncbi;
using namespace std;

void ChkReqId(CCgiContext& ctx, QmReqID reqId)
{

     CCgiResponse& response = ctx.GetResponse();
     static string strtmp;

     if (!reqId) {

          PrtMes::PrintErrorHeader(&response, "VastSrv", 1);
          PrtMes::PrintMsgWithTail(&response, "ReqId = 0");

          ERR_POST("ReqId=0");
          throw exception();

     }

     QmReqStatus status = sqmGetStatus(reqId);

     if (status != qmReqStatusUnknown) return;
     else {

        strtmp = "Invalid ReqId = " + ReqId2Str(reqId) + ".<br>"
                + "The search results are kept for only 10 days, please resubmit your query. <br> If the ReqId is wrong, please try again.";

        PrtMes::PrintErrorHeader(&response, "VastSrv", 1);
        PrtMes::PrintMsgWithTail(&response, strtmp);

        ERR_POST(strtmp);
        throw exception();

     }


}       // ChkReqId();



void ChkGrpId(CCgiContext& ctx, QmReqID ori_gid, QmReqID new_gid)
{
	CCgiResponse& resp = ctx.GetResponse();
	string 	strtmp;
	
	if (!new_gid || (ori_gid && ori_gid != new_gid)) {

 	    if (!new_gid) strtmp = "Group Id = 0";
	    else 
	       strtmp = (string)"Got wrong group id when associating with a new reqid.\n"
			+ "<br>Original gid = " + ReqId2Str(ori_gid)
			+ "<br>New gid = " + ReqId2Str(new_gid);
	    
	    PrtMes::PrintErrorHeader(&resp, "VastSrch", 1);
	    PrtMes::PrintMsgWithTail(&resp, "Group Id = 0");
	
            ERR_POST(strtmp);
            throw exception();

	}
	else return;

} 	// ChkGrpId()
