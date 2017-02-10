/*
 * $Id: cVAcmd.cpp,v 1.1 2005/07/26 16:25:22 chenj Exp $
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
 *
 * $Log: cVAcmd.cpp,v $
 * Revision 1.1  2005/07/26 16:25:22  chenj
 * source for vastalign.cgi
 *
 *
 */


#include <corelib/ncbiexec.hpp>
#include <corelib/ncbifile.hpp>
#include <serial/objostr.hpp>

#include "hVAapp.hpp"
#include "hVAcmd.hpp"
#include "hUtilib.hpp"
#include "PubVastApi.hpp"

#include <string.h>
#include <sys/time.h>
#include <math.h>


using namespace ncbi;


CVACommand::CVACommand(CNcbiResource& resource) : CNcbiCommand(resource)
{};

string CVACommand::GetEntry() const
{
	return string("cmdVA"); 
}

CNcbiCommand* CVACommand :: Clone(void) const
{
	return new CVACommand(GetVAResource());
}

string CVACommand :: GetName() const
{
	return string("VastAlign");
}


void CVACommand :: Execute(CCgiContext& ctx)
{

	CCgiRequest&  req = ctx.GetRequest();
   	CCgiResponse& resp = ctx.GetResponse();
	
	string mst = req.GetEntry("master").GetValue();
	if ( mst.empty() ) {
	
		PrtMes::PrintErrorHeader(&resp, "VastAlign", 1);
      		PrtMes::PrintMsgWithTail(&resp,
            		"<h3>No master privided. </h3>n");
 	}

	if ( 4 == mst.size() ) mst += ' ';
	if ( 6 <= mst.size() ) {
		mst[4] = mst[5];
		if ( 6 == mst.size() ) mst[5] = '\0';
	}
	if ( 'x' == mst[4] ) mst[4] = ' ';
	unsigned i;
	for (i=0; i< mst.size(); i++) mst[i]  = toupper(mst[i]);

	string slv = req.GetEntry("slave").GetValue();
	if ( slv.empty() ) {

	     PrtMes::PrintErrorHeader(&resp, "VastAlign", 1);
	     PrtMes::PrintMsgWithTail(&resp, "<h3>No slave privided. </h3>n");
        }

	if ( 4 == slv.size() ) slv += ' ';
	if ( 6 <= slv.size() ) {
		slv[4] = slv[5];
		if ( 6 == slv.size()) slv[5] = '\0';
	}
	
        if ( 'x' == slv[4] ) slv[4] = ' ';
	for (i=0; i< slv.size(); i++) slv[i] = toupper(slv[i]);

	
	string value = req.GetEntry("from").GetValue();
	unsigned from = 0;
	if ( !value.empty() ) from = atoi(value.c_str());

	value = req.GetEntry("to").GetValue();
	unsigned to = 0;
	if ( !value.empty() ) to = atoi(value.c_str());

	VastPageDataPtr vpp = NewDataType <VastPageData> (1);
  	unsigned row_count =
      		constructTopNBestAlignedVastRows(vpp,1, (char *)mst.c_str(),
					(char *)slv.c_str(), from, to, 0, 0,1);

  	if (!row_count) vpp[0].IpdpLen = 0;
  	BiostrucAnnotSetPtr pbsa = constructBASPFromVastPagePtr(vpp, 1);

{AsnIoPtr aip;
aip = AsnIoOpen("pbsa.out", "w");
BiostrucAnnotSetAsnWrite(pbsa, aip, NULL);
AsnIoClose(aip);
}

	auto_ptr <CObjectOStream> oos
	      ( CObjectOStream::Open( eSerial_AsnBinary, resp.out() ) ); //Cn3D
	CObjectOStream::AsnIo aip(*oos, "Biostruc-annot-set");

	//  printf("Content-type: chemical/ncbi-asn1-binary\n\n"); launch Cn3D 
	//  printf("Content-type: text/html\n\n"); for test

/*
   auto_ptr <CObjectOStream> oos
              ( CObjectOStream::Open( eSerial_AsnText, resp.out() ) );
   aip = AsnIoOpen("pbsa.out", "w");  
*/

	resp.SetContentType("application/octet-stream"); 	// Cn3D
	resp.WriteHeader();
  	BiostrucAnnotSetAsnWrite(pbsa, aip, NULL);
	aip.End();

  	BiostrucAnnotSetFree(pbsa);
	delete [] vpp;

  	return;
}
	
