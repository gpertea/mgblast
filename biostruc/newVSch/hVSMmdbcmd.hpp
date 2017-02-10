#ifndef _VSMmdbCMD_HPP
#define _VSMmdbCMD_HPP

/*
 * $Id: hVSMmdbcmd.hpp,v 1.1 2005/07/26 17:11:46 chenj Exp $
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
 * $Log: hVSMmdbcmd.hpp,v $
 * Revision 1.1  2005/07/26 17:11:46  chenj
 * Making linux VSMmdb.cgi
 *
 *
 *
 */


#include <cgi/cgictx.hpp>
#include <sstream>
#include "hVSMmdbres.hpp"
#include "qman.hpp"

#define PDBMAXLEN 8000000  // 8Mb

enum FileDispOp {

	LAUNCH_VIEWER= 0,
	SEE_FILE,
	SAVE_FILE
};


using namespace ncbi;
using namespace std;

BEGIN_NCBI_SCOPE


class DataInfo {

  public: 
        DataInfo () {};
        ~DataInfo() {};

        static 	string   	JobID;
        static 	string   	subset;
        static 	stringstream   	PDB_ios;
	static 	stringstream   	PDB_chk; 
        static 	QmReqID 	ReqID;
	static 	QmGrpID 	GrpID;
	static 	string		JobType;
};


class CVSMmdbCommand : public CNcbiCommand
{

 public:
      
	CVSMmdbCommand(CNcbiResource& resource);
	virtual ~CVSMmdbCommand(void) {};

	virtual void Execute(CCgiContext& ctx) = 0;
	virtual string GetLink(CCgiContext&) const
		{ return NcbiEmptyString; }

 protected:

	CVSMmdbResource& GetVSMmdbResource() const
	  { return dynamic_cast<CVSMmdbResource&> (GetResource()); }

	CVSMmdbCommand& operator = (const CVSMmdbCommand& rhs);

	virtual string GetEntry() const;
};



// CVSMmdbSubmitCommand

class CVSMmdbSubmitCommand : public CVSMmdbCommand
{
 public:

	CVSMmdbSubmitCommand(CNcbiResource& resource);
	virtual ~CVSMmdbSubmitCommand(void) {};

	virtual CNcbiCommand* Clone(void) const;
	virtual string GetName() const;
        virtual void Execute(CCgiContext& ctx);

 protected:

	string 	ReadReqAndPDBSubmission(CCgiContext& ctx);

};



class CVSMmdbStrTextCommand : public CVSMmdbCommand
{
 public:

        CVSMmdbStrTextCommand(CNcbiResource& resource);
        virtual ~CVSMmdbStrTextCommand(void) {};

        virtual CNcbiCommand* Clone(void) const;
        virtual string GetName() const;
        virtual void Execute(CCgiContext& ctx);
};


class CVSMmdbStrImgCommand : public CVSMmdbCommand
{
 public:

   	CVSMmdbStrImgCommand(CNcbiResource& resource);
	virtual ~CVSMmdbStrImgCommand(void) {};

	virtual CNcbiCommand* Clone(void) const;
	virtual string GetName() const;
        virtual void Execute(CCgiContext& ctx);

};




class CVSMmdbViewCommand : public CVSMmdbCommand
{
  public:

        CVSMmdbViewCommand(CNcbiResource& resource);
        virtual ~CVSMmdbViewCommand(void) {};

	virtual CNcbiCommand* Clone(void) const;
	virtual string GetName() const;
	virtual void Execute(CCgiContext& ctx); 

  protected:
	
	void SendStructureMIME(char Filetype,  int Mime, int Complexity,
        					int Models, CCgiContext& ctx);
};


class CVSMmdbSearchCommand : public CVSMmdbCommand
{
 public:
	
        CVSMmdbSearchCommand(CNcbiResource& resource);
        virtual ~CVSMmdbSearchCommand(void) {};

	virtual CNcbiCommand* Clone(void) const;
	virtual string GetName() const;
	virtual void Execute(CCgiContext& ctx);

};



class CVSMmdbChkErrCommand : public CVSMmdbCommand
{
 public:

        CVSMmdbChkErrCommand(CNcbiResource& resource);
        virtual ~CVSMmdbChkErrCommand(void) {};

        virtual CNcbiCommand* Clone(void) const;
        virtual string GetName() const;
        virtual void Execute(CCgiContext& ctx);

};

END_NCBI_SCOPE

#endif
