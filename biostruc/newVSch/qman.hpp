#pragma once
#ifndef qman_h_included
#define qman_h_included


#ifdef __cplusplus
extern "C" {
#endif


//  QmReqID is the main 64bit identifier for accounting job Request 
#ifdef WIN32 
    #define I64Format "%I64i" // use this for QmReqID printfs 
	typedef __int64 QmReqID;
	typedef __int64 QmJobID;
	typedef __int64 QmGrpID;
#else 
    #define I64Format "%llu" // use this for QmReqID printfs 
    typedef long long QmReqID;
	typedef long long QmJobID;
	typedef long long QmGrpID;
#endif

		// this Handler defines the interface with database : like FILE * for file stream functions
typedef void * QmHandle;
extern char qmDebugLog[];

		// initialize the DB connection with server 
		// if(server==0) default will be used : DDDSQL5
QmHandle qmOpen( const char * server); 
		// free the database resources associated with QmHandle 
void qmClose(QmHandle handle); 

	//  submit the job
QmReqID qmSubmit(QmHandle handle, // the handle you got while opening DB connection
		// the name of the service you want to run: pcngbr, omssa, ...
	const char * service, 
		// set of CGI style parameters to form 
		// command line parameters for service executeable 
	const char * params, 
		// the blobs to be associated with this job request
		// pay attention that this variable amount argument function 
		// allows you to associate multiple blobs with a single job request
		// to end the list terminate it with dataname=0
	//const char * datanam, // name of the blob 
	//const void * datablob, // data of the blob 
	//unsigned int datasize, // the size of the blob
	... ); // terminate by 0 after lastblob definition

// qmSubmit usage samples :
//
// TO submit a job without blobs associated 
//		QmReqID req=qmSubmit(hdl, "zip", "par1=value1&par2=value2", 0 ); 
//
// TO submit a job with two blobs associated
//		QmReqID req=qmSubmit(hdl, "zip", "par=value"
//			, "firstfile" , file1Buf, file1Size 
//			, "secondfile" , file2Buf, file2Size, 
//			0 ); // zero is the termination
//
// if a blob is associated with a job during the submission 
// the request will be put to the status qmReqActionRun, which tells the job 
// handler to execute this job. If job has no blobs associated during submission
// the request stays in Wait(qmReqActionNone) status expecting you to add data.
// You can use qmSetIOData to add data to that job later. After you are 
// finished adding IOData  or you want to run the job wihtout associated data:
// use qmSetStatus(qmReqActionRun)
//


// associate IO data with request ID 
bool qmSetIOData( QmHandle handle, 
	QmReqID req, // the request Identifier
	const char * dataname,  // data blob name 
	const void * datablob , // data blob itself
	unsigned int datasize );// the size of the blob

// get associated data blob 
// NOTE: the caller function is responsible for freeing the pointer 
void * qmGetIOData( QmHandle handle, 
	QmReqID req, // the request Id
	const char * dataname, // the name of datablob requested
	unsigned int * datasize ); // pointer to return the size of the data blob, may be 0

// service specific infromation
typedef struct tagQmSrv{
	char name[64]; // the name of the service 
	char cmdLine [1024] ; // service command line
	int maxJobs; // the maximum number of jobs for this type of service
	unsigned int jobType; // the default type of job for this service
	int srvStatus; // service status 
	unsigned int ncpu; // the number of days data and requests are alive for this kind of job
//	int cntJobs; // the number of jobs running now 
	char queues[1024];// the computers running this service 
	unsigned int maxtry;  // the maximum number of tries to launch jobs of this service
}QmSrv;
unsigned int qmGetSrv(QmHandle handle, const char *wheresrvname, QmSrv * srv); // get the service information structure for the service satisfying the where clause 
void qmSetSrvStatus(QmHandle handle, const char * srvname, int srvStatus); // set the service status for the given service


// reqId-group specific functions
QmGrpID qmSetGroup(QmHandle handle, QmReqID req, QmGrpID grp ); // creates/associates group with the given req Id
unsigned int qmGetGroupGrps(QmHandle handle, QmReqID req, QmGrpID * grps, unsigned int maxCnt); // retrieves the group ID byt req ID
unsigned int qmGetGroupReqs(QmHandle handle, QmGrpID grp, QmReqID * reqs, unsigned int maxCnt); // retrieves the list of reqIDs by the group ID



// job specific 
QmReqID qmGrabAJobToRun(QmHandle handle, const char * whereclause); // grab a job satisfying the whereclause
bool qmRegisterServiceJob(QmHandle handle, const char * srvname, const char * hostname, int jobid); // regiseter a jobid

// manipulate errors 
void qmSetError( QmHandle handle, QmReqID req, const char * error); 
char * qmGetError( QmHandle handle, QmReqID req); // NOTE: the caller function is responsible for freeing the pointer 

// manipulate logs 
void qmSetLog( QmHandle handle, QmReqID req, const char * log); 
char * qmGetLog( QmHandle handle, QmReqID req); // NOTE: the caller function is responsible for freeing the pointer 

// manipulate config parameters
void qmSetConfig( QmHandle handle, const char * nam, const char * value);
char * qmGetConfig( QmHandle handle, const char * nam);

// jobids 
void qmSetJobId( QmHandle handle, QmReqID req, QmJobID jobid);
QmJobID qmGetJobId( QmHandle handle, QmReqID req);

//emails per job 
void qmSetEmail( QmHandle handle, QmReqID req, const char * email);
char * qmGetEmail( QmHandle handle, QmReqID req);

// ips per job 
void qmSetIP( QmHandle handle, QmReqID req, const char * email);
char * qmGetIP( QmHandle handle, QmReqID req);


typedef enum {
        qmReqTypeNormal=0, // normal job to be executed whenever possible 
		qmReqTypeFast, // fast job to be executed as soon as possible
		qmReqTypeFastService, // constantly running service : to be kept alive 
		qmReqTypePullService, // constantly running service : to be kept alive on an idividual machine
}QmReqType;
void qmSetType( QmHandle handle, QmReqID req, QmReqType type);
QmReqType qmGetType( QmHandle handle, QmReqID req);


// manipulate the status 
typedef enum {
        qmReqStatusWaiting=0, // the inital status  ... nothing yet happened to this request
		qmReqStatusProcessing, // job was grabbed but not yet running 
        qmReqStatusRunning, // the job was submitted 
        qmReqStatusSuspended, // the job was submitted either by the user or by the sysyem 
        qmReqStatusDone, // the job has finished 
		qmReqStatusKilled, // job was killed by the request of user 
        qmReqStatusProgError, // programatic error, job executeable reported an error
        qmReqStatusSysError, // system failure, infrastructure has had a problem
		qmReqStatusUnknown=255, // unknown status: probably wrong ReqID
}QmReqStatus;
void qmSetStatus( QmHandle handle, QmReqID req, QmReqStatus status );
QmReqStatus qmGetStatus( QmHandle handle, QmReqID req);
const char * qmGetStatusText (QmReqStatus stat);


// manipulate the action to be performed on this Request
typedef enum { // types of actions 
	qmReqActionNone=0,  // initial state ... job is in the table but not to be submitted
	qmReqActionRun, // tell job handler to submit this job
	qmReqActionKill,  // tell job handler to kill this job
	qmReqActionSuspend, //  ... to suspend 
	qmReqActionResume  // ... to resumed previously suspended job 
} QmReqAction; 
void qmSetAction( QmHandle handle, QmReqID req, QmReqAction action);
QmReqAction qmGetAction( QmHandle handle, QmReqID req);
const char * qmGetActionText(QmReqAction  action);


// information
const char * qmGetReqTime( QmHandle handle, QmReqID req, char * timeBuf, unsigned int sizeTimeBuf);
void qmTransactionStart(QmHandle handle, const char * tName);
void qmTransactionCommit(QmHandle handle, const char * tName);




/////////////////////////////////////////////////////////////
// 
// Static DB handle using version of the same functions 
// User doesn't need to pass qmHandle for these functions
//
bool sqmInit( void ); 
bool sqmFree(void ); 

QmReqID sqmSubmit(const char * srvname, const char * params, ... );

bool sqmSetIOData (QmReqID req,	const char * dataname, const void * datablob , unsigned int datasize );
void * sqmGetIOData (QmReqID req, const char * dataname, unsigned int * datasize );

QmGrpID sqmSetGroup(QmReqID req, QmGrpID grp);
unsigned int sqmGetGroupGrps(QmReqID req, QmGrpID * grp, unsigned int maxCnt);
unsigned int sqmGetGroupReqs(QmGrpID grp, QmReqID * req, unsigned int maxCnt);

void sqmSetError( QmReqID req, const char * error);
char * sqmGetError ( QmReqID req);
void sqmSetLog( QmReqID req, const char * log) ;
char * sqmGetLog( QmReqID req);
void sqmSetEmail( QmReqID req, const char * log) ;
char * sqmGetEmail( QmReqID req);
void sqmSetIP( QmReqID req, const char * ip) ;
char * sqmGetIP( QmReqID req);
void sqmSetConfig( const char * nam, const char * value);
char * sqmGetConfig( const char * nam);

void sqmSetType( QmReqID req, QmReqType type);
QmReqType sqmGetType( QmReqID req);
void sqmSetStatus (QmReqID req, QmReqStatus status ) ;
void sqmSetJobId( QmReqID req, QmJobID jobid);
QmJobID sqmGetJobId( QmReqID req);
QmReqStatus sqmGetStatus( QmReqID req);
void sqmSetAction( QmReqID req, QmReqAction action) ;
QmReqAction sqmGetAction ( QmReqID req ) ;

const char * sqmGetReqTime( QmReqID req);
void sqmTransactionStart(const char * tName);
void sqmTransactionCommit(const char * tName);




#ifdef __cplusplus
}

#endif


#endif





