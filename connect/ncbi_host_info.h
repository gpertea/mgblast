#ifndef CONNECT___NCBI_HOST_INFO__H
#define CONNECT___NCBI_HOST_INFO__H

/*  $Id: ncbi_host_info.h,v 6.7 2006/03/06 20:23:59 lavr Exp $
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
 * Author:  Anton Lavrentiev
 *
 * File Description:
 *   NCBI host info getters
 *
 *   Host information handle becomes available from SERV_Get[Next]InfoEx()
 *   calls of the service mapper (ncbi_service.c) and remains valid until
 *   destructed by passing into free(). All API functions declared below
 *   accept NULL as 'host_info' parameter, and as the result return a failure
 *   status as described individually for each API call.
 *
 */

#include <connect/connect_export.h>


/** @addtogroup ServiceSupport
 *
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


struct SHostInfoTag;
typedef struct SHostInfoTag* HOST_INFO;


/* Return CPU count or -1 if error occurred.
 */
extern NCBI_XCONNECT_EXPORT int HINFO_CpuCount(const HOST_INFO host_info);


/* Return task count or -1 if error occurred.
 */
extern NCBI_XCONNECT_EXPORT int HINFO_TaskCount(const HOST_INFO host_info);


/* Return non-zero on success and store load averages in the
 * provided array "lavg", with the standard load average for last
 * minute stored at index [0], and instant load average
 * (aka BLAST) stored at index [1]. Return 0 on error.
 */
extern NCBI_XCONNECT_EXPORT int/*bool*/ HINFO_LoadAverage
(const HOST_INFO host_info,
 double          lavg[2]
 );


/* Return non-zero on success and store host status coefficients in
 * the provided array "status", with status based on the standard
 * load average stored at index [0], and that based on instant load
 * average stored at index [1]. Status may return as 0 if the host
 * does not provide such information. Return 0 on error.
 */
extern NCBI_XCONNECT_EXPORT int/*bool*/ HINFO_Status
(const HOST_INFO host_info,
 double          status[2]
 );


/* Obsolete.  Always returns 0 and does not touch its "blast" argument.
 */
extern NCBI_XCONNECT_EXPORT int/*bool*/ HINFO_BLASTParams
(const HOST_INFO host_info,
 unsigned int    blast[8]
 );


/* Obtain and return host environment. The host environment is the
 * sequence of lines (separated by \n), all having the form "name=value",
 * which are provided to load-balancing service mapping daemon (LBSMD)
 * in the configuration file on that host. Return 0 if the host
 * environment cannot be obtained. If completed successfully, the
 * host environment remains valid until the handle 'host_info' deleted
 * in the application program.
 */
extern NCBI_XCONNECT_EXPORT const char* HINFO_Environment
(const HOST_INFO host_info);


/* Obtain affinity argument and value that has keyed the service
 * selection (if affinities have been used at all).  NULL gets returned
 * as argument if no affinity has been found (in this case value
 * will be returned 0 as well).  Otherwise, NULL gets returned as
 * value if there was no particular value matched but the argument
 * played alone; "" is the value has been used empty, or any other
 * substring from the host environment that has keyed the decision.
 */
extern NCBI_XCONNECT_EXPORT const char* HINFO_AffinityArgument
(const HOST_INFO host_info);

extern NCBI_XCONNECT_EXPORT const char* HINFO_AffinityArgvalue
(const HOST_INFO host_info);


#ifdef __cplusplus
}  /* extern "C" */
#endif


/* @} */


/*
 * --------------------------------------------------------------------------
 * $Log: ncbi_host_info.h,v $
 * Revision 6.7  2006/03/06 20:23:59  lavr
 * Added "const" qualifier to all host-infos when passed to getters
 *
 * Revision 6.6  2006/03/05 17:33:15  lavr
 * +HINFO_AffinityArgument, +HINFO_AffinityArgvalue
 *
 * Revision 6.5  2003/04/09 19:05:42  siyan
 * Added doxygen support
 *
 * Revision 6.4  2003/02/08 21:03:51  lavr
 * Unimportant change in comments
 *
 * Revision 6.3  2003/01/08 01:59:32  lavr
 * DLL-ize CONNECT library for MSVC (add NCBI_XCONNECT_EXPORT)
 *
 * Revision 6.2  2002/11/08 17:16:11  lavr
 * NULL parameter acceptance explicitly stated
 *
 * Revision 6.1  2002/10/28 20:12:02  lavr
 * Initial revision
 *
 * ==========================================================================
 */

#endif /* CONNECT___NCBI_HOST_INFO__H */
