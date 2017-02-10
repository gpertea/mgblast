/*  $Id: ncbi_priv.c,v 6.9 2005/10/25 18:53:10 lavr Exp $
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
 * Author:  Denis Vakatov
 *
 * File Description:
 *   Private aux. code for the "ncbi_*.[ch]"
 *
 */

#include "ncbi_priv.h"
#if defined(NCBI_OS_UNIX)
#  include <unistd.h>
#elif defined(NCBI_OS_MSWIN)
#  include <windows.h>
#else
#  include <connect/ncbi_socket.h>
#endif /*NCBI_OS_...*/
#include <stdarg.h>
#include <string.h>


int     g_NCBI_ConnectRandomSeed = 0;

MT_LOCK g_CORE_MT_Lock           = 0;
LOG     g_CORE_Log               = 0;
REG     g_CORE_Registry          = 0;


extern int g_NCBI_ConnectSrandAddend(void)
{
#if defined(NCBI_OS_UNIX)
    return (int) getpid(); 
#elif defined(NCBI_OS_MSWIN)
    return (int) GetCurrentProcessId();
#else
    return SOCK_gethostbyname(0);
#endif /*NCBI_OS_...*/ 
}


extern const char* g_CORE_Sprintf(const char* fmt, ...)
{
  va_list args;
  static char str[4096];

  va_start(args, fmt);
  *str = '\0';
  vsprintf(str, fmt, args);
  assert(strlen(str) < sizeof(str));
  va_end(args);
  return str;
}


extern char* g_CORE_RegistryGET
(const char* section,
 const char* name,
 char*       value,
 size_t      value_size,
 const char* def_value)
{
    char* ret_value;
    CORE_LOCK_READ;
    ret_value = REG_Get(g_CORE_Registry,
                        section, name, value, value_size, def_value);
    CORE_UNLOCK;
    return ret_value;
}


/*
 * ---------------------------------------------------------------------------
 * $Log: ncbi_priv.c,v $
 * Revision 6.9  2005/10/25 18:53:10  lavr
 * Fix ADDEND spelling (finally, hopefully)
 *
 * Revision 6.8  2005/07/11 18:14:14  lavr
 * Fix ADDEND spelling
 *
 * Revision 6.7  2005/05/03 13:56:40  lavr
 * +<connect/ncbi_socket.h> for non-UNIX, non-Windows platforms
 *
 * Revision 6.6  2005/05/03 11:50:19  ivanov
 * Added MS Win specific for NCBI_CONNECT_SRAND_ADDEND, removing dependency
 * from socket library.
 *
 * Revision 6.5  2005/05/02 16:04:20  lavr
 * Use global random seed
 *
 * Revision 6.4  2002/09/24 15:06:40  lavr
 * Log moved to end
 *
 * Revision 6.3  2002/06/18 18:39:38  ucko
 * Explicitly initialize global variables to avoid a MacOS X linker bug.
 *
 * Revision 6.2  2002/03/22 19:52:17  lavr
 * Do not include <stdio.h>: included from ncbi_util.h or ncbi_priv.h
 *
 * Revision 6.1  2000/03/24 22:53:35  vakatov
 * Initial revision
 *
 * ===========================================================================
 */
