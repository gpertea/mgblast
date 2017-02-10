#ifndef TEST_ASSERT_IMPL__H
#define TEST_ASSERT_IMPL__H

/*  $Id: test_assert_impl.h,v 1.1 2006/03/07 18:22:09 lavr Exp $
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
 *   Setup #NDEBUG and #_DEBUG preprocessor macro in a way that ASSERTs
 *   will be active even in the "Release" mode (it's useful for test apps).
 *
 */

#ifndef TEST_ASSERT__H
#  error "Must not use this header alone, but from a proper wrapper."
#endif /*TEST_ASSERT__H*/

#if defined(NCBI_OS_MAC) || \
   (defined(NCBI_OS_DARWIN) && defined(NCBI_COMPILER_METROWERKS))
#  include <stdio.h>
#  include <stdlib.h>

#elif defined(NCBI_OS_MSWIN)
#  ifdef   _ASSERT
#    undef _ASSERT
#  endif
#  define  Type aType
#  include <crtdbg.h>
#  include <stdio.h>
#  include <windows.h>
#  undef   Type

/* Suppress popup messages on execution errors.
 * NOTE: Windows-specific, suppresses all error message boxes in both runtime
 * and in debug libraries, as well as all General Protection Fault messages.
 * Environment variable DIAG_SILENT_ABORT must be set to "Y" or "y".
 */

/* Handler for "Unhandled" exceptions */
static LONG CALLBACK _SEH_Handler(EXCEPTION_POINTERS* ep)
{
    /* Always terminate a program */
    return EXCEPTION_EXECUTE_HANDLER;
}

static int _SuppressDiagPopupMessages(void)
{
    /* Check environment variable for silent abort app at error */
    const char* value = getenv("DIAG_SILENT_ABORT");
    if (value  &&  (*value == 'Y'  ||  *value == 'y')) {
        /* Windows GPF errors */
        SetErrorMode(SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX |
                     SEM_NOOPENFILEERRORBOX);

        /* Runtime library */
        _set_error_mode(_OUT_TO_STDERR);

        /* Debug library */
        _CrtSetReportFile(_CRT_WARN,   _CRTDBG_FILE_STDERR);
        _CrtSetReportMode(_CRT_WARN,   _CRTDBG_MODE_FILE);
        _CrtSetReportFile(_CRT_ERROR,  _CRTDBG_FILE_STDERR);
        _CrtSetReportMode(_CRT_ERROR,  _CRTDBG_MODE_FILE);
        _CrtSetReportFile(_CRT_ASSERT, _CRTDBG_FILE_STDERR);
        _CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE);

        /* Exceptions */
        SetUnhandledExceptionFilter(_SEH_Handler);
    }
    return 0;
}

/* Put this function at startup init level 'V', far enough not to mess up with
 * base RTL init, which happens at preceding levels in alphabetical order.
 */
#  if _MSC_VER >= 1400
#    pragma section( ".CRT$XIV", read)
#  endif
#  pragma data_seg(".CRT$XIV")
static int (*_SDPM)(void) = _SuppressDiagPopupMessages;
#  pragma data_seg()

#endif /*defined(NCBI_OS_...)*/


#ifdef   NDEBUG
#  undef NDEBUG
#endif
#ifdef   assert
#  undef assert
#endif

/* IRIX stdlib fix (MIPSpro compiler tested): assert.h already included above*/
#ifdef NCBI_OS_IRIX
#  ifdef   __ASSERT_H__
#    undef __ASSERT_H__
#  endif
#endif

/* Likewise on OSF/1 (at least with GCC 3, but this never hurts) */
#ifdef NCBI_OS_OSF1
#  ifdef   _ASSERT_H_
#    undef _ASSERT_H_
#  endif
#endif

/* ...and on Darwin (at least with GCC 3, but this never hurts) */
#ifdef NCBI_OS_DARWIN
#  ifdef   FIXINC_BROKEN_ASSERT_STDLIB_CHECK
#    undef FIXINC_BROKEN_ASSERT_STDLIB_CHECK
#  endif
#endif

#include <assert.h>

#ifdef   _ASSERT
#  undef _ASSERT
#endif
#define  _ASSERT assert


/*
 * --------------------------------------------------------------------------
 * $Log: test_assert_impl.h,v $
 * Revision 1.1  2006/03/07 18:22:09  lavr
 * Initial revision
 *
 * ===========================================================================
 */

#endif  /* TEST_ASSERT_IMPL__H */
