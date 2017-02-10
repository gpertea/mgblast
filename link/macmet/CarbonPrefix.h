/*   carbonmacros.h
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
* File Name:  carbonmacros.h
*
* Author:  Joshua Juran
*
* Version Creation Date:   2001-03-21
*
* $Revision: 1.3 $
*
* File Description:
*		header for controlling the level of enforcement of Carbon-compliance
*		Mac only
*
* ==========================================================================
*/
#ifndef _CARBONPREFIX_
#define _CARBONPREFIX_

// Carbonization
// -------------

// Phase 1, 2, 3
# define ACCESSOR_CALLS_ARE_FUNCTIONS 1
# define OPAQUE_TOOLBOX_STRUCTS 1
# define TARGET_API_MAC_CARBON 1

#endif
